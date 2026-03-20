import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
# import requests
# import json
import sys
import time
import random
import itertools
"""
This script is for testing whether we have reached saturation of phosphosites - adaptation of GSB script
"""

# usage: python Site_based_format_GSB_counts_pipeline.py file_list.txt meta.csv PRIDE/SDRFs/ NA 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]
# OR
# python Site_based_format_GSB_counts_pipeline.py file_list.txt NA NA simple_meta/ 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]
# OR
# python Site_based_format_GSB_counts_pipeline.py file_list.txt NA NA NA 5 2 [optional: DECOY_Prefix ie. DECOY] [optional: CONTAM_prefix ie. CONTAM] [optional: modification:target:decoy ie. phospho:STY:A]

# read in csv from FLR pipeline
# DECOY and CONTAM removed, filtered for contains Phospho, exploded for site based, binomial adjustment and collapsed by peptidoform

# TO DO: probably need some error handling here, to check using the right params
folder_list_file = open(sys.argv[1], "r")
folder_list_file = folder_list_file.read()
folder_list = folder_list_file.replace('\n', ';').split(";")

dataset_list = []
for i in folder_list:
    index = [idx for idx, s in enumerate(i.split("/")) if 'PXD' in s][0]
    if i.split("/")[index] not in dataset_list:
        dataset_list.append(i.split("/")[index])


gold_count = int(sys.argv[5])
silver_count = int(sys.argv[6])

if len(sys.argv) > 7:
    decoy_prefix = sys.argv[7]
    contam_prefix = sys.argv[8]
else:
    decoy_prefix = "DECOY"
    contam_prefix = "CONTAM"
if len(sys.argv) > 9:
    search_mod = sys.argv[9].split(":")[0]
    target_aa = sys.argv[9].split(":")[1]
    decoy_aa = sys.argv[9].split(":")[2]
else:
    decoy_aa = "A"
    search_mod = "Phospho"
    target_aa = "STY"

print(f"Using decoy prefix: {decoy_prefix}")
print(f"Using contam prefix: {contam_prefix}")
print(f"Using modification: {search_mod}")
print(f"Using target amino acids: {target_aa}")
print(f"Using decoy amino acid: {decoy_aa}")


def ratio(df, targets, decoy):
    STY_count = 0
    for target in list(targets):
        T_count = df['Peptide'].str.count(target).sum()
        STY_count += T_count
    A_count = df['Peptide'].str.count(decoy).sum()
    STY_A_ratio = STY_count / A_count
    return STY_A_ratio


print("Gold threshold: " + str(gold_count) + "\nSilver threshold: " + str(silver_count))
flr_filter = 0.05

df = pd.read_csv("/mnt/hc-storage/users/hleboswe/Yeast_GSB_Nov2024/All_site_formats_Updated/all_datasets_merged_Site_Peptidoform_centric.tsv", sep="\t")
# Calculate STY:A ratio from the PSM centric format - for calculating FLR at GSB levels
STY_ratio = ratio(df, target_aa, decoy_aa)
print(f"Using target:decoy ratio of: {STY_ratio}")


# generate random sample of data sets 3x 7, 6 and 5 datasets
random.seed(1234)
seven_sample = random.sample(dataset_list, 3)
print(seven_sample)

dataset_list_7_1 = [i for i in dataset_list if i != seven_sample[0]]
print(dataset_list_7_1)
dataset_list_7_2 = [i for i in dataset_list if i != seven_sample[1]]
print(dataset_list_7_2)
dataset_list_7_3 = [i for i in dataset_list if i != seven_sample[2]]
print(dataset_list_7_3)
print(dataset_list)

six_sample = random.sample(dataset_list, 6)
print(six_sample)
dataset_list_6_1 = [i for i in dataset_list if i not in six_sample[0:2]]
print(dataset_list_6_1)
dataset_list_6_2 = [i for i in dataset_list if i not in six_sample[2:4]]
print(dataset_list_6_2)
dataset_list_6_3 = [i for i in dataset_list if i not in six_sample[4:6]]
print(dataset_list_6_3)

random.shuffle(dataset_list)
print(dataset_list)
dataset_list_5_1 = dataset_list[0:5]
print(dataset_list_5_1)
dataset_list_5_2 = dataset_list[1:3] + dataset_list[5:8]
print(dataset_list_5_2)
dataset_list_5_3 = dataset_list[3:8]
print(dataset_list_5_3)

overall_dataset_samples = [dataset_list_7_1, dataset_list_7_2, dataset_list_7_3, dataset_list_6_1, dataset_list_6_2, dataset_list_6_3, dataset_list_5_1, dataset_list_5_2, dataset_list_5_3]
overall_dataset_name = ["sample_7_1", "sample_7_2", "sample_7_3", "sample_6_1", "sample_6_2", "sample_6_3", "sample_5_1", "sample_5_2", "sample_5_3"]


### GSB counts
# filter for sites seen multiple times across dataset studies
for random_sample, sample_name in zip(overall_dataset_samples, overall_dataset_name):
    print(random_sample)
    for m in ["single", "all"]:
        print(f"Creating GSB format mapped to {m} proteins")
        for dataset in random_sample:
            loc_full = "/mnt/hc-storage/users/hleboswe/Yeast_GSB_Nov2024/All_site_formats_Updated/" + dataset + "_merged_Site_Peptidoform_centric.tsv"
            df = pd.read_csv(loc_full, sep="\t")
            df = df.loc[df['p' + decoy_aa + '_q_value_BA'] <= flr_filter]
            PSM_threshold = dict(df.groupby('Peptide_mod_pos')['0.05FLR_threshold_count'].sum())

            df = df.sort_values(['Peptide_mod_pos', 'p' + decoy_aa + '_q_value_BA', 'Binomial_final_score'],
                                     ascending=[True, False, True])
            df = df.drop_duplicates(subset=('Peptide_mod_pos'), keep="last", inplace=False)

            df['0.05FLR_threshold_count'] = df['Peptide_mod_pos'].map(PSM_threshold)

            if m == "all":
                if len(df) >= 1:
                    df['PTM_residue'] = df.apply(lambda x: x['Peptide'][x['PTM positions'] - 1], axis=1)
                    df['All_PTM_protein_positions'] = df['All_PTM_protein_positions'].str.split(':')
                    df['All_Proteins'] = df['All_Proteins'].str.split(':')
                    df = df.explode(['All_PTM_protein_positions', 'All_Proteins'])
                    df["Protein-pos"] = df['All_Proteins'] + "-" + df['All_PTM_protein_positions']
                else:
                    df['PTM_residue'] = ""
            else:
                if len(df) >= 1:
                    df['PTM_residue'] = df.apply(lambda x: x['Peptide'][x['PTM positions'] - 1], axis=1)
                    df['Protein-pos'] = df['Protein'] + "-" + df['Protein position'].astype(str)
                else:
                    df['PTM_residue'] = ""
            df['PTM_residue'] = np.where(df['PTM positions'] == 0, "N-term", df['PTM_residue'])
            df = df.sort_values(['Protein-pos', 'p' + decoy_aa + '_q_value_BA', 'Binomial_final_score'],
                                          ascending=[True, False, True])

            df['Protein_pos_res'] = df['Protein-pos'] + "_" + df['PTM_residue']

            PSM_threshold_2 = dict(df.groupby('Protein_pos_res')['0.05FLR_threshold_count'].sum())

            df = df.drop_duplicates(subset=('Protein_pos_res'), keep="last", inplace=False)

            # column for all PSM counts at 5%FLR
            df['Sum_of_PSM_counts(5%FLR)'] = df['Protein_pos_res'].map(PSM_threshold_2)

            df['PXD'] = dataset
            if dataset == random_sample[0]:
                df_counts = df
            else:
                df_counts = pd.concat([df_counts, df])

            df_temp = df[
                ['Peptide_mod_pos', 'p' + decoy_aa + '_q_value_BA', 'Binomial_final_score', 'Protein_pos_res',
                 'Sum_of_PSM_counts(5%FLR)', '0.05FLR_threshold_count']]
            df_temp.rename(columns={'p' + decoy_aa + '_q_value_BA': dataset + "_FLR",
                                    'Binomial_final_score': dataset + '_BinomialScore',
                                    'Peptide_mod_pos': dataset + "_peptide_mod_pos",
                                    'Sum_of_PSM_counts(5%FLR)': dataset + "_Sum_of_PSM_counts(5%FLR)",
                                    '0.05FLR_threshold_count': dataset + "_peptidoform_PSMcount(5%FLR)"}, inplace=True)
            df_temp = df_temp.set_index('Protein_pos_res')

            if dataset == random_sample[0]:
                df_final = df_temp
            else:
                df_final = pd.concat([df_final, df_temp], axis=1)

        column_list = [PXD + "_Sum_of_PSM_counts(5%FLR)" for PXD in random_sample]
        df_final["Sum_of_PSM_counts(5%FLR)"] = df_final[column_list].sum(axis=1)

        cols = df_final.columns.tolist()
        cols = [cols[-1]] + cols[:-1]
        df_final = df_final[cols]
        df_final = df_final.drop(column_list, axis=1)

        df_final = df_final[~df_final.index.str.contains(decoy_prefix)]
        df_final = df_final[~df_final.index.str.contains(contam_prefix)]

        df_counts['PTM_residue'] = df_counts['Protein_pos_res'].str.rsplit("_", n=1).str[-1]
        # added 130326 due to error in cross tab (seems to be from the index in concat step)
        df_counts = df_counts.reset_index(drop=True)

        counts_res = pd.crosstab(df_counts['PTM_residue'], df_counts['PXD']).replace(0,
                                                                                     np.nan).stack().reset_index().rename(
            columns={0: 'Count'})
        print(counts_res)
        counts_res.to_csv("/mnt/hc-storage/users/hleboswe/yeast_GSB_per_build_size/outputs/G" + str(gold_count) + "S" + str(silver_count) + "B_" + str(
            flr_filter) + "_Residue_counts_" + m+ "_" + sample_name + ".csv", index=False)

        for i in df_final.index.values.tolist():
            df_final.loc[i, 'Protein'] = i.rsplit("-", 1)[0]
            df_final.loc[i, 'Protein_pos'] = i.rsplit("-", 1)[1].split("_")[0]
            df_final.loc[i, 'PTM_residue'] = i.rsplit("_", 1)[-1]
            flr1_count = 0
            for dataset in random_sample:
                FLR_count_all = 0
                FLR_cols = [x for x in random_sample if dataset in x]
                for FLR_col in FLR_cols:
                    if df_final.loc[i, FLR_col + "_FLR"] != "N/A":
                        if float(df_final.loc[i, FLR_col + "_FLR"]) <= 0.01:
                            FLR_count_all += 1
                if FLR_count_all != 0:
                    flr1_count += 1
            if flr1_count >= gold_count:
                df_final.loc[i, 'PTM_FLR_category'] = "Gold"
            elif flr1_count >= silver_count:
                df_final.loc[i, 'PTM_FLR_category'] = "Silver"
            else:
                df_final.loc[i, 'PTM_FLR_category'] = "Bronze"

        cols = list(df_final.columns.values)
        cols.pop(cols.index('Protein'))
        cols.pop(cols.index('Protein_pos'))
        cols.pop(cols.index('PTM_residue'))
        cols.pop(cols.index('PTM_FLR_category'))

        df_final['Decoy_mod'] = np.where(df_final['PTM_residue'] == decoy_aa, 1, 0)

        df_final = df_final[['Protein', 'Protein_pos', 'PTM_residue', 'Decoy_mod', 'PTM_FLR_category'] + cols]

        # replace 0 PSM counts -> 1
        column_list = [PXD + "_peptidoform_PSMcount(5%FLR)" for PXD in random_sample]
        column_list.append("Sum_of_PSM_counts(5%FLR)")
        df_final[column_list] = df_final[column_list].replace(0, 1)

        df_final.to_csv("/mnt/hc-storage/users/hleboswe/yeast_GSB_per_build_size/outputs/G" + str(gold_count) + "S" + str(silver_count) + "B_" + str(
            flr_filter) + "_protein_pos" +"_"+ m +"_prot_mapping_" + sample_name + ".csv", index=False)

        counts = pd.crosstab(df_final.PTM_residue, df_final.PTM_FLR_category).replace(0,
                                                                                      np.nan).stack().reset_index().rename(
            columns={0: 'Count'})

        # calculate FLR estimate using STY:A ratio from PSM centric format
        f = lambda x, y, z: np.nan if x != decoy_aa else str(
            round(z * STY_ratio / counts.loc[counts['PTM_FLR_category'] == y, 'Count'].sum() * 100, 2)) + "%"

        counts['FLR'] = counts.apply(lambda x: f(x.PTM_residue, x.PTM_FLR_category, x.Count), axis=1)
        counts['FLR'] = counts.groupby('PTM_FLR_category')['FLR'].ffill()
        counts['PTM_FLR_category_FLR'] = counts['PTM_FLR_category'] + " - " + counts['FLR']
        counts['PTM_FLR_category'] = pd.Categorical(counts['PTM_FLR_category'], ["Bronze", "Silver", "Gold"])
        counts = counts.sort_values("PTM_FLR_category")
        print(counts)

        fig, axes = plt.subplots(ncols=3)
        for i, (name, group) in enumerate(counts.groupby("PTM_FLR_category_FLR", sort=False)):
            axes[i].set_title(name)
            group.plot(kind="bar", x="PTM_residue", y="Count", ax=axes[i], legend=False)
            axes[i].set_ylabel("count")
            axes[i].set_xlabel("")
        plt.tight_layout()
        # plt.show()
        plt.savefig("/mnt/hc-storage/users/hleboswe/yeast_GSB_per_build_size/outputs/G" + str(gold_count) + "S" + str(silver_count) + "B_" + str(
            flr_filter) + "_protein_pos_categories" +"_"+ m +"_prot_mapping_" + sample_name + ".png", dpi=300)
