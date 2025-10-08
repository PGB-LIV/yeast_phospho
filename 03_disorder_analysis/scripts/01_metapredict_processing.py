########################
# Disorder analysis    #
# YEAST PHOSPHO        #  
########################
import pandas as pd
import Bio.SeqIO as SeqIO
import os
##########################################
# 1. Get max seq lenghth                 #
# Needed for reading in df correctly     # 
# with column header for index of score  #
##########################################

# create results dir
os.makedirs("metapredict_disorder_processed_v3", exist_ok=True)
# get max seq length
seq_len = pd.read_csv("yeast_metapredict_v3_270625.csv", usecols=[0,1], header=None)
seq_len.columns =["protein","seq"]
# rm spaces
seq_len['seq'] = seq_len['seq'].str.strip()
# col with seq length
seq_len['seq_len'] = seq_len['seq'].str.len()
seq_len.to_csv("metapredict_disorder_processed_v3/seqlen.csv")
# Longest seq (ie column header needs to go up to this number)
max_seq_length = seq_len['seq_len'].max()
print("Max sequence length:", max_seq_length)

# Column names
col_list = ["Protein","Sequence"]
# col names of index pos
pos_list = list(range(1,max_seq_length+1))
col_list = col_list + pos_list

###############################
# 2. Now read in full df with #
# column names as above       #   
###############################

# Read in data with specificied col names
df = pd.read_csv("yeast_metapredict_v3_270625.csv",header=None,names=col_list) # From metapredict (on cluster), need to gives names due to different lengths of rows
# remove spaces from seq column
df['Sequence'] = df['Sequence'].str.strip()
df[["Protein","Descrip"]] = df["Protein"].str.split(" ", n=1,expand=True)
df = df.drop('Descrip', axis=1)
# Remove decoy and contm
df = df[~df.Protein.str.contains("rev_")]
contaminants_file = open("cRAP_contaminants.txt","r")
contaminants_list = []
# Remove target species from contam list
for line in contaminants_file:
    row = line.strip()
    contaminants_list.append(row)
# Filter out contaminants from the target species
non_target_contaminants = [contaminants for contaminants in contaminants_list if "YEAST" not in contaminants]
print(non_target_contaminants)
# Contaminant prefix
contam = "CONTAM_" 
df = df.reset_index(drop=True)
# add contam prefix
df["Protein"] = df.apply(
    lambda row: contam + row["Protein"] if any(contaminant in row["Protein"] for contaminant in non_target_contaminants) else row["Protein"],
    axis=1)
df.to_csv("metapredict_disorder_processed_v3/contam.csv")
# Remove contam
df=df[~df.Protein.str.contains(contam)] 
df.to_csv("metapredict_disorder_processed_v3/testing_contam.csv")

# Wide -> long
df = pd.melt(df, id_vars=["Protein","Sequence"])
df.to_csv("metapredict_disorder_processed_v3/melt.csv")
# remove blanks (ie proteins with shorter seq)
df = df.dropna(subset=["value"])  
df = df.rename(columns={"variable":"Protein_pos","value":"disorder"})

# get aa corresponding to row
df["PTM_residue"] = df.apply(
    lambda row: row["Sequence"][row["Protein_pos"] - 1] if pd.notnull(row["Sequence"]) else None,
    axis=1
)
# STY only
df = df[(df["PTM_residue"]=="S")|(df["PTM_residue"]=="T")|(df["PTM_residue"]=="Y")]
         
# add col to say if aa is a phosphosite in gsb
gsb = pd.read_csv("/mnt/hc-storage/users/hleboswe/Yeast_GSB_Nov2024/All_site_formats_Updated/G2S1B_0.05_protein_pos_single_prot_mapping.csv")
gsb = gsb[(gsb["PTM_residue"]=="S")|(gsb["PTM_residue"]=="T")|(gsb["PTM_residue"]=="Y")]
gsb = gsb[["Protein","Protein_pos", "PTM_residue", "PTM_FLR_category"]]
gsb["Phospho"] = True
df = df.merge(gsb,how="left", on=["Protein", "Protein_pos","PTM_residue"])
df=df[["Protein","Protein_pos", "disorder", "PTM_residue","PTM_FLR_category","Phospho"]]
df.to_csv("metapredict_disorder_processed_v3/disorder_processed.csv",index=False)

