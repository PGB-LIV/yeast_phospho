import pandas as pd
import os
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
#https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ

###########################################################################
# This script is to extract secondary structural info from AF3 files      #  
# after running dssp.                                                     # 
# The aim is to compare models with Gold PTM sites versus no PTM sites.   #
###########################################################################

os.makedirs("AF_model_dssp_processed_results", exist_ok=True)
os.makedirs("dssp_summary_files_v3", exist_ok=True)
############################# 
# multi GSB mapping file    #
#                           #  
#############################

GSB = pd.read_csv("/mnt/hc-storage/users/hleboswe/Yeast_GSB_Nov2024/All_site_formats_Updated/G2S1B_0.05_protein_pos_all_prot_mapping.csv")
# In the AF models we only use Gold sites in the PTM version
GSB_gold = GSB[GSB['PTM_FLR_category'] == 'Gold']
# Remove decoy hits
GSB_gold = GSB_gold[GSB_gold['PTM_residue'] != "A"]
# Protein ID (same format as used elsewhere for later merging)
GSB_gold["Protein"] = GSB_gold["Protein"].str.split("|").str[1]
# Dictionart with protein ID as key and PTM positions as values
gold_phosphosites = GSB_gold.groupby('Protein')['Protein_pos'].apply(list).to_dict()
#print(gold_phosphosites)

#######################################
# .cif files of AF models from dssp   #
#                                     # 
#######################################
counter = 0
overall_protein_structure = []
protein_ptm_results = []
protein_confidence_df_list = []
# loop through each file (modelled with and without PTMs for each protein)
for file_name in os.listdir(): 
    # get dssp results only
    if file_name.endswith("_dssp.cif"):
        counter += 1
        # obtain protein id
        #file = file_path.split("/")[-1]
        #protein = file_path.split("_sp_")[1]
        protein_id = file_name.split("_")[0].upper()
        # extract gold PTM positions for protein
        ptm_pos = gold_phosphosites[protein_id]
        # parse dssp file into dict
        mmcif_dict = MMCIF2Dict(file_name)
        mmcif_dict_struct_filtered = {}
        mmcif_dict_confidence_filtered = {}
        
        # loop through dict
        for key, value in mmcif_dict.items():
            # if the key is related to structure
            if "_struct_conf." in key:
                # add to new dict
                mmcif_dict_struct_filtered[key] = value
            elif "_atom_site." in key:
                # add to new dict
                mmcif_dict_confidence_filtered[key] = value 

        ########################
        # confidence score df  #
        ########################
        
        # convert to df with confidence score           
        confidence_df = pd.DataFrame.from_dict(mmcif_dict_confidence_filtered)
        # convert columns to numeric type
        confidence_df['_atom_site.auth_seq_id'] = pd.to_numeric(confidence_df['_atom_site.auth_seq_id'], errors='coerce')
        # Keep STY only            
        confidence_df = confidence_df[confidence_df['_atom_site.label_comp_id'].isin(['SER', 'SEP', 'THR', 'TPO', 'TYR', 'PTR'])]
        # PTM positions - check it if it a match
        confidence_df["phospho"] = confidence_df["_atom_site.auth_seq_id"].isin(ptm_pos)
        # Column with the list of PTM positions for protein
        confidence_df["PTM_position_all"] = ", ".join(map(str,ptm_pos))
        # to numeric (needed to calc mean)
        confidence_df['_atom_site.B_iso_or_equiv'] = pd.to_numeric(confidence_df['_atom_site.B_iso_or_equiv'], errors='coerce')
        # calculate mean confidence score for each amino acid
        confidence_mean = confidence_df.groupby(["_atom_site.label_seq_id"])['_atom_site.B_iso_or_equiv'].mean()
        confidence_df = confidence_df.set_index(['_atom_site.label_seq_id'])
        # add the score back to df
        confidence_df["mean_confidence_score"] = confidence_mean
        # Add column with protein id
        confidence_df["Protein"] = protein_id
        # add column for PTM model
        confidence_df['PTM_model'] = "no_PTM" if "non_phospho" in file_name else "PTM"
        # keep 1 row per aa - average per aa rather than atom
        confidence_df = confidence_df.drop_duplicates('_atom_site.auth_seq_id', keep='last')
        # collect results per file
        protein_confidence_df_list.append(confidence_df)

        #################
        # structure df  # 
        #################
        
        # convert to df - seoncdary structure          
        df=pd.DataFrame.from_dict(mmcif_dict_struct_filtered)
        # convert columns to numeric type
        df['_struct_conf.beg_label_seq_id'] = pd.to_numeric(df['_struct_conf.beg_label_seq_id'], errors='coerce')
        df['_struct_conf.end_label_seq_id'] = pd.to_numeric(df['_struct_conf.end_label_seq_id'], errors='coerce')
        # loop through ptm positions in list for protein. If the position occurs within the range of the struct feature, put the position in the col
        df['PTM_position'] = df.apply(
        lambda row: [position for position in ptm_pos 
                     if position in range(row['_struct_conf.beg_label_seq_id'], row['_struct_conf.end_label_seq_id'] + 1)],axis=1)
        # Convert to string
        df['PTM_position'] = df['PTM_position'].apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)
        # Column with the list of PTM positions for protein
        df["PTM_position_all"] = ", ".join(map(str,ptm_pos))
        if counter < 5:
            # filtered dssp file per protein model
            df.to_csv("AF_model_dssp_processed_results/"+(file_name[:-4])+".csv") 
        
        
        #####################################################
        # Extract the struct info for downstream analysis   # 
        #####################################################
        
        # Keep rows that map to a PTM position only
        #df2 = df[df['PTM_position'] != '']
        # Column with protein id (for merging all protein models later)
        df["Protein"] = protein_id
        # Column to distinguish whether PTMs where included in AF model (PTM AF model should always have "phospho" in file name)
        df['PTM_model'] = "no_PTM" if "non_phospho" in file_name else "PTM"
        # Filter columns
        df2 = df[["Protein","PTM_position","PTM_position_all", "PTM_model", "_struct_conf.conf_type_id","_struct_conf.beg_label_seq_id","_struct_conf.end_label_seq_id"]]
        # collect df for each protein in a list (to combine after loop) - struct at site
        protein_ptm_results.append(df2)
       
        
        # Extract number of aa contained in row for structural feature
        df["aa_total"] = df.apply(
        lambda row: row['_struct_conf.end_label_seq_id'] - row['_struct_conf.beg_label_seq_id'] + 1,
        axis=1)
        # count of aa in each structural feature category
        structural_groups = df.groupby('_struct_conf.conf_type_id')['aa_total'].sum().reset_index()
        # protein id
        structural_groups["Protein"] = protein_id
        # rename col
        structural_groups = structural_groups.rename(columns={"_struct_conf.conf_type_id" : "structure"})
        # ptm model
        structural_groups['PTM_model'] = "no_PTM" if "non_phospho" in file_name else "PTM"
        # collect df for each protein in a list (to combine after loop) - overall struct
        overall_protein_structure.append(structural_groups)



########################################################
# Overall output 1                                     # 
# PTM position level                                   #  
# Is there a change between PTM and no PTM model?      # 
########################################################

# combine all outputs from models 
final_ptm_pos = pd.concat(protein_ptm_results, ignore_index=True)
# Combine start and end amino acid number into a single column
final_ptm_pos["structure_start_end"] = final_ptm_pos["_struct_conf.beg_label_seq_id"].astype(str)+"-"+final_ptm_pos["_struct_conf.end_label_seq_id"].astype(str)
# filter for phospho
final_ptm_pos_phospho = final_ptm_pos[final_ptm_pos['PTM_position'] != ""]

# String to list
final_ptm_pos_phospho["PTM_position"] = final_ptm_pos_phospho["PTM_position"].str.split(", ")
#final_ptm_pos = final_ptm_pos.drop(columns=["_struct_conf.beg_label_seq_id", "_struct_conf.end_label_seq_id"])
# expand so each position is on its own row
final_ptm_pos_exploded = final_ptm_pos_phospho.explode("PTM_position")
final_ptm_pos_exploded = final_ptm_pos_exploded.sort_values(by=["Protein", "_struct_conf.beg_label_seq_id"])
final_ptm_pos_exploded.to_csv("dssp_summary_files_v3/summary_per_ptm_pos_protein_secondary_struct_features_dssp.csv",index=False)


###########################################
# Overall output 2                        #
# Confidence score and structure at STY   #
###########################################

# df with confidence
confidence_final = pd.concat(protein_confidence_df_list, ignore_index=True)
confidence_final = confidence_final[["Protein", "PTM_model", "PTM_position_all", "phospho", "mean_confidence_score", "_atom_site.label_comp_id", "_atom_site.auth_seq_id"]]
# df with struct (before explode)
final_ptm_pos = final_ptm_pos[["Protein", "PTM_model", "PTM_position", "PTM_position_all", "_struct_conf.conf_type_id", "_struct_conf.beg_label_seq_id",
 "_struct_conf.end_label_seq_id"]]
#list to str
#final_ptm_pos["PTM_position"] = final_ptm_pos["PTM_position"].str.join(", ")

#merge 
conf_merge = confidence_final.merge(final_ptm_pos, on =["Protein", "PTM_model", "PTM_position_all"])
# filter to keep only rows where the position of confidence score falls within range of struct
conf_merge = conf_merge[
(conf_merge['_atom_site.auth_seq_id'] >= conf_merge['_struct_conf.beg_label_seq_id']) &
(conf_merge['_atom_site.auth_seq_id'] <= conf_merge['_struct_conf.end_label_seq_id'])]
conf_merge.to_csv("dssp_summary_files_v3/conf_merge.csv",index=False)


##################################################################
# Overall output 3                                               #
# At the protein level, does the amount of                       #   
# structural features change between the PTM and no PTM models?  # 
##################################################################

protein_model_overall_df = pd.concat(overall_protein_structure, ignore_index=True)
protein_model_overall_df.to_csv("dssp_summary_files_v3/summary_per_protein_secondary_struct_features_dssp.csv",index=False)


##################################################
# Overall output 4                               #
# Extract confidence score per aa in protein     #
# combine into one file                          #   
##################################################
      
confidence_final = pd.concat(protein_confidence_df_list, ignore_index=True)
confidence_final.to_csv("dssp_summary_files_v3/confidence_score_AF.csv", index=False)