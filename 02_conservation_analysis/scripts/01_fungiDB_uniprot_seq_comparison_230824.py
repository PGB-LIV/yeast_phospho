#################################################
# This script is for checking if the sequence   #
# in FungiDB is the same as in UniProt          #
#################################################
from Bio import SeqIO
import pandas as pd
import os
from itertools import chain
if not os.path.exists("outputs"):
    os.mkdir("outputs")
if not os.path.exists("outputs/processing_files"):
    os.mkdir("outputs/processing_files")

###################################
# 1 GSB Mappings                  # 
#                                 #
# Need to convert to FungiDB ID   #  
# using the uniprot mappings file #
###################################

# Multiple mappings GSB
yeast_build = pd.read_csv("/mnt/hc-storage/users/hleboswe/Yeast_GSB_Nov2024/All_site_formats_Updated/G2S1B_0.05_protein_pos_all_prot_mapping.csv")
yeast_build[["sp","name"]] = yeast_build["Protein"].str.split("|", n=1,expand=True)
yeast_build[["uniprot_transcript_id","entry_name"]] = yeast_build["name"].str.split("|", n=1,expand=True)
yeast_build[["uniprot_accession_id","transcript"]] = yeast_build["uniprot_transcript_id"].str.split("-", n=1,expand=True)

#counts
stya_val = yeast_build["PTM_residue"].value_counts()
print("STYA count", stya_val)
# remove ala
yeast_build_filtered = yeast_build[yeast_build["PTM_residue"].str.contains("S|T|Y")]
#counts
sty_val = yeast_build_filtered["PTM_residue"].value_counts()
print("STY count", sty_val)

# Mappings table - need to convert uniprot id -> other id
id_convert = pd.read_excel("inputs/uniprotkb_proteome_UP000002311_2024_07_19.xlsx") # Downloaded from UniProt (~6000 entries) as table
# dictionary
df_dict = id_convert.set_index("Entry")["Gene Names (ordered locus)"].to_dict() # better to use ordered locus - get complete mappings (alt. method was to loop through names and check to see if in FungiDB fasta but there are cases whereby there is a ordered locus name but it is not in FungiDB- check_ids.py)
#https://www.statology.org/pandas-convert-dictionary-to-dataframe/
df_convert = pd.DataFrame(list(df_dict.items()), columns=["uniprot_accession_id","ordered_locus_id"])
df_convert.to_csv("outputs/processing_files/uniprot_id_to_ordered_locus.csv")
# split on ; (multiple ids)
df_convert["ordered_locus_id"] = df_convert["ordered_locus_id"].str.split("; ")
# each mapping on own row
df_convert = df_convert.explode("ordered_locus_id")
df_convert.to_csv("outputs/processing_files/uniprot_id_to_ordered_locus_exploded.csv")
# keep first mapping
df_convert = df_convert.drop_duplicates(subset=["uniprot_accession_id"])
df_convert.to_csv("outputs/processing_files/uniprot_id_to_ordered_locus_exploded_unique.csv")


# Add the ordered locus to the multiple mappings
final = pd.merge(yeast_build_filtered, df_convert, how="left", on="uniprot_accession_id")
final = final[["Protein", "Protein_pos", "PTM_residue", "uniprot_transcript_id","uniprot_accession_id","ordered_locus_id"]]
final.to_csv("outputs/processing_files/GSB_converted.csv", index=False)

############################
# 2. FungiDB               #
#                          # 
############################

# K: seq, val= ID
fungi_db_seq_to_id = {}
fungi_db_ids = []
# K=ID, val=Seq
fungi_db_id_to_seq = {}
# FungiDB fasta file
for record in SeqIO.parse("inputs/FungiDB-68_ScerevisiaeS288C_AnnotatedProteins.fasta", "fasta"):
    # sequence
    seq = str(record.seq)
    protein_id = record.id.split("-t")[0]
    fungi_db_ids.append(protein_id) # store the FungiDB IDs
    # If the seq is already in dict, add the ID to the list
    if seq in fungi_db_seq_to_id:
        fungi_db_seq_to_id[seq].append(protein_id)
    # If seq not in dict, add ID
    else:
        fungi_db_seq_to_id[seq] = [protein_id]
        
    # If the protein_id is already in dict, add the seq to the list
    if protein_id in fungi_db_id_to_seq:
        fungi_db_id_to_seq[protein_id].append(seq)
    # If protein_id not in dict, add seq
    else:
        fungi_db_id_to_seq[protein_id] = [seq]

for seq, prot_id in fungi_db_seq_to_id.items():
    prot_id_str = ",".join(prot_id)
    fungi_db_seq_to_id[seq] = prot_id_str
    
# remove brackets        
for prot_id, seq in fungi_db_id_to_seq.items():
    seq_str = "".join(seq)
    fungi_db_id_to_seq[prot_id] = seq_str


###################################
# 2.2 UniProt                     # 
#                                 #
# Dict                            #
###################################
uniprot_db_id_to_seq = {}
uniprot_db_seq_to_id = {}

# UniProt fasta file (used in GSB mappings)
for record in SeqIO.parse("inputs/uniprotkb_proteome_UP000002311_2023_11_27.fasta", "fasta"):
    seq = str(record.seq)
    protin_record_id = record.id
    protein_id_accession = protin_record_id.split("|")[1] # this is if using entry
    # K=id, val=seq
    # if the id is already in dict, add seq to list
    if protein_id_accession in uniprot_db_id_to_seq:
        uniprot_db_id_to_seq[protein_id_accession].append(seq)
    else:
        uniprot_db_id_to_seq[protein_id_accession] = [seq]
   
   # check if seq in dict and add id as value
    if seq in uniprot_db_seq_to_id:           
        uniprot_db_seq_to_id[seq].append(protein_id)
    else:
        #print(seq)
        uniprot_db_seq_to_id[seq] = [protein_id]

for prot_id, seq in uniprot_db_id_to_seq.items():
    seq_str = "".join(seq)
    uniprot_db_id_to_seq[prot_id] = seq_str

for seq, prot_id in uniprot_db_seq_to_id.items():
    prot_id_str = ",".join(prot_id)
    uniprot_db_seq_to_id[seq] = prot_id_str

######################################
#                                    #
# 3. Check - FungiDB vs UniProt      #
#                                    # 
######################################

# uniprot id col - look up seq in the uniprot K= id, V=seq (id includes transcript)
final["uniprot_seq"] = final["uniprot_transcript_id"].map(uniprot_db_id_to_seq)
# ordered locus - look up seq in fungi_db2 dict K= id
final["fungi_seq"] = final["ordered_locus_id"].map(fungi_db_id_to_seq)
# use uniprot seq - loop up id in fungi_db K= seq
final["fungi_id_from_uni_seq"] = final["uniprot_seq"].map(fungi_db_seq_to_id)
# use uniprot seq - look up id in uniprot_seq_key K=seq
final["uniprot_seq_check"] = final["uniprot_seq"].map(uniprot_db_seq_to_id)

# Is the ordered locus ID in the FungiDB FASTA?
final["id_in_db"] = final["ordered_locus_id"].isin(fungi_db_ids)
# Check if the ordered locusID is the same as FungiDB ID from the seq
final["fungi_id_comparison"] = final["ordered_locus_id"] == final["fungi_id_from_uni_seq"]
final["uniprot_id_comparison"] = final["uniprot_transcript_id"] == final["uniprot_seq_check"]
final["comparison"] = final["uniprot_seq"] == final["fungi_seq"]


final.to_csv("outputs/processing_files/final.csv", index=False)
final_non_dup = final.drop_duplicates("Protein")
final_non_dup.to_csv("outputs/processing_files/final_non_dup.csv", index=False)

######################################
#                                    #
# 4. Filter                          #
#                                    # 
######################################

# Remove IDs whereby uniprot seq and fungi db seq not the same 
final = final[final["comparison"] == True]
final = final[["Protein","Protein_pos","PTM_residue","uniprot_accession_id","ordered_locus_id"]]
# to csv
final.to_csv("outputs/processing_files/final_filtered.csv", index=False)
# to tsv
final.to_csv("outputs/GSB_STY_conservation_yeast.tsv", sep="\t",index=None)
