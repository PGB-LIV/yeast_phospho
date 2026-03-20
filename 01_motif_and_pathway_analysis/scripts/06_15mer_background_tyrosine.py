######################################################
# This script is for generating 15mer background     # 
# for tyrosine                                       #   
######################################################
import sys
import pandas as pd
from Bio import SeqIO
import os
import re

# Database used in analysis
database = "/home/hleboswe/hc-storage/PXD028028/2023-11-27-decoys-contam-uniprotkb_proteome_UP000002311_2023_11_27.fasta"
# into dictionary
seq_dict = SeqIO.to_dict(SeqIO.parse(database, "fasta"))

# input and outputs for each exp
fdr_input = sys.argv[1]
background_output= sys.argv[2]

overall_background = []

# For each exp...
# FDR_output.csv from mzidFLR
all_sites = fdr_input + "/FDR_0.01/FDR_output.csv"
df = pd.read_csv(all_sites)
# <1% FDR
df = df.loc[df['FDR'] <= float(0.01)]
# Drop duplicates on protein - peptide
df = df.drop_duplicates(subset=["Protein","Peptide"], keep="first").reset_index(drop=True)
protein_list = df['Protein'].to_list()
peptide_list = df['Peptide'].to_list()
# Loop through protein/peptide
for peptide, protein in zip(peptide_list, protein_list):
    # if peptide maps to multiple proteins, keep first only
    protein = protein.split(":")[0]
    # Look up full protein seq in dict
    record = seq_dict[protein]
    seq_temp = str(record.seq)
    # index of Y in peptide
    pep_index = [y.start() for y in re.finditer("Y", peptide)]
    # If there is at least 1 Y
    if len(pep_index) > 0:
        # find index of peptide in protein
        prot_index = [peptide_pos.start() for peptide_pos in re.finditer(peptide, seq_temp)]
        # for each instance of the peptide in the protein (I guess this wont occur very often)
        for index in prot_index:
            # create new list with index of Y in the protein (index of Y in peptide + index of peptide in protein)
            aa_index_in_prot = [y_index + index for y_index in pep_index]
            # For each index of Y in the protein (from identified peptides)
            for aa in aa_index_in_prot:
                #if len(prot_index)>1:
                    #print(aa,index,prot_index, peptide, seq_temp) 
                # clean protein seq
                seq_temp = str(record.seq)
                if seq_temp[aa] != "Y":
                    print("15mer fail from", peptide, "/", seq_temp[aa])
                # needs to be at least 7 aa after Y, so add "_" to end of seq
                if aa + 8 > len(seq_temp):
                    while aa + 8 > len(seq_temp):
                        seq_temp += "_"
                # needs to be at least 7 aa before Y, so add "_" before seq
                if aa < 7:
                    seq_temp = ("_" * (7 - aa)) + seq_temp
                    aa += (7 - aa) # update index of Y in protein since the dashes will shift the position
                    if seq_temp[aa] != "Y":
                        print("15mer fail 2", seq_temp[aa - 7:aa + 8])
                # string of 15mer and protein
                protein_str = record.id + ":" + seq_temp[aa - 7:aa  + 8]    
                if protein_str not in overall_background:
                    overall_background.append(protein_str)    
    
print(len(overall_background))
df2=pd.DataFrame(list(zip(overall_background)),columns=["Proteins"])
# Split protein and 15mer into columns   
df2[["Proteins","Sequence"]]=df2["Proteins"].str.split(":", n=1,expand=True)                               
updated_protein_list=[]

#df2=df2.reset_index(drop=True)
# UniProt accession
for a in range(len(df2)):
    protein = df2.loc[a,'Proteins']
    try:
        uni_id = protein.split("|")[1]
    except:
        print("error")
        uni_id = "-"
    updated_protein_list.append(uni_id)

# ID = UniProt accession      
df2['ID'] = updated_protein_list
df2 = df2[["Sequence","Proteins","ID"]]
df2.to_csv(background_output, index=False)
     


