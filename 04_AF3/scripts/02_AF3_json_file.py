####################################################
# This py script is for generating json files      #
# for AF3 (on the HPC).                            #
####################################################

# Specifically, for every yeast protein with >=10 gold phosphosites, generate two AF3 json input files (one with phosphosites and one without)

import pandas as pd
import Bio.SeqIO as SeqIO
import json
import os
# Create output directory
os.makedirs("af_input/", exist_ok=True)

# Input file requires that CCD codes are given for mods
ccd_ptm_code = {
"S" : "SEP",
"T" : "TPO",
"Y" : "PTR"
}

# File with yeast proteins that have >=10 gold phosphosites (site_count_per_protein.R)
gold_proteins = pd.read_csv("multiple_gold_sites.csv")
# List of proteins
protein_list = gold_proteins["Protein"].tolist()

# GSB multiple mapping file (needed to obtain PTM pos and residue)
gsb = pd.read_csv("G2S1B_0.05_protein_pos_all_prot_mapping.csv")
# Filter to contain gold sites only
gsb_gold = gsb[gsb['PTM_FLR_category'] == 'Gold']
# Filter to contain STY sites only (ie remove the A decoys)
gsb_gold = gsb_gold[gsb_gold['PTM_residue'].isin(["S", "T", "Y"])]
gsb_dict = {}

# GSB DF -> dictionary with protein ID as the key and the PTM pos and residue as a list of values
for protein, ptm_pos, ptm_res in zip(gsb_gold.Protein, gsb_gold.Protein_pos, gsb_gold.PTM_residue):
    # If the protein is already in the dict, add a list of the next PTM pos and res
    if protein in gsb_dict:
        gsb_dict[protein].append([ptm_pos, ptm_res])
    else:
    # If the protein is not already in the dict, add it with the PTM pos and residue as a list of values
        gsb_dict[protein] = [[ptm_pos, ptm_res]]

# check
#print(gsb_dict)


# Database fasta file - we need to obtain the sequence of the proteins
yeast_proteins = "2023-11-27-decoys-contam-uniprotkb_proteome_UP000002311_2023_11_27.fasta"
# parse into dict with protein ID as key
records = SeqIO.to_dict(SeqIO.parse(yeast_proteins, "fasta"))
#print(protein_list)

file_versions = ["non_phospho","phospho"]
# Now, loop through each protein and create two json files
for protein in protein_list:
    #af_protein = {}
    # Use protein ID as the key to get the protein record
    protein_record = records[protein]
    protein_seq = str(protein_record.seq)
    protein_name = protein.split("|")[1]
    for version in file_versions:
        # ID and seq
        seq_dict = {
        "id" : "A",
        "sequence" : protein_seq,
        "templates" : []}
        # For phospho file, we need to add mods
        if version == "phospho":
            # Look up PTM pos and residue for protein
            ptm_pos = gsb_dict[protein]
            mods_ls = []
            # Each PTM is specified in its own dict
            for ptm in ptm_pos:
                # Convert from residue to CCD code
                ccd = ccd_ptm_code[ptm[1]]
                mod_dict = {
                "ptmType" : ccd,
                "ptmPosition" : ptm[0]
                }
                # Add PTM info (dict) to list
                mods_ls.append(mod_dict)
                # Dict that now contains ID, seq and mods
            seq_dict["modifications"] = mods_ls
        # put ID, seq (and mods) into dict
        prot_dict = {
        "protein" : seq_dict}
        # Job name - specify PTM version
        id_name = protein_name + "_" + version
        
        # Add to final dict
        af_protein = {
        "name" : id_name,
        "sequences" : [prot_dict], 
        "modelSeeds" : [1],
        "dialect" : "alphafold3",
        "version" : 1
        }
       
        # json output file
        with open("af_input/" + protein_name + "_" + version + ".json", "w") as f:
            json.dump(af_protein, f, indent=4)
          