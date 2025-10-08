#######################################################################
# This script is for generating pre-alignment fasta files of target   # 
# and ortholog                                                        #   
#Plan:                                                                #
# Read tsv file downloaded from FungiDB                               #
# Grab all orthologs that are the same as the target input protein    #
# Create a fasta file                                                 #
# Call muscle                                                         #
# Write protein-level MSA                                             #   
#######################################################################
import os
import Bio.SeqIO as SeqIO

INPUT_FILE = "inputs/GenesByOrthologs_Summary_all_yeast_syntenic_orthologs.txt"
DATA_FOLDER = "outputs/fastas_pre_aligment_syntenic_Sept24/"

yeast_proteins = "inputs/FungiDB-68_ScerevisiaeS288C_AnnotatedProteins.fasta"
records = SeqIO.parse(yeast_proteins, "fasta")

# DB of FungiDB Fasta, K: ID, V: seq
yeast_id_to_seq = {}
for record in records:
    #print(record.id)
    gene_id = record.id.split("-t")[0]
    yeast_id_to_seq[gene_id] = record.seq
    print("inserting",gene_id)

# Ortholog file
f = open(INPUT_FILE,"r")
# K: yeast protein, V: orthologs
input_protein_to_orthologs = {}
counter = 0
for line in f:
    if counter != 0:
        line = line.rstrip()
        cells = line.split("\t")
        input_proteins = cells[5]
        # Some have two entries
        proteins = input_proteins.split(",")

        for input_protein in proteins:
            #print(input_protein)
            mapped_ortholog_protein = cells[0]
            mapped_ortho_seq = cells[8]

            ortho_proteins = []
            if input_protein in input_protein_to_orthologs:#if already in dict, get back the values
                ortho_proteins = input_protein_to_orthologs[input_protein]
                #print(ortho_proteins)
            # check to see if the ortholog is already in dict for input protein (ie we only want to kee the first one- the first transcript - ortholog need to be unique)
            ortho_present = False
            for ortho in ortho_proteins:
                #ortholog id
                ortho_split = ortho.split(":")[0]
                #if the mapped ortholog is the same as an ortholog in the dict,
                if mapped_ortholog_protein == ortho_split:
                    #then set ortho_present to True
                    ortho_present = True
                    print(mapped_ortholog_protein, "is already in dictionary for", input_protein)
            # If the mapped ortholog did not match any orthologs in dict (so ortho_present is still False)
            if ortho_present is False:
                # add ortholog to dict
                ortho_proteins.append(mapped_ortholog_protein + ":" + mapped_ortho_seq)
                input_protein_to_orthologs[input_protein] = ortho_proteins # key is input protein, val is the gene id
    counter += 1

#print(input_protein_to_orthologs.keys())
#print(yeast_id_to_seq.keys())
# output dir
if os.path.exists(DATA_FOLDER) == False:
    os.mkdir(DATA_FOLDER)
for input_protein in input_protein_to_orthologs:
    # for each yeast protein, get the orthologs
    orthologs = input_protein_to_orthologs[input_protein]
    # if the input protein is in the FungiDB fasta dict, get the seq
    if input_protein in yeast_id_to_seq:
        yeast_seq = yeast_id_to_seq[input_protein]
        # new file per input protein
        f_out = open(DATA_FOLDER + input_protein + "_orthos_notaligned.fasta","w")
        # yeast id
        f_out.write(">" + input_protein + "\n")
        # yeast seq
        f_out.write(str(yeast_seq) + "\n")
        # add each ortho id and seq for the yeast protein
        for ortholog in orthologs: #orthologs for input protein
            cells = ortholog.split(":")
            ortho_id = cells[0]
            ortho_seq = cells[1]
            f_out.write(">" + ortho_id + "\n")
            f_out.write(ortho_seq + "\n")
        f_out.close()
    else:
        print("protein not found in yeast input set",input_protein)
