############################################################ 
# This script is for exporting mapped yeast phosphosites   #
# in ortholog proteins from the conservation analysis      #
############################################################
geneprefix_to_species = {}
# dict of ortho protein -> species, seq
f = open("inputs/GenesByOrthologs_Summary_all_yeast_syntenic_orthologs.txt","r")
for line in f:
    line = line.rstrip()
    cells = line.split("\t")
    geneprefix_to_species[cells[0]] = [cells[2],cells[8]]
f.close()

# file containing conservation results
f = open("outputs/conservation_stats_yeast/summary_stats.tsv","r")

f_out = open("outputs/Supplementary_File_8_mapped_yeast_sites.tsv","w")
f_out.write("Protein\tPosition\tGenome\tPredicted_sequence\tMapped_yeast_site_FungiDB_ID\tMapped_yeast_site_UniProt_ID\n")

for line in f:
    line = line.rstrip()
    cells = line.split("\t")
    # yeast protein, PTM pos, residue
    yeast_prot_pos_res = cells[0] + "_" + cells[2] + "_" + cells[3]
    yeast_prot_pos_res_uniprot = cells[1] + "_" + cells[2] + "_" + cells[3]

    # mapped sites  
    if cells[12] != "":
        mapped_prots = cells[12].split(";")
        for mapped_prot in mapped_prots:
            #print(mapped_prot)
            mapped_prot_splitup = mapped_prot.rsplit("_",1)
            mapped_prot_id = mapped_prot_splitup[0]
            mapped_prot_pos = mapped_prot_splitup[1]
            species_prefix = mapped_prot.split("_")[0] #Species prefix is before first underscore
            
            genome = ""
            if mapped_prot_id in geneprefix_to_species:
                genome = geneprefix_to_species[mapped_prot_id][0]
                seq = geneprefix_to_species[mapped_prot_id][1]
            f_out.write(mapped_prot_id + "\t" + mapped_prot_pos + "\t" + genome + "\t" + seq + "\t" + yeast_prot_pos_res + "\t" + yeast_prot_pos_res_uniprot + "\n" )
