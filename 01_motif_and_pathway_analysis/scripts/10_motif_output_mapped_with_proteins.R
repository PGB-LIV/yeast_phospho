#################################################### 
# This script is to add columns to rmotifx output. #
# This is to add proteins the motif is present in. #
####################################################
library(dplyr)
library(readxl)
library(stringr)

# file with IDs from uniprot
uniprot_conversion <- read_xlsx("01_motif_and_pathway_analysis/inputs/uniprotkb_proteome_UP000002311_2026_02_19.xlsx")
# if no gene symbol, add entry to column
uniprot_conversion$`Gene Names (primary)`<-ifelse(is.na(uniprot_conversion$`Gene Names (primary)`), uniprot_conversion$Entry, uniprot_conversion$`Gene Names (primary)`)


for (flr_category in c("gold", "gsb")){
  # rmotifx foreground
  foreground <- read.csv(file = paste("01_motif_and_pathway_analysis/inputs/", flr_category, "_motif_seqs.txt", sep = ""), header = FALSE)
  # rmotifx results
  motif_df <- read.csv(file = paste("01_motif_and_pathway_analysis/outputs/", flr_category,"_motif_seqs/All_motifs_", flr_category, "_motif_seqs.csv", sep = ""))
  # column with motifs
  motifs <- motif_df$motif
  #print(flr_category)
  #print(foreground)
  counter <- 0
  overall <- list()
  # loop through motifs
  for (motif in motifs){
    counter <- counter + 1
    print(motif)
    # filter to get only proteins in foreground that contain the motif
    protein_with_motif <- filter(foreground, grepl(motif, V1))
    # protein uniprot id
    protein_with_motif_uniprot <- protein_with_motif$V2
    # keep unique
    protein_with_motif_uniprot <- unique(protein_with_motif_uniprot)
    # split on "|"
    protein_with_motif_uniprot_entry <- str_split(protein_with_motif_uniprot, "\\|")
    # keep accession
    accession <- sapply(protein_with_motif_uniprot_entry, `[`, 2)
    # match accession in lookup df
    match_pos <- match(accession, uniprot_conversion$Entry)
    # if the id was not in lookup table (NA), input accession, else put gene name primary from lookup
    gene_name <- ifelse(is.na(match_pos), accession, uniprot_conversion$`Gene Names (primary)`[match_pos])
    # join into list - uniprot id and primary gene names
    all_protein_with_motif_uniprot_gene_names <- paste(gene_name, collapse = ":")
    all_protein_with_motif_uniprot <- paste(protein_with_motif_uniprot, collapse = ":")
    print(any(grepl("[\r\n]", all_protein_with_motif_uniprot)))
    # collect proteins for each motif in df
    overall[[counter]] <- data.frame(motif = motif, proteins_with_motif_UniProt_ID = all_protein_with_motif_uniprot, proteins_with_motif_gene_names = all_protein_with_motif_uniprot_gene_names)
  }
  
  # join results for all motifs in FLR cat
  motif_proteins_df <- do.call(rbind, overall)
  # add the primary protein name back to rmotifx results
  flr_motif_proteins <- left_join(motif_df, motif_proteins_df, by = "motif")
  # rename col
  names(flr_motif_proteins)[names(flr_motif_proteins) == "proteins_with_motif_UniProt_ID"] <- "proteins_with_motif (UniProt ID)"
  names(flr_motif_proteins)[names(flr_motif_proteins) == "proteins_with_motif_gene_names"] <- "proteins_with_motif (gene names)"
  
  # output
  write.csv(flr_motif_proteins, paste("01_motif_and_pathway_analysis/outputs/", flr_category, "_motif_seqs/All_motifs_", flr_category, "_motif_seqs_with_proteins.csv", sep = ""), row.names = F)

  }
