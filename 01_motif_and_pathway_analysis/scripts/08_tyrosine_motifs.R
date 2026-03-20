##############################################
# this script is for the motif analysis      #
# central residue: Y                         #  
#                                            #  
#             ---- Gold ----                 #
#                   &                        #   
#             ---- GSB ----                  #     
##############################################
require(rmotifx)

######################
# 1. Motif analysis  #
#                    #
# rmotifx            #
######################


# rmotifx background - 15mers from peptides in FDR_output
background <- read.csv(file="01_motif_and_pathway_analysis/inputs/motif_Y_background.txt", header = FALSE)
bg.seqs <- as.character(background$V1)
bg_list <- unique(as.character(background$V3))
write.table(bg.seqs,  file = paste0("01_motif_and_pathway_analysis/outputs/background_seq_tyrosine.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# Gold and GSB
file_list  <- c("01_motif_and_pathway_analysis/inputs/gold_motif_seqs.txt", "01_motif_and_pathway_analysis/inputs/gsb_motif_seqs.txt")
for (f in file_list){
  foreground <- read.csv(file = f, header = FALSE)
  # central residue in 15mer
  foreground$central <- substr(foreground$V1, 8, 8)
  # filter for y 15mers
  foreground <- subset(foreground, central == "Y")
  fg.seqs <- as.character(foreground$V1)
  
  
  
  f <- strsplit(f,"/")[[1]][3]
  f <- substring(f, 1, nchar(f) - 4)
  print(f)
  #find enriched motifs - y
  mot <- motifx(fg.seqs, bg.seqs, central.res = 'Y', min.seqs = 20, pval.cutoff =  1e-6)
  print(mot)
  write.csv(mot, file = paste0("01_motif_and_pathway_analysis/outputs/tyrosine/", f, "/All_motifs_", f, ".csv"), row.names = FALSE) 
  
}
