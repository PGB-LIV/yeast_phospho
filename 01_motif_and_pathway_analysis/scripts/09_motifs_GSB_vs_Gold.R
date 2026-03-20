#################################################
# This script is for finding the                #
# overlap in motifs from the rmotifx analysis   #
# with the Gold set versus the GSB set          #
#################################################
library(VennDiagram)
# load files
motif_gold <- read.csv(file="01_motif_and_pathway_analysis/outputs/gold_motif_seqs/All_motifs_gold_motif_seqs.csv")
motif_gsb <- read.csv(file="01_motif_and_pathway_analysis/outputs/gsb_motif_seqs/All_motifs_gsb_motif_seqs.csv")

# comparison of motifs
only_in_gold <- motif_gold$motif[!(motif_gold$motif %in% motif_gsb$motif)]
only_in_gsb <- motif_gsb$motif[!(motif_gsb$motif %in% motif_gold$motif)]
overlap <- motif_gsb$motif[motif_gsb$motif %in% motif_gold$motif]

#plot
venn.diagram(
  x = list(motif_gold$motif, motif_gsb$motif),
  category.names = c("Gold motifs" , "GSB motifs"),
  filename = "01_motif_and_pathway_analysis/outputs/motifs_comparison/motif_GSB_versus_Gold.png", width = 4800, height = 2800, resolution = 330, 
  main = "Comparison of significantly enriched motifs", main.cex = 2, lwd = 3, cex = 2, cat.cex = 2, fill = "#6699CC"
)



#####################
# top 30 fg matches #
#####################
motif_gsb <- read.csv(file="01_motif_and_pathway_analysis/outputs/gsb_motif_seqs/All_motifs_gsb_motif_seqs.csv")
motif_gsb <- motif_gsb [order(motif_gsb$fg.matches, decreasing = T), ]
motif_gsb <- motif_gsb[1:30,]

motif_gold <- read.csv(file="01_motif_and_pathway_analysis/outputs/gold_motif_seqs/All_motifs_gold_motif_seqs.csv")
motif_gold <- motif_gold[order(motif_gold$fg.matches, decreasing = T), ]
motif_gold <- motif_gold[1:30,]

venn.diagram(
  x = list(motif_gold$motif, motif_gsb$motif),
  category.names = c("Gold motifs" , "GSB motifs"),
  filename = "01_motif_and_pathway_analysis/outputs/motifs_comparison/motif_GSB_versus_Gold_fg_matches.png", width = 4800, height = 2800, resolution = 330, 
  main = "Comparison of significantly enriched motifs", main.cex = 2, lwd = 3, cex = 2, cat.cex = 2, fill = "#6699CC"
)

