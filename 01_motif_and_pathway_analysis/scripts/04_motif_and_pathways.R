##############################################
# this script is for the motif analysis      #
# central residue= ST (Y not major in yeast) #  
#                                            #  
#             ---- Gold ----                 #
#                   &                        #   
#             ---- GSB ----                  #     
##############################################
require(rmotifx)
require(clusterProfiler)
library(stringr)
require(ggplot2)
require(ggseqlogo)
require(gridExtra)
library(dplyr)
library(org.Sc.sgd.db)
library(tidyverse)
library(pheatmap)
library(svglite)
######################
# 1. Motif analysis  #
#                    #
# rmotifx            #
######################


# rmotifx background - 15mers from peptides in FDR_output
background <- read.csv(file ="01_motif_and_pathway_analysis/inputs/motif_ST_background.txt", header = FALSE)
bg.seqs <- as.character(background$V1)
bg_list <- unique(as.character(background$V3))
write.table(bg.seqs, file = paste0("01_motif_and_pathway_analysis/outputs/background_seq.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# Gold (motif and enrichment analysis) and GSB (motif only)
file_list<-c("01_motif_and_pathway_analysis/inputs/gold_motif_seqs.txt", "01_motif_and_pathway_analysis/inputs/gsb_motif_seqs.txt")
for (f in file_list){
  foreground <- read.csv(file = f, header = FALSE)
  # central residue in 15mer
  foreground$central <- substr(foreground$V1, 8, 8)
  # filter for ST 15mers
  foreground <- subset(foreground, central == "S" | central == "T")
  fg.seqs <- as.character(foreground$V1)
  
  f <- strsplit(f, "/")[[1]][3]
  f <- substring(f, 1, nchar(f) - 4)
  print(f)
  
  #find enriched motifs
  mot <- motifx(fg.seqs, bg.seqs, central.res = 'ST', min.seqs = 20, pval.cutoff = 1e-6)
  write.csv(mot, file = paste0("01_motif_and_pathway_analysis/outputs/", f, "/All_motifs_", f, ".csv"), row.names = FALSE) 
  
  
  ################################################################################
  # 2. Group motifs:                                                             #      
  # https://www.pnas.org/doi/full/10.1073/pnas.0609836104                        #
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2846625/                        # 
  # https://www.sciencedirect.com/science/article/pii/S0022030222005483#bib49    #
  # May need to change/adjust depending on motifs (particularly the other group) #
  ################################################################################ 
  # Group 1: Proline directed
  mot$pro <- str_sub(mot$motif, 12, 12) == "P"
  
  # Acidic: ≥5 Glu/Asp at +1 to +6 (acidic)
  mot[c("minus","pos")] <- str_split_fixed(mot$motif, "ST", 2)
  mot$pos_loc <- substr(mot$pos, 2, 7)
  mot$de <- str_count(mot$pos_loc,"D") + str_count(mot$pos_loc, "E")
  mot$acidic_1 <- ifelse(mot$de >= 5, TRUE, FALSE)
  
  # Group 2: Basophilic  Arg/Lys at −3 (basic)
  mot$basic <- (str_sub(mot$motif, 5, 5) == "R") |
    (str_sub(mot$motif, 5, 5) == "K")
  
  # Group 3: acidophilic Glu/Asp at +1/+2 or +3 (acidic),
  mot$acidic_2 <- 
    (str_sub(mot$motif, 12, 12) %in% c("E", "D"))| 
    (str_sub(mot$motif, 13, 13) %in% c("E", "D"))| 
    (str_sub(mot$motif, 14, 14) %in% c("E", "D"))
  
  # Group 2: Basophilic ≥2 Arg/Lys at −6 to −1 (basic)
  mot$minus_loc <- substr(mot$minus, 2, 7)
  mot$rk <- str_count(mot$minus_loc,"R") + str_count(mot$minus_loc, "K")
  mot$basic_2 <- ifelse(mot$rk > 1,TRUE, FALSE)
  
  # Group 3: Other
  
  # Divide up the Other group
  # Group 3.1: Any E or D residues - acidic (but less strict than the acidophilic definition)
  mot$DE_residues <- str_detect(mot$motif, "E") | str_detect(mot$motif, "D")
  # Group 3.2: Any N, Q, S or T (more than 1 for ST as this is the [ST]) - hydrophilic neutral aa
  mot$NQST_residues <- str_detect(mot$motif, "N") | str_detect(mot$motif, "Q") | str_count(mot$motif, fixed("S")) > 1 | (str_count(mot$motif,fixed("T")) > 1)
  # Group 3.4: Any KR residues  (basophilic)
  mot$KR_residues <- str_detect(mot$motif, "K") | str_detect(mot$motif, "R")
  
  
  
  # Decision tree
  # is it proline directed?, is it acidophilic?, is it basophilic?, is it acidophilic?, is it basophilic?, does it contain aa with polar uncharged side chains?, does it contain DE residues (acidic)? does it contain basic aa?
  mot$group <- ifelse(mot$pro,"Proline-directed",
                    ifelse(mot$acidic_1,"Acidophilic",
                           ifelse(mot$basic,"Basophilic",
                                  ifelse(mot$acidic_2,"Acidophilic",
                                         ifelse(mot$basic_2,"Basophilic",
                                                ifelse(mot$NQST_residues,"NQST_residues",
                                                       ifelse(mot$DE_residues,"DE_residues",
                                                              ifelse(mot$KR_residues,"KR_residues",
                                                                     "Other"))))))))
  
  
  groups <- c("Acidophilic", "Basophilic", "Proline-directed", "NQST_residues", "DE_residues", "KR_residues", "Other") 

  #ClusterProfiler background - all gold (ST) phosphoproteins
  cluster_background <- read.csv(file= "01_motif_and_pathway_analysis/inputs/All_phosphoproteins_background.txt", header= FALSE)
  cluster_bg_list <- unique(as.character(cluster_background$V1))
  
  # seqlogo
  seq_plots <- list()
  # cp
  cp_plots <- list()
  #create GO dataframe
  GO_full <- data.frame()
  GO_simplified <- data.frame()
  counter <- 1
  #groups <-c("Acidophilic","Basophilic", "Proline-directed", "Other")
  # Gold foreground - clusterprofiler
  if (f == "gold_motif_seqs"){
    for (motif_group in groups){
      
      ################################################
      #1. motif seqlogo of enriched motif sequences  #
      ################################################
      pos <- c("-7", "-6", "-5", "-4", "-3", "-2", "-1", "p", "1", "2", "3", "4", "5", "6", "7")
      break_list <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
      print(motif_group)
      
      # Get the motifs that fall into a group
      mot_filtered <- subset(mot, group == motif_group)
      motifs_in_group <- mot_filtered$motif
      # Filter foreground by motifs in the motif group
      foreground_motif <- filter(foreground, grepl(paste(motifs_in_group, collapse='|'),V1))
      # seq
      foreground_motif_seq <- foreground_motif$V1
      
      # csv
      write.table(foreground_motif_seq, file = paste0("01_motif_and_pathway_analysis/outputs/",f,"/",motif_group, "_15mer_", f, ".txt"), quote = FALSE, row.names = FALSE,col.names = FALSE)
      # seq logo
      motif_plot <- ggplot() + geom_logo(foreground_motif_seq, method = "probability") + scale_x_continuous(labels = pos, breaks = break_list) + 
        ggtitle(motif_group) + theme(legend.position = "none")
      seq_plots[[counter]] <- motif_plot
      
      ####################################################
      # 2. Pathway enrichment analysis - ClusterProfiler #
      ####################################################
      #https://hbctraining.github.io/DGE_workshop_salmon/lessons/AnnotationDbi_lesson.html
      #https://rpubs.com/mano2991/GO
      #https://rpubs.com/nsato/gsea_dsigdb
      #https://mirrors.dotsrc.org/bioconductor/help/course-materials/2019/CSAMA/L1.5-bioc-annotation.html
      #https://mirrors.dotsrc.org/bioconductor/help/course-materials/2019/CSAMA/materials/lectures/lecture-08a-annotation.html#34
      
      
      # background = all gold ST phosphosites
      # Define foreground and background gene lists.
      # The foreground list should be contained within the background list.
      foreground_motif$motif <- motif_group
      colnames(foreground_motif) <- c("15mer", "Protein", "PTM_position", "UniProt_accession", "Central_residue", "Motif_group")
      write.csv(foreground_motif, paste0("01_motif_and_pathway_analysis/outputs/", f, "/", motif_group, "_", f, ".csv"), row.names = F)
      
      # UniProt accession - unique
      fg_list <- unique(as.character(foreground_motif$UniProt_accession))
      
      #  GO over-representation test
      enrich_result <- enrichGO(gene = fg_list,
                                OrgDb = org.Sc.sgd.db,
                                universe = cluster_bg_list,
                                keyType = "UNIPROT",
                                ont = "ALL",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 1, 
                                qvalueCutoff  = 1)
      
      
      # To a dataframe
      GO_term_dataframe <- as.data.frame(enrich_result@result)
      #add motif
      GO_term_dataframe$motif_group <- motif_group
      # Write the data to a CSV file
      output <- paste("01_motif_and_pathway_analysis/outputs/", f, "/","GO_full_", motif_group, "_", f, ".csv", sep="")
      write.csv(GO_term_dataframe, output, row.names = F)
      
      # add to an overall dataframe
      GO_full <- rbind(GO_full, GO_term_dataframe)

      # filter for significance (p.adjust < 0.01)
      results_filtered <- filter(enrich_result, p.adjust < 0.01, qvalue < 0.2)
      
      # Simplify using 0.6 similarity cutoff
      result_simplified <- clusterProfiler::simplify(results_filtered, 
                                                     by = "p.adjust",
                                                     cutoff = 0.6,
                                                     select_fun = min)
      
      
      # filter on pval andq val
      #https://yulab-smu.top/biomedical-knowledge-mining-book/clusterProfiler-dplyr.html
      #https://support.bioconductor.org/p/9139934/
      # https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html - most signif
      
      # dot plot by ontology
      dotplot_fig <- clusterProfiler::dotplot(result_simplified, title = motif_group, split = "ONTOLOGY",showCategory = 5) + facet_grid(ONTOLOGY ~., scale = "free")
      ggsave(file = paste0("01_motif_and_pathway_analysis/outputs/", f, "/cp_", motif_group, "_", f, ".svg"), dpi = 330, width = 6.5, height = 6.2)
      ggsave(file = paste0("01_motif_and_pathway_analysis/outputs/", f, "/cp_", motif_group, "_", f, ".png"), dpi = 330, width = 6.5, height = 6.2)
      
      cp_plots[[counter]] <- dotplot_fig
      
      # to dataframe
      result_simplified <- as.data.frame(result_simplified@result)
      # add motif group
      result_simplified$motif <- motif_group
      # add to overall dataframe for enrichment analysis after simplify
      GO_simplified <- rbind(GO_simplified, result_simplified)
      cp_output <- paste("01_motif_and_pathway_analysis/outputs/", f, "/GO_simplified_", motif_group, "_", f, ".csv", sep = "")
      write.csv(result_simplified, cp_output, row.names=F)
      
      counter <- counter+1
    }
    # combined results across motifs
    write.csv(GO_full, file = paste0("01_motif_and_pathway_analysis/outputs/",f,"/GO_full_",f,".csv"),row.names=F)
    write.csv(GO_simplified, file = paste0("01_motif_and_pathway_analysis/outputs/",f,"/GO_simplified_",f,".csv"),row.names=F)
    
    overall_seq_plot <- grid.arrange(grobs = seq_plots, ncol = 2)
    overall_cp_plot <- grid.arrange(grobs = cp_plots, ncol = 2)
    
    ggsave(file = paste0("01_motif_and_pathway_analysis/outputs/", f, "/rmotifx_all_", f, ".png"), overall_seq_plot, height = 10, width = 7, dpi = 330)
    ggsave(file = paste0("01_motif_and_pathway_analysis/outputs/", f, "/cp_all_", f, ".png"), overall_cp_plot, dpi = 330, height = 30, width = 16)
    
    #########################
    # 3. Ontology heatmap   #
    #                       #
    #########################
    # heatmap of each ont class
    GO_ont <- c("BP", "MF", "CC")
    for (ontology_group in GO_ont){
      GO_simplified_ont <- subset(GO_simplified, ONTOLOGY == ontology_group)
      # order by p.adjust 
      GO_simplified_ont <- GO_simplified_ont[order(GO_simplified_ont$p.adjust, decreasing = FALSE),]
      # unique Go ids from simplified df - top 25 most signif terms
      unique_go_id_values <- unique(GO_simplified_ont[['ID']])[1:25]
      # filter full df by the unique values in the simplified df
      data_filtered <- GO_full[GO_full$ID %in% unique_go_id_values, ]
      # filter for 0.01 p val
      data_filtered <- data_filtered %>% group_by(ID) %>% filter(any(p.adjust <= 0.01)) %>% ungroup()
      # negative log
      data_filtered_pval <- data_filtered %>% mutate(neg_log_p.adjust = -log10(p.adjust)) %>% dplyr::select(Description, motif_group, neg_log_p.adjust) %>% spread(motif_group, neg_log_p.adjust) %>% column_to_rownames('Description')
      # Counts
      heatmap_counts <- data_filtered %>% dplyr::select(Description, motif_group, Count) %>% spread(motif_group, Count) %>% column_to_rownames('Description')
      # heatmap_counts
      heatmap_counts <- as.matrix(heatmap_counts)
      # asterisks = significant values
      heatmap_counts <- ifelse(data_filtered_pval >= -log10(0.01), paste0("n= ", heatmap_counts, "*"), paste0("n= ", heatmap_counts))
      
      
      # rows and col number
      row_count <- nrow(data_filtered_pval)
      col_count <- ncol(data_filtered_pval)
      # start with black
      overall_colour <- matrix("black", nrow = row_count, ncol = col_count)
      
      # heatmap plot
      output <- paste0("01_motif_and_pathway_analysis/outputs/", f, "/heatmap_", ontology_group, "_", f, ".png")
      ont_heatmap <- pheatmap(data_filtered_pval, 
                              scale = 'none', 
                              clustering_distance_rows = 'euclidean', 
                              clustering_method = 'complete', 
                              show_rownames = TRUE, 
                              display_numbers = heatmap_counts,  
                              cellnote.fontcolor = 'blue',  
                              fontsize = 16, 
                              fontsize_number = 12,
                              angle_col = 0, 
                              color = colorRampPalette(c("#6699CC", "white", "#CC6677"),space ="rgb")(100),
                              filename = output,  
                              width = 18, 
                              height = 14, 
                              number_color =  overall_colour,
                              units = "in")
    }
  } else{
    # for non Gold foreground, motif plot only
    for (motif_group in groups){
      
      ################################################
      #1. motif seqlogo of enriched motif sequences  #
      ################################################
      pos <- c("-7", "-6", "-5", "-4", "-3", "-2", "-1", "p", "1", "2", "3", "4", "5", "6", "7")
      break_list <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
      print(motif_group)
      
      # Get the motifs that fall into a group
      mot_filtered <- subset(mot, group == motif_group)
      motifs_in_group <- mot_filtered$motif
      # Filter foreground by motifs in the motif group
      foreground_motif <- filter(foreground, grepl(paste(motifs_in_group, collapse='|'), V1))
      # seq
      foreground_motif_seq <- foreground_motif$V1
      
      # csv
      write.table(foreground_motif_seq, file = paste0("01_motif_and_pathway_analysis/outputs/", f, "/", motif_group, "_15mer_", f, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
      # seq logo
      motif_plot <- ggplot() + geom_logo(foreground_motif_seq, method ="probability") + scale_x_continuous(labels = pos, breaks = break_list) + ggtitle(motif_group) + theme(legend.position = "none")
      seq_plots[[counter]] <- motif_plot
      counter <- counter + 1

  
    motif_plot_overall <- grid.arrange(grobs = seq_plots, ncol = 2)
    ggsave(file = paste0("01_motif_and_pathway_analysis/outputs/", f, "/rmotifx_all_", f, ".png"), motif_plot_overall, height = 10, width = 7, dpi = 330)
    }
  }
}
