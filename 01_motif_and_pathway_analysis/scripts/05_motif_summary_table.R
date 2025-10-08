############################
# Motif summary table      #
#                          #
############################


# rmotifx output
mot <- read.csv(file = paste0("01_motif_and_pathway_analysis/outputs/gold_motif_seqs/All_motifs_gold_motif_seqs.csv")) 
################################################################################
# 2. Group motifs:                                                             #      
# https://www.pnas.org/doi/full/10.1073/pnas.0609836104                        #
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2846625/                        # 
# https://www.sciencedirect.com/science/article/pii/S0022030222005483#bib49    #
# May need to change/adjust depending on motifs (particularly the other group) #
################################################################################ 
# Group 1: Proline directed
mot$pro<- str_sub(mot$motif, 12, 12) == "P"

# Acidic: ≥5 Glu/Asp at +1 to +6 (acidic)
mot[c("minus","pos")]<-str_split_fixed(mot$motif,"ST",2)
mot$pos_loc<-substr(mot$pos,2,7)
mot$de<-str_count(mot$pos_loc,"D")+str_count(mot$pos_loc,"E")
mot$acidic_1<-ifelse(mot$de>=5,TRUE,FALSE)

# Group 2: Basophilic  Arg/Lys at −3 (basic)
mot$basic<-(str_sub(mot$motif, 5, 5) =="R") |
  (str_sub(mot$motif, 5, 5)=="K")

# Group 3: acidophilic Glu/Asp at +1/+2 or +3 (acidic),
mot$acidic_2<- 
  (str_sub(mot$motif, 12, 12) %in% c("E", "D"))| 
  (str_sub(mot$motif, 13, 13) %in% c("E", "D"))| 
  (str_sub(mot$motif, 14, 14) %in% c("E", "D"))

# Group 2: Basophilic ≥2 Arg/Lys at −6 to −1 (basic)
mot$minus_loc<-substr(mot$minus,2,7)
mot$rk<-str_count(mot$minus_loc,"R")+str_count(mot$minus_loc,"K")
mot$basic_2<-ifelse(mot$rk>1,TRUE,FALSE)

# Group 3: Other

# Divide up the Other group
# Group 3.1: Any E or D residues - acidic (but less strict than the acidophilic definition)
mot$DE_residues<-str_detect(mot$motif,"E")|str_detect(mot$motif,"D")
# Group 3.2: Any N, Q, S or T (more than 1 for ST as this is the [ST]) - hydrophilic neutral aa
mot$NQST_residues<-str_detect(mot$motif,"N")|str_detect(mot$motif,"Q")|str_count(mot$motif,fixed("S"))>1| (str_count(mot$motif,fixed("T"))>1)
# Group 3.4: Any KR residues  (basophilic)
mot$KR_residues<-str_detect(mot$motif,"K")|str_detect(mot$motif,"R")



# Decision tree
# is it proline directed?, is it acidophilic?, is it basophilic?, is it acidophilic?, is it basophilic?, does it contain aa with polar uncharged side chains?, does it contain DE residues (acidic)? does it contain basic aa?
mot$group<-ifelse(mot$pro,"Proline-directed",
                  ifelse(mot$acidic_1,"Acidophilic",
                         ifelse(mot$basic,"Basophilic",
                                ifelse(mot$acidic_2,"Acidophilic",
                                       ifelse(mot$basic_2,"Basophilic",
                                              ifelse(mot$NQST_residues,"NQST_residues",
                                                     ifelse(mot$DE_residues,"DE_residues",
                                                            ifelse(mot$KR_residues,"KR_residues",
                                                                   "Other"))))))))





#################
# 3. summary    #
#               #
#################

groups<-c("Acidophilic","Basophilic","Proline-directed","NQST_residues","DE_residues","KR_residues", "Other") 
table(mot$group)

# for each class
for (motif_group in groups){
  # filter rmotifx result
  mot_class<-subset(mot, group == motif_group)
  print(motif_group)
  print(mot_class$motif)
  # in each motif class, get the proteins in the foreground
  protein_per_mot_class<-read.csv(file=paste0("01_motif_and_pathway_analysis/outputs/gold_motif_seqs/",motif_group,"_gold_motif_seqs.csv"))
  # unique proteins that have motif
  print(protein_per_mot_class %>% 
    distinct(Protein) %>% 
    count())
  
}
