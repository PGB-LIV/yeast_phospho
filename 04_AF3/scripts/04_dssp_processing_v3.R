#############################################################
# This script is for DSSP outputs                           #
#                   YEAST PHOSPHO                           #
#                                                           #
#############################################################
# confidence score - STY only, mean confidence score
# ptm - based on struc (filtered for phospho). exploded on PTM column (so structure row might not be unique
# if more than one ptm site occurs within a given structral range)
# conf_merge - ptm (before explode so contains non-ptm sites andf one row per struct range) plus confidence - atom id in conf has to be within struc range 
# protein - per protein aa count per struct element
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(corrplot)
library(chisq.posthoc.test)
library(stringr)
#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
safe_colorblind_palette <- c("#6699CC","#117733", "#CC6677", "#888888")
#############################
# Plot 1                    # 
# AA in struct (whole prot) #
#############################  

df<-read.csv(file=paste0("04_AF3/inputs/summary_per_protein_secondary_struct_features_dssp.csv"))

# adding zeros
df <- df%>%
  complete(
    nesting(Protein,PTM_model),
    structure= c("BEND", "HELX_LH_PP_P", "HELX_RH_3T_P", "HELX_RH_AL_P",
                        "HELX_RH_PI_P", "OTHER", "STRN", "TURN_TY1_P" ), fill = list(aa_total = 0))

p1<-ggplot(df,aes( y=aa_total, x=factor(structure),fill=PTM_model))+ geom_boxplot()+theme_bw()+
  scale_fill_manual(values=safe_colorblind_palette[2:3])+xlab("Structure")+ylab("Amino acid count")
p1
ggsave("04_AF3/outputs/secondary_structural_boxplot.png",
       dpi=330,width=14, height=10)


###########################################
# Plot 2                                  #
# proportion of aa in struct of protein   #
###########################################

# aa count for each protein - group by protein and PTM model
structure_sum_prot<-df%>% 
  group_by(Protein,PTM_model) %>% 
  summarise(total_aa_in_protein = sum(aa_total))

# add overall count
df_overall<-left_join(df, structure_sum_prot)
# calc proportion
df_overall$proportion<-df_overall$aa_total/df_overall$total_aa_in_protein

p2<-ggplot(df_overall,aes( y=proportion, x=factor(structure),fill=PTM_model))+ geom_boxplot()+theme_bw()+
  scale_fill_manual(values=safe_colorblind_palette[2:3])+ylab("Proportion of amino acids")+xlab("Structure")

p2
ggsave("04_AF3/outputs/secondary_structural_proportion_boxplot.png",
       dpi=330,width=14,height=10)


###################
# Plot 3          #
#                 #
###################

# add up for  structure
structure_sum<-df %>% 
  group_by(PTM_model,structure) %>% 
  summarise(aa_per_struct = sum(aa_total))


p3<-ggplot(structure_sum,aes( y=aa_per_struct, x=factor(structure),fill=PTM_model))+ 
  geom_bar(position='dodge', stat='identity')+
  ylab("Amino acid count")+
  xlab("Structure")+
  labs(fill="Model")+theme(text = element_text(size=10))+
  scale_fill_manual(values=safe_colorblind_palette)+ggtitle("Secondary structural elements (protein level)")+
  theme_bw()
p3

ggsave("04_AF3/outputs/secondary_structural_barchart_count_per_structure_prot_level.png",
       dpi=330,width=14,height=8)

#################
# Chi-sq test   #
#################
#https://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
#https://cran.r-project.org/web/packages/chisq.posthoc.test/vignettes/chisq_posthoc_test.html
#https://www.geeksforgeeks.org/r-language/chi-square-test-in-r/
#https://stackoverflow.com/questions/73496335/how-to-add-the-stars-to-my-barplot-in-ggplot-based-on-the-values-in-the-padj
# long-> wide
df_wide <- pivot_wider(structure_sum, names_from = PTM_model, values_from = aa_per_struct)
#col1 -> rowname
df_wide <- data.frame(df_wide, row.names = 1)
chisq_res<-chisq.test((df_wide))
chisq_posthoc<-chisq.posthoc.test((df_wide),
                   method = "bonferroni")

chisq_pval<-subset(chisq_posthoc, Value=="p values")
chisq_df <- chisq_pval %>%
  mutate(
    signif_label = case_when(
      PTM < 0.001 ~ "***",
      PTM < 0.01  ~ "**",
      PTM < 0.05  ~ "*",
      PTM >= 0.05 ~ ""
    )
  )

chisq_df$y_pos<-c(50000,50000,50000,220000,50000,255000,50000,60000)

safe_colorblind_palette<-c("#332288","#88CCEE")
p4<-ggplot(structure_sum,aes(y=aa_per_struct, x=factor(structure),fill=PTM_model))+ 
  geom_bar(position='dodge', stat='identity')+
  ylab("Amino acid count")+
  xlab("Structure")+
  labs(fill="Model")+theme(text = element_text(size=18))+
  scale_fill_manual(values=safe_colorblind_palette)+ggtitle("Secondary structural elements (protein level)")+
  theme_bw()+
  geom_text(data=chisq_df,aes(x=Dimension, y=y_pos,label=signif_label), inherit.aes = FALSE, size = 3.5) + 
  scale_x_discrete(labels=c("Bend","Left-handed alpha-helix","3-10 helix", "Alpha-helix", "Pi-helix", "Loop", "Beta-bridge/Strand", "Turn")) +
  scale_y_continuous(labels=scales::comma)
p4

ggsave("04_AF3/outputs/chi_sq_barplot.png",
       dpi=330,width=12,height=6)
##########################
# Plot 5                 #   
# 2.aa per ptm position  #
##########################
df<-read.csv("04_AF3/inputs/summary_per_ptm_pos_protein_secondary_struct_features_dssp.csv")

# Frequency of secondary structure at PTM site in PTM model versus no PTM model  
df<-df %>%
  group_by(PTM_model, X_struct_conf.conf_type_id) %>%
  mutate(structure_count = n()) %>% ungroup()
#t<-subset(df,PTM_model=="PTM")
#t<-subset(df,PTM_model=="no_PTM")
#table(t$X_struct_conf.conf_type_id)

ptm_struc<-subset(df, select=c(PTM_model,X_struct_conf.conf_type_id,structure_count))
ptm_struc<-distinct(ptm_struc)
p5<-ggplot(ptm_struc,aes(x=factor(X_struct_conf.conf_type_id),fill=PTM_model,y=structure_count))+ 
  geom_bar(position="dodge", stat='identity')+theme_bw()+xlab("Structure")+ylab("Frequency at PTM site")+
  scale_fill_manual(values=safe_colorblind_palette)
p5
ggsave("04_AF3/outputs/structure_at_PTM_site.png",dpi=330, width=12,height=9)

# -> proportion
sum_prot<-ptm_struc%>% 
  group_by(PTM_model) %>% 
  summarise(ptm_sites_total = sum(structure_count))

# add overall count
df_overall<-left_join(ptm_struc, sum_prot)
df_overall$proportion<-df_overall$structure_count/df_overall$ptm_sites_total
p5<-ggplot(df_overall,aes(x=factor(X_struct_conf.conf_type_id),y=proportion, fill=PTM_model))+ 
  geom_bar(position="dodge", stat='identity')+theme_bw()+xlab("Structure")+ylab("Proportion at PTM site")+
  scale_fill_manual(values=safe_colorblind_palette)+ylim(0,1)
p5
####################################################
# Plot 6                                           #   
# Are we confident in the predictions being made?  #
#                                                  #
####################################################

# df containing plddt scores per residues -PTM/no PTM model and phospho/non-phospho sites
confidence_score<-read.csv("04_AF3/inputs/confidence_score_AF.csv")
confidence_score<-subset(confidence_score,select=c(Protein,X_atom_site.auth_seq_id,X_atom_site.label_comp_id,phospho,PTM_model,mean_confidence_score))

##################
# a. Serine plot #
##################

# Filter for serine (SER) and phosphoserine (SEP)
serine<-subset(confidence_score, X_atom_site.label_comp_id=="SER"|X_atom_site.label_comp_id=="SEP")
# col -> PTM model/ phosphosite?
serine$group<-paste0(serine$PTM_model,"_",serine$phospho)


# Re-order
order<-c("no_PTM_False","PTM_False","no_PTM_True","PTM_True")
p6<-ggplot(serine,aes( y=mean_confidence_score, x=factor(group, level=order),fill=PTM_model))+ geom_boxplot()+
  theme_bw()+ylab("pLDDT mean atom score")+xlab("Residue")+ggtitle("Serine")+
  scale_fill_manual(values=safe_colorblind_palette)+scale_x_discrete(labels=c("S", "S", "pS", "pS"))
p6
ggsave("04_AF3/outputs/confidence_Ser_boxplot.png",dpi=330, height=5, width=5)


#####################
# b. Threonine plot #
#####################

# Filter for threonine (THR) and phosphothreonine (TPO)
thr<-subset(confidence_score, X_atom_site.label_comp_id=="THR"|X_atom_site.label_comp_id=="TPO")
# column -> PTM model/phosphosite?
thr$group<-paste0(thr$PTM_model,"_",thr$phospho)

p7<-ggplot(thr,aes( y=mean_confidence_score, x=factor(group,level=order),fill=PTM_model))+ geom_boxplot()+
  theme_bw()+ylab("pLDDT mean atom score")+xlab("Residue")+ggtitle("Threonine")+
  scale_fill_manual(values=safe_colorblind_palette)+scale_x_discrete(labels=c("T", "T", "pT", "pT"))
p7
ggsave("04_AF3/outputs/confidence_Thr_boxplot.png",dpi=330, height=5, width=5)


#####################
# c. Tyr plot       #
#####################

# Filter for tyrosine and phospho tyr
tyr<-subset(confidence_score, X_atom_site.label_comp_id=="TYR"|X_atom_site.label_comp_id=="PTR")
# column -> PTM model/phosphosite?
tyr$group<-paste0(tyr$PTM_model,"_",tyr$phospho)

p8<-ggplot(tyr,aes( y=mean_confidence_score, x=factor(group,level=order),fill=PTM_model))+ geom_boxplot()+
  theme_bw()+ylab("pLDDT mean atom score")+xlab("Residue")+ggtitle("Tyrosine")+
  scale_fill_manual(values=safe_colorblind_palette)+scale_x_discrete(labels=c("Y", "Y", "pY", "pY"))
p8
ggsave("04_AF3/outputs/confidence_Tyr_boxplot.png",dpi=330, height=5, width=5)





##########################
# Plot 9                 #     
# confidence per struct  #
#                        #  
##########################
df<-read.csv("04_AF3/inputs/conf_merge.csv")
df$group<-paste0(df$PTM_model,"_",df$phospho)

ser<-subset(df, X_atom_site.label_comp_id=="SER"|X_atom_site.label_comp_id=="SEP")
safe_colorblind_palette <- c("#888888","#44AA99", "#882255","#88CCEE")
p9<-ggplot(ser,aes( y=mean_confidence_score, x=factor(X_struct_conf.conf_type_id),fill=group))+ geom_boxplot()+
  theme_bw()+ylab("pLDDT mean atom score")+xlab("Structure")+ggtitle("A. Serine")+
  scale_fill_manual(values=safe_colorblind_palette,name = "Model", 
                    labels = c("No PTM model/ not a phosphosite", "No PTM model/ phosphosite","PTM model/ not a phosphosite",
                               "PTM model/ phosphosite")) + scale_x_discrete(labels=c("Bend","Left-handed alpha-helix","3-10 helix",
                                                                                      "Alpha-helix", "Pi-helix", "Loop", "Beta-bridge/Strand",
                                                                                      "Turn"))+theme(text = element_text(size=18))
p9
ggsave("04_AF3/outputs/confidence_ser_per_struct.png",dpi=330,width=14, height=5)

thr<-subset(df, X_atom_site.label_comp_id=="THR"|X_atom_site.label_comp_id=="TPO")
p10<-ggplot(thr,aes( y=mean_confidence_score, x=factor(X_struct_conf.conf_type_id),fill=group))+ geom_boxplot()+
  theme_bw()+ylab("pLDDT mean atom score")+xlab("Structure")+ggtitle("B. Threonine")+
  scale_fill_manual(values=safe_colorblind_palette,name = "Model", 
                    labels = c("No PTM model/ not a phosphosite","No PTM model/ phosphosite","PTM model/ not a phosphosite",
                               "PTM model/ phosphosite")) + scale_x_discrete(labels=c("Bend","Left-handed alpha-helix","3-10 helix",
                                                                                     "Alpha-helix", "Pi-helix", "Loop", "Beta-bridge/Strand",
                                                                                     "Turn"))+theme(text = element_text(size=18))

p10
ggsave("04_AF3/outputs/confidence_thr_per_struct.png",dpi=330,width=14,height=5)

tyr<-subset(df, X_atom_site.label_comp_id=="TYR"|X_atom_site.label_comp_id=="PTR")
p11<-ggplot(tyr,aes( y=mean_confidence_score, x=factor(X_struct_conf.conf_type_id),fill=group))+ geom_boxplot()+
  theme_bw()+ylab("pLDDT mean atom score")+xlab("Structure")+ggtitle("Tyrosine")+
  scale_fill_manual(values=safe_colorblind_palette,name = "", labels = c("No PTM model/ not a phosphosite","No PTM model/ phosphosite","PTM model/ not a phosphosite","PTM model/ phosphosite"))

p11
ggsave("04_AF3/outputs/confidence_tyr_per_struct.png",dpi=330,width=14,height=5)

# combine - ST only- we are exploring the disorder to order transition (from metapredict we know Tyrosine residues tend to be ordered in the non-phospho state)
overall_plot <- p9 / p10 + plot_layout(axis_titles= "collect", guides= "collect")
ggsave("04_AF3/outputs/confidence_per_structure_STY.png",plot= overall_plot,dpi=330, width=18,height=10)



####################################
# confidence diff - ptm vs no ptm  #
#                                  #
####################################

# compare confidence of all STY in PTM model versus no PTM model
confidence_score<-read.csv("04_AF3/inputs/confidence_score_AF.csv")
confidence_score<-subset(confidence_score,select=c(Protein,X_atom_site.auth_seq_id,X_atom_site.label_comp_id,phospho,PTM_model,mean_confidence_score  ))

confidence_score$aa<-confidence_score$X_atom_site.label_comp_id
confidence_score$aa[confidence_score$aa=="SEP"] <- "SER"
confidence_score$aa[confidence_score$aa=="TPO"] <- "THR"
confidence_score$aa[confidence_score$aa=="PTR"] <- "TYR"


# remove the column containing the phospho representation of aa as we need the results to match for PTM and no PTM models
confidence_score<-subset(confidence_score, select=-c(X_atom_site.label_comp_id))
df_wide<-confidence_score %>% 
  pivot_wider(names_from = PTM_model, values_from = mean_confidence_score)

# 1. plot of change in STY confidence scores
df_wide$diff<- df_wide$PTM - df_wide$no_PTM
#df_wide<-subset(df_wide,phospho=="True")
p12<-ggplot(df_wide,aes( y=diff, x=factor(aa),fill=aa))+ geom_boxplot()+
  theme_bw()+ylab("Difference in pLDDT mean atom score")+xlab("Residue")+ggtitle("PTM versus no PTM model - all STY")+scale_fill_manual(values=safe_colorblind_palette)
p12
ggsave("04_AF3/outputs/confidence_diff_PTM_versus_no_PTM_model_all.png",dpi=330,width=12)


# filter for phospho
df_wide_phospho<-subset(df_wide,phospho=="True")
p13<-ggplot(df_wide_phospho,aes( y=diff, x=factor(aa),fill=aa))+ geom_boxplot()+
  theme_bw()+ylab("Difference in pLDDT mean atom score")+xlab("Residue")+ggtitle("PTM versus no PTM model - phosphosites")+scale_fill_manual(values=safe_colorblind_palette)
p13
ggsave("04_AF3/outputs/confidence_diff_PTM_versus_no_PTM_model_phosphosite.png",dpi=330,width=12)

# filter for phospho
df_wide_no_phospho<-subset(df_wide,phospho=="False")
p14<-ggplot(df_wide_no_phospho,aes( y=diff, x=factor(aa),fill=aa))+ geom_boxplot()+
  theme_bw()+ylab("Difference in pLDDT mean atom score")+xlab("Residue")+ggtitle("PTM versus no PTM model - not phosphosites")+scale_fill_manual(values=safe_colorblind_palette)
p14


#############################
# Plot 15                   #     
#length - ptm vs no ptm     #
#                           # 
#############################
df<-read.csv("04_AF3/inputs/summary_per_ptm_pos_protein_secondary_struct_features_dssp.csv")
df$phospho <- ifelse(is.na(df$PTM_position),"non-phospho", "phospho")
df<-subset(df, select=c(Protein, PTM_position_all, PTM_model, X_struct_conf.conf_type_id, structure_start_end, phospho))
df <- df[!duplicated(df), ]
df[c("start", "end")] <- str_split_fixed(df$structure_start_end, "-", 2)
df$start<-as.numeric(df$start)
df$end<-as.numeric(df$end)
df$aa_count<- df$end - df$start + 1
df$group<-paste0(df$PTM_model,"_",df$phospho)

# stats
alpha_helix<-subset(df, X_struct_conf.conf_type_id == "HELX_RH_AL_P")
alpha_helix_ptm<- subset(alpha_helix, group == "PTM_phospho")
quantile(alpha_helix_ptm$aa_count)
alpha_helix_no_ptm<- subset(alpha_helix, group == "no_PTM_phospho")
quantile(alpha_helix_no_ptm$aa_count)


safe_colorblind_palette<-c("#44AA99", "#88CCEE")
p15<-ggplot(df,aes(y=aa_count, x=factor(X_struct_conf.conf_type_id),fill=group))+ geom_boxplot(outliers = F)+
  theme_bw()+ylab("Length of structure (amino acids)")+xlab("Structure")+
  scale_fill_manual("Model",values=safe_colorblind_palette)+
  scale_x_discrete(labels=c("Bend","Left-handed alpha-helix","3-10 helix", "Alpha-helix", "Pi-helix", "Loop", "Beta-bridge/Strand", "Turn"))+theme(text = element_text(size=18))
p15
ggsave("04_AF3/outputs/length_struct_PTM_site.png",dpi=330, width=18, height=8)

# remove other
df<-subset(df, X_struct_conf.conf_type_id!="OTHER")
p16<-ggplot(df,aes( y=aa_count, x=factor(X_struct_conf.conf_type_id),fill=group))+ geom_boxplot(outliers = F)+
  theme_bw()+ylab("Length of structure (amino acids)")+xlab("Structure")+
  scale_fill_manual("Model", values=safe_colorblind_palette)+
  scale_x_discrete(labels=c("Bend","Left-handed alpha-helix","3-10 helix", "Alpha-helix", "Pi-helix", "Loop", "Beta-bridge/Strand", "Turn"))
p16
ggsave("04_AF3/outputs/length_struct_PTM_site_without_other.png",dpi=330,width=12)



