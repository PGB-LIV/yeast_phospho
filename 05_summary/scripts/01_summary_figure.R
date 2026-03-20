########################
# Yeast manuscript     #  
# Figure 1             #
######################## 
# https://github.com/PGB-LIV/Rice_Phospho_Manuscript/blob/main/Figs/Rice_fig1.R

library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(dplyr)
library(stringr)
library(patchwork)

safe_colorblind_palette <- c("#6699CC","#117733", "#AA4499","#CC6677","#DDCC77")

########################################################################
#a)bar chart - count non-redundant and redundant phosphopep and sites  #
########################################################################


#############################################
# Collapse experiments to PXD level (sum?)  #
#############################################

##################
# A. PSM counts  #
##################
# Aggregate per PXD
FLR_counts<-read.csv("05_summary/inputs/FLRcounts_noA_decoy_methods.csv",check.names=FALSE,header=F)
#rm file path
FLR_counts$V1<-stringr::str_remove(FLR_counts$V1,"/mnt/hc-storage/users/hleboswe/")
colnames(FLR_counts)<-FLR_counts[1,]
FLR_counts<-FLR_counts[-1,]
FLR_counts<-select(FLR_counts, -c("No. RAW files","Spectral Count","FLR0.01 combined probability cutoff","FLR0.05 combined probability cutoff" ,"FLR0.1 combined probability cutoff","NA",
                                  "Peptidoform-site count", "FLR pA0.01 Peptidoform-site Count", "FLR pA0.05 Peptidoform-site Count",
                                  "FLR pA0.1 Peptidoform-site Count" ))
# wide->long
FLR_counts<-melt(FLR_counts, id.vars="Dataset", variable.name="Count")
# numeric
FLR_counts$value<-as.numeric(FLR_counts$value)
# separate col for PXD and exp
FLR_counts[c('PXD', 'exp')]  <- str_split_fixed(FLR_counts$Dataset, '/', 2)
# Add counts per PXD for each category
df2 <- FLR_counts %>% group_by(PXD,Count) %>%
  summarize(sum = sum(value))
# Rename columns for later merging
colnames(df2)<-c("Source.Dataset.Identifier", "Count","sum")
# Reorder columns for later merging
df2<-df2[,c("Source.Dataset.Identifier","sum","Count")]
# Add facet column
df2$Group<-"PSM"
# Rename for legend
levels(df2$Count)
levels(df2$Count)<-c("FDR 0.01 PSM count", "FDR 0.01 phosphopeptide count")


##########################
# B. Peptidoform counts  #  
##########################
FLR<-read.csv("05_summary/inputs/all_datasets_merged_Site_Peptidoform_centric_Uniprot.tsv",sep="\t")
table(FLR$Source.Dataset.Identifier)

# DF of all peptidoforms per PXD (including decoy sites)
overall<-as.data.frame(table(FLR$Source.Dataset.Identifier))
colnames(overall)<-c("Source.Dataset.Identifier", "sum")
overall$Count<-"Peptidoform-site count"

# Remove decoy sites
FLR<-subset(FLR,Decoy.Modification.Site==0)

# DF of counts of peptidoforms at 5% FLR
peptido_5<-FLR%>% group_by(Source.Dataset.Identifier)%>%
  summarize(sum=sum(Site.Passes.Threshold..0.05.))
peptido_5$Count<-"FLR 0.05 peptidoform-site count"

# DF of counts of peptidoforms at 1% FLR
peptido_1<-FLR%>% group_by(Source.Dataset.Identifier)%>%
  summarize(sum=sum(Site.Passes.Threshold..0.01.))
peptido_1$Count<-"FLR 0.01 peptidoform-site count"

# Combine peptidoform counts
peptidoform_overall<-rbind(overall,peptido_5,peptido_1)
# col for facet
peptidoform_overall$Group<-"Peptidoform"
# Add PSM and peptidoform counts
final<-rbind(df2,peptidoform_overall)


psm<-ggplot2::ggplot(df2, aes(fill=Count, y=as.numeric(sum), x=Source.Dataset.Identifier)) + geom_bar(position='dodge', stat='identity')+
  theme(axis.text.x = element_text(angle = 45 , hjust=1))+ggtitle("A.")+facet_grid(factor(Group, levels=c("PSM","Peptidoform")) ~ ., scales="free_y")+
  scale_fill_manual(values=safe_colorblind_palette)+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(), # https://stackoverflow.com/questions/54583624/combining-plots-that-share-common-x-axis-but-not-y (remove PXDs)
    panel.grid.minor = element_blank())+xlab("Data set")+ylab("Count")+theme(text = element_text(size=11))


peptido<-ggplot2::ggplot(peptidoform_overall, aes(fill=Count, y=as.numeric(sum), x=Source.Dataset.Identifier)) + geom_bar(position='dodge', stat='identity')+
  theme(axis.text.x = element_text(angle = 45 , hjust=1))+
  ggtitle("B.")+
  scale_fill_manual(values=safe_colorblind_palette[3:7], name="Peptidoform FLR category")+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+xlab("Data set")+ylab("Count")+theme(text = element_text(size=18))+
    scale_y_continuous(labels=scales::comma)


# patchwork plot
b<-psm/peptido+plot_layout(axis_titles="collect", axes="collect_x")

#https://tidytales.ca/snippets/2022-12-22_patchwork-shared-axis-labels/#shared-y-axis-labels
#https://www.data-imaginist.com/posts/2024-01-05-patchwork-1-2-0/

######################################################
#b)Gold silver bronze counts per site, with alanines #
######################################################
counts_A<-read.csv("05_summary/inputs/G2S1B_0.05_protein_pos_all_prot_mapping.csv")
safe_colorblind_palette <- c("#888888","#6699CC","#117733", "#CC6677")
cat_counts<-subset(counts_A,select=c("PTM_FLR_category","PTM_residue"))
cat_counts<-melt(table(cat_counts))

level_order<-c("Bronze","Silver","Gold")

gsb_summary<-ggplot(cat_counts,aes(fill=PTM_residue, y=value, x=factor(PTM_FLR_category,level=level_order)))+ 
  geom_bar(position='dodge', stat='identity')+
  geom_text(aes(label = value),size = 4, vjust = -0.4, position = position_dodge(0.9))+
  ylab("Count of sites")+
  xlab("Category")+
  ggtitle("C.")+theme(text = element_text(size=18))+
  scale_fill_manual(values=safe_colorblind_palette,name="Phosphosite residue")+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+scale_y_continuous(labels=scales::comma)



#################################
# log10(PSM count) by PTM res   #
# GSB                           #
#################################
# remove ala
counts_noA<-filter(counts_A,PTM_residue!="A")
counts_noA$log_PSM<-log10(counts_noA$Sum_of_PSM_counts.5.FLR.)
level_order<-c("Bronze","Silver","Gold")
safe_colorblind_palette <- c("#6699CC","#117733", "#CC6677")

psm_counts<-ggplot(data=counts_noA,aes(y=log_PSM, x=factor(PTM_residue), fill=PTM_residue))+ 
  geom_boxplot()+facet_wrap(~factor(PTM_FLR_category,level=level_order))+theme_bw()+
  scale_fill_manual(values=safe_colorblind_palette,name="Phosphosite residue")+ylab("log10(PSM count)")+ xlab("Residue")+ggtitle("D.")+theme(text = element_text(size=18))


##############
# final plot #
##############

#patchwork
#patch_plot<-b/(c|e)
#ggsave("05_summary/outputs/yeast_phosphobuild_summary.png",plot=patch_plot,dpi=330, width=16,height=12)
patch_plot<-peptido/(gsb_summary|psm_counts)
ggsave("05_summary/outputs/yeast_phosphobuild_summary_JPR.png",plot=patch_plot,dpi=330, width=18,height=10)

#####################
# STY distribution  #
#####################

# Distribution of STY in GSB
table(counts_noA$PTM_FLR_category)
table(counts_noA$PTM_FLR_category)/56694

# gold 
gold <- subset(counts_noA, PTM_FLR_category  == "Gold")
table(gold$PTM_residue)
table(gold$PTM_residue)/22311 #total in gold  

# silver 
silver <- subset(counts_noA, PTM_FLR_category  == "Silver")
table(silver$PTM_residue)
table(silver$PTM_residue)/9890  #total in silver

# bronze 
bronze <- subset(counts_noA, PTM_FLR_category  == "Bronze")
table(bronze$PTM_residue)
table(bronze$PTM_residue)/24493    #total in bronze
