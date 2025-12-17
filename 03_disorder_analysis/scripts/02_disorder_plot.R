#############################################
# Disorder plots                            #
# input: from 01_metapredict_processing.py  # 
#############################################
library(ggplot2)
library(ggpubr)
library(ggsignif)
safe_colorblind_palette <- c("#6699CC","#117733", "#CC6677")


df <- read.csv(paste("03_disorder_analysis/inputs/disorder_processed.csv",sep=""))

# GSB sites
phospho <- subset(df,Phospho=="True")
# Column with the amino acid (not with "pS")
phospho$aa<-phospho$PTM_residue
# Add "p" to show phospho
phospho$PTM_residue<-paste("p",phospho$PTM_residue,sep="")
phospho<-subset(phospho,select=c(Protein,PTM_residue,disorder,PTM_FLR_category,aa))

# All other sites not phosphorylated
other_STY<-subset(df, Phospho!= "True")
# STY
other_STY<-subset(other_STY, PTM_residue=="S"|PTM_residue=="T"|PTM_residue=="Y")
other_STY<-subset(other_STY,select=c(Protein,PTM_residue,disorder,PTM_FLR_category))
other_STY$aa <- other_STY$PTM_residue

######################## 
# disorder score stats #
########################

# median
pY<-subset(phospho,PTM_residue=="pY")
quantile(pY$disorder)
pS<-subset(phospho,PTM_residue=="pS")
quantile(pS$disorder)
pT<-subset(phospho,PTM_residue=="pT")
quantile(pT$disorder)

Y_results<-subset(other_STY,PTM_residue=="Y")
quantile(Y_results$disorder)
S_results<-subset(other_STY,PTM_residue=="S")
quantile(S_results$disorder)
T_results<-subset(other_STY,PTM_residue=="T")
quantile(T_results$disorder)


#Kolmogorov-Smirnov 
ks.test(pY$disorder, Y_results$disorder, exact = TRUE)
ks.test(pS$disorder, S_results$disorder, exact = TRUE)
ks.test(pT$disorder, T_results$disorder, exact = TRUE)


##############
# Plots      #
##############

# combine
overall<-rbind(phospho,other_STY)
p1<-ggplot(data=overall, mapping=aes(x=PTM_residue, y=disorder,fill=aa))+geom_boxplot()+theme_bw()+labs(title="A. pSTY versus STY")+
  scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score") + xlab("Residue")+ geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)+
  geom_signif(comparisons = list(c("S", "pS"),c("T", "pT"), c("Y", "pY")),annotation = c("D= 0.368***","D= 0.327***","D= 0.213***"), y_position=c(1.0,1.04,1.08))


# Gold-Silver-Bronze
phospho$PTM_FLR_category<-factor(phospho$PTM_FLR_category,levels=c("Gold","Silver","Bronze"))
p2<-ggplot(data=phospho, mapping=aes(x=PTM_residue, y=disorder,fill=aa))+geom_boxplot()+theme_bw()+labs(title="B. pSTY in Gold-Silver-Bronze")+facet_wrap(~PTM_FLR_category)+
  scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score")+xlab("Residue")+geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)

ggarrange(p1,p2,ncol=1, common.legend = TRUE)
ggsave(paste("03_disorder_analysis/outputs/yeast_disorder_pSTY.png",sep=""),dpi=330, height=14, width=14)


# violin plot
p3<-ggplot(data=overall, mapping=aes(x=PTM_residue, y=disorder,fill=aa))+geom_violin()+theme_bw()+labs(title="A. pSTY versus STY")+
  scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score") + xlab("Residue")+ geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)
p4<-ggplot(data=phospho, mapping=aes(x=PTM_residue, y=disorder,fill=aa))+geom_violin()+theme_bw()+labs(title="B. pSTY in Gold-Silver-Bronze")+facet_wrap(~PTM_FLR_category)+
  scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score")+xlab("Residue")+geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)
ggarrange(p3,p4,ncol=1, common.legend = TRUE)
ggsave(paste("03_disorder_analysis/outputs/yeast_disorder_pSTY_violin_plot.png",sep=""),dpi=330, height=6, width=6)


# abstract - gold
gold <-subset(phospho,PTM_FLR_category  == "Gold")
abstract_plot<-ggplot(data=gold, mapping=aes(x=PTM_residue, y=disorder,fill=aa))+geom_boxplot()+theme_bw()+facet_wrap(~PTM_FLR_category)+
  scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score")+xlab("Residue")+geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)
ggsave(paste("03_disorder_analysis/outputs/yeast_disorder_abstract.png",sep=""),dpi=330, height=3, width=3)
