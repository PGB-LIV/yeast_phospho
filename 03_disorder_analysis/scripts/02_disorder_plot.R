#############################################
# Disorder plots                            #
# input: from 01_metapredict_processing.py  # 
#############################################
library(ggplot2)
library(ggpubr)
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

# combine
overall<-rbind(phospho,other_STY)
p1<-ggplot(data=overall, mapping=aes(x=PTM_residue, y=disorder,fill=aa))+geom_boxplot()+theme_bw()+labs(title="A. pSTY versus STY")+
  scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score") + xlab("Residue")+ geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)

# Gold-Silver-Bronze
phospho$PTM_FLR_category<-factor(phospho$PTM_FLR_category,levels=c("Gold","Silver","Bronze"))
p2<-ggplot(data=phospho, mapping=aes(x=PTM_residue, y=disorder,fill=aa))+geom_boxplot()+theme_bw()+labs(title="B. pSTY in Gold-Silver-Bronze")+facet_wrap(~PTM_FLR_category)+
  scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score")+xlab("Residue")+geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)


ggarrange(p1,p2,ncol=1, common.legend = TRUE)
ggsave(paste("03_disorder_analysis/outputs/yeast_disorder_pSTY.png",sep=""),dpi=330, height=6, width=6)

# disorder score stats
pY<-subset(phospho,PTM_residue=="pY")
quantile(pY$disorder)
pS<-subset(phospho,PTM_residue=="pS")
quantile(pS$disorder)
pT<-subset(phospho,PTM_residue=="pT")
quantile(pT$disorder)

Y<-subset(other_STY,PTM_residue=="Y")
quantile(Y$disorder)
S<-subset(other_STY,PTM_residue=="S")
quantile(S$disorder)
T<-subset(other_STY,PTM_residue=="T")
quantile(T$disorder)

