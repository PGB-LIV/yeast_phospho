############################## 
# Count of sites per prot    #
##############################
library(dplyr)
library(ggplot2)
safe_colorblind_palette <- c("#888888","#6699CC","#117733", "#CC6677")


counts_A<-read.csv("04_AF3/inputs/G2S1B_0.05_protein_pos_all_prot_mapping.csv")
# remove alanine
counts_noA<-filter(counts_A, PTM_residue!="A")
# Group per FLR category
sum <- counts_noA %>% group_by(PTM_FLR_category, Protein) %>% 
  summarise(prot_count = n())
sum$prot_count[sum$prot_count>=10] = ">=10"
sum$prot_count<-factor(sum$prot_count, levels=c("1","2","3","4","5","6","7","8","9",">=10"))

sum$PTM_FLR_category<-factor(sum$PTM_FLR_category,levels=c("Gold","Silver","Bronze"))
# Site count per protein in GSB
x<-ggplot(sum, aes(x=prot_count)) + 
  geom_bar()+facet_wrap(~PTM_FLR_category)+theme_bw()+xlab("Phosphosite count")+ylab("Number of proteins")
ggsave("phosphosite_freq_GSB.png",plot=x,dpi=330,width=14, height=5)

multiple_gold_site_protein <- subset(sum, PTM_FLR_category=="Gold"&prot_count==">=10")
write.csv(multiple_gold_site_protein,"04_AF3/inputs/multiple_gold_sites.csv",row.names=F)  

# Overall count per protein
sum_overall <- counts_noA %>% group_by(Protein) %>% 
  summarise(prot_count = n())
sum_overall$prot_count[sum_overall$prot_count>=10] = ">=10"
sum_overall$prot_count<-factor(sum_overall$prot_count, levels=c("1","2","3","4","5","6","7","8","9",">=10"))

x<-ggplot(sum_overall, aes(x=prot_count)) + 
  geom_bar()+theme_bw()+xlab("Phosphosite count")+ylab("Number of proteins")
ggsave("04_AF3/outputs/phosphosite_freq_overall.png",plot=x,dpi=330,width=6)

# Compare 1 site versus >1 site
sum_overall <- counts_noA %>% group_by(Protein) %>% 
  summarise(prot_count = n())
sum_overall$prot_count[sum_overall$prot_count>1] = ">1"
sum_overall$prot_count<-factor(sum_overall$prot_count, levels=c("1",">1"))
x<-ggplot(sum_overall, aes(x=prot_count)) + 
  geom_bar()+theme_bw()+xlab("Phosphosite count")+ylab("Number of proteins")
ggsave("04_AF3/outputs/phosphosite_freq_binary.png",plot=x,dpi=330,width=6)
