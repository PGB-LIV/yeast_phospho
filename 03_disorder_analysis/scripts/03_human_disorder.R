###################################
# human disorder analysis         #
# files from Anton                #
###################################
library(ggplot2)
safe_colorblind_palette <- c("#6699CC","#117733", "#CC6677")


# All sites
human_sty <- read.csv("03_disorder_analysis/inputs/all_disorder_human_psites_STY_only_v2.csv")
human_sty$aa <- human_sty$res
human_sty<-subset(human_sty, select=c(res,disorder,aa))
p1<-ggplot(data=human_sty, mapping=aes(x=res, y=disorder,fill=res))+geom_boxplot()+theme_bw()+labs(title="pSTY")+scale_fill_manual(values=safe_colorblind_palette)


# Gold
human_psty <- read.csv("03_disorder_analysis/inputs/Order_analysis_input.csv")
human_psty <- subset(human_psty, select=c(Res,disorder))
p1 <- ggplot(data=human_psty, mapping=aes(x=Res, y=disorder,fill=Res))+geom_boxplot()+theme_bw()+labs(title="STY")+scale_fill_manual(values=safe_colorblind_palette)

#rename for merge
names(human_psty)[1] <- "res"
human_psty$aa <- human_psty$res
# Add "p" to show phospho
human_psty$res<-paste("p",human_psty$res,sep="")

# join
combined<-rbind(human_psty,human_sty)
p1<-ggplot(data=combined, mapping=aes(x=res, y=disorder,fill=aa))+geom_boxplot()+theme_bw()+labs(title="pSTY versus STY in humans")+scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score") + xlab("Residue")+ 
  geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)

ggsave("03_disorder_analysis/outputs/human_disorder.png",dpi=330, height=4, width=6)

# violin plot
p2<-ggplot(data=combined, mapping=aes(x=res, y=disorder,fill=aa))+geom_violin()+theme_bw()+labs(title="pSTY versus STY in humans")+scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score") + xlab("Residue")+ 
  geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)
ggsave("03_disorder_analysis/outputs/human_disorder_violin_plot.png",dpi=330, height=4, width=6)
