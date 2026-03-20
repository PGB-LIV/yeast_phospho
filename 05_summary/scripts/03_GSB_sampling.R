##############################################
# This script is for testing the             #
# saturation of phosphosites                 #
# by randomly sampling 5, 6 or 7 datasets    #
# and comparing Gold sites                   #
##############################################
library(reshape2)
library(dplyr)
library(ggplot2)

samples <- c("sample_5_1","sample_5_2","sample_5_3","sample_6_1", "sample_6_2","sample_6_3",
             "sample_7_1", "sample_7_2", "sample_7_3")

# All 8 datasets
full <- read.csv("05_summary/inputs/G2S1B_0.05_protein_pos_all_prot_mapping.csv")
full <- select(full, PTM_FLR_category, PTM_residue)
# Gold sites
full <- filter(full, PTM_FLR_category == "Gold")
overall_df <- melt(table(full))
overall_df$group <-"Full_sample"

# Loop through the random samples generated in the py script
for (sample in samples){
  sample_df <- read.csv(paste("05_summary/inputs/G2S1B_0.05_protein_pos_all_prot_mapping_", sample, ".csv", sep = ""))
  sample_df <- select(sample_df, PTM_FLR_category,PTM_residue)
  # Gold sites
  sample_df <- filter(sample_df, PTM_FLR_category == "Gold")
  sample_df <- melt(table(sample_df))
  sample_df$group <- sample
  # Join counts
  overall_df <- rbind(overall_df, sample_df)
  
}

# order of group for plot
level_order <- c("sample_5_1", "sample_5_2", "sample_5_3", "sample_6_1", "sample_6_2", "sample_6_3",
               "sample_7_1", "sample_7_2", "sample_7_3", "Full_sample")
###########################
# plot - STYA comparison  #
###########################
safe_colorblind_palette <- c("#888888", "#6699CC", "#117733", "#CC6677")

stya_comparison <- ggplot(overall_df, aes(fill = PTM_residue, y = value, x = factor(group, level = level_order))) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  geom_text(aes(label = value), size = 4, vjust = -0.4, position = position_dodge(0.9)) +
  ylab("Count of Gold sites") +
  xlab("Sample") + theme(text = element_text(size = 18)) +
  scale_fill_manual(values = safe_colorblind_palette, name = "Phosphosite residue") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) + scale_y_continuous(labels = scales::comma)
stya_comparison
ggsave("05_summary/outputs/GSB_sampling_STYA.png", plot = stya_comparison, dpi = 330, width = 18)


#################################
# plot - Gold sites comparison  #
#################################
# Filter STY sites
sty_sites <- filter(overall_df, PTM_residue == "S" | PTM_residue == "T" | PTM_residue == "Y")
sty_sites <- select(sty_sites, group, value)

# count STY per sample
gold_count <- sty_sites %>% 
  group_by(group) %>% 
  summarise(sum_group = sum(value))

gold_comparison <- ggplot(gold_count, aes(y = sum_group, x = factor(group, level = level_order))) + 
  geom_bar(position = 'dodge', stat = 'identity', fill = "#882255") +
  geom_text(aes(label = sum_group), size = 4, vjust = -0.4, position = position_dodge(0.9)) +
  ylab("Count of Gold STY sites") +
  xlab("Sample") + theme(text = element_text(size = 18)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) + scale_y_continuous(labels = scales::comma)


gold_comparison
ggsave("05_summary/outputs/GSB_sampling_gold.png", plot = gold_comparison, dpi = 330, width = 18)
