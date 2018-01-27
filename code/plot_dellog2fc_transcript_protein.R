source("code/calc_del_log2_fc.R")

library(ggpmisc)

# gathering all of the dellog2fc data into a single dataframe by joining weekes and gurczynski datasets
protein_vs_transcript_coord <- full_join(weekes_data_host[, c("Symbol", "del_log2_fc_pm2", "del_log2_fc_wcl2")],
                                         sam_probe_data[, c("SYMBOL", "del_log2_fc_towne", "del_log2_fc_ad169")],
                                         by = c("Symbol" = "SYMBOL"))

# importing gene list from GO analysis of transcriptional data
transcript_gene_list <- read.csv("data/transcript_filtered_gene_list.txt", header = F)
transcript_gene_list <- transcript_gene_list$V1
length(unique(transcript_gene_list))

# importing gene list from GO analysis of protein data
protein_gene_list <- read.csv("data/protein_filtered_gene_list.txt", header = F)
protein_gene_list <- protein_gene_list$V1
length(unique(protein_gene_list))

# filtering all rows of dellog2fc dataframe that correspond to genes in transcript
# or protein gene lists
protein_vs_transcript_coord <- protein_vs_transcript_coord %>% 
  filter(Symbol %in% protein_gene_list | Symbol %in% transcript_gene_list)

# number of proteins/transcripts for
# pm2 and towne
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_pm2)) # 631
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_towne)) # 1027
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_pm2) &
      !is.na(protein_vs_transcript_coord$del_log2_fc_towne)) #474


# wcl2 and towne
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_wcl2)) # 1057
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_towne)) # 1027
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_wcl2) &
      !is.na(protein_vs_transcript_coord$del_log2_fc_towne)) # 806


# pm2 and ad169
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_pm2)) # 631
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_ad169)) # 1027
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_pm2) &
      !is.na(protein_vs_transcript_coord$del_log2_fc_ad169)) # 474


# wcl2 and ad169
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_wcl2)) # 1057
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_ad169)) # 1027
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_wcl2) &
      !is.na(protein_vs_transcript_coord$del_log2_fc_ad169)) # 806


# towne and ad169
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_towne)) # 1027
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_ad169)) # 1027
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_towne) &
      !is.na(protein_vs_transcript_coord$del_log2_fc_ad169)) # 1027


# pm and wcl
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_wcl2)) # 1057
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_pm2)) # 631
sum(!is.na(protein_vs_transcript_coord$del_log2_fc_wcl2) &
      !is.na(protein_vs_transcript_coord$del_log2_fc_pm2)) # 580


# data is ready to plot

length(setdiff(unique(protein_gene_list), unique(transcript_gene_list)))

formula <- y ~ x

# protein v towne
ggplot(data = protein_vs_transcript_coord) +
  geom_hline(yintercept = seq(-6,6,1), color = "grey", size = 0.1) +
  geom_vline(xintercept = seq(-6,6,1), color = "grey", size = 0.1) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  geom_point(mapping = aes(x = del_log2_fc_pm2, y = del_log2_fc_towne), color = "firebrick3", alpha = 0.35) +
  geom_point(mapping = aes(x = del_log2_fc_wcl2, y = del_log2_fc_towne), color = "cyan4", alpha = 0.35) +
  geom_smooth(mapping = aes(x = del_log2_fc_pm2, y = del_log2_fc_towne), method = lm, fill = "firebrick3", color = "black", alpha = 0.1, size = 0.5) +
  geom_smooth(mapping = aes(x = del_log2_fc_wcl2, y = del_log2_fc_towne), method = lm, fill = "cyan4", color = "black", alpha = 0.1, size = 0.5) +
  stat_poly_eq(aes(x = del_log2_fc_pm2, y = del_log2_fc_towne, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = formula, parse = TRUE, size = 3) +
  stat_poly_eq(aes(x = del_log2_fc_wcl2, y = del_log2_fc_towne, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.3,
               formula = formula, parse = TRUE, size = 3) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 1), expand = c(0,0)) +
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        text = element_text(face = "bold", family = "serif", size = 12),
        legend.position = "top")

#ggsave("figures/protein_vs_towne_plot.png", dpi = 600)
#ggsave("figures/protein_vs_towne_plot.pdf", dpi = 600)
# pm = 0.16, wcl = 0.34

# protein v ad169
ggplot(data = protein_vs_transcript_coord) +
  geom_hline(yintercept = seq(-6,6,1), color = "grey", size = 0.1) +
  geom_vline(xintercept = seq(-6,6,1), color = "grey", size = 0.1) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  geom_point(mapping = aes(x = del_log2_fc_pm2, y = del_log2_fc_ad169), color = "firebrick3", alpha = 0.35) +
  geom_point(mapping = aes(x = del_log2_fc_wcl2, y = del_log2_fc_ad169), color = "cyan4", alpha = 0.35) +
  geom_smooth(mapping = aes(x = del_log2_fc_pm2, y = del_log2_fc_ad169), method = lm, fill = "firebrick3", color = "black", alpha = 0.1, size = 0.5) +
  geom_smooth(mapping = aes(x = del_log2_fc_wcl2, y = del_log2_fc_ad169), method = lm, fill = "cyan4", color = "black", alpha = 0.1, size = 0.5) +
  stat_poly_eq(aes(x = del_log2_fc_pm2, y = del_log2_fc_ad169, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = formula, parse = TRUE, size = 3) +
  stat_poly_eq(aes(x = del_log2_fc_wcl2, y = del_log2_fc_ad169, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.3,
               formula = formula, parse = TRUE, size = 3) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 1), expand = c(0,0)) +
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        text = element_text(face = "bold", family = "serif", size = 12),
        legend.position = "top")

#ggsave("figures/protein_vs_ad169_plot.png", dpi = 600)
#ggsave("figures/protein_vs_ad169_plot.pdf", dpi = 600)
# pm = 0.13, wcl = 0.26

# towne v ad169
ggplot(data = protein_vs_transcript_coord) +
  geom_hline(yintercept = seq(-6,6,1), color = "grey", size = 0.1) +
  geom_vline(xintercept = seq(-6,6,1), color = "grey", size = 0.1) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  geom_point(mapping = aes(x = del_log2_fc_towne, y = del_log2_fc_ad169), color = "firebrick3", alpha = 0.35) +
  geom_smooth(mapping = aes(x = del_log2_fc_towne, y = del_log2_fc_ad169), method = lm, fill = "firebrick3", color = "black", alpha = 0.1, size = 0.5) +
  stat_poly_eq(aes(x = del_log2_fc_towne, y = del_log2_fc_ad169, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = formula, parse = TRUE, size = 3) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 1), expand = c(0,0)) +
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        text = element_text(face = "bold", family = "serif", size = 12),
        legend.position = "top")

#ggsave("figures/towne_vs_ad169_plot.png", dpi = 600)
#ggsave("figures/towne_vs_ad169_plot.pdf", dpi = 600)
# r = 0.87

# pm vs wcl
ggplot(data = protein_vs_transcript_coord) +
  geom_hline(yintercept = seq(-6,6,1), color = "grey", size = 0.1) +
  geom_vline(xintercept = seq(-6,6,1), color = "grey", size = 0.1) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  geom_point(mapping = aes(x = del_log2_fc_wcl2, y = del_log2_fc_pm2), color = "firebrick3", alpha = 0.35) +
  geom_smooth(mapping = aes(x = del_log2_fc_wcl2, y = del_log2_fc_pm2), method = lm, fill = "firebrick3", color = "black", alpha = 0.1, size = 0.5) +
  stat_poly_eq(aes(x = del_log2_fc_wcl2, y = del_log2_fc_pm2, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = 0.15,
               formula = formula, parse = TRUE, size = 3) +
  scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 1), expand = c(0,0)) +
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 1), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        text = element_text(face = "bold", family = "serif", size = 12),
        legend.position = "top")

#ggsave("figures/pm_vs_wcl_plot.png", dpi = 600)
#ggsave("figures/pm_vs_wcl_plot.pdf", dpi = 600)
# r = 0.39

