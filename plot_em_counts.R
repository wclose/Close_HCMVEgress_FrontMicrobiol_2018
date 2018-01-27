library(tidyverse)
library(extrafont)
library(viridis)

em_count <- tibble(Virus = c("Parent", "UL103+", "UL103-"), Inside = c(28, 62, 51), Outside = c(61, 97, 158))
em_count <- em_count %>% 
  mutate(Total = Inside + Outside)
em_count_percent <- em_count %>%
  mutate(Inside = Inside / Total * 100) %>% 
  mutate(Outside = Outside / Total * 100)
em_count_percent

em_count_v2 <- em_count %>%
  gather(key = "Sample", value = "Count", 2:4) %>% 
  spread(Virus, value = "Count")

em_count_v3 <- em_count_percent %>%
  gather(key = "Sample", value = "Percent", 2:4) %>% 
  filter(Sample != "Total")

chisq.test(em_count_v2[1:2, 2:3])
chisq.test(em_count_v2[1:2, 3:4])
chisq.test(em_count_v2[1:2, c(2,4)])

secretion_proteins_plot_data$Gene <- factor(secretion_proteins_plot_data$Gene, levels = levels(secretion_transcripts_plot_data$Gene))

em_count_v3$Virus <- factor(em_count_v3$Virus, levels = c("Parent", "UL103+", "UL103-"))

ggplot(em_count_v3) +
  geom_bar(mapping = aes(x = Virus, y = Percent, group = Sample, fill = Sample), stat = "identity",
           position = "stack", width = 0.85, color = "black") +
  scale_fill_grey() +
  theme(panel.background = element_rect(fill = "white"),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      text = element_text(face = "bold", family = "Times New Roman", size = 12),
      legend.position = "top") +
  scale_y_continuous(expand = c(0,0))

ggsave("plot_em_counts.png", dpi = 600)
ggsave("plot_em_counts.pdf", dpi = 600)


