source("code/import_pellett_sam_probe.R")
source("code/import_weekes_protein.R")

library(viridis)
library(extrafont)

# list of proteins in SV schematic
secretion_schematic_proteins <- c('SYTL2', 'RAB3GAP1', 'RAB3IL1', 'RAB3IP', 'TBC1D10A',
                                  'RAB3A', 'RAB3B', 'RAB8A', 'RAB8B', 'RAB27A',
                                  'SYTL4', 'STXBP1', 'STXBP2', 'SNAP25',
                                  'STX2', 'VAMP2', 'SYT3', 'STX3', 'VAMP8')



# curating the pellett transcript data --------------------------------------------

# selecting out all of the rows corresponding to the schematic
secretion_transcripts <- pellett_probe[pellett_probe$SYMBOL %in% secretion_schematic_proteins,]

# removing columns that aren't being used for heatmap
col_to_remove <- c("\\w3\\..*", "\\w4\\..*", "\\w5\\..*", "\\w7\\..*", "\\w9\\..*", "\\w10\\..*", "\\w12\\..*")
remove_col <- grep(paste(col_to_remove, collapse = "|"), colnames(secretion_transcripts))
secretion_transcripts[colnames(secretion_transcripts)[c(remove_col)]] <- NULL

# removing rows without any detected data
secretion_transcripts <- secretion_transcripts[!is.na(secretion_transcripts$A1.AVG_Signal), ]

# selecting cols with data
secretion_transcripts_data <- filter(secretion_transcripts[,c(2,grep("AVG_Signal", colnames(secretion_transcripts)))])

# combining the data from probes matching the same transcript/protein
secretion_transcripts_data <- as_data_frame(aggregate(. ~ SYMBOL, data = secretion_transcripts_data, FUN = mean))

# tidying the data and adding the replicate numbers
secretion_transcripts_data <- secretion_transcripts_data %>%
  gather(key = "Sample", value = "Mean", 2:16) %>%
  spread(key = "SYMBOL", value = "Mean") %>%
  add_column(Replicate = rep(1:3, each = 5), .before = 1) %>% 
  mutate(Sample = gsub("[ABC]1\\..*", "NA.Mock", Sample),
         Sample = gsub("[ABC]2\\..*", "Towne.12h", Sample),
         Sample = gsub("[ABC]6\\..*", "AD169.12h", Sample),
         Sample = gsub("[ABC]8\\..*", "Towne.96h", Sample),
         Sample = gsub("[ABC]11\\..*", "AD169.96h", Sample)) %>% 
  separate(col = Sample, into = c("Virus", "Sample"))
        
# log2 transforming the values
secretion_transcripts_log2 <- secretion_transcripts_data %>% 
  mutate_at(vars(-c(Replicate, Virus, Sample)), funs(log2))

# subtracting the mock values from the corresponding samples
secretion_transcripts_log2_fc <- secretion_transcripts_log2 %>% 
  group_by(Replicate) %>% 
  mutate_at(vars(3:21), funs(ifelse(Sample != "Mock", .-.[Sample == "Mock"], .))) %>%
  ungroup() %>% 
  filter(Sample != "Mock") %>%
  select(-Replicate)

# calculating the mean for each gene
secretion_transcripts_log2_fc_mean <- secretion_transcripts_log2_fc %>%
  group_by(Virus, Sample) %>% 
  summarise_all(mean) %>% 
  ungroup()

# calculating the std deviation of each gene
secretion_transcripts_log2_fc_sd <- secretion_transcripts_log2_fc %>% 
  group_by(Virus, Sample) %>% 
  summarise_all(sd) %>% 
  ungroup()



# curating the weekes protein data ----------------------------------------

# selecting out all of the rows corresponding to the schematic
secretion_proteins <- weekes_data_host[weekes_data_host$Symbol %in% secretion_schematic_proteins,]

# removing any proteins that don't have data for both pm and wcl columns
secretion_proteins <- secretion_proteins[!is.na(secretion_proteins$pm2_Mock_1) | !is.na(secretion_proteins$wcl2_Mock_1), ]

# Curating data columns of interest and tidying data
secretion_proteins_data <- secretion_proteins %>% 
  select(c(3,4,5,7,12,13,14,16,21)) %>% 
  gather(key = "Sample", value = "Abund", 2:9) %>% 
  spread(key = "Symbol", value = "Abund") %>% 
  mutate(Sample = gsub("Mock\\_[12]", "Mock", Sample)) %>% 
  separate(col = Sample, into = c("Fraction", "Sample"))

# averaging the mock samples and log2 transforming all of the data
secretion_proteins_log2 <- secretion_proteins_data %>% 
  group_by(Fraction, Sample) %>% 
  summarise_all(mean) %>% 
  ungroup %>% 
  mutate_at(vars(-c(Fraction, Sample)), funs(log2))

# calculating the fold change relative to mock
secretion_proteins_log2_fc <- secretion_proteins_log2 %>% 
  group_by(Fraction) %>% 
  mutate_at(vars(2:16), funs(ifelse(Sample != "Mock", .-.[Sample == "Mock"], .))) %>%
  ungroup() %>% 
  filter(Sample != "Mock")


# plotting the transcriptional data -------------------------------------------------------

# reformatting the transcriptional data for plotting
secretion_transcripts_plot_sd <- secretion_transcripts_log2_fc_sd %>% 
  gather(key = "Gene", value = "sd", 3:21) %>% 
  unite("Virus", Virus, Sample, sep = "_")

secretion_transcripts_plot_data <- secretion_transcripts_log2_fc_mean %>% 
  gather(key = "Gene", value = "log2_fc", 3:21) %>%
  unite("Virus", Virus, Sample, sep = "_") %>% 
  mutate(Label = ifelse(is.na(log2_fc), "ND", NA)) %>% 
  mutate(log2_fc = ifelse(is.na(log2_fc), 0, log2_fc))

# appending the sd data to the plot_data tibble
secretion_transcripts_plot_data <- left_join(secretion_transcripts_plot_data, secretion_transcripts_plot_sd, by = c("Virus", "Gene"))

# setting order of plot based on log2_fc data
secretion_transcripts_plot_data$Gene <- reorder(secretion_transcripts_plot_data$Gene, secretion_transcripts_plot_data$log2_fc)

# checking the plotting conditions
levels(secretion_transcripts_plot_data$Gene)
summary(secretion_transcripts_plot_data)

# plotting the transcriptional data
ggplot(secretion_transcripts_plot_data) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_hline(yintercept = c(-2, -1, 1, 2, 3, 4), color = "grey", size = 0.15) +
  geom_bar(mapping = aes(x = Gene, y = log2_fc, fill = Virus), stat = "identity", width = 0.75,
           position = position_dodge(width = .75), color = "black") +
  geom_errorbar(mapping = aes(x = Gene, ymin = log2_fc - sd, ymax = log2_fc + sd, group = Virus),
                position = position_dodge(width = 0.75), width = 0.5) +
  geom_text(mapping = aes(x = Gene, y = log2_fc, label = Label, group = Virus), angle = 90,
            hjust = -0.25, na.rm = T, size = 2.5, position = position_dodge(width = .75)) +
  scale_fill_viridis(discrete = T, begin = 1, end = 0, option = "magma") +
  scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1), expand = c(0,0)) +
  coord_fixed(ratio = 1) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x = element_line(color = "grey"), axis.ticks.y = element_line(color = "black"),
        axis.line.x = element_line(color = "grey"), axis.line.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(face = "bold", family = "Times New Roman", size = 12),
        legend.position = "top")

# save the plots
#ggsave("figures/secretion_transcripts_plot.png", dpi = 600)
#ggsave("figures/secretion_transcripts_plot.pdf", dpi = 600)


# plotting the protein data -----------------------------------------------

# finding which proteins aren't included in the data so columns can be added for plotting uniformity
# compared to the transcriptional data
# setdiff(secretion_schematic_proteins, unique(secretion_proteins_plot_data$Gene))

# reformatting the protein data for plotting
secretion_proteins_plot_data <- secretion_proteins_log2_fc %>%
  add_column("RAB3IL1" = rep(NA, 4), "RAB3IP" = rep(NA, 4), "SNAP25" = rep(NA, 4), "VAMP8" =  rep(NA, 4)) %>% 
  gather(key = "Gene", value = "log2_fc", 3:21) %>% 
  unite("Fraction", Fraction, Sample, sep = "_") %>%
  mutate(Label = ifelse(is.na(log2_fc), "ND", NA)) %>% 
  mutate(log2_fc = ifelse(is.na(log2_fc), 0, log2_fc))

# setting the order of plotting to be similar to transcriptional plot
secretion_proteins_plot_data$Gene <- factor(secretion_proteins_plot_data$Gene, levels = levels(secretion_transcripts_plot_data$Gene))

# checking the plotting conditions
levels(secretion_proteins_plot_data$Gene)
summary(secretion_proteins_plot_data)

# plotting the protein data
ggplot(secretion_proteins_plot_data) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_hline(yintercept = c(-2, -1, 1, 2, 3, 4), color = "grey", size = 0.15) +
  geom_bar(mapping = aes(x = Gene, y = log2_fc, fill = Fraction), stat = "identity", width = 0.75,
           position = position_dodge(width = .75), color = "black") +
  geom_text(mapping = aes(x = Gene, y = log2_fc, label = Label, group = Fraction), angle = 90,
            hjust = -0.25, na.rm = T, size = 2.5, position = position_dodge(width = .75)) +
  scale_fill_viridis(discrete = T, begin = 1, end = 0, option = "magma") +
  scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1), expand = c(0,0)) +
  coord_fixed(ratio = 1) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x = element_line(color = "grey"), axis.ticks.y = element_line(color = "black"),
        axis.line.x = element_line(color = "grey"), axis.line.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(face = "bold", family = "serif", size = 12),
        legend.position = "top")

# save the plots
#ggsave("figures/secretion_proteins_plot.png", dpi = 600)
#ggsave("figures/secretion_proteins_plot.pdf", dpi = 600)

