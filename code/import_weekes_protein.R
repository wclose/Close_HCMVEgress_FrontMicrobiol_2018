library(tidyverse)

# read in the weekes data file for experiments pm2 and wcl2
weekes_data <- read_tsv("data/weekes_raw.txt")

# filter out any HCMV genes because only care about cellular genes
weekes_data_host <- weekes_data %>%
  filter(weekes_data$Species != "HCMVM")

# ready to calculate the dellog2fc for the data (uses del_log2_fc.R)