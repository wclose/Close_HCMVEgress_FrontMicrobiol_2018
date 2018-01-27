source("code/calc_del_log2_fc.R")

# filter transcript data based on significance
transcript_map_data <- sam_probe_data %>% 
  filter(log2_ad169_12_mock >=1 | log2_ad169_12_mock <= -1 |
           log2_ad169_96_mock >=1 | log2_ad169_96_mock <= -1 |
           log2_towne_12_mock >=1 | log2_towne_12_mock <= -1 |
           log2_towne_96_mock >=1 | log2_towne_96_mock <= -1 |
           del_log2_fc_ad169 >= 1 | del_log2_fc_ad169 <= -1 |
           del_log2_fc_towne >= 1 | del_log2_fc_towne <= -1)

transcript_map_data_list <- unique(c(transcript_map_data$alt_Symbol, transcript_map_data$SYMBOL))

#write(transcript_map_data_list, "data/transcript_gene_ontology_input_list.txt")

# filter protein data based on significance
protein_map_data <- weekes_data_host %>% 
  filter(log2_pm2_12_mock >= 1 | log2_pm2_12_mock <= -1 |
           log2_pm2_96_mock >= 1 | log2_pm2_96_mock <= -1 |
           log2_wcl2_12_mock >= 1 | log2_wcl2_12_mock <= -1 |
           log2_wcl2_96_mock >= 1 | log2_wcl2_96_mock <= -1 |
           del_log2_fc_pm2 >= 1 | del_log2_fc_pm2 <= -1 |
           del_log2_fc_wcl2 >= 1 | del_log2_fc_wcl2 <= -1)

protein_map_data_list <- unique(protein_map_data$Symbol)

#write(protein_map_data_list, "data/protein_gene_ontology_input_list.txt")
