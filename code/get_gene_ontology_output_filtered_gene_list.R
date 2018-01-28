library(tidyverse)

# start with running get_gene_ontology_input_list and use the output in cytoscape to generate a lists of terms
# annotated by vesicle-mediated transport, secretion, and intracellular transport

# importing transcript data
transcript_vesicle_mediated_transport_16192 <- read_tsv("data/transcript_vesicle_mediated_transport_16192.txt", col_names = F)
transcript_secretion_46903 <- read_tsv("data/transcript_secretion_46903.txt", col_names = F)
transcript_intracellular_transport_46907 <- read_tsv("data/transcript_intracellular_transport_46907.txt", col_names = F)

transcript_filtered_gene_list <- full_join(full_join(transcript_vesicle_mediated_transport_16192,
                                                     transcript_secretion_46903, by = "X1"),
                                           transcript_intracellular_transport_46907, by = "X1")

#write(transcript_filtered_gene_list$X1, "data/transcript_filtered_gene_list.txt")

# importing transcript data
protein_vesicle_mediated_transport_16192 <- read_tsv("data/protein_vesicle_mediated_transport_16192.txt", col_names = F)
protein_secretion_46903 <- read_tsv("data/protein_secretion_46903.txt", col_names = F)
protein_intracellular_transport_46907 <- read_tsv("data/protein_intracellular_transport_46907.txt", col_names = F)

protein_filtered_gene_list <- full_join(full_join(protein_vesicle_mediated_transport_16192,
                                                     protein_secretion_46903, by = "X1"),
                                           protein_intracellular_transport_46907, by = "X1")

#write(protein_filtered_gene_list$X1, "data/protein_filtered_gene_list.txt")

