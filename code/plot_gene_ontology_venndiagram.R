library(VennDiagram)
library(tidyverse)

# importing transcript data
transcript_transport_6810 <- read_tsv("data/transcript_transport_6810.txt", col_names = F)
transcript_vesicle_mediated_transport_16192 <- read_tsv("data/transcript_vesicle_mediated_transport_16192.txt", col_names = F)
transcript_secretion_46903 <- read_tsv("data/transcript_secretion_46903.txt", col_names = F)
transcript_intracellular_transport_46907 <- read_tsv("data/transcript_intracellular_transport_46907.txt", col_names = F)
transcript_neutrophil_degranulation_43312 <- read_tsv("data/transcript_neutrophil_degranulation_43312.txt", col_names = F)

# importing protein data
protein_transport_6810 <- read_tsv("data/protein_transport_6810.txt", col_names = F)
protein_vesicle_mediated_transport_16192 <- read_tsv("data/protein_vesicle_mediated_transport_16192.txt", col_names = F)
protein_secretion_46903 <- read_tsv("data/protein_secretion_46903.txt", col_names = F)
protein_intracellular_transport_46907 <- read_tsv("data/protein_intracellular_transport_46907.txt", col_names = F)
protein_neutrophil_degranulation_43312 <- read_tsv("data/protein_neutrophil_degranulation_43312.txt", col_names = F)


# making 'transport' venn diagram 
#pdf(file="figures/vp_transport.pdf")
vp_transport <- venn.diagram(list(t=transcript_transport_6810$X1,p=protein_transport_6810$X1), 
                   fill = c("firebrick3", "cyan4"), alpha = 0.3, filename = NULL);
grid.draw(vp_transport)
#dev.off()


# vesicle mediated transport diagram
#pdf(file="figures/vp_vesicle_mediated_transport.pdf")
vp_vesicle_mediated_transport <- venn.diagram(list(t=transcript_vesicle_mediated_transport_16192$X1,p=protein_vesicle_mediated_transport_16192$X1),
                   fill = c("firebrick3", "cyan4"), alpha = 0.3, filename = NULL);
grid.draw(vp_vesicle_mediated_transport)
#dev.off()


# secretion diagram
#pdf(file="figures/vp_secretion.pdf")
vp_secretion <- venn.diagram(list(t=transcript_secretion_46903$X1, p=protein_secretion_46903$X1),
                   fill = c("firebrick3", "cyan4"), alpha = 0.3, filename = NULL);
grid.draw(vp_secretion)
#dev.off()


# intracellular transport diagram
#pdf(file="figures/vp_intracellular_transport.pdf")
vp_intracellular_transport <- venn.diagram(list(t=transcript_intracellular_transport_46907$X1, p=protein_intracellular_transport_46907$X1), 
                   fill = c("firebrick3", "cyan4"), alpha = 0.3, filename = NULL);
grid.draw(vp_intracellular_transport)
#dev.off()


# neutrophil degranulation diagram
#pdf(file="figures/vp_neutrophil_degranulation.pdf")
vp_neutrophil_degranulation <- venn.diagram(list(t=transcript_neutrophil_degranulation_43312$X1, p=protein_neutrophil_degranulation_43312$X1), 
                                           fill = c("firebrick3", "cyan4"), alpha = 0.3, filename = NULL);
grid.draw(vp_neutrophil_degranulation)
#dev.off()
