library(VennDiagram)
library(tidyverse)

transcript_transport_6810 <- read_tsv("transcript_transport_6810", col_names = F)
transcript_vesicle_mediated_transport_16192 <- read_tsv("transcript_vesicle_mediated_transport_16192", col_names = F)
transcript_secretion_46903 <- read_tsv("transcript_secretion_46903", col_names = F)
transcript_intracellular_transport_46907 <- read_tsv("transcript_intracellular_transport_46907", col_names = F)


protein_transport_6810 <- read_tsv("protein_transport_6810", col_names = F)
protein_vesicle_mediated_transport_16192 <- read_tsv("protein_vesicle_mediated_transport_16192", col_names = F)
protein_secretion_46903 <- read_tsv("protein_secretion_46903", col_names = F)
protein_intracellular_transport_46907 <- read_tsv("protein_intracellular_transport_46907", col_names = F)



pdf(file="vp_transport.pdf")
vp_transport <- venn.diagram(list(t=transcript_transport_6810$X1,p=protein_transport_6810$X1), 
                   fill = c("firebrick3", "cyan4"), alpha = 0.3, filename = NULL);
grid.draw(vp_transport)
dev.off()



pdf(file="vp_vesicle_mediated_transport.pdf")
vp_vesicle_mediated_transport <- venn.diagram(list(t=transcript_vesicle_mediated_transport_16192$X1,p=protein_vesicle_mediated_transport_16192$X1),
                   fill = c("firebrick3", "cyan4"), alpha = 0.3, filename = NULL);
grid.draw(vp_vesicle_mediated_transport)
dev.off()



pdf(file="vp_secretion.pdf")
vp_secretion <- venn.diagram(list(t=transcript_secretion_46903$X1, p=protein_secretion_46903$X1),
                   fill = c("firebrick3", "cyan4"), alpha = 0.3, filename = NULL);
grid.draw(vp_secretion)
dev.off()



pdf(file="vp_intracellular_transport.pdf")
vp_intracellular_transport <- venn.diagram(list(t=transcript_intracellular_transport_46907$X1, p=protein_intracellular_transport_46907$X1), 
                   fill = c("firebrick3", "cyan4"), alpha = 0.3, filename = NULL);
grid.draw(vp_intracellular_transport)
dev.off()
