source("import_pellett_sam_probe.R")

library(limma)
library(gplots)
library(dendextend)

# import data sets needed
secretory_genes_of_interest <- read.table("secretory_genes_of_interest.txt", header = T, sep = "\t")
gene_ontology_secretory_genes <- read.table("transcript_filtered_gene_list.txt", header = F)

# adjust dataframes to be characters for joining functions later
secretory_genes_of_interest$PROBE_ID <- as.character(secretory_genes_of_interest$PROBE_ID)
pellett_probe$PROBE_ID <- as.character(pellett_probe$PROBE_ID)

secretory_genes_of_interest$SYMBOL <- as.character(secretory_genes_of_interest$SYMBOL)
gene_ontology_secretory_genes$V1 <- as.character(gene_ontology_secretory_genes$V1)
pellett_probe$SYMBOL <- as.character(pellett_probe$SYMBOL)

# finding the gene ontology related genes not represented in the initial curation of vesicle-mediated
# transport related transcripts
genes_not_matched <- anti_join(gene_ontology_secretory_genes, secretory_genes_of_interest, by = c("V1" = "SYMBOL"))

# accessing raw data set to attach data to initial genes of interest and missing genes
secretory_genes_of_interest_data <- left_join(secretory_genes_of_interest, sam_probe_data, by = c("PROBE_ID", "SYMBOL"))
genes_not_matched_data <- sam_probe_data %>% 
  filter(SYMBOL %in% genes_not_matched$V1)

# combining all of the genes of interest with their data
heatmap_genes_of_interest <- full_join(secretory_genes_of_interest_data, genes_not_matched_data, rm.na = T)

# removing columns that aren't being used for heatmap
col_to_remove <- c("\\w3\\..*", "\\w4\\..*", "\\w5\\..*", "\\w7\\..*", "\\w9\\..*", "\\w10\\..*", "\\w12\\..*")
remove_col <- grep(paste(col_to_remove, collapse = "|"), colnames(heatmap_genes_of_interest))
heatmap_genes_of_interest[colnames(heatmap_genes_of_interest)[c(remove_col)]] <- NULL

# removing rows without any detected data
heatmap_genes_of_interest <- heatmap_genes_of_interest[!is.na(heatmap_genes_of_interest$A1.AVG_Signal), ]

# changing colnames to be compatible for assignment when importing and generating heatmap
colnames(heatmap_genes_of_interest) <- gsub("A1\\.", "1\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("A2\\.", "2\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("A6\\.", "3\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("A8\\.", "4\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("A11\\.", "5\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("B1\\.", "6\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("B2\\.", "7\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("B6\\.", "8\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("B8\\.", "9\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("B11\\.", "10\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("C1\\.", "11\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("C2\\.", "12\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("C6\\.", "13\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("C8\\.", "14\\.", colnames(heatmap_genes_of_interest))
colnames(heatmap_genes_of_interest) <- gsub("C11\\.", "15\\.", colnames(heatmap_genes_of_interest))

heatmap_genes_of_interest <- heatmap_genes_of_interest %>% 
  select(-alt_Symbol)

# exporting as a table for future use
write_tsv(heatmap_genes_of_interest, "transcript_heatmap_genes_of_interest.txt")



# setting up the data for plotting ----------------------------------------

targets <- readTargets()
eset <- read.ilmn(files = "transcript_heatmap_genes_of_interest.txt", probeid = "SYMBOL")

design <- cbind("Mock 1"=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                "Mock 2"=c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
                "Mock 3"=c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
                "Towne.12hpi 1"=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                "Towne.12hpi 2"=c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
                "Towne.12hpi 3"=c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
                "Towne.96hpi 1"=c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
                "Towne.96hpi 2"=c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
                "Towne.96hpi 3"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
                "AD169.12hpi 1"=c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                "AD169.12hpi 2"=c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
                "AD169.12hpi 3"=c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
                "AD169.96hpi 1"=c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
                "AD169.96hpi 2"=c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
                "AD169.96hpi 3"=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1))

fit <- lmFit(log2(eset$E),design)

esetSel <- as.matrix(fit)

# calculates distance using pearson's correlation coefficient
dist.pear <- function (x) as.dist(1-cor(t(x)))

# uses average linkage clustering
hclust.ave <- function(x) hclust(x, method = "average")

# sets 
lmat = rbind(4:3,2:1)

# size of output plot
dev.new(width=25, height=30)

lwid = c(6,12.5)

lhei = c(4,21)

# make the first version of the graph to generate a col dendrogram
transcript_secrete <- heatmap.2(esetSel, col=redgreen(100), scale="row", key=TRUE,
                                symkey=FALSE, density.info="none", trace="none",
                                cexCol=1.0, margins=c(6,1), srtCol=45,
                                distfun=dist.pear, hclustfun=hclust.ave, lmat=lmat,
                                lhei=lhei, lwid=lwid)

# changing the order of the col dend
col_dend_rotate <- rotate(transcript_secrete$colDendrogram, c(7:15,1:6))

# now for final plotting
#pdf(file="transcript_heatmap.pdf")

heatmap.2(esetSel, col=redgreen(100), scale="row", key=TRUE,
          symkey=FALSE, density.info="none", trace="none",
          cexCol=1.0, margins=c(6,1), srtCol=45,
          distfun=dist.pear, hclustfun=hclust.ave, lmat=lmat,
          lhei=lhei, lwid=lwid, Colv = col_dend_rotate)

#dev.off()


# plotting row mean log2 information --------------------------------------

# creating matrix form of heatmap plot
transcript_secrete_heatmap_matrix <- esetSel[rev(transcript_secrete$rowInd), c(1:6, 10:12, 7:9, 13:15)]

# pulling out mean log2, adding col of symbol names (creating unique names if duplicated),
# and creating col for ordering as desired
transcript_heatmap_row_means <- data_frame(rowMeans = rowMeans(transcript_secrete_heatmap_matrix, na.rm = T)) %>% 
  add_column(Symbol = make.names(rownames(transcript_secrete_heatmap_matrix), unique = T), .before = 1) %>% 
  add_column(Order = 1:nrow(transcript_secrete_heatmap_matrix), .before = 1)

# plotting the row means
ggplot(transcript_heatmap_row_means) +
  geom_bar(mapping = aes(x = reorder(Symbol, rev(Order)), y = rowMeans), stat = "identity", fill = "grey50") +
  geom_hline(yintercept = c(5,10,15), color = "black", size = .25) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5), expand = c(0,0), position = "right") +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        text = element_text(face = "bold", family = "serif", size = 12, color = "black")) +
  coord_flip()

# save the plots
#ggsave("transcript_heatmap_row_means.png", dpi = 600)
#ggsave("transcript_heatmap_row_means.pdf", dpi = 600)

