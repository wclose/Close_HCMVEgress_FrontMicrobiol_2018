source("code/import_weekes_protein.R")

library(dendextend)
library(gplots)

# Importing data ----------------------------------------------------------

# importing list of genes of interest
# genes from literature
secretory_genes_of_interest <- read.table("data/secretory_genes_of_interest.txt", header = T, sep = "\t")
secretory_genes_of_interest$PROBE_ID <- as.character(secretory_genes_of_interest$PROBE_ID)
secretory_genes_of_interest$SYMBOL <- as.character(secretory_genes_of_interest$SYMBOL)

#genes from gene ontology network
gene_ontology_secretory_genes <- read.table("data/protein_filtered_gene_list.txt", header = F)
gene_ontology_secretory_genes$V1 <- as.character(gene_ontology_secretory_genes$V1)

# combining lists
genes_not_matched <- anti_join(gene_ontology_secretory_genes, secretory_genes_of_interest, by = c("V1" = "SYMBOL"))
heatmap_genes_of_interest <- c(secretory_genes_of_interest$SYMBOL, genes_not_matched$V1)



# Curating weekes data ----------------------------------------------------

# filtering weekes data using each column of gene symbols to make sure protein aliases are detected
weekes_heatmap_data <- filter(weekes_data_host, Symbol %in% heatmap_genes_of_interest)

# removing any proteins that don't have data in either pm or wcl columns (removes 8 rows)
weekes_heatmap_data <- weekes_heatmap_data[!is.na(weekes_heatmap_data$pm2_Mock_1) | !is.na(weekes_heatmap_data$wcl2_Mock_1), ]

# pulling out the columns for the heatmap
weekes_heatmap_data_col <- weekes_heatmap_data[,c("Symbol", "pm2_Mock_1", "pm2_Mock_2", "pm2_12h", "pm2_96h",
                                                  "wcl2_Mock_1", "wcl2_Mock_2", "wcl2_12h", "wcl2_96h")]



# Creating matrix for plotting --------------------------------------------

# creating a naming vector for rows of the matrix
rnames <- weekes_heatmap_data_col[,1]

# turning the curated weekes data into a matrix
mat_data <- data.matrix(weekes_heatmap_data_col[,2:ncol(weekes_heatmap_data_col)])

# naming the matrix rows
rownames(mat_data) <- rnames$Symbol

# log2 transorming the data
mat_data_log2 <- log(mat_data, 2)

# replacing any NA/Inf values with 0 for generation of clustering dendrogram
# (if this isn't done, the distance matrix/clustering errors and can't generate a plot)
mat_data_log2_fin <- ifelse(!is.finite(mat_data_log2), 0, mat_data_log2)
length(mat_data_log2_fin[,1])


# Plotting the data -------------------------------------------------------

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

# make the first version of the graph to generate row and col dendrograms using matrix with non-finite
# values as 0 for calculation of distance matrices and clustering
protein_secrete <- heatmap.2(mat_data_log2_fin, col=redgreen(100), scale="row", key=TRUE,
                                symkey=FALSE, density.info="none", trace="none",
                                cexCol=1.0, margins=c(6,1), srtCol=45,
                                distfun=dist.pear, hclustfun=hclust.ave, lmat=lmat,
                                lhei=lhei, lwid=lwid)

# changing the order of the col dend to put in desired order
protein_col_dend_rotate <- rotate(protein_secrete$colDendrogram, c(8,7,5,6,3,4,1,2))

# now for final plotting
#pdf(file="figures/protein_heatmap.pdf")

heatmap.2(mat_data_log2, col=redgreen(100), scale="row", key=TRUE,
          symkey=FALSE, density.info="none", trace="none",
          cexCol=1.0, margins=c(6,1), srtCol=45, lmat=lmat,
          lhei=lhei, lwid=lwid, Colv = protein_col_dend_rotate,
          Rowv = protein_secrete$rowDendrogram, na.color = "white")

#dev.off()



# plotting row mean log2 information --------------------------------------

# creating matrix form of heatmap plot
protein_secrete_heatmap_matrix <- mat_data_log2[rev(protein_secrete$rowInd), c(5,6,7,8,1,2,3,4)]

# pulling out mean log2, adding col of symbol names (creating unique names if duplicated),
# and creating col for ordering as desired
protein_heatmap_row_means <- data_frame(rowMeans = rowMeans(protein_secrete_heatmap_matrix, na.rm = T)) %>% 
  add_column(Symbol = make.names(rownames(protein_secrete_heatmap_matrix), unique = T), .before = 1) %>% 
  add_column(Order = 1:nrow(protein_secrete_heatmap_matrix), .before = 1)

# plotting the row means
ggplot(protein_heatmap_row_means) +
  geom_bar(mapping = aes(x = reorder(Symbol, rev(Order)), y = rowMeans), stat = "identity", fill = "grey50") +
  geom_hline(yintercept = c(5,10,15), color = "black", size = .25) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5), expand = c(0,0), position = "right") +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        text = element_text(face = "bold", family = "serif", size = 12, color = "black")) +
  coord_flip()

# save the plots
#ggsave("figures/protein_heatmap_rows_means.png", dpi = 600)
#ggsave("figures/protein_heatmap_rows_means.pdf", dpi = 600)



###########################################################################################
# this section was used as a test to generate a dendrogram based solely on wcl data

# pulling out only the wcl data
#mat_data_wcl <- data.matrix(weekes_heatmap_data_col[,6:ncol(weekes_heatmap_data_col)])

# assigning row names
#rownames(mat_data_wcl) <- rnames$Symbol

# log2 transforming the data
#mat_data_wcl_log2 <- log(mat_data_wcl, 2)

# had to include random number generation to avoid error due to lack of variance if all 0's were used
# to replace the Na/Inf values
#mat_data_wcl_log2_fin <- ifelse(!is.finite(mat_data_wcl_log2), rnorm(10, mean = 0.03, sd = 0.01), mat_data_wcl_log2)

# size of output plot
#dev.new(width=25, height=30)

# sets size and position of matrix
#lmat = rbind(4:3,2:1)

#lwid = c(6,12.5)

#lhei = c(4,21)

# creating a dendrogram based on clustering of the wcl samples alone being that there's more data
#protein_secrete_wcl <- heatmap.2(mat_data_wcl_log2_fin, col=redgreen(100), scale="row", key=TRUE,
                                 #symkey=FALSE, density.info="none", trace="none",
                                 #cexCol=1.0, margins=c(6,1), srtCol=45,
                                 #distfun=dist.pear, hclustfun=hclust.ave, lmat=lmat,
                                 #lhei=lhei, lwid=lwid)

# changing the order of the col dend to put in desired order
#protein_col_dend_rotate <- rotate(protein_secrete$colDendrogram, c(8,7,5,6,3,4,1,2))

# now for final plotting
#png(file="figures/protein_heatmap_v2.png", width = 2718, height = 4000, res = 300)

#heatmap.2(mat_data_log2, col=redgreen(100), scale="row", key=TRUE,
          #symkey=FALSE, density.info="none", trace="none",
          #cexCol=1.0, margins=c(6,1), srtCol=45, lmat=lmat,
          #lhei=lhei, lwid=lwid, Colv = protein_col_dend_rotate,
          #Rowv = protein_secrete_wcl$rowDendrogram, na.color = "gray35")

#dev.off()

###########################################################################################


