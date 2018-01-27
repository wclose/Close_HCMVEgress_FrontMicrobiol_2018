source("heatmap_probe.R")

library(limma)
library(gplots)
library(dendextend)

targets <- readTargets()
eset <- read.ilmn(files = "heatmap_genes_of_interest.txt", probeid = "SYMBOL")

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
secrete <- heatmap.2(esetSel, col=redgreen(100), scale="row", key=TRUE,
                     symkey=FALSE, density.info="none", trace="none",
                     cexCol=1.0, margins=c(6,1), srtCol=45,
                     distfun=dist.pear, hclustfun=hclust.ave, lmat=lmat,
                     lhei=lhei, lwid=lwid)

# changing the order of the col dend
col_dend_rotate <- rotate(secrete$colDendrogram, c(7:15,1:6))

# now for final plotting

png(file="heatmap2.png", width = 2718, height = 4000, res = 300)

heatmap.2(esetSel, col=redgreen(100), scale="row", key=TRUE,
                     symkey=FALSE, density.info="none", trace="none",
                     cexCol=1.0, margins=c(6,1), srtCol=45,
                     distfun=dist.pear, hclustfun=hclust.ave, lmat=lmat,
                     lhei=lhei, lwid=lwid, Colv = col_dend_rotate)

dev.off()

summary(secrete)
