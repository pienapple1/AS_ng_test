###R4
library(Seurat)
library(stringr)
library(ggplot2)
sp <- readRDS("W10.RDS")
sp <- NormalizeData(sp,normalization.method = "LogNormalize")
gene <- read.csv("plaque_gene.csv")
gene1 <- gene$plaque_gene_mm
sp1 <- AddModuleScore(sp, features = list(gene1))
colnames(sp1@meta.data) <- str_replace(colnames(sp1@meta.data),"Cluster1","plaque_score")
ggplot(sp1@meta.data, aes(y=y,x=x,colour=plaque_score)) + geom_point(size=0.01) +labs(x="",y="")+scale_colour_gradient(low ="#54BAD3",high = "#ce1256",limits=c(0,2))+theme_bw()