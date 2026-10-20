library(Seurat)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(stringr)
gss <- list(SMC_like=c("MYH11","ACTA2","CNN1","TAGLN","LMOD1","PLN","MYOCD","SRF"),
     Macrophage_like=c("LGALS3","CD68","LAMP2","CD36","SPP1","CCL2","KLF4","IRF8"),
     Chondro_like=c("RUNX2","ALPL","SPP1","SOX9","ACAN","BMP2"),
     Fibromyocyte=c("TCF21","FN1","COL1A1","LUM","DCN","PDGFRB","AHR")
    )
sc <- readRDS("sc.RDS")
sc <- NormalizeData(sc)
sc <- AddModuleScore(sc,gss[["SMC_like"]])
colnames(sc@meta.data) <- str_replace(colnames(sc@meta.data),"Cluster1","SMC-like")####The same applies to other genesets.
FeaturePlot(sc,features = "SMC-like")