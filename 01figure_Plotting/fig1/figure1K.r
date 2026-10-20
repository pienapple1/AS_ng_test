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
sp <- readRDS("Q40018K2_1_fn.RDS")
sp <- NormalizeData(sp)
sp <- AddModuleScore(sp,gss[["SMC_like"]])
colnames(sp@meta.data) <- str_replace(colnames(sp@meta.data),"Cluster1","SMC-like")####The same applies to other genesets.
SpatialFeaturePlot(sp,features = "SMC-like",pt.size.factor = 0.6,max.cutoff = 1.2)+ggplot2::scale_fill_continuous(low ="#54BAD3",high = "#ce1256")+
coord_cartesian(xlim = c(0, 22000), ylim = c(0,22000))