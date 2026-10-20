library(Seurat)
library(GSVA)
library(limma)
library(pheatmap)
library(ggplot2)
library(stringr)
library(Seurat)
library(stringr)
library(Seurat)
library(stringr)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggplotify)
library(ggpubr)

drds <- "/DATA/User/liqian7/pienapple/24AS_fig_add/00DATA/05mm_data/" ####seurat sample data (finsh base all ST deal Data processing workflow)
i=0
for(s in Sys.glob(file.path(drds,"*_mm.RDS"))){
    i=i+1
    key <- gsub("_mm.RDS","",basename(s))
    sp <- readRDS(paste0(drds,key,"_mm.RDS"))
    sp@meta.data$Celltype <- factor(sp@meta.data$Celltype,levels = c('mSMCss','Macro.','FBs','EC','SMCs',"T cells"))
    n=0
    for(k in unique(sp@meta.data$cd)){
        n=n+1
        pbmc <- subset(sp,subset = cd==k)###cd is sample * six area
        week <- unique(pbmc@meta.data$Week)
        sample <- unique(pbmc@meta.data$sample)
        counts <- as.data.frame(rowSums(as.matrix(pbmc@assays$RNA@counts)))
        colnames(counts) <- k
        ct <- as.data.frame(prop.table(table(pbmc$Celltype)))
        colnames(ct) <- c("celltype","percent")
        ct$percent <- as.numeric(ct$percent)
        ct <- as.data.frame(t(ct))
        colnames(ct) <- ct[1,]
        ct <- ct[-1,]
        rownames(ct) <- k
        ct$Mac_Msmc_SMC <- paste0(as.numeric(ct$Macro.)+as.numeric(ct$mSMCss)-as.numeric(ct$SMCs))
        ct$sample <- sample
        ct$week <- week
        if(n==1){
            m <- counts
            ct1 <- ct
        }
        else{
            m <- merge(m,counts,all=TRUE,by=0)
            rownames(m) <- m$Row.names
            m <- m[,-1]
            ct1 <- rbind(ct1,ct)
        }
    }
    if(i==1){
        mer <- m
        meta <- ct1
    }
    else{
        mer <- merge(mer,m,all=TRUE,by=0)
        rownames(mer) <- mer$Row.names
        mer <- mer[,-1]
        meta <- rbind(meta,ct1)
    }
}

mer[is.na(mer)] <- 0

bulk <- CreateSeuratObject(counts = mer,meta.data = meta)###Create area pseudobulk seurat data

bulk <- NormalizeData(bulk,normalization.method = "LogNormalize",scale.factor = 400000)

expr <- bulk@assays$RNA@data
gene <- read.csv("1.7hotspot_plaque_related_gene.csv")
g1 <- list(unique(gene$plaque_gene_mm))
names(g1) <- "plaque_geneset"

gsva.res1 <- gsva(
    as.matrix(expr),
    as.list(g1),
    method="gsva",
    kcdf = "Gaussian",
    mx.diff = F,
    verbose = F)

bulk@meta.data <- merge(bulk@meta.data,t(gsva.res1),all = TRUE,by = 0)
colnames(bulk@meta.data) <- str_replace(colnames(bulk@meta.data),"plaque_geneset","plaque_score")

df <- bulk@meta.data
df$stage_ud <- "NA"
df$Mac_Msmc_SMC <- as.numeric(df$Mac_Msmc_SMC)
df[which(df$plaque_score <= -0.5&df$Mac_Msmc_SMC<0),]$stage_ud <- "stage1"
df[which(df$plaque_score<=0&df$plaque_score> -0.5&df$Mac_Msmc_SMC> -0.5&df$Mac_Msmc_SMC<0),]$stage_ud <- "stage2"
df[which(df$plaque_score<= 0.5&df$plaque_score>0&df$Mac_Msmc_SMC> -0.5&df$Mac_Msmc_SMC<0.5),]$stage_ud <- "stage3"
df[which(df$plaque_score>0.5&df$Mac_Msmc_SMC> 0),]$stage_ud <- "stage4"

df <- df[which(df$stage_ud!="NA"&df$sample!="filter"),]
saveRDS(df,"1.7result_data.RDS")





