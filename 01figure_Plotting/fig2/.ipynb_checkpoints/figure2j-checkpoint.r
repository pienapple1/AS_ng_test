library(Seurat)
library(pheatmap)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggplotify)
library(ggpubr)
options(warn=-1)

####just select layer*0-25*
sp <- readRDS("2I_raw.RDS")
sp@meta.data$stage_cell <- factor(sp@meta.data$stage_cell,levels=c('stage4_adventitia','stage3_adventitia','stage2_adventitia','stage1_adventitia',
                                                                   "stage1_intima","stage2_intima","stage3_intima","stage4_intima"))
Idents(sp) <- as.factor(sp@meta.data$stage_cell)

####dotplot 
p <- DotPlot(object = sp4, cols =c("#89ABE3FF","#F55665FF"),scale.max = 40,scale.min = 0,
        features=c("Acta2","Hp","Treml4","Ly6c2","Nes","Lgals3","Mmp12","Acp5","Gngt2","Cd9","F13a1",
                   "Lyve1","Folr2","Pf4","Mrc1","C4b","Maf","Gm26917","Fcgr3","Ctsd","Gpnmb","Spp1","Fn1"))+theme(axis.text.x = element_text(face="italic", hjust=1,angle = 90), axis.text.y = element_text(face="bold")) +theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="")
p

# compute the gene set score base AddModuleScore
gs <- readRDS("2I_geneset.RDS")

obj_merged = AddModuleScore(sp, gs)
for (i in 1:ntype) {
  obj_merged@meta.data[[subtype[i]]] = obj_merged@meta.data[[glue('Cluster{i}')]]
}

