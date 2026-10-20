library(hdf5r)
library(SeuratDisk)
library(Seurat)
library(ggplot2)

drds <- "/00data/02cellsegmentation/"
dout <- "/00data/04seurat_data/"

name <- list.files(path =drds,full.names = FALSE)
for(chip in name){
    adata_data.loom <-  Connect(filename = paste0(drds,chip,"/",chip,"_data.loom"), mode = 'r')
    metadata <- read.csv(paste0(drds,chip,"/",chip,"_metadata.csv"), row.names = 1)
    #read matrix 
    full.matrix <- adata_data.loom[["matrix"]][,]
    ###read gene
    gene <- adata_data.loom[["row_attrs"]][["Gene"]][]
    #check cell
    barcode <- adata_data.loom[["col_attrs"]][["CellID"]][]
    #rownames and colnames
    rownames(full.matrix) <- barcode
    colnames(full.matrix) <- gene
    #Always remember to close loom files when done
    adata_data.loom$close_all()
    # coor flip
    full.matrix <- t(full.matrix)
    #Create seurat file
    object  <- CreateSeuratObject(full.matrix, meta.data = metadata)
    colnames(object@meta.data)
    p <- ggplot(object@meta.data,aes(x=x,y =y,color=orig.ident))+geom_point()+theme_classic ()
    print(p)
    saveRDS(object,paste0(dout,chip,"_cellsegmentation_seurat.RDS"))
}