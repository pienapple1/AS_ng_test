###Author Jinpei Lin
###Times 2024.07

library(spacexr)
library(Matrix)
library(Seurat)
library(stringr)
###based on single cell be a refence
# sc <- readRDS("/DATA/User/liqian7/pienapple/00AS_project/22RCTD/bin/sc_AS.RDS")
sc <- readRDS("./sc_mm_ann2_V2.RDS")
counts <- as.data.frame(sc@assays$RNA@counts)
sc@meta.data$barcode <- rownames(sc@meta.data)
meta_data  <- subset(sc@meta.data,select = c('barcode',"celltype_level1","nUMI"))
# create cell_types named list
cell_types <- meta_data$celltype_level1
names(cell_types) <- meta_data$barcode 
# convert to factor data type
cell_types <- as.factor(cell_types) 
# create nUMI named list
nUMI <- meta_data$nUMI
names(nUMI) <- meta_data$barcode 
reference <- Reference(counts, cell_types, nUMI)
str(reference)


## Examine reference object (optional)
#observe Digital Gene Expression matrix
print(dim(reference@counts)) 

###and then make a cycle to ST_data run RCTD


drds <- "/seurat_data/"
dout <- "/RCTD/"
for(s in Sys.glob(file.path(drds,"*_cellsegmentation_seurat.RDS"))){
    key <- gsub("_cellsegmentation_seurat.RDS","",basename(s))
    rds <- paste0(drds,key,"_cellsegmentation_seurat.RDS")
    dir.create(paste0(dout,key,"/"))
    setwd(paste0(dout,key,"/"))
    sp <- readRDS(rds)
    counts <- as.data.frame(sp@assays$RNA@counts)
    ##ST coordante
    #sp@meta.data$col <- sp@images$slice1@coordinates$col
    #sp@meta.data$row <- sp@images$slice1@coordinates$row
    coords <- subset(sp@meta.data,select = c("x","y"))
    # In this case, total counts per pixel is nUMI
    nUMI <- colSums(counts) 
    ### Create SpatialRNA object
    puck <- SpatialRNA(coords, counts, nUMI)
    str(puck)
    
    ## Examine SpatialRNA object (optional)
    print(dim(puck@counts))

    # histogram of log_2 nUMI
    hist(log(puck@nUMI,2))
    
    print(head(puck@coords)) # start of coordinate data.frame
    barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 
    
    # This list can be restricted if you want to crop the puck e.g. 
    # puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
    # on the plot:
    p <- plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')
    
    myRCTD <- create.RCTD(puck, reference, max_cores = 10,UMI_min = 50)
    myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
    str(myRCTD)
    results <- myRCTD@results
    
    # normalize the cell type proportions to sum to 1.
    norm_weights = normalize_weights(results$weights) 
    
    #list of cell type names
    cell_type_names <- myRCTD@cell_type_info$info[[2]] 
    
    spatialRNA <- myRCTD@spatialRNA
    
    ## you may change this to a more accessible directory on your computer.
    resultsdir <- 'RCTD_Plots' 
    dir.create(resultsdir)
    
    # make the plots 
    plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
    plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
    

    plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

    plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
    
    saveRDS(results,paste0(key,"_RCTD_result_list.RDS"))
    
}