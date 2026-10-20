library(spacexr)
library(Matrix)
library(Seurat)

###based on single cell be a refence
sc <- readRDS("./sc_mm_ann2_V2.RDS")
counts <- as.data.frame(sc@assays$SCT@counts)
sc@meta.data$barcode <- rownames(sc@meta.data)
meta_data  <- subset(sc@meta.data,select = c('barcode',"celltype_level2","nUMI"))
# create cell_types named list
cell_types <- meta_data$celltype_level2
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

drds <- "/DATA/User/liqian7/pienapple/10collating_ASdata/04fndata_RDS/"
dout <- "/DATA/User/liqian7/pienapple/17other_test/02rctd_V2/"

for(s in Sys.glob(file.path(drds,"*.RDS"))){
    key <- gsub(".RDS","",basename(s))
    rds <- paste0(drds,key,".RDS")
    dir.create(paste0(dout,key,"/"), showWarnings = FALSE)
    setwd(paste0(dout,key,"/"))
    sp <- readRDS(rds)
    counts <- as.data.frame(sp@assays$RNA@counts)
    coords <- subset(sp@meta.data,select = c("x","y"))
    nUMI <- colSums(counts) 
    ### Create SpatialRNA object
    puck <- SpatialRNA(coords, counts, nUMI)
    str(puck)
    
    ## Examine SpatialRNA object (optional)
    print(dim(puck@counts))

    # histogram of log_2 nUMI
    hist(log(puck@nUMI,2))
    
    print(head(puck@coords))
    barcodes <- colnames(puck@counts)
    
    p <- plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')
    
    # UMI_min_sigma = 50
    cat("Processing", key, ": Trying UMI_min_sigma = 50\n")
    
    tryCatch({

        myRCTD <- create.RCTD(puck, reference, max_cores = 10, UMI_min = 50)
        myRCTD@config$UMI_min_sigma <- 50
        

        myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
        
        cat("Success with UMI_min_sigma = 50 for", key, "\n")
        

        sigma_used <- 50
        
    }, error = function(e) {

        cat("Error with UMI_min_sigma = 50 for", key, ":", conditionMessage(e), "\n")
        cat("Retrying with UMI_min_sigma = 100\n")

        myRCTD <<- create.RCTD(puck, reference, max_cores = 10, UMI_min = 50)
        myRCTD@config$UMI_min_sigma <<- 100
        

        myRCTD <<- run.RCTD(myRCTD, doublet_mode = 'doublet')
        
        cat("Success with UMI_min_sigma = 100 for", key, "\n")
        

        sigma_used <<- 100
    })
    

    cat("For sample", key, "used UMI_min_sigma =", sigma_used, "\n")
    
    str(myRCTD)
    results <- myRCTD@results
    
    # normalize the cell type proportions to sum to 1.
    norm_weights = normalize_weights(results$weights) 
    
    #list of cell type names
    cell_type_names <- myRCTD@cell_type_info$info[[2]] 
    
    spatialRNA <- myRCTD@spatialRNA
    
    ## you may change this to a more accessible directory on your computer.
    resultsdir <- 'RCTD_Plots' 
    dir.create(resultsdir, showWarnings = FALSE)
    
    # make the plots 
    plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
    plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
    plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
    plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
    

    saveRDS(results, paste0(key, "_RCTD_result_list_sigma", sigma_used, ".RDS"))
    
    ##merge_RCTD dataframe and sp@meta.data to RAW ST RDS
    rowname <- rownames(results$results_df)
    sp@meta.data$rown <- rownames(sp@meta.data)
    sp1 <- subset(sp, subset = rown %in% rowname)
    
    sp1@meta.data <- merge(sp1@meta.data, results$results_df, by = 0, all = TRUE)
    rownames(sp1@meta.data) <- sp1@meta.data$Row.names
    sp1@meta.data <- sp1@meta.data[,-1]
    

    saveRDS(sp1, paste0(key, "_merge_RCTD_result_sigma", sigma_used, ".RDS"))
    

    write.table(
        data.frame(sample = key, UMI_min_sigma_used = sigma_used),
        file = paste0(dout, "sigma_usage_log.txt"),
        append = TRUE,
        quote = FALSE,
        row.names = FALSE,
        col.names = !file.exists(paste0(dout, "sigma_usage_log.txt"))
    )
}