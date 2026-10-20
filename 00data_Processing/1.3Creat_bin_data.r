library(raster)
library(Seurat)
dgz <- "/00gem/"
dout <- "/02RDS_raw/"
for(s in Sys.glob(file.path(dgz,"*.gem.gz"))){
    name <- gsub(".gem.gz","",basename(s))
    input <- paste0(dgz,name,".gem.gz")
    binsize <- 50
    output <- dout
    source(file.path("1.3Create_bin_data.R"))
    seurat_spatialObj <- LoadBGI_Spatial(input, outdir = output, bin_data = T, bin_size = binsize, cell_mask = F, area_mask = F, pro_name = name)
}