###3o ECM deg
library(scRNAtoolVis)
m <- readRDS("3o.RDS")###deg base FindAllmarkers()
jjVolcano(diffData = m,
          tile.col = c("#DCFFEB","#88C7ED","#ED7A90"),
          size  = 3.5,topGeneN = 5,
          fontface = 'italic',pSize=2.5,legend.position=c(0.9,0.97))