library(vcd)
library(stringr)
m <- readRDS("2g.RDS")

order <- c('W20_N','W0','W2','W6','W10','W20')
m$week <- factor(m$week,levels = order)
order1 <- c('GC-AA','GC-AR','GC-DA','LC-AA','LC-AR','LC-DA')
m$AD <- factor(m$AD,levels = order1)
mosaic(~week + AD+satge_p_all, data = m)