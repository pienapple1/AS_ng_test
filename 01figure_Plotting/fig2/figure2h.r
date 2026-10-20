####2h
###R4

#####sp==merge all sample
library(ggplot2)
library(ggridges)
library(RColorBrewer)
library(scales)

cols <- c("Macrophage" = '#DF65B0', "Modulated_SMC" = "#A1CFFA","SMC" = '#1F78B4', "EC" = '#238B45', "Fibroblast" = "#8B4A4B")

sp@meta.data$week <- factor(sp@meta.data$week, levels = c("W0", "W2", "W6", "W10", "W20", "W20_N"))
order <- c("Fibroblast", "SMC", "Modulated_SMC", "Macrophage", "EC")
sp@meta.data$first_type <- factor(sp@meta.data$first_type, levels = order)

####plot
ggplot(drds@meta.data, aes(x = layer1, y = first_type, fill = first_type, color = first_type)) +
geom_density_ridges() +facet_wrap(~ week, scales = "free_x", ncol = length(unique(drds@meta.data$week)))