library(affy) 
library(hgu133plus2.db) 
library(hgu133plus2cdf)
library(limma) 
library(readr)

cel_file <- "your/file/path" 
data <- ReadAffy(filenames = cel_file) 
data_rma <- rma(data) 
exprs_data_rma <- exprs(data_rma) 
head(exprs_data_rma)

#Convert the probe ID to the gene ID
process_and_save_cel_file <- function(file, output_dir) {
  data <- ReadAffy(filenames = file)
  data@cdfName <- "hgu133plus2"
  data_rma <- rma(data)
  exprs_data_rma <- exprs(data_rma)
  probe_ids <- rownames(exprs_data_rma)
  gene_symbols <- mapIds(hgu133plus2.db, keys = probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  exprs_data_rma_with_genes <- data.frame(GeneSymbol = gene_symbols, exprs_data_rma)
  file_name <- tools::file_path_sans_ext(basename(file))
  output_file <- file.path(output_dir, paste0(file_name, ".csv"))
  write.csv(exprs_data_rma_with_genes, file = output_file, row.names = FALSE) 
  print(paste("Results saved to", output_file))
}

data_dir <- "/Users/yangbolin/Downloads/03.BGI/08.AS/01.LPCA/00.rawdata/GSE41571_RAW"
output_dir <- "/Users/yangbolin/Downloads/03.BGI/08.AS/01.LPCA/02.2gene_name/GSE41571"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cel_files <- list.files(data_dir, pattern = "\\.CEL$", full.names = TRUE) 
print(cel_files)
if (length(cel_files) > 0) {
  sapply(cel_files, function(file) process_and_save_cel_file(file, output_dir))
} else {
  print("No CEL files found in the specified directory.")
}
input_dir <- "/Users/yangbolin/Downloads/03.BGI/08.AS/01.LPCA/02.2gene_name/GSE41571"
output_dir <- "/Users/yangbolin/Downloads/03.BGI/08.AS/01.LPCA/03.same_gene_intergrate/GSE41571"

if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

file_paths <- list.files(input_dir, pattern = ".csv", full.names = TRUE)

for (file_path in file_paths) {
  data <- read.csv(file_path, header = TRUE)
  summed_data <- aggregate(data[, 2], by = list(data$GeneSymbol), FUN = sum)
  output_file <- file.path(output_dir, basename(file_path))
  write.csv(summed_data, file = output_file, row.names = FALSE)
}

input_dir <- "/Users/yangbolin/Downloads/03.BGI/08.AS/01.LPCA/02.2gene_name/GSE41571"
output_dir <- "/Users/yangbolin/Downloads/03.BGI/08.AS/01.LPCA/03.same_gene_intergrate/GSE41571"
file_paths <- list.files(input_dir, pattern = ".csv", full.names = TRUE)

for (file_path in file_paths) {
  data <- read.csv(file_path, header = TRUE)
  summed_data <- aggregate(data[, 2], by = list(data$GeneSymbol), FUN = sum)
  output_file <- file.path(output_dir, basename(file_path))
  write.csv(summed_data, file = output_file, row.names = FALSE)
}

###Integrate to a matrix
input_dir <- "/Users/yangbolin/Downloads/03.BGI/08.AS/01.LPCA/03.same_gene_intergrate/GSE41571"
file_paths <- list.files(input_dir, pattern = ".csv", full.names = TRUE)
data_frames <- list()
for (file_path in file_paths) {
  data <- read.csv(file_path, header = TRUE)
  file_name <- tools::file_path_sans_ext(basename(file_path))
  row_names <- data[["Group.1"]]
  values <- data[[2]]
  data_frame <- data.frame(row_names, values, row.names = NULL)
  colnames(data_frame) <- c("Group.1", file_name)
  data_frames[[file_name]] <- data_frame
}
result_matrix <- Reduce(function(x, y) merge(x, y, by = "Group.1", all = TRUE), data_frames)
output_file <- "/Users/yangbolin/Downloads/03.BGI/08.AS/01.LPCA/04.matrix/GSE41571/merged_matrix.csv"
write.csv(result_matrix, file = output_file, row.names = FALSE)


library(pheatmap)
library(ggplot2)
gene <- read.csv("/Users/yangbolin/Downloads/03.BGI/08.AS/01.LPCA/MSMC_MAC_hotspot_with_human_genes.csv")
g1 <- list(unique(gene$Human.gene))
names(g1) <- "plaque_geneset"

result_matrix <- as.matrix(result_matrix)
params <- gsvaParam(exprData = result_matrix,geneSets = g1,kcdf = "Gaussian",maxDiff = F)
gsva.res1 <- gsva(params)
write.csv(gsva.res1, file = "gsva_results.csv", row.names = FALSE)

GSVA_matrix <- t(GSVA_matrix)
colnames(GSVA_matrix) <- GSVA_matrix[1, ]
GSVA_matrix <- GSVA_matrix[-1, ]


######Visualization

GSVA_matrix$atherosclerotic_plaque <- factor(GSVA_matrix$atherosclerotic_plaque, levels = c('Ruptured_Plaque', 'Stable_Plaque'))
p1 <- ggplot(data, aes(x = atherosclerotic_plaque, y = gsva_score, fill = atherosclerotic_plaque)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange", color = "red") +
  labs(title = "Gene Set Variation Analysis (GSVA)", x = "", y = "gsva_score") +
  geom_signif(comparisons = list(c('Ruptured_Plaque', 'Stable_Plaque')), test = “t.test”) + 

  scale_fill_manual(values = c('advanced' = '#66c2a4', 'early' = '#ce1256')) + 
  theme_classic()


GSVA_matrix$atherosclerotic_plaque <- factor(GSVA_matrix$atherosclerotic_plaque, levels = c('Ruptured_Plaque', 'Stable_Plaque'))
p2 <- ggplot(GSVA_matrix, aes(x = gsva_score, fill = atherosclerotic_plaque)) +
  geom_density(alpha = 0.5) +
  labs(title = "Gene Set Variation Analysis (GSVA) Density Plot", x = "GSVA Score", y = "Density") +
  scale_fill_manual(values = c('Ruptured_Plaque' = '#66c2a4', 'Stable_Plaque' = '#ce1256')) +
  theme_classic()
print(p2)


GSVA_matrix$atherosclerotic_plaque <- factor(GSVA_matrix$atherosclerotic_plaque, levels = c('Ruptured_Plaque', 'Stable_Plaque'))

p3 <- ggplot(GSVA_matrix, aes(x = atherosclerotic_plaque, y = gsva_score, fill = atherosclerotic_plaque)) +
  geom_boxplot() +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange", color = "red") +
  labs(title = "Gene Set Variation Analysis (GSVA) Box Plot", x = "", y = "GSVA Score") +
  geom_signif(comparisons = list(c('Ruptured_Plaque', 'Stable_Plaque')), test = "t.test") +
  scale_fill_manual(values = c('Ruptured_Plaque' = '#66c2a4', 'Stable_Plaque' = '#ce1256')) +
  theme_classic()
print(p3)
