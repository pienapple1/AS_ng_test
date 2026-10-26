#!/usr/bin/env Rscript
## =============================================================================
##  sc_qc_and_anno.R
##
##  Quality control, integration, clustering and level-1 annotation of the
##  mouse single-cell RNA-seq atlas.
##
##    Section 1-2  <- cell / gene QC
##    Section 3-4  <- normalisation, cell cycle, SCTransform, PCA, Harmony
##    Section 5-7  <- neighbours, UMAP/tSNE, resolution sweep, markers
##    Section 8-9  <- annotation review figures, level-1 cell type labels
##
##  INPUTS
##    <RAW_C4_DIR>                              the per-sample C4 BGI count
##                                              directories that were merged into
##                                              the atlas. Read when the directory
##                                              is reachable; otherwise Section 1
##                                              falls back to the pre-QC checkpoint
##                                              below.
##    <INPUT_QC_RDATA>                          -> object `merged_seurat`
##    <CELLCYCLE_RDATA>                         -> `g2m_genes`, `s_genes`
##    <FILTER_FUN_R>                            -> function filter_genes()
##
##  OUTPUTS (all written to OUT_DIR)
##    all_seurat_filtered.RData                   QC-passed object
##    filtered_seurat_cellcycle_score.RData       after SCTransform
##    seurat_integrated_seurat.RData              after Harmony
##    all_integrated_seurat_res.Rdata             all resolutions, pre-selection
##    all_integrated_seurat_res0.3_cluster.Rdata  resolution 0.3 selected
##    annotation_level1.Rdata                     FINAL annotated object
##    all_integrated_top{5,10,20,50}_DEG.csv      marker tables
##    all_integrated_DEG.csv                      full marker table
##    qc_summary.csv, celltype_level1_cluster_map.csv
##    *.pdf                                       figures
##    sc_QC_and_anno_sessionInfo.txt              software versions
##
##  ENVIRONMENT
##    R 4.1.1; Seurat 4.3.0 (SeuratObject 4.1.3); harmony 1.2.0;
##    clustree 0.5.0; ggsci 2.9; tidyverse 1.3.2; ggplot2 3.4.2
##
##  USAGE
##    conda activate R411
##    Rscript sc_QC_and_anno.R                       # run everything
##    Rscript sc_QC_and_anno.R --from 5              # resume at Section 5
##                                                   # (requires checkpoints)
##  Expected runtime on a compute node: several hours, peak RAM ~60 GB
## =============================================================================


## -----------------------------------------------------------------------------
## 0. Global setup
## -----------------------------------------------------------------------------

rm(list = ls())

.libPaths(c("/DATA/User/hemingmin/packages/R/x86_64-conda-linux-gnu-library/4.1",
            "/DATA/User/hemingmin/.conda/envs/R411/lib/R/library"))

suppressPackageStartupMessages({
  library(Seurat)      # 4.3.0
  library(harmony)     # 1.2.0
  library(clustree)    # 0.5.0
  library(ggsci)       # pal_igv()
  library(dplyr)       # %>%, group_by(), top_n()
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

## ---- paths ------------------------------------------------------------------
BASE_DIR         <- "/DATA/User/hemingmin/AS/NG_results/0.base"
SOURCE_DIR       <- "/DATA/User/hemingmin/AS/NG_results/source"
## The raw count directory can be pointed elsewhere with SC_QC_ANNO_RAW, so that
## the script does not have to be edited when the data are mounted under a
## different path.
RAW_C4_DIR       <- Sys.getenv("SC_QC_ANNO_RAW",
                               unset = "/DATA/User/liuyujia/MouseArtery/2.single_cell_rnaseq/mouse_sample/data")
INPUT_QC_RDATA   <- file.path(SOURCE_DIR, "00_QCpre_71235_cells.RData")
CELLCYCLE_RDATA  <- "/DATA/User/hemingmin/AS/merge5/data/mouse_cell_cycle_markers.RData"
FILTER_FUN_R     <- "/DATA/User/hemingmin/00.scripts/filter_genes.r"

## OUT_DIR defaults to the directory the original scripts wrote to.
OUT_DIR          <- Sys.getenv("SC_QC_ANNO_OUT", unset = BASE_DIR)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)

## ---- run control ------------------------------------------------------------
## Sections can be run individually with `--from N` / `--to M`. Skipping a
## section requires the checkpoint it would have written (see SAVE_CHECKPOINTS).
.cli <- commandArgs(trailingOnly = TRUE)
.cli_int <- function(flag, default) {
  i <- match(flag, .cli)
  if (is.na(i) || i == length(.cli)) return(default)
  as.integer(.cli[i + 1])
}
RUN_FROM          <- .cli_int("--from", 1)     # first section to execute (1-9)
RUN_UNTIL         <- .cli_int("--to",   9)     # last section to execute
SAVE_CHECKPOINTS  <- !("--no-checkpoints" %in% .cli)
LOG_TO_FILE       <- TRUE   # tee console output to a timestamped log
PARAM_SEED        <- 123    # seed for clustering / UMAP / tSNE

## ---- convention -------------------------------------------------------------
## Every analysis parameter is declared at the top of the section that uses it,
## immediately above the code it controls (see `## parameters` blocks below).

## ---- helpers ----------------------------------------------------------------
## active(n): TRUE when section n lies inside the requested --from / --to range.
active <- function(n) RUN_FROM <= n && n <= RUN_UNTIL
section <- function(n, title) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat(sprintf(" SECTION %s | %s\n", n, title))
  cat(strrep("=", 78), "\n\n", sep = "")
}
tick <- function(...) cat(sprintf("[ok] %s\n", paste0(...)))
checkpoint <- function(obj, filename, name = deparse(substitute(obj))) {
  if (!SAVE_CHECKPOINTS) return(invisible(NULL))
  assign(name, obj)
  save(list = name, file = file.path(OUT_DIR, filename), envir = environment())
  tick("checkpoint -> ", filename, sprintf(" (%.2f GB)", file.size(file.path(OUT_DIR, filename)) / 1e9))
  invisible(NULL)
}

if (LOG_TO_FILE) {
  log_path <- file.path(OUT_DIR, sprintf("sc_QC_and_anno_%s.log",
                                         format(Sys.time(), "%Y%m%d_%H%M%S")))
  log_con  <- file(log_path, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  cat("Log file:", log_path, "\n")
}

set.seed(PARAM_SEED)
cat("Working directory:", OUT_DIR, "\n")
cat("Started:", format(Sys.time()), "\n")


## =============================================================================
## 1. Assemble the pre-QC expression matrix and apply the sample design
## =============================================================================

## ---- parameters --------------------------------------------------------------
QC_MIN_CELLS_PER_SAMPLE <- 800   # discard libraries with fewer cells than this

## ---- 1.1 library-level assembly from the per-sample count matrices ----------
## Each sample was delivered as a count directory (BGI format, hence
## gene.column = 1). Libraries with < QC_MIN_CELLS_PER_SAMPLE cells were
## discarded before merging, and all remaining libraries were merged into a
## single object.
build_pre_qc_object <- function(raw_dir, min_cells = QC_MIN_CELLS_PER_SAMPLE) {
  samples <- list.files(raw_dir)
  cat(sprintf("Reading %d count directories from %s\n", length(samples), raw_dir))
  stopifnot(length(samples) > 1)

  sce_list <- lapply(samples, function(s) {
    counts <- Read10X(data.dir = file.path(raw_dir, s), gene.column = 1)
    obj <- CreateSeuratObject(counts = counts, project = s)
    obj$sample <- s
    obj
  })

  n_before <- length(sce_list)
  sce_list <- sce_list[vapply(sce_list, ncol, integer(1)) > min_cells]
  cat(sprintf("Libraries: %d -> %d retained (> %d cells)\n",
              n_before, length(sce_list), min_cells))
  obj <- merge(x = sce_list[[1]], y = sce_list[-1])
  cat(sprintf("Merged object: %d genes x %d cells\n", nrow(obj), ncol(obj)))

  ## Per-cell complexity and mitochondrial content. "^mt-" is the mouse
  ## mitochondrial prefix;
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  obj$mitoRatio <- PercentageFeatureSet(object = obj, pattern = "^mt-") / 100

  md <- obj@meta.data
  md$cells <- rownames(md)
  md <- md %>%
    dplyr::rename(seq_folder = orig.ident, nUMI = nCount_RNA, nGene = nFeature_RNA)
  obj@meta.data <- md
  obj
}

## ---- 1.2 sample design ------------------------------------------------------
## Library-to-group table. `type` separates the atherosclerotic (AS) animals
## from their littermate controls; `timepoint_2` additionally separates the two
## 20-weeks cohorts, which were collected together and are compared separately
## in the composition analyses.
SAMPLE_TYPE_CON <- c("LN20_2", "LN20_3", "LW0_2", "LW0_3", "SSLW0_1", "SSLW0_2",
                     "fsLW0_2", "LLW0_1", "LLW0_2")
SAMPLE_TYPE_AS  <- c("LW2_1", "LW2_2", "LW2_3", "LW6_1", "LW6_2", "LW6_3",
                     "LW10_2", "LW10_3", "LW10_3_1", "LW10_4", "LW20_4", "LW20_5",
                     "SLW10_2", "SLW10_3", "SLW20_1", "SSLW2_1",
                     "fsLW10_1", "fsLW10_2", "fsLW10_3", "fsLW10_4", "fsLW10_5",
                     "fsLW10_7", "fsLW10_8", "fsLW10_9", "fsLW20_1", "fsLW20_2",
                     "tsLW2_1", "tsLW2_2", "tsLW2_3", "tsLW6_4", "tsLW6_5")

SAMPLE_TIME <- list(
  `0`  = c("LW0_2", "LW0_3", "SSLW0_1", "SSLW0_2", "fsLW0_2", "LLW0_1", "LLW0_2"),
  `2`  = c("LW2_1", "LW2_2", "LW2_3", "SSLW2_1", "tsLW2_1", "tsLW2_2", "tsLW2_3"),
  `6`  = c("LW6_1", "LW6_2", "LW6_3", "tsLW6_4", "tsLW6_5"),
  `10` = c("LW10_2", "LW10_3", "LW10_3_1", "LW10_4", "SLW10_2", "SLW10_3",
           "fsLW10_1", "fsLW10_2", "fsLW10_3", "fsLW10_4", "fsLW10_5",
           "fsLW10_7", "fsLW10_8", "fsLW10_9"),
  `20` = c("LN20_2", "LN20_3", "LW20_4", "LW20_5", "SLW20_1", "fsLW20_1", "fsLW20_2")
)
SAMPLE_TIME2 <- list(
  `0`      = SAMPLE_TIME[["0"]],
  `2`      = SAMPLE_TIME[["2"]],
  `6`      = SAMPLE_TIME[["6"]],
  `10`     = SAMPLE_TIME[["10"]],
  `20_As`  = c("LW20_4", "LW20_5", "SLW20_1", "fsLW20_1", "fsLW20_2"),
  `20_Con` = c("LN20_2", "LN20_3")
)
SAMPLE_BATCH <- list(
  `1` = c("LN20_2", "LN20_3", "LW0_2", "LW0_3", "LW2_1", "LW2_2", "LW2_3",
          "LW6_1", "LW6_2", "LW6_3", "LW10_2", "LW10_3", "LW10_3_1", "LW10_4",
          "LW20_4", "LW20_5"),
  `2` = c("SLW10_2", "SLW10_3", "SLW20_1"),
  `3` = c("SSLW0_1", "SSLW0_2", "SSLW2_1"),
  `4` = c("tsLW2_1", "tsLW2_2", "tsLW2_3", "tsLW6_4", "tsLW6_5"),
  `5` = c("fsLW0_2", "fsLW10_1", "fsLW10_2", "fsLW10_3", "fsLW10_4", "fsLW10_5",
          "fsLW20_1", "fsLW20_2"),
  `6` = c("fsLW10_7", "fsLW10_8", "fsLW10_9"),
  `7` = c("LLW0_1", "LLW0_2")
)

## The design columns are rebuilt from `sample` so that the tables above are the
## single source of truth. This is a no-op for the pre-QC checkpoint, which
## already carries identical columns, and it makes the script independent of
## which of the two inputs was used.
annotate_design <- function(md) {
  if (!"sample" %in% colnames(md)) {
    md$sample <- if ("seq_folder" %in% colnames(md)) md$seq_folder else md$orig.ident
  }

  md$type <- NA_character_
  md$type[md$sample %in% SAMPLE_TYPE_CON] <- "Con"
  md$type[md$sample %in% SAMPLE_TYPE_AS]  <- "AS"

  md$timepoint <- NA_character_
  for (g in names(SAMPLE_TIME)) {
    md$timepoint[md$sample %in% SAMPLE_TIME[[g]]] <- g
  }
  md$timepoint <- factor(md$timepoint, levels = c("0", "2", "6", "10", "20"))

  md$timepoint_2 <- NA_character_
  for (g in names(SAMPLE_TIME2)) {
    md$timepoint_2[md$sample %in% SAMPLE_TIME2[[g]]] <- g
  }
  ## Level order follows the order used in the composition figures.
  md$timepoint_2 <- factor(md$timepoint_2,
                           levels = c("20_Con", "0", "2", "6", "10", "20_As"))

  md$batch <- NA_character_
  for (b in names(SAMPLE_BATCH)) {
    md$batch[md$sample %in% SAMPLE_BATCH[[b]]] <- b
  }
  md
}

if (active(1)) {
  section(1, "Assemble the pre-QC object and apply the sample design")

  ## Prefer the raw count directories; fall back to the pre-QC checkpoint when
  ## they are not mounted (the original location is a shared project directory
  ## that is no longer available on this node).
  raw_available <- dir.exists(RAW_C4_DIR) && length(list.files(RAW_C4_DIR)) > 1

  if (raw_available) {
    input_used <- "raw C4 count matrices"
    merged_seurat <- build_pre_qc_object(RAW_C4_DIR)
  } else {
    input_used <- "pre-QC checkpoint (00_QCpre_71235_cells.RData)"
    stopifnot(file.exists(INPUT_QC_RDATA))
    message("Raw count directory not found: ", RAW_C4_DIR)
    message("Falling back to the pre-QC checkpoint ", basename(INPUT_QC_RDATA))
    load(INPUT_QC_RDATA)            # -> merged_seurat (71,235 cells)
  }
  cat("Input used:", input_used, "\n")

  ## The design columns (including `sample`) are applied before any of them are
  ## read, because build_pre_qc_object() only produces `seq_folder`.
  merged_seurat@meta.data <- annotate_design(merged_seurat@meta.data)
  cat("Object:", nrow(merged_seurat), "genes x", ncol(merged_seurat), "cells,",
      length(unique(merged_seurat$sample)), "libraries\n")

  ## Sample bookkeeping: any library not covered by SAMPLE_TIME2 would silently
  ## drop out of the composition figures, and an empty group would break them.
  covered   <- unlist(SAMPLE_TIME2, use.names = FALSE)
  uncovered <- setdiff(unique(as.character(merged_seurat$sample)), covered)
  if (length(uncovered)) {
    message("NOTE: libraries not assigned to a time point: ",
            paste(uncovered, collapse = ", "))
  }
  empty <- names(SAMPLE_TIME2)[vapply(SAMPLE_TIME2, function(g)
    !any(as.character(merged_seurat$sample) %in% g), logical(1))]
  if (length(empty)) {
    message("NOTE: groups with no cells in the object: ", paste(empty, collapse = ", "))
  }

  cat("Cells per type:\n");        print(table(merged_seurat$type))
  cat("Cells per timepoint:\n");   print(table(merged_seurat$timepoint))
  cat("Cells per timepoint_2:\n"); print(table(merged_seurat$timepoint_2, useNA = "ifany"))
  cat("Cells per batch:\n");       print(table(merged_seurat$batch))
  tick("sample design applied")
}


## =============================================================================
## 2. Gene filtering and cell quality control
## =============================================================================
if (active(2)) {
  section(2, "Gene filtering and cell QC")

  ## ---- parameters ------------------------------------------------------------
  QC_MIN_UMI       <- 500     # lower bound on nUMI (counts per cell)
  QC_MAX_UMI_QUANT <- 0.95    # upper bound on nUMI, as a quantile of the data
  QC_MIN_GENE      <- 1000    # lower bound on nGene (detected genes per cell)
  QC_MIN_LOG10_GENES_PER_UMI <- 0.80   # complexity floor, removes empty droplets
  QC_MAX_MITO      <- 0.05    # mitochondrial fraction ceiling (mouse, "^mt-")

  ## Section 1 has no checkpoint of its own; `merged_seurat` must already exist.
  if (!exists("merged_seurat")) {
    stop("Section 2 requires `merged_seurat`. Section 1 has no checkpoint, ",
         "so start with --from 1.")
  }
  if (!exists("filter_genes")) {
    stopifnot(file.exists(FILTER_FUN_R))
    source(FILTER_FUN_R)            # -> filter_genes()
  }

  n_cells_pre  <- ncol(merged_seurat)
  n_genes_pre  <- nrow(merged_seurat)

  ## Remove mitochondrial, ribosomal, heat-shock, predicted and Rik genes that
  ## inflate per-cell complexity without carrying cell-type information.
  filtered_seurat <- filter_genes(merged_seurat)
  cat(sprintf("Genes: %d -> %d after filter_genes()\n",
              n_genes_pre, nrow(filtered_seurat)))

  ## Cell-level thresholds. The upper bound on nUMI is the 95th percentile of
  ## the gene-filtered data, i.e. it is computed on the same matrix that is
  ## subsequently subset; it removes putative doublets.
  data      <- filtered_seurat@meta.data
  max_UMI   <- quantile(data$nUMI, QC_MAX_UMI_QUANT)
  cat(sprintf("nUMI 95th percentile (upper bound): %.0f\n", max_UMI))

  filtered_seurat <- subset(
    x = filtered_seurat,
    subset = (nUMI >= QC_MIN_UMI) &
             (nUMI <  max_UMI) &
             (nGene >= QC_MIN_GENE) &
             (log10GenesPerUMI > QC_MIN_LOG10_GENES_PER_UMI) &
             (mitoRatio < QC_MAX_MITO))

  n_cells_post <- ncol(filtered_seurat)
  cat(sprintf("Cells: %d -> %d retained (%.1f%%)\n",
              n_cells_pre, n_cells_post, 100 * n_cells_post / n_cells_pre))
  print(table(filtered_seurat$timepoint_2, useNA = "ifany"))

  ## Summary table for the methods section / point-by-point response.
  qc_summary <- data.frame(
    step      = c("pre-QC", "after filter_genes() + cell thresholds"),
    n_genes   = c(n_genes_pre, nrow(filtered_seurat)),
    n_cells   = c(n_cells_pre, n_cells_post),
    retained  = c(1, n_cells_post / n_cells_pre),
    stringsAsFactors = FALSE)
  write.csv(qc_summary, file.path(OUT_DIR, "qc_summary.csv"), row.names = FALSE)
  print(qc_summary, row.names = FALSE)

  checkpoint(filtered_seurat, "all_seurat_filtered.RData")
}


## =============================================================================
## 3. Normalisation, cell-cycle scoring and variance stabilisation
## =============================================================================
if (active(3)) {
  section(3, "Normalise, score cell cycle, SCTransform")

  if (RUN_FROM > 2) load(file.path(OUT_DIR, "all_seurat_filtered.RData"))

  stopifnot(file.exists(CELLCYCLE_RDATA))
  load(CELLCYCLE_RDATA)             # -> g2m_genes, s_genes

  filtered_seurat <- NormalizeData(filtered_seurat, verbose = TRUE)

  ## Cell-cycle scores are regressed out in SCTransform below, so that
  ## proliferating and quiescent cells of the same identity are not split.
  filtered_seurat <- CellCycleScoring(filtered_seurat,
                                      g2m.features = g2m_genes,
                                      s.features   = s_genes)
  print(table(filtered_seurat$Phase))

  filtered_seurat <- SCTransform(filtered_seurat,
                                 vars.to.regress = c("S.Score", "G2M.Score"),
                                 verbose = TRUE)

  checkpoint(filtered_seurat, "filtered_seurat_cellcycle_score.RData")
}


## =============================================================================
## 4. PCA and Harmony integration across samples
## =============================================================================
if (active(4)) {
  section(4, "PCA and Harmony batch correction")

  if (RUN_FROM > 3) load(file.path(OUT_DIR, "filtered_seurat_cellcycle_score.RData"))

  ## `sample` (33 libraries) is used as the batch variable.
  seurat_integrated <- RunPCA(filtered_seurat, verbose = FALSE)

  set.seed(PARAM_SEED)              # reproducible Harmony embedding
  seurat_integrated <- seurat_integrated %>%
    RunHarmony("sample", plot_convergence = TRUE)

  cat("Harmony components:", ncol(Embeddings(seurat_integrated, "harmony")), "\n")

  checkpoint(seurat_integrated, "seurat_integrated_seurat.RData")
}


## =============================================================================
## 5. Neighbours, UMAP / tSNE and the clustering resolution sweep
## =============================================================================
if (active(5)) {
  section(5, "Graph construction, embedding and resolution sweep")

  ## ---- parameters ------------------------------------------------------------
  N_DIMS      <- 1:50                    # Harmony components used downstream
  CLUSTER_RES <- seq(0.1, 1, by = 0.1)   # resolutions evaluated for cluster stability

  if (RUN_FROM > 4) load(file.path(OUT_DIR, "seurat_integrated_seurat.RData"))

  seurat_integrated <- FindNeighbors(seurat_integrated,
                                     reduction = "harmony",
                                     dims = N_DIMS,
                                     verbose = FALSE)
  set.seed(PARAM_SEED)
  seurat_integrated <- RunUMAP(seurat_integrated,
                               reduction = "harmony",
                               dims = N_DIMS,
                               verbose = FALSE)

  set.seed(PARAM_SEED)
  seurat_integrated <- RunTSNE(seurat_integrated,
                               reduction = "harmony",
                               dims = N_DIMS)

  ## Resolution sweep, used to justify the choice of resolution 0.3.
  set.seed(PARAM_SEED)
  for (res in CLUSTER_RES) {
    seurat_integrated <- FindClusters(seurat_integrated,
                                      resolution = res,
                                      verbose = FALSE)
  }

  checkpoint(seurat_integrated, "all_integrated_seurat_res.Rdata")

  ## Cluster stability across resolutions, saved as a named figure.
  p_clustree <- clustree(seurat_integrated)
  ggsave(p_clustree, filename = "clustree_resolution_sweep.pdf",
         path = OUT_DIR, width = 12, height = 10)
  tick("clustree figure written")

  for (res in CLUSTER_RES) {
    column <- paste0("SCT_snn_res.", res)
    cat(sprintf("resolution %.1f -> %d clusters\n", res,
                length(unique(seurat_integrated@meta.data[[column]]))))
  }
}


## =============================================================================
## 6. Select resolution 0.3 and fix factor levels for plotting
## =============================================================================
if (active(6)) {
  section(6, "Select resolution 0.3")

  ## ---- parameters ------------------------------------------------------------
  SELECTED_RES <- "SCT_snn_res.0.3"   # column produced by the Section 5 sweep

  if (RUN_FROM > 5) load(file.path(OUT_DIR, "all_integrated_seurat_res.Rdata"))

  stopifnot(SELECTED_RES %in% colnames(seurat_integrated@meta.data))

  Idents(seurat_integrated) <- seurat_integrated@meta.data[[SELECTED_RES]]
  seurat_integrated$seurat_clusters <- seurat_integrated@meta.data[[SELECTED_RES]]

  ## Plotting order: control first, then the time series, then the aged AS group.
  seurat_integrated$timepoint <- factor(seurat_integrated$timepoint,
                                        levels = c("0", "2", "6", "10", "20"))
  seurat_integrated$timepoint_2 <- factor(seurat_integrated$timepoint_2,
                                          levels = c("20_Con", "0", "2", "6",
                                                     "10", "20_As"))

  cat("Clusters at resolution 0.3:", length(unique(seurat_integrated$seurat_clusters)), "\n")
  print(table(seurat_integrated$seurat_clusters))
  cat("Cells per time point:\n")
  print(table(seurat_integrated$timepoint_2))

  checkpoint(seurat_integrated, "all_integrated_seurat_res0.3_cluster.Rdata")
}


## =============================================================================
## 7. Cluster marker genes
## =============================================================================
if (active(7)) {
  section(7, "Marker gene discovery")

  ## ---- parameters ------------------------------------------------------------
  MARKER_TEST    <- "MAST"    # test used by the original analysis
  MARKER_MIN_PCT <- 0.1       # gene must be detected in >= 10% of a cluster
  MARKER_LOGFC   <- 0.25      # minimum log2 fold change
  TOP_N_TABLES   <- c(5, 10, 20, 50)   # marker tables written to disk

  if (RUN_FROM > 6) load(file.path(OUT_DIR, "all_integrated_seurat_res0.3_cluster.Rdata"))

  ## Markers are called on the normalised RNA assay rather than on the
  ## SCTransform residuals, so that log2 fold changes are on a conventional
  ## scale and the tables are comparable with the published ones.
  DefaultAssay(seurat_integrated) <- "RNA"
  seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

  ## ScaleData on all genes is retained from the original script.
  all.genes <- rownames(seurat_integrated)
  seurat_integrated <- ScaleData(seurat_integrated, features = all.genes,
                                 verbose = FALSE)

  seurat_integrated.markers <- FindAllMarkers(seurat_integrated,
                                              test.use       = MARKER_TEST,
                                              only.pos       = TRUE,
                                              min.pct        = MARKER_MIN_PCT,
                                              logfc.threshold = MARKER_LOGFC,
                                              verbose        = TRUE)
  cat("Marker rows:", nrow(seurat_integrated.markers), "\n")
  for (n in TOP_N_TABLES) {
    tab <- seurat_integrated.markers %>% group_by(cluster) %>% top_n(n = n, wt = avg_log2FC)
    write.csv(tab, file.path(OUT_DIR, sprintf("all_integrated_top%d_DEG.csv", n)),
              row.names = FALSE)
    assign(sprintf("top%d", n), tab)
    tick(sprintf("top%d marker table: %d rows", n, nrow(tab)))
  }
  write.csv(seurat_integrated.markers, file.path(OUT_DIR, "all_integrated_DEG.csv"),
            row.names = FALSE)
}


## =============================================================================
## 8. Annotation review figures
## =============================================================================
if (active(8)) {
  section(8, "Annotation review figures")

  ## ---- parameters ------------------------------------------------------------
  CLUSTER_COLOUR_N <- 18      # clusters at resolution 0.3 (18 colours, pal_igv)
  PALETTE_CLUSTERS <- pal_igv("default")(CLUSTER_COLOUR_N)

  if (RUN_FROM > 7) load(file.path(OUT_DIR, "all_integrated_seurat_res0.3_cluster.Rdata"))
  n_clusters <- length(unique(seurat_integrated$seurat_clusters))
  if (n_clusters != length(PALETTE_CLUSTERS)) {
    warning(sprintf("Palette holds %d colours but the object has %d clusters; regenerating.",
                    length(PALETTE_CLUSTERS), n_clusters), call. = FALSE)
    PALETTE_CLUSTERS <- pal_igv("default")(n_clusters)
  }

  ## Marker panel used for manual review of the clusters, grouped by lineage.
  ## Gene symbols were checked against their mouse standard names.
  features <- c(
    "Pecam1", "Vwf", "Edn1",                                              # endothelial
    "Pi16", "Clec3b", "Gpx3", "Serping1", "Dcn", "Col1a2", "Col1a1", "Wif1",  # fibroblast
    "Acta2", "Myh11", "Mylk", "Myl9", "Cnn1", "Tagln", "Tpm2", "Flna",    # SMC
    "Lum", "Ly6a", "Vcam1", "Ly6c1", "Spp1", "Fn1", "Lcn2", "Timp1",      # modulated SMC
    "Cd14", "Fcgr3", "Cx3cr1", "Ly6c2", "Ly6g", "Ccr2",                   # monocyte
    "Cd68", "Adgre1",                                                     # macrophage
    "F13a1", "Mrc1", "Lyve1", "Gas6", "Folr2", "Cbr2", "Selenop", "Ccl9", # resident-like macrophage
    "Il1b", "Nlrp3", "Tlr2", "Tnf", "Ccl3", "Ccl4", "Ccl5", "Cebpb",      # inflammatory macrophage
    "Trem2", "Cd9", "Itgax", "Lgals3", "Ctsb", "Ctsd", "Ctsz",            # Trem2 macrophage
    "Ndufa4l2", "Vtn", "Notch3", "Colec11", "Higd1b",                     # pericyte
    "Upk3b", "Msln", "Fxyd3", "Gpm6a", "Nkain4", "Krt19", "Krt7",         # epithelial
    "Plp1", "Kcna1", "Gfra3", "Egfl8", "Mbp", "Cadm4",                    # neuron
    "Cd79a", "Ighm", "Cd19", "Cd79b",                                     # B cell
    "Cd2", "Cd3d", "Cd3e", "Cd3g"                                         # T cell
  )
  DefaultAssay(seurat_integrated) <- "RNA"
  missing_features <- setdiff(features, rownames(seurat_integrated))
  if (length(missing_features)) {
    message("Marker panel genes absent from the object: ",
            paste(missing_features, collapse = ", "))
  }
  features <- intersect(features, rownames(seurat_integrated))

  ## 8.1 lineage marker panel
  p_panel <- DotPlot(seurat_integrated, features = features,
                     col.min = 0, cols = c("lightgrey", "red")) +
    coord_flip() +
    theme(axis.text.x = element_text(size = 10, hjust = 1),
          axis.text.y = element_text(size = 15))
  ggsave(p_panel, filename = "all_seurat_integrated_dotplot_markers.pdf",
         path = OUT_DIR, width = 10, height = 15)
  tick("lineage marker dot plot written")

  ## 8.2 top-5 marker genes per cluster (dot plots, both orientations)
  top5 <- read.csv(file.path(OUT_DIR, "all_integrated_top5_DEG.csv"))

  p_top5_wide <- DotPlot(seurat_integrated, features = unique(top5$gene),
                         col.min = 0, cols = c("lightgrey", "red")) +
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 12))
  ggsave(p_top5_wide, filename = "top5_dotplot_harmony.pdf",
         path = OUT_DIR, width = 16, height = 6)

  p_top5_tall <- DotPlot(seurat_integrated, features = unique(top5$gene),
                         col.min = 0, cols = c("lightgrey", "red")) +
    coord_flip() +
    theme(axis.text.x = element_text(size = 12, hjust = 1),
          axis.text.y = element_text(size = 10, face = "bold"))
  ggsave(p_top5_tall, filename = "top5_dotplot.pdf",
         path = OUT_DIR, width = 6, height = 16)
  tick("top-5 marker dot plots written")

  ## 8.3 cluster embeddings
  ggsave(DimPlot(seurat_integrated, group.by = "seurat_clusters", label = TRUE,
                 label.size = 5, cols = PALETTE_CLUSTERS, reduction = "tsne"),
         filename = "tsne.pdf", path = OUT_DIR, width = 8, height = 8)

  ggsave(DimPlot(seurat_integrated, group.by = "seurat_clusters", label = TRUE,
                 label.size = 5, cols = PALETTE_CLUSTERS, reduction = "umap"),
         filename = "umap.pdf", path = OUT_DIR, width = 8, height = 8)

  ggsave(DimPlot(seurat_integrated, group.by = "seurat_clusters",
                 split.by = "timepoint_2", label = TRUE, label.size = 5,
                 cols = PALETTE_CLUSTERS, reduction = "tsne"),
         filename = "tsne_split_timepoint2.pdf", path = OUT_DIR, width = 18, height = 4)

  ggsave(DimPlot(seurat_integrated, group.by = "seurat_clusters",
                 split.by = "timepoint_2", label = TRUE, label.size = 5,
                 cols = PALETTE_CLUSTERS, reduction = "umap"),
         filename = "umap_split_timepoint2.pdf", path = OUT_DIR, width = 18, height = 4)
  tick("cluster embeddings written")

  ## 8.4 cluster composition across time points
  p_composition <- seurat_integrated@meta.data %>%
    ggplot(aes(x = timepoint_2, fill = seurat_clusters)) +
    geom_bar(aes(fill = seurat_clusters), position = "fill") +
    scale_fill_manual(values = PALETTE_CLUSTERS) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("NCells")
  ggsave(p_composition, filename = "seurat_clusters_nCell_timepoint2.pdf",
         path = OUT_DIR, width = 10, height = 8)
  tick("cluster composition plot written")
}


## =============================================================================
## 9. Level-1 cell type annotation
## =============================================================================
if (active(9)) {
  section(9, "Level-1 annotation")

  ## ---- parameters ------------------------------------------------------------
  ## Level order used for every downstream table and figure.
  CELLTYPE_LEVELS <- c("SMC", "Modulated_SMC", "Fibroblast",
                       "Macrophage", "EC", "T cell")

  ## Cluster-to-cell-type assignment for level-1 annotation. The mapping is
  ## stated as an explicit table so that every label can be traced to one cluster.
  CLUSTER_TO_LEVEL1 <- c(
    `0`  = "SMC",           `1`  = "SMC",           `2`  = "Modulated_SMC",
    `3`  = "Fibroblast",    `4`  = "SMC",           `5`  = "SMC",
    `6`  = "Macrophage",    `7`  = "SMC",           `8`  = "SMC",
    `9`  = "EC",            `10` = "SMC",           `11` = "T cell",
    `12` = "SMC",           `13` = "SMC",           `14` = "SMC",
    `15` = "SMC",           `16` = "SMC",           `17` = "SMC"
  )

  if (RUN_FROM > 8) load(file.path(OUT_DIR, "all_integrated_seurat_res0.3_cluster.Rdata"))

  ## The original script cross-tabulated the new clusters against a
  ## `celltype_level1` column that was already present in the R object it had
  ## in memory. That column is NOT part of the resolution-0.3 checkpoint, so the
  ## cross-tabulation is now conditional: when the column is absent, the mapping
  ## below defines the labels de novo.
  if ("celltype_level1" %in% colnames(seurat_integrated@meta.data)) {
    cat("Cross-tabulation of existing celltype_level1 vs clusters:\n")
    print(table(seurat_integrated@meta.data$celltype_level1,
                seurat_integrated@meta.data$seurat_clusters))
  } else {
    message("No pre-existing `celltype_level1` column in the object: ",
            "labels are assigned from CLUSTER_TO_LEVEL1 (Section 9).")
  }

  cluster_ids <- levels(factor(seurat_integrated$seurat_clusters))
  unmapped    <- setdiff(cluster_ids, names(CLUSTER_TO_LEVEL1))
  if (length(unmapped)) {
    stop("No cell type assigned to cluster(s): ", paste(unmapped, collapse = ", "),
         ". Update CLUSTER_TO_LEVEL1 in Section 9.")
  }

  ## An explicit vector lookup replaces RenameIdents() so that the
  ## cluster-to-label assignment is a single auditable table.
  label_vec <- unname(CLUSTER_TO_LEVEL1[cluster_ids])
  names(label_vec) <- cluster_ids

  seurat_integrated$celltype_level1 <- factor(
    label_vec[as.character(seurat_integrated$seurat_clusters)],
    levels = CELLTYPE_LEVELS)
  Idents(seurat_integrated) <- "celltype_level1"

  cat("Cells per level-1 cell type:\n")
  print(table(seurat_integrated$celltype_level1))

  ## Audit table: cluster, cell type and number of cells.
  cluster_map <- data.frame(
    cluster   = cluster_ids,
    n_cells   = as.integer(table(seurat_integrated$seurat_clusters)[cluster_ids]),
    celltype_level1 = label_vec,
    stringsAsFactors = FALSE)
  write.csv(cluster_map, file.path(OUT_DIR, "celltype_level1_cluster_map.csv"),
            row.names = FALSE)
  print(cluster_map, row.names = FALSE)

  ## Release the RNA scale.data matrix, which is not used downstream and would
  ## add ~8.6 GB to the saved object. The original workflow avoided this by
  ## reloading a checkpoint that predates ScaleData(); doing it explicitly keeps
  ## the saved object identical in content while staying within one session.
  rna_assay <- seurat_integrated[["RNA"]]
  rna_assay@scale.data <- matrix()
  seurat_integrated[["RNA"]] <- rna_assay
  DefaultAssay(seurat_integrated) <- "RNA"

  ## Write the annotated object under the name used by every downstream script
  final_path <- file.path(OUT_DIR, "annotation_level1.Rdata")
  if (file.exists(final_path)) {
    final_path <- file.path(OUT_DIR, sprintf("annotation_level1_rerun_%s.Rdata",
                                             format(Sys.time(), "%Y%m%d_%H%M%S")))
    message("annotation_level1.Rdata already exists; writing to ",
            basename(final_path), " instead. Rename it manually if you intend ",
            "to replace the existing annotated object.")
  }
  save(seurat_integrated, file = final_path)
  tick("annotated object -> ", basename(final_path),
       sprintf(" (%.2f GB)", file.size(final_path) / 1e9))
}


## =============================================================================
## 10. Provenance
## =============================================================================
section(10, "Session information")
si <- capture.output(sessionInfo())
cat(paste(si, collapse = "\n"), "\n")
writeLines(si, file.path(OUT_DIR, "sc_QC_and_anno_sessionInfo.txt"))
cat("\nFinished:", format(Sys.time()), "\n")

if (LOG_TO_FILE) {
  sink(type = "message")
  sink()
}
