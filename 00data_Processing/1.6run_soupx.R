#!/usr/bin/env Rscript
# =============================================================================
# run_soupx.R
#
#
# Input (produced by the upstream export script; one <slug> sub-directory per lineage):
#   <export_dir>/soup_profile_SMC.csv           # gene,est,counts — pooled SMC soup profile
#   <export_dir>/<slug>/matrix.mtx              # lineage counts (MatrixMarket)
#   <export_dir>/<slug>/genes.tsv, barcodes.tsv
#   <export_dir>/<slug>/obs_<slug>.csv          # cell metadata, with seurat_clusters
#   <export_dir>/<slug>/umap_<slug>.csv         # UMAP_1, UMAP_2 (optional, attached
#                                               #   to the SoupChannel object)
#   <export_dir>/<slug>/soup_profile_<slug>.csv # self-lineage soup profile
#
# Output (<out_dir>/<slug>/):
#   matrix_adjusted.mtx / genes.tsv / barcodes.tsv   # corrected counts (status=OK)
#   rho_per_cell.csv, qc_delta.csv, gene_delta.csv, rho_summary.csv, run_<slug>.log
#   <out_dir>/rho_summary_all.csv                    # summary across all lineages
#
# Requirements: R >= 4.1; SoupX (>= 1.6); Matrix
# Usage:
#   Rscript run_soupx_perLineage.R <export_dir> <out_dir>
# Example:
#   Rscript run_soupx_perLineage.R 01_export 02_soupx_out
# =============================================================================
suppressPackageStartupMessages({
  library(SoupX)
  library(Matrix)
})

## ---------------- Arguments ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_soupx_perLineage.R <export_dir> <out_dir>\n",
       "  <export_dir>: upstream export directory (per-lineage sub-directories ",
       "plus soup_profile_SMC.csv)\n",
       "  <out_dir>   : output directory (created automatically)", call. = FALSE)
}
EXP  <- normalizePath(args[[1]], mustWork = TRUE)
OUTD <- args[[2]]

## ---------------- Configuration ----------------
SMC_GENES <- c("Acta2", "Myh11", "Mylk", "Myl9", "Cnn1", "Tagln", "Tpm2", "Flna")
SLUGS     <- c("mSMC", "Fib", "Macro", "EC", "Tcell")
RHO_LO    <- 0.01
RHO_HI    <- 0.50

dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
cat("[env] R", as.character(getRversion()), "| SoupX", as.character(packageVersion("SoupX")),
    "| Matrix", as.character(packageVersion("Matrix")), "\n")

read_profile <- function(path, genes) {
  p <- read.csv(path)
  rownames(p) <- p$gene
  p <- p[genes, c("est", "counts")]
  p[is.na(p)] <- 0
  p
}

estimate_rho <- function(ch, gl) {
  ute <- estimateNonExpressingCells(ch, gl, clusters = FALSE)
  ne  <- sum(as.logical(ute))
  r <- tryCatch({
    ch2 <- calculateContaminationFraction(ch, gl, ute)
    list(rho = median(ch2$metaData$rho, na.rm = TRUE), ch = ch2, ne = ne, err = NA_character_)
  }, error = function(e) {
    list(rho = NA_real_, ch = NULL, ne = ne, err = conditionMessage(e))
  })
  r
}

summary_rows <- list()

for (slug in SLUGS) {
  cat("\n==================== ", slug, " ====================\n")
  d  <- file.path(EXP, slug)
  od <- file.path(OUTD, slug)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  con <- file(file.path(od, paste0("run_", slug, ".log")), open = "wt")
  say <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, con) }

  tryCatch({
    cnt <- as(readMM(file.path(d, "matrix.mtx")), "CsparseMatrix")
    g   <- readLines(file.path(d, "genes.tsv"))
    bc  <- readLines(file.path(d, "barcodes.tsv"))
    dimnames(cnt) <- list(g, bc)
    say("[read] ", nrow(cnt), " genes x ", ncol(cnt), " cells, nnz = ", length(cnt@x))

    soup <- SoupChannel(cnt, cnt, calcSoupProfile = FALSE)
    meta <- read.csv(file.path(d, paste0("obs_", slug, ".csv")), row.names = 1)
    cl   <- if ("seurat_clusters" %in% colnames(meta))
              as.character(meta[colnames(cnt), "seurat_clusters"]) else rep("c1", ncol(cnt))
    cl[is.na(cl)] <- "c1"
    soup <- setClusters(soup, setNames(cl, colnames(cnt)))
    um   <- read.csv(file.path(d, paste0("umap_", slug, ".csv")), row.names = 1)
    soup <- setDR(soup, as.matrix(um[colnames(cnt), c("UMAP_1", "UMAP_2")]))

    gl <- list(smc = SMC_GENES)
    p_smc  <- read_profile(file.path(EXP, "soup_profile_SMC.csv"), rownames(cnt))
    p_self <- read_profile(file.path(d, paste0("soup_profile_", slug, ".csv")), rownames(cnt))
    say("[profile] smc sum(est)=", round(sum(p_smc$est), 6),
        " self sum(est)=", round(sum(p_self$est), 6))

    ch_smc  <- setSoupProfile(soup, p_smc)
    ch_self <- setSoupProfile(soup, p_self)
    r_smc  <- estimate_rho(ch_smc,  gl)
    r_self <- estimate_rho(ch_self, gl)
    say("[rho] smc-profile  = ", round(r_smc$rho, 4),  " (non-expressing cells ", r_smc$ne, ")",
        ifelse(is.na(r_smc$err), "", paste0("  ERR: ", r_smc$err)))
    say("[rho] self-profile = ", round(r_self$rho, 4), " (non-expressing cells ", r_self$ne, ")",
        ifelse(is.na(r_self$err), "", paste0("  ERR: ", r_self$err)))

    ## ---------- Decision rule ----------
    ok <- function(x) !is.na(x) && x > RHO_LO && x < RHO_HI
    if (ok(r_smc$rho)) {
      mode <- "smc";  use <- r_smc$ch
    } else if (ok(r_self$rho)) {
      mode <- "self"; use <- r_self$ch
    } else {
      mode <- "-";    use <- NULL
    }
    status <- if (is.null(use)) "SKIPPED" else "OK"
    say("[decide] mode = ", mode, "  status = ", status)
    if (!is.na(r_smc$rho) && !is.na(r_self$rho) && abs(r_smc$rho - r_self$rho) > 0.20)
      say("[warn] rho differs by > 0.20 between the two soup profiles (",
          round(r_smc$rho, 3), " vs ", round(r_self$rho, 3),
          "): results are sensitive to the ambient profile composition")

    ## ---------- adjustCounts ----------
    if (status == "OK") {
      out <- tryCatch(adjustCounts(use, roundToInt = FALSE), error = function(e) {
        say("[adjust] failed: ", conditionMessage(e)); NULL })
      if (is.null(out)) {
        status <- "SKIPPED"; mode <- "-"
      } else {
        writeMM(out, file.path(od, "matrix_adjusted.mtx"))
        writeLines(rownames(out), file.path(od, "genes.tsv"))
        writeLines(colnames(out), file.path(od, "barcodes.tsv"))
        b <- Matrix::colSums(cnt); a <- Matrix::colSums(out)
        qc <- data.frame(barcode = names(b), total_before = as.numeric(b),
                         total_after = as.numeric(a[names(b)]))
        qc$frac_removed <- 1 - qc$total_after / qc$total_before
        write.csv(qc, file.path(od, "qc_delta.csv"), row.names = FALSE)
        gb <- Matrix::rowSums(cnt); ga <- Matrix::rowSums(out)
        gd <- data.frame(gene = names(gb), before = as.numeric(gb), after = as.numeric(ga[names(gb)]))
        gd$delta <- gd$after - gd$before
        gd$delta_pct <- 100 * gd$delta / gd$before
        write.csv(gd[order(gd$delta), ], file.path(od, "gene_delta.csv"), row.names = FALSE)
        say("[adjust] total UMI before ", sum(b), " after ", round(sum(a)),
            " ; median fraction removed = ", round(median(qc$frac_removed), 4))
        say("[adjust] top-10 most-depleted genes: ",
            paste(head(gd$gene[order(gd$delta)], 10), collapse = ", "))
        ## Key sanity check: the share of SMC markers in the new matrix must drop
        f_before <- sum(gb[SMC_GENES]) / sum(gb) * 100
        f_after  <- sum(ga[SMC_GENES]) / sum(ga) * 100
        say("[sanity] SMC 8 markers as % of total UMI: before ", round(f_before, 2),
            "% -> after ", round(f_after, 2), "%")
      }
    }

    if (status == "OK" && !is.null(use)) {
      write.csv(data.frame(barcode = colnames(cnt), rho = use$metaData$rho),
                file.path(od, "rho_per_cell.csv"), row.names = FALSE)
    }
    write.csv(data.frame(slug = slug, n_cells = ncol(cnt), mode = mode,
                         rho_smc_profile = r_smc$rho, rho_self_profile = r_self$rho,
                         status = status),
              file.path(od, "rho_summary.csv"), row.names = FALSE)

    summary_rows[[slug]] <- data.frame(
      slug = slug, n_cells = ncol(cnt), mode = mode,
      rho_smc_profile = round(r_smc$rho, 4), rho_self_profile = round(r_self$rho, 4),
      status = status)
  }, error = function(e) {
    say("[FATAL] ", conditionMessage(e))
    summary_rows[[slug]] <<- data.frame(slug = slug, n_cells = NA, mode = "-",
                                        rho_smc_profile = NA, rho_self_profile = NA,
                                        status = "ERROR")
  })
  close(con)
}

all <- do.call(rbind, summary_rows)
write.csv(all, file.path(OUTD, "rho_summary_all.csv"), row.names = FALSE)
cat("\n==================== Summary ====================\n")
print(all, row.names = FALSE)
cat("\nDONE_SOUPX_PER_LINEAGE\n")
