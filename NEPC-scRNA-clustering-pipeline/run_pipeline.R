# ============================================================================
# NEPC vs adenocarcinoma scRNA-seq: per-dataset clustering + annotation,
#     group-specific cell sets (populations recurring across datasets of a
#     group), cell composition per group, and correlation of the NEPC-specific
#     cell-set signatures with the NEPC (Beltran, custom UP) signature in bulk
#     RNA-seq from PRAD TCGA and prad_su2c_2019 (cBioPortal)
# ============================================================================
#
# INPUT — two dataset groups, one folder per dataset (config.R):
#   NEPC  (7):  <NEPC_ROOT>/GSE137829, GSE210358, GSE210358_TKO (mouse), GSE235036_TKO (mouse),
#               GSE264573, GSE292074, GSE296986_TKO (mouse)
#   Adeno (9):  <ADENO_ROOT>/GSE137829, GSE141445, GSE176031, GSE181294, GSE210358, GSE264573,
#               GSE268307, GSE292074, GSE296986 (mouse)
#
# SAMPLE FORMATS (one sample = one of these inside a dataset folder)
#   <sample>_barcodes.tsv.gz + <sample>_features.tsv.gz + <sample>_matrix.mtx.gz
#   <sample>.txt.gz | <sample>.csv.gz | <sample>.zip | <sample>.matrix.zip | <sample>.tar.gz
#   All other files are ignored.
#
# PART 1  coarse clustering + annotation per dataset;  PART 1b stromal tier per group (ASPC)
# PART 2  populations in > cfg$recurrence_min_datasets datasets of a group = group-specific
#         cell sets; mean cell composition per group (per sample, coarse and stromal tiers)
# PART 3  signatures of the NEPC-specific cell sets vs the NEPC signature in PRAD TCGA and
#         prad_su2c_2019 (cBioPortal)
#
# FILE LAYOUT
#   config.R                    paths, parameters, dataset list, signatures, marker panels
#   R/helpers.R                 logging + Seurat v4/v5 compatibility helpers
#   R/readers.R                 format-specific matrix readers
#   R/discovery.R               sample discovery + per-sample loading and QC
#   R/annotation.R              cluster annotation (marker panels + SingleR)
#   R/clustering.R              PCA -> Harmony -> clustering at several resolutions -> UMAP
#   R/process_dataset.R         PART 1: coarse tier, one dataset end to end
#   R/stromal_aspc.R            PART 1b: integrated stromal tier (ASPC) across datasets
#   R/recurrence.R              PART 2: collect Part 1 / 1b outputs
#   R/signatures_cbioportal.R   PART 3: signatures, cBioPortal, scoring, correlation
#   run_pipeline.R              this file: runs Parts 1-3
#
# USAGE
#   1. Edit config.R (NEPC_ROOT, datasets, cfg).
#   2. Rscript run_pipeline.R        (or source() it from RStudio)
#   See README.md for the method description and the list of outputs.
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)
options(stringsAsFactors = FALSE, future.globals.maxSize = 8 * 1024^3)

# ---- load configuration and modules ----------------------------------------
# The folder that holds THIS file is used for config.R and R/, whether the
# script is run with Rscript, source()d, or executed from the RStudio editor.
PIPELINE_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f[1])))
  for (fr in sys.frames()) {                       # source("…/run_pipeline.R")
    of <- fr$ofile
    if (!is.null(of) && nzchar(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))   # "Run" / "Source" button in RStudio
  }
  getwd()
})
if (!file.exists(file.path(PIPELINE_DIR, "config.R")) || !dir.exists(file.path(PIPELINE_DIR, "R"))) {
  stop("Cannot find config.R and R/ next to run_pipeline.R (looked in ", PIPELINE_DIR,
       "). Run with Rscript run_pipeline.R from the pipeline folder, or setwd() to it first.")
}
source(file.path(PIPELINE_DIR, "config.R"))
for (mod in c("helpers", "readers", "discovery", "clustering", "annotation", "process_dataset",
              "stromal_aspc", "recurrence", "composition", "signatures_cbioportal")) {
  source(file.path(PIPELINE_DIR, "R", paste0(mod, ".R")))
}
cat("\n[INFO] Pipeline version", SCRIPT_VERSION, "loaded from", PIPELINE_DIR, "\n")
cat("[INFO] Seurat", as.character(packageVersion("Seurat")), "| R", R.version$major, R.version$minor, "\n\n")
# guard against a stale copy of the modules (e.g. an older unzipped folder)
if (!exists("plausible_matrix_file") || is.null(cfg$min_cells_dataset)) {
  stop("The files under ", PIPELINE_DIR, " are an older version of the pipeline. ",
       "Replace the whole folder (config.R, R/, run_pipeline.R) with the current one and re-run.")
}

dir.create(OUT_CROSS, showWarnings = FALSE, recursive = TRUE)
run_log <- list()

if (isTRUE(cfg$run_part1)) {
  hdr("PART 1 — clustering + annotation per dataset")
  for (ds in datasets) {
    t0 <- Sys.time()
    res <- tryCatch(process_dataset(ds), error = function(e) {
      warn_msg("DATASET FAILED [%s]: %s", ds$name, clean_msg(conditionMessage(e)))
      data.frame(dataset = ds$name, species = ds$species, dir = ds$dir,
                 status = paste("failed:", clean_msg(conditionMessage(e))),
                 n_samples = NA_integer_, n_cells = NA_integer_, n_clusters = NA_integer_)
    })
    res$minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
    run_log[[length(run_log) + 1]] <- res
    gc(verbose = FALSE)
  }
  run_log_df <- bind_rows(run_log)
  print(run_log_df)
  write.csv(run_log_df, file.path(OUT_CROSS, "nepc_run_log.csv"), row.names = FALSE)
}

# ============================================================================
# PART 1b — integrated stromal tier (ASPC), once per group
# ============================================================================
if (isTRUE(cfg$run_part1) && isTRUE(cfg$run_stromal_tier)) {
  for (g in GROUPS) {
    tryCatch(run_stromal_tier(g), error = function(e) {
      warn_msg("STROMAL TIER FAILED [%s]: %s", g, clean_msg(conditionMessage(e))) })
  }
}

# ============================================================================
# PART 2 — group-specific cell sets + cell composition
# ============================================================================
rec <- list(); comp <- list()
if (isTRUE(cfg$run_part2) || isTRUE(cfg$run_part3)) {
  for (g in GROUPS) {
    rec[[g]] <- tryCatch(run_recurrence(g), error = function(e) {
      warn_msg("PART 2 FAILED [%s]: %s", g, clean_msg(conditionMessage(e))); NULL })
  }
  # flag cell sets that are recurrent in the other group as well
  for (g in names(rec)) {
    if (is.null(rec[[g]])) next
    other <- setdiff(names(rec), g)
    other_sets <- unlist(lapply(other, function(o) if (!is.null(rec[[o]])) paste(rec[[o]]$recurrent$tier, rec[[o]]$recurrent$population) else character(0)))
    r <- rec[[g]]$recurrent
    r$also_recurrent_in_other_group <- paste(r$tier, r$population) %in% other_sets
    rec[[g]]$recurrent <- r
    write.csv(r, file.path(OUT_CROSS, paste0(g, "_specific_cell_sets.csv")), row.names = FALSE)
  }
  for (g in GROUPS) {
    if (is.null(rec[[g]])) next
    comp[[g]] <- tryCatch(run_composition(g, rec[[g]]$recurrent), error = function(e) {
      warn_msg("COMPOSITION FAILED [%s]: %s", g, clean_msg(conditionMessage(e))); NULL })
  }
  if (length(comp) >= 2) tryCatch(compare_composition(comp), error = function(e) {
    warn_msg("COMPOSITION COMPARISON FAILED: %s", clean_msg(conditionMessage(e))) })
}

# ============================================================================
# PART 3 — cell-set signatures vs NEPC signature in PRAD TCGA and SU2C 2019
# ============================================================================
if (isTRUE(cfg$run_part3)) {
  for (g in intersect(cfg$correlate_groups, names(rec))) {
    if (is.null(rec[[g]])) next
    tryCatch(run_part3(g, rec[[g]]$recurrent, read_all_markers(g)), error = function(e) {
      warn_msg("PART 3 FAILED [%s]: %s", g, clean_msg(conditionMessage(e))) })
  }
}

hdr("PIPELINE COMPLETE — cross-dataset outputs in %s", OUT_CROSS)
