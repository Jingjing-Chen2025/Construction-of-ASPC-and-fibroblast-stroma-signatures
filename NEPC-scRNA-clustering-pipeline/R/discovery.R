# ============================================================================
# discovery.R — sample discovery inside a dataset folder + per-sample loading/QC
# ============================================================================
ANNOTATION_ONLY_PATTERNS <- paste(
  "rename_cluster", "TCRresult", "TCR_result", "tcr", "cell_annotation", "cellannotation",
  "_metadata", "metadata_", "clonotype", "barcode_annotation", sep = "|")

PIPELINE_OUTPUT_PATTERNS <- paste0(
  "^(metadata_merged|sample_summary|cluster_counts|sample_manifest|batch_run_log|",
  "matrix_file_comparison|cluster\\.markers|seurat_merged|umap_|nepc_|07_|tau_|ratio_|target_markers)")

scan_dataset <- function(dir_path) {
  files <- list.files(dir_path, full.names = FALSE, recursive = FALSE)
  files <- files[!grepl("^\\.", files)]
  samples <- list(); claimed <- character(0)
  add <- function(s) samples[[length(samples) + 1]] <<- s

  # --- 10x triplets: <pfx>_barcodes.tsv(.gz) + _features/_genes + _matrix.mtx ---
  bc10x <- grep("_barcodes\\.tsv(\\.gz)?$", files, value = TRUE)
  for (b in bc10x) {
    pfx <- sub("_barcodes\\.tsv(\\.gz)?$", "", b)
    f <- grep(paste0("^", pfx, "_(features|genes)\\.tsv(\\.gz)?$"), files, value = TRUE)
    m <- grep(paste0("^", pfx, "_matrix\\.mtx(\\.gz)?$"), files, value = TRUE)
    if (length(f) && length(m)) {
      add(list(sample_id = pfx, type = "triplet_10x", barcodes = file.path(dir_path, b),
               features = file.path(dir_path, f[1]), matrix = file.path(dir_path, m[1])))
      claimed <- c(claimed, b, f[1], m[1])
    }
  }
  # --- custom triplets: .barcode.csv / .genes.csv / .counts.mtx ---
  bccus <- grep("\\.barcode\\.csv(\\.gz)?$", files, value = TRUE)
  for (b in bccus) {
    pfx <- sub("\\.barcode\\.csv(\\.gz)?$", "", b)
    g <- grep(paste0("^", pfx, "\\.genes\\.csv(\\.gz)?$"), files, value = TRUE)
    m <- grep(paste0("^", pfx, "\\.counts\\.mtx(\\.gz)?$"), files, value = TRUE)
    if (length(g) && length(m)) {
      add(list(sample_id = pfx, type = "triplet_custom", barcodes = file.path(dir_path, b),
               features = file.path(dir_path, g[1]), matrix = file.path(dir_path, m[1])))
      claimed <- c(claimed, b, g[1], m[1])
    }
  }
  # --- single files, archives and 10x sub-directories ---
  for (f in setdiff(files, claimed)) {
    fp <- file.path(dir_path, f)
    if (grepl(PIPELINE_OUTPUT_PATTERNS, f)) next
    if (dir.exists(fp)) {
      has_mtx <- length(list.files(fp, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE)) > 0
      if (has_mtx) add(list(sample_id = f, type = "dir_10x", path = fp))
      next
    }
    if (grepl("\\.(rds|rda|RData)(\\.gz)?$", f, ignore.case = TRUE)) next
    if (grepl("\\.(pdf|log|md|json|png|jpg|xlsx|xls|html|R|r|py)$", f)) next
    if (grepl(ANNOTATION_ONLY_PATTERNS, f, ignore.case = TRUE)) next

    ty <- if (grepl("\\.h5$", f))                                         "h5_10x"
    else if (grepl("_dense\\.csv\\.gz$", f))                              "dense_csv"
    else if (grepl("_gene_cell_exprs_table\\.txt\\.gz$", f))              "gene_cell_table"
    else if (grepl("(_dge|data\\.raw\\.matrix|data\\.matrix|counts?|matrix|expr)\\.(txt|tsv)(\\.gz)?$", f, ignore.case = TRUE)) "matrix_txt"
    else if (grepl("\\.count(s)?\\.csv(\\.gz)?$", f))                     "matrix_csv"
    else if (grepl("\\.tar\\.gz$|\\.tgz$", f))                            "archive_tar"
    else if (grepl("\\.zip$", f))                                         "archive_zip"
    else if (grepl("\\.csv(\\.gz)?$", f))                                 "labelled_csv"
    else NA_character_
    if (is.na(ty)) next

    sid <- f
    sid <- sub("_dense\\.csv\\.gz$", "", sid)
    sid <- sub("_gene_cell_exprs_table\\.txt\\.gz$", "", sid)
    sid <- sub("_filtered_feature_bc_matrix", "", sid)
    sid <- sub("_raw_feature_bc_matrix", "", sid)
    sid <- sub("\\.count(s)?\\.csv(\\.gz)?$", "", sid)
    sid <- sub("(_dge)?\\.(txt|tsv)(\\.gz)?$", "", sid)
    sid <- sub("\\.csv(\\.gz)?$", "", sid)
    sid <- sub("\\.tar\\.gz$|\\.tgz$|\\.zip$|\\.h5$", "", sid)
    add(list(sample_id = sid, type = ty, path = fp))
  }
  samples
}

load_sample <- function(s, ds) {
  m <- switch(s$type,
    triplet_10x     = read_triplet_10x(s$barcodes, s$features, s$matrix, s$sample_id),
    triplet_custom  = read_triplet_custom(s$barcodes, s$features, s$matrix, s$sample_id),
    dir_10x         = read_dir_10x(s$path, s$sample_id),
    h5_10x          = read_h5_10x(s$path, s$sample_id),
    dense_csv       = read_dense_csv(s$path, s$sample_id),
    gene_cell_table = read_gene_cell_table(s$path, s$sample_id),
    matrix_txt      = read_matrix_flat(s$path, s$sample_id),
    matrix_csv      = read_matrix_flat(s$path, s$sample_id),
    labelled_csv    = read_labelled_csv(s$path, s$sample_id),
    archive_tar     = read_archive(s$path, s$sample_id, "tar"),
    archive_zip     = read_archive(s$path, s$sample_id, "zip"),
    stop("Unknown sample type: ", s$type)
  )
  if (is.finite(cfg$max_cells_per_sample) && ncol(m) > cfg$max_cells_per_sample) {
    msg("  %s: subsampling %d -> %d cells.", s$sample_id, ncol(m), cfg$max_cells_per_sample)
    m <- m[, sample(colnames(m), cfg$max_cells_per_sample), drop = FALSE]
  }
  colnames(m) <- paste0(s$sample_id, "__", colnames(m))

  seu <- CreateSeuratObject(counts = m, project = s$sample_id,
                            min.cells = cfg$min_cells_gene, min.features = cfg$min_features_cell)
  rm(m); gc(verbose = FALSE)

  # ---- per-cell QC ----
  mt_pattern <- if (identical(ds$species, "mouse")) "^mt-" else "^MT-"
  n_mt <- sum(grepl(mt_pattern, rownames(seu)))
  seu$percent.mt <- if (n_mt > 0) PercentageFeatureSet(seu, pattern = mt_pattern)[, 1] else rep(0, ncol(seu))
  n_before <- ncol(seu)
  keep <- seu$percent.mt <= cfg$max_mito_pct & seu$nFeature_RNA <= cfg$max_features_cell
  if (any(!keep)) seu <- subset(seu, cells = colnames(seu)[keep])
  msg("  QC %s: %d -> %d cells (mito > %g%% or nFeature > %s removed; %d MT genes)",
      s$sample_id, n_before, ncol(seu), cfg$max_mito_pct, format(cfg$max_features_cell), n_mt)

  seu$sample  <- s$sample_id
  seu$dataset <- ds$name
  seu$species <- ds$species
  seu
}
