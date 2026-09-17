# ============================================================================
# discovery.R — sample discovery inside a dataset folder + per-sample loading/QC
# ============================================================================
# Files that are never expression matrices even when they are compressed.
ANNOTATION_ONLY_PATTERNS <- paste(
  "rename_cluster", "TCRresult", "TCR_result", "clonotype", "cell_annotation", "cellannotation",
  "_metadata", "metadata_", "barcode_annotation", sep = "|")

# A text file is only treated as an expression matrix when it has enough columns.
plausible_matrix_file <- function(path) {
  lines <- tryCatch(peek_lines(path, n = 2), error = function(e) character(0))
  if (length(lines) < 2) return(FALSE)
  sep <- detect_sep(lines)
  max(lengths(strsplit(lines, sep, fixed = TRUE))) >= cfg$min_fields_expression
}

# Strip the format suffix to obtain the sample id.
sample_id_from_file <- function(f) {
  sid <- f
  sid <- sub("\\.matrix\\.zip$|\\.zip$|\\.tar\\.gz$|\\.tgz$", "", sid)
  sid <- sub("\\.(txt|csv|tsv)\\.gz$", "", sid)
  sid <- sub("_filtered_feature_bc_matrix$|_raw_feature_bc_matrix$", "", sid)
  sid <- sub("_gene_cell_exprs_table$|_dense$|_dge$|\\.counts?$|_counts?$|_matrix$|_expression_matrix$", "", sid)
  sid
}

# ---- sample discovery --------------------------------------------------------
# One sample is EITHER
#   (i)  a 10x triplet  <sample>_barcodes.tsv.gz + <sample>_features.tsv.gz (or _genes.tsv.gz)
#                       + <sample>_matrix.mtx.gz, OR
#   (ii) one file       <sample>.txt.gz | <sample>.csv.gz | <sample>.zip | <sample>.matrix.zip
#                       | <sample>.tar.gz
# Every other file in the folder (uncompressed .csv/.txt, .rds, .pdf, analysis
# outputs of the companion scripts, sub-directories) is ignored.
scan_dataset <- function(dir_path) {
  files <- list.files(dir_path, full.names = FALSE, recursive = FALSE)
  files <- files[!grepl("^\\.", files)]
  files <- files[!dir.exists(file.path(dir_path, files))]
  samples <- list(); claimed <- character(0)
  add <- function(s) samples[[length(samples) + 1]] <<- s

  # (i) 10x triplets
  for (b in grep("_barcodes\\.tsv\\.gz$", files, value = TRUE)) {
    pfx <- sub("_barcodes\\.tsv\\.gz$", "", b)
    f <- files[files %in% paste0(pfx, c("_features.tsv.gz", "_genes.tsv.gz"))]
    m <- files[files %in% paste0(pfx, "_matrix.mtx.gz")]
    if (length(f) && length(m)) {
      add(list(sample_id = sample_id_from_file(pfx), type = "triplet_10x",
               barcodes = file.path(dir_path, b), features = file.path(dir_path, f[1]),
               matrix = file.path(dir_path, m[1])))
      claimed <- c(claimed, b, f[1], m[1])
    } else {
      warn_msg("  %s has no matching _features.tsv.gz / _matrix.mtx.gz; triplet ignored", b)
      claimed <- c(claimed, b)
    }
  }
  triplet_parts <- "_(barcodes|features|genes)\\.tsv\\.gz$|_matrix\\.mtx\\.gz$"

  # (ii) single-file samples
  for (f in setdiff(files, claimed)) {
    if (grepl(triplet_parts, f)) next                       # stray part of an incomplete triplet
    ty <- if (grepl("\\.tar\\.gz$|\\.tgz$", f))             "archive_tar"
      else if (grepl("\\.zip$", f))                         "archive_zip"
      else if (grepl("\\.(txt|csv|tsv)\\.gz$", f))          "text_gz"
      else NA_character_
    if (is.na(ty)) next                                     # not an accepted sample format
    fp <- file.path(dir_path, f)
    if (grepl(ANNOTATION_ONLY_PATTERNS, f, ignore.case = TRUE)) {
      msg("  Skipping %s (annotation table, not expression)", f); next
    }
    if (!is.null(cfg$sample_name_regex) && !grepl(cfg$sample_name_regex, f)) {
      msg("  Skipping %s (name does not match cfg$sample_name_regex '%s')", f, cfg$sample_name_regex); next
    }
    if (ty == "text_gz" && !plausible_matrix_file(fp)) {
      msg("  Skipping %s (fewer than %d columns; not an expression matrix)", f, cfg$min_fields_expression); next
    }
    add(list(sample_id = sample_id_from_file(f), type = ty, path = fp))
  }
  samples
}

load_sample <- function(s, ds) {
  m <- switch(s$type,
    triplet_10x = read_triplet_10x(s$barcodes, s$features, s$matrix, s$sample_id),
    text_gz     = read_text_matrix(s$path, s$sample_id),
    archive_zip = read_archive(s$path, s$sample_id, "zip"),
    archive_tar = read_archive(s$path, s$sample_id, "tar"),
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
  # case-insensitive so "MT-CO1", "mt-Co1" and "Mt-co1" are all recognised
  mt_genes <- grep("^MT-", rownames(seu), ignore.case = TRUE, value = TRUE)
  n_mt <- length(mt_genes)
  if (n_mt > 0) {
    # col.name= works in Seurat v4 and v5 (v5 returns a bare vector otherwise)
    seu <- PercentageFeatureSet(seu, features = mt_genes, col.name = "percent.mt")
  } else {
    seu$percent.mt <- rep(0, ncol(seu))
  }
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
