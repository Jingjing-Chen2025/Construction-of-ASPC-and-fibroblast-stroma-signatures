# ============================================================================
# helpers.R — logging, small utilities, Seurat v4/v5 compatibility
# ============================================================================
SCRIPT_VERSION <- "nepc_scrna_pipeline_v2.0"

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

msg      <- function(...) cat("[INFO]", sprintf(...), "\n")
warn_msg <- function(...) cat("[WARN]", sprintf(...), "\n")
hdr      <- function(...) cat("\n========================================\n",
                              sprintf(...), "\n========================================\n", sep = "")

# ---- raise the macOS vector memory ceiling if the R build supports it -----
if (exists("mem.maxVSize")) {
  cur <- try(mem.maxVSize(), silent = TRUE)
  if (!inherits(cur, "try-error")) {
    try(mem.maxVSize(vsize = 64 * 1024), silent = TRUE)
  }
}


# ============================================================================
# GENERIC HELPERS
# ============================================================================
sanitize_gene_names <- function(x) {
  x <- str_trim(as.character(x))
  x <- gsub("﻿", "", x, fixed = TRUE)
  x <- gsub("\r", "", x, fixed = TRUE)
  x <- gsub('"', "", x, fixed = TRUE)
  x[x == "" | is.na(x)] <- "NA_FEATURE"
  x
}

sanitize_cell_names <- function(x) {
  x <- str_trim(as.character(x))
  x <- gsub("﻿", "", x, fixed = TRUE)
  x <- gsub("\r", "", x, fixed = TRUE)
  x <- gsub('"', "", x, fixed = TRUE)
  x[x == "" | is.na(x)] <- "NA_CELL"
  x
}

collapse_duplicate_features <- function(mat) {
  rn <- rownames(mat)
  if (is.null(rn) || !any(duplicated(rn))) return(mat)
  msg("  Collapsing %d duplicated feature name(s).", sum(duplicated(rn)))
  f <- factor(rn, levels = unique(rn))
  G <- Matrix::sparse.model.matrix(~ 0 + f)
  agg <- Matrix::t(G) %*% mat
  rownames(agg) <- levels(f)
  as(agg, "dgCMatrix")
}

finalize_matrix <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) mat <- as(as(mat, "CsparseMatrix"), "dgCMatrix")
  rownames(mat) <- sanitize_gene_names(rownames(mat))
  colnames(mat) <- make.unique(sanitize_cell_names(colnames(mat)))
  collapse_duplicate_features(mat)
}

open_maybe_gz <- function(path) if (grepl("\\.gz$", path)) gzfile(path) else file(path)

peek_lines <- function(path, n = 5) {
  con <- open_maybe_gz(path)
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

detect_sep <- function(lines) {
  if (!length(lines)) return("\t")
  counts <- c(
    "\t" = max(lengths(strsplit(lines, "\t", fixed = TRUE))) - 1,
    ","  = max(lengths(strsplit(lines, ",",  fixed = TRUE))) - 1,
    " "  = max(lengths(strsplit(lines, " ",  fixed = TRUE))) - 1
  )
  names(which.max(counts))
}

looks_like_barcodes <- function(x, n = 200) {
  x <- head(as.character(x), n)
  if (!length(x)) return(FALSE)
  frac_acgt <- mean(grepl("^[ACGTN]{8,}(-[0-9]+)?$", toupper(x)))
  frac_long <- mean(nchar(x) >= 14)
  frac_acgt > 0.5 || frac_long > 0.7
}

sniff_table <- function(path, sample_name) {
  lines <- peek_lines(path, n = 3)
  sep <- detect_sep(lines)
  hdr_n  <- length(strsplit(lines[1], sep, fixed = TRUE)[[1]])
  data_n <- if (length(lines) >= 2) length(strsplit(lines[2], sep, fixed = TRUE)[[1]]) else hdr_n
  list(sep = sep, header_fields = hdr_n, data_fields = data_n, offset = (data_n == hdr_n + 1))
}

df_to_sparse_by_chunks <- function(df, chunk = cfg$chunk_cols, label = "") {
  nc <- ncol(df)
  if (nc == 0) stop("No numeric columns to convert (", label, ")")
  starts <- seq(1, nc, by = chunk)
  blocks <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    a <- starts[i]; b <- min(a + chunk - 1, nc)
    blk <- as.matrix(df[, a:b, drop = FALSE])
    suppressWarnings(storage.mode(blk) <- "double")
    blk[is.na(blk)] <- 0
    blocks[[i]] <- as(as(blk, "CsparseMatrix"), "dgCMatrix")
    rm(blk)
    if (i %% 5 == 0 || i == length(starts)) gc(verbose = FALSE)
  }
  out <- if (length(blocks) == 1) blocks[[1]] else do.call(cbind, blocks)
  rm(blocks); gc(verbose = FALSE)
  out
}

# Seurat v4 / v5 compatibility -------------------------------------------------
# Seurat v5 keeps one counts/data layer per merged sample; join them so that
# FindAllMarkers, GetAssayData(layer=) and SingleR see one matrix. No-op on v4.
join_layers_safe <- function(obj) {
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    obj <- tryCatch(SeuratObject::JoinLayers(obj), error = function(e) {
      warn_msg("JoinLayers failed: %s", conditionMessage(e)); obj })
  }
  obj
}

seurat_v5 <- function() utils::packageVersion("SeuratObject") >= "5.0.0"

get_data_layer <- function(obj, assay = "RNA") {
  if (seurat_v5()) GetAssayData(obj, assay = assay, layer = "data")
  else GetAssayData(obj, assay = assay, slot = "data")
}

get_counts_layer <- function(obj, assay = "RNA") {
  if (seurat_v5()) GetAssayData(obj, assay = assay, layer = "counts")
  else GetAssayData(obj, assay = assay, slot = "counts")
}

# strip ANSI colour codes from rlang/cli error messages before logging
clean_msg <- function(x) gsub("\033\\[[0-9;]*m", "", x)

