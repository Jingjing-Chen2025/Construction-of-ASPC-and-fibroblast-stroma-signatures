# ============================================================================
# readers.R — every reader returns a dgCMatrix (genes x cells)
# ============================================================================
read_triplet_10x <- function(barcodes, features, matrix_mtx, sample_name) {
  bc <- readLines(open_maybe_gz(barcodes))
  bc <- make.unique(sanitize_cell_names(bc))
  ft <- data.table::fread(features, header = FALSE)
  genes <- if (ncol(ft) >= 2) ft[[2]] else ft[[1]]
  genes <- sanitize_gene_names(genes)
  m <- as(as(Matrix::readMM(open_maybe_gz(matrix_mtx)), "CsparseMatrix"), "dgCMatrix")
  if (nrow(m) != length(genes) || ncol(m) != length(bc)) {
    if (ncol(m) == length(genes) && nrow(m) == length(bc)) {
      msg("  %s: matrix transposed; correcting.", sample_name)
      m <- Matrix::t(m)
    } else {
      stop(sprintf("Dimension mismatch for %s: matrix %dx%d, features %d, barcodes %d",
                   sample_name, nrow(m), ncol(m), length(genes), length(bc)))
    }
  }
  rownames(m) <- genes; colnames(m) <- bc
  finalize_matrix(m)
}

read_dir_10x <- function(dir_path, sample_name) {
  mtx <- list.files(dir_path, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  if (!length(mtx)) stop("No matrix.mtx in directory: ", dir_path)
  mat_dir <- dirname(mtx[1])
  bc <- list.files(mat_dir, pattern = "barcodes\\.tsv(\\.gz)?$", full.names = TRUE)
  ft <- list.files(mat_dir, pattern = "(features|genes)\\.tsv(\\.gz)?$", full.names = TRUE)
  if (!length(bc) || !length(ft)) stop("Incomplete 10x directory: ", mat_dir)
  read_triplet_10x(bc[1], ft[1], mtx[1], sample_name)
}

read_h5_10x <- function(path, sample_name) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) stop("hdf5r is required to read ", basename(path))
  m <- Seurat::Read10X_h5(path)
  if (is.list(m)) {
    nm <- names(m)
    pick <- if ("Gene Expression" %in% nm) "Gene Expression" else nm[1]
    msg("  %s: multi-modal h5, using '%s'", sample_name, pick)
    m <- m[[pick]]
  }
  finalize_matrix(m)
}

read_triplet_custom <- function(barcode_csv, genes_csv, counts_mtx, sample_name) {
  bc_df <- data.table::fread(cmd = paste("gzip -cdf", shQuote(barcode_csv)), data.table = FALSE)
  bc <- if (ncol(bc_df) == 1) bc_df[[1]] else bc_df[[ncol(bc_df)]]
  bc <- make.unique(sanitize_cell_names(bc))
  g_df <- data.table::fread(cmd = paste("gzip -cdf", shQuote(genes_csv)), data.table = FALSE)
  genes <- if (ncol(g_df) == 1) g_df[[1]] else g_df[[ncol(g_df)]]
  genes <- sanitize_gene_names(genes)
  m <- as(as(Matrix::readMM(open_maybe_gz(counts_mtx)), "CsparseMatrix"), "dgCMatrix")
  if (nrow(m) == length(genes) && ncol(m) == length(bc)) {
    rownames(m) <- genes; colnames(m) <- bc
  } else if (ncol(m) == length(genes) && nrow(m) == length(bc)) {
    msg("  %s: matrix transposed; correcting.", sample_name)
    m <- Matrix::t(m); rownames(m) <- genes; colnames(m) <- bc
  } else {
    stop(sprintf("Dimension mismatch for %s: matrix %dx%d, genes %d, barcodes %d",
                 sample_name, nrow(m), ncol(m), length(genes), length(bc)))
  }
  finalize_matrix(m)
}

read_dense_csv <- function(path, sample_name) {
  dt <- data.table::fread(cmd = paste("gzip -cdf", shQuote(path)), data.table = FALSE)
  if (!("V1" %in% colnames(dt))) colnames(dt)[1] <- "V1"
  cells <- make.unique(sanitize_cell_names(dt[["V1"]]))
  drop  <- intersect(c("V1", "CLUSTER", "Cluster", "cluster"), colnames(dt))
  keep  <- setdiff(colnames(dt), drop)
  genes <- sanitize_gene_names(keep)
  sp <- df_to_sparse_by_chunks(dt[, keep, drop = FALSE], label = sample_name)
  rm(dt); gc(verbose = FALSE)
  rownames(sp) <- cells
  colnames(sp) <- genes
  msg("  %s: dense_csv %d cells x %d genes -> transposing.", sample_name, nrow(sp), ncol(sp))
  finalize_matrix(Matrix::t(sp))
}

read_gene_cell_table <- function(path, sample_name) {
  dt <- data.table::fread(cmd = paste("gzip -cdf", shQuote(path)), data.table = FALSE)
  nms <- colnames(dt)
  if (!all(c("Gene_ID", "Symbol") %in% nms)) stop("Unrecognised gene_cell_exprs_table format: ", basename(path))
  genes <- sanitize_gene_names(dt[["Symbol"]])
  keep  <- setdiff(nms, c("Gene_ID", "Symbol"))
  sp <- df_to_sparse_by_chunks(dt[, keep, drop = FALSE], label = sample_name)
  rm(dt); gc(verbose = FALSE)
  rownames(sp) <- genes
  colnames(sp) <- make.unique(sanitize_cell_names(keep))
  finalize_matrix(sp)
}

read_matrix_flat <- function(path, sample_name) {
  sn <- sniff_table(path, sample_name)
  read_try <- function(sep, header_offset) {
    args <- list(cmd = paste("gzip -cdf", shQuote(path)), sep = sep,
                 data.table = FALSE, showProgress = FALSE)
    if (header_offset) {
      args$header <- FALSE
      d <- do.call(data.table::fread, args)
      hd <- as.character(unlist(d[1, ]))
      d <- d[-1, , drop = FALSE]
      hd <- hd[!is.na(hd) & hd != ""]
      colnames(d) <- c("GENE", hd)[seq_len(ncol(d))]
      d
    } else {
      args$header <- TRUE
      do.call(data.table::fread, args)
    }
  }
  dt <- NULL
  attempts <- list(
    list(sep = sn$sep, off = sn$offset), list(sep = sn$sep, off = !sn$offset),
    list(sep = "\t", off = FALSE), list(sep = ",", off = FALSE),
    list(sep = " ", off = TRUE), list(sep = " ", off = FALSE)
  )
  for (a in attempts) {
    dt <- tryCatch(read_try(a$sep, a$off), error = function(e) NULL)
    if (!is.null(dt) && ncol(dt) >= 2 && nrow(dt) >= 2) break
    dt <- NULL
  }
  if (is.null(dt)) stop("Could not parse flat matrix: ", basename(path))
  ids  <- sanitize_gene_names(dt[[1]])
  keep <- colnames(dt)[-1]
  sp <- df_to_sparse_by_chunks(dt[, keep, drop = FALSE], label = sample_name)
  rm(dt); gc(verbose = FALSE)
  rownames(sp) <- ids
  colnames(sp) <- make.unique(sanitize_cell_names(keep))
  if (looks_like_barcodes(rownames(sp)) && !looks_like_barcodes(colnames(sp))) {
    msg("  %s: row IDs look like barcodes -> transposing.", sample_name)
    sp <- Matrix::t(sp)
  }
  finalize_matrix(sp)
}

read_labelled_csv <- function(path, sample_name) {
  sn <- sniff_table(path, sample_name)
  dt <- data.table::fread(cmd = paste("gzip -cdf", shQuote(path)), sep = sn$sep,
                          data.table = FALSE, showProgress = FALSE)
  if (ncol(dt) < 2) stop("Labelled CSV has < 2 columns: ", basename(path))
  if (!nzchar(colnames(dt)[1]) || colnames(dt)[1] %in% c("V1", "")) colnames(dt)[1] <- "ROWID"
  rowid <- make.unique(sanitize_cell_names(dt[[1]]))
  is_num <- vapply(dt[-1], function(v) {
    if (is.numeric(v)) return(TRUE)
    mean(!is.na(suppressWarnings(as.numeric(head(v, 500))))) > 0.9
  }, logical(1))
  keep <- names(is_num)[is_num]
  if (!length(keep)) stop("No numeric columns found in ", basename(path), " (annotation table, not expression)")
  sp <- df_to_sparse_by_chunks(dt[, keep, drop = FALSE], label = sample_name)
  rm(dt); gc(verbose = FALSE)
  rownames(sp) <- rowid
  colnames(sp) <- sanitize_gene_names(keep)
  row_bc <- looks_like_barcodes(rownames(sp)); col_bc <- looks_like_barcodes(colnames(sp))
  if (row_bc && !col_bc) {
    sp <- Matrix::t(sp)
  } else if (!row_bc && !col_bc && ncol(sp) > nrow(sp)) {
    msg("  %s: more columns than rows -> transposing (shape heuristic).", sample_name)
    sp <- Matrix::t(sp)
  }
  finalize_matrix(sp)
}

read_archive <- function(path, sample_name, kind = c("tar", "zip")) {
  kind <- match.arg(kind)
  exdir <- file.path(tempdir(), paste0("ex_", gsub("[^A-Za-z0-9_]", "_", sample_name), "_", as.integer(Sys.time())))
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(exdir, recursive = TRUE, force = TRUE), add = TRUE)
  if (kind == "tar") utils::untar(path, exdir = exdir) else utils::unzip(path, exdir = exdir)
  read_dir_10x(exdir, sample_name)
}
