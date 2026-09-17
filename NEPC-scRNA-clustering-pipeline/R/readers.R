# ============================================================================
# readers.R — every reader returns a dgCMatrix (genes x cells)
# ============================================================================
# Accepted sample formats (discovered by scan_dataset() in discovery.R):
#   triplet_10x   <sample>_barcodes.tsv.gz + <sample>_features.tsv.gz (or _genes.tsv.gz)
#                 + <sample>_matrix.mtx.gz
#   text_gz       <sample>.txt.gz or <sample>.csv.gz — ONE delimited matrix per sample:
#                   * genes x cells (first column = gene symbol), or
#                   * cells x genes (first column = barcode; e.g. *_dense.csv.gz,
#                     annotation columns such as CLUSTER are dropped), or
#                   * Gene_ID + Symbol + one column per cell (*_gene_cell_exprs_table.txt.gz)
#                 Orientation is detected and the matrix is returned as genes x cells.
#   archive_zip   <sample>.zip / <sample>.matrix.zip
#   archive_tar   <sample>.tar.gz
#                 An archive wraps either a 10x directory (barcodes/features/matrix.mtx)
#                 or one text matrix as above.
# ============================================================================

# ---- 10x triplet ------------------------------------------------------------
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
  msg("  %s: 10x triplet, %d genes x %d cells", sample_name, nrow(m), ncol(m))
  finalize_matrix(m)
}

# ---- 10x directory (inside an extracted archive) ----------------------------
read_10x_dir <- function(dir_path, sample_name) {
  mtx <- list.files(dir_path, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  mtx <- mtx[!grepl("__MACOSX", mtx, fixed = TRUE)]
  if (!length(mtx)) stop("No matrix.mtx in directory: ", dir_path)
  mat_dir <- dirname(mtx[1])
  bc <- list.files(mat_dir, pattern = "barcodes\\.tsv(\\.gz)?$", full.names = TRUE)
  ft <- list.files(mat_dir, pattern = "(features|genes)\\.tsv(\\.gz)?$", full.names = TRUE)
  if (!length(bc) || !length(ft)) stop("Incomplete 10x directory: ", mat_dir)
  read_triplet_10x(bc[1], ft[1], mtx[1], sample_name)
}

# ---- one delimited text matrix ---------------------------------------------
# Column names that are cell annotations, never genes (cells x genes files).
NON_GENE_COLUMNS <- "^(cluster|clusters|cell_type|celltype|sample|batch|orig\\.ident|seurat_clusters|barcode|cell|cells|x)$"

read_text_matrix <- function(path, sample_name) {
  sn  <- sniff_table(path, sample_name)
  cmd <- paste("gzip -cdf", shQuote(path))

  read_try <- function(sep, header_offset) {
    args <- list(cmd = cmd, sep = sep, data.table = FALSE, showProgress = FALSE)
    if (header_offset) {                 # header has one field fewer than the data rows
      # fread would silently drop the short header line, so read it separately
      hd <- strsplit(peek_lines(path, n = 1), sep, fixed = TRUE)[[1]]
      hd <- gsub('"', "", hd, fixed = TRUE)
      hd <- hd[!is.na(hd) & hd != ""]
      args$header <- FALSE; args$skip <- 1
      d <- do.call(data.table::fread, args)
      if (ncol(d) == length(hd) + 1) colnames(d) <- c("ROWID", hd)
      else if (ncol(d) == length(hd)) colnames(d) <- hd
      else stop("header/data field mismatch")
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
  if (is.null(dt)) stop("Could not parse text matrix: ", basename(path))
  nms <- colnames(dt)

  # (a) Gene_ID + Symbol + one column per cell
  if (all(c("Gene_ID", "Symbol") %in% nms)) {
    genes <- sanitize_gene_names(dt[["Symbol"]])
    keep  <- setdiff(nms, c("Gene_ID", "Symbol"))
    sp <- df_to_sparse_by_chunks(dt[, keep, drop = FALSE], label = sample_name)
    rm(dt); gc(verbose = FALSE)
    rownames(sp) <- genes
    colnames(sp) <- make.unique(sanitize_cell_names(keep))
    msg("  %s: Gene_ID/Symbol table, %d genes x %d cells", sample_name, nrow(sp), ncol(sp))
    return(finalize_matrix(sp))
  }

  # (b) generic: first column = row IDs; keep numeric columns that are not annotations
  if (!nzchar(nms[1]) || nms[1] %in% c("V1", "")) colnames(dt)[1] <- "ROWID"
  rowid <- make.unique(sanitize_cell_names(dt[[1]]))
  is_num <- vapply(dt[-1], function(v) {
    if (is.numeric(v)) return(TRUE)
    mean(!is.na(suppressWarnings(as.numeric(head(v, 500))))) > 0.9
  }, logical(1))
  is_ann <- grepl(NON_GENE_COLUMNS, names(is_num), ignore.case = TRUE)
  drop <- names(is_num)[!is_num | is_ann]
  if (length(drop)) {
    msg("  %s: dropping %d annotation column(s): %s", sample_name, length(drop),
        paste(head(drop, 6), collapse = ", "))
  }
  keep <- names(is_num)[is_num & !is_ann]
  if (length(keep) < 2) stop("No numeric expression columns in ", basename(path))
  sp <- df_to_sparse_by_chunks(dt[, keep, drop = FALSE], label = sample_name)
  rm(dt); gc(verbose = FALSE)
  rownames(sp) <- rowid
  colnames(sp) <- sanitize_gene_names(keep)

  # orientation -> genes x cells
  row_bc <- looks_like_barcodes(rownames(sp)); col_bc <- looks_like_barcodes(colnames(sp))
  if (row_bc && !col_bc) {
    msg("  %s: rows are cells -> transposing to genes x cells.", sample_name)
    sp <- Matrix::t(sp)
  } else if (!row_bc && !col_bc && ncol(sp) > nrow(sp)) {
    msg("  %s: more columns than rows -> transposing (shape heuristic).", sample_name)
    sp <- Matrix::t(sp)
  }
  msg("  %s: text matrix, %d genes x %d cells", sample_name, nrow(sp), ncol(sp))
  finalize_matrix(sp)
}

# ---- archives (.zip / .matrix.zip / .tar.gz) --------------------------------
read_archive <- function(path, sample_name, kind = c("tar", "zip")) {
  kind  <- match.arg(kind)
  exdir <- file.path(tempdir(), paste0("ex_", gsub("[^A-Za-z0-9_]", "_", sample_name), "_", as.integer(Sys.time())))
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(exdir, recursive = TRUE, force = TRUE), add = TRUE)
  if (kind == "tar") utils::untar(path, exdir = exdir) else utils::unzip(path, exdir = exdir)

  inside <- list.files(exdir, recursive = TRUE, full.names = TRUE)
  inside <- inside[!grepl("__MACOSX|/\\._", inside)]
  if (any(grepl("matrix\\.mtx(\\.gz)?$", inside))) {
    msg("  %s: archive contains a 10x matrix directory", sample_name)
    return(read_10x_dir(exdir, sample_name))
  }
  txt <- inside[grepl("\\.(txt|csv|tsv)(\\.gz)?$", inside, ignore.case = TRUE)]
  if (!length(txt)) stop("Archive holds neither a 10x matrix directory nor a text matrix: ", basename(path))
  if (length(txt) > 1) {
    txt <- txt[which.max(file.info(txt)$size)]
    msg("  %s: several text files in archive; using the largest (%s)", sample_name, basename(txt))
  }
  read_text_matrix(txt, sample_name)
}
