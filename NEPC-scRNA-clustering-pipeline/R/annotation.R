# ============================================================================
# annotation.R — cluster annotation: canonical marker panels + SingleR
# ============================================================================
# cluster x panel score matrix: mean log-normalised expression per cluster,
# z-scored per gene across clusters, averaged over the panel genes.
cluster_panel_scores <- function(seu, panels, cluster_col = "seurat_clusters") {
  dat <- get_data_layer(seu)
  cl  <- as.character(seu@meta.data[[cluster_col]])
  lv  <- unique(cl)
  num <- suppressWarnings(as.numeric(lv))
  lv  <- if (!anyNA(num)) lv[order(num)] else sort(lv)
  cl  <- factor(cl, levels = lv)
  # cell x cluster indicator matrix (works for a single cluster too)
  G   <- Matrix::sparseMatrix(i = seq_along(cl), j = as.integer(cl), x = 1,
                              dims = c(length(cl), nlevels(cl)))
  n_per <- Matrix::colSums(G)

  genes_up    <- toupper(rownames(dat))
  panel_genes <- unique(unlist(panels))
  keep <- which(genes_up %in% panel_genes)
  keep <- keep[!duplicated(genes_up[keep])]
  if (!length(keep)) stop("None of the panel genes are present in the expression matrix.")

  means <- as.matrix(dat[keep, , drop = FALSE] %*% G)
  means <- sweep(means, 2, n_per, "/")
  rownames(means) <- genes_up[keep]
  colnames(means) <- levels(cl)

  z <- if (ncol(means) > 1) t(scale(t(means))) else means
  z[!is.finite(z)] <- 0

  scores <- vapply(names(panels), function(p) {
    g <- intersect(panels[[p]], rownames(z))
    if (length(g) < cfg$ann_min_genes) return(rep(NA_real_, ncol(z)))
    colMeans(z[g, , drop = FALSE])
  }, numeric(ncol(z)))
  scores <- matrix(scores, nrow = ncol(z), dimnames = list(colnames(z), names(panels)))

  detected <- vapply(names(panels), function(p) length(intersect(panels[[p]], rownames(z))), integer(1))
  list(scores = scores, cluster_means = means, detected = detected)
}

assign_population <- function(scores) {
  rows <- lapply(seq_len(nrow(scores)), function(i) {
    s <- scores[i, ]; s <- s[!is.na(s)]
    if (!length(s)) {
      return(data.frame(cluster = rownames(scores)[i], population = "Unassigned",
                        best_panel = NA_character_, panel_score = NA_real_,
                        runner_up = NA_character_, runner_up_score = NA_real_, margin = NA_real_))
    }
    o <- order(s, decreasing = TRUE)
    top <- unname(s[o[1]])
    second <- if (length(s) > 1) unname(s[o[2]]) else NA_real_
    margin <- if (is.na(second)) top else top - second
    lab <- if (top >= cfg$ann_min_score && margin >= cfg$ann_min_margin) names(s)[o[1]] else "Unassigned"
    data.frame(cluster = rownames(scores)[i], population = lab,
               best_panel = names(s)[o[1]], panel_score = top,
               runner_up = if (length(s) > 1) names(s)[o[2]] else NA_character_,
               runner_up_score = second, margin = margin)
  })
  bind_rows(rows)
}

# Map SingleR labels of either reference onto the shared vocabulary.
harmonize_singler_label <- function(x) {
  x <- as.character(x)
  out <- unname(SINGLER_LABEL_MAP[x])
  miss <- is.na(out) & !is.na(x)
  if (any(miss)) {
    y <- gsub("[ .-]+", "_", trimws(x[miss]))
    y <- sub("_?[Cc]ells$", "_cell", y)
    y <- sub("^(.*[^s])s$", "\\1", y)         # simple plural -> singular
    out[miss] <- y
  }
  out
}

get_singler_ref <- function(species) {
  if (identical(species, "mouse")) celldex::MouseRNAseqData() else celldex::HumanPrimaryCellAtlasData()
}

run_singler <- function(seu, species) {
  if (!requireNamespace("SingleR", quietly = TRUE) || !requireNamespace("celldex", quietly = TRUE) ||
      !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    warn_msg("SingleR / celldex / SingleCellExperiment not installed; skipping reference annotation.")
    return(NULL)
  }
  tryCatch({
    obj <- seu
    if (ncol(obj) > cfg$singler_max_cells) {
      msg("  Downsampling to %d cells for SingleR.", cfg$singler_max_cells)
      obj <- subset(obj, cells = sample(colnames(obj), cfg$singler_max_cells))
    }
    sce <- suppressWarnings(Seurat::as.SingleCellExperiment(obj))
    ref <- get_singler_ref(species)
    msg("  SingleR reference: %s", if (species == "mouse") "MouseRNAseqData" else "HumanPrimaryCellAtlasData")
    sr  <- SingleR::SingleR(test = sce, ref = ref,
                            labels = SummarizedExperiment::colData(ref)[[cfg$singler_labels %||% "label.main"]],
                            assay.type.test = cfg$singler_assay %||% "counts")
    out <- setNames(sr$labels, colnames(obj))
    rm(sce, obj, sr); gc(verbose = FALSE)
    out
  }, error = function(e) { warn_msg("  SingleR failed: %s", conditionMessage(e)); NULL })
}

