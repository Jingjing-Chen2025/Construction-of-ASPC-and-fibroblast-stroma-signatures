# ============================================================================
# clustering.R — PCA -> (Harmony) -> SNN -> Louvain at several resolutions -> UMAP
# Shared by the coarse tier (per dataset, integrated by sample) and the stromal
# tier (all datasets, integrated by dataset and sample).
# ============================================================================
run_harmony_safe <- function(seu, group_by, npcs) {
  group_by <- intersect(group_by, colnames(seu@meta.data))
  if (!length(group_by)) return(list(seu = seu, reduction = "pca"))
  n_groups <- length(unique(do.call(paste, c(seu@meta.data[group_by], sep = "|"))))
  if (n_groups < 2) return(list(seu = seu, reduction = "pca"))
  if (!requireNamespace("harmony", quietly = TRUE)) {
    warn_msg("harmony is not installed; clustering on the unintegrated PCA (install.packages('harmony'))")
    return(list(seu = seu, reduction = "pca"))
  }
  res <- tryCatch(
    harmony::RunHarmony(seu, group.by.vars = group_by, dims.use = 1:npcs,
                        reduction.save = "harmony", verbose = FALSE),
    error = function(e) { warn_msg("RunHarmony failed (%s); using unintegrated PCA", conditionMessage(e)); NULL })
  if (is.null(res)) return(list(seu = seu, reduction = "pca"))
  msg("  Harmony integration by %s (%d groups)", paste(group_by, collapse = " + "), n_groups)
  list(seu = res, reduction = "harmony")
}

embed_and_cluster <- function(seu, npcs, resolution, extra_resolutions = NULL, integrate_by = NULL, label = "") {
  seu <- RunPCA(seu, npcs = npcs, verbose = FALSE)
  red <- "pca"
  if (!is.null(integrate_by)) {
    h <- run_harmony_safe(seu, integrate_by, npcs); seu <- h$seu; red <- h$reduction
  }
  seu <- FindNeighbors(seu, reduction = red, dims = 1:npcs, verbose = FALSE)
  # extra resolutions first; the requested one last so that seurat_clusters holds it
  for (r in unique(c(setdiff(extra_resolutions, resolution), resolution))) {
    seu <- FindClusters(seu, resolution = r, verbose = FALSE)
  }
  seu <- RunUMAP(seu, reduction = red, dims = 1:npcs, verbose = FALSE)
  Idents(seu) <- seu$seurat_clusters
  seu@misc$reduction_used <- red
  seu@misc$resolution_used <- resolution
  msg("%sClusters: %d (resolution %g, %d dims, reduction '%s')", label,
      length(levels(seu$seurat_clusters)), resolution, npcs, red)
  seu
}

# clusters found at every stored resolution (stability check, clustree-style)
resolution_table <- function(seu) {
  cols <- grep("_snn_res\\.", colnames(seu@meta.data), value = TRUE)
  if (!length(cols)) return(NULL)
  data.frame(resolution = as.numeric(sub(".*_snn_res\\.", "", cols)),
             n_clusters = vapply(cols, function(cn) length(unique(seu@meta.data[[cn]])), integer(1)),
             used = as.numeric(sub(".*_snn_res\\.", "", cols)) == seu@misc$resolution_used)
}
