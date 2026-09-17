# ============================================================================
# process_dataset.R — PART 1: merge, cluster, annotate one dataset
# ============================================================================
process_dataset <- function(ds) {
  hdr("DATASET: %s (%s)", ds$name, ds$species)
  if (!dir.exists(ds$dir)) stop("Dataset directory not found: ", ds$dir)
  out_dir <- file.path(ds$dir, "nepc_clustering")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  rds_path <- file.path(out_dir, "nepc_seurat_annotated.rds")

  if (isTRUE(cfg$reuse_existing) && file.exists(rds_path) &&
      file.exists(file.path(out_dir, "nepc_cluster_annotation.csv"))) {
    msg("Existing outputs found in %s; skipping (cfg$reuse_existing = TRUE).", out_dir)
    ann <- read.csv(file.path(out_dir, "nepc_cluster_annotation.csv"))
    return(data.frame(dataset = ds$name, species = ds$species, dir = out_dir, status = "reused",
                      n_samples = NA_integer_, n_cells = sum(ann$n_cells), n_clusters = nrow(ann)))
  }

  # ---- STEP 1: discover + load samples ----
  samples <- scan_dataset(ds$dir)
  if (!length(samples)) stop("No recognisable samples found in ", ds$dir)
  msg("Detected %d sample(s): %s", length(samples),
      paste(vapply(samples, function(s) sprintf("%s[%s]", s$sample_id, s$type), character(1)), collapse = ", "))

  seurat_list <- list(); load_log <- list()
  for (s in samples) {
    err <- NA_character_
    res <- tryCatch(load_sample(s, ds), error = function(e) {
      err <<- clean_msg(conditionMessage(e)); warn_msg("  FAILED %s: %s", s$sample_id, err); NULL })
    load_log[[length(load_log) + 1]] <- data.frame(
      dataset = ds$name, sample_id = s$sample_id, type = s$type,
      file = basename(s$path %||% s$matrix), status = if (is.null(res)) "failed" else "ok", error = err,
      n_cells = if (is.null(res)) NA_integer_ else ncol(res),
      n_genes = if (is.null(res)) NA_integer_ else nrow(res))
    if (!is.null(res)) seurat_list[[s$sample_id]] <- res
    gc(verbose = FALSE)
  }
  write.csv(bind_rows(load_log), file.path(out_dir, "nepc_sample_manifest.csv"), row.names = FALSE)
  if (!length(seurat_list)) stop("All samples failed to load in ", ds$name, " (see nepc_sample_manifest.csv)")

  # ---- STEP 2: merge + unsupervised clustering ----
  seu <- if (length(seurat_list) == 1) seurat_list[[1]] else
    merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = NULL)
  rm(seurat_list); gc(verbose = FALSE)
  msg("Merged: %d genes x %d cells from %d sample(s)", nrow(seu), ncol(seu), length(unique(seu$sample)))

  if (ncol(seu) < cfg$min_cells_dataset) {
    stop(sprintf("only %d cells after QC (< cfg$min_cells_dataset = %d); check nepc_sample_manifest.csv",
                 ncol(seu), cfg$min_cells_dataset))
  }
  seu <- join_layers_safe(seu)   # one counts layer (Seurat v5) before normalisation
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = cfg$nfeatures, verbose = FALSE)
  vf <- VariableFeatures(seu)
  if (!length(vf)) {
    warn_msg("vst found no variable features; retrying with selection.method = 'dispersion'")
    seu <- FindVariableFeatures(seu, selection.method = "dispersion", nfeatures = cfg$nfeatures, verbose = FALSE)
    vf <- VariableFeatures(seu)
  }
  if (!length(vf)) stop("No variable features could be determined")
  seu <- ScaleData(seu, features = vf, verbose = FALSE)
  npcs <- min(cfg$npcs, ncol(seu) - 1, length(vf) - 1)
  seu <- RunPCA(seu, npcs = npcs, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = cfg$resolution, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:npcs, min.dist = cfg$umap_min_dist, verbose = FALSE)
  seu <- join_layers_safe(seu)
  Idents(seu) <- seu$seurat_clusters
  msg("Clusters found: %d (resolution %g, %d PCs)", length(levels(seu$seurat_clusters)), cfg$resolution, npcs)
  gc(verbose = FALSE)

  # ---- STEP 3: markers ----
  markers <- NULL
  if (length(levels(seu$seurat_clusters)) >= 2) {
    markers <- tryCatch(
      FindAllMarkers(seu, only.pos = TRUE, min.pct = cfg$marker_min_pct,
                     logfc.threshold = cfg$marker_logfc, verbose = FALSE),
      error = function(e) { warn_msg("FindAllMarkers failed: %s", conditionMessage(e)); NULL })
  } else warn_msg("Only one cluster; skipping markers.")
  if (!is.null(markers) && nrow(markers)) {
    markers$cluster <- as.character(markers$cluster)
    markers$dataset <- ds$name
    markers$gene_upper <- toupper(markers$gene)
    p_col <- if ("p_val_adj" %in% colnames(markers)) "p_val_adj" else "p_val"
    markers <- markers[markers[[p_col]] < cfg$marker_p_adj, , drop = FALSE]
  }

  # ---- STEP 4: annotation (a) canonical panels ----
  ps  <- cluster_panel_scores(seu, POPULATION_MARKERS)
  ann <- assign_population(ps$scores)
  write.csv(cbind(cluster = rownames(ps$scores), as.data.frame(round(ps$scores, 4))),
            file.path(out_dir, "nepc_panel_scores_by_cluster.csv"), row.names = FALSE)
  msg("Panel genes detected: %s",
      paste(sprintf("%s=%d", names(ps$detected), ps$detected), collapse = ", "))

  # ---- STEP 4: annotation (b) SingleR ----
  lab <- rep(NA_character_, ncol(seu)); names(lab) <- colnames(seu)
  if (isTRUE(cfg$run_singler)) {
    sr <- run_singler(seu, ds$species)
    if (!is.null(sr)) lab[names(sr)] <- sr
  }
  seu$singler_label <- lab
  sr_major <- seu@meta.data %>%
    mutate(cluster = as.character(seurat_clusters)) %>%
    filter(!is.na(singler_label)) %>%
    count(cluster, singler_label, name = "n_lab") %>%
    group_by(cluster) %>% mutate(frac = n_lab / sum(n_lab)) %>%
    slice_max(order_by = n_lab, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(cluster, singler_majority = singler_label, singler_fraction = round(frac, 3))

  # ---- STEP 5: assemble the cluster annotation table ----
  sizes <- seu@meta.data %>% mutate(cluster = as.character(seurat_clusters)) %>%
    count(cluster, name = "n_cells") %>% mutate(pct_cells = round(100 * n_cells / sum(n_cells), 2))
  top_mk <- if (!is.null(markers) && nrow(markers)) {
    markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>%
      slice_head(n = 10) %>% summarise(top_markers = paste(gene, collapse = ";"), .groups = "drop")
  } else data.frame(cluster = character(0), top_markers = character(0))

  ann <- ann %>% left_join(sizes, by = "cluster") %>% left_join(sr_major, by = "cluster") %>%
    left_join(top_mk, by = "cluster") %>%
    mutate(dataset = ds$name, species = ds$species,
           population_label = ifelse(population == "Unassigned" & !is.na(singler_majority),
                                     paste0("Unassigned/", singler_majority), population)) %>%
    select(dataset, species, cluster, n_cells, pct_cells, population, population_label,
           best_panel, panel_score, runner_up, runner_up_score, margin,
           singler_majority, singler_fraction, top_markers) %>%
    arrange(suppressWarnings(as.numeric(cluster)))
  write.csv(ann, file.path(out_dir, "nepc_cluster_annotation.csv"), row.names = FALSE)

  pop_map <- setNames(ann$population, ann$cluster)
  seu$population <- unname(pop_map[as.character(seu$seurat_clusters)])
  seu$cluster_population <- paste0("C", seu$seurat_clusters, ":", seu$population)

  if (!is.null(markers) && nrow(markers)) {
    markers$population <- unname(pop_map[markers$cluster])
    write.csv(markers, file.path(out_dir, "nepc_cluster_markers_sig.csv"), row.names = FALSE)
    top_tbl <- markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>%
      slice_head(n = cfg$top_markers_report) %>% ungroup()
    write.csv(top_tbl, file.path(out_dir, "nepc_cluster_markers_top.csv"), row.names = FALSE)
  }

  pop_counts <- seu@meta.data %>% count(population, name = "n_cells") %>%
    mutate(pct_cells = round(100 * n_cells / sum(n_cells), 2), dataset = ds$name)
  write.csv(pop_counts, file.path(out_dir, "nepc_population_counts.csv"), row.names = FALSE)
  write.csv(seu@meta.data, file.path(out_dir, "nepc_cell_metadata.csv"), row.names = TRUE)

  # ---- STEP 6: plots ----
  p1 <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle(paste0(ds$name, " - clusters"))
  p2 <- DimPlot(seu, group.by = "cluster_population", label = TRUE, repel = TRUE, label.size = 3) +
    ggtitle(paste0(ds$name, " - annotated populations")) + theme(legend.text = element_text(size = 7))
  p3 <- DimPlot(seu, group.by = "sample") + ggtitle(paste0(ds$name, " - sample"))
  pdf(file.path(out_dir, "nepc_umap_clusters.pdf"), width = 8, height = 6); print(p1); dev.off()
  pdf(file.path(out_dir, "nepc_umap_populations.pdf"), width = 11, height = 6); print(p2); dev.off()
  pdf(file.path(out_dir, "nepc_umap_sample.pdf"), width = 9, height = 6); print(p3); dev.off()
  if (any(!is.na(seu$singler_label))) {
    p4 <- DimPlot(seu, group.by = "singler_label", label = TRUE, repel = TRUE, label.size = 3) +
      ggtitle(paste0(ds$name, " - SingleR"))
    pdf(file.path(out_dir, "nepc_umap_singler.pdf"), width = 10, height = 6); print(p4); dev.off()
  }
  pdf(file.path(out_dir, "nepc_umap_clusters_and_populations.pdf"), width = 18, height = 6); print(p1 + p2); dev.off()

  if (isTRUE(cfg$save_rds)) saveRDS(seu, rds_path)
  msg("Outputs written to %s", out_dir)
  print(ann[, c("cluster", "n_cells", "population", "panel_score", "singler_majority")])

  res <- data.frame(dataset = ds$name, species = ds$species, dir = out_dir, status = "ok",
                    n_samples = length(unique(seu$sample)), n_cells = ncol(seu),
                    n_clusters = length(levels(seu$seurat_clusters)))
  rm(seu); gc(verbose = FALSE)
  res
}
