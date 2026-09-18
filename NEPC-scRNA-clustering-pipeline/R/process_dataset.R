# ============================================================================
# process_dataset.R — PART 1: merge, cluster, annotate one dataset
# ============================================================================
process_dataset <- function(ds) {
  hdr("DATASET: %s (%s, %s)", ds$name, ds$group %||% "?", ds$species)
  if (!dir.exists(ds$dir)) stop("Dataset directory not found: ", ds$dir)
  out_dir <- file.path(ds$dir, "nepc_clustering")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  rds_path <- file.path(out_dir, "nepc_seurat_annotated.rds")

  if (isTRUE(cfg$reuse_existing) && file.exists(rds_path) &&
      file.exists(file.path(out_dir, "nepc_cluster_annotation.csv"))) {
    msg("Existing outputs found in %s; skipping (cfg$reuse_existing = TRUE).", out_dir)
    ann <- read.csv(file.path(out_dir, "nepc_cluster_annotation.csv"))
    return(data.frame(dataset = ds$name, group = ds$group, species = ds$species, dir = out_dir, status = "reused",
                      n_samples = NA_integer_, n_cells = sum(ann$n_cells), n_clusters = nrow(ann)))
  }

  # ---- STEP 1: discover + load samples ----
  samples <- scan_dataset(ds$dir)
  if (!length(samples)) stop("No recognisable samples found in ", ds$dir)
  samples <- triage_duplicate_gsm(samples)
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
  integ <- if (identical(cfg$integrate_samples, "harmony")) "sample" else NULL
  res_used <- cfg$resolution_by_group[[ds$group %||% ""]] %||% cfg$resolution
  seu <- embed_and_cluster(seu, npcs, res_used, cfg$extra_resolutions, integ, label = "")
  rt <- resolution_table(seu)
  if (!is.null(rt)) write.csv(rt, file.path(out_dir, "nepc_clusters_by_resolution.csv"), row.names = FALSE)
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
  keys <- gene_keys(rownames(seu), ds$species)
  score_panels <- c(POPULATION_MARKERS, list(ASPC_specific = ASPC_SPECIFIC, Fibroblast_matrix = FIBROBLAST_MATRIX))
  cell_scores <- tryCatch(score_cells(seu, score_panels, keys = keys), error = function(e) {
    warn_msg("per-cell scoring failed: %s", conditionMessage(e)); NULL })
  if (!is.null(cell_scores)) {
    for (p in colnames(cell_scores)) seu[[paste0("score_", p)]] <- unname(cell_scores[colnames(seu), p])
    msg("  per-cell panel scores: %s (%d panels)", attr(cell_scores, "method"), ncol(cell_scores))
  }
  if (identical(cfg$cluster_label_method, "ucell") && !is.null(cell_scores)) {
    ps <- cluster_panel_scores_ucell(cell_scores, seu$seurat_clusters, POPULATION_MARKERS)
    write.csv(cbind(cluster = rownames(ps$cluster_means), as.data.frame(round(ps$cluster_means, 4))),
              file.path(out_dir, "nepc_panel_ucell_mean_by_cluster.csv"), row.names = FALSE)
  } else {
    ps <- cluster_panel_scores(seu, POPULATION_MARKERS, keys = keys)
  }
  if (!is.null(ps$frac_top)) {
    ann <- assign_population_ucell(ps)
    write.csv(cbind(cluster = rownames(ps$frac_top), as.data.frame(round(ps$frac_top, 3))),
              file.path(out_dir, "nepc_panel_best_fraction_by_cluster.csv"), row.names = FALSE)
  } else {
    ann <- assign_population(ps$scores); ann$label_cell_fraction <- NA_real_
  }
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
  seu$singler_label <- lab                                   # reference label as returned
  seu$singler_population <- harmonize_singler_label(lab)     # shared human/mouse vocabulary
  # optional prostate-specific reference (Song 2022 / Henry 2018 ...) as a third label
  ref_lab <- run_singler_custom(seu, cfg$prostate_reference_rds, cfg$prostate_reference_label)
  seu$reference_label <- NA_character_
  if (!is.null(ref_lab)) seu$reference_label[match(names(ref_lab), colnames(seu))] <- unname(ref_lab)
  ref_major <- if (!is.null(ref_lab)) seu@meta.data %>%
    mutate(cluster = as.character(seurat_clusters)) %>% filter(!is.na(reference_label)) %>%
    count(cluster, reference_label, name = "n_lab") %>% group_by(cluster) %>%
    mutate(frac = n_lab / sum(n_lab)) %>% slice_max(order_by = n_lab, n = 1, with_ties = FALSE) %>%
    ungroup() %>% transmute(cluster, reference_majority = reference_label, reference_fraction = round(frac, 3))
  else data.frame(cluster = character(0), reference_majority = character(0), reference_fraction = numeric(0))
  # ASPC-high cells: ASPC-specific score above threshold and above the fibroblast matrix score
  aspc_tbl <- if (!is.null(cell_scores) && all(c("ASPC_specific", "Fibroblast_matrix") %in% colnames(cell_scores))) {
    seu$aspc_high <- cell_scores[colnames(seu), "ASPC_specific"] >= cfg$aspc_ucell_min &
                     cell_scores[colnames(seu), "ASPC_specific"] > cell_scores[colnames(seu), "Fibroblast_matrix"]
    seu@meta.data %>% mutate(cluster = as.character(seurat_clusters)) %>% group_by(cluster) %>%
      summarise(frac_cells_aspc_high = round(mean(aspc_high), 3),
                mean_aspc_specific = round(mean(score_ASPC_specific), 4),
                mean_fibroblast_matrix = round(mean(score_Fibroblast_matrix), 4), .groups = "drop")
  } else data.frame(cluster = character(0), frac_cells_aspc_high = numeric(0), mean_aspc_specific = numeric(0), mean_fibroblast_matrix = numeric(0))
  sr_major <- seu@meta.data %>%
    mutate(cluster = as.character(seurat_clusters)) %>%
    filter(!is.na(singler_population)) %>%
    count(cluster, singler_population, name = "n_lab") %>%
    group_by(cluster) %>% mutate(frac = n_lab / sum(n_lab)) %>%
    slice_max(order_by = n_lab, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(cluster, singler_majority = singler_population, singler_fraction = round(frac, 3))

  # ---- STEP 5: assemble the cluster annotation table ----
  sizes <- seu@meta.data %>% mutate(cluster = as.character(seurat_clusters)) %>%
    count(cluster, name = "n_cells") %>% mutate(pct_cells = round(100 * n_cells / sum(n_cells), 2))
  top_mk <- if (!is.null(markers) && nrow(markers)) {
    markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>%
      slice_head(n = 10) %>% summarise(top_markers = paste(gene, collapse = ";"), .groups = "drop")
  } else data.frame(cluster = character(0), top_markers = character(0))

  use_singler <- identical(cfg$annotation_method, "singler")
  ann <- ann %>% rename(population_panel = population) %>%
    left_join(sizes, by = "cluster") %>% left_join(sr_major, by = "cluster") %>%
    left_join(top_mk, by = "cluster") %>% left_join(ref_major, by = "cluster") %>% left_join(aspc_tbl, by = "cluster") %>%
    mutate(dataset = ds$name, group = ds$group, species = ds$species,
           cluster_name = paste0("Cluster_", cluster),
           annotation_method = if (use_singler) "singler_majority" else "marker_panel",
           population = if (use_singler) ifelse(is.na(singler_majority), "Unassigned", singler_majority)
                        else population_panel,
           compartment_panel = compartment_of(population_panel),
           compartment_singler = compartment_of(singler_majority),
           compartment_agreement = ifelse(is.na(compartment_singler) | population_panel == "Unassigned", NA,
                                          compartment_panel == compartment_singler)) %>%
    select(dataset, group, species, cluster, cluster_name, n_cells, pct_cells, population, annotation_method,
           population_panel, label_cell_fraction, singler_majority, singler_fraction, compartment_agreement,
           reference_majority, reference_fraction, frac_cells_aspc_high, mean_aspc_specific, mean_fibroblast_matrix,
           best_panel, panel_score, runner_up, runner_up_score, margin, top_markers) %>%
    arrange(suppressWarnings(as.numeric(cluster)))
  # ---- coarse-tier ASPC rule ----
  ann$aspc_rule_applied <- FALSE
  if (isTRUE(cfg$aspc_coarse_rule)) {
    stromal_like <- ann$population %in% c("Fibroblast", "ASPC_adipose_progenitor", "Smooth_muscle_myofibroblast", "Pericyte", "Unassigned") |
                    ann$singler_majority %in% c("Fibroblast", "Tissue_stem_cell", "Smooth_muscle_cell")
    msc <- ann$singler_majority %in% "Tissue_stem_cell" & ann$singler_fraction >= 0.5
    genes_win <- !is.na(ann$mean_aspc_specific) & ann$mean_aspc_specific > ann$mean_fibroblast_matrix &
                 ann$frac_cells_aspc_high >= 0.5
    hit <- stromal_like & (msc | genes_win) & ann$population != "ASPC_adipose_progenitor" &
           !ann$population %in% c("Luminal_epithelial", "Basal_epithelial", "Club_Hillock_epithelial", "Neuroendocrine",
                                  "T_cell", "NK_cell", "B_cell", "Plasma_cell", "Macrophage_myeloid", "Dendritic_cell",
                                  "Mast_cell", "Neutrophil", "Endothelial", "Lymphatic_endothelial", "Erythroid")
    if (any(hit)) {
      msg("  coarse ASPC rule: cluster(s) %s relabelled ASPC_adipose_progenitor (was %s)",
          paste(ann$cluster[hit], collapse = ","), paste(ann$population[hit], collapse = ","))
      ann$population[hit] <- "ASPC_adipose_progenitor"; ann$aspc_rule_applied[hit] <- TRUE
    }
  }
  n_dis <- sum(ann$compartment_agreement %in% FALSE)
  if (n_dis > 0) warn_msg("%d cluster(s) where the marker panel and SingleR disagree on the compartment; see nepc_cluster_annotation.csv", n_dis)
  if (use_singler && all(is.na(ann$singler_majority))) {
    warn_msg("SingleR produced no labels; clusters are Unassigned (set cfg$annotation_method = 'panel' to use marker panels)")
  }
  write.csv(ann, file.path(out_dir, "nepc_cluster_annotation.csv"), row.names = FALSE)

  pop_map <- setNames(ann$population, ann$cluster)
  seu$population <- unname(pop_map[as.character(seu$seurat_clusters)])
  panel_map <- setNames(ann$population_panel, ann$cluster)
  seu$population_panel <- unname(panel_map[as.character(seu$seurat_clusters)])
  seu$cluster_name <- paste0("Cluster_", seu$seurat_clusters)
  seu$cluster_population <- paste0("C", seu$seurat_clusters, ":", seu$population)
  seu$celltype <- seu$singler_label                       # batch-pipeline column name

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
  p1 <- DimPlot(seu, group.by = "cluster_name", label = TRUE, repel = TRUE, label.size = 4) +
    ggtitle(paste0(ds$name, " - Named Clusters")) + labs(color = "Cluster")
  p2 <- DimPlot(seu, group.by = "cluster_population", label = TRUE, repel = TRUE, label.size = 3) +
    ggtitle(paste0(ds$name, " - annotated populations")) + theme(legend.text = element_text(size = 7))
  p3 <- DimPlot(seu, group.by = "sample") + ggtitle(paste0(ds$name, " - Sample")) + labs(color = "Sample")
  pdf(file.path(out_dir, "nepc_umap_clusters_named.pdf"), width = 8, height = 6); print(p1); dev.off()
  pdf(file.path(out_dir, "nepc_umap_populations.pdf"), width = 11, height = 6); print(p2); dev.off()
  pdf(file.path(out_dir, "nepc_umap_sample.pdf"), width = 9, height = 6); print(p3); dev.off()
  if (any(!is.na(seu$celltype))) {
    p4 <- DimPlot(seu, group.by = "celltype", label = TRUE, repel = TRUE, label.size = 3) +
      ggtitle(paste0(ds$name, " - Cell Type (SingleR)")) + labs(color = "Cell Type")
    pdf(file.path(out_dir, "nepc_umap_celltypes.pdf"), width = 9, height = 6); print(p4); dev.off()
    pdf(file.path(out_dir, "nepc_umap_clusters_named_and_celltypes.pdf"), width = 14, height = 6); print(p1 + p4); dev.off()
  }
  pdf(file.path(out_dir, "nepc_umap_clusters_and_populations.pdf"), width = 18, height = 6); print(p1 + p2); dev.off()

  if (isTRUE(cfg$save_rds)) saveRDS(seu, rds_path)
  msg("Outputs written to %s", out_dir)
  print(ann[, c("cluster", "n_cells", "population", "singler_majority", "frac_cells_aspc_high", "aspc_rule_applied")])

  res <- data.frame(dataset = ds$name, group = ds$group, species = ds$species, dir = out_dir, status = "ok",
                    n_samples = length(unique(seu$sample)), n_cells = ncol(seu),
                    n_clusters = length(levels(seu$seurat_clusters)))
  rm(seu); gc(verbose = FALSE)
  res
}
