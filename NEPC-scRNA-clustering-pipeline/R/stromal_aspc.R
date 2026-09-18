# ============================================================================
# stromal_aspc.R — PART 1b: integrated stromal tier for ASPC identification
# ============================================================================
# 1. From every dataset's coarse result take the mesenchymal cells (clusters
#    labelled fibroblast / ASPC / smooth muscle / pericyte by the marker panel,
#    or fibroblast / smooth muscle / tissue stem cell by SingleR).
# 2. Map mouse genes to human orthologs (babelgene), merge all datasets, and
#    integrate with Harmony by dataset and sample.
# 3. Recluster at cfg$stromal_resolution and label clusters with the stromal
#    panels (ASPC, committed preadipocyte, matrix fibroblast, myCAF, iCAF,
#    smooth muscle, pericyte, Schwann, cycling, contaminants).
# 4. Call ASPC per cell by consensus: stromal cluster label + per-cell ASPC
#    score (UCell) above threshold and above the competing panels (+ optional
#    reference label, e.g. Emont 2022 adipose atlas).
# 5. Write the integrated object, per-cluster composition across datasets,
#    per-dataset population counts (used by Part 2) and per-dataset markers
#    (used by Part 3).
# ============================================================================

map_to_human_symbols <- function(mat, species) {
  if (!identical(species, "mouse")) return(mat)
  if (isTRUE(cfg$map_mouse_orthologs) && requireNamespace("babelgene", quietly = TRUE)) {
    orth <- tryCatch(babelgene::orthologs(genes = rownames(mat), species = "mouse", human = FALSE),
                     error = function(e) NULL)
    if (!is.null(orth) && nrow(orth)) {
      orth <- orth[!duplicated(orth$symbol), c("symbol", "human_symbol")]
      hs   <- orth$human_symbol[match(rownames(mat), orth$symbol)]
      keep <- !is.na(hs) & hs != ""
      msg("  ortholog mapping (babelgene): %d of %d mouse genes -> human symbols", sum(keep), nrow(mat))
      mat <- mat[keep, , drop = FALSE]
      rownames(mat) <- hs[keep]
      return(collapse_duplicate_features(mat))
    }
    warn_msg("  babelgene mapping failed; falling back to toupper()")
  } else if (isTRUE(cfg$map_mouse_orthologs)) {
    warn_msg("  babelgene not installed; mouse genes mapped by toupper() (install.packages('babelgene'))")
  }
  rownames(mat) <- toupper(rownames(mat))
  collapse_duplicate_features(mat)
}

extract_stromal_cells <- function(ds) {   # returns list(obj = Seurat, by_sample = data.frame)
  rds <- file.path(dataset_out_dir(ds), "nepc_seurat_annotated.rds")
  if (!file.exists(rds)) { warn_msg("  %s: no coarse result (%s)", ds$name, rds); return(NULL) }
  seu <- readRDS(rds)
  md  <- seu@meta.data
  pick <- (md$population %in% cfg$stromal_populations) |
          (if ("population_panel" %in% colnames(md)) md$population_panel %in% cfg$stromal_populations else FALSE) |
          (if ("singler_population" %in% colnames(md)) md$singler_population %in% cfg$stromal_singler_populations else FALSE)
  pick[is.na(pick)] <- FALSE
  cells <- colnames(seu)[pick]
  msg("  %s: %d of %d cells are mesenchymal", ds$name, length(cells), ncol(seu))
  by_sample <- data.frame(dataset = ds$name, group = ds$group, sample = md$sample, stromal = pick) %>%
    group_by(dataset, group, sample) %>%
    summarise(n_total_cells = n(), n_stromal_cells = sum(stromal), .groups = "drop")
  if (length(cells) < cfg$min_cells_population) { warn_msg("  %s: too few stromal cells (%d); skipped", ds$name, length(cells)); return(NULL) }
  if (length(cells) > cfg$stromal_max_cells_per_dataset) {
    cells <- sample(cells, cfg$stromal_max_cells_per_dataset)
    msg("  %s: subsampled to %d stromal cells", ds$name, length(cells))
  }
  counts <- get_counts_layer(seu)[, cells, drop = FALSE]
  counts <- map_to_human_symbols(counts, ds$species)
  keep_cols <- intersect(c("sample", "dataset", "species", "seurat_clusters", "population", "population_panel",
                           "singler_population", "singler_label"), colnames(md))
  meta <- md[cells, keep_cols, drop = FALSE]
  names(meta)[names(meta) == "seurat_clusters"] <- "coarse_cluster"
  names(meta)[names(meta) == "population"]      <- "coarse_population"
  meta$coarse_cluster <- as.character(meta$coarse_cluster)
  meta$group <- ds$group
  used <- as.data.frame(table(meta$sample), stringsAsFactors = FALSE); names(used) <- c("sample", "n_stromal_used")
  by_sample <- by_sample %>% left_join(used, by = "sample") %>% mutate(n_stromal_used = ifelse(is.na(n_stromal_used), 0L, n_stromal_used))
  rm(seu, md); gc(verbose = FALSE)
  list(obj = CreateSeuratObject(counts = counts, project = ds$name, meta.data = meta, min.cells = 0, min.features = 0),
       by_sample = by_sample)
}

# per-cell ASPC decision from the score matrix and the cluster label
call_aspc_cells <- function(scores, cluster_population) {
  aspc <- scores[, "ASPC_adipose_progenitor"]
  competitors <- intersect(c("myCAF", "iCAF", "Smooth_muscle", "Pericyte", "Committed_preadipocyte",
                             "Fibroblast_matrix", "Schwann_neural"), colnames(scores))
  best_other <- if (length(competitors)) apply(scores[, competitors, drop = FALSE], 1, max) else -Inf
  in_cluster <- cluster_population == "ASPC_adipose_progenitor"
  thr <- cfg$aspc_cell_min_score
  if (is.null(thr)) {
    bg <- aspc[!in_cluster]
    thr <- if (length(bg) >= 10) stats::median(bg) + stats::mad(bg) else stats::median(aspc)
  }
  pass <- aspc >= thr & aspc >= best_other
  call <- ifelse(in_cluster & pass, "ASPC_consensus",
          ifelse(in_cluster, "ASPC_cluster_only",
          ifelse(pass, "ASPC_cell_only", "non_ASPC")))
  list(call = call, pass = pass, threshold = thr)
}

run_stromal_tier <- function(group) {
  hdr("PART 1b — integrated stromal tier (ASPC): %s datasets", group)
  out_dir <- file.path(OUT_CROSS, paste0("stromal_", group))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  objs <- list(); extraction <- list()
  for (ds in group_datasets(group)) {
    o <- tryCatch(extract_stromal_cells(ds), error = function(e) {
      warn_msg("  %s: stromal extraction failed: %s", ds$name, clean_msg(conditionMessage(e))); NULL })
    if (!is.null(o)) { objs[[ds$name]] <- o$obj; extraction[[ds$name]] <- o$by_sample }
  }
  if (length(extraction)) write.csv(bind_rows(extraction), file.path(out_dir, "nepc_stromal_extraction_by_sample.csv"), row.names = FALSE)
  if (length(objs) < 2) stop("fewer than two ", group, " datasets contribute stromal cells; nothing to integrate")

  seu <- merge(objs[[1]], y = objs[-1]); rm(objs); gc(verbose = FALSE)
  seu <- join_layers_safe(seu)
  msg("Stromal object: %d genes x %d cells from %d datasets", nrow(seu), ncol(seu), length(unique(seu$dataset)))

  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = cfg$stromal_nfeatures, verbose = FALSE)
  vf  <- VariableFeatures(seu)
  if (!length(vf)) stop("no variable features in the stromal object")
  seu <- ScaleData(seu, features = vf, verbose = FALSE)
  npcs <- min(cfg$stromal_npcs, ncol(seu) - 1, length(vf) - 1)
  integ <- if (identical(cfg$stromal_integrate, "harmony")) c("dataset", "sample") else NULL
  seu <- embed_and_cluster(seu, npcs, cfg$stromal_resolution, NULL, integ, label = "stromal ")

  # ---- cluster labels from the stromal panels ----
  ps  <- cluster_panel_scores(seu, STROMAL_MARKERS)
  ann <- assign_population(ps$scores, cfg$stromal_ann_min_score, cfg$stromal_ann_min_margin)
  ann$population[grepl("^Contaminant", ann$population)] <- "Contaminant"
  write.csv(cbind(cluster = rownames(ps$scores), as.data.frame(round(ps$scores, 4))),
            file.path(out_dir, "nepc_stromal_panel_scores_by_cluster.csv"), row.names = FALSE)

  seu$stromal_cluster <- as.character(seu$seurat_clusters)
  pop_map <- setNames(ann$population, ann$cluster)
  seu$stromal_population <- unname(pop_map[seu$stromal_cluster])

  # ---- per-cell scores + ASPC consensus call ----
  sc <- score_cells(seu, STROMAL_MARKERS)
  sc <- sc[colnames(seu), , drop = FALSE]
  for (p in colnames(sc)) seu[[paste0("score_", p)]] <- unname(sc[, p])
  cl <- call_aspc_cells(sc, seu$stromal_population)
  seu$aspc_cell_score_pass <- cl$pass
  seu$aspc_call <- cl$call
  msg("Per-cell ASPC threshold (%s): %.3f; calls: %s", attr(sc, "method"), cl$threshold,
      paste(names(table(cl$call)), table(cl$call), sep = "=", collapse = ", "))

  # ---- optional reference vote (e.g. Emont 2022 adipose atlas) ----
  ref_lab <- run_singler_custom(seu, cfg$stromal_reference_rds, cfg$stromal_reference_label)
  seu$stromal_reference_label <- NA_character_
  if (!is.null(ref_lab)) {
    seu$stromal_reference_label[match(names(ref_lab), colnames(seu))] <- unname(ref_lab)
    ref_major <- seu@meta.data %>% filter(!is.na(stromal_reference_label)) %>%
      count(stromal_cluster, stromal_reference_label, name = "n_lab") %>% group_by(stromal_cluster) %>%
      mutate(frac = n_lab / sum(n_lab)) %>% slice_max(order_by = n_lab, n = 1, with_ties = FALSE) %>% ungroup() %>%
      transmute(cluster = stromal_cluster, reference_majority = stromal_reference_label, reference_fraction = round(frac, 3))
    ann <- ann %>% left_join(ref_major, by = "cluster")
  } else {
    ann$reference_majority <- NA_character_; ann$reference_fraction <- NA_real_
  }

  # ---- composition of every stromal cluster across datasets ----
  md <- seu@meta.data
  comp <- md %>% count(stromal_cluster, dataset, name = "n") %>%
    pivot_wider(names_from = dataset, values_from = n, values_fill = 0)
  ds_cols <- setdiff(colnames(comp), "stromal_cluster")
  comp$n_datasets_present <- rowSums(comp[, ds_cols, drop = FALSE] >= cfg$min_cells_population_stromal)
  sizes <- md %>% count(stromal_cluster, name = "n_cells") %>%
    mutate(pct_cells = round(100 * n_cells / sum(n_cells), 2))
  aspc_frac <- md %>% group_by(stromal_cluster) %>%
    summarise(frac_cells_pass_aspc = round(mean(aspc_cell_score_pass), 3), .groups = "drop")
  ann_out <- ann %>% left_join(sizes, by = c("cluster" = "stromal_cluster")) %>%
    left_join(aspc_frac, by = c("cluster" = "stromal_cluster")) %>%
    left_join(comp, by = c("cluster" = "stromal_cluster")) %>%
    arrange(suppressWarnings(as.numeric(cluster)))
  write.csv(ann_out, file.path(out_dir, "nepc_stromal_cluster_annotation.csv"), row.names = FALSE)

  # ---- per-dataset population table (Part 2 reads this) ----
  by_ds <- md %>% filter(!stromal_population %in% c("Contaminant", "Unassigned")) %>%
    count(dataset, group, species, stromal_population, name = "n_cells") %>%
    group_by(dataset) %>% mutate(pct_of_stromal = round(100 * n_cells / sum(n_cells), 2)) %>% ungroup() %>%
    transmute(dataset, group, species, cluster = paste0("S_", stromal_population), population = stromal_population,
              n_cells, pct_of_stromal, panel_score = NA_real_, tier = "stromal")
  write.csv(by_ds, file.path(out_dir, "nepc_stromal_population_by_dataset.csv"), row.names = FALSE)

  aspc_by_ds <- md %>% count(dataset, species, aspc_call, name = "n_cells") %>%
    pivot_wider(names_from = aspc_call, values_from = n_cells, values_fill = 0)
  write.csv(aspc_by_ds, file.path(out_dir, "nepc_aspc_calls_by_dataset.csv"), row.names = FALSE)
  write.csv(md, file.path(out_dir, "nepc_stromal_cell_metadata.csv"), row.names = TRUE)

  # ---- markers: integrated (by stromal population) and per dataset (Part 3 consensus) ----
  Idents(seu) <- seu$stromal_population
  mk_int <- tryCatch(FindAllMarkers(seu, only.pos = TRUE, min.pct = cfg$marker_min_pct,
                                    logfc.threshold = cfg$marker_logfc, verbose = FALSE),
                     error = function(e) { warn_msg("stromal FindAllMarkers failed: %s", conditionMessage(e)); NULL })
  if (!is.null(mk_int) && nrow(mk_int)) {
    mk_int <- mk_int[mk_int$p_val_adj < cfg$marker_p_adj, , drop = FALSE]
    mk_int$population <- as.character(mk_int$cluster)
    write.csv(mk_int, file.path(out_dir, "nepc_stromal_markers_integrated.csv"), row.names = FALSE)
  }
  per_ds <- list()
  pops <- setdiff(unique(seu$stromal_population), c("Contaminant", "Unassigned", "Cycling"))
  for (d in unique(seu$dataset)) {
    sub <- subset(seu, cells = colnames(seu)[seu$dataset == d])
    for (p in pops) {
      n_in <- sum(sub$stromal_population == p); n_out <- sum(sub$stromal_population != p)
      if (n_in < cfg$min_cells_population_stromal || n_out < cfg$min_cells_population_stromal) next
      m <- tryCatch(FindMarkers(sub, ident.1 = p, group.by = "stromal_population", only.pos = TRUE,
                                min.pct = cfg$marker_min_pct, logfc.threshold = cfg$marker_logfc, verbose = FALSE),
                    error = function(e) NULL)
      if (is.null(m) || !nrow(m)) next
      m <- m[m$p_val_adj < cfg$marker_p_adj, , drop = FALSE]
      if (!nrow(m)) next
      per_ds[[length(per_ds) + 1]] <- data.frame(dataset = d, cluster = paste0("S_", p), population = p,
                                                 gene = rownames(m), gene_upper = toupper(rownames(m)),
                                                 avg_log2FC = m$avg_log2FC, p_val_adj = m$p_val_adj,
                                                 pct.1 = m$pct.1, pct.2 = m$pct.2, tier = "stromal")
    }
    rm(sub); gc(verbose = FALSE)
  }
  if (length(per_ds)) write.csv(bind_rows(per_ds), file.path(out_dir, "nepc_stromal_markers_by_dataset.csv"), row.names = FALSE)

  # ---- plots ----
  p1 <- DimPlot(seu, group.by = "stromal_cluster", label = TRUE, repel = TRUE) + ggtitle(paste0(group, " stromal tier - clusters"))
  p2 <- DimPlot(seu, group.by = "stromal_population", label = TRUE, repel = TRUE, label.size = 3) + ggtitle("Stromal tier - populations")
  p3 <- DimPlot(seu, group.by = "dataset") + ggtitle("Stromal tier - dataset")
  p4 <- DimPlot(seu, group.by = "aspc_call") + ggtitle("ASPC call")
  pdf(file.path(out_dir, "nepc_stromal_umap.pdf"), width = 18, height = 12); print((p1 + p2) / (p3 + p4)); dev.off()
  feats <- intersect(paste0("score_", c("ASPC_adipose_progenitor", "myCAF", "Smooth_muscle", "Pericyte", "Committed_preadipocyte", "iCAF")),
                     colnames(seu@meta.data))
  pdf(file.path(out_dir, "nepc_stromal_umap_scores.pdf"), width = 15, height = 10)
  print(FeaturePlot(seu, features = feats, ncol = 3)); dev.off()
  if (any(!is.na(seu$stromal_reference_label))) {
    pdf(file.path(out_dir, "nepc_stromal_umap_reference.pdf"), width = 10, height = 7)
    print(DimPlot(seu, group.by = "stromal_reference_label", label = TRUE, repel = TRUE, label.size = 3)); dev.off()
  }
  if (isTRUE(cfg$save_rds)) saveRDS(seu, file.path(out_dir, "nepc_stromal_integrated.rds"))

  cat("\nStromal clusters (cells, population, datasets present):\n")
  print(as.data.frame(ann_out[, c("cluster", "n_cells", "population", "panel_score", "frac_cells_pass_aspc", "n_datasets_present")]))
  msg("Stromal outputs written to %s", out_dir)
  invisible(ann_out)
}
