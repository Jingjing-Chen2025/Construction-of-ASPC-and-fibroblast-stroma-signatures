# ============================================================================
# report.R — compact ASPC report from the result folders (no Seurat needed)
# ============================================================================
# Summarises, for every dataset and group: run status, coarse clusters called
# ASPC or ASPC-rich, mesenchymal-stem-cell clusters by SingleR, the stromal
# tier (present? ASPC clusters, ASPC calls per dataset, top ASPC markers), the
# group-specific cell sets, ASPC composition and the Part 3 correlations.
# Run at the end of run_pipeline.R, or standalone: Rscript report_aspc.R
# ============================================================================
.rd <- function(f, ...) if (file.exists(f)) tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, ...), error = function(e) NULL) else NULL

write_aspc_report <- function(out_file = file.path(OUT_CROSS, "nepc_ASPC_report.txt")) {
  L <- character(0); add <- function(...) L <<- c(L, sprintf(...))
  tbl <- function(df, cols = NULL, n = 50) {
    if (is.null(df) || !nrow(df)) { add("  (none)"); return(invisible()) }
    if (!is.null(cols)) df <- df[, intersect(cols, colnames(df)), drop = FALSE]
    txt <- capture.output(print(head(as.data.frame(df), n), row.names = FALSE)); L <<- c(L, paste0("  ", txt))
  }
  add("ASPC REPORT  %s  (%s)", format(Sys.time(), "%Y-%m-%d %H:%M"), SCRIPT_VERSION)
  add("results folder: %s", OUT_CROSS)

  add("\n== Run log =="); tbl(.rd(file.path(OUT_CROSS, "nepc_run_log.csv")), c("dataset", "group", "status", "n_samples", "n_cells", "n_clusters", "minutes"))

  for (g in GROUPS) {
    add("\n================ %s ================", g)
    for (ds in group_datasets(g)) {
      ann <- .rd(file.path(dataset_out_dir(ds), "nepc_cluster_annotation.csv"))
      if (is.null(ann)) { add("\n-- %s: no coarse result", ds$name); next }
      add("\n-- %s: %d clusters, %d cells; populations: %s", ds$name, nrow(ann), sum(ann$n_cells),
          paste(sort(unique(ann$population)), collapse = ", "))
      aspc <- ann[ann$population == "ASPC_adipose_progenitor" | (if ("aspc_rule_applied" %in% names(ann)) ann$aspc_rule_applied %in% TRUE else FALSE), ]
      add("   coarse ASPC clusters: %s", if (nrow(aspc)) paste(sprintf("C%s (%d cells)", aspc$cluster, aspc$n_cells), collapse = ", ") else "none")
      if ("frac_cells_aspc_high" %in% names(ann)) {
        rich <- ann[!is.na(ann$frac_cells_aspc_high) & ann$frac_cells_aspc_high >= 0.3, ]
        add("   ASPC-rich clusters (>= 30%% ASPC-high cells): %s",
            if (nrow(rich)) paste(sprintf("C%s %s %.0f%%", rich$cluster, rich$population, 100 * rich$frac_cells_aspc_high), collapse = "; ") else "none")
      }
      msc <- ann[ann$singler_majority %in% c("Tissue_stem_cell", "Stem_cell"), ]
      add("   SingleR mesenchymal-stem-cell clusters: %s",
          if (nrow(msc)) paste(sprintf("C%s %s (%d cells, %.0f%%)", msc$cluster, msc$population, msc$n_cells, 100 * msc$singler_fraction), collapse = "; ") else "none")
    }
    sd <- file.path(OUT_CROSS, paste0("stromal_", g))
    add("\n-- stromal tier %s: %s", g, if (dir.exists(sd)) "folder present" else "NOT RUN (folder missing)")
    if (dir.exists(sd)) {
      sa <- .rd(file.path(sd, "nepc_stromal_cluster_annotation.csv"))
      if (is.null(sa)) add("   nepc_stromal_cluster_annotation.csv missing -> the stromal tier failed; see the log for 'STROMAL TIER FAILED'")
      else {
        add("   stromal clusters:"); tbl(sa, c("cluster", "n_cells", "population", "panel_score", "frac_cells_pass_aspc", "n_datasets_present", "reference_majority"))
        add("   ASPC calls per dataset:"); tbl(.rd(file.path(sd, "nepc_aspc_calls_by_dataset.csv")))
        mk <- .rd(file.path(sd, "nepc_stromal_markers_integrated.csv"))
        if (!is.null(mk) && "population" %in% names(mk)) {
          mk <- mk[mk$population == "ASPC_adipose_progenitor", ]; mk <- mk[order(-mk$avg_log2FC), ]
          add("   top ASPC markers (integrated): %s", if (nrow(mk)) paste(head(mk$gene, 20), collapse = ", ") else "none")
        }
      }
    }
    add("\n-- %s-specific cell sets:", g); tbl(.rd(file.path(OUT_CROSS, paste0(g, "_specific_cell_sets.csv"))), c("tier", "population", "n_datasets", "n_clusters", "total_cells", "also_recurrent_in_other_group"))
    cm <- .rd(file.path(OUT_CROSS, paste0("nepc_cell_composition_mean_", g, ".csv")))
    add("-- %s composition of ASPC / fibroblast populations:", g)
    if (!is.null(cm)) tbl(cm[grepl("ASPC|Fibroblast|myCAF|iCAF|Tissue_stem", cm$population), ], c("tier", "population", "n_samples", "n_datasets_present", "mean_pct", "sd_pct", "mean_pct_of_all", "in_specific_cell_set")) else add("  (none)")
  }
  for (g in GROUPS) {
    cs <- .rd(file.path(OUT_CROSS, paste0("nepc_nepc_correlation_summary_", g, ".csv")))
    if (!is.null(cs)) { add("\n== Part 3 (%s): correlation with the NEPC signature (Pearson) ==", g); tbl(cs[cs$method == "pearson", ], c("study_id", "tier", "population", "signature_type", "r", "p", "p_adj_BH", "n_samples")) }
  }
  writeLines(L, out_file)
  cat(L, sep = "\n"); cat("\n[INFO] report written to", out_file, "\n")
  invisible(L)
}
