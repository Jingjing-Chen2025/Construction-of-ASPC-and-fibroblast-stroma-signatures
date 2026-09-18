# ============================================================================
# NEPC scRNA-seq: per-dataset clustering + cell-population annotation,
#     recurrent populations across datasets, and correlation of their
#     signatures with the NEPC (Beltran, custom UP) signature in bulk
#     RNA-seq from PRAD TCGA and prad_su2c_2019 (cBioPortal)
# ============================================================================
#
# INPUT  — seven NEPC dataset folders (one folder = one dataset):
#   /Volumes/Jingjing_Chen/NEPC scRNA dataset/GSE137829         human
#   /Volumes/Jingjing_Chen/NEPC scRNA dataset/GSE210358         human
#   /Volumes/Jingjing_Chen/NEPC scRNA dataset/GSE210358_TKO     mouse
#   /Volumes/Jingjing_Chen/NEPC scRNA dataset/GSE235036_TKO     mouse
#   /Volumes/Jingjing_Chen/NEPC scRNA dataset/GSE264573         human
#   /Volumes/Jingjing_Chen/NEPC scRNA dataset/GSE292074         human
#   /Volumes/Jingjing_Chen/NEPC scRNA dataset/GSE296986_TKO     mouse
#
# FILE LAYOUT
#   config.R                    paths, parameters, dataset list, signatures, marker panels
#   R/helpers.R                 logging + Seurat v4/v5 compatibility helpers
#   R/readers.R                 format-specific matrix readers
#   R/discovery.R               sample discovery + per-sample loading and QC
#   R/annotation.R              cluster annotation (marker panels + SingleR)
#   R/clustering.R              PCA -> Harmony -> clustering at several resolutions -> UMAP
#   R/process_dataset.R         PART 1: coarse tier, one dataset end to end
#   R/stromal_aspc.R            PART 1b: integrated stromal tier (ASPC) across datasets
#   R/recurrence.R              PART 2: collect Part 1 / 1b outputs
#   R/signatures_cbioportal.R   PART 3: signatures, cBioPortal, scoring, correlation
#   run_pipeline.R              this file: runs Parts 1-3
#
# USAGE
#   1. Edit config.R (NEPC_ROOT, datasets, cfg).
#   2. Rscript run_pipeline.R        (or source() it from RStudio)
#   See README.md for the method description and the list of outputs.
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)
options(stringsAsFactors = FALSE, future.globals.maxSize = 8 * 1024^3)

# ---- load configuration and modules ----------------------------------------
# The folder that holds THIS file is used for config.R and R/, whether the
# script is run with Rscript, source()d, or executed from the RStudio editor.
PIPELINE_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f[1])))
  for (fr in sys.frames()) {                       # source("…/run_pipeline.R")
    of <- fr$ofile
    if (!is.null(of) && nzchar(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))   # "Run" / "Source" button in RStudio
  }
  getwd()
})
if (!file.exists(file.path(PIPELINE_DIR, "config.R")) || !dir.exists(file.path(PIPELINE_DIR, "R"))) {
  stop("Cannot find config.R and R/ next to run_pipeline.R (looked in ", PIPELINE_DIR,
       "). Run with Rscript run_pipeline.R from the pipeline folder, or setwd() to it first.")
}
source(file.path(PIPELINE_DIR, "config.R"))
for (mod in c("helpers", "readers", "discovery", "clustering", "annotation", "process_dataset",
              "stromal_aspc", "recurrence", "signatures_cbioportal")) {
  source(file.path(PIPELINE_DIR, "R", paste0(mod, ".R")))
}
cat("\n[INFO] Pipeline version", SCRIPT_VERSION, "loaded from", PIPELINE_DIR, "\n")
cat("[INFO] Seurat", as.character(packageVersion("Seurat")), "| R", R.version$major, R.version$minor, "\n\n")
# guard against a stale copy of the modules (e.g. an older unzipped folder)
if (!exists("plausible_matrix_file") || is.null(cfg$min_cells_dataset)) {
  stop("The files under ", PIPELINE_DIR, " are an older version of the pipeline. ",
       "Replace the whole folder (config.R, R/, run_pipeline.R) with the current one and re-run.")
}

dir.create(OUT_CROSS, showWarnings = FALSE, recursive = TRUE)
run_log <- list()

if (isTRUE(cfg$run_part1)) {
  hdr("PART 1 — clustering + annotation per dataset")
  for (ds in datasets) {
    t0 <- Sys.time()
    res <- tryCatch(process_dataset(ds), error = function(e) {
      warn_msg("DATASET FAILED [%s]: %s", ds$name, clean_msg(conditionMessage(e)))
      data.frame(dataset = ds$name, species = ds$species, dir = ds$dir,
                 status = paste("failed:", clean_msg(conditionMessage(e))),
                 n_samples = NA_integer_, n_cells = NA_integer_, n_clusters = NA_integer_)
    })
    res$minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
    run_log[[length(run_log) + 1]] <- res
    gc(verbose = FALSE)
  }
  run_log_df <- bind_rows(run_log)
  print(run_log_df)
  write.csv(run_log_df, file.path(OUT_CROSS, "nepc_run_log.csv"), row.names = FALSE)
}

# ============================================================================
# PART 1b — integrated stromal tier (ASPC) across all datasets
# ============================================================================
if (isTRUE(cfg$run_part1) && isTRUE(cfg$run_stromal_tier)) {
  stromal_res <- tryCatch(run_stromal_tier(), error = function(e) {
    warn_msg("STROMAL TIER FAILED: %s", clean_msg(conditionMessage(e))); NULL })
}

# ============================================================================
# PART 2 — recurrent populations across datasets
# ============================================================================
recurrent <- NULL
if (isTRUE(cfg$run_part2) || isTRUE(cfg$run_part3)) {
  hdr("PART 2 — recurrent populations (present in > %d %s)",
      cfg$recurrence_min_datasets, cfg$recurrence_count_by)
  all_ann <- read_all_annotations()
  if (!nrow(all_ann)) stop("No Part 1 annotation tables found; run Part 1 first.")
  write.csv(all_ann, file.path(OUT_CROSS, "nepc_all_cluster_annotations.csv"), row.names = FALSE)

  all_ann <- all_ann %>% filter(n_cells >= cfg$min_cells_population)
  recurrence <- all_ann %>%
    filter(!population %in% c("Unassigned", "Contaminant")) %>%
    group_by(tier, population) %>%
    summarise(n_datasets = n_distinct(dataset),
              n_clusters = n(),
              total_cells = sum(n_cells),
              datasets = paste(sort(unique(dataset)), collapse = ";"),
              clusters = paste(paste0(dataset, ":C", cluster), collapse = ";"),
              mean_panel_score = round(mean(panel_score, na.rm = TRUE), 3),
              .groups = "drop") %>%
    mutate(n_datasets_total = length(unique(all_ann$dataset)),
           count_used = if (cfg$recurrence_count_by == "clusters") n_clusters else n_datasets,
           recurrent = count_used > cfg$recurrence_min_datasets) %>%
    arrange(desc(n_datasets), desc(n_clusters))
  write.csv(recurrence, file.path(OUT_CROSS, "nepc_population_recurrence.csv"), row.names = FALSE)

  presence <- all_ann %>% filter(!population %in% c("Unassigned", "Contaminant")) %>%
    count(tier, population, dataset, name = "n_clusters") %>%
    pivot_wider(names_from = dataset, values_from = n_clusters, values_fill = 0)
  write.csv(presence, file.path(OUT_CROSS, "nepc_population_presence_matrix.csv"), row.names = FALSE)

  recurrent <- recurrence %>% filter(recurrent)
  write.csv(recurrent, file.path(OUT_CROSS, "nepc_recurrent_populations.csv"), row.names = FALSE)

  cat("\nPopulations appearing in more than", cfg$recurrence_min_datasets, cfg$recurrence_count_by, ":\n")
  if (nrow(recurrent)) {
    print(as.data.frame(recurrent[, c("tier", "population", "n_datasets", "n_clusters", "total_cells", "datasets")]))
  } else cat("  (none)\n")
}

# ============================================================================
# PART 3 — signatures of recurrent populations vs NEPC signature in bulk
# ============================================================================
if (isTRUE(cfg$run_part3)) {
  hdr("PART 3 — recurrent-population signatures vs NEPC signature (cBioPortal)")
  if (is.null(recurrent) || !nrow(recurrent)) {
    warn_msg("No recurrent populations; Part 3 has nothing to correlate.")
  } else {
    all_markers <- read_all_markers()
    signatures  <- build_signatures(recurrent, all_markers)
    if (!nrow(signatures)) stop("No signatures could be built for the recurrent populations.")
    write.csv(signatures, file.path(OUT_CROSS, "nepc_recurrent_population_signatures.csv"), row.names = FALSE)
    sig_n <- signatures %>% filter(used_for_scoring) %>% count(tier, population, signature_type, name = "n_genes")
    msg("Signatures (genes used for scoring): %s",
        paste(sprintf("%s:%s/%s=%d", sig_n$tier, sig_n$population, sig_n$signature_type, sig_n$n_genes), collapse = ", "))

    genes_needed <- unique(c(NEPC_BELTRAN_CUSTOM_UP, signatures$gene[signatures$used_for_scoring]))
    sig_sets <- signatures %>% filter(used_for_scoring) %>%
      group_by(tier, population, signature_type) %>% summarise(genes = list(gene), .groups = "drop")

    summary_rows <- list(); gene_rows <- list(); profiles_used <- list(); plots <- list()

    for (key in names(cfg$cbio_studies)) {
      st <- cfg$cbio_studies[[key]]
      hdr("Bulk cohort: %s (%s)", st$label, st$study_id)
      bulk <- tryCatch(get_bulk_matrix(key, st, genes_needed), error = function(e) {
        warn_msg("Could not obtain expression for %s: %s", st$study_id, conditionMessage(e)); NULL })
      if (is.null(bulk)) next

      logmat <- to_log2(bulk$mat)
      msg("  %d genes x %d samples (of %d requested genes)", nrow(logmat), ncol(logmat), length(genes_needed))
      write.csv(cbind(gene = rownames(logmat), as.data.frame(round(logmat, 4))),
                file.path(OUT_CROSS, sprintf("nepc_bulk_log2_expression_%s.csv", st$study_id)), row.names = FALSE)
      profiles_used[[key]] <- data.frame(study_id = st$study_id, label = st$label, profile_id = bulk$profile_id,
                                         profile_name = bulk$profile_name, sample_list = bulk$sample_list,
                                         n_samples = ncol(logmat), n_genes = nrow(logmat))

      nepc <- score_signature(logmat, NEPC_BELTRAN_CUSTOM_UP)
      msg("  NEPC signature: %d/%d genes present", length(nepc$genes), length(NEPC_BELTRAN_CUSTOM_UP))
      scores <- data.frame(sample_id = colnames(logmat), NEPC_score = unname(nepc$score))

      for (i in seq_len(nrow(sig_sets))) {
        pop <- sig_sets$population[i]; typ <- sig_sets$signature_type[i]; tr <- sig_sets$tier[i]
        gs  <- sig_sets$genes[[i]]
        sc  <- score_signature(logmat, gs)
        col <- paste0(tr, ":", pop, "__", typ)
        scores[[col]] <- unname(sc$score)
        for (m in c("pearson", "spearman")) {
          cp <- cor_pair(sc$score, nepc$score, m)
          summary_rows[[length(summary_rows) + 1]] <- data.frame(
            study_id = st$study_id, cohort = st$label, tier = tr, population = pop, signature_type = typ,
            n_genes_signature = length(gs), n_genes_present = length(sc$genes),
            method = m, r = unname(cp["r"]), p = unname(cp["p"]), n_samples = unname(cp["n"]))
        }
        for (g in sc$genes) {
          cp <- cor_pair(logmat[g, ], nepc$score, "spearman")
          gene_rows[[length(gene_rows) + 1]] <- data.frame(
            study_id = st$study_id, tier = tr, population = pop, signature_type = typ, gene = g,
            rho = unname(cp["r"]), p = unname(cp["p"]), n_samples = unname(cp["n"]))
        }
      }
      write.csv(scores, file.path(OUT_CROSS, sprintf("nepc_bulk_signature_scores_%s.csv", st$study_id)), row.names = FALSE)

      # scatter plots: population score vs NEPC score
      long <- scores %>% pivot_longer(-c(sample_id, NEPC_score), names_to = "signature", values_to = "score") %>%
        filter(is.finite(score), is.finite(NEPC_score))
      if (nrow(long)) {
        p <- ggplot(long, aes(x = NEPC_score, y = score)) +
          geom_point(alpha = 0.5, size = 1) + geom_smooth(method = "lm", se = TRUE, colour = "firebrick") +
          facet_wrap(~ signature, scales = "free_y") + theme_bw(base_size = 9) +
          labs(title = sprintf("%s: population signature score vs NEPC (Beltran custom UP) score", st$label),
               x = "NEPC signature score (mean z)", y = "population signature score (mean z)")
        pdf(file.path(OUT_CROSS, sprintf("nepc_nepc_correlation_scatter_%s.pdf", st$study_id)), width = 12, height = 9)
        print(p); dev.off()
      }
    }

    if (length(summary_rows)) {
      cor_summary <- bind_rows(summary_rows) %>%
        group_by(study_id, method) %>% mutate(p_adj_BH = p.adjust(p, method = "BH")) %>% ungroup() %>%
        mutate(r = round(r, 4), p = signif(p, 4), p_adj_BH = signif(p_adj_BH, 4)) %>%
        arrange(study_id, method, desc(abs(r)))
      write.csv(cor_summary, file.path(OUT_CROSS, "nepc_nepc_correlation_summary.csv"), row.names = FALSE)
      cat("\nCorrelation of recurrent-population signatures with the NEPC signature:\n")
      print(as.data.frame(cor_summary %>% select(study_id, tier, population, signature_type, method, r, p, p_adj_BH, n_samples)))

      cor_wide <- cor_summary %>% filter(method == "pearson") %>%
        select(tier, population, signature_type, study_id, r) %>%
        pivot_wider(names_from = study_id, values_from = r, names_prefix = "pearson_r_")
      write.csv(cor_wide, file.path(OUT_CROSS, "nepc_nepc_correlation_wide.csv"), row.names = FALSE)

      # heatmap of correlation coefficients
      hm <- cor_summary %>% mutate(sig = paste0(tier, ": ", population, "\n", signature_type))
      p_hm <- ggplot(hm, aes(x = study_id, y = sig, fill = r)) + geom_tile() +
        geom_text(aes(label = sprintf("%.2f\np=%.1e", r, p)), size = 2.5) +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) +
        facet_wrap(~ method) + theme_bw(base_size = 9) + labs(x = NULL, y = NULL, fill = "r")
      pdf(file.path(OUT_CROSS, "nepc_nepc_correlation_heatmap.pdf"), width = 8, height = 2 + 0.5 * length(unique(hm$sig)))
      print(p_hm); dev.off()
    }
    if (length(gene_rows)) {
      gene_cor <- bind_rows(gene_rows) %>%
        group_by(study_id) %>% mutate(p_adj_BH = p.adjust(p, method = "BH")) %>% ungroup() %>%
        mutate(rho = round(rho, 4), p = signif(p, 4), p_adj_BH = signif(p_adj_BH, 4)) %>%
        arrange(study_id, population, signature_type, desc(rho))
      write.csv(gene_cor, file.path(OUT_CROSS, "nepc_nepc_correlation_per_gene.csv"), row.names = FALSE)
    }
    if (length(profiles_used)) {
      write.csv(bind_rows(profiles_used), file.path(OUT_CROSS, "nepc_cbioportal_profiles_used.csv"), row.names = FALSE)
    }

    prov <- c(
      sprintf("script_version: %s", SCRIPT_VERSION),
      sprintf("run_date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      sprintf("datasets: %s", paste(vapply(datasets, function(d) d$name, character(1)), collapse = ", ")),
      sprintf("clustering: nfeatures=%d npcs=%d resolution=%g min_features=%d max_mito=%g",
              cfg$nfeatures, cfg$npcs, cfg$resolution, cfg$min_features_cell, cfg$max_mito_pct),
      sprintf("annotation: panel z-score, min_score=%g min_margin=%g; SingleR=%s",
              cfg$ann_min_score, cfg$ann_min_margin, cfg$run_singler),
      sprintf("recurrence: population in > %d %s", cfg$recurrence_min_datasets, cfg$recurrence_count_by),
      sprintf("recurrent populations: %s", paste(recurrent$population, collapse = ", ")),
      sprintf("consensus signature: markers in >= %d datasets, top %d by (n_datasets, mean log2FC); excluded %s",
              cfg$consensus_min_datasets, cfg$consensus_top_n, cfg$exclude_gene_regex),
      sprintf("NEPC-overlapping genes excluded from population signatures: %s", cfg$exclude_nepc_genes_from_signatures),
      sprintf("NEPC signature (%d genes): %s", length(NEPC_BELTRAN_CUSTOM_UP), paste(NEPC_BELTRAN_CUSTOM_UP, collapse = ",")),
      "bulk scoring: log2(x+1) if linear -> gene-wise z across samples -> mean z over signature genes",
      "correlation: Pearson and Spearman of population score vs NEPC score; BH-adjusted within study and method",
      if (length(profiles_used)) paste0("cBioPortal profiles: ",
        paste(vapply(profiles_used, function(p) sprintf("%s=%s", p$study_id, p$profile_id), character(1)), collapse = "; ")) else "cBioPortal profiles: none"
    )
    writeLines(prov, file.path(OUT_CROSS, "nepc_PROVENANCE.txt"))
  }
}

hdr("nepc_ COMPLETE — cross-dataset outputs in %s", OUT_CROSS)
