# ============================================================================
# composition.R — cell composition of every population per sample and group
# ============================================================================
# Coarse tier: % of all cells of a sample in each population.
# Stromal tier: % of the sample's stromal cells (stromal_population), and the
#   same rescaled to % of all cells using the stromal fraction of the sample
#   recorded at extraction (nepc_stromal_extraction_by_sample.csv).
# Means are reported across samples and as the mean of dataset means, so that a
# dataset with many samples does not dominate the group estimate.
# ============================================================================

# cell-metadata CSVs are written with row names (empty first header field): header = TRUE is required
read_cell_metadata <- function(f, cols) {
  d <- data.table::fread(f, data.table = FALSE, header = TRUE)
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop("columns missing in ", basename(f), ": ", paste(miss, collapse = ", "))
  d[, cols, drop = FALSE]
}

composition_by_sample <- function(group) {
  rows <- list()
  # ---- coarse tier ----
  for (ds in group_datasets(group)) {
    f <- file.path(dataset_out_dir(ds), "nepc_cell_metadata.csv")
    if (!file.exists(f)) { warn_msg("  %s: no cell metadata (%s)", ds$name, f); next }
    md <- read_cell_metadata(f, c("sample", "population"))
    md$population[is.na(md$population) | md$population == ""] <- "Unassigned"
    cnt <- md %>% count(sample, population, name = "n") %>%
      group_by(sample) %>% mutate(n_total = sum(n)) %>% ungroup()
    rows[[length(rows) + 1]] <- cnt %>%
      transmute(group = group, dataset = ds$name, species = ds$species, tier = "coarse", sample, population,
                n_cells = n, n_denominator = n_total, pct = 100 * n / n_total, pct_of_all = 100 * n / n_total)
  }
  # ---- stromal tier ----
  sdir <- file.path(OUT_CROSS, paste0("stromal_", group))
  f_md <- file.path(sdir, "nepc_stromal_cell_metadata.csv"); f_ex <- file.path(sdir, "nepc_stromal_extraction_by_sample.csv")
  if (file.exists(f_md) && file.exists(f_ex)) {
    smd <- read_cell_metadata(f_md, c("dataset", "sample", "stromal_population"))
    ex  <- read.csv(f_ex)
    sp  <- vapply(group_datasets(group), function(d) d$species, character(1)); names(sp) <- vapply(group_datasets(group), function(d) d$name, character(1))
    cnt <- smd %>% filter(!stromal_population %in% c("Contaminant", "Unassigned")) %>%
      count(dataset, sample, stromal_population, name = "n") %>%
      group_by(dataset, sample) %>% mutate(n_stromal_clean = sum(n)) %>% ungroup() %>%
      left_join(ex, by = c("dataset", "sample")) %>%
      mutate(frac_used_clean = n_stromal_clean / pmax(n_stromal_used, 1),
             pct = 100 * n / n_stromal_clean,
             pct_of_all = pct / 100 * (n_stromal_cells * frac_used_clean) / n_total_cells * 100)
    rows[[length(rows) + 1]] <- cnt %>%
      transmute(group = group, dataset, species = unname(sp[dataset]), tier = "stromal", sample,
                population = stromal_population, n_cells = n, n_denominator = n_stromal_clean, pct, pct_of_all)
  } else msg("  %s: no stromal-tier metadata; composition from the coarse tier only", group)
  if (!length(rows)) return(NULL)
  bs <- bind_rows(rows)
  # zero-fill populations absent from a sample (within tier)
  bs <- bs %>% group_by(tier) %>%
    tidyr::complete(tidyr::nesting(group, dataset, species, sample, n_denominator), population,
                    fill = list(n_cells = 0L, pct = 0, pct_of_all = 0)) %>% ungroup()
  bs
}

run_composition <- function(group, recurrent) {
  hdr("Cell composition — %s samples", group)
  bs <- composition_by_sample(group)
  if (is.null(bs) || !nrow(bs)) { warn_msg("No composition data for %s", group); return(NULL) }
  small <- bs %>% filter(tier == "coarse", n_denominator < cfg$composition_min_cells_sample) %>% distinct(dataset, sample)
  if (nrow(small)) msg("  %d sample(s) with < %d cells excluded from the means", nrow(small), cfg$composition_min_cells_sample)
  bs <- bs %>% anti_join(small, by = c("dataset", "sample"))
  sets <- if (!is.null(recurrent) && nrow(recurrent)) paste(recurrent$tier, recurrent$population) else character(0)
  bs$in_specific_cell_set <- paste(bs$tier, bs$population) %in% sets
  write.csv(bs %>% mutate(pct = round(pct, 3), pct_of_all = round(pct_of_all, 3)),
            file.path(OUT_CROSS, paste0("nepc_cell_composition_by_sample_", group, ".csv")), row.names = FALSE)

  ds_means <- bs %>% group_by(tier, population, dataset) %>% summarise(pct = mean(pct), pct_of_all = mean(pct_of_all), .groups = "drop")
  mean_tbl <- bs %>% group_by(tier, population, in_specific_cell_set) %>%
    summarise(n_samples = n(), n_samples_present = sum(n_cells > 0), n_datasets_present = n_distinct(dataset[n_cells > 0]),
              mean_pct = mean(pct), sd_pct = stats::sd(pct), median_pct = stats::median(pct),
              mean_pct_of_all = mean(pct_of_all), sd_pct_of_all = stats::sd(pct_of_all), .groups = "drop") %>%
    left_join(ds_means %>% group_by(tier, population) %>%
                summarise(mean_of_dataset_means_pct = mean(pct), mean_of_dataset_means_pct_of_all = mean(pct_of_all), .groups = "drop"),
              by = c("tier", "population")) %>%
    mutate(group = group, denominator = ifelse(tier == "coarse", "all cells of the sample", "stromal cells of the sample")) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    arrange(tier, desc(in_specific_cell_set), desc(mean_pct))
  write.csv(mean_tbl, file.path(OUT_CROSS, paste0("nepc_cell_composition_mean_", group, ".csv")), row.names = FALSE)
  cell_sets <- mean_tbl %>% filter(in_specific_cell_set)
  write.csv(cell_sets, file.path(OUT_CROSS, paste0(group, "_specific_cell_sets_mean_composition.csv")), row.names = FALSE)

  cat("\nMean composition of the ", group, "-specific cell sets (% of sample, mean across samples):\n", sep = "")
  print(as.data.frame(cell_sets[, c("tier", "population", "n_samples", "mean_pct", "sd_pct", "mean_pct_of_all", "mean_of_dataset_means_pct")]))

  # plots: one panel per tier, samples as points
  pl <- bs %>% mutate(population = factor(population, levels = mean_tbl$population[order(mean_tbl$tier, -mean_tbl$mean_pct)] %>% unique()))
  p <- ggplot(pl, aes(x = population, y = pct, fill = in_specific_cell_set)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_jitter(aes(colour = dataset), width = 0.2, size = 0.8) +
    facet_wrap(~ tier, scales = "free", ncol = 1) + theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(title = paste0(group, ": cell composition per sample"), x = NULL, y = "% of cells (coarse: all cells; stromal: stromal cells)",
         fill = paste0(group, "-specific cell set"))
  pdf(file.path(OUT_CROSS, paste0("nepc_cell_composition_", group, ".pdf")), width = 14, height = 10); print(p); dev.off()
  list(by_sample = bs, mean = mean_tbl)
}

# NEPC vs adenocarcinoma: same population, sample-level % (Wilcoxon), both tiers
compare_composition <- function(comp) {
  hdr("Cell composition — %s", paste(names(comp), collapse = " vs "))
  bs <- bind_rows(lapply(comp, function(x) x$by_sample))
  groups <- unique(bs$group)
  if (length(groups) < 2) return(NULL)
  g1 <- groups[1]; g2 <- groups[2]
  out <- bs %>% group_by(tier, population) %>%
    summarise(n_samples_1 = sum(group == g1), n_samples_2 = sum(group == g2),
              mean_pct_1 = mean(pct[group == g1]), mean_pct_2 = mean(pct[group == g2]),
              mean_pct_of_all_1 = mean(pct_of_all[group == g1]), mean_pct_of_all_2 = mean(pct_of_all[group == g2]),
              p_wilcoxon = if (n_samples_1 >= 3 && n_samples_2 >= 3)
                suppressWarnings(stats::wilcox.test(pct[group == g1], pct[group == g2])$p.value) else NA_real_,
              .groups = "drop") %>%
    filter(n_samples_1 > 0 & n_samples_2 > 0) %>%
    mutate(p_adj_BH = p.adjust(p_wilcoxon, method = "BH"),
           log2_ratio = log2((mean_pct_1 + 0.1) / (mean_pct_2 + 0.1))) %>%
    mutate(across(where(is.numeric), ~ signif(.x, 4))) %>% arrange(tier, p_wilcoxon)
  names(out) <- sub("_1$", paste0("_", g1), names(out)); names(out) <- sub("_2$", paste0("_", g2), names(out))
  write.csv(out, file.path(OUT_CROSS, sprintf("nepc_cell_composition_%s_vs_%s.csv", g1, g2)), row.names = FALSE)
  p <- ggplot(bs, aes(x = population, y = pct, fill = group)) + geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    facet_wrap(~ tier, scales = "free", ncol = 1) + theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(title = sprintf("Cell composition per sample: %s vs %s", g1, g2), x = NULL, y = "% of cells")
  pdf(file.path(OUT_CROSS, sprintf("nepc_cell_composition_%s_vs_%s.pdf", g1, g2)), width = 14, height = 10); print(p); dev.off()
  cat("\nTop composition differences (", g1, " vs ", g2, "):\n", sep = "")
  print(as.data.frame(head(out, 15)))
  out
}
