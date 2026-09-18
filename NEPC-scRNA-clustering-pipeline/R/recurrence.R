# ============================================================================
# recurrence.R — PART 2 helpers: collect Part 1 outputs across datasets
# ============================================================================
read_all_annotations <- function(group) {
  tbls <- lapply(group_datasets(group), function(ds) {
    f <- file.path(dataset_out_dir(ds), "nepc_cluster_annotation.csv")
    if (!file.exists(f)) { warn_msg("Missing annotation for %s (%s)", ds$name, f); return(NULL) }
    d <- read.csv(f, colClasses = c(cluster = "character")); d$tier <- "coarse"; d$group <- group; d
  })
  coarse <- bind_rows(tbls)
  # stromal tier (Part 1b): one row per dataset x stromal population
  f <- file.path(OUT_CROSS, paste0("stromal_", group), "nepc_stromal_population_by_dataset.csv")
  if (file.exists(f)) {
    st <- read.csv(f, colClasses = c(cluster = "character"))
    msg("Stromal-tier populations added to the recurrence table (%d rows)", nrow(st))
    coarse <- bind_rows(coarse, st)
  }
  coarse
}

read_all_markers <- function(group) {
  tbls <- lapply(group_datasets(group), function(ds) {
    f <- file.path(dataset_out_dir(ds), "nepc_cluster_markers_sig.csv")
    if (!file.exists(f)) return(NULL)
    d <- read.csv(f, colClasses = c(cluster = "character")); d$tier <- "coarse"; d
  })
  mk <- bind_rows(tbls)
  f <- file.path(OUT_CROSS, paste0("stromal_", group), "nepc_stromal_markers_by_dataset.csv")
  if (file.exists(f)) mk <- bind_rows(mk, read.csv(f, colClasses = c(cluster = "character")))
  mk
}

# ============================================================================
# PART 2 — recurrent populations of one group = group-specific cell sets
# ============================================================================
run_recurrence <- function(group) {
  hdr("PART 2 — %s: populations present in > %d %s = %s-specific cell sets",
      group, cfg$recurrence_min_datasets, cfg$recurrence_count_by, group)
  all_ann <- read_all_annotations(group)
  if (!nrow(all_ann)) stop("No Part 1 annotation tables found; run Part 1 first.")
  write.csv(all_ann, file.path(OUT_CROSS, paste0("nepc_all_cluster_annotations_", group, ".csv")), row.names = FALSE)
  
  all_ann <- all_ann %>% filter(n_cells >= ifelse(tier == "stromal", cfg$min_cells_population_stromal, cfg$min_cells_population))
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
  write.csv(recurrence, file.path(OUT_CROSS, paste0("nepc_population_recurrence_", group, ".csv")), row.names = FALSE)
  
  presence <- all_ann %>% filter(!population %in% c("Unassigned", "Contaminant")) %>%
    count(tier, population, dataset, name = "n_clusters") %>%
    pivot_wider(names_from = dataset, values_from = n_clusters, values_fill = 0)
  write.csv(presence, file.path(OUT_CROSS, paste0("nepc_population_presence_matrix_", group, ".csv")), row.names = FALSE)
  
  recurrent <- recurrence %>% filter(recurrent)
  recurrent$group <- group
    recurrent$cell_set <- paste0(group, "_specific:", recurrent$tier, ":", recurrent$population)
    write.csv(recurrent, file.path(OUT_CROSS, paste0(group, "_specific_cell_sets.csv")), row.names = FALSE)
  
  cat("\n", group, "-specific cell sets (populations in more than ", cfg$recurrence_min_datasets, " ", cfg$recurrence_count_by, "):\n", sep = "")
  if (nrow(recurrent)) {
    print(as.data.frame(recurrent[, c("tier", "population", "n_datasets", "n_clusters", "total_cells", "datasets")]))
  } else cat("  (none)\n")
  
  list(recurrent = recurrent, recurrence = recurrence, all_ann = all_ann)
}
