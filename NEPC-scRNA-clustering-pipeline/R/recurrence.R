# ============================================================================
# recurrence.R — PART 2 helpers: collect Part 1 outputs across datasets
# ============================================================================
read_all_annotations <- function() {
  tbls <- lapply(datasets, function(ds) {
    f <- file.path(ds$dir, "nepc_clustering", "nepc_cluster_annotation.csv")
    if (!file.exists(f)) { warn_msg("Missing annotation for %s (%s)", ds$name, f); return(NULL) }
    d <- read.csv(f, colClasses = c(cluster = "character")); d$tier <- "coarse"; d
  })
  coarse <- bind_rows(tbls)
  # stromal tier (Part 1b): one row per dataset x stromal population
  f <- file.path(OUT_CROSS, "stromal", "nepc_stromal_population_by_dataset.csv")
  if (file.exists(f)) {
    st <- read.csv(f, colClasses = c(cluster = "character"))
    msg("Stromal-tier populations added to the recurrence table (%d rows)", nrow(st))
    coarse <- bind_rows(coarse, st)
  }
  coarse
}

read_all_markers <- function() {
  tbls <- lapply(datasets, function(ds) {
    f <- file.path(ds$dir, "nepc_clustering", "nepc_cluster_markers_sig.csv")
    if (!file.exists(f)) return(NULL)
    d <- read.csv(f, colClasses = c(cluster = "character")); d$tier <- "coarse"; d
  })
  mk <- bind_rows(tbls)
  f <- file.path(OUT_CROSS, "stromal", "nepc_stromal_markers_by_dataset.csv")
  if (file.exists(f)) mk <- bind_rows(mk, read.csv(f, colClasses = c(cluster = "character")))
  mk
}
