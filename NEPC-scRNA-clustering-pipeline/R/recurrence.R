# ============================================================================
# recurrence.R — PART 2 helpers: collect Part 1 outputs across datasets
# ============================================================================
read_all_annotations <- function() {
  tbls <- lapply(datasets, function(ds) {
    f <- file.path(ds$dir, "nepc_clustering", "nepc_cluster_annotation.csv")
    if (!file.exists(f)) { warn_msg("Missing annotation for %s (%s)", ds$name, f); return(NULL) }
    read.csv(f, colClasses = c(cluster = "character"))
  })
  bind_rows(tbls)
}

read_all_markers <- function() {
  tbls <- lapply(datasets, function(ds) {
    f <- file.path(ds$dir, "nepc_clustering", "nepc_cluster_markers_sig.csv")
    if (!file.exists(f)) return(NULL)
    read.csv(f, colClasses = c(cluster = "character"))
  })
  bind_rows(tbls)
}
