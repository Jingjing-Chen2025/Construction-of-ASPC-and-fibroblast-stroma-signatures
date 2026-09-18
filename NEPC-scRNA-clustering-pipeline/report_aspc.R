# Standalone ASPC report: summarises the result folders without loading Seurat.
#   Rscript report_aspc.R          (run from the pipeline folder, or setwd() to it)
suppressPackageStartupMessages({ library(dplyr) })
PIPELINE_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f[1])))
  for (fr in sys.frames()) { of <- fr$ofile; if (!is.null(of) && nzchar(of)) return(dirname(normalizePath(of))) }
  getwd()
})
source(file.path(PIPELINE_DIR, "config.R")); source(file.path(PIPELINE_DIR, "R", "helpers.R")); source(file.path(PIPELINE_DIR, "R", "report.R"))
write_aspc_report()
