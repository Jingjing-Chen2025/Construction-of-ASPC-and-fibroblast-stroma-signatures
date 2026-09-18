# NEPC scRNA clustering pipeline

R pipeline for neuroendocrine prostate cancer (NEPC) single-cell RNA-seq datasets:

1. **Per-dataset clustering and cell-population annotation** — every GEO dataset folder is processed on its own: samples are auto-detected, merged, QC-filtered, clustered (Seurat) and every cluster is annotated.
2. **Recurrent populations** — annotated populations are counted across datasets; those present in more than a configurable number of datasets (default: more than 3) are listed.
3. **Correlation with the NEPC signature in bulk cohorts** — signatures of the recurrent populations are scored in PRAD TCGA and SU2C/PCF 2019 (both fetched from cBioPortal) and correlated with the Beltran custom NEPC UP signature.

## Datasets

One folder = one dataset. The default configuration (`config.R`) expects:

| Dataset | Species | Folder |
| --- | --- | --- |
| GSE137829 | human | `<NEPC_ROOT>/GSE137829` |
| GSE210358 | human | `<NEPC_ROOT>/GSE210358` |
| GSE210358_TKO | mouse | `<NEPC_ROOT>/GSE210358_TKO` |
| GSE235036_TKO | mouse | `<NEPC_ROOT>/GSE235036_TKO` |
| GSE264573 | human | `<NEPC_ROOT>/GSE264573` |
| GSE292074 | human | `<NEPC_ROOT>/GSE292074` |
| GSE296986_TKO | mouse | `<NEPC_ROOT>/GSE296986_TKO` |

`NEPC_ROOT` defaults to `/Volumes/Jingjing_Chen/NEPC scRNA dataset`.

**Sample formats.** Inside a dataset folder, one sample is either

* a 10x triplet: `<sample>_barcodes.tsv.gz` + `<sample>_features.tsv.gz` (or `_genes.tsv.gz`) + `<sample>_matrix.mtx.gz`, or
* one file: `<sample>.txt.gz`, `<sample>.csv.gz`, `<sample>.zip` / `<sample>.matrix.zip`, or `<sample>.tar.gz`.

A `.txt.gz` / `.csv.gz` file holds one delimited matrix, either genes x cells (first column = gene), cells x genes (first column = barcode; annotation columns such as `CLUSTER` are dropped) or a `Gene_ID` + `Symbol` + one-column-per-cell table; the orientation is detected automatically. An archive wraps either a 10x directory or one such text matrix. Every other file in the folder is ignored: uncompressed `.csv` / `.txt` files, `.rds`, `.pdf`, sub-directories and the outputs of the companion analysis scripts. A compressed text file is skipped when it has fewer than `cfg$min_fields_expression` columns or matches a known annotation-table name; `cfg$sample_name_regex` (default `NULL`) can additionally restrict single-file samples by name, e.g. `"^GS[ME][0-9]+"`.

## Requirements

R >= 4.1 with: Seurat (v4 or v5), SingleR, celldex, SingleCellExperiment, Matrix, data.table, dplyr, tidyr, tibble, stringr, ggplot2, patchwork, httr, jsonlite.

Recommended (used automatically when installed, skipped with a warning otherwise):

| Package | Used for |
| --- | --- |
| `harmony` | integration by sample (coarse tier) and by dataset + sample (stromal tier) |
| `scDblFinder` | doublet removal per sample |
| `UCell` | per-cell marker-panel scores (fallback: mean of per-gene z-scores) |
| `babelgene` | mouse-to-human ortholog mapping for the stromal tier (fallback: upper-casing) |
| `scuttle` | log-normalising a user-supplied reference that only has counts |

```r
install.packages(c("Seurat", "Matrix", "data.table", "dplyr", "tidyr", "tibble",
                   "stringr", "ggplot2", "patchwork", "httr", "jsonlite", "harmony", "babelgene"))
BiocManager::install(c("SingleR", "celldex", "SingleCellExperiment", "scDblFinder", "UCell", "scuttle"))
```

## Usage

```bash
# 1. edit config.R: NEPC_ROOT, the dataset list and any parameter in cfg
# 2. run everything
Rscript run_pipeline.R
```

Parts can be toggled with `cfg$run_part1`, `cfg$run_part2`, `cfg$run_part3`. Part 1 is skipped for a dataset whose `nepc_clustering/nepc_seurat_annotated.rds` already exists when `cfg$reuse_existing = TRUE`; Parts 2 and 3 only read Part 1 outputs, so they can be re-run quickly.

## Method

### Part 1 — coarse tier: clustering and annotation (per dataset)

* Per-cell QC: minimum detected genes (`cfg$min_features_cell`), mitochondrial percentage (`cfg$max_mito_pct`, default 20), doublets removed per sample with scDblFinder (`cfg$remove_doublets`).
* All samples of a dataset are merged and layers joined, then `NormalizeData` → `FindVariableFeatures(nfeatures = 2000)` → `ScaleData` → `RunPCA(npcs = 30)` → Harmony by sample (`cfg$integrate_samples`) → `FindNeighbors` → `FindClusters` at `cfg$resolution` (0.3) → `RunUMAP`. Clusters are also computed at `cfg$extra_resolutions` and the cluster counts are written to `nepc_clusters_by_resolution.csv` as a stability check.
* `FindAllMarkers` (positive markers, `min.pct` 0.25, `logfc.threshold` 0.25, BH-adjusted p < 0.05).
* Cluster annotation combines three sources:
  * **Marker panels** (`POPULATION_MARKERS` in `config.R`: luminal, basal, club/hillock, neuroendocrine, cycling, fibroblast, ASPC-like adipose progenitor, smooth muscle/myofibroblast, pericyte, endothelial, lymphatic, T, NK, B, plasma, macrophage/myeloid, dendritic, mast, neutrophil, Schwann, erythroid, adipocyte). Per cluster, mean log-normalised expression of every panel gene is z-scored across clusters and averaged per panel; per cell, the same panels are scored with UCell (`score_<panel>` columns in the cell metadata).
  * **SingleR** with celldex `HumanPrimaryCellAtlasData` / `MouseRNAseqData` (`label.main`, on counts), labels harmonised across the two references (`SINGLER_LABEL_MAP`), majority per cluster. Its compartment (epithelial / stromal / endothelial / immune) is compared with the panel label; disagreements are flagged in `compartment_agreement`.
  * **Optional prostate reference** (`cfg$prostate_reference_rds`, e.g. Song et al. 2022 or Henry et al. 2018 saved as a Seurat/SingleCellExperiment object with a label column) run through SingleR and reported as `reference_majority`.
  * `cfg$annotation_method` chooses which label names the cluster: `"panel"` (default) or `"singler"`. The batch-pipeline settings (1000 genes, 20 PCs, resolution 0.05, no mito filter, SingleR majority) are documented in `config.R` and can be restored there.

### Part 1b — stromal tier: integrated ASPC identification (all datasets)

Adipose stem and progenitor cells (ASPC, PDGFRA+ CD34+ DPP4+ PI16+) are a minority of the mesenchymal compartment and are not a class in the celldex references, so the coarse tier labels them fibroblasts. The stromal tier resolves them the way the adipose atlases do:

1. From every dataset's coarse result, mesenchymal cells are taken (clusters labelled fibroblast, ASPC, smooth muscle or pericyte by the panels, or fibroblast / smooth muscle / tissue stem cell by SingleR; `cfg$stromal_populations`, `cfg$stromal_singler_populations`), at most `cfg$stromal_max_cells_per_dataset` per dataset.
2. Mouse genes are mapped to human orthologs (babelgene), all datasets are merged and integrated with Harmony by dataset and sample, and reclustered at `cfg$stromal_resolution` (0.5).
3. Clusters are labelled with the stromal panels (`STROMAL_MARKERS`): ASPC, committed preadipocyte, matrix fibroblast, myCAF, iCAF, smooth muscle, pericyte, Schwann, cycling, and contaminant panels (epithelial, immune, endothelial) that catch doublets and mis-assigned cells.
4. ASPC is called per cell by consensus: the cell sits in an ASPC-labelled cluster **and** its per-cell ASPC score is above threshold (median + MAD of non-ASPC clusters, or `cfg$aspc_cell_min_score`) and above every competing stromal panel. `aspc_call` records `ASPC_consensus`, `ASPC_cluster_only`, `ASPC_cell_only` or `non_ASPC`. An optional adipose reference with ASPC labels (`cfg$stromal_reference_rds`, e.g. Emont et al. 2022 via Azimuth) adds a third, reference-based vote (`reference_majority`).
5. Outputs in `<NEPC_ROOT>/nepc_cross_dataset/stromal/`: the integrated object, cluster annotation with the number of cells per dataset and `n_datasets_present`, per-dataset stromal population counts (fed to Part 2), ASPC calls per dataset, integrated and per-dataset markers (fed to Part 3), and UMAPs by cluster, population, dataset, ASPC call and panel scores.

### Part 2 — recurrent populations

Every non-`Unassigned` population of both tiers (coarse and stromal, `tier` column) with at least `cfg$min_cells_population` cells in a dataset is counted across datasets. A population is *recurrent* when it is present in **more than** `cfg$recurrence_min_datasets` datasets (default 3, i.e. at least 4 of 7). Set `cfg$recurrence_count_by = "clusters"` to count cluster occurrences instead of datasets.

### Part 3 — correlation with the NEPC signature

* Two signatures per recurrent population: the **canonical** marker panel (for SingleR-derived populations via `SINGLER_TO_PANEL`, e.g. `Neuron` → neuroendocrine panel, `Tissue_stem_cell` → ASPC-like progenitor panel) and a **consensus-marker** signature (positive cluster markers of that population found in at least `cfg$consensus_min_datasets` datasets, top `cfg$consensus_top_n` by recurrence and mean log2FC; mitochondrial/ribosomal/haemoglobin genes excluded). Genes that overlap the NEPC signature are excluded from population signatures by default (`cfg$exclude_nepc_genes_from_signatures`) to avoid circular correlations.
* Bulk mRNA expression for `prad_tcga` and `prad_su2c_2019` is fetched from the cBioPortal REST API. The mRNA profile is chosen at run time from the study's non-z-score `MRNA_EXPRESSION` profiles using the preference patterns in `cfg$cbio_studies`, and the choice is logged and written to `nepc_cbioportal_profiles_used.csv`. A profile can be pinned with `profile_id`. Without API access, point `cfg$cbio_local_files` at downloaded cBioPortal datahub expression files.
* Scoring: linear values are log2(x + 1) transformed, every gene is z-scored across samples, a signature score is the mean z over its genes. The NEPC score uses `NEPC_BELTRAN_CUSTOM_UP` (29 genes).
* Correlation: Pearson and Spearman correlation of every population score with the NEPC score across samples (BH-adjusted within study and method), plus Spearman correlation of every signature gene with the NEPC score.

## Outputs

Per dataset, in `<dataset folder>/nepc_clustering/`:

| File | Content |
| --- | --- |
| `nepc_sample_manifest.csv` | detected samples, format, load status, cell/gene counts |
| `nepc_cluster_annotation.csv` | one row per cluster: size, population (SingleR majority by default), majority share, panel label and score, top markers |
| `nepc_panel_scores_by_cluster.csv` | full cluster x panel score matrix |
| `nepc_cluster_markers_sig.csv`, `nepc_cluster_markers_top.csv` | FindAllMarkers results (all significant / top N per cluster) |
| `nepc_population_counts.csv` | cells per population |
| `nepc_cell_metadata.csv` | per-cell metadata incl. cluster, population, SingleR label |
| `nepc_umap_*.pdf` | UMAPs: named clusters, populations, sample, SingleR cell types, clusters + cell types |
| `nepc_seurat_annotated.rds` | annotated Seurat object |

Cross-dataset, in `<NEPC_ROOT>/nepc_cross_dataset/`:

| File | Content |
| --- | --- |
| `nepc_run_log.csv` | per-dataset status, cell and cluster counts, run time |
| `nepc_all_cluster_annotations.csv` | all cluster annotations stacked |
| `nepc_population_recurrence.csv` | every population with number of datasets/clusters and where it occurs |
| `nepc_population_presence_matrix.csv` | population x dataset matrix of cluster counts |
| `nepc_recurrent_populations.csv` | populations (coarse and stromal tier) present in more than `cfg$recurrence_min_datasets` datasets |
| `stromal/nepc_stromal_cluster_annotation.csv` | integrated stromal clusters: population, panel score, fraction of cells passing the ASPC score, cells per dataset |
| `stromal/nepc_aspc_calls_by_dataset.csv` | ASPC consensus / cluster-only / cell-only / non-ASPC counts per dataset |
| `stromal/nepc_stromal_markers_integrated.csv`, `stromal/nepc_stromal_markers_by_dataset.csv` | stromal population markers (integrated, and per dataset for consensus signatures) |
| `stromal/nepc_stromal_umap*.pdf`, `stromal/nepc_stromal_integrated.rds` | stromal-tier plots and object |
| `nepc_recurrent_population_signatures.csv` | canonical and consensus signatures per recurrent population |
| `nepc_cbioportal_profiles_used.csv` | study, molecular profile and sample list used |
| `nepc_bulk_log2_expression_<study>.csv` | fetched expression (log2) |
| `nepc_bulk_signature_scores_<study>.csv` | per-sample NEPC and population scores |
| `nepc_nepc_correlation_summary.csv`, `nepc_nepc_correlation_wide.csv` | signature-level correlations with the NEPC score |
| `nepc_nepc_correlation_per_gene.csv` | gene-level Spearman correlations with the NEPC score |
| `nepc_nepc_correlation_scatter_<study>.pdf`, `nepc_nepc_correlation_heatmap.pdf` | plots |
| `nepc_PROVENANCE.txt` | parameters and profiles used for the run |

## Notes

* Malignant versus normal epithelium is not resolved by this pipeline; run inferCNV or copyKAT on the epithelial cells of the coarse result, using immune and stromal cells as the diploid reference.
* Mouse datasets are annotated with the same human panels via uppercase symbol matching; consensus signatures pool human and mouse markers the same way.
