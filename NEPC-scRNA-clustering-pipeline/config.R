# ============================================================================
# config.R — paths, parameters, dataset list, signatures and marker panels
# Edit this file to point the pipeline at your data.
# ============================================================================
# CONFIG
# ============================================================================
NEPC_ROOT <- "/Volumes/Jingjing_Chen/NEPC scRNA dataset"
OUT_CROSS <- file.path(NEPC_ROOT, "nepc_cross_dataset")

cfg <- list(
  run_part1 = TRUE,
  run_part2 = TRUE,
  run_part3 = TRUE,
  reuse_existing = TRUE,          # skip Part 1 for datasets already processed

  # ---- Part 1: QC + clustering ----
  min_cells_gene       = 3,
  min_features_cell    = 200,
  max_features_cell    = Inf,     # e.g. 8000 to drop likely doublets
  max_mito_pct         = Inf,     # batch-pipeline setting: no mitochondrial filter (e.g. 20 to enable)
  nfeatures            = 1000,    # batch-pipeline setting
  npcs                 = 20,      # batch-pipeline setting
  resolution           = 0.05,    # batch-pipeline setting
  max_cells_per_sample = Inf,     # e.g. 10000 to subsample huge samples
  chunk_cols           = 2000,    # dense -> sparse conversion chunk size
  min_cells_dataset    = 100,     # abort a dataset with fewer cells after QC
  min_fields_expression = 10,     # a text matrix needs at least this many columns
  sample_name_regex    = NULL,    # e.g. "^GS[ME][0-9]+" to accept only GEO-named files

  # ---- Part 1: markers ----
  marker_min_pct   = 0.25,
  marker_logfc     = 0.25,
  marker_p_adj     = 0.05,
  top_markers_report = 30,

  # ---- Part 1: annotation ----
  # "singler": each cluster is named by the majority SingleR label of its cells
  #            (celldex HumanPrimaryCellAtlasData / MouseRNAseqData, label.main,
  #            run on counts) — the batch-pipeline setting. The marker-panel
  #            label is still computed and reported as population_panel.
  # "panel":   name clusters by the canonical marker panels instead.
  annotation_method  = "singler",
  run_singler        = TRUE,
  singler_max_cells  = 20000,
  singler_assay      = "counts",  # batch-pipeline setting ("logcounts" also works)
  singler_labels     = "label.main",
  ann_min_genes  = 3,             # panel genes that must be detected
  ann_min_score  = 0.30,          # mean z of the winning panel
  ann_min_margin = 0.10,          # winner minus runner-up
  save_rds           = TRUE,

  # ---- Part 2: recurrence ----
  recurrence_min_datasets = 3,    # populations in MORE THAN this many datasets
  recurrence_count_by     = "datasets",   # "datasets" | "clusters"

  # ---- Part 3: signatures ----
  consensus_min_datasets = 2,
  consensus_top_n        = 50,
  exclude_gene_regex     = "^(MT-|RP[SL][0-9]|RPLP|MRP[SL]|HB[AB])",
  exclude_nepc_genes_from_signatures = TRUE,   # avoid circular overlap

  # ---- Part 3: cBioPortal ----
  cbio_base      = "https://www.cbioportal.org/api",
  cbio_timeout   = 180,
  cbio_retries   = 4,
  cbio_use_cache = TRUE,
  cbio_gene_batch = 200,
  cbio_studies = list(
    prad_tcga = list(
      study_id = "prad_tcga",
      label    = "PRAD TCGA",
      # explicit profile wins; otherwise these patterns are tried in order,
      # then any non-z-score MRNA_EXPRESSION profile
      profile_id       = NULL,
      profile_patterns = c("rna_seq_v2_mrna$", "rna_seq_mrna$", "mrna_seq_fpkm", "mrna$")
    ),
    prad_su2c_2019 = list(
      study_id = "prad_su2c_2019",
      label    = "SU2C/PCF 2019 (mCRPC)",
      profile_id       = NULL,
      profile_patterns = c("fpkm_polya$", "polya$", "fpkm_capture$", "capture$",
                           "rna_seq_mrna$", "mrna$")
    )
  ),
  # Optional offline fallback: named list study_id -> path to a cBioPortal
  # datahub expression file (data_mrna_seq_v2_rsem.txt / data_mrna_seq_fpkm*.txt)
  cbio_local_files = list(),

  min_samples_cor = 10
)

datasets <- list(
  list(name = "GSE137829",     dir = file.path(NEPC_ROOT, "GSE137829"),     species = "human"),
  list(name = "GSE210358",     dir = file.path(NEPC_ROOT, "GSE210358"),     species = "human"),
  list(name = "GSE210358_TKO", dir = file.path(NEPC_ROOT, "GSE210358_TKO"), species = "mouse"),
  list(name = "GSE235036_TKO", dir = file.path(NEPC_ROOT, "GSE235036_TKO"), species = "mouse"),
  list(name = "GSE264573",     dir = file.path(NEPC_ROOT, "GSE264573"),     species = "human"),
  list(name = "GSE292074",     dir = file.path(NEPC_ROOT, "GSE292074"),     species = "human"),
  list(name = "GSE296986_TKO", dir = file.path(NEPC_ROOT, "GSE296986_TKO"), species = "mouse")
)

# ============================================================================
# SIGNATURES
# ============================================================================
NEPC_BELTRAN_CUSTOM_UP <- toupper(c(
  "ASXL3","CAND2","ETV5","GPX2","JAKMIP2","KIAA0408","SOGA3","TRIM9","BRINP1",
  "C7ORF76","GNAO1","KCNB2","KCND2","LRRC16B","MAP10","NRSN1",
  "PCSK1","PROX1","RGS7","SCG3","SEC11C","SEZ6","ST8SIA3","SVOP","SYT11",
  "AURKA","DNMT1","EZH2","MYCN"
))

# Symbols that cBioPortal may know under a different HUGO name. Tried only
# when the primary symbol is not found.
GENE_ALIASES <- c(
  LRRC16B  = "CARMIL3",
  C7ORF76  = "C7orf76",
  KIAA0408 = "KIAA0408",
  ACPP     = "ACP3",
  FAM64A   = "PIMREG",
  WISP2    = "CCN5",
  CTGF     = "CCN2",
  CYR61    = "CCN1"
)

# Canonical marker panels used for cluster annotation (human symbols; mouse
# symbols are matched after toupper()). Each panel = one candidate population.
POPULATION_MARKERS <- list(
  Luminal_epithelial          = c("KLK3","KLK2","AR","NKX3-1","TMPRSS2","KRT8","KRT18","EPCAM","FOLH1","MSMB","ACPP"),
  Basal_epithelial            = c("KRT5","KRT14","KRT15","TP63","DST","LGALS7","KRT17"),
  Club_Hillock_epithelial     = c("SCGB1A1","SCGB3A1","PIGR","KRT13","KRT4","LTF","MMP7","WFDC2"),
  Neuroendocrine              = c("CHGA","CHGB","SYP","ENO2","ASCL1","INSM1","NCAM1","SCG2","SCG3","NEUROD1","FOXA2","PCSK1","SOX2","STMN2","TUBB3"),
  Cycling                     = c("MKI67","TOP2A","CENPF","BIRC5","UBE2C","CCNB1","STMN1","TYMS"),
  Fibroblast                  = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","FBLN1","C7","APOD","SFRP2"),
  ASPC_adipose_progenitor     = c("PDGFRA","DPP4","PI16","CD34","CD55","ANXA3","WNT2","CD248","SEMA3C","MFAP5","CLEC3B","IGFBP6"),
  Smooth_muscle_myofibroblast = c("ACTA2","MYH11","TAGLN","CNN1","DES","MYL9","ACTG2","TPM2"),
  Pericyte                    = c("RGS5","PDGFRB","MCAM","NOTCH3","KCNJ8","ABCC9","HIGD1B"),
  Endothelial                 = c("PECAM1","VWF","CDH5","CLDN5","FLT1","KDR","EMCN","PLVAP"),
  Lymphatic_endothelial       = c("PROX1","LYVE1","PDPN","CCL21","FLT4","MMRN1"),
  T_cell                      = c("CD3D","CD3E","CD3G","CD2","TRAC","CD8A","CD4","IL7R","CD7"),
  NK_cell                     = c("NKG7","GNLY","KLRD1","KLRF1","NCR1","PRF1","GZMB"),
  B_cell                      = c("MS4A1","CD79A","CD79B","CD19","BANK1","PAX5"),
  Plasma_cell                 = c("MZB1","JCHAIN","SDC1","XBP1","IGKC","TNFRSF17","PRDM1"),
  Macrophage_myeloid          = c("CD68","CD14","LYZ","C1QA","C1QB","CD163","AIF1","FCGR3A","CSF1R","MS4A7"),
  Dendritic_cell              = c("CLEC9A","CD1C","FCER1A","LAMP3","IRF8","CLEC10A","XCR1"),
  Mast_cell                   = c("TPSAB1","TPSB2","CPA3","KIT","MS4A2","HDC"),
  Neutrophil                  = c("S100A8","S100A9","CSF3R","FCGR3B","CXCR2","G0S2","RETNLG"),
  Schwann_neural              = c("PLP1","S100B","MPZ","SOX10","NGFR","CDH19"),
  Erythroid                   = c("HBB","HBA1","HBA2","ALAS2","HBA-A1","HBA-A2","HBB-BS"),
  Adipocyte                   = c("ADIPOQ","LEP","PLIN1","CIDEA","FABP4","CIDEC")
)
POPULATION_MARKERS <- lapply(POPULATION_MARKERS, toupper)

# SingleR label.main vocabularies differ between the human (HPCA) and mouse
# (MouseRNAseq) references; both are mapped onto one shared vocabulary so that a
# population can be counted across human and mouse datasets. Unlisted labels are
# normalised automatically (spaces -> "_", plural "cells" -> "cell").
SINGLER_LABEL_MAP <- c(
  "T_cells" = "T_cell", "T cells" = "T_cell",
  "B_cell" = "B_cell", "B cells" = "B_cell", "Pre-B_cell_CD34-" = "B_cell", "Pro-B_cell_CD34+" = "B_cell",
  "NK_cell" = "NK_cell", "NK cells" = "NK_cell",
  "Macrophage" = "Macrophage", "Macrophages" = "Macrophage", "Microglia" = "Macrophage",
  "Monocyte" = "Monocyte", "Monocytes" = "Monocyte",
  "DC" = "Dendritic_cell", "Dendritic cells" = "Dendritic_cell",
  "Granulocytes" = "Granulocyte", "Neutrophils" = "Granulocyte", "Myelocyte" = "Granulocyte",
  "Pro-Myelocyte" = "Granulocyte", "GMP" = "Myeloid_progenitor", "CMP" = "Myeloid_progenitor",
  "MEP" = "Myeloid_progenitor", "HSC_-G-CSF" = "Myeloid_progenitor", "HSC_CD34+" = "Myeloid_progenitor",
  "Erythroblast" = "Erythroid", "Erythrocytes" = "Erythroid", "Platelets" = "Platelet",
  "Fibroblasts" = "Fibroblast",
  "Endothelial_cells" = "Endothelial_cell", "Endothelial cells" = "Endothelial_cell",
  "Epithelial_cells" = "Epithelial_cell", "Epithelial cells" = "Epithelial_cell",
  "Keratinocytes" = "Epithelial_cell", "Hepatocytes" = "Hepatocyte",
  "Smooth_muscle_cells" = "Smooth_muscle_cell", "Cardiomyocytes" = "Smooth_muscle_cell",
  "Neurons" = "Neuron", "Neuroepithelial_cell" = "Neuron",
  "Astrocyte" = "Glia", "Astrocytes" = "Glia", "Oligodendrocytes" = "Glia", "OPCs" = "Glia",
  "Tissue_stem_cells" = "Tissue_stem_cell", "MSC" = "Tissue_stem_cell",
  "Chondrocytes" = "Chondrocyte", "Osteoblasts" = "Osteoblast",
  "Adipocytes" = "Adipocyte",
  "iPS_cells" = "Stem_cell", "Embryonic_stem_cells" = "Stem_cell"
)

# Canonical marker panel used for a SingleR-derived population in Part 3.
SINGLER_TO_PANEL <- c(
  T_cell = "T_cell", B_cell = "B_cell", NK_cell = "NK_cell",
  Macrophage = "Macrophage_myeloid", Monocyte = "Macrophage_myeloid", Dendritic_cell = "Dendritic_cell",
  Granulocyte = "Neutrophil", Erythroid = "Erythroid",
  Fibroblast = "Fibroblast", Endothelial_cell = "Endothelial",
  Smooth_muscle_cell = "Smooth_muscle_myofibroblast",
  Neuron = "Neuroendocrine", Tissue_stem_cell = "ASPC_adipose_progenitor",
  Adipocyte = "Adipocyte", Epithelial_cell = "Epithelial_cell", Plasma_cell = "Plasma_cell"
)

# Extra canonical signatures for populations without a panel of their own.
EXTRA_SIGNATURES <- list(
  Epithelial_cell = toupper(c("EPCAM","KRT8","KRT18","KRT19","CDH1","CLDN4","CLDN7","KRT5","KRT14","KLK3","AR"))
)
