# ============================================================================
# signatures_cbioportal.R — PART 3: signatures, cBioPortal access, scoring, correlation
# ============================================================================
# ---- 3a. signature construction -------------------------------------------
build_consensus_signature <- function(pop, all_markers) {
  if (is.null(all_markers) || !nrow(all_markers)) return(NULL)
  df <- all_markers %>%
    filter(population == pop, avg_log2FC > 0) %>%
    mutate(gene = toupper(gene)) %>%
    filter(!grepl(cfg$exclude_gene_regex, gene)) %>%
    group_by(dataset, gene) %>%
    summarise(avg_log2FC = max(avg_log2FC), .groups = "drop") %>%
    group_by(gene) %>%
    summarise(n_datasets = n_distinct(dataset), mean_log2FC = mean(avg_log2FC), .groups = "drop") %>%
    filter(n_datasets >= cfg$consensus_min_datasets) %>%
    arrange(desc(n_datasets), desc(mean_log2FC)) %>%
    slice_head(n = cfg$consensus_top_n)
  if (!nrow(df)) return(NULL)
  df
}

# canonical panel for a population name from either annotation method
canonical_signature_for <- function(pop) {
  if (!is.null(POPULATION_MARKERS[[pop]])) return(POPULATION_MARKERS[[pop]])
  if (!is.null(EXTRA_SIGNATURES[[pop]]))   return(EXTRA_SIGNATURES[[pop]])
  panel <- unname(SINGLER_TO_PANEL[pop])
  if (!is.na(panel)) {
    if (!is.null(POPULATION_MARKERS[[panel]])) return(POPULATION_MARKERS[[panel]])
    if (!is.null(EXTRA_SIGNATURES[[panel]]))   return(EXTRA_SIGNATURES[[panel]])
  }
  NULL
}

build_signatures <- function(recurrent, all_markers) {
  out <- list()
  for (pop in recurrent$population) {
    canon <- canonical_signature_for(pop)
    if (is.null(canon)) warn_msg("No canonical panel for population '%s'; consensus markers only", pop)
    if (!is.null(canon)) {
      out[[length(out) + 1]] <- data.frame(population = pop, signature_type = "canonical", gene = canon,
                                           n_datasets = NA_integer_, mean_log2FC = NA_real_)
    }
    cons <- build_consensus_signature(pop, all_markers)
    if (!is.null(cons)) {
      out[[length(out) + 1]] <- data.frame(population = pop, signature_type = "consensus_markers",
                                           gene = cons$gene, n_datasets = cons$n_datasets,
                                           mean_log2FC = round(cons$mean_log2FC, 3))
    } else warn_msg("No consensus markers for %s (need >= %d datasets)", pop, cfg$consensus_min_datasets)
  }
  sig <- bind_rows(out)
  if (!nrow(sig)) return(sig)
  sig$in_nepc_signature <- sig$gene %in% NEPC_BELTRAN_CUSTOM_UP
  sig$used_for_scoring  <- !(isTRUE(cfg$exclude_nepc_genes_from_signatures) & sig$in_nepc_signature)
  sig
}

# ---- 3b. cBioPortal access ---------------------------------------------------
cbio_cache_dir <- file.path(OUT_CROSS, "cbioportal_cache")

cbio_request <- function(path, query = NULL, body_json = NULL, cache_key = NULL) {
  if (!requireNamespace("httr", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE))
    stop("Packages 'httr' and 'jsonlite' are required for cBioPortal access.")
  cf <- if (!is.null(cache_key)) file.path(cbio_cache_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", cache_key), ".rds")) else NULL
  if (!is.null(cf) && isTRUE(cfg$cbio_use_cache) && file.exists(cf)) return(readRDS(cf))
  url <- paste0(cfg$cbio_base, path)
  last_err <- "unknown"
  for (attempt in seq_len(cfg$cbio_retries)) {
    resp <- tryCatch({
      if (is.null(body_json)) {
        httr::GET(url, query = query, httr::accept_json(), httr::timeout(cfg$cbio_timeout))
      } else {
        httr::POST(url, query = query, body = body_json, httr::content_type_json(),
                   httr::accept_json(), httr::timeout(cfg$cbio_timeout))
      }
    }, error = function(e) { last_err <<- conditionMessage(e); NULL })
    if (!is.null(resp)) {
      if (httr::status_code(resp) == 200) {
        txt <- httr::content(resp, as = "text", encoding = "UTF-8")
        out <- jsonlite::fromJSON(txt, simplifyVector = TRUE)
        if (!is.null(cf)) { dir.create(cbio_cache_dir, showWarnings = FALSE, recursive = TRUE); saveRDS(out, cf) }
        return(out)
      }
      last_err <- sprintf("HTTP %d: %s", httr::status_code(resp),
                          substr(httr::content(resp, as = "text", encoding = "UTF-8"), 1, 200))
    }
    warn_msg("  cBioPortal request failed (attempt %d/%d) %s -> %s", attempt, cfg$cbio_retries, path, last_err)
    Sys.sleep(2^attempt)
  }
  stop("cBioPortal request failed: ", path, " (", last_err, ")")
}

json_array <- function(x) jsonlite::toJSON(as.list(x), auto_unbox = TRUE)

cbio_pick_profile <- function(st) {
  profs <- cbio_request(sprintf("/studies/%s/molecular-profiles", st$study_id),
                        cache_key = paste0(st$study_id, "_profiles"))
  if (!is.data.frame(profs) || !nrow(profs)) stop("No molecular profiles returned for ", st$study_id)
  msg("  %s molecular profiles:", st$study_id)
  for (i in seq_len(nrow(profs)))
    msg("    %-55s %-18s %s", profs$molecularProfileId[i], profs$molecularAlterationType[i], profs$name[i])
  if (!is.null(st$profile_id)) {
    if (!st$profile_id %in% profs$molecularProfileId) stop("Requested profile not found: ", st$profile_id)
    return(profs[profs$molecularProfileId == st$profile_id, , drop = FALSE])
  }
  cand <- profs[profs$molecularAlterationType == "MRNA_EXPRESSION" &
                !grepl("zscore", profs$molecularProfileId, ignore.case = TRUE), , drop = FALSE]
  if (!nrow(cand)) stop("No non-z-score MRNA_EXPRESSION profile for ", st$study_id)
  for (pat in st$profile_patterns) {
    hit <- grep(pat, cand$molecularProfileId, ignore.case = TRUE)
    if (length(hit)) return(cand[hit[1], , drop = FALSE])
  }
  cand[1, , drop = FALSE]
}

cbio_pick_sample_list <- function(study_id) {
  sl <- cbio_request(sprintf("/studies/%s/sample-lists", study_id), cache_key = paste0(study_id, "_sample_lists"))
  if (!is.data.frame(sl) || !nrow(sl)) stop("No sample lists for ", study_id)
  pref <- c(paste0(study_id, "_all"), grep("mrna|rna_seq", sl$sampleListId, value = TRUE, ignore.case = TRUE))
  pref <- pref[pref %in% sl$sampleListId]
  if (length(pref)) return(pref[1])
  sl$sampleListId[which.max(sl$sampleCount)]
}

cbio_gene_lookup <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  fetch <- function(syms) {
    if (!length(syms)) return(NULL)
    res <- cbio_request("/genes/fetch", query = list(geneIdType = "HUGO_GENE_SYMBOL", projection = "SUMMARY"),
                        body_json = json_array(syms),
                        cache_key = paste0("genes_", substr(digest_string(paste(sort(syms), collapse = ",")), 1, 12)))
    if (is.data.frame(res) && nrow(res)) res[, c("entrezGeneId", "hugoGeneSymbol")] else NULL
  }
  found <- list()
  for (i in seq(1, length(symbols), by = 500)) {
    found[[length(found) + 1]] <- fetch(symbols[i:min(i + 499, length(symbols))])
  }
  found <- bind_rows(found)
  found$query_symbol <- if (nrow(found)) toupper(found$hugoGeneSymbol) else character(0)

  # retry missing symbols via aliases (and the C#orf# capitalisation rule)
  missing <- setdiff(symbols, found$query_symbol)
  if (length(missing)) {
    alt <- vapply(missing, function(s) {
      if (s %in% names(GENE_ALIASES)) GENE_ALIASES[[s]]
      else if (grepl("^C[0-9XY]+ORF[0-9]+$", s)) sub("ORF", "orf", s)
      else NA_character_
    }, character(1))
    alt <- alt[!is.na(alt) & alt != names(alt)]
    if (length(alt)) {
      res <- fetch(unname(alt))
      if (!is.null(res) && nrow(res)) {
        res$query_symbol <- names(alt)[match(res$hugoGeneSymbol, alt)]
        res$query_symbol[is.na(res$query_symbol)] <- names(alt)[match(toupper(res$hugoGeneSymbol), toupper(alt))][is.na(res$query_symbol)]
        found <- bind_rows(found, res[!is.na(res$query_symbol), ])
      }
    }
  }
  found <- found[!duplicated(found$query_symbol), , drop = FALSE]
  still_missing <- setdiff(symbols, found$query_symbol)
  if (length(still_missing)) msg("  %d symbol(s) not resolvable in cBioPortal: %s",
                                 length(still_missing), paste(still_missing, collapse = ", "))
  found
}

digest_string <- function(x) {
  # lightweight stable hash without extra dependencies
  v <- utf8ToInt(x)
  h <- 0
  for (b in v) h <- (h * 31 + b) %% 2147483647
  sprintf("%d_%d", h, nchar(x))
}

cbio_fetch_expression <- function(profile_id, sample_list_id, gene_tbl) {
  ids <- unique(gene_tbl$entrezGeneId)
  chunks <- split(ids, ceiling(seq_along(ids) / cfg$cbio_gene_batch))
  parts <- list()
  for (k in seq_along(chunks)) {
    body <- jsonlite::toJSON(list(entrezGeneIds = as.list(chunks[[k]]), sampleListId = sample_list_id),
                             auto_unbox = TRUE)
    res <- cbio_request(sprintf("/molecular-profiles/%s/molecular-data/fetch", profile_id),
                        query = list(projection = "SUMMARY"), body_json = body,
                        cache_key = sprintf("%s_%s_expr_%d_of_%d_%s", profile_id, sample_list_id, k,
                                            length(chunks), digest_string(paste(chunks[[k]], collapse = ","))))
    if (is.data.frame(res) && nrow(res)) parts[[k]] <- res[, c("sampleId", "entrezGeneId", "value")]
    msg("  fetched gene batch %d/%d (%d rows)", k, length(chunks), if (is.data.frame(res)) nrow(res) else 0)
  }
  df <- bind_rows(parts)
  if (!nrow(df)) stop("No expression values returned for profile ", profile_id)
  df$value <- suppressWarnings(as.numeric(df$value))
  mat <- tapply(df$value, list(df$entrezGeneId, df$sampleId), function(v) v[1])
  mat <- as.matrix(mat)
  sym <- gene_tbl$query_symbol[match(as.integer(rownames(mat)), gene_tbl$entrezGeneId)]
  keep <- !is.na(sym)
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- sym[keep]
  mat
}

load_bulk_local <- function(path) {
  dt <- data.table::fread(path, data.table = FALSE)
  sym_col <- intersect(c("Hugo_Symbol", "gene", "Gene", "symbol"), colnames(dt))[1]
  if (is.na(sym_col)) stop("Cannot find a gene-symbol column in ", path)
  drop <- intersect(c(sym_col, "Entrez_Gene_Id"), colnames(dt))
  sym  <- toupper(dt[[sym_col]])
  m <- as.matrix(dt[, setdiff(colnames(dt), drop), drop = FALSE])
  suppressWarnings(storage.mode(m) <- "double")
  keep <- !is.na(sym) & sym != "" & !duplicated(sym)
  m <- m[keep, , drop = FALSE]; rownames(m) <- sym[keep]
  m
}

get_bulk_matrix <- function(key, st, genes_needed) {
  local <- cfg$cbio_local_files[[key]]
  if (!is.null(local) && file.exists(local)) {
    msg("  Using local expression file for %s: %s", st$study_id, local)
    m <- load_bulk_local(local)
    return(list(mat = m[rownames(m) %in% genes_needed, , drop = FALSE],
                profile_id = paste0("local:", basename(local)), sample_list = NA_character_,
                profile_name = "local file"))
  }
  prof <- cbio_pick_profile(st)
  sl   <- cbio_pick_sample_list(st$study_id)
  msg("  %s -> profile '%s' (%s), sample list '%s'", st$study_id, prof$molecularProfileId, prof$name, sl)
  genes <- cbio_gene_lookup(genes_needed)
  if (!nrow(genes)) stop("No genes resolved for ", st$study_id)
  m <- cbio_fetch_expression(prof$molecularProfileId, sl, genes)
  list(mat = m, profile_id = prof$molecularProfileId, sample_list = sl, profile_name = prof$name)
}

# ---- 3c. scoring + correlation ---------------------------------------------
to_log2 <- function(mat) {
  mx <- suppressWarnings(max(mat, na.rm = TRUE))
  mn <- suppressWarnings(min(mat, na.rm = TRUE))
  if (is.finite(mx) && mx > 50 && is.finite(mn) && mn >= 0) {
    msg("  values look linear (max %.1f) -> log2(x + 1)", mx); log2(mat + 1)
  } else { msg("  values look already log-scaled (max %.2f); used as is", mx); mat }
}

score_signature <- function(logmat, genes) {
  g <- intersect(unique(toupper(genes)), rownames(logmat))
  if (length(g) < 2) return(list(score = setNames(rep(NA_real_, ncol(logmat)), colnames(logmat)), genes = g))
  sub <- logmat[g, , drop = FALSE]
  sds  <- apply(sub, 1, stats::sd, na.rm = TRUE)
  keep <- rowSums(is.finite(sub)) >= cfg$min_samples_cor & !is.na(sds) & sds > 0
  sub  <- sub[keep, , drop = FALSE]
  if (nrow(sub) < 2) return(list(score = setNames(rep(NA_real_, ncol(logmat)), colnames(logmat)), genes = rownames(sub)))
  z <- t(scale(t(sub)))
  list(score = colMeans(z, na.rm = TRUE), genes = rownames(sub))
}

cor_pair <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < cfg$min_samples_cor) return(c(r = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method))
  c(r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
