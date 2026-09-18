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
  if (!is.null(EXTRA_SIGNATURES[[pop]]))   return(EXTRA_SIGNATURES[[pop]])   # stromal-tier panels first
  if (!is.null(POPULATION_MARKERS[[pop]])) return(POPULATION_MARKERS[[pop]])
  panel <- unname(SINGLER_TO_PANEL[pop])
  if (!is.na(panel)) {
    if (!is.null(POPULATION_MARKERS[[panel]])) return(POPULATION_MARKERS[[panel]])
    if (!is.null(EXTRA_SIGNATURES[[panel]]))   return(EXTRA_SIGNATURES[[panel]])
  }
  NULL
}

build_signatures <- function(recurrent, all_markers) {
  out <- list()
  if (!"tier" %in% colnames(recurrent)) recurrent$tier <- "coarse"
  if (!"tier" %in% colnames(all_markers) && nrow(all_markers)) all_markers$tier <- "coarse"
  for (i in seq_len(nrow(recurrent))) {
    pop <- recurrent$population[i]; tier <- recurrent$tier[i]
    canon <- canonical_signature_for(pop)
    if (is.null(canon)) warn_msg("No canonical panel for population '%s'; consensus markers only", pop)
    if (!is.null(canon)) {
      out[[length(out) + 1]] <- data.frame(population = pop, tier = tier, signature_type = "canonical", gene = canon,
                                           n_datasets = NA_integer_, mean_log2FC = NA_real_)
    }
    cons <- build_consensus_signature(pop, all_markers[all_markers$tier == tier, , drop = FALSE])
    if (!is.null(cons)) {
      out[[length(out) + 1]] <- data.frame(population = pop, tier = tier, signature_type = "consensus_markers",
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

# ============================================================================
# PART 3 — cell-set signatures of one group vs the NEPC signature in bulk cohorts
# ============================================================================
run_part3 <- function(group, recurrent, all_markers) {
  hdr("PART 3 — %s-specific cell-set signatures vs NEPC signature (cBioPortal)", group)
  suffix <- paste0("_", group)
  if (is.null(recurrent) || !nrow(recurrent)) {
    warn_msg("No %s-specific cell sets; nothing to correlate.", group); return(invisible(NULL))
  }
  signatures  <- build_signatures(recurrent, all_markers)
  if (!nrow(signatures)) stop("No signatures could be built for the recurrent populations.")
  write.csv(signatures, file.path(OUT_CROSS, paste0("nepc_recurrent_population_signatures", suffix, ".csv")), row.names = FALSE)
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
              file.path(OUT_CROSS, sprintf("nepc_bulk_log2_expression_%s%s.csv", st$study_id, suffix)), row.names = FALSE)
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
    write.csv(scores, file.path(OUT_CROSS, sprintf("nepc_bulk_signature_scores_%s%s.csv", st$study_id, suffix)), row.names = FALSE)
  
    # scatter plots: population score vs NEPC score
    long <- scores %>% pivot_longer(-c(sample_id, NEPC_score), names_to = "signature", values_to = "score") %>%
      filter(is.finite(score), is.finite(NEPC_score))
    if (nrow(long)) {
      p <- ggplot(long, aes(x = NEPC_score, y = score)) +
        geom_point(alpha = 0.5, size = 1) + geom_smooth(method = "lm", se = TRUE, colour = "firebrick") +
        facet_wrap(~ signature, scales = "free_y") + theme_bw(base_size = 9) +
        labs(title = sprintf("%s: %s-specific cell-set signatures vs NEPC (Beltran custom UP) score", st$label, group),
             x = "NEPC signature score (mean z)", y = "population signature score (mean z)")
      pdf(file.path(OUT_CROSS, sprintf("nepc_nepc_correlation_scatter_%s%s.pdf", st$study_id, suffix)), width = 12, height = 9)
      print(p); dev.off()
    }
  }
  
  if (length(summary_rows)) {
    cor_summary <- bind_rows(summary_rows) %>%
      group_by(study_id, method) %>% mutate(p_adj_BH = p.adjust(p, method = "BH")) %>% ungroup() %>%
      mutate(r = round(r, 4), p = signif(p, 4), p_adj_BH = signif(p_adj_BH, 4)) %>%
      arrange(study_id, method, desc(abs(r)))
    write.csv(cor_summary, file.path(OUT_CROSS, paste0("nepc_nepc_correlation_summary", suffix, ".csv")), row.names = FALSE)
    cat("\nCorrelation of recurrent-population signatures with the NEPC signature:\n")
    print(as.data.frame(cor_summary %>% select(study_id, tier, population, signature_type, method, r, p, p_adj_BH, n_samples)))
  
    cor_wide <- cor_summary %>% filter(method == "pearson") %>%
      select(tier, population, signature_type, study_id, r) %>%
      pivot_wider(names_from = study_id, values_from = r, names_prefix = "pearson_r_")
    write.csv(cor_wide, file.path(OUT_CROSS, paste0("nepc_nepc_correlation_wide", suffix, ".csv")), row.names = FALSE)
  
    # heatmap of correlation coefficients
    hm <- cor_summary %>% mutate(sig = paste0(tier, ": ", population, "\n", signature_type))
    p_hm <- ggplot(hm, aes(x = study_id, y = sig, fill = r)) + geom_tile() +
      geom_text(aes(label = sprintf("%.2f\np=%.1e", r, p)), size = 2.5) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", limits = c(-1, 1)) +
      facet_wrap(~ method) + theme_bw(base_size = 9) + labs(x = NULL, y = NULL, fill = "r")
    pdf(file.path(OUT_CROSS, paste0("nepc_nepc_correlation_heatmap", suffix, ".pdf")), width = 8, height = 2 + 0.5 * length(unique(hm$sig)))
    print(p_hm); dev.off()
  }
  if (length(gene_rows)) {
    gene_cor <- bind_rows(gene_rows) %>%
      group_by(study_id) %>% mutate(p_adj_BH = p.adjust(p, method = "BH")) %>% ungroup() %>%
      mutate(rho = round(rho, 4), p = signif(p, 4), p_adj_BH = signif(p_adj_BH, 4)) %>%
      arrange(study_id, population, signature_type, desc(rho))
    write.csv(gene_cor, file.path(OUT_CROSS, paste0("nepc_nepc_correlation_per_gene", suffix, ".csv")), row.names = FALSE)
  }
  if (length(profiles_used)) {
    write.csv(bind_rows(profiles_used), file.path(OUT_CROSS, paste0("nepc_cbioportal_profiles_used", suffix, ".csv")), row.names = FALSE)
  }
  
  prov <- c(
    sprintf("script_version: %s", SCRIPT_VERSION),
    sprintf("run_date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("group: %s; datasets: %s", group, paste(vapply(group_datasets(group), function(d) d$name, character(1)), collapse = ", ")),
    sprintf("clustering: nfeatures=%d npcs=%d resolution=%g min_features=%d max_mito=%g",
            cfg$nfeatures, cfg$npcs, cfg$resolution, cfg$min_features_cell, cfg$max_mito_pct),
    sprintf("annotation: panel z-score, min_score=%g min_margin=%g; SingleR=%s",
            cfg$ann_min_score, cfg$ann_min_margin, cfg$run_singler),
    sprintf("recurrence: population in > %d %s", cfg$recurrence_min_datasets, cfg$recurrence_count_by),
    sprintf("%s-specific cell sets: %s", group, paste(paste(recurrent$tier, recurrent$population, sep = ":"), collapse = ", ")),
    sprintf("consensus signature: markers in >= %d datasets, top %d by (n_datasets, mean log2FC); excluded %s",
            cfg$consensus_min_datasets, cfg$consensus_top_n, cfg$exclude_gene_regex),
    sprintf("NEPC-overlapping genes excluded from population signatures: %s", cfg$exclude_nepc_genes_from_signatures),
    sprintf("NEPC signature (%d genes): %s", length(NEPC_BELTRAN_CUSTOM_UP), paste(NEPC_BELTRAN_CUSTOM_UP, collapse = ",")),
    "bulk scoring: log2(x+1) if linear -> gene-wise z across samples -> mean z over signature genes",
    "correlation: Pearson and Spearman of population score vs NEPC score; BH-adjusted within study and method",
    if (length(profiles_used)) paste0("cBioPortal profiles: ",
      paste(vapply(profiles_used, function(p) sprintf("%s=%s", p$study_id, p$profile_id), character(1)), collapse = "; ")) else "cBioPortal profiles: none"
  )
  writeLines(prov, file.path(OUT_CROSS, paste0("nepc_PROVENANCE", suffix, ".txt")))
  invisible(if (length(summary_rows)) cor_summary else NULL)
}
