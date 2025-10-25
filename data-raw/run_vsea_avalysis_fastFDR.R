library(EdgeCount)
library(data.table)


run_vsea_analysis_fastFDR <- function(object, element_ranks, scoring_statistic = "log2_Anscombe_ratio", n_permutations = 1000, seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- STEP 0 & 1: Conditioning and Core Scoring (Validated) ---
  message("--- Step 0 & 1: Conditioning and calculating real scores ---")
  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) {
    stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
  }
  analysis_object <- reduce_universe_by_elements(object, valid_elements)
  element_ranks <- rank(element_ranks[valid_elements])

  # Pre-calculate constant values
  all_element_degrees <- unlist(analysis_object@ecprob@degrees)
  bipartite_edges <- as.data.table(to_dataframe(analysis_object))
  setnames(bipartite_edges, c("term", "element"), c("term_id", "element_id"))
  bipartite_edges[, term_id := as.character(term_id)]

  term_sizes <- lengths(analysis_object@ecprob@adj[analysis_object@terms])
  term_summary <- data.table(
    term_id = names(term_sizes),
    term_size = term_sizes,
    term_degree = term_sizes
  )

  # Internal core scoring engine
  score_core <- function(obj, ranks_dt, summary_only = FALSE) {
    ranks_dt[, degree := all_element_degrees[element_id]]
    ranks_dt[, cumsum_degrees := cumsum(degree)]
    setkey(bipartite_edges, element_id)
    setkey(ranks_dt, element_id)
    final_dt <- ranks_dt[bipartite_edges, on = "element_id"]
    setorder(final_dt, term_id, global_rank)
    final_dt[term_summary, on = "term_id", `:=`(term_size = i.term_size, term_degree = i.term_degree)]
    final_dt[, `:=`(observed_ec = 1:.N, max_ec = pmin(term_size, global_rank)), by = term_id]
    final_dt[, lambda := (term_degree / (2 * obj@ecprob@graph_size)) * cumsum_degrees]
    final_dt[, `:=`(
      p_value = calculate_p_value(obj@ecprob, observed_ec, max_ec, lambda),
      log2_Anscombe_ratio = 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))
    )]
    if (summary_only) {
      return(final_dt[, .(
        min_score = min(get(scoring_statistic), na.rm = TRUE),
        max_score = max(get(scoring_statistic), na.rm = TRUE),
        median_score = median(get(scoring_statistic), na.rm = TRUE)
      ), by = term_id])
    } else {
      return(final_dt)
    }
  }

  # --- Helper for rich summary (equivalent to table_terms_ecranks_statistics) ---
  summarize_ranks_full <- function(flat_dt) {
    flat_dt[, {
      scores <- get(scoring_statistic)
      idx_min <- which.min(scores)
      idx_max <- which.max(scores)
      idx_median <- which.min(abs(scores - median(scores, na.rm = TRUE)))
      .(
        term_size = term_size[1],
        min_score = scores[idx_min], element_at_min = element_id[idx_min], rank_at_min = global_rank[idx_min],
        median_score = scores[idx_median], element_at_median = element_id[idx_median], rank_at_median = global_rank[idx_median],
        max_score = scores[idx_max], element_at_max = element_id[idx_max], rank_at_max = global_rank[idx_max]
      )
    }, by = term_id]
  }

  # --- STEP 1: Calculate "Real" Scores (ONCE) ---
  message("--- Step 1: Calculating real enrichment scores ---")

  real_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
  setorder(real_ranks_dt, global_rank)
  real_scores_flat_dt <- score_core(analysis_object, real_ranks_dt, summary_only = FALSE)

  real_summary_rich <- summarize_ranks_full(real_scores_flat_dt)
  real_summary_lean <- real_scores_flat_dt[, .(
    min_score = min(get(scoring_statistic), na.rm = TRUE),
    max_score = max(get(scoring_statistic), na.rm = TRUE),
    median_score = median(get(scoring_statistic), na.rm = TRUE)
  ), by = term_id]

  # --- STEP 2: Run Permutations using the lean summary mode ---
  message(paste("--- Step 3: Running", n_permutations, "permutations ---"))

  perm_results_list <- replicate(n_permutations, {
    shuffled_ranks <- sample(element_ranks)
    shuffled_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = shuffled_ranks)
    setorder(shuffled_ranks_dt, global_rank)
    score_core(analysis_object, shuffled_ranks_dt, summary_only = TRUE)
  }, simplify = FALSE)

  null_scores_long <- rbindlist(perm_results_list, idcol = "perm_id")

  # --- STEP 4: Calculate NES and FDR (New High-Performance Version) ---
  message("--- Step 4: Calculating NES and FDR for each summary metric ---")

  calculate_nes_fdr <- function(real_summary_dt, null_scores_dt, score_col) {

    # 1. Calculate normalization factors (mean of nulls)
    mean_nulls <- null_scores_dt[, .(
      mean_pos = mean(get(score_col)[get(score_col) > 0], na.rm = TRUE),
      mean_neg = mean(abs(get(score_col)[get(score_col) < 0]), na.rm = TRUE)
    ), by = term_id]
    mean_nulls[is.nan(mean_pos), mean_pos := 1]
    mean_nulls[is.nan(mean_neg), mean_neg := 1]

    # 2. Calculate NES for real scores
    real_summary_dt[mean_nulls, on="term_id", `:=`(mean_null_pos = i.mean_pos, mean_null_neg = i.mean_neg)]
    real_summary_dt[, nes := ifelse(get(score_col) > 0, get(score_col) / mean_null_pos, get(score_col) / mean_null_neg)]

    # 3. Calculate NES for null scores
    null_scores_dt[mean_nulls, on="term_id", `:=`(mean_null_pos = i.mean_pos, mean_null_neg = i.mean_neg)]
    null_scores_dt[, null_nes := ifelse(get(score_col) > 0, get(score_col) / mean_null_pos, get(score_col) / mean_null_neg)]

    # --- 4. High-Performance FDR Calculation (using findInterval) ---
    all_null_nes <- na.omit(null_scores_dt$null_nes)

    # Pre-sort the null distributions
    null_nes_pos <- sort(all_null_nes[all_null_nes > 0])
    null_nes_neg <- sort(all_null_nes[all_null_nes < 0])

    N_null_pos <- length(null_nes_pos)
    N_null_neg <- length(null_nes_neg)

    # 4a. Calculate frac_real for all scores at once
    real_summary_dt[nes > 0, rank_real := frank(-nes, ties.method = "max")]
    real_summary_dt[nes < 0, rank_real := frank(nes, ties.method = "max")]
    real_summary_dt[nes > 0, frac_real := rank_real / sum(nes > 0)]
    real_summary_dt[nes < 0, frac_real := rank_real / sum(nes < 0)]

    # 4b. Calculate frac_null for all scores at once using fast binary search
    if (N_null_pos > 0) {
      count_null_pos <- N_null_pos - findInterval(real_summary_dt[nes > 0]$nes, null_nes_pos, left.open = TRUE)
      real_summary_dt[nes > 0, frac_null := (count_null_pos + 1) / (N_null_pos + 1)]
    } else {
      real_summary_dt[nes > 0, frac_null := 1.0]
    }

    if (N_null_neg > 0) {
      count_null_neg <- findInterval(real_summary_dt[nes < 0]$nes, null_nes_neg)
      real_summary_dt[nes < 0, frac_null := (count_null_neg + 1) / (N_null_neg + 1)]
    } else {
      real_summary_dt[nes < 0, frac_null := 1.0]
    }

    # 4c. Calculate final FDR and clean up
    real_summary_dt[is.na(frac_real), frac_real := 1.0] # Handle scores of 0
    real_summary_dt[is.na(frac_null), frac_null := 1.0] # Handle scores of 0

    real_summary_dt[, fdr_q_value := frac_null / frac_real]
    real_summary_dt[fdr_q_value > 1, fdr_q_value := 1]

    return(real_summary_dt)
  }

  results_max <- calculate_nes_fdr(real_summary_lean[, .(term_id, max_score)], null_scores_long, "max_score")
  results_min <- calculate_nes_fdr(real_summary_lean[, .(term_id, min_score)], null_scores_long, "min_score")
  results_median <- calculate_nes_fdr(real_summary_lean[, .(term_id, median_score)], null_scores_long, "median_score")

  # --- Step 5: Combine rich summaries with final stats ---
  final_max <- real_summary_rich[results_max, on = "term_id"]
  setorder(final_max, fdr_q_value)

  final_min <- real_summary_rich[results_min, on = "term_id"]
  setorder(final_min, fdr_q_value)

  final_median <- real_summary_rich[results_median, on = "term_id"]
  setorder(final_median, fdr_q_value)

  return(list(
    max_score_summary = final_max,
    min_score_summary = final_min,
    median_score_summary = final_median
  ))
}

run_vsea_analysis_old <-function(object,
                                 element_ranks,
                                 scoring_statistic = "log2_Anscombe_ratio",
                                 n_permutations = 1000,
                                 seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) {
    stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
  }
  analysis_object <- reduce_universe_by_elements(object, valid_elements)
  element_ranks <- rank(element_ranks[valid_elements])

  all_element_degrees <- unlist(analysis_object@ecprob@degrees)
  bipartite_edges <- as.data.table(to_dataframe(analysis_object))
  setnames(bipartite_edges, c("term", "element"), c("term_id", "element_id"))
  bipartite_edges[, term_id := as.character(term_id)]

  term_sizes <- lengths(analysis_object@ecprob@adj[analysis_object@terms])
  term_summary <- data.table(
    term_id = names(term_sizes),
    term_size = term_sizes,
    term_degree = term_sizes
  )

  # Internal scoring engine
  score_core <- function(obj, ranks_dt, summary_only = FALSE) {
    ranks_dt[, degree := all_element_degrees[element_id]]
    ranks_dt[, cumsum_degrees := cumsum(degree)]
    setkey(bipartite_edges, element_id)
    setkey(ranks_dt, element_id)
    final_dt <- ranks_dt[bipartite_edges, on = "element_id"]
    setorder(final_dt, term_id, global_rank)
    final_dt[term_summary, on = "term_id", `:=`(term_size = i.term_size, term_degree = i.term_degree)]
    final_dt[, `:=`(observed_ec = 1:.N, max_ec = pmin(term_size, global_rank)), by = term_id]
    final_dt[, lambda := (term_degree / (2 * obj@ecprob@graph_size)) * cumsum_degrees]
    final_dt[, `:=`(
      p_value = calculate_p_value(obj@ecprob, observed_ec, max_ec, lambda),
      log2_Anscombe_ratio = 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))
    )]
    if (summary_only) {
      return(final_dt[, .(
        min_score = min(get(scoring_statistic), na.rm = TRUE),
        max_score = max(get(scoring_statistic), na.rm = TRUE),
        median_score = median(get(scoring_statistic), na.rm = TRUE)
      ), by = term_id])
    } else {
      return(final_dt)
    }
  }

  # Real scores
  real_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
  setorder(real_ranks_dt, global_rank)
  real_scores_flat_dt <- score_core(analysis_object, real_ranks_dt, summary_only = FALSE)

  real_summary_lean <- real_scores_flat_dt[, .(
    min_score = min(get(scoring_statistic), na.rm = TRUE),
    max_score = max(get(scoring_statistic), na.rm = TRUE),
    median_score = median(get(scoring_statistic), na.rm = TRUE)
  ), by = term_id]

  # Permutations
  perm_results_list <- replicate(n_permutations, {
    shuffled_ranks <- sample(element_ranks)
    shuffled_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = shuffled_ranks)
    setorder(shuffled_ranks_dt, global_rank)
    score_core(analysis_object, shuffled_ranks_dt, summary_only = TRUE)
  }, simplify = FALSE)

  null_scores_long <- rbindlist(perm_results_list, idcol = "perm_id")

  # NES and FDR
  calculate_nes_fdr <- function(real_summary_dt, null_scores_dt, score_col) {
    mean_nulls <- null_scores_dt[, .(
      mean_pos = mean(get(score_col)[get(score_col) > 0], na.rm = TRUE),
      mean_neg = mean(abs(get(score_col)[get(score_col) < 0]), na.rm = TRUE)
    ), by = term_id]
    mean_nulls[is.nan(mean_pos), mean_pos := 1]
    mean_nulls[is.nan(mean_neg), mean_neg := 1]

    real_summary_dt[mean_nulls, on="term_id", `:=`(mean_null_pos = i.mean_pos, mean_null_neg = i.mean_neg)]
    real_summary_dt[, nes := ifelse(get(score_col) > 0, get(score_col) / mean_null_pos, get(score_col) / mean_null_neg)]

    null_scores_dt[mean_nulls, on="term_id", `:=`(mean_null_pos = i.mean_pos, mean_null_neg = i.mean_neg)]
    null_scores_dt[, null_nes := ifelse(get(score_col) > 0, get(score_col) / mean_null_pos, get(score_col) / mean_null_neg)]

    all_null_nes <- na.omit(null_scores_dt$null_nes)
    null_nes_pos <- all_null_nes[all_null_nes > 0]
    null_nes_neg <- all_null_nes[all_null_nes < 0]

    nes_real_vec <- real_summary_dt$nes

    fdr_values <- vapply(nes_real_vec, function(score) {
      if (score > 0) {
        if (length(null_nes_pos) == 0) return(1)
        frac_null <- (sum(null_nes_pos >= score) + 1) / (length(null_nes_pos) + 1)
        frac_real <- sum(nes_real_vec > 0 & nes_real_vec >= score) / sum(nes_real_vec > 0)
        return(if (frac_real == 0) 1 else frac_null / frac_real)
      } else if (score < 0) {
        if (length(null_nes_neg) == 0) return(1)
        frac_null <- (sum(null_nes_neg <= score) + 1) / (length(null_nes_neg) + 1)
        frac_real <- sum(nes_real_vec < 0 & nes_real_vec <= score) / sum(nes_real_vec < 0)
        return(if (frac_real == 0) 1 else frac_null / frac_real)
      } else {
        return(1)
      }
    }, FUN.VALUE = numeric(1))

    fdr_values[fdr_values > 1] <- 1
    real_summary_dt[, fdr_q_value := fdr_values]

    return(real_summary_dt)
  }

  results_max <- calculate_nes_fdr(real_summary_lean[, .(term_id, max_score)], null_scores_long, "max_score")
  results_min <- calculate_nes_fdr(real_summary_lean[, .(term_id, min_score)], null_scores_long, "min_score")
  results_median <- calculate_nes_fdr(real_summary_lean[, .(term_id, median_score)], null_scores_long, "median_score")

  summarize_ranks_full <- function(flat_dt) {
    flat_dt[, {
      scores <- get(scoring_statistic)
      idx_min <- which.min(scores)
      idx_max <- which.max(scores)
      idx_median <- which.min(abs(scores - median(scores, na.rm = TRUE)))
      .(
        term_size = term_size[1],
        min_score = scores[idx_min],
        element_at_min = element_id[idx_min],
        rank_at_min = global_rank[idx_min],
        median_score = scores[idx_median],
        element_at_median = element_id[idx_median],
        rank_at_median = global_rank[idx_median],
        max_score = scores[idx_max],
        element_at_max = element_id[idx_max],
        rank_at_max = global_rank[idx_max]
      )
    }, by = term_id]
  }
  real_summary_rich <- summarize_ranks_full(real_scores_flat_dt)

  final_max <- real_summary_rich[results_max, on = "term_id"]
  setorder(final_max, fdr_q_value)

  final_min <- real_summary_rich[results_min, on = "term_id"]
  setorder(final_min, fdr_q_value)

  final_median <- real_summary_rich[results_median, on = "term_id"]
  setorder(final_median, fdr_q_value)

  return(list(max_score_summary = final_max,
              min_score_summary = final_min,
              median_score_summary = final_median))
}

# --- NEW: Function to run without the FDR calculation ---
run_vsea_analysis_noFDR <- function(object, element_ranks, scoring_statistic = "log2_Anscombe_ratio", n_permutations = 1000, seed = NULL) {

  # This function is identical to run_vsea_analysis_fastFDR,
  # but stops before the NES/FDR calculation.

  if (!is.null(seed)) {
    set.seed(seed)
  }

  message("--- Step 0 & 1: Conditioning and calculating real scores (noFDR) ---")
  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) {
    stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
  }
  analysis_object <- reduce_universe_by_elements(object, valid_elements)
  element_ranks <- rank(element_ranks[valid_elements])

  all_element_degrees <- unlist(analysis_object@ecprob@degrees)
  bipartite_edges <- as.data.table(to_dataframe(analysis_object))
  setnames(bipartite_edges, c("term", "element"), c("term_id", "element_id"))
  bipartite_edges[, term_id := as.character(term_id)]

  term_sizes <- lengths(analysis_object@ecprob@adj[analysis_object@terms])
  term_summary <- data.table(
    term_id = names(term_sizes),
    term_size = term_sizes,
    term_degree = term_sizes
  )

  score_core <- function(obj, ranks_dt, summary_only = FALSE) {
    ranks_dt[, degree := all_element_degrees[element_id]]
    ranks_dt[, cumsum_degrees := cumsum(degree)]
    setkey(bipartite_edges, element_id)
    setkey(ranks_dt, element_id)
    final_dt <- ranks_dt[bipartite_edges, on = "element_id"]
    setorder(final_dt, term_id, global_rank)
    final_dt[term_summary, on = "term_id", `:=`(term_size = i.term_size, term_degree = i.term_degree)]
    final_dt[, `:=`(observed_ec = 1:.N, max_ec = pmin(term_size, global_rank)), by = term_id]
    final_dt[, lambda := (term_degree / (2 * obj@ecprob@graph_size)) * cumsum_degrees]
    final_dt[, `:=`(
      p_value = calculate_p_value(obj@ecprob, observed_ec, max_ec, lambda),
      log2_Anscombe_ratio = 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))
    )]
    if (summary_only) {
      return(final_dt[, .(
        min_score = min(get(scoring_statistic), na.rm = TRUE),
        max_score = max(get(scoring_statistic), na.rm = TRUE),
        median_score = median(get(scoring_statistic), na.rm = TRUE)
      ), by = term_id])
    } else {
      return(final_dt)
    }
  }

  message("--- Step 1: Calculating real enrichment scores (noFDR) ---")
  real_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
  setorder(real_ranks_dt, global_rank)
  real_scores_flat_dt <- score_core(analysis_object, real_ranks_dt, summary_only = FALSE)

  message(paste("--- Step 3: Running", n_permutations, "permutations (noFDR) ---"))
  perm_results_list <- replicate(n_permutations, {
    shuffled_ranks <- sample(element_ranks)
    shuffled_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = shuffled_ranks)
    setorder(shuffled_ranks_dt, global_rank)
    score_core(analysis_object, shuffled_ranks_dt, summary_only = TRUE)
  }, simplify = FALSE)

  null_scores_long <- rbindlist(perm_results_list, idcol = "perm_id")

  message("--- Skipping NES/FDR calculation ---")

  return(null_scores_long) # Return the permutation results
}

run_vsea_analysis_prealloc <- function(object, element_ranks, scoring_statistic = "log2_Anscombe_ratio", n_permutations = 1000, seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- STEP 0 & 1: Conditioning and Core Scoring (Validated) ---
  message("--- Step 0 & 1: Conditioning and calculating real scores ---")
  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) {
    stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
  }
  analysis_object <- reduce_universe_by_elements(object, valid_elements)
  element_ranks <- rank(element_ranks[valid_elements])

  # Pre-calculate constant values
  all_element_degrees <- unlist(analysis_object@ecprob@degrees)
  bipartite_edges <- as.data.table(to_dataframe(analysis_object))
  setnames(bipartite_edges, c("term", "element"), c("term_id", "element_id"))
  bipartite_edges[, term_id := as.character(term_id)]

  term_sizes <- lengths(analysis_object@ecprob@adj[analysis_object@terms])
  term_summary <- data.table(
    term_id = names(term_sizes),
    term_size = term_sizes,
    term_degree = term_sizes
  )

  # Internal core scoring engine (identical to the _fastFDR version)
  score_core <- function(obj, ranks_dt, summary_only = FALSE) {
    ranks_dt[, degree := all_element_degrees[element_id]]
    ranks_dt[, cumsum_degrees := cumsum(degree)]
    setkey(bipartite_edges, element_id)
    setkey(ranks_dt, element_id)
    final_dt <- ranks_dt[bipartite_edges, on = "element_id"]
    setorder(final_dt, term_id, global_rank)
    final_dt[term_summary, on = "term_id", `:=`(term_size = i.term_size, term_degree = i.term_degree)]
    final_dt[, `:=`(observed_ec = 1:.N, max_ec = pmin(term_size, global_rank)), by = term_id]
    final_dt[, lambda := (term_degree / (2 * obj@ecprob@graph_size)) * cumsum_degrees]
    final_dt[, `:=`(
      p_value = calculate_p_value(obj@ecprob, observed_ec, max_ec, lambda),
      log2_Anscombe_ratio = 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))
    )]
    if (summary_only) {
      return(final_dt[, .(
        min_score = min(get(scoring_statistic), na.rm = TRUE),
        max_score = max(get(scoring_statistic), na.rm = TRUE),
        median_score = median(get(scoring_statistic), na.rm = TRUE)
      ), by = term_id])
    } else {
      return(final_dt)
    }
  }

  # --- Helper for rich summary (identical to the _fastFDR version) ---
  summarize_ranks_full <- function(flat_dt) {
    flat_dt[, {
      scores <- get(scoring_statistic)
      idx_min <- which.min(scores)
      idx_max <- which.max(scores)
      idx_median <- which.min(abs(scores - median(scores, na.rm = TRUE)))
      .(
        term_size = term_size[1],
        min_score = scores[idx_min], element_at_min = element_id[idx_min], rank_at_min = global_rank[idx_min],
        median_score = scores[idx_median], element_at_median = element_id[idx_median], rank_at_median = global_rank[idx_median],
        max_score = scores[idx_max], element_at_max = element_id[idx_max], rank_at_max = global_rank[idx_max]
      )
    }, by = term_id]
  }

  # --- STEP 1: Calculate "Real" Scores (ONCE) ---
  message("--- Step 1: Calculating real enrichment scores ---")

  real_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
  setorder(real_ranks_dt, global_rank)
  real_scores_flat_dt <- score_core(analysis_object, real_ranks_dt, summary_only = FALSE)

  real_summary_rich <- summarize_ranks_full(real_scores_flat_dt)
  real_summary_lean <- real_scores_flat_dt[, .(
    min_score = min(get(scoring_statistic), na.rm = TRUE),
    max_score = max(get(scoring_statistic), na.rm = TRUE),
    median_score = median(get(scoring_statistic), na.rm = TRUE)
  ), by = term_id]

  # --- STEP 2: Run Permutations with Pre-Allocation ---
  message(paste("--- Step 3: Running", n_permutations, "permutations (Pre-Allocation) ---"))

  # Pre-allocate the final null_scores table
  all_term_ids <- real_summary_lean$term_id
  n_terms <- length(all_term_ids)
  null_scores_long <- data.table(
    perm_id = rep(1:n_permutations, each = n_terms),
    term_id = rep(all_term_ids, times = n_permutations),
    min_score = NA_real_,
    max_score = NA_real_,
    median_score = NA_real_
  )
  setkey(null_scores_long, perm_id, term_id)

  # Use a for loop instead of replicate
  for (i in 1:n_permutations) {
    shuffled_ranks <- sample(element_ranks)
    shuffled_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = shuffled_ranks)
    setorder(shuffled_ranks_dt, global_rank)

    # Run the core scoring
    perm_summary_dt <- score_core(analysis_object, shuffled_ranks_dt, summary_only = TRUE)

    # Fill the pre-allocated table
    # This join is very fast on keyed tables
    null_scores_long[.(i, perm_summary_dt$term_id), `:=`(
      min_score = perm_summary_dt$min_score,
      max_score = perm_summary_dt$max_score,
      median_score = perm_summary_dt$median_score
    )]
  }

  # --- STEP 4: Calculate NES and FDR (Identical to _fastFDR version) ---
  message("--- Step 4: Calculating NES and FDR for each summary metric ---")

  message("--- Step 4: Calculating NES and FDR for each summary metric ---")

  calculate_nes_fdr <- function(real_summary_dt, null_scores_dt, score_col) {

    # 1. Calculate normalization factors (mean of nulls)
    mean_nulls <- null_scores_dt[, .(
      mean_pos = mean(get(score_col)[get(score_col) > 0], na.rm = TRUE),
      mean_neg = mean(abs(get(score_col)[get(score_col) < 0]), na.rm = TRUE)
    ), by = term_id]
    mean_nulls[is.nan(mean_pos), mean_pos := 1]
    mean_nulls[is.nan(mean_neg), mean_neg := 1]

    # 2. Calculate NES for real scores
    real_summary_dt[mean_nulls, on="term_id", `:=`(mean_null_pos = i.mean_pos, mean_null_neg = i.mean_neg)]
    real_summary_dt[, nes := ifelse(get(score_col) > 0, get(score_col) / mean_null_pos, get(score_col) / mean_null_neg)]

    # 3. Calculate NES for null scores
    null_scores_dt[mean_nulls, on="term_id", `:=`(mean_null_pos = i.mean_pos, mean_null_neg = i.mean_neg)]
    null_scores_dt[, null_nes := ifelse(get(score_col) > 0, get(score_col) / mean_null_pos, get(score_col) / mean_null_neg)]

    # --- 4. High-Performance FDR Calculation (using findInterval) ---
    all_null_nes <- na.omit(null_scores_dt$null_nes)

    # Pre-sort the null distributions
    null_nes_pos <- sort(all_null_nes[all_null_nes > 0])
    null_nes_neg <- sort(all_null_nes[all_null_nes < 0])

    N_null_pos <- length(null_nes_pos)
    N_null_neg <- length(null_nes_neg)

    # 4a. Calculate frac_real for all scores at once
    real_summary_dt[nes > 0, rank_real := frank(-nes, ties.method = "max")]
    real_summary_dt[nes < 0, rank_real := frank(nes, ties.method = "max")]
    real_summary_dt[nes > 0, frac_real := rank_real / sum(nes > 0)]
    real_summary_dt[nes < 0, frac_real := rank_real / sum(nes < 0)]

    # 4b. Calculate frac_null for all scores at once using fast binary search
    if (N_null_pos > 0) {
      count_null_pos <- N_null_pos - findInterval(real_summary_dt[nes > 0]$nes, null_nes_pos, left.open = TRUE)
      real_summary_dt[nes > 0, frac_null := (count_null_pos + 1) / (N_null_pos + 1)]
    } else {
      real_summary_dt[nes > 0, frac_null := 1.0]
    }

    if (N_null_neg > 0) {
      count_null_neg <- findInterval(real_summary_dt[nes < 0]$nes, null_nes_neg)
      real_summary_dt[nes < 0, frac_null := (count_null_neg + 1) / (N_null_neg + 1)]
    } else {
      real_summary_dt[nes < 0, frac_null := 1.0]
    }

    # 4c. Calculate final FDR and clean up
    real_summary_dt[is.na(frac_real), frac_real := 1.0] # Handle scores of 0
    real_summary_dt[is.na(frac_null), frac_null := 1.0] # Handle scores of 0

    real_summary_dt[, fdr_q_value := frac_null / frac_real]
    real_summary_dt[fdr_q_value > 1, fdr_q_value := 1]

    return(real_summary_dt)
  }


  results_max <- calculate_nes_fdr(real_summary_lean[, .(term_id, max_score)], null_scores_long, "max_score")
  results_min <- calculate_nes_fdr(real_summary_lean[, .(term_id, min_score)], null_scores_long, "min_score")
  results_median <- calculate_nes_fdr(real_summary_lean[, .(term_id, median_score)], null_scores_long, "median_score")

  # --- Step 5: Combine rich summaries with final stats ---
  final_max <- real_summary_rich[results_max, on = "term_id"]
  setorder(final_max, fdr_q_value)

  final_min <- real_summary_rich[results_min, on = "term_id"]
  setorder(final_min, fdr_q_value)

  final_median <- real_summary_rich[results_median, on = "term_id"]
  setorder(final_median, fdr_q_value)

  return(list(
    max_score_summary = final_max,
    min_score_summary = final_min,
    median_score_summary = final_median
  ))
}

# ---
# MAIN BENCHMARK SCRIPT
# ---

data("sample_ects")
term_size <- 10
term_dt <- data.table(term = sample_ects@terms,
                      term_degree = unlist(sample_ects@ecprob@degrees[sample_ects@terms]))
term_selection_dt <- term_dt[term_degree >= term_size]
selected_terms <- unlist(term_selection_dt$term)
ects <- reduce_universe_by_terms(sample_ects, selected_terms)

elements <- ects@elements
element_ranks <- setNames(sample(1:length(elements)), elements)

message("-- running fast FDR version (replicate) --")
start_time <- Sys.time()
results_fast <- run_vsea_analysis_fastFDR(ects,
                                          element_ranks,
                                          scoring_statistic = "log2_Anscombe_ratio",
                                          n_permutations = 1000,
                                          seed = NULL)
end_time <- Sys.time()
taken_time <- as.numeric(end_time - start_time, units = "secs")
message("-- taken time (replicate) = ", taken_time)

message("-- running fast FDR version (pre-alloc) --")
start_time <- Sys.time()
results_prealloc <- run_vsea_analysis_prealloc(ects,
                                               element_ranks,
                                               scoring_statistic = "log2_Anscombe_ratio",
                                               n_permutations = 1000,
                                               seed = NULL)
end_time <- Sys.time()
taken_time <- as.numeric(end_time - start_time, units = "secs")
message("-- taken time (pre-alloc) = ", taken_time)

message("-- running old FDR version (vapply) --")
start_time <- Sys.time()
results_old <- run_vsea_analysis_old(ects,
                                     element_ranks,
                                     scoring_statistic = "log2_Anscombe_ratio",
                                     n_permutations = 1000,
                                     seed = NULL)
end_time <- Sys.time()
taken_time <- as.numeric(end_time - start_time, units = "secs")
message("-- taken time (vapply) = ", taken_time)

message("-- running no FDR version --")
start_time <- Sys.time()
results_nofdr <- run_vsea_analysis_noFDR(ects,
                                         element_ranks,
                                         scoring_statistic = "log2_Anscombe_ratio",
                                         n_permutations = 1000,
                                         seed = NULL)
end_time <- Sys.time()
taken_time <- as.numeric(end_time - start_time, units = "secs")
message("-- taken time (no FDR) = ", taken_time)
