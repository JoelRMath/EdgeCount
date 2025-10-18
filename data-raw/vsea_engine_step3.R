library(EdgeCount)
library(data.table)

run_vsea_analysis <- function(object, element_ranks, scoring_statistic = "log2_Anscombe_ratio", n_permutations = 1000, seed = NULL) {

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

  real_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
  setorder(real_ranks_dt, global_rank)
  real_scores_flat_dt <- score_core(analysis_object, real_ranks_dt, summary_only = FALSE)

  real_summary_lean <- real_scores_flat_dt[, .(
    min_score = min(get(scoring_statistic), na.rm = TRUE),
    max_score = max(get(scoring_statistic), na.rm = TRUE),
    median_score = median(get(scoring_statistic), na.rm = TRUE)
  ), by = term_id]

  # --- STEP 2: Run Permutations (Validated) ---
  message(paste("--- Step 2: Running", n_permutations, "permutations ---"))

  perm_results_list <- replicate(n_permutations, {
    shuffled_ranks <- sample(element_ranks)
    shuffled_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = shuffled_ranks)
    setorder(shuffled_ranks_dt, global_rank)
    score_core(analysis_object, shuffled_ranks_dt, summary_only = TRUE)
  }, simplify = FALSE)

  null_scores_long <- rbindlist(perm_results_list, idcol = "perm_id")

  # --- STEP 3: Calculate NES and FDR ---
  message("--- Step 3: Calculating NES and FDR ---")

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


library(EdgeCount)
library(data.table)

start_time <- Sys.time()
# 1. Load the sample data included with the package
data(sample_ects)

# 2. Create a small, reproducible example object
#    We'll select 3 terms of size 5 and 3 terms of size 10.
set.seed(1) # for reproducibility
term_sizes <- lengths(sample_ects@ecprob@adj[sample_ects@terms])

terms_size_5 <- names(term_sizes[term_sizes == 5])
terms_size_10 <- names(term_sizes[term_sizes == 10])

selected_terms <- c(
  sample(terms_size_5, 3),
  sample(terms_size_10, 3)
)

ects_example <- reduce_universe_by_terms(sample_ects, selected_terms)

# 3. Create a reproducible random ranked list from the example object's universe
element_ranks <- setNames(
  sample(seq_along(ects_example@elements)),
  ects_example@elements
)

# 4. Run the full VSEA analysis
#    We use a small number of permutations to ensure the example runs quickly.
vsea_results <- run_vsea_analysis(
  ects_example,
  element_ranks,
  scoring_statistic = "log2_Anscombe_ratio",
  n_permutations = 100,
  seed = 1
)

# 5. View the top results for enrichment at the top of the list
print(vsea_results$max_score_summary)

end_time <- Sys.time()
time_diff <- as.numeric(end_time - start_time, units = "secs")
print(paste("Time:", round(time_diff, 4), "seconds"))




