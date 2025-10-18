library(EdgeCount)
library(data.table)

run_vsea_analysis <- function(object, element_ranks, scoring_statistic = "log2_Anscombe_ratio", n_permutations = 1000, seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- STEP 0: Conditioning analysis on the input universe (ONCE) ---
  message("--- Step 0: Conditioning analysis on the input universe ---")
  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) {
    stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
  }
  analysis_object <- reduce_universe_by_elements(object, valid_elements)
  element_ranks <- rank(element_ranks[valid_elements])

  # --- Pre-calculate constant values once ---
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

  # --- This is the core scoring engine, now as an internal helper ---
  score_core <- function(obj, ranks_dt, summary_only = FALSE) {
    ranks_dt[, degree := all_element_degrees[element_id]]
    ranks_dt[, cumsum_degrees := cumsum(degree)]

    setkey(bipartite_edges, element_id)
    setkey(ranks_dt, element_id)
    final_dt <- ranks_dt[bipartite_edges, on = "element_id"]

    setorder(final_dt, term_id, global_rank)
    final_dt[term_summary, on = "term_id", `:=`(term_size = i.term_size, term_degree = i.term_degree)]

    final_dt[, `:=`(
      observed_ec = 1:.N,
      max_ec = pmin(term_size, global_rank)
    ), by = term_id]

    final_dt[, lambda := (term_degree / (2 * obj@ecprob@graph_size)) * cumsum_degrees]

    final_dt[, `:=`(
      p_value = calculate_p_value(obj@ecprob, observed_ec, max_ec, lambda),
      log2_Anscombe_ratio = 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))
    )]

    if (summary_only) {
      summary_dt <- final_dt[, .(
        min_score = min(get(scoring_statistic), na.rm = TRUE),
        max_score = max(get(scoring_statistic), na.rm = TRUE),
        median_score = median(get(scoring_statistic), na.rm = TRUE)
      ), by = term_id]
      return(summary_dt)
    } else {
      return(final_dt)
    }
  }

  # --- STEP 1: Calculate "Real" Scores (ONCE) ---
  message("--- Step 1: Calculating real enrichment scores ---")
  real_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
  setorder(real_ranks_dt, global_rank)
  real_scores_flat_dt <- score_core(analysis_object, real_ranks_dt, summary_only = FALSE)

  real_summary_lean <- real_scores_flat_dt[, .(
    min_score = min(get(scoring_statistic), na.rm = TRUE),
    max_score = max(get(scoring_statistic), na.rm = TRUE),
    median_score = median(get(scoring_statistic), na.rm = TRUE)
  ), by = term_id]

  # --- STEP 2: Run Permutations using the lean summary mode ---
  message(paste("--- Step 2: Running", n_permutations, "permutations ---"))

  perm_results_list <- replicate(n_permutations, {
    shuffled_ranks <- sample(element_ranks)
    shuffled_ranks_dt <- data.table(element_id = names(element_ranks), global_rank = shuffled_ranks)
    setorder(shuffled_ranks_dt, global_rank)
    score_core(analysis_object, shuffled_ranks_dt, summary_only = TRUE)
  }, simplify = FALSE)

  null_scores_long <- rbindlist(perm_results_list, idcol = "perm_id")

  # --- NEW STEP 3: Calculate NES values for real and null scores ---
  message("--- Step 3: Calculating Normalized Enrichment Scores (NES) ---")

  # This helper function performs the normalization
  calculate_nes <- function(real_summary_dt, null_scores_dt, score_col) {

    # Calculate mean of positive and negative nulls for each term
    mean_nulls <- null_scores_dt[, .(
      mean_pos = mean(get(score_col)[get(score_col) > 0], na.rm = TRUE),
      mean_neg = mean(abs(get(score_col)[get(score_col) < 0]), na.rm = TRUE)
    ), by = term_id]
    mean_nulls[is.nan(mean_pos), mean_pos := 1]
    mean_nulls[is.nan(mean_neg), mean_neg := 1]

    # Calculate NES for real scores
    real_summary_dt[mean_nulls, on="term_id", `:=`(
      mean_null_pos = i.mean_pos,
      mean_null_neg = i.mean_neg
    )]
    real_summary_dt[, nes := ifelse(get(score_col) > 0, get(score_col) / mean_null_pos, get(score_col) / mean_null_neg)]

    # Calculate NES for null scores
    null_scores_dt[mean_nulls, on="term_id", `:=`(
      mean_null_pos = i.mean_pos,
      mean_null_neg = i.mean_neg
    )]
    null_scores_dt[, null_nes := ifelse(get(score_col) > 0, get(score_col) / mean_null_pos, get(score_col) / mean_null_neg)]

    return(list(
      real_nes = real_summary_dt$nes,
      null_nes = null_scores_dt$null_nes
    ))
  }

  nes_max <- calculate_nes(real_summary_lean[, .(term_id, max_score)], null_scores_long, "max_score")
  nes_min <- calculate_nes(real_summary_lean[, .(term_id, min_score)], null_scores_long, "min_score")
  nes_median <- calculate_nes(real_summary_lean[, .(term_id, median_score)], null_scores_long, "median_score")

  # --- DEBUGGING: Return intermediate NES results for inspection ---
  message("--- Pausing before FDR calculation. Returning intermediate NES scores. ---")
  return(list(
    real_nes_max = nes_max$real_nes,
    null_nes_max = nes_max$null_nes,
    real_nes_min = nes_min$real_nes,
    null_nes_min = nes_min$null_nes,
    real_nes_median = nes_median$real_nes,
    null_nes_median = nes_median$null_nes
  ))
}


data(sample_ects)
set.seed(1)

terms <- sample_ects@terms
term_selection_dt <- data.table(term = sample_ects@terms, term_degree = unlist(sample_ects@ecprob@degrees[sample_ects@terms]))
term_selection_dt <- term_selection_dt[term_degree >= 13]
selected_terms <- unlist(term_selection_dt$term)
ects <- reduce_universe_by_terms(sample_ects, selected_terms)

element_ranks <- setNames(sample(seq_along(ects@elements)), ects@elements)

# Run the modified VSEA function
intermediate_results <- run_vsea_analysis(
  ects_test,
  element_ranks,
  n_permutations = 1000,
  seed = 1
)

# --- Verification for Step 2 (NES Calculation) ---
message("\n\n--- Verifying NES Distributions ---")

real_nes <- intermediate_results$real_nes_max
null_nes <- intermediate_results$null_nes_max

message("\n-- Summary of REAL NES (max_score) --")
print(summary(real_nes))

message("\n-- Summary of NULL NES (max_score) --")
print(summary(null_nes))

# Perform a t-test to see if the distributions are significantly different
t_test_result <- t.test(real_nes, null_nes)

message("\n-- T-test comparing real vs. null NES (max_score) --")
print(t_test_result)

if (t_test_result$p.value < 0.05) {
  message("\nFAILURE: The distribution of 'real' NES is significantly different from the null distribution.")
} else {
  message("\nSUCCESS: The 'real' NES distribution is statistically indistinguishable from the null.")
  message("The NES normalization step appears to be working correctly.")
}
