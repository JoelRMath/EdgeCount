library(EdgeCount)
library(data.table)

# terms_ecranks_statistics_core_with_summary <- function(object, element_ranks, scoring_statistic = "log2_Anscombe_ratio") {
#
#   all_element_degrees <- unlist(object@ecprob@degrees)
#
#   ## From here ....
#   ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
#   setorder(ranks_dt, global_rank)
#   ranks_dt[, degree := all_element_degrees[element_id]]
#   ranks_dt[, cumsum_degrees := cumsum(degree)]
#
#   bipartite_edges <- as.data.table(to_dataframe(object))
#   setnames(bipartite_edges, c("term", "element"), c("term_id", "element_id"))
#   bipartite_edges[, term_id := as.character(term_id)]
#
#   setkey(bipartite_edges, element_id)
#   setkey(ranks_dt, element_id)
#   final_dt <- ranks_dt[bipartite_edges, on = "element_id"]
#
#   setorder(final_dt, term_id, global_rank)
#   final_dt[, rank_in_term := 1:.N, by = term_id]
#
#   term_sizes <- lengths(object@ecprob@adj[object@terms])
#   term_summary <- data.table(
#     term_id = names(term_sizes),
#     term_size = term_sizes,
#     term_degree = term_sizes
#   )
#   final_dt[term_summary, on = "term_id", `:=`(term_size = i.term_size, term_degree = i.term_degree)]
#
#   final_dt[, `:=`(
#     observed_ec = rank_in_term,
#     max_ec = pmin(term_size, global_rank),
#     lambda = (term_degree / (2 * object@ecprob@graph_size)) * cumsum_degrees
#   )]
#
#   final_dt[, `:=`(
#     p_value = calculate_p_value(object@ecprob, observed_ec, max_ec, lambda),
#     log2_Anscombe_ratio = 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))
#   )]
#
#   # --- STEP 2: Perform the lean summary (Direct, high-performance method) ---
#   summary_new <- final_dt[, .(
#     min_score = min(get(scoring_statistic), na.rm = TRUE),
#     max_score = max(get(scoring_statistic), na.rm = TRUE),
#     median_score = median(get(scoring_statistic), na.rm = TRUE)
#   ), by = term_id]
#
#   ## ... To here
#
#   # --- STEP 3: Perform summary via the split-then-summarize method for comparison ---
#   results_list <- split(final_dt, by = "term_id")
#
#   summary_old_list <- lapply(results_list, function(dt) {
#     scores <- dt[[scoring_statistic]]
#     # Return a list for rbindlist to combine
#     list(
#       term_id = dt$term_id[1],
#       min_score = min(scores, na.rm = TRUE),
#       max_score = max(scores, na.rm = TRUE),
#       median_score = median(scores, na.rm = TRUE)
#     )
#   })
#   summary_old <- rbindlist(summary_old_list)
#
#   # Return both for external comparison
#   return(list(summary_new = summary_new, summary_old = summary_old))
# }
#
#
# ---
# SCRIPT TO RUN AND TEST THE NEW FUNCTION
# ---
# data(sample_ects)
# set.seed(1)
#
# # To test this function, ensure the element_ranks have the same universe
# # as the ects object.
# ects_test <- reduce_universe_by_terms(sample_ects, sample(sample_ects@terms, 1000))
# element_ranks <- setNames(seq_along(ects_test@elements), ects_test@elements)
#
# # Run the new core summary function
# vsea_summaries <- terms_ecranks_statistics_core_with_summary(ects_test, element_ranks, "p_value")
#
# # Compare the two outputs
# setorder(vsea_summaries$summary_new, term_id)
# setorder(vsea_summaries$summary_old, term_id)
# if (isTRUE(all.equal(vsea_summaries$summary_new, vsea_summaries$summary_old))){
#   print("Succes")
# } else {
#   print("Failure")
# }


run_vsea_analysis <- function(object, element_ranks, scoring_statistic = "log2_Anscombe_ratio", n_permutations = 1000, seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- STEP 0: Condition the analysis on the input universe (ONCE) ---
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
  # It now assumes the input `ranks` are already correctly ordered.
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
  # Create and sort the data.table for the real ranks ONCE
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
    # Shuffle the ranks, preserving the names
    shuffled_ranks <- sample(element_ranks)
    names(shuffled_ranks) <- names(element_ranks)

    # Create the data.table for this permutation. NOTE: It is NOT sorted.
    shuffled_ranks_dt <- data.table(element_id = names(shuffled_ranks), global_rank = shuffled_ranks)

    score_core(analysis_object, shuffled_ranks_dt, summary_only = TRUE)
  }, simplify = FALSE)

  null_scores_long <- rbindlist(perm_results_list, idcol = "perm_id")

  # --- DEBUGGING: Return intermediate results for inspection ---
  message("--- Pausing before NES/FDR calculation. Returning intermediate scores. ---")
  return(list(
    real_summary = real_summary_lean,
    null_scores = null_scores_long
  ))
}

data(sample_ects)
set.seed(1)

terms <- sample_ects@terms
term_selection_dt <- data.table(term = sample_ects@terms, term_degree = unlist(sample_ects@ecprob@degrees[sample_ects@terms]))
term_selection_dt <- term_selection_dt[term_degree >= 13]
print(term_selection_dt)
selected_terms <- unlist(term_selection_dt$term)
ects <- reduce_universe_by_terms(sample_ects, selected_terms)
elements <- ects@elements
element_ranks <- setNames(seq_along(elements), elements)
intermediate_results <- run_vsea_analysis(ects, element_ranks, scoring_statistic = "log2_Anscombe_ratio", n_permutations = 10, seed = 2)
real_scores <- intermediate_results$real_summary
null_scores <- intermediate_results$null_scores

message("\n-- Summary of REAL max_scores --")
print(summary(real_scores$max_score))

message("\n-- Summary of NULL max_scores --")
print(summary(null_scores$max_score))

# Perform a t-test to see if the distributions are significantly different
t_test_result <- t.test(real_scores$max_score, null_scores$max_score)

message("\n-- T-test comparing real vs. null max_scores --")
print(t_test_result)

if (t_test_result$p.value < 0.05) {
  message("\nFAILURE: The distribution of 'real' scores is significantly different from the null distribution.")
} else {
  message("\nSUCCESS: The 'real' score distribution is statistically similar to the null.")
}

message("\n-- Summary of REAL min_scores --")
print(summary(real_scores$min_score))

message("\n-- Summary of NULL min_scores --")
print(summary(null_scores$min_score))

# Perform a t-test to see if the distributions are significantly different
t_test_result <- t.test(real_scores$min_score, null_scores$min_score)

message("\n-- T-test comparing real vs. null min_scores --")
print(t_test_result)

if (t_test_result$p.value < 0.05) {
  message("\nFAILURE: The distribution of 'real' scores is significantly different from the null distribution.")
} else {
  message("\nSUCCESS: The 'real' score distribution is statistically similar to the null.")
}

message("\n-- Summary of REAL median_scores --")
print(summary(real_scores$median_score))

message("\n-- Summary of NULL median_scores --")
print(summary(null_scores$median_score))

# Perform a t-test to see if the distributions are significantly different
t_test_result <- t.test(real_scores$median_score, null_scores$median_score)

message("\n-- T-test comparing real vs. null median_scores --")
print(t_test_result)

if (t_test_result$p.value < 0.05) {
  message("\nFAILURE: The distribution of 'real' scores is significantly different from the null distribution.")
} else {
  message("\nSUCCESS: The 'real' score distribution is statistically similar to the null.")
}
