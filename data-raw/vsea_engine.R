library(EdgeCount)
library(data.table)

terms_ecranks_statistics_core_with_summary <- function(object, element_ranks, scoring_statistic = "log2_Anscombe_ratio") {

  # --- STEP 1: Core vectorized calculation to get the "flat" table ---
  # This version assumes the input object and element_ranks have the same universe.
  all_element_degrees <- unlist(object@ecprob@degrees)
  ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
  setorder(ranks_dt, global_rank)
  ranks_dt[, degree := all_element_degrees[element_id]]
  ranks_dt[, cumsum_degrees := cumsum(degree)]

  bipartite_edges <- as.data.table(to_dataframe(object))
  setnames(bipartite_edges, c("term", "element"), c("term_id", "element_id"))
  bipartite_edges[, term_id := as.character(term_id)]

  setkey(bipartite_edges, element_id)
  setkey(ranks_dt, element_id)
  final_dt <- ranks_dt[bipartite_edges, on = "element_id"]

  setorder(final_dt, term_id, global_rank)
  final_dt[, rank_in_term := 1:.N, by = term_id]

  term_sizes <- lengths(object@ecprob@adj[object@terms])
  term_summary <- data.table(
    term_id = names(term_sizes),
    term_size = term_sizes,
    term_degree = term_sizes
  )
  final_dt[term_summary, on = "term_id", `:=`(term_size = i.term_size, term_degree = i.term_degree)]

  final_dt[, `:=`(
    observed_ec = rank_in_term,
    max_ec = pmin(term_size, global_rank),
    lambda = (term_degree / (2 * object@ecprob@graph_size)) * cumsum_degrees
  )]

  final_dt[, `:=`(
    p_value = calculate_p_value(object@ecprob, observed_ec, max_ec, lambda),
    log2_Anscombe_ratio = 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))
  )]

  # --- STEP 2: Perform the lean summary (Direct, high-performance method) ---
  summary_new <- final_dt[, .(
    min_score = min(get(scoring_statistic), na.rm = TRUE),
    max_score = max(get(scoring_statistic), na.rm = TRUE),
    median_score = median(get(scoring_statistic), na.rm = TRUE)
  ), by = term_id]

  # --- STEP 3: Perform summary via the split-then-summarize method for comparison ---
  results_list <- split(final_dt, by = "term_id")

  summary_old_list <- lapply(results_list, function(dt) {
    scores <- dt[[scoring_statistic]]
    # Return a list for rbindlist to combine
    list(
      term_id = dt$term_id[1],
      min_score = min(scores, na.rm = TRUE),
      max_score = max(scores, na.rm = TRUE),
      median_score = median(scores, na.rm = TRUE)
    )
  })
  summary_old <- rbindlist(summary_old_list)

  # Return both for external comparison
  return(list(summary_new = summary_new, summary_old = summary_old))
}


# ---
# SCRIPT TO RUN AND TEST THE NEW FUNCTION
# ---
data(sample_ects)
set.seed(1)

# To test this function, ensure the element_ranks have the same universe
# as the ects object.
ects_test <- reduce_universe_by_terms(sample_ects, sample(sample_ects@terms, 1000))
element_ranks <- setNames(seq_along(ects_test@elements), ects_test@elements)

# Run the new core summary function
vsea_summaries <- terms_ecranks_statistics_core_with_summary(ects_test, element_ranks, "p_value")

# Compare the two outputs
setorder(vsea_summaries$summary_new, term_id)
setorder(vsea_summaries$summary_old, term_id)
if (isTRUE(all.equal(vsea_summaries$summary_new, vsea_summaries$summary_old))){
  print("Succes")
} else {
  print("Failure")
}

