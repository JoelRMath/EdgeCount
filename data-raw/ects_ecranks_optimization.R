library(EdgeCount)
library(data.table)

# ---
# GOLD STANDARD: A "flat" data.table implementation for easy verification
# ---
terms_ecranks_statistics_gs <- function(object, element_ranks) {

  # --- Condition the analysis on the provided ranked list ---
  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) {
    stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
  }
  object <- reduce_universe_by_elements(object, valid_elements)
  element_ranks <- rank(element_ranks[valid_elements])

  # --- Pre-computation of core data structures ---

  # 1. Create a data.table of ranked elements with their degrees & cumsum
  all_element_degrees <- unlist(object@ecprob@degrees)
  ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
  setorder(ranks_dt, global_rank)
  ranks_dt[, degree := all_element_degrees[element_id]]
  ranks_dt[, cumsum_degrees := cumsum(degree)]

  # 2. Create the long-form term-element table
  bipartite_edges <- as.data.table(to_dataframe(object))
  setnames(bipartite_edges, c("term", "element"), c("term_id", "element_id"))
  bipartite_edges[, term_id := as.character(term_id)]

  # 3. Join all information into a single flat table
  setkey(bipartite_edges, element_id)
  setkey(ranks_dt, element_id)
  final_dt <- ranks_dt[bipartite_edges, on = "element_id"]

  # 4. Sort by term, then by global rank to calculate rank_in_term
  setorder(final_dt, term_id, global_rank)
  final_dt[, rank_in_term := 1:.N, by = term_id]

  # 5. Add term-level information (size and degree) using joins
  #    In this bipartite context, a term's degree IS its size.
  term_sizes <- lengths(object@ecprob@adj[object@terms])
  term_summary <- data.table(
    term_id = names(term_sizes),
    term_size = term_sizes,
    term_degree = term_sizes  # Assign the same value to both
  )
  final_dt[term_summary, on = "term_id", `:=`(term_size = i.term_size, term_degree = i.term_degree)]

  # 6. Final, fully vectorized statistics calculation
  final_dt[, `:=`(
    observed_ec = rank_in_term,
    max_ec = pmin(term_size, global_rank),
    lambda = (term_degree / (2 * object@ecprob@graph_size)) * cumsum_degrees
  )]

  # NEW: Add the final statistical columns
  final_dt[, `:=`(
    p_value = calculate_p_value(object@ecprob, observed_ec, max_ec, lambda),
    log2_Anscombe_ratio = 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))
  )]

  # 7. Reorder for final output
  setcolorder(final_dt, c("p_value", "log2_Anscombe_ratio", "lambda", "observed_ec",
                          "max_ec", "term_id", "element_id", "term_size", "term_degree",
                          "global_rank", "cumsum_degrees"))
  setorder(final_dt, p_value)
  return(final_dt)
}


# ---
# SCRIPT TO RUN THE FUNCTION
# ---
data(sample_ects)
set.seed(2)

ects <- sample_ects
term_selection_dt <- data.table(term = ects@terms,
                                term_degree = unlist(ects@ecprob@degrees[ects@terms]))
term_selection_dt <- term_selection_dt[term_degree >= 2]
selected_terms <- sample(unlist(term_selection_dt$term), 10)
ects <- reduce_universe_by_terms(ects, selected_terms)
elements <- ects@elements
element_ranks <- setNames(seq_along(elements), elements)

# Run the gold standard function
gs_data <- terms_ecranks_statistics_gs(ects, element_ranks)

# Print the head of the final data.table for inspection
print(head(gs_data, 20))

#' @title Summarize a "flat" gold standard data table
#'
#' @description This function takes the "flat" data table from
#' `terms_ecranks_statistics_gs` and computes summary statistics (min, median, max)
#' for the `log2_Anscombe_ratio` for each term.
#'
#' @param gs_data The data.table output from `terms_ecranks_statistics_gs`.
#'
#' @return A `data.table` with one row per term, containing the summary statistics.
summarize_gs_statistics <- function(gs_data) {

  if (!is.data.table(gs_data) || nrow(gs_data) == 0) {
    return(data.table())
  }

  # This single, vectorized operation groups by term_id and calculates all
  # summary statistics at once.
  summary_dt <- gs_data[, {
    # Find the indices of the min, max, and median score
    idx_min <- which.min(log2_Anscombe_ratio)
    idx_max <- which.max(log2_Anscombe_ratio)
    idx_median <- which.min(abs(log2_Anscombe_ratio - median(log2_Anscombe_ratio, na.rm = TRUE)))

    # Create a list of the summary stats for this term
    .(
      term_size = term_size[1],
      min_score = log2_Anscombe_ratio[idx_min],
      element_at_min = element_id[idx_min],
      observed_ec_at_min = observed_ec[idx_min],
      median_score = log2_Anscombe_ratio[idx_median],
      element_at_median = element_id[idx_median],
      observed_ec_at_median = observed_ec[idx_median],
      max_score = log2_Anscombe_ratio[idx_max],
      element_at_max = element_id[idx_max],
      observed_ec_at_max = observed_ec[idx_max]
    )
  }, by = term_id]

  return(summary_dt)
}


# ---
# SCRIPT TO RUN THE FUNCTION
# ---
# data(sample_ects)
# set.seed(2)
#
# ects_full <- sample_ects
# term_selection_dt <- data.table(term = ects_full@terms, term_degree = unlist(ects_full@ecprob@degrees[ects_full@terms]))
# term_selection_dt <- term_selection_dt[term_degree >= 2]
# selected_terms <- sample(unlist(term_selection_dt$term), 10)
# ects_reduced <- reduce_universe_by_terms(ects_full, selected_terms)
# elements <- ects_reduced@elements
# element_ranks <- setNames(seq_along(elements), elements)
#
# # Run the gold standard function to get the flat table
# gs_data <- terms_ecranks_statistics_gs(ects_reduced, element_ranks)
#
# # Run the new summary function on the result
# gs_summary <- summarize_gs_statistics(gs_data)
#
# # Print the summary table for inspection
# print(gs_summary)

