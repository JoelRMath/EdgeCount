library(EdgeCount)
library(data.table)

# ---
# GOLD STANDARD: A "slow but safe" implementation using loops.
# ---
terms_ecranks_statistics_gs_loop <- function(object, element_ranks) {

  # --- Step 1: Condition the analysis on the provided ranked list ---
  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) {
    stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
  }
  object <- reduce_universe_by_elements(object, valid_elements)
  element_ranks <- rank(element_ranks[valid_elements]) # Re-rank from 1 to N_valid

  # --- Step 2: Pre-compute globals from the REDUCED object ---
  ecprob <- object@ecprob
  M_g <- ecprob@graph_size
  N_e <- length(object@elements)
  all_terms <- object@terms

  # --- Step 3: Pre-compute the cumulative degree sum ---
  # Create a data.frame of ranks, sort it, and get the sorted element list
  df <- data.frame(element = names(element_ranks),
                   rank = as.numeric(element_ranks),
                   stringsAsFactors = FALSE)
  df <- df[order(df$rank),]
  ranked_elements_sorted <- df$element

  # Get the degrees for these sorted elements
  K <- unlist(ecprob@degrees[ranked_elements_sorted])
  cumul_sum_K <- cumsum(K)


  # --- Step 4: Loop through each term and calculate its stats ---

  # This is the inner helper function, similar to the original 'score_one_term'
  score_one_term_safe <- function(term_id) {

    # Get elements for this term (from the reduced object)
    elements_in_term <- get_neighbors(ecprob, term_id)

    # This intersect is now redundant but harmless, as both are from the same universe
    elements_in_term <- intersect(elements_in_term, names(element_ranks))

    if (length(elements_in_term) == 0) return(NULL)

    # Get the ranks for these specific elements
    ranks_for_term <- element_ranks[elements_in_term]

    # Sort them to get the running order
    ranks_for_term_sorted <- sort(ranks_for_term)
    elements_term_sorted <- names(ranks_for_term_sorted)

    sz <- length(elements_term_sorted)

    # --- Calculate all vectors ---
    observed_ec_term <- 1:sz
    term_size <- rep(sz, sz)
    max_ec_term <- pmin(term_size, ranks_for_term_sorted)
    term_degree <- rep(ecprob@degrees[[term_id]], sz)

    # Use the ranks to look up the correct cumulative sum
    lambda_term <- (term_degree / (2 * M_g)) * cumul_sum_K[ranks_for_term_sorted]

    p_value_term <- calculate_p_value(ecprob, observed_ec_term, max_ec_term, lambda_term)
    log2_anscombe_term <- 0.5 * (log2(observed_ec_term + 3/8) - log2(lambda_term + 3/8))
    log2_rel_change_term <- log2(observed_ec_term) - log2(lambda_term)

    # Assemble the final data.frame for this term
    df <- data.frame(
      element_id = elements_term_sorted,
      global_rank = element_ranks[elements_term_sorted],
      rank_in_term = rank(ranks_for_term_sorted),
      observed_ec = observed_ec_term,
      max_ec = max_ec_term,
      term_size = term_size,
      lambda = lambda_term,
      p_value = p_value_term,
      log2_Anscombe_ratio = log2_anscombe_term,
      stringsAsFactors = FALSE
    )
    row.names(df) <- NULL
    return(as.data.table(df))
  }

  # Run the loop
  all_results_list <- lapply(all_terms, score_one_term_safe)
  names(all_results_list) <- all_terms

  # Filter out any terms that were empty in the reduced universe
  all_results_list <- all_results_list[!sapply(all_results_list, is.null)]

  return(all_results_list)
}

# ---
# SCRIPT TO RUN THE FUNCTION
# ---
set.seed(3)
te_df <- data.frame(
  term = c("TA", "TA", "TB", "TB", "TB", "TC", "TC", "TC", "TC", "TD", "TD", "TD"),
  element = c("E1", "E2", "E2", "E3", "E4", "E2", "E3", "E5", "E6", "E7", "E8", "E5")
)
ects <- ECTermScoring(te_df)
elements <- c(ects@elements, "E9")
elements <- sample(elements)
element_ranks <- setNames(sample(1:length(elements)), elements)

result_gs <- terms_ecranks_statistics_gs_loop(ects, element_ranks)
result_s4 <- terms_ecranks_statistics(ects, element_ranks)

message("\n--- Verifying correctness ---")

# 1. Check that they produced results for the same terms
if (!isTRUE(all.equal(sort(names(result_gs)), sort(names(result_s4))))) {
  comparison <- "The two methods produced results for different sets of terms."
} else {
  # 2. Sort the lists by name before comparing
  sorted_slow <- result_gs[sort(names(result_gs))]
  sorted_fast <- result_s4[sort(names(result_s4))]

  # 3. Compare the lists directly.
  #    Both are now lists of data.tables with identical columns and order.
  comparison <- all.equal(sorted_slow, sorted_fast)
}

# 4. Print the final result
if (isTRUE(comparison)) {
  message("SUCCESS: The fast S4 method produces identical results to the 'gs_loop' gold standard.")
} else {
  message("FAILURE: The outputs of the slow and fast methods are different.")
  print(comparison)
}

