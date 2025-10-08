library(EdgeCount)
library(data.table)

# ---
# GOLD STANDARD: The original implementation, corrected for universe reduction
# ---
terms_ecranks_statistics_old <- function(object, element_ranks) {

  # --- Condition the analysis on the provided ranked list ---
  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) {
    stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
  }
  object <- reduce_universe(object, valid_elements)
  element_ranks <- element_ranks[valid_elements]

  input_terms <- object@terms

  # --- All subsequent calculations now use the CORRECT, reduced object ---
  ecprob <- object@ecprob
  M_g <- ecprob@graph_size
  N_e <- length(object@elements)

  # Prepare the ranked elements and cumulative sum once
  sorted_ranks <- sort(element_ranks)
  ranked_elements <- names(sorted_ranks)

  K <- unlist(ecprob@degrees[ranked_elements])
  K[is.na(K)] <- 0 # Ensure elements not in degrees get a 0
  cumul_sum_K <- cumsum(K)
  # Create a lookup for rank positions for the cumsum
  cumsum_lookup <- setNames(cumul_sum_K, ranked_elements)

  # --- Inner function to score one term ---
  score_one_term <- function(obj, term, element_to_ranks, cumul_sum_lookup, M, N){
    elements_term <- get_neighbors(obj, term)
    if(length(elements_term) == 0) return(NULL)

    # Sort the term's elements by their rank
    ranks_term <- unlist(element_to_ranks[elements_term])
    elements_term_sorted <- names(sort(ranks_term))
    ranks_term_sorted <- sort(ranks_term)

    sz <- length(elements_term_sorted)

    K_term <- obj@degrees[[term]]
    lambda_term <- (K_term / (2*M)) * cumul_sum_lookup[elements_term_sorted]

    observed_ec_term <- 1:sz
    max_ec_term <- ranks_term_sorted

    p_value_term <- calculate_p_value(obj, observed_ec_term, max_ec_term, lambda_term)
    log2_anscombe_ratio_term <- 0.5 * (log2(observed_ec_term + 3/8) - log2(lambda_term + 3/8))
    log2_relative_change_term <- log2(observed_ec_term) - log2(lambda_term)

    df <- data.frame(
      element = elements_term_sorted,
      element_relative_rank = ranks_term_sorted / N,
      lambda = lambda_term,
      observed_edge_count = observed_ec_term,
      max_ec = max_ec_term,
      term_size = sz,
      log2_Anscombe_ratio = log2_anscombe_ratio_term,
      log2_relative_change = log2_relative_change_term,
      p_value = p_value_term,
      stringsAsFactors = FALSE
    )
    return(df)
  }

  all_results_list <- lapply(input_terms, function(term_id_iter) {
    score_one_term(
      obj = ecprob,
      term = term_id_iter,
      element_to_ranks = element_ranks,
      cumul_sum_lookup = cumsum_lookup,
      M = M_g,
      N = N_e
    )
  })
  names(all_results_list) <- input_terms
  all_results_list <- all_results_list[!sapply(all_results_list, is.null)]
  return(all_results_list)
}


# ---
# MAIN DEBUGGING SCRIPT
# ---

# 1. Create a more complex toy example designed to fail
message("--- Creating complex toy example to test universe reduction ---")
te_df <- data.frame(
  term = c("TermA", "TermA", "TermB", "TermB", "TermB", "TermC", "TermD", "TermE", "TermE"),
  element = c("E1", "E3", "E3", "E4", "E5", "E5", "E6", "E7", "E8")
)
ects <- ECTermScoring(te_df)

# Ranked list is a SUBSET of the elements, designed to trigger edge cases.
# It omits E3 and E5, which will change the composition of TermA, TermB, and TermC,
# and it omits E6, which should cause TermD to be dropped entirely.
element_ranks <- c("E1" = 1, "E2" = 2, "E4" = 3, "E7" = 4, "E8" = 5)


# 2. Run the "slow but safe" gold standard
message("\n--- Running 'slow' (original) version ---")
results_slow <- suppressWarnings(terms_ecranks_statistics_old(ects, element_ranks))
message("Output from slow version:")
print(results_slow)


# 3. Run the new, fast S4 method
message("\n--- Running 'fast' vectorized version ---")
results_fast <- suppressWarnings(terms_ecranks_statistics(ects, element_ranks))
message("Output from fast version:")
print(results_fast)


# --- Verification ---
message("\n--- Verifying correctness ---")

# Harmonize the objects for a fair comparison
results_fast_df <- lapply(results_fast, as.data.frame)
results_slow_clean <- lapply(results_slow, function(df) { row.names(df) <- NULL; df })
results_fast_clean <- lapply(results_fast_df, function(df) { row.names(df) <- NULL; df })

# Sort the lists by name before comparing
if (!isTRUE(all.equal(sort(names(results_slow_clean)), sort(names(results_fast_clean))))) {
  comparison <- "The two methods produced results for different sets of terms."
} else {
  sorted_slow <- results_slow_clean[sort(names(results_slow_clean))]
  sorted_fast <- results_fast_clean[sort(names(results_fast_clean))]
  comparison <- all.equal(sorted_slow, sorted_fast)
}


if (isTRUE(comparison)) {
  message("SUCCESS: The new, fast method produces identical results to the original.")
} else {
  message("FAILURE: The outputs of the slow and fast methods are different.")
  print(comparison)
}
