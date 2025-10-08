library(EdgeCount)
library(data.table)

# ---
# GOLD STANDARD: The original implementation, corrected for universe reduction
# ---
terms_ecranks_statistics_old <- function(object, element_ranks, terms = NULL) {

  # This initial filtering is correct, but not sufficient on its own.
  input_elements <- names(element_ranks)
  input_elements <- input_elements[input_elements %in% object@elements]
  if (length(input_elements) < 1){
    stop("no input element in the ECTErmScoring element universe")
  }

  input_terms <- if (is.null(terms)) object@terms else terms[terms %in% object@terms]

  ecprob <- object@ecprob
  M_g <- ecprob@graph_size
  N_e <- length(object@elements)

  df <- data.frame(elements = names(element_ranks), ranks = as.numeric(unlist(element_ranks)), stringsAsFactors = FALSE)
  df <- df[order(df$ranks),]
  ranked_elements <- unlist(df$elements)

  # --- THE FIX: Robustly handle degrees for elements not in the graph ---
  K <- unlist(ecprob@degrees[ranked_elements])
  K[is.na(K)] <- 0 # Assign degree of 0 to elements not in the object
  # --- END FIX ---

  cumul_sum_K <- cumsum(K)

  score_one_term <- function(obj, term, element_to_ranks, cumul_K, M, N){

    elements_term <- get_neighbors(obj, term)

    # Filter for elements that are actually in the ranked list
    elements_term <- intersect(elements_term, names(element_to_ranks))
    sz <- length(elements_term)
    if(sz == 0) return(NULL)

    size_term <- rep(sz, sz)
    ranks_term <- unlist(element_to_ranks[elements_term])
    df_term <- data.frame(ranks = ranks_term, elements = elements_term)
    df_term <- df_term[order(df_term$ranks),]
    elements_term <- df_term$elements
    ranks_term <- df_term$ranks
    K_term <- rep(obj@degrees[[term]], length(elements_term))
    one_over_2M <- rep((1/(2*M)), length(elements_term))
    lambda_term <- K_term * one_over_2M * cumul_sum_K[ranks_term]
    observed_ec_term <- 1:length(elements_term)
    max_ec_term <- ranks_term
    log2_anscombe_ratio_term <- 0.5 * (log2(observed_ec_term + 3/8) - log2(lambda_term + 3/8))
    log2_relative_change_term <- log2(observed_ec_term) - log2(lambda_term)
    p_value_term <- mapply(calculate_p_value,
                           z = observed_ec_term,
                           m = max_ec_term,
                           lambda = lambda_term,
                           MoreArgs = list(object = obj),
                           SIMPLIFY = TRUE)
    df <- data.frame(element = elements_term,
                     element_relative_rank = ranks_term/N,
                     lambda = lambda_term,
                     observed_edge_count = observed_ec_term,
                     max_ec = max_ec_term,
                     term_size = size_term,
                     log2_Anscombe_ratio = log2_anscombe_ratio_term,
                     log2_relative_change = log2_relative_change_term,
                     p_value = p_value_term)
    return(df)
  }

  all_results_list <- lapply(input_terms, function(term_id_iter) {
    score_one_term(
      obj = ecprob,
      term = term_id_iter,
      element_to_ranks = element_ranks,
      cumul_K = cumul_sum_K,
      M = M_g,
      N = N_e
    )
  })
  names(all_results_list) <- input_terms

  all_results_list <- all_results_list[!sapply(all_results_list, is.null)]

  return(all_results_list)
}

# ---
# MAIN BENCHMARK SCRIPT
# ---

message("--- Loading sample data ---")
data(sample_ecg)
data(sample_ects)

# 1. Create a large, reproducible ranked list for the benchmark
message("--- Creating benchmark ranked list using a subset of sample_ects ---")
set.seed(123)
# Use a subset of elements to properly test the universe reduction logic
ranked_elements_subset <- sample(sample_ects@elements, 5000)
element_ranks <- setNames(seq_along(ranked_elements_subset), ranked_elements_subset)


# 2. Run the "slow but safe" gold standard
message("\n--- Running 'slow' (original) version ---")
start_time_slow <- Sys.time()
results_slow <- suppressWarnings(terms_ecranks_statistics_old(sample_ects, element_ranks))
end_time_slow <- Sys.time()
time_diff_slow <- as.numeric(end_time_slow - start_time_slow, units = "secs")
print(paste("Slow method time:", round(time_diff_slow, 4), "seconds"))


# 3. Run the new, fast S4 method
message("\n--- Running 'fast' vectorized version ---")
start_time_fast <- Sys.time()
results_fast <- suppressWarnings(terms_ecranks_statistics(sample_ects, element_ranks))
end_time_fast <- Sys.time()
time_diff_fast <- as.numeric(end_time_fast - start_time_fast, units = "secs")
print(paste("Fast method time:", round(time_diff_fast, 4), "seconds"))


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

