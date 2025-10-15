# This script benchmarks the old ("build-then-filter") vs. new ("filter-then-build")
# implementations of the reduce_universe methods for ECTermScoring objects.

library(EdgeCount)
library(data.table)

# ---
# METHOD DEFINITIONS
# ---

# --- "Old" Slow Versions (Build-then-filter) ---
reduce_universe_by_elements_old <- function(object, elements_to_keep) {
  bipartite_edges <- as.data.table(to_dataframe(object))
  reduced_edges <- bipartite_edges[element %in% elements_to_keep]
  return(ECTermScoring(reduced_edges))
}

reduce_universe_by_terms_old <- function(object, terms_to_keep) {
  bipartite_edges <- as.data.table(to_dataframe(object))
  reduced_edges <- bipartite_edges[term %in% terms_to_keep]
  return(ECTermScoring(reduced_edges))
}


# --- "New" Fast Versions (Filter-then-build) ---
reduce_universe_by_elements_fast <- function(object, elements_to_keep) {
  term_adj <- object@ecprob@adj[object@terms]
  filtered_adj <- lapply(term_adj, function(neighbors) {
    intersect(neighbors, elements_to_keep)
  })
  filtered_adj <- filtered_adj[lengths(filtered_adj) > 0]
  if (length(filtered_adj) == 0) {
    return(ECTermScoring(data.frame(term=character(), element=character())))
  }
  reduced_edges <- data.table(
    term = rep(names(filtered_adj), lengths(filtered_adj)),
    element = unlist(filtered_adj, use.names = FALSE)
  )
  return(ECTermScoring(reduced_edges))
}

reduce_universe_by_terms_fast <- function(object, terms_to_keep) {
  terms_to_keep_valid <- intersect(object@terms, terms_to_keep)
  filtered_adj <- object@ecprob@adj[terms_to_keep_valid]
  if (length(filtered_adj) == 0) {
    return(ECTermScoring(data.frame(term=character(), element=character())))
  }
  reduced_edges <- data.table(
    term = rep(names(filtered_adj), lengths(filtered_adj)),
    element = unlist(filtered_adj, use.names = FALSE)
  )
  return(ECTermScoring(reduced_edges))
}


# ---
# MAIN BENCHMARK SCRIPT
# ---

message("--- Loading sample data ---")
data(sample_ects)

# 1. Create data for the benchmark
set.seed(1)

n_simulations <- 100
terms_old <- vector("numeric", n_simulations)
terms_new <- vector("numeric", n_simulations)
terms_identical <- TRUE
elements_old <- vector("numeric", n_simulations)
elements_new <- vector("numeric", n_simulations)
elements_identical <- TRUE


for (i in 1:n_simulations){

  print(i)

  elements_to_keep <- sample(sample_ects@elements, 5000)
  terms_to_keep <- sample(sample_ects@terms, 5000)

  # --- Benchmark for reduce_universe_by_elements ---

  start_time_slow_elem <- Sys.time()
  result_slow_elem <- reduce_universe_by_elements_old(sample_ects, elements_to_keep)
  end_time_slow_elem <- Sys.time()
  time_diff_slow_elem <- as.numeric(end_time_slow_elem - start_time_slow_elem, units = "secs")
  elements_old[i] <- time_diff_slow_elem

  start_time_fast_elem <- Sys.time()
  result_fast_elem <- reduce_universe_by_elements_fast(sample_ects, elements_to_keep)
  end_time_fast_elem <- Sys.time()
  time_diff_fast_elem <- as.numeric(end_time_fast_elem - start_time_fast_elem, units = "secs")
  elements_new[i] <- time_diff_fast_elem

  # Verification
  slow_df_elem <- as.data.table(to_dataframe(result_slow_elem))
  fast_df_elem <- as.data.table(to_dataframe(result_fast_elem))
  setorder(slow_df_elem, term, element)
  setorder(fast_df_elem, term, element)
  comparison_elem <- all.equal(slow_df_elem, fast_df_elem)
  if (!isTRUE(comparison_elem)) {
    elements_identical <- FALSE
  }


  # --- Benchmark for reduce_universe_by_terms ---

  start_time_slow_term <- Sys.time()
  result_slow_term <- reduce_universe_by_terms_old(sample_ects, terms_to_keep)
  end_time_slow_term <- Sys.time()
  time_diff_slow_term <- as.numeric(end_time_slow_term - start_time_slow_term, units = "secs")
  terms_old[i] <- time_diff_slow_term

  start_time_fast_term <- Sys.time()
  result_fast_term <- reduce_universe_by_terms_fast(sample_ects, terms_to_keep)
  end_time_fast_term <- Sys.time()
  time_diff_fast_term <- as.numeric(end_time_fast_term - start_time_fast_term, units = "secs")
  terms_new[i] <- time_diff_fast_term

  # Verification
  slow_df_term <- as.data.table(to_dataframe(result_slow_term))
  fast_df_term <- as.data.table(to_dataframe(result_fast_term))
  setorder(slow_df_term, term, element)
  setorder(fast_df_term, term, element)
  comparison_term <- all.equal(slow_df_term, fast_df_term)
  if (!isTRUE(comparison_term)) {
    terms_identical <- FALSE
  }
}

message("--- t-test for elements --- ")
t <- t.test(elements_new, elements_old)
print(t)
message("--- t-test for terms --- ")
t <- t.test(terms_new, terms_old)
print(t)

message("--- checking identical elements ---")
print(paste("Identical elements:", elements_identical))
message("--- checking identical terms ---")
print(paste("Identical terms:", terms_identical))

