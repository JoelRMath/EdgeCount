library(EdgeCount)
library(data.table)

# ---
# SETUP: Create the benchmark data
# ---

create_benchmark_sets <- function(ects_object, n_sets = 100, set_size = 50) {
  set_size <- min(set_size, length(ects_object@elements))
  set_list <- lapply(1:n_sets, function(i) {
    data.table(
      set_id = paste0("Set_", i),
      element = sample(ects_object@elements, set_size)
    )
  })
  return(rbindlist(set_list))
}


# ---
# GOLD STANDARD: Slow but safe implementation, now returning a list
# ---

terms_ecset_statistics_slow <- function(object, input_sets_dt, lambda_method = "fast") {
  all_set_ids <- unique(input_sets_dt$set_id)
  results_list <- vector("list", length(all_set_ids))

  for (i in seq_along(all_set_ids)) {
    current_set_id <- all_set_ids[i]
    current_element_set <- input_sets_dt[set_id == current_set_id, element]

    single_result_df <- terms_ecset_statistics(
      object,
      element_set = current_element_set,
      lambda_method = lambda_method
    )

    if (!is.null(single_result_df) && nrow(single_result_df) > 0) {
      single_result_dt <- as.data.table(single_result_df)
      single_result_dt[, input_set_id := current_set_id]
      results_list[[i]] <- single_result_dt
    }
  }

  # Combine all results into a single table first
  combined_dt <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

  # Split the combined table into a named list
  if (nrow(combined_dt) > 0) {
    results_as_list <- split(combined_dt, by = "input_set_id")
    return(results_as_list)
  } else {
    return(list())
  }
}


# ---
# HIGH-PERFORMANCE VERSION: With flexible input handling
# ---

terms_ecset_statistics_fast_vectorized <- function(object, input_sets_dt, lambda_method = "fast") {

  if (is.list(input_sets_dt) && !is.data.frame(input_sets_dt)) {
    input_sets_dt <- as.data.table(utils::stack(input_sets_dt))
    setnames(input_sets_dt, c("values", "ind"), c("element", "set_id"))
  }

  bipartite_edges <- as.data.table(to_dataframe(object))
  setnames(bipartite_edges, "term", "term_id")
  bipartite_edges[, term_id := as.character(term_id)]

  input_sets_dt_unique <- unique(input_sets_dt, by = c("set_id", "element"))

  setkey(input_sets_dt_unique, element)
  setkey(bipartite_edges, element)
  all_connections <- bipartite_edges[input_sets_dt_unique, on = "element", nomatch = 0, allow.cartesian = TRUE]

  observed_edges_dt <- all_connections[, .(observed_edges = .N), by = .(input_set_id = set_id, term_id)]

  all_element_degrees <- unlist(object@ecprob@degrees)

  term_degrees <- all_element_degrees[bipartite_edges[, unique(term_id)]]
  term_summary <- data.table(term_id = names(term_degrees), term_degree = term_degrees)

  input_set_summary <- input_sets_dt_unique[,
                                            .(sum_degrees_set = sum(all_element_degrees[element], na.rm = TRUE)),
                                            by = .(input_set_id = set_id)
  ]

  final_dt <- copy(observed_edges_dt)

  final_dt[term_summary, on = "term_id", term_degree := i.term_degree]
  final_dt[input_set_summary, on = "input_set_id", sum_degrees_set := i.sum_degrees_set]

  input_set_sizes <- input_sets_dt_unique[, .(set_size = .N), by = .(input_set_id = set_id)]
  final_dt[input_set_sizes, on = "input_set_id", max_possible_edges := i.set_size]

  final_dt[, lambda := (term_degree * sum_degrees_set) / (2 * object@ecprob@graph_size)]
  final_dt[, p_value := calculate_p_value(object@ecprob, observed_edges, max_possible_edges, lambda)]
  final_dt[, log2_Anscombe_ratio := 0.5 * (log2(observed_edges + 3/8) - log2(lambda + 3/8))]

  results_list <- split(final_dt, by = "input_set_id")

  return(results_list)
}


# ---
# MAIN BENCHMARK SCRIPT
# ---

message("--- Loading sample data ---")
data(sample_ecg)
data(sample_ects)

# 1. Create the random sets to test against
message("--- Creating benchmark data: 100 random sets of 50 elements ---")
benchmark_sets_dt <- create_benchmark_sets(sample_ects, n_sets = 100, set_size = 50)

# 2. Run the "slow but safe" gold standard
message("\n--- Running 'slow but safe' version ---")
start_time_slow <- Sys.time()
results_slow <- terms_ecset_statistics_slow(sample_ects, benchmark_sets_dt)
end_time_slow <- Sys.time()
time_diff_slow <- as.numeric(end_time_slow - start_time_slow, units = "secs")
print(paste("Slow method time:", round(time_diff_slow, 4), "seconds"))

# 3. Run the fast version using a NAMED LIST as input to test the new feature
message("\n--- Running 'fast' vectorized version with list input ---")
# Convert the benchmark data.table to a named list
benchmark_sets_list <- split(benchmark_sets_dt$element, benchmark_sets_dt$set_id)
start_time_fast <- Sys.time()
results_fast <- terms_ecset_statistics_fast_vectorized(sample_ects, benchmark_sets_list)
end_time_fast <- Sys.time()
time_diff_fast <- as.numeric(end_time_fast - start_time_fast, units = "secs")
print(paste("Fast method time:", round(time_diff_fast, 4), "seconds"))


# --- Verification ---
message("\n--- Verifying correctness ---")

# Harmonize column names and order for comparison
harmonize_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if ("log2_Anscombe" %in% names(df)) setnames(df, "log2_Anscombe", "log2_Anscombe_ratio")
  if ("observed_edge_count" %in% names(df)) setnames(df, "observed_edge_count", "observed_edges")

  cols_to_keep <- c("term_id", "observed_edges", "lambda", "p_value", "log2_Anscombe_ratio")

  return(df[, ..cols_to_keep])
}

# Perform a rigorous comparison of the two lists
all_identical <- TRUE
all_set_ids <- names(results_slow)

for(set_id in all_set_ids) {
  slow_dt <- harmonize_df(results_slow[[set_id]])
  fast_dt <- harmonize_df(results_fast[[set_id]])

  setorder(slow_dt, term_id)
  setorder(fast_dt, term_id)

  comparison <- all.equal(slow_dt, fast_dt)
  if (!isTRUE(comparison)) {
    all_identical <- FALSE
    message(paste("FAILURE: Results for set_id", set_id, "are different."))
    print(comparison)
    break
  }
}

if (isTRUE(all_identical)) {
  message("SUCCESS: The fast vectorized method produces identical results to the slow method.")
}

