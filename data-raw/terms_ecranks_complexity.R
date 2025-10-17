library(EdgeCount)
library(data.table)



terms_ecranks_statistics_core <- function(object, element_ranks) {

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

  return(final_dt)
}


# --- Time complexity testing
# ---

data(sample_ects)

# Pre-filter
bipartite_dt <- as.data.table(to_dataframe(sample_ects))
term_sizes <- bipartite_dt[, .N, by = term]
terms_to_keep <- term_sizes[N >= 3, term]
ects_base <- reduce_universe_by_terms(sample_ects, terms_to_keep)
all_available_terms <- ects_base@terms

# Benchmark parameters
term_sample_sizes <- seq(from = 250, to = 5000, by = 250) # (3) Outer loop: Varies T
n_term_sets <- 10       # (2) Intermediate loop: Number of random sets for each T
n_replicates <- 100     # (1) Inner loop: Number of times to time each specific set

# A list to store run results
benchmark_results <- list()

# Simulation
message("\n--- Starting benchmark simulations ---")
set.seed(1) # for reproducibility

for (T in term_sample_sizes) {

  message(paste("\nRunning simulations for T =", T, "terms..."))

  for (j in 1:n_term_sets) {

    # --- Preprocessing (outside the timed block) ---
    terms_subset <- sample(all_available_terms, T)
    ects_subset <- reduce_universe_by_terms(ects_base, terms_subset)
    elements_subset <- ects_subset@elements

    # --- Complexity metrics for this specific run ---
    term_sizes_for_run <- lengths(ects_subset@ecprob@adj[ects_subset@terms])
    term_sizes_loggable <- term_sizes_for_run[term_sizes_for_run > 1]
    complexity_metric_terms <- sum(term_sizes_loggable * log(term_sizes_loggable))

    N_e <- length(elements_subset)
    complexity_metric_elements <- N_e * log(N_e)

    new_complexity_metric <- complexity_metric_terms + complexity_metric_elements

    # --- Inner Loop for Timing ---
    run_times <- replicate(n_replicates, {

      # Simulate a real use case
      shuffled_ranks <- sample(seq_along(elements_subset))
      element_ranks_shuffled <- setNames(shuffled_ranks, elements_subset)

      # Measure the execution time of the core scoring function
      system.time({
        suppressWarnings(terms_ecranks_statistics_core(ects_subset, element_ranks_shuffled))
      })[["elapsed"]]
    })

    # Stores results for this specific set
    benchmark_results <- c(benchmark_results, list(data.table(
      num_terms_T = T,
      num_elements_N = N_e,
      complexity_metric_old = complexity_metric_terms,
      complexity_metric_elements = complexity_metric_elements,
      new_complexity_metric = new_complexity_metric,
      med_runtime_secs = median(run_times)
    )))
  }
}

final_summary <- rbindlist(benchmark_results)


message("\n\n--- Benchmark Summary ---")
print(final_summary)

output_file <- "data-raw/ects_ranks_benchmark_results.tsv"
fwrite(final_summary, output_file, sep = "\t")
message(paste("\nSaved detailed benchmark results to", output_file))

message("\n\n--- Linear Model Summary (Runtime vs. OLD Complexity Metric) ---")
model_old <- lm(med_runtime_secs ~ complexity_metric_old, data = final_summary)
print(summary(model_old))

message("\n\n--- Linear Model Summary (Runtime vs. NEW Complexity Metric) ---")
model_new <- lm(med_runtime_secs ~ new_complexity_metric, data = final_summary)
print(summary(model_new))

message("\n\n--- Multiple Linear Regression Model ---")
model_multiple <- lm(
  med_runtime_secs ~ complexity_metric_old + complexity_metric_elements,
  data = final_summary
)
print(summary(model_multiple))
