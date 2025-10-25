library(EdgeCount)
library(data.table)

# ---
# SCRIPT TO RUN THE BENCHMARK
# ---

# 1. Load and prepare a FIXED base dataset
message("--- Loading and preparing base data ---")
data(sample_ects)

# We will use a fixed, reasonably large subset of terms for all runs
set.seed(42)
terms_to_keep <- sample(sample_ects@terms, 3000) # Use 3000 terms
ects_benchmark <- reduce_universe_by_terms(sample_ects, terms_to_keep)
elements_benchmark <- ects_benchmark@elements
element_ranks_benchmark <- setNames(seq_along(elements_benchmark), elements_benchmark)

message(paste("Benchmark will run on a fixed dataset with",
              length(ects_benchmark@terms), "terms and",
              length(elements_benchmark), "elements."))


# 2. Define the benchmark parameters
# We will vary the number of permutations
permutation_steps <- c(100, 300, 500, 700, 1000)

# A list to store the results
benchmark_results <- list()

# 3. Run the simulation loop
message("\n--- Starting benchmark simulations ---")

for (n_perms in permutation_steps) {

  message(paste("Running VSEA with n_permutations =", n_perms, "..."))

  # Measure the execution time of the full function
  time_taken <- system.time({
    suppressWarnings(
      run_vsea_analysis_fastFDR( # Use the fastest function
        object = ects_benchmark,
        element_ranks = element_ranks_benchmark,
        n_permutations = n_perms,
        seed = 42
      )
    )
  })[["elapsed"]]

  # Store the results for this step
  benchmark_results <- c(benchmark_results, list(data.table(
    n_permutations = n_perms,
    runtime_secs = time_taken
  )))
}

# 4. Combine and display the final results
final_summary <- rbindlist(benchmark_results)

message("\n\n--- Benchmark Summary ---")
print(final_summary)

# 5. Perform and print the linear model summary
message("\n\n--- Linear Model Summary (Runtime vs. Number of Permutations) ---")
model <- lm(runtime_secs ~ n_permutations, data = final_summary)
print(summary(model))
