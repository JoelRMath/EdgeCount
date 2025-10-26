library(EdgeCount)
library(data.table)

# ---
# Helper function to get complexity metrics for a given object
# ---
get_complexity_metric <- function(ects, n_permutations) {

  ne <- length(ects@elements)
  ne_log_ne <- ne * log(ne)

  # Get term sizes (K_t) from the reduced object
  kt <- lengths(ects@ecprob@adj[ects@terms])
  kt_log_kt <- kt[kt > 0] * log(kt[kt > 0]) # Ensure log(0) is not computed
  sum_kt_log_kt <- sum(kt_log_kt)

  # The full hypothesized complexity of the *core* loop
  core_complexity <- ne_log_ne + sum_kt_log_kt

  return(list(
    ne_log_ne = ne_log_ne,
    sum_kt_log_kt = sum_kt_log_kt,
    # The total theoretical complexity
    full_complexity = n_permutations * core_complexity
  ))
}

# ---
# SCRIPT TO RUN THE BENCHMARK
# ---

message("--- Loading and preparing base data ---")
data(sample_ects)
set.seed(1)

# 1. Pre-filter the ECTS object to have a more interesting term distribution
all_terms <- sample_ects@terms
term_dt <- data.table(term = all_terms,
                      term_degree = unlist(sample_ects@ecprob@degrees[all_terms]))
term_selection_dt <- term_dt[term_degree >= 3]
selected_terms <- unlist(term_selection_dt$term)
full_ects <- reduce_universe_by_terms(sample_ects, selected_terms)

message(paste("Base object created with", length(full_ects@terms), "terms (size >= 3)."))


# 2. Define the benchmark parameters
term_sample_sizes <- seq(from = 200, to = 3000, by = 200) # (3) Outer loop: Varies T
n_repeats <- 5          # Number of random element rank permutations
n_perms_for_test <- 10  # Fixed number of permutations for the VSEA run

# Pre-sample the term lists
n_to_terms <-lapply(term_sample_sizes, function(n){
  sample(full_ects@terms, n)
})
names(n_to_terms) <- as.character(term_sample_sizes)

# A list to store the aggregated results
benchmark_summary <- list()


# 3. Run the nested simulation loop
message("\n--- Starting benchmark simulations ---")
message(paste("Running with fixed n_permutations =", n_perms_for_test))

for (T_size in term_sample_sizes) {

  message(paste("\nRunning benchmark for T =", T_size, "terms..."))

  # --- Preprocessing (outside the timed block) ---
  # Get the pre-sampled list of terms for this size
  selected_terms <- n_to_terms[[as.character(T_size)]]
  ects_subset <- reduce_universe_by_terms(full_ects, selected_terms)

  # Calculate complexity metrics for this ects_subset
  lst <- get_complexity_metric(ects_subset, n_perms_for_test)

  # --- Inner Loop for Timing ---
  replicate_times <- replicate(n_repeats, {

    # Create a new random element ranking for this repeat
    elements <- ects_subset@elements
    element_ranks <- setNames(sample(1:length(elements)), elements)

    # Time the FULL run_vsea_analysis function
    ptm <- proc.time()
    suppressWarnings(
      run_vsea_analysis(
        ects_subset,
        element_ranks,
        n_permutations = n_perms_for_test,
        seed = NULL # Use a different seed for each run
      )
    )
    time_taken <- (proc.time() - ptm)[["user.self"]] # Get CPU time

    time_taken
  })

  # Store the 'min' (best-case) runtime to be robust to system noise
  benchmark_summary <- c(benchmark_summary, list(data.table(
    num_terms_T = T_size,
    ne_log_ne = lst$ne_log_ne,
    sum_kt_log_kt = lst$sum_kt_log_kt,
    min_run_time = min(replicate_times)
  )))
}

# 4. Combine all individual runs into a single, detailed data.table
final_summary <- rbindlist(benchmark_summary)

# 5. Save the detailed results to a file and print
output_file <- "data-raw/res/vsea_full_complexity.tsv"
fwrite(final_summary, output_file, sep = "\t")
message(paste("\nSaved detailed benchmark results to", output_file))


# 6. Perform and print the linear model summary
message("\n\n--- Linear Model Summary (Run Time vs. Complexity Metrics) ---")
model_multiple <- lm(
  min_run_time ~ ne_log_ne + sum_kt_log_kt,
  data = final_summary
)
print(summary(model_multiple))

# Plot the combined complexity vs. run time
plot(
  final_summary$sum_kt_log_kt + final_summary$ne_log_ne,
  final_summary$min_run_time,
  xlab = "Combined Complexity (N_e*log(N_e) + sum(K_t*log(K_t)))",
  ylab = "Min CPU Time (seconds)",
  main = "VSEA Full Runtime vs. Core Complexity"
)
abline(lm(min_run_time ~ I(ne_log_ne + sum_kt_log_kt), data = final_summary), col = "red")
