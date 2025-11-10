library(EdgeCount)
library(data.table)

# ---------------------
# FUNCTION DEFINITIONS
# ---------------------

calculate_in_stats_slow_old <- function(object, set_membership_dt){

  set_membership_dt <- unique(set_membership_dt, by = c("set_id", "element"))

  all_sets_to_test <- unique(set_membership_dt$set_id)
  results_list <- vector("list", length(all_sets_to_test))

  for (i in seq_along(all_sets_to_test)) {

    current_set_id <- all_sets_to_test[i]

    # Get all elements for the current term that are in the graph
    all_elements_for_set <- set_membership_dt[set_id == current_set_id, element]

    valid_elements <- unique(all_elements_for_set[all_elements_for_set %in% object@names])
    set_size <- length(valid_elements)
    max_possible_edges <- set_size * (set_size - 1) / 2

    # --- THE FIX: Ensure both branches create a list with the same 5 items ---
    if (length(valid_elements) < 2) {
      # For sets of size 0 or 1, create a default, 5-item result.
      current_list <- list(
        set_id = current_set_id,
        observed_edge_count = 0L,
        lambda = 0.0,
        p_value = 1.0,
        log2_Anscombe_ratio = 0.0,
        set_size = as.integer(set_size),
        max_possible_edges = max_possible_edges
      )
    } else {
      # For larger sets, calculate the stats and construct a 5-item list.
      observed_ec <- get_edge_count_in(object, valid_elements)
      stats <- edge_count_statistics(object, valid_elements, NULL, observed_ec, "fast")
      stats[["set_size"]] <- set_size
      stats[["max_possible_edges"]] <- max_possible_edges
      current_list <- list(
        set_id = current_set_id,
        observed_edge_count = stats$observed_edge_count,
        lambda = stats$lambda,
        p_value = stats$p_value,
        log2_Anscombe_ratio = stats$log2_Anscombe_ratio,
        set_size = stats$set_size,
        max_possible_edges = stats$max_possible_edges
      )
    }
    results_list[[i]] <- current_list
  }

  results_dt <- rbindlist(results_list)

  # Select and reorder columns to match the fast version's output
  final_dt <- results_dt[, .(
    set_id,
    observed_edge_count = observed_edge_count,
    lambda,
    p_value,
    log2_Anscombe_ratio,
    set_size,
    max_possible_edges
  )]

  return(final_dt)
}
data(sample_ecg)
data(sample_ects)
ecprob <- ECProb(sample_ecg)
set_membership_dt <- as.data.table(to_dataframe(sample_ects))
setnames(set_membership_dt, c("term", "element"), c("set_id", "element"))

sets_to_test <- data.table(set_id = unique(set_membership_dt$set_id))

message("\n--- Benchmarking calculate_in_stats ---")

message("Running 'slow but safe' loop version...")
start_time_slow <- Sys.time()
results_slow <- calculate_in_stats_slow_old(ecprob, set_membership_dt)
end_time_slow <- Sys.time()
time_diff_slow <- as.numeric(end_time_slow - start_time_slow, units = "secs")
print(paste("Slow method time:", round(time_diff_slow, 4), "seconds"))

message("Running fast vectorized S4 method...")
start_time_fast <- Sys.time()
results_fast <- calculate_in_stats_fast_vectorized(ecprob, sets_to_test, set_membership_dt)
end_time_fast <- Sys.time()
time_diff_fast <- as.numeric(end_time_fast - start_time_fast, units = "secs")
print(paste("Fast method time:", round(time_diff_fast, 4), "seconds"))

message("\n--- Comparing outputs for identity ---")
setorder(results_slow, set_id)
setorder(results_fast, set_id)

comparison_result <- all.equal(results_slow, results_fast)

if (isTRUE(comparison_result)) {
  message("SUCCESS: The outputs of the slow and fast methods are identical.")
} else {
  message("FAILURE: The outputs are different. Details below:")
  print(comparison_result)
}

message("\n--- Annotating and saving results ---")

# 1. Convert the named vector lookup to a data.table for joining
lookup_dt <- data.table(
  set_id = names(sample_term_lookup),
  term_name = sample_term_lookup
)

# 2. Join the results with the term names
setkey(results_fast, set_id)
setkey(lookup_dt, set_id)
annotated_results <- lookup_dt[results_fast, on = "set_id"]

# 3. Reorder for a clean output
setcolorder(annotated_results, c("set_id", "term_name", "observed_edge_count", "lambda", "p_value", "log2_Anscombe_ratio", "set_size", "max_possible_edges"))

# 4. Create the directory if it doesn't exist
output_dir <- "data-raw/res"

# 5. Save the file
output_file <- file.path(output_dir, "go_term_in_stats.tsv")
fwrite(annotated_results, output_file, sep = "\t")

message(paste("Successfully saved annotated results to:", output_file))
print(head(annotated_results[order(p_value)]))
