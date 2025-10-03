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

    # --- THE FIX: Ensure both branches create a list with the same 5 items ---
    if (length(valid_elements) < 2) {
      # For sets of size 0 or 1, create a default, 5-item result.
      current_list <- list(
        set_id = current_set_id,
        observed_edge_count = 0L,
        lambda = 0.0,
        p_value = 1.0,
        log2_Anscombe_ratio = 0.0
      )
    } else {
      # For larger sets, calculate the stats and construct a 5-item list.
      observed_ec <- get_edge_count_in(object, valid_elements)
      stats <- edge_count_statistics(object, valid_elements, NULL, observed_ec, "fast")

      current_list <- list(
        set_id = current_set_id,
        observed_edge_count = stats$observed_edge_count,
        lambda = stats$lambda,
        p_value = stats$p_value,
        log2_Anscombe_ratio = stats$log2_Anscombe_ratio
      )
    }
    results_list[[i]] <- current_list
  }

  results_dt <- rbindlist(results_list)

  # Select and reorder columns to match the fast version's output
  final_dt <- results_dt[, .(
    set_id,
    observed_edges = observed_edge_count,
    lambda,
    p_value,
    log2_Anscombe_ratio
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
slow_results <- calculate_in_stats_slow(ecprob, set_membership_dt)
end_time_slow <- Sys.time()
time_diff_slow <- as.numeric(end_time_slow - start_time_slow, units = "secs")
print(paste("Slow method time:", round(time_diff_slow, 4), "seconds"))

message("Running fast vectorized S4 method...")
start_time_fast <- Sys.time()
results_fast <- calculate_in_stats_fast_vectorized(ecp, sets_to_test, set_membership_dt)
end_time_fast <- Sys.time()
time_diff_fast <- as.numeric(end_time_fast - start_time_fast, units = "secs")
print(paste("Fast method time:", round(time_diff_fast, 4), "seconds"))

message("\n--- Comparing outputs for identity ---")
# setorder(results_slow, set_id)
# setorder(results_fast, set_id)

comparison_result <- all.equal(results_slow, results_fast)

if (isTRUE(comparison_result)) {
  message("SUCCESS: The outputs of the slow and fast methods are identical.")
} else {
  message("FAILURE: The outputs are different. Details below:")
  print(comparison_result)
}
