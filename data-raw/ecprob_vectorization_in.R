# library(EdgeCount)
# library(data.table)
#
# # --- A small, self-contained toy example for step-by-step validation ---
#
# # 1. The main graph
# edges_dt <- data.table(p1 = c("E1", "E3", "E6", "E2", "E1"), p2 = c("E2", "E4", "E7", "E4", "E6"))
# ecg <- ECGraph(edges_dt)
# ecp <- ECProb(ecg)
#
# # 2. The list of sets we want to test
# sets_dt <- data.table(set_id = c("SA", "SB", "SC"))
#
# # 3. The membership of elements in those sets
# set_membership_dt <- data.table(
#   set_id = c("SA", "SA", "SA", "SB", "SB", "SB", "SC"),
#   element = c("E1", "E2", "E3", "E3", "E4", "E5", "E6") # E5 is not in the graph
# )
#
# message("--- STEP 1: Calculate Observed Edge Counts ---")
#
# # 1a. Get the canonical edge list from the main graph
# ecg_edges <- data.table(to_dataframe(ecp))
# setnames(ecg_edges, c("from", "to"), c("e1", "e2"))
# ecg_edges[, `:=`(canon1 = pmin(e1, e2), canon2 = pmax(e1, e2))]
# ecg_edges <- unique(ecg_edges, by = c("canon1", "canon2"))
# setkey(ecg_edges, canon1, canon2)
#
# message("\n-- Canonical Edges in Graph --")
# print(ecg_edges)
#
# # 1b. Self-join the set memberships to get all possible intra-set element pairs
# set_membership_dt_unique <- unique(set_membership_dt, by = c("set_id", "element"))
# setkey(set_membership_dt_unique, set_id)
# possible_edges <- set_membership_dt_unique[set_membership_dt_unique, on = "set_id", allow.cartesian = TRUE]
# possible_edges <- possible_edges[element < i.element] # Keep unique pairs only
#
# message("\n-- All Possible Intra-Set Edges --")
# print(possible_edges)
#
# # 1c. Create canonical representation and join with real edges
# possible_edges[, `:=`(canon1 = pmin(element, i.element), canon2 = pmax(element, i.element))]
# observed_edges_long <- ecg_edges[possible_edges, on = .(canon1, canon2), nomatch = 0]
#
# message("\n-- Observed Intra-Set Edges (the intersection) --")
# print(observed_edges_long)
#
# # 1d. Aggregate to get the final count for each set
# observed_edges_dt <- observed_edges_long[, .(observed_edges = .N), by = set_id]
#
# message("\n-- Final Observed Edge Counts per Set --")
# print(observed_edges_dt)
#
#
# message("\n\n--- STEP 2: Calculate Lambda Components ---")
#
# # 2a. Create a lookup table for element degrees
# all_element_degrees_dt <- data.table(
#   element = ecp@names,
#   degree = unlist(ecp@degrees)
# )
# setkey(all_element_degrees_dt, element)
#
# # 2b. Join to get the degree for each element in each set (filtering out unknown elements)
# setkey(set_membership_dt_unique, element)
# sets_with_degrees <- all_element_degrees_dt[set_membership_dt_unique, on = "element", nomatch = 0]
#
# message("\n-- Set Memberships with Degrees (elements not in graph are dropped) --")
# print(sets_with_degrees)
#
# # 2c. Aggregate to get sum of degrees and sum of squared degrees
# term_summary <- sets_with_degrees[, .(
#   sum_of_degrees = sum(degree, na.rm = TRUE),
#   sum_of_sq_degrees = sum(degree^2, na.rm = TRUE),
#   set_size = .N
# ), by = set_id]
#
# message("\n-- Final Summary of Degrees per Set --")
# print(term_summary)
#
#
# message("\n\n--- STEP 3: Final Combination ---")
# # This demonstrates the final join and calculation logic from the main function
#
# final_dt <- copy(sets_dt)
# final_dt[observed_edges_dt, on = "set_id", observed_edges := i.observed_edges]
# final_dt[term_summary, on = "set_id", `:=`(
#   sum_of_degrees = i.sum_of_degrees,
#   sum_of_sq_degrees = i.sum_of_sq_degrees,
#   set_size = i.set_size
# )]
# final_dt[is.na(observed_edges), observed_edges := 0L]
#
# final_dt[, lambda := (sum_of_degrees^2 - sum_of_sq_degrees) / (4 * ecp@graph_size)]
#
# message("\n-- Final Combined Table (before stats calculation) --")
# print(final_dt)
#
# calculate_in_stats_slow <- function(object, sets_to_test, set_membership_dt){
#
#   results_list <- vector("list", length(sets_to_test))
#
#   for (i in seq_along(sets_to_test)) {
#     current_set_id <- sets_to_test[i]
#
#     # Get all elements for the current term and filter for valid ones
#     all_elements_for_set <- set_membership_dt[set_id == current_set_id, element]
#     valid_elements <- unique(all_elements_for_set[all_elements_for_set %in% object@names])
#
#     if (length(valid_elements) < 2) {
#       # Handle sets that are too small for an "in" calculation
#       current_list <- list(
#         set_id = current_set_id, observed_edges = 0L,
#         sum_of_degrees = sum(unlist(object@degrees[valid_elements])),
#         sum_of_sq_degrees = sum(unlist(object@degrees[valid_elements])^2),
#         set_size = length(valid_elements)
#       )
#     } else {
#       observed_ec <- get_edge_count_in(object, valid_elements)
#
#       # Manually calculate lambda components for the "fast" model
#       k <- unlist(object@degrees[valid_elements])
#       sum_k <- sum(k)
#       sum_k_sq <- sum(k^2)
#
#       current_list <- list(
#         set_id = current_set_id,
#         observed_edges = as.integer(observed_ec),
#         sum_of_degrees = sum_k,
#         sum_of_sq_degrees = sum_k_sq,
#         set_size = length(valid_elements)
#       )
#     }
#     results_list[[i]] <- current_list
#   }
#
#   slow_results_dt <- rbindlist(results_list)
#   slow_results_dt[, lambda := (sum_of_degrees^2 - sum_of_sq_degrees) / (4 * object@graph_size)]
#
#   return(slow_results_dt)
# }
#
# # Run the slow version
# slow_results <- calculate_in_stats_slow(ecp, sets_dt$set_id, set_membership_dt_unique)
# message("\n-- Slow Method Result --")
# print(slow_results)
#
#
# # --- Final Comparison ---
# message("\n\n--- Comparing Fast and Slow Results ---")
# # Order both tables to ensure a fair comparison
# setorder(final_dt, set_id)
# setorder(slow_results, set_id)
#
# comparison <- all.equal(
#   final_dt[, .(set_id, observed_edges, lambda)],
#   slow_results[, .(set_id, observed_edges, lambda)]
# )
#
# if (isTRUE(comparison)) {
#   message("SUCCESS: The outputs for observed_edges and lambda are identical.")
# } else {
#   message("FAILURE: The outputs are different.")
#   print(comparison)
# }

library(EdgeCount)
library(data.table)

# ---
# This script validates the final, high-performance vectorized method for
# calculating within-set statistics against a slow, safe, for-loop version.
# ---


# --- Gold Standard: A slow but safe for-loop implementation ---
calculate_in_stats_slow <- function(object, sets_to_test, set_membership_dt){

  results_list <- vector("list", length(sets_to_test))

  for (i in seq_along(sets_to_test)) {
    current_set_id <- sets_to_test[i]

    # Get all elements for the current term and filter for valid ones
    all_elements_for_set <- set_membership_dt[set_id == current_set_id, element]
    valid_elements <- unique(all_elements_for_set[all_elements_for_set %in% object@names])

    # Manually calculate all required statistics
    observed_ec <- if (length(valid_elements) < 2) 0L else get_edge_count_in(object, valid_elements)

    k <- unlist(object@degrees[valid_elements])
    sum_k <- sum(k)
    sum_k_sq <- sum(k^2)
    set_size <- length(valid_elements)

    lambda <- (sum_k^2 - sum_k_sq) / (4 * object@graph_size)
    max_ec <- set_size * (set_size - 1) / 2

    p_val <- calculate_p_value(object, observed_ec, max_ec, lambda)
    anscombe <- 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))

    results_list[[i]] <- list(
      set_id = current_set_id,
      observed_edges = as.integer(observed_ec),
      lambda = lambda,
      p_value = p_val,
      log2_Anscombe_ratio = anscombe
    )
  }

  slow_results_dt <- rbindlist(results_list)
  return(slow_results_dt)
}


# --- Final, High-Performance Vectorized Function ---
calculate_in_stats_fast_vectorized <- function(object, sets_dt, set_membership_dt) {

  # Ensure input is unique
  set_membership_dt <- unique(set_membership_dt, by = c("set_id", "element"))

  # --- Step 1: Calculate Observed Edge Counts ---
  ecg_edges <- data.table(to_dataframe(object))
  setnames(ecg_edges, c("from", "to"), c("e1", "e2"))
  ecg_edges[, `:=`(canon1 = pmin(e1, e2), canon2 = pmax(e1, e2))]
  ecg_edges <- unique(ecg_edges, by = c("canon1", "canon2"))
  setkey(ecg_edges, canon1, canon2)

  setkey(set_membership_dt, set_id)
  possible_edges <- set_membership_dt[set_membership_dt, on = "set_id", allow.cartesian = TRUE]
  possible_edges <- possible_edges[element < i.element]

  if (nrow(possible_edges) > 0) {
    possible_edges[, `:=`(canon1 = pmin(element, i.element), canon2 = pmax(element, i.element))]
    observed_edges_long <- ecg_edges[possible_edges, on = .(canon1, canon2), nomatch = 0]
    observed_edges_dt <- observed_edges_long[, .(observed_edges = .N), by = set_id]
  } else {
    observed_edges_dt <- data.table(set_id = character(), observed_edges = integer())
  }

  # --- Step 2: Calculate Lambda Components ---
  all_element_degrees_dt <- data.table(element = object@names, degree = unlist(object@degrees))
  setkey(all_element_degrees_dt, element)

  setkey(set_membership_dt, element)
  sets_with_degrees <- all_element_degrees_dt[set_membership_dt, on = "element", nomatch = 0]

  term_summary <- sets_with_degrees[, .(
    sum_of_degrees = sum(degree, na.rm = TRUE),
    sum_of_sq_degrees = sum(degree^2, na.rm = TRUE),
    set_size = .N
  ), by = set_id]

  # --- Step 3: Join all results and perform final calculations ---
  final_dt <- copy(sets_dt)

  final_dt[observed_edges_dt, on = "set_id", observed_edges := i.observed_edges]
  final_dt[term_summary, on = "set_id", `:=`(
    sum_of_degrees = i.sum_of_degrees,
    sum_of_sq_degrees = i.sum_of_sq_degrees,
    set_size = i.set_size
  )]

  # Clean NAs from joins for all columns
  cols_to_clean <- c("observed_edges", "sum_of_degrees", "sum_of_sq_degrees", "set_size")
  for (col in cols_to_clean) {
    final_dt[is.na(get(col)), (col) := 0]
  }

  # Calculate final statistics
  final_dt[, lambda := (sum_of_degrees^2 - sum_of_sq_degrees) / (4 * object@graph_size)]
  final_dt[, max_possible_edges := set_size * (set_size - 1) / 2]
  final_dt[, p_value := calculate_p_value(object, observed_edges, max_possible_edges, lambda)]
  final_dt[, log2_Anscombe_ratio := 0.5 * (log2(observed_edges + 3/8) - log2(lambda + 3/8))]

  return(final_dt[, .(set_id, observed_edges, lambda, p_value, log2_Anscombe_ratio)])
}

calculate_in_stats_slow_v0 <- function(object, set_membership_dt){

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

# --- THE BENCHMARK SCRIPT ---
message("--- Loading sample data ---")
data(sample_ecg)
data(sample_ects)

# Prepare the inputs
ecp <- ECProb(sample_ecg)
set_membership_dt <- as.data.table(to_dataframe(sample_ects))
setnames(set_membership_dt, c("term", "element"), c("set_id", "element"))
sets_to_test <- data.table(set_id = unique(set_membership_dt$set_id))

message("\n--- Benchmarking calculate_in_stats ---")

message("Running 'slow but safe' loop version...")
start_time_slow <- Sys.time()
results_slow <- calculate_in_stats_slow(ecp, sets_to_test$set_id, set_membership_dt)
end_time_slow <- Sys.time()
time_diff_slow <- as.numeric(end_time_slow - start_time_slow, units = "secs")
print(paste("Slow method time:", round(time_diff_slow, 4), "seconds"))

message("Running 'slow but safe' v0 loop version...")
start_time_slow <- Sys.time()
results_slow_v0 <- calculate_in_stats_slow_v0(ecp, set_membership_dt)
end_time_slow <- Sys.time()
time_diff_slow <- as.numeric(end_time_slow - start_time_slow, units = "secs")
print(paste("Slow v0 method time:", round(time_diff_slow, 4), "seconds"))

message("Running fast vectorized method...")
start_time_fast <- Sys.time()
results_fast <- calculate_in_stats_fast_vectorized(ecp, sets_to_test, set_membership_dt)
end_time_fast <- Sys.time()
time_diff_fast <- as.numeric(end_time_fast - start_time_fast, units = "secs")
print(paste("Fast method time:", round(time_diff_fast, 4), "seconds"))

message("\n--- Comparing slow and fast outputs for identity ---")
setorder(results_slow, set_id)
setorder(results_fast, set_id)

comparison_result <- all.equal(results_slow, results_fast)

if (isTRUE(comparison_result)) {
  message("SUCCESS: The outputs of the slow and fast methods are identical.")
} else {
  message("FAILURE: The outputs are different. Details below:")
  print(comparison_result)
}

message("\n--- Comparing slow and slow_v0 outputs for identity ---")
setorder(results_slow_v0, set_id)
setorder(results_slow, set_id)

comparison_result <- all.equal(results_slow, results_fast)

if (isTRUE(comparison_result)) {
  message("SUCCESS: The outputs of the slow and slow_v0 methods are identical.")
} else {
  message("FAILURE: The outputs are different. Details below:")
  print(comparison_result)
}

