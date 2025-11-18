library(EdgeCount)
library(data.table)

# ---------------------
# FUNCTION DEFINITIONS
# ---------------------

get_candidate_pairs <- function(ecg, ects) {

  network_edges <- data.table(to_dataframe(ecg))
  setnames(network_edges, c("from", "to"), c("element1", "element2"))
  bipartite_edges <- as.data.table(to_dataframe(ects))
  bipartite_edges[, `:=`(term = as.character(term), element = as.character(element))]
  merged1 <- network_edges[bipartite_edges, on = .(element1 = element), nomatch = 0, allow.cartesian = TRUE]
  setnames(merged1, "term", "term1")
  merged2 <- merged1[bipartite_edges, on = .(element2 = element), nomatch = 0, allow.cartesian = TRUE]
  setnames(merged2, "term", "term2")
  unique_pairs_dt <- merged2[term1 != term2,
                             .(term1 = pmin(term1, term2), term2 = pmax(term1, term2))
  ][, unique(.SD)]
  data.table::setnames(unique_pairs_dt, c("term1", "term2"), c("set1", "set2"))
  return(unique_pairs_dt)
}

#' Gold Standard: Slow but safe for-loop implementation
calculate_between_stats_slow <- function(object, pairs_dt, set_membership_dt) {
  term_to_elements_list <- split(set_membership_dt$element, set_membership_dt$set_id)
  all_element_degrees <- unlist(object@degrees)
  results_list <- vector("list", nrow(pairs_dt))

  for (i in 1:nrow(pairs_dt)) {
    set1_id <- pairs_dt$set1[i]
    set2_id <- pairs_dt$set2[i]
    elements1 <- term_to_elements_list[[set1_id]]
    elements2 <- term_to_elements_list[[set2_id]]
    elements1_d <- setdiff(elements1, elements2)
    size1 <- length(elements1_d)
    elements2_d <- setdiff(elements2, elements1)
    size2 <- length(elements2_d)
    sum_deg1 <- sum(all_element_degrees[elements1_d], na.rm = TRUE)
    sum_deg2 <- sum(all_element_degrees[elements2_d], na.rm = TRUE)
    observed <- get_edge_count_between(object, elements1, elements2)
    lambda <- (sum_deg1 * sum_deg2) / (2 * object@graph_size)
    max_ec <- as.numeric(length(elements1_d)) * as.numeric(length(elements2_d))
    p_val <- calculate_p_value(object, observed, max_ec, lambda)
    anscombe <- 0.5 * (log2(observed + 3/8) - log2(lambda + 3/8))

    results_list[[i]] <- list(
      observed_edges = as.integer(observed),
      lambda = lambda,
      p_value = p_val,
      log2_Anscombe_ratio = anscombe,
      size1 = size1,
      size2 = size2
    )
  }
  results_dt <- rbindlist(results_list)
  return(cbind(pairs_dt, results_dt))
}

run_benchmark <- function(){

  data(sample_ecg)
  data(sample_ects)

  sample_ecp <- ECProb(sample_ecg)
  # candidate_pairs <- get_candidate_pairs(sample_ecg, sample_ects)

  set_membership <- as.data.table(to_dataframe(sample_ects))
  setnames(set_membership, c("term", "element"), c("set_id", "element"))
  candidate_pairs <- get_disjoint_connected_pairs(sample_ecp, set_membership)

  print(nrow(candidate_pairs))
  test_pairs <- head(candidate_pairs, 10000)
  test_pairs[, observed_edges := NULL]

  message("Running 'slow but safe' version...")
  start_time_slow <- Sys.time()
  results_slow <- calculate_between_stats_slow(sample_ecp, test_pairs, set_membership)
  end_time_slow <- Sys.time()
  time_diff_slow <- as.numeric(end_time_slow - start_time_slow, units = "secs")
  print(paste("Slow method time:", round(time_diff_slow, 4), "seconds"))

  message("Running fast vectorized S4 method...")
  start_time_fast <- Sys.time()
  results_fast <- calculate_between_stats_fast_vectorized(sample_ecp, test_pairs, set_membership)
  end_time_fast <- Sys.time()
  time_diff_fast <- as.numeric(end_time_fast - start_time_fast, units = "secs")
  print(paste("Fast method time:", round(time_diff_fast, 4), "seconds"))

  message("\n--- Comparing outputs for identity ---")
  setorder(results_slow, set1, set2)
  setcolorder(results_slow, sort(names(results_slow)))
  setorder(results_fast, set1, set2)
  setcolorder(results_fast, sort(names(results_fast)))

  comparison_result <- all.equal(results_slow, results_fast)

  if (isTRUE(comparison_result)) {
    message("SUCCESS: The outputs of the slow and fast methods are identical.")
  } else {
    message("FAILURE: The outputs are different. Details below:")
    print(comparison_result)
  }
}

save_pairs <- function(){

  data(sample_ecg)
  data(sample_ects)

  sample_ecp <- ECProb(sample_ecg)


  set_membership <- as.data.table(to_dataframe(sample_ects))
  setnames(set_membership, c("term", "element"), c("set_id", "element"))
  print(system.time(
    candidate_pairs <- get_disjoint_connected_pairs(sample_ecp, set_membership))
  )
  candidate_pairs[, observed_edges := NULL]

  print(nrow(candidate_pairs))
  test_pairs <- candidate_pairs

  message("Running fast vectorized S4 method...")
  print(system.time(
    results <- calculate_between_stats_fast_vectorized(sample_ecp, test_pairs, set_membership)
  ))
  setorder(results, -log2_Anscombe_ratio)
  results_red <- results[size1 > 3 & size2 > 3]
  select_dt <- results[(size1 == 6 & size2 == 5 & observed_edges == 21) | (size1 == 5 & size2 == 6 & observed_edges == 21)]
  results_red <- results[1:5000]

  data("sample_term_lookup")
  lookup_dt <- data.table(
    set_id = names(sample_term_lookup),
    term_name = sample_term_lookup
  )

  annotated_results <- copy(results_red)
  annotated_results[lookup_dt, on = .(set1 = set_id), term_name_1 := i.term_name]
  annotated_results[lookup_dt, on = .(set2 = set_id), term_name_2 := i.term_name]

  annotated_select <- copy(select_dt)
  annotated_select[lookup_dt, on = .(set1 = set_id), term_name_1 := i.term_name]
  annotated_select[lookup_dt, on = .(set2 = set_id), term_name_2 := i.term_name]

  fwrite(annotated_results, "data-raw/res/between_examples.tsv", sep = "\t", quote = FALSE)
  fwrite(annotated_select, "data-raw/res/between_select.tsv", sep = "\t", quote = FALSE)
  # print(annotated_results)
}

# --- MAIN ---

# run_benchmark()

save_pairs()
