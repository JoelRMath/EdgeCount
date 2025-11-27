library(EdgeCount)
library(data.table)


#- added giant-join version 11/26/2025
run_vsea_analysis_gj <- function(object, element_ranks, scoring_statistic = "log2_Anscombe_ratio", n_permutations = 1000, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # --- 1. Setup and Validations ---
  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) stop("No valid elements found in the ECTermScoring object.")

  analysis_object <- reduce_universe_by_elements(object, valid_elements)

  # --- 2. Pre-computation (Mapping to Integer Indices) ---
  # We map term_id and element_id to integers 1..N to save RAM and speed up joins
  dt_edges <- data.table::as.data.table(to_dataframe(analysis_object))
  data.table::setnames(dt_edges, c("term", "element"), c("term_id", "element_id"))

  # Create Integer Maps
  u_terms <- unique(dt_edges$term_id)
  u_elems <- unique(dt_edges$element_id)

  term_map <- setNames(seq_along(u_terms), u_terms)
  elem_map <- setNames(seq_along(u_elems), u_elems)

  # Convert Edges to Integer
  dt_edges[, `:=`(
    t_idx = term_map[term_id],
    e_idx = elem_map[element_id]
  )]

  # Get degrees (vectorized lookup)
  deg_vec <- unlist(analysis_object@ecprob@degrees[u_elems]) # Ordered by u_elems map

  # Term Info
  term_sizes <- lengths(analysis_object@ecprob@adj[u_terms])
  term_k_vec <- term_sizes

  # --- 3. Build the GIANT Permutation Table ---
  # Rank Vector for Real Data (perm = 0)
  rank_vec_real <- rank(element_ranks[u_elems])

  n_elem <- length(u_elems)

  # Generate random ranks efficiently (Matrix approach)
  perm_matrix <- replicate(n_permutations, sample.int(n_elem), simplify = "matrix")

  # Melt to Long Format: (perm_id, e_idx, rank)
  dt_perms <- data.table::data.table(
    perm_id = rep(1:n_permutations, each = n_elem),
    e_idx = rep(1:n_elem, n_permutations),
    rank = as.vector(perm_matrix)
  )

  # Add Real Data as perm_id = 0
  dt_real <- data.table::data.table(perm_id = 0L, e_idx = 1:n_elem, rank = as.integer(rank_vec_real))
  dt_perms <- data.table::rbindlist(list(dt_real, dt_perms))

  # Cleanup to save RAM
  rm(perm_matrix); gc()

  # --- 4. The GIANT Join ---
  data.table::setkey(dt_edges, e_idx)
  data.table::setkey(dt_perms, e_idx)

  # Join Edges with Perms
  # Result: M_edges * (N_perms + 1) rows
  dt_joined <- dt_perms[dt_edges, on = "e_idx", allow.cartesian = TRUE]

  # Add element degrees
  dt_joined[, e_deg := deg_vec[e_idx]]

  # --- 5. The Vectorized Calculation ---
  # Sort by Permutation -> Term -> Rank
  data.table::setorder(dt_joined, perm_id, t_idx, rank)

  global_M2 <- 2 * analysis_object@ecprob@graph_size

  # Running Score Calculation
  dt_joined[, `:=`(
    t_deg = term_k_vec[t_idx],
    observed_ec = 1:.N,
    cumsum_deg = cumsum(e_deg)
  ), by = .(perm_id, t_idx)]

  dt_joined[, lambda := (t_deg / global_M2) * cumsum_deg]

  # Score Formula
  dt_joined[, max_ec := pmin(t_deg, rank)]

  # We calculate the specific statistic requested, but here hardcoded to Anscombe for speed
  dt_joined[, score := 0.5 * (log2(observed_ec + 0.375) - log2(lambda + 0.375))]

  # --- 6. Aggregation ---
  # Calculate Max, Min, Median for every Permutation/Term pair
  summary_dt <- dt_joined[, .(
    min_score = min(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    # Only need extra metadata for the Real run (perm_id == 0)
    rank_at_max = if(perm_id[1] == 0) rank[which.max(score)] else NA_integer_,
    rank_at_min = if(perm_id[1] == 0) rank[which.min(score)] else NA_integer_,
    term_size = t_deg[1]
  ), by = .(perm_id, t_idx)]

  # Cleanup big table
  rm(dt_joined); gc()

  # --- 7. FDR Calculation (Compare Perm 0 to Perms 1..N) ---

  # Helper to calculate stats for a specific metric (Max, Min, or Median)
  calc_metric_stats <- function(metric_col, rank_col = NULL) {
    real_stats <- summary_dt[perm_id == 0, .(t_idx, score = get(metric_col))]
    null_stats <- summary_dt[perm_id > 0, .(t_idx, score = get(metric_col))]

    # Calculate Mean Nulls
    mean_nulls <- null_stats[, .(
      mean_pos = mean(score[score >= 0], na.rm=TRUE),
      mean_neg = mean(abs(score[score < 0]), na.rm=TRUE)
    ), by = t_idx]

    mean_nulls[is.na(mean_pos), mean_pos := 1]
    mean_nulls[is.na(mean_neg), mean_neg := 1]

    # Calculate NES
    res <- merge(real_stats, mean_nulls, by = "t_idx", all.x = TRUE)
    res[, nes := ifelse(score >= 0, score / mean_pos, score / mean_neg)]

    # Add Metadata (Term ID, Size, Ranks)
    res[, term_id := u_terms[t_idx]]

    # Join back original metadata from perm=0
    meta_real <- summary_dt[perm_id == 0]
    res[meta_real, on = "t_idx", term_size := i.term_size]
    if (!is.null(rank_col)) {
      res[meta_real, on = "t_idx", (rank_col) := get(paste0("i.", rank_col))]
    }

    # Calculate Empirical FDR (Non-equi join)
    data.table::setkey(null_stats, t_idx, score)

    # Count nulls strictly better (or equal)
    # Logic:
    # If Score > 0, count nulls >= score
    # If Score < 0, count nulls <= score

    # We split positive and negative handling
    res_pos <- res[score >= 0]
    if (nrow(res_pos) > 0) {
      res_pos[, n_better := null_stats[res_pos, on=.(t_idx, score >= score), .N, by=.EACHI]$N]
    } else {
      res_pos[, n_better := integer(0)]
    }

    res_neg <- res[score < 0]
    if (nrow(res_neg) > 0) {
      res_neg[, n_better := null_stats[res_neg, on=.(t_idx, score <= score), .N, by=.EACHI]$N]
    } else {
      res_neg[, n_better := integer(0)]
    }

    res_final <- rbind(res_pos, res_neg)
    res_final[, fdr_q_value := (n_better + 1) / (n_permutations + 1)]

    # Rename score column back to specific metric name
    data.table::setnames(res_final, "score", metric_col)

    # Sort
    data.table::setorder(res_final, fdr_q_value, -nes)

    # Final Cleanup
    res_final[, c("t_idx", "mean_pos", "mean_neg", "n_better") := NULL]

    # Reorder columns for readability
    cols <- c("term_id", metric_col, "nes", "fdr_q_value", "term_size")
    if (!is.null(rank_col)) cols <- c(cols, rank_col)

    return(res_final[, ..cols])
  }

  # Generate the 3 Summary Tables
  final_max <- calc_metric_stats("max_score", "rank_at_max")
  final_min <- calc_metric_stats("min_score", "rank_at_min")
  final_median <- calc_metric_stats("median_score", NULL)

  return(list(max_score_summary = final_max,
              min_score_summary = final_min,
              median_score_summary = final_median))
}


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
# complexity BENCHMARK
# ---

complexity_benchmark <- function(){
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
}

data("sample_ects")
element_ranks <- setNames(sample(1:length(sample_ects@elements)), sample_ects@elements)
n_perms <- c(10, 100, 1000)
for (n_p in n_perms){
  print(n_p)
  print(system.time(
    run_vsea_analysis_gj(sample_ects, element_ranks, n_permutations = n_p)
  ))
}
