library(EdgeCount)
library(data.table)

# Internal helper: calculates stats but keeps result as a flat data.table
# to_dataframe/edges extraction happens OUTSIDE this function.
calc_ecset_stats_core_proto <- function(object, input_sets_dt, bipartite_edges_dt) {

  # 1. Set Keys for fast joining
  # input_sets_dt must have: set_id, element
  # bipartite_edges_dt must have: term_id, element

  # Ensure uniqueness of input pairs
  input_sets_dt_unique <- unique(input_sets_dt, by = c("set_id", "element"))
  setkey(input_sets_dt_unique, element)
  setkey(bipartite_edges_dt, element)

  # 2. The Big Join (Cartesian product of Sets * Terms via Elements)
  all_connections <- bipartite_edges_dt[input_sets_dt_unique, on = "element", nomatch = 0, allow.cartesian = TRUE]

  # 3. Aggregation (Observed Edge Count)
  # Note: input_set_id comes from the 'i' table (input_sets_dt)
  observed_edges_dt <- all_connections[, .(observed_edge_count = .N), by = .(input_set_id = set_id, term_id)]

  # 4. Pre-computation of Degrees (Global properties)
  all_element_degrees <- unlist(object@ecprob@degrees)

  # Get degrees only for terms actually involved in connections
  term_degrees <- all_element_degrees[bipartite_edges_dt[, unique(term_id)]]
  term_summary <- data.table(term_id = names(term_degrees), term_degree = term_degrees)

  # 5. Set Summaries (Sum of degrees of elements in the set)
  valid_input_elements <- all_connections[, .(input_set_id = set_id, element)] |> unique()
  input_set_summary <- valid_input_elements[,
                                            .(sum_degrees_set = sum(all_element_degrees[element], na.rm = TRUE),
                                              set_size = .N),
                                            by = input_set_id
  ]

  # 6. Merge and Calculate Statistics
  final_dt <- copy(observed_edges_dt)
  final_dt[term_summary, on = "term_id", term_degree := i.term_degree]
  final_dt[input_set_summary, on = "input_set_id", `:=`(sum_degrees_set = i.sum_degrees_set, max_possible_edges = i.set_size)]

  # lambda calculation (Fast Approximation)
  final_dt[, lambda := (term_degree * sum_degrees_set) / (2 * object@ecprob@graph_size)]

  # P-value and Effect Size
  final_dt[, p_value := calculate_p_value(object@ecprob, observed_edge_count, max_possible_edges, lambda)]
  final_dt[, log2_Anscombe_ratio := 0.5 * (log2(observed_edge_count + 3/8) - log2(lambda + 3/8))]

  # Rename for consistency with package standards
  setnames(final_dt, "input_set_id", "set1")
  setnames(final_dt, "term_id", "set2")

  return(final_dt)
}

terms_ecset_statistics_vectorized_proto <- function(object, input_sets, lambda_method = "fast") {

  # 1. Standardize Input Sets
  if (is.list(input_sets) && !is.data.frame(input_sets)) {
    input_sets_dt <- as.data.table(utils::stack(input_sets))
    setnames(input_sets_dt, c("values", "ind"), c("element", "set_id"))
  } else {
    input_sets_dt <- as.data.table(input_sets)
  }

  # 2. Prepare Bipartite Edges (Once)
  bipartite_edges <- as.data.table(to_dataframe(object))
  setnames(bipartite_edges, "term", "term_id")
  bipartite_edges[, term_id := as.character(term_id)]

  # 3. Call the Core
  final_dt <- calc_ecset_stats_core_proto(object, input_sets_dt, bipartite_edges)

  # 4. Return List (Splitting is the expensive part here)
  return(split(final_dt, by = "set1"))
}


calculate_empirical_fdr_proto <- function(object, real_results_dt, n_permutations = 1000, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # --- STEP 0: Copy Input ---
  # data.table modifies by reference. We must copy to avoid side effects
  # on the user's original object (renaming cols, adding fdr cols).
  real_results_dt <- copy(real_results_dt)

  # --- PRE-STEP: Standardize Column Names ---
  # The core engine uses 'set2' for terms, but single-set methods return 'term_id'.
  # We normalize to 'set2' for this calculation.
  if ("term_id" %in% names(real_results_dt) && !"set2" %in% names(real_results_dt)) {
    setnames(real_results_dt, "term_id", "set2")
  }

  # --- STEP 1: Setup Inputs ---
  if (!"max_possible_edges" %in% names(real_results_dt)) {
    set_size <- max(real_results_dt$observed_edge_count)
  } else {
    set_size <- max(real_results_dt$max_possible_edges)
  }

  bipartite_edges <- as.data.table(to_dataframe(object))
  setnames(bipartite_edges, "term", "term_id")
  bipartite_edges[, term_id := as.character(term_id)]

  # --- STEP 2: Generate "Big Table" of Random Sets ---
  random_elements <- unlist(replicate(n_permutations, sample(object@elements, set_size), simplify = FALSE))

  null_inputs_dt <- data.table(
    set_id = rep(paste0("R", 1:n_permutations), each = set_size),
    element = random_elements
  )

  # --- STEP 3: Run Core Calculation ---
  null_stats_dt <- calc_ecset_stats_core_proto(object, null_inputs_dt, bipartite_edges)

  # --- STEP 4: Statistical Correction ---

  # 4a. Calculate Mean Null Scores per Term
  score_col <- "log2_Anscombe_ratio"

  mean_nulls <- null_stats_dt[, .(
    mean_pos = mean(get(score_col)[get(score_col) >= 0], na.rm = TRUE),
    mean_neg = mean(abs(get(score_col)[get(score_col) < 0]), na.rm = TRUE)
  ), by = .(term_id = set2)]

  mean_nulls[is.nan(mean_pos), mean_pos := 1]
  mean_nulls[is.nan(mean_neg), mean_neg := 1]

  # 4b. Calculate NES for Real Results
  real_results_dt[mean_nulls, on = .(set2 = term_id), `:=`(mean_null_pos = i.mean_pos, mean_null_neg = i.mean_neg)]

  real_results_dt[is.na(mean_null_pos), mean_null_pos := 1]
  real_results_dt[is.na(mean_null_neg), mean_null_neg := 1]

  real_results_dt[, nes := ifelse(get(score_col) >= 0,
                                  get(score_col) / mean_null_pos,
                                  get(score_col) / mean_null_neg)]

  # 4c. Calculate NES for Null Results
  null_stats_dt[mean_nulls, on = .(set2 = term_id), `:=`(mean_null_pos = i.mean_pos, mean_null_neg = i.mean_neg)]
  null_stats_dt[is.na(mean_null_pos), mean_null_pos := 1]
  null_stats_dt[is.na(mean_null_neg), mean_null_neg := 1]

  null_stats_dt[, null_nes := ifelse(get(score_col) >= 0,
                                     get(score_col) / mean_null_pos,
                                     get(score_col) / mean_null_neg)]

  # 4d. Calculate FDR (GSEA Style)
  all_null_nes <- na.omit(null_stats_dt$null_nes)
  null_nes_pos <- all_null_nes[all_null_nes >= 0]
  null_nes_neg <- all_null_nes[all_null_nes < 0]

  nes_real_vec <- real_results_dt$nes

  if(any(is.na(nes_real_vec))) {
    warning("NA values found in NES. Replacing with 0 for FDR calculation.")
    nes_real_vec[is.na(nes_real_vec)] <- 0
  }

  fdr_values <- vapply(nes_real_vec, function(score) {
    if (score >= 0) {
      if (length(null_nes_pos) == 0) return(1)
      frac_null <- (sum(null_nes_pos >= score, na.rm = TRUE) + 1) / (length(null_nes_pos) + 1)
      frac_real <- sum(nes_real_vec >= 0 & nes_real_vec >= score, na.rm = TRUE) / sum(nes_real_vec >= 0, na.rm = TRUE)

      return(if (frac_real == 0) 1 else frac_null / frac_real)
    } else {
      if (length(null_nes_neg) == 0) return(1)
      frac_null <- (sum(null_nes_neg <= score, na.rm = TRUE) + 1) / (length(null_nes_neg) + 1)
      frac_real <- sum(nes_real_vec < 0 & nes_real_vec <= score, na.rm = TRUE) / sum(nes_real_vec < 0, na.rm = TRUE)

      return(if (frac_real == 0) 1 else frac_null / frac_real)
    }
  }, FUN.VALUE = numeric(1))

  fdr_values[fdr_values > 1] <- 1
  real_results_dt[, fdr_q_value := fdr_values]

  # Cleanup
  real_results_dt[, c("mean_null_pos", "mean_null_neg") := NULL]
  setorder(real_results_dt, fdr_q_value, -nes)

  # --- POST-STEP: Restore Column Names ---
  # Return the table with 'term_id' if that's what came in, to match User API
  if ("set2" %in% names(real_results_dt)) {
    setnames(real_results_dt, "set2", "term_id")
  }

  return(real_results_dt)
}
# --- VALIDATION SECTION ---

validate_core <- function(n_trials){

  global_success <- TRUE
  set.seed(1)

  # 1. Setup Data
  data("sample_ects")

  for (trial in 1:n_trials){
    target_term <- sample(sample_ects@terms, 1)
    set_1 <- sample_ects@ecprob@adj[[target_term]] # The element set

    # 2. Run Ground Truth (Existing Package Method)
    # We convert to data.table for fair comparison
    cat("Running Ground Truth (terms_ecset_statistics)...\n")
    truth_dt <- as.data.table(terms_ecset_statistics(sample_ects, set_1, lambda_method = "fast"))
    setorder(truth_dt, term_id)

    # 3. Run Proto Core (New Vectorized Core)
    cat("Running Proto Core (calc_ecset_stats_core_proto)...\n")

    # Prepare inputs for the core
    # A) Input Sets (Single set 'S1')
    input_sets_dt <- data.table(
      set_id = "S1",
      element = set_1
    )

    # B) Bipartite Edges (External to function)
    bipartite_edges <- as.data.table(to_dataframe(sample_ects))
    setnames(bipartite_edges, "term", "term_id")
    bipartite_edges[, term_id := as.character(term_id)]

    # C) Execute
    proto_dt <- calc_ecset_stats_core_proto(sample_ects, input_sets_dt, bipartite_edges)

    # Clean up Proto result to match Truth structure for comparison
    # (Remove extra columns returned by core like degrees, sum_degrees, input_set_id)
    proto_clean <- proto_dt[, .(term_id = set2,
                                p_value,
                                lambda,
                                observed_edge_count,
                                log2_Anscombe_ratio)]
    setorder(proto_clean, term_id)

    # 4. Comparisons
    cat("\n--- VALIDATION REPORT ---\n")

    # A. Dimension Check
    if (nrow(truth_dt) == nrow(proto_clean)) {
      cat("[PASS] Row counts match:", nrow(truth_dt), "\n")
    } else {
      stop("[FAIL] Row counts mismatch! Truth:", nrow(truth_dt), " Proto:", nrow(proto_clean))
    }

    # B. Column Equality Check (Tolerance for floating point)
    cols_to_check <- c("lambda", "p_value", "log2_Anscombe_ratio", "observed_edge_count")
    all_pass <- TRUE

    for (col in cols_to_check) {
      # Handle NA comparisons safely
      truth_vec <- truth_dt[[col]]
      proto_vec <- proto_clean[[col]]

      # Check numeric difference
      # Use a small tolerance (e.g., 1e-14)
      max_diff <- max(abs(truth_vec - proto_vec), na.rm = TRUE)

      if (max_diff < 1e-14) {
        cat(sprintf("[PASS] Column '%s' matches (Max Diff: %g)\n", col, max_diff))
      } else {
        cat(sprintf("[FAIL] Column '%s' mismatch! (Max Diff: %g)\n", col, max_diff))
        all_pass <- FALSE
      }
    }

    if (!all_pass) {
      global_success <- FALSE
    }
  }
  if (global_success){
    message("--- Global Success ---")
  } else {
    message("--- FAILURE ---")
  }
}

validate_vectorized_wrapper <- function() {
  message("\n--- VALIDATING WRAPPER C (Vectorized List Return) ---")

  data("sample_ects")
  set.seed(123)

  # 1. Create a multi-set input (3 distinct sets)
  # We pick 3 random terms to generate 3 sets of elements
  target_terms <- sample(sample_ects@terms, 3)

  # Construct the input data.table manually to mimic user input
  input_list <- lapply(target_terms, function(t) sample_ects@ecprob@adj[[t]])
  names(input_list) <- paste0("Set_", 1:3)

  # Convert to DT for the Core comparison
  input_dt <- as.data.table(stack(input_list))
  setnames(input_dt, c("values", "ind"), c("element", "set_id"))

  # 2. Run the Wrapper (Proto C)
  # Pass the LIST to test the internal conversion logic
  cat("Running Vectorized Wrapper (Input = List)...\n")
  wrapper_result <- terms_ecset_statistics_vectorized_proto(sample_ects, input_list)

  # 3. Run the Core (Proto D) for comparison
  # Pass the DT and external edges
  cat("Running Core Engine for comparison...\n")
  bipartite_edges <- as.data.table(to_dataframe(sample_ects))
  setnames(bipartite_edges, "term", "term_id")
  bipartite_edges[, term_id := as.character(term_id)]

  core_result <- calc_ecset_stats_core_proto(sample_ects, input_dt, bipartite_edges)

  # 4. CHECKS

  # Check 1: Is the output a list?
  if (is.list(wrapper_result) && !is.data.table(wrapper_result)) {
    cat("[PASS] Output is a list.\n")
  } else {
    stop("[FAIL] Output is not a list.")
  }

  # Check 2: Length matches input sets
  if (length(wrapper_result) == 3) {
    cat("[PASS] Output list length is 3.\n")
  } else {
    stop("[FAIL] Output list length mismatch.")
  }

  # Check 3: Data Integrity (Re-bind list and compare to Core result)
  wrapper_flat <- rbindlist(wrapper_result)

  # Sort both for comparison
  setorder(wrapper_flat, set1, set2)
  setorder(core_result, set1, set2)

  # Check row counts
  if (nrow(wrapper_flat) == nrow(core_result)) {
    cat("[PASS] Re-bound list matches Core row count.\n")
  } else {
    stop("[FAIL] Data loss detected in wrapper splitting.")
  }

  # Check a value sum (hash-like check)
  diff_sum <- sum(abs(wrapper_flat$log2_Anscombe_ratio - core_result$log2_Anscombe_ratio))
  if (diff_sum < 1e-10) {
    cat("[PASS] Numerical values match Core exactly.\n")
  } else {
    stop("[FAIL] Numerical mismatch between Wrapper and Core.")
  }

  message("--- Wrapper C Validation: SUCCESS ---")
}

validate_fdr_wrapper <- function() {
  message("\n--- VALIDATING WRAPPER B (Empirical FDR) ---")

  data("sample_ects")
  set.seed(1)

  # 1. Setup a "Real" Result with a strong signal
  # Use a known term so we expect a low FDR
  target_term <- "GO:0005787"
  set_1 <- sample_ects@ecprob@adj[[target_term]]

  cat("Generating 'Real' results for a known strong term...\n")
  # Run single set analysis (Ground Truth style)
  real_results <- as.data.table(terms_ecset_statistics(sample_ects, set_1, lambda_method = "fast"))

  # Add max_possible_edges column (Required input for FDR function)
  # In this context, max edges = set size
  real_results[, max_possible_edges := length(set_1)]

  # 2. Run FDR Proto
  # We use n=100 for speed in dev testing
  cat("Running FDR Calculation (n=100 permutations)...\n")
  fdr_results <- calculate_empirical_fdr_proto(sample_ects, real_results, n_permutations = 100, seed = 1)

  # 3. CHECKS

  # Check 1: Columns exist
  if (all(c("nes", "fdr_q_value") %in% names(fdr_results))) {
    cat("[PASS] NES and fdr_q_value columns created.\n")
  } else {
    stop("[FAIL] Missing output columns.")
  }

  # Check 2: Value Bounds
  if (min(fdr_results$fdr_q_value) >= 0 && max(fdr_results$fdr_q_value) <= 1) {
    cat("[PASS] FDR q-values are within [0, 1].\n")
  } else {
    stop("[FAIL] FDR q-values out of bounds.")
  }

  # Check 3: Logical Consistency (Signal Detection)
  # The target term "GO:0005787" should be at the top with low FDR
  top_hit <- fdr_results[term_id == target_term]

  cat(sprintf("Target Term Stats -> NES: %.2f, FDR: %.4f\n", top_hit$nes, top_hit$fdr_q_value))

  if (top_hit$fdr_q_value < 0.05) {
    cat("[PASS] Target term correctly identified as significant (FDR < 0.05).\n")
  } else {
    warning("[WARN] Target term FDR is high. This might happen with low n_permutations, but check logic.")
  }

  # Check 4: Sorting
  # The function should return sorted by FDR then NES
  is_sorted <- !is.unsorted(fdr_results$fdr_q_value)
  if (is_sorted) {
    cat("[PASS] Result table is sorted by FDR.\n")
  } else {
    stop("[FAIL] Result table is not sorted.")
  }

  message("--- Wrapper B Validation: SUCCESS ---")
}

###### Main testing ##########

# validate_core(100)
validate_vectorized_wrapper()
validate_fdr_wrapper()
