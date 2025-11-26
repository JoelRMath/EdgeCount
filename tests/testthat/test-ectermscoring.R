library(testthat)
library(EdgeCount)

create_temp_edge_file <- function(data, file_path = tempfile(fileext = ".txt")) {
  write.table(data, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  return(file_path)
}

test_that("ECTermScoring constructor works with file input and default columns", {
  term_element_data <- data.frame(
    TermCol = c("T1", "T1", "T2", "T2", "T3", "T5", "T5"),
    ElementCol = c("E1", "E2", "E2", "E3", "E4", "E1", "E5"),
    stringsAsFactors = FALSE
  )
  temp_file <- create_temp_edge_file(term_element_data)

  ets_obj <- ECTermScoring(temp_file)

  expected_terms <- sort(c("T1", "T2", "T3", "T5"))
  expected_elements <- sort(c("E1", "E2", "E3", "E4", "E5"))
  expected_all_vertices_in_graph <- sort(c(expected_terms, expected_elements))

  expect_s4_class(ets_obj, "ECTermScoring")
  expect_s4_class(ets_obj@ecprob, "ECProb")

  expect_equal(sort(ets_obj@terms), expected_terms)
  expect_equal(sort(ets_obj@elements), expected_elements)

  expect_equal(sort(ets_obj@ecprob@names), expected_all_vertices_in_graph)
  expect_equal(ets_obj@ecprob@graph_size, 7)

  expect_equal(ets_obj@ecprob@degrees[["T1"]], 2)
  expect_equal(ets_obj@ecprob@degrees[["E2"]], 2)
  expect_equal(ets_obj@ecprob@degrees[["T3"]], 1)
  expect_equal(ets_obj@ecprob@degrees[["E5"]], 1)

  unlink(temp_file)
})

test_that("ECTermScoring constructor works with data frame input and named columns", {
  term_element_df <- data.frame(
    MyCategory = c("TC1", "TC1", "TC2", "TC2"),
    MyItem = c("ItemX", "ItemY", "ItemY", "ItemZ"),
    stringsAsFactors = FALSE
  )

  ets_obj <- ECTermScoring(term_element_df, col_term = "MyCategory", col_element = "MyItem")

  expected_terms_df <- sort(c("TC1", "TC2"))
  expected_elements_df <- sort(c("ItemX", "ItemY", "ItemZ"))

  expect_s4_class(ets_obj, "ECTermScoring")
  expect_equal(sort(ets_obj@terms), expected_terms_df)
  expect_equal(sort(ets_obj@elements), expected_elements_df)
  expect_equal(ets_obj@ecprob@graph_size, 4) # 4 unique links
  expect_equal(ets_obj@ecprob@degrees[["TC1"]], 2)
  expect_equal(ets_obj@ecprob@degrees[["ItemY"]], 2)
})

test_that("ECTermScoring constructor handles empty input gracefully", {
  empty_df <- data.frame(Term = character(0), Element = character(0))

  ets_obj_empty <- ECTermScoring(empty_df)

  expect_s4_class(ets_obj_empty, "ECTermScoring")
  expect_length(ets_obj_empty@terms, 0)
  expect_length(ets_obj_empty@elements, 0)
  expect_s4_class(ets_obj_empty@ecprob, "ECProb")
  expect_length(ets_obj_empty@ecprob@names, 0)
  expect_equal(ets_obj_empty@ecprob@graph_size, 0)
  expect_equal(length(ets_obj_empty@ecprob@adj), 0)
  expect_equal(length(ets_obj_empty@ecprob@degrees), 0)
})

test_that("ECTermScoring constructor stops if vertices are identified as both terms and elements", {
  overlapping_df <- data.frame(
    TermColumn = c("vertexA", "vertexB", "vertexC"), # vertexC is also an element
    ElementColumn = c("vertexX", "vertexC", "vertexY"), # vertexC is also a term
    stringsAsFactors = FALSE
  )

  expect_error(
    ECTermScoring(overlapping_df, col_term = "TermColumn", col_element = "ElementColumn"),
    "Some vertex IDs are identified as both elements and terms within the constructed graph: vertexC"
  )
})

test_that("ECTermScoring constructor handles input with terms/elements not forming edges", {

  edge_data_with_isolates_in_input <- data.frame(
    T = c("T1", "T1", "T2"), # T_Isolated is a potential term, E_Isolated a potential element
    E = c("E1", "E2", "E1"),
    stringsAsFactors = FALSE
  )

  ets_obj <- ECTermScoring(edge_data_with_isolates_in_input, col_term = "T", col_element = "E")

  expect_equal(sort(ets_obj@terms), sort(c("T1", "T2")))
  expect_equal(sort(ets_obj@elements), sort(c("E1", "E2")))
  expect_equal(ets_obj@ecprob@graph_size, 3) # T1-E1, T1-E2, T2-E1
})

test_that("ECTermScoring constructor handles ... for file reading options", {
  semicolon_data <- data.frame(
    c("TermX", "TermX"),
    c("ElemA", "ElemB")
  )
  temp_file_sc <- tempfile(fileext = ".txt")
  write.table(semicolon_data, temp_file_sc, sep = ";", row.names = FALSE, quote = FALSE, col.names = FALSE)

  expect_error(ECTermScoring(temp_file_sc))

  ets_obj_sc <- ECTermScoring(temp_file_sc, col_term = 1, col_element = 2, header = FALSE, sep = ";")

  expect_s4_class(ets_obj_sc, "ECTermScoring")
  expect_equal(sort(ets_obj_sc@terms), "TermX")
  expect_equal(sort(ets_obj_sc@elements), sort(c("ElemA", "ElemB")))
  expect_equal(ets_obj_sc@ecprob@graph_size, 2)

  unlink(temp_file_sc)
})

test_that("ECTermScoring one set", {
  ect <- sample_ects
  element_set <- sample(ect@elements, 10)
  t <- terms_ecset_statistics(ect, element_set = element_set, lambda_method = "fast")
})

create_vsea_test_object <- function(full_ects, n_terms) {

  if (n_terms >= length(full_ects@terms)) {
    selected_terms <- full_ects@terms
  } else {
    selected_terms <- sample(full_ects@terms, n_terms)
  }

  adj_list <- full_ects@ecprob@adj

  term_neighbors_list <- adj_list[selected_terms]

  subset_edges_df <- utils::stack(term_neighbors_list)
  names(subset_edges_df) <- c("element", "term") # Rename for clarity

  new_ects <- ECTermScoring(subset_edges_df, col_term = "term", col_element = "element")

  return(new_ects)
}

test_that("terms_ecset_statistics_vectorized is correct", {

  # --- Gold Standard: A slow but safe wrapper around the original function ---
  terms_ecset_statistics_slow_wrapper <- function(object, input_sets_dt, lambda_method = "fast") {
    all_set_ids <- unique(input_sets_dt$set_id)
    results_list <- list()

    for (current_set_id in all_set_ids) {
      current_element_set <- input_sets_dt[set_id == current_set_id, element]

      single_result_df <- terms_ecset_statistics(
        object,
        element_set = current_element_set,
        lambda_method = lambda_method
      )

      if (!is.null(single_result_df) && nrow(single_result_df) > 0) {
        single_result_dt <- as.data.table(single_result_df)
        # Add the input_set_id to the results table
        single_result_dt[, input_set_id := current_set_id]
        results_list[[current_set_id]] <- single_result_dt
      }
    }
    return(results_list)
  }

  # --- 1. Define the test data ---
  te_df <- data.frame(
    term = c("T1", "T1", "T2", "T3"),
    element = c("E1", "E2", "E2", "E3")
  )
  ects <- ECTermScoring(te_df)

  input_sets_dt <- data.table(
    set_id = c("SetA", "SetA", "SetB", "SetB"),
    element = c("E1", "E3", "E2", "E99") # E99 is not in the graph
  )

  # --- 2. Run both the slow and fast versions ---
  results_slow <- terms_ecset_statistics_slow_wrapper(ects, input_sets_dt)

  # This now calls the formal S4 method from the package
  results_fast <- terms_ecset_statistics_vectorized(ects, input_sets_dt)

  # --- 3. Perform the automated comparison ---
  expect_setequal(names(results_slow), names(results_fast))

  harmonize_df <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(data.table())
    # The original function has slightly different column names
    if ("log2_Anscombe" %in% names(df)) setnames(df, "log2_Anscombe", "log2_Anscombe_ratio")
    if ("observed_edge_count" %in% names(df)) setnames(df, "observed_edge_count", "observed_edges")

    # Rename for consistency
    if("input_set_id" %in% names(df)) setnames(df, "input_set_id", "set1")
    if("term_id" %in% names(df)) setnames(df, "term_id", "set2")

    cols_to_keep <- c("set1", "set2", "observed_edges", "lambda", "p_value", "log2_Anscombe_ratio")

    # Ensure all columns to keep actually exist before subsetting
    cols_that_exist <- intersect(cols_to_keep, names(df))

    return(df[, ..cols_that_exist])
  }

  for(set_id in names(results_slow)) {
    slow_dt <- harmonize_df(results_slow[[set_id]])
    fast_dt <- harmonize_df(results_fast[[set_id]])

    setorder(slow_dt, set2)
    setorder(fast_dt, set2)

    expect_equal(slow_dt, fast_dt, info = paste("Results for set_id", set_id, "do not match."))
  }
})

test_that("ECTermScoring vsea", {

  toggle <- FALSE
  if (toggle){
    ect <- sample_ects
    n_terms <- 100
    ects <- create_vsea_test_object(ect, n_terms)
    element_to_ranks <- as.list(1:length(ects@elements))
    names(element_to_ranks) <- sample(ects@elements, length(ects@elements))
    t <- run_vsea_analysis(ects, element_to_ranks,
                           scoring_statistic = "log2_Anscombe_ratio",
                           n_permutations = 1000,
                           seed = NULL)
    write.table(t$max_score_summary, test_path("res/vsea_max_summary.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(t$min_score_summary, test_path("res/vsea_min_summary.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(t$median_score_summary, test_path("res/vsea_median_summary.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }

})

test_that("ECTermScoring terms_ecranks_statistics S4 method vs gold standard", {

  terms_ecranks_statistics_gs_loop <- function(object, element_ranks) {

    # 1
    valid_elements <- intersect(names(element_ranks), object@elements)
    if (length(valid_elements) < 1) {
      stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
    }
    object <- reduce_universe_by_elements(object, valid_elements)
    element_ranks <- rank(element_ranks[valid_elements]) # Re-rank from 1 to N_valid

    # 2
    ecprob <- object@ecprob
    M_g <- ecprob@graph_size
    N_e <- length(object@elements)
    all_terms <- object@terms

    # 3
    df <- data.frame(element = names(element_ranks),
                     rank = as.numeric(element_ranks),
                     stringsAsFactors = FALSE)
    df <- df[order(df$rank),]
    ranked_elements_sorted <- df$element

    K <- unlist(ecprob@degrees[ranked_elements_sorted])
    cumul_sum_K <- cumsum(K)

    score_one_term_safe <- function(term_id) {

      elements_in_term <- get_neighbors(ecprob, term_id)

      elements_in_term <- intersect(elements_in_term, names(element_ranks))

      if (length(elements_in_term) == 0) return(NULL)

      ranks_for_term <- element_ranks[elements_in_term]

      ranks_for_term_sorted <- sort(ranks_for_term)
      elements_term_sorted <- names(ranks_for_term_sorted)

      sz <- length(elements_term_sorted)

      observed_ec_term <- 1:sz
      term_size <- rep(sz, sz)
      max_ec_term <- pmin(term_size, ranks_for_term_sorted)
      term_degree <- rep(ecprob@degrees[[term_id]], sz)

      lambda_term <- (term_degree / (2 * M_g)) * cumul_sum_K[ranks_for_term_sorted]

      p_value_term <- calculate_p_value(ecprob, observed_ec_term, max_ec_term, lambda_term)
      log2_anscombe_term <- 0.5 * (log2(observed_ec_term + 3/8) - log2(lambda_term + 3/8))
      log2_rel_change_term <- log2(observed_ec_term) - log2(lambda_term)

      df <- data.frame(
        element_id = elements_term_sorted,
        global_rank = element_ranks[elements_term_sorted],
        rank_in_term = rank(ranks_for_term_sorted),
        observed_ec = observed_ec_term,
        max_ec = max_ec_term,
        term_size = term_size,
        lambda = lambda_term,
        p_value = p_value_term,
        log2_Anscombe_ratio = log2_anscombe_term,
        stringsAsFactors = FALSE
      )
      row.names(df) <- NULL
      return(as.data.table(df))
    }

    all_results_list <- lapply(all_terms, score_one_term_safe)
    names(all_results_list) <- all_terms

    all_results_list <- all_results_list[!sapply(all_results_list, is.null)]

    return(all_results_list)
  }

  te_df <- data.frame(
    term = c("TA", "TA", "TB", "TB", "TB", "TC", "TC", "TC", "TC", "TD", "TD", "TD"),
    element = c("E1", "E2", "E2", "E3", "E4", "E2", "E3", "E5", "E6", "E7", "E8", "E5")
  )
  ects <- ECTermScoring(te_df)

  element_ranks <- c(E1 = 2, E2 = 1, E3 = 3, E5 = 4, E6 = 7, E7 = 6, E8 = 8)

  results_gs <- terms_ecranks_statistics_gs_loop(ects, element_ranks)
  results_s4 <- terms_ecranks_statistics(ects, element_ranks)

  expect_setequal(names(results_gs), names(results_s4))

  sorted_gs <- results_gs[sort(names(results_gs))]
  sorted_s4 <- results_s4[sort(names(results_s4))]

  expect_equal(sorted_gs, sorted_s4)

})

test_that("run_vsea_analysis correctly finds enrichment", {

  te_df <- data.frame(
    term = c("TermA", "TermA", "TermB", "TermB", "TermC", "TermC"),
    element = c("E1", "E2", "E3", "E4", "E5", "E6")
  )
  ects <- ECTermScoring(te_df)

  element_ranks <- c(E1 = 1, E2 = 2, E3 = 3, E4 = 4, E5 = 5, E6 = 6)

 vsea_results <- suppressWarnings(run_vsea_analysis(
    ects,
    element_ranks,
    scoring_statistic = "log2_Anscombe_ratio", # Use the default, robust statistic
    n_permutations = 100,
    seed = 1
  ))

  expect_true(is.list(vsea_results))
  expect_equal(sort(names(vsea_results)),
               sort(c("max_score_summary", "min_score_summary", "median_score_summary")))

  max_summary <- vsea_results$max_score_summary
  expect_equal(max_summary$term_id[1], "TermA")
  expect_gt(max_summary$nes[1], 0)

  min_summary <- vsea_results$min_score_summary
  median_summary <- vsea_results$median_score_summary

  nes_TermA <- min_summary[term_id == "TermA", nes]
  nes_TermB <- min_summary[term_id == "TermB", nes]

  expect_lt(nes_TermB, nes_TermA)
})

test_that("terms_ecset_statistics_fdr works as a wrapper", {

  # 1. Setup a small graph and known set
  data("sample_ects")
  target_term <- "GO:0005787"
  set_1 <- sample_ects@ecprob@adj[[target_term]]

  # 2. Run the wrapper
  # We use n=50 for speed in testing
  results <- terms_ecset_statistics_fdr(sample_ects, set_1, n_permutations = 50, seed = 123)

  # 3. Check Return Type
  expect_true(is.data.table(results))

  # 4. Check Columns
  expected_cols <- c("term_id", "observed_edge_count", "log2_Anscombe_ratio",
                     "p_value", "nes", "fdr_q_value")
  expect_true(all(expected_cols %in% names(results)))

  # 5. Check Signal Logic
  # The target term itself should be highly significant
  top_hit <- results[term_id == target_term]
  expect_equal(nrow(top_hit), 1)
  expect_lt(top_hit$fdr_q_value, 0.05) # Should be significant

  # 6. Check robustness with empty input
  # An empty set should return NULL (handled by the wrapper)
  empty_res <- suppressWarnings(terms_ecset_statistics_fdr(sample_ects, character(0), n_permutations = 10))
  expect_null(empty_res)

  # 7. Check robustness with disjoint input (no connections)
  # Create a dummy element that doesn't exist in the graph
  disjoint_res <- suppressWarnings(terms_ecset_statistics_fdr(sample_ects, c("GHOST_ELEMENT"), n_permutations = 10))
  expect_null(disjoint_res)
})
