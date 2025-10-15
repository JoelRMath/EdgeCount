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



test_that("bipartite degrees", {
  simul <- FALSE
  if (simul){
    ect <- sample_ects
    degrees <- unlist(ect@ecprob@degrees[ect@elements])
    nsim <- 10000
    m <- sapply(1:nsim, function(x){
        sample(degrees, size=length(degrees), replace = FALSE)
      }
    )
    print(dim(m))
    saveRDS(m,test_path("res/DegreeMat.rds"))
  }
  toggle <- FALSE
  if (toggle){
    m <- readRDS(test_path("res/DegreeMat.rds"))
    n <- length(m[,1])
    sk <- sd(m[,1])
    i <- c(1:n)
    N <- rep(n, n)
    N_1 <- rep(n-1, n)
    sk2 <- rep(sk*sk, n)
    sdp <- sqrt(i*sk2*(N-i)/(N_1))
    m <- apply(m, 2, cumsum)
    ind <- c(1:length(m[,1]))
    mn <- apply(m, 1, mean)
    s <- apply(m, 1, sd)
    df <- data.frame(index = ind, mean = mn, sd = s)
    plot(df$index, df$mean)
    plot(df$index, df$sd)
    print(cor.test(df$index, df$mean))
    print(cor.test(df$index, df$sd))
    print(length(i))
    print(length(N))
    print(length(N_1))
    print(length(N-i))
    print(length(sk2))
    print(length(sdp))
    plot(df$sd, sdp)
    print(cor.test(df$sd, sdp))
  }
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
