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
  ect <- ECTermScoring(test_path("network/terms_edges.txt"))
  element_set <- sample(ect@elements, 10)
  t <- terms_ecset_statistics(ect, element_set = element_set, lambda_method = "fast")
})

test_that("ECTermScoring rank scoring", {
  ect <- ECTermScoring(test_path("network/terms_edges.txt"))
  element_to_ranks <- as.list(1:length(ect@elements))
  names(element_to_ranks) <- sample(ect@elements, length(ect@elements))
  n_terms <- 10
  terms <- sample(ect@terms, n_terms)
  time_start <- Sys.time()
  terms_ecranks_statistics(ect, element_to_ranks, terms)
  time_diff <- Sys.time() - time_start
  # print(length(ect@terms))
  # print(as.numeric(time_diff, units = "secs"))
  
})

test_that("ECTermScoring table rank scoring", {
  ect <- ECTermScoring(test_path("network/terms_edges.txt"))
  element_to_ranks <- as.list(1:length(ect@elements))
  names(element_to_ranks) <- sample(ect@elements, length(ect@elements))
  n_terms <- length(ect@terms)
  n_terms <- 100
  terms <- sample(ect@terms, n_terms)
  lst <- terms_ecranks_statistics(ect, element_to_ranks, terms)
  allowed = c("lambda", "z_score", "log2_anscombe_ratio", 
              "log2_relative_change", "p_value")
  
  expect_error(
    table_terms_ecranks_statistics(2,
                                   "z_score",
                                   "max"))

  scoring_statistic <- "zscore"
  expect_warning(
    table_terms_ecranks_statistics(lst,
                                   scoring_statistic,
                                   "max"),
    paste0("Invalid 'scoring_statistic': ", scoring_statistic, "Allowed choices are: ", toString(allowed)))
  
  scoring_statistic <- "z_score"
  res <- table_terms_ecranks_statistics(lst,
                                        scoring_statistic,
                                        "max")
  write.table(res, test_path("res/z_score_max.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  res <- table_terms_ecranks_statistics(lst,
                                        scoring_statistic,
                                        "median")
  write.table(res, test_path("Res/z_score_median.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  res <- table_terms_ecranks_statistics(lst,
                                        scoring_statistic,
                                        "min")
  write.table(res, test_path("Res/z_score_min.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
})
