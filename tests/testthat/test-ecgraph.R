library(testthat)
library(EdgeCount)

test_that("ECGraph constructor", {

  toggle <- TRUE
  if (toggle){
    ec_graph <- ECGraph(test_path("network/toy.txt"))
    # adj list
    expect_setequal(ec_graph@adj$A, c("E"))
    expect_setequal(ec_graph@adj$B, c("E","C"))
    expect_setequal(ec_graph@adj$C, c("B","E","F"))
    expect_setequal(ec_graph@adj$D, c("G"))
    expect_setequal(ec_graph@adj$E, c("A", "B", "C", "F"))
    expect_setequal(ec_graph@adj$F, c("E", "C"))
    expect_setequal(ec_graph@adj$G, c("D", "H"))
    expect_setequal(ec_graph@adj$H, c("G"))
    # degrees
    expect_equal(ec_graph@degrees[["A"]], 1)
    expect_equal(ec_graph@degrees[["B"]], 2)
    expect_equal(ec_graph@degrees[["C"]], 3)
    expect_equal(ec_graph@degrees[["D"]], 1)
    expect_equal(ec_graph@degrees[["E"]], 4)
    expect_equal(ec_graph@degrees[["F"]], 2)
    expect_equal(ec_graph@degrees[["G"]], 2)
    expect_equal(ec_graph@degrees[["H"]], 1)
    # names
    expect_setequal(ec_graph@names, c("A", "B", "C", "D", "E", "F", "G", "H"))
  }
})

test_that("EC in", {
  
  toggle <- TRUE
  if (toggle){
    ec <- ECGraph(test_path("network/toy.txt"))
    expect_equal(get_edge_count_in(ec,NULL), 0)
    expect_equal(get_edge_count_in(ec,c("E")), 0)
    expect_equal(get_edge_count_in(ec,c("E","I")), 0)
    expect_equal(get_edge_count_in(ec,c("E", "D")), 0)
    expect_equal(get_edge_count_in(ec,c("E", "D", "G")), 1)
    expect_equal(get_edge_count_in(ec,c("E", "D", "G", "C", "B")), 4)
  }
})

test_that("EC between", {
  
  toggle <- TRUE
  if (toggle){
    ec <- ECGraph(test_path("network/toy.txt"))
    expect_equal(get_edge_count_between(ec,NULL,NULL), 0)
    expect_equal(get_edge_count_between(ec,c("A"),NULL), 0)
    expect_equal(get_edge_count_between(ec,c("A"),c("G")), 0)
    expect_equal(get_edge_count_between(ec,c("E","A"),c("D","G")), 0)
    expect_equal(get_edge_count_between(ec,c("E","A"),c("B","C","F")), 3)
    expect_equal(get_edge_count_between(ec,c("H","D","A","B","C","F"),c("E","G")), 6)
  }
})

test_that("EC in max fds", {
  
  toggle <- TRUE
  if (toggle){
    ec <- ECGraph(test_path("network/toy.txt"))
    expect_equal(get_edge_count_in_max_fds(ec,NULL), 0)
    expect_equal(get_edge_count_in_max_fds(ec,c("E")), 0)
    expect_equal(get_edge_count_in_max_fds(ec,c("B","C","E")), 3)
    expect_equal(get_edge_count_in_max_fds(ec,c("B","C","E","G")), 5)
  }
})

test_that("EC between max fds", {
  toggle <- TRUE
  if (toggle){
    ec <- ECGraph(test_path("network/toy.txt")) 
    expect_equal(get_edge_count_between_max_fds(ec,NULL,NULL), 0)
    expect_equal(get_edge_count_between_max_fds(ec,c("E"),NULL), 0)
    expect_equal(get_edge_count_between_max_fds(ec,c("E"),c("A")), 1)
    expect_equal(get_edge_count_between_max_fds(ec,c("E"),c("D","H")), 2)
    # Corrected expectation for the last test case
    expect_equal(get_edge_count_between_max_fds(ec,c("E","F"),c("D","G","H","B","C")), 6) 
  }
})

test_that("EC to_dataframe", {
  
  toggle <- TRUE
  if (toggle){
    ec1 <- ECGraph(test_path("network/toy.txt"))
    df1 <- to_dataframe(ec1)
    ec2 <- ECGraph(df1)
    expect_equal(ec1@adj, ec2@adj)
    expect_equal(ec1@degrees, ec2@degrees)
    expect_equal(ec1@names, ec2@names)
  }
})

