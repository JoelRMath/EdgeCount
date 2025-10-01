library(EdgeCount)
library(data.table)

get_candidate_pairs_v1 <- function(){

  candidate_pairs <- list()

  network_edges <- data.table(to_dataframe(sample_ecg))
  setnames(network_edges, c("from", "to"), c("element1", "element2"))

  bipartite_edges <- as.data.table(to_dataframe(sample_ects))
  bipartite_edges[, `:=`(term = as.character(term), element = as.character(element))]

  all_terms <- unique(bipartite_edges[,term])

  for (term_a in all_terms){
    elements_a <- bipartite_edges[term == term_a, element]
    neighbors <- get_neighbors(sample_ecg, elements_a)
    connected_terms <- unlist(sapply(neighbors, function (vertex){
      term <- bipartite_edges[element == vertex, term]
      term
    }))
    connected_terms <- unique(connected_terms)
    valid_terms_b <- NULL
    if (length(connected_terms > 0)){
      for (term_b in connected_terms){
        if (term_b != term_a){
          valid_terms_b <- c(valid_terms_b, term_b)
        }
      }
      if (!is.null(valid_terms_b)){
        valid_terms_b <- unique(valid_terms_b)
        candidate_pairs[[term_a]] <- valid_terms_b
      }
    }
  }
  candidate_dt <- utils::stack(candidate_pairs)
  candidate_dt <- as.data.table(candidate_dt)
  setnames(candidate_dt, c("values", "ind"), c("term2", "term1"))
  candidate_dt <- candidate_dt[, .(term1, term2)]

  candidate_dt[, `:=`(term1 = as.character(term1), term2 = as.character(term2))]

  canonical_pairs <- candidate_dt[, .(
    term1 = pmin(term1, term2),
    term2 = pmax(term1, term2)
  )]
  unique_pairs_dt <- unique(canonical_pairs)

  return(unique_pairs_dt)
}

get_candidate_pairs_v2 <- function(){

  # We'll collect small data.tables in a list
  pair_list <- list()

  network_edges <- data.table(to_dataframe(sample_ecg))
  setnames(network_edges, c("from", "to"), c("element1", "element2"))

  bipartite_edges <- as.data.table(to_dataframe(sample_ects))
  bipartite_edges[, `:=`(term = as.character(term), element = as.character(element))]

  all_terms <- unique(bipartite_edges[,term])

  for (term_a in all_terms){
    elements_a <- bipartite_edges[term == term_a, element]
    neighbors <- get_neighbors(sample_ecg, elements_a)

    connected_terms_list <- sapply(neighbors, function (vertex){
      bipartite_edges[element == vertex, term]
    })
    connected_terms <- unique(unlist(connected_terms_list, use.names = FALSE))

    valid_terms_b <- NULL
    if (length(connected_terms) > 0){
      # Filter out self-loops
      valid_terms_b <- connected_terms[connected_terms != term_a]

      if (length(valid_terms_b) > 0){
        # Create a small data.table for this iteration's pairs
        # and add it to our list
        pairs_found <- data.table(term1 = term_a, term2 = valid_terms_b)
        pair_list <- c(pair_list, list(pairs_found))
      }
    }
  }

  # --- STEP 1: Combine all found pairs into one data.table ---
  # This will contain duplicates like (A, B) and (B, A)
  all_pairs_dt <- rbindlist(pair_list)
  all_pairs_dt[, `:=`(term1 = as.character(term1), term2 = as.character(term2))]

  # --- STEP 2: Create a canonical representation to find unique edges ---
  # pmin/pmax sorts each row alphabetically, ensuring (B, A) becomes (A, B)
  canonical_pairs <- all_pairs_dt[, .(
    term1 = pmin(term1, term2),
    term2 = pmax(term1, term2)
  )]

  # --- STEP 3: Get the final, unique set of undirected edges ---
  unique_pairs_dt <- unique(canonical_pairs)

  return(unique_pairs_dt)
}

get_candidate_pairs_v3 <- function() {

  # 1. Get the two fundamental edge lists as data.tables
  network_edges <- data.table(to_dataframe(sample_ecg))
  setnames(network_edges, c("from", "to"), c("element1", "element2"))

  bipartite_edges <- as.data.table(to_dataframe(sample_ects))
  bipartite_edges[, `:=`(term = as.character(term), element = as.character(element))]

  # 2. Annotate the network edges with their corresponding terms using joins
  #    This is the core of the fast algorithm.
  merged1 <- network_edges[bipartite_edges, on = .(element1 = element), nomatch = 0, allow.cartesian = TRUE]
  setnames(merged1, "term", "term1")

  merged2 <- merged1[bipartite_edges, on = .(element2 = element), nomatch = 0, allow.cartesian = TRUE]
  setnames(merged2, "term", "term2")

  # 3. Get the unique, canonical pairs from the merged result
  unique_pairs_dt <- merged2[term1 != term2, # Exclude self-connections
                             .(term1 = pmin(term1, term2), term2 = pmax(term1, term2)) # Create canonical pairs
  ][, unique(.SD)] # Get the unique rows

  return(unique_pairs_dt)
}

score_pairs_v1 <- function(candidate_pairs_dt) {

  # --- Setup: Create objects needed for calculations once ---
  # This is more efficient than creating them inside the loop.
  ecp <- ECProb(sample_ecg)

  bipartite_edges <- as.data.table(to_dataframe(sample_ects))
  bipartite_edges[, term := as.character(term)]

  # Create a fast lookup list (hash map) for term -> elements
  term_to_elements_list <- split(bipartite_edges$element, bipartite_edges$term)

  # Pre-allocate a list to store the results. This is more memory-efficient
  # than repeatedly growing a list with c().
  results_list <- vector("list", nrow(candidate_pairs_dt))


  # --- The "Slow and Safe" Loop ---
  # We iterate through each row of the candidate pairs.
  for (i in 1:nrow(candidate_pairs_dt)) {
    # Print progress
    # if (i %% 10000 == 0) {
    #   message(paste("Scoring pair", i, "of", nrow(candidate_pairs_dt)))
    # }

    # Get the terms for the current row
    term1 <- candidate_pairs_dt$term1[i]
    term2 <- candidate_pairs_dt$term2[i]

    # Get the element sets using the fast lookup list
    elements1 <- term_to_elements_list[[term1]]
    elements2 <- term_to_elements_list[[term2]]

    # Calculate the observed edge count
    observed <- get_edge_count_between(sample_ecg, elements1, elements2)

    # Get the full statistics
    stats <- edge_count_statistics(ecp, elements1, elements2, observed, lambda_method = "fast")

    # Store the results for this row in our pre-allocated list
    results_list[[i]] <- list(
      observed_edges = as.integer(observed),
      p_value = stats$p_value,
      lambda = stats$lambda,
      log2_anscombe = stats$log2_Anscombe_ratio
    )
  }

  # --- Finalization: Combine results ---
  # Convert the list of results into a data.table
  results_dt <- rbindlist(results_list)

  # Combine the results with the original pairs to create the final table
  final_dt <- cbind(candidate_pairs_dt, results_dt)

  return(final_dt)
}




# # Load the lookup tables
# data(sample_term_lookup)
# data(sample_gene_symbol_lookup) # Not used here, but for annotating element lists
#
# # Add new columns with the term names using the fast lookup
# scored_pairs_v1[, term1_name := sample_term_lookup[term1]]
# scored_pairs_v1[, term2_name := sample_term_lookup[term2]]
#
# # Save the annotated data frame to a file
# final_scored_pairs_v1 <- scored_pairs_v1[p_value < 0.0001]
# data.table::fwrite(final_scored_pairs_v1, "data-raw/scored_pairs_annotated.tsv", sep = "\t")
#
# # View the result
# print(head(scored_pairs_v1))
#

# bipartite_edges <- as.data.table(to_dataframe(sample_ects))
# disjoint_sets_dt <- get_disjoint_sets(sample_ects@ecprob,
#                                        candidate_pairs_v3,
#                                        bipartite_edges)

get_disjoint_sets_v1_slow <- function(pairs_dt, bipartite_edges) {

  # Setup: Create a fast term -> element lookup list
  bipartite_edges[, term := as.character(term)]
  term_to_elements_list <- split(bipartite_edges$element, bipartite_edges$term)

  # Pre-allocate lists to store the results
  results1 <- vector("list", nrow(pairs_dt))
  results2 <- vector("list", nrow(pairs_dt))

  for (i in 1:nrow(pairs_dt)) {
    term1 <- pairs_dt$term1[i]
    term2 <- pairs_dt$term2[i]

    elements1 <- term_to_elements_list[[term1]]
    elements2 <- term_to_elements_list[[term2]]

    results1[[i]] <- setdiff(elements1, elements2)
    results2[[i]] <- setdiff(elements2, elements1)
  }

  # Combine results into the final data.table
  final_dt <- copy(pairs_dt)
  final_dt[, `:=`(
    elements1_disjoint = results1,
    elements2_disjoint = results2
  )]

  return(final_dt)
}


# Fast version: A wrapper around the new S4 method
get_disjoint_sets_v2_fast <- function(pairs_dt, bipartite_edges) {
  # This directly calls the new S4 method we are testing
  get_disjoint_sets(sample_ects@ecprob, pairs_dt, bipartite_edges)
}

get_disjoint_sets_v3_fast <- function(pairs_dt, bipartite_edges) {

  # Create a unique ID for each pair for grouping
  pairs_with_id <- pairs_dt[, .(pair_id = .I, term1, term2)]

  # Melt the pairs table into a long format
  long_pairs <- melt(pairs_with_id,
                     id.vars = "pair_id",
                     measure.vars = c("term1", "term2"),
                     value.name = "term")

  # Join to get all elements for all terms in all pairs
  setkey(long_pairs, term)
  setkey(bipartite_edges, term)
  all_elements_by_pair <- bipartite_edges[long_pairs, on = "term", allow.cartesian = TRUE]

  # Find the INTERSECTING elements for each pair
  intersecting_elements <- all_elements_by_pair[, .N, by = .(pair_id, element)][N == 2]
  setkey(intersecting_elements, pair_id, element)

  # Use an ANTI-JOIN to find the DISJOINT elements
  setkey(all_elements_by_pair, pair_id, element)
  disjoint_elements <- all_elements_by_pair[!intersecting_elements]

  # --- NEW: Robustly aggregate the results back ---
  # Aggregate the disjoint elements for term1
  disjoint1 <- disjoint_elements[variable == "term1", .(elements1_disjoint = list(element)), by = pair_id]
  # Aggregate the disjoint elements for term2
  disjoint2 <- disjoint_elements[variable == "term2", .(elements2_disjoint = list(element)), by = pair_id]

  # Perform a full join to merge the two sets of results
  setkey(disjoint1, pair_id)
  setkey(disjoint2, pair_id)
  merged_disjoint <- merge(disjoint1, disjoint2, by = "pair_id", all = TRUE)

  # Join back with the original pairs table to get the final output
  setkey(pairs_with_id, pair_id)
  final_output <- merged_disjoint[pairs_with_id, on = "pair_id"]

  # Replace any NA list-columns (from pairs with no disjoint sets) with an empty list
  final_output[is.na(elements1_disjoint), elements1_disjoint := list(list(character(0)))]
  final_output[is.na(elements2_disjoint), elements2_disjoint := list(list(character(0)))]

  return(final_output[, .(term1, term2, elements1_disjoint, elements2_disjoint)])
}

#' @title Vectorized Setdiff (Step 1: Reshape Data)
#'
#' @description This is the first step in a high-performance, vectorized
#' setdiff operation. It takes a "wide" data.table of term pairs and a
#' term-element mapping and reshapes them into a single "long" or "tidy"
#' data.table. This format is the prerequisite for fast, join-based set
#' operations.
#'
#' @param pairs_dt A `data.table` with two columns ("term1", "term2").
#' @param bipartite_edges A `data.table` with "term" and "element" columns.
#'
#' @return A "long" data.table with the columns: `pair_id`, `variable` (which
#'   indicates if the term was from the original "term1" or "term2" column),
#'   `term`, and `element`.
vectorized_setdiff_step1 <- function(pairs_dt, bipartite_edges) {

  # --- 1a. Create a unique ID for each pair for grouping ---
  # The .I special symbol gives the row number, creating a unique ID.
  pairs_with_id <- pairs_dt[, .(pair_id = .I, term1, term2)]
  print("-- 1 -- pairs_with_id")
  print(pairs_with_id)

  # --- 1b. Melt the "wide" pairs table into a "long" format ---
  # This converts the term1 and term2 columns into key-value rows.
  long_pairs <- data.table::melt(
    pairs_with_id,
    id.vars = "pair_id",
    measure.vars = c("term1", "term2"),
    value.name = "term"
  )
  print("-- 2 -- long_pairs")
  print(long_pairs)
  # The 'variable' column now helpfully tracks whether the term was
  # originally a 'term1' or a 'term2'.

  # --- 1c. Join to get all elements for all terms in all pairs ---
  # This is the final step that brings in the element information.
  data.table::setkey(long_pairs, term)
  data.table::setkey(bipartite_edges, term)

  all_elements_by_pair <- bipartite_edges[long_pairs,
                                          on = "term",
                                          allow.cartesian = TRUE]
  print("-- 3 -- all_elements_by_pair")
  print(all_elements_by_pair)

  return(all_elements_by_pair)
}

vectorized_setdiff_step2 <- function(long_dt) {

  # Group by the pair_id and the element. If a group has a size of 2,
  # it means that element was present for both term1 and term2 for that pair.
  # This single, highly optimized command finds all intersections at once.
  with_element_count_in_long_dt <- long_dt[, .N, by = .(pair_id, element)]
  print("-- 4 -- with_element_count_in_long_dt")
  print(with_element_count_in_long_dt)
  intersecting_elements <- with_element_count_in_long_dt[,][N == 2]
  print("-- 5 -- intersecting_elements")
  print(intersecting_elements)
  # We only need the pair_id and element columns for the result.
  print("-- 6 -- only pair_id and element")
  print(intersecting_elements[, .(pair_id, element)])
  return(intersecting_elements[, .(pair_id, element)])
}

vectorized_setdiff_step3 <- function(long_dt, intersecting_dt) {

  # To perform a fast anti-join, both tables must have keys set
  # on the columns we are joining by.
  data.table::setkey(long_dt, pair_id, element)
  data.table::setkey(intersecting_dt, pair_id, element)

  # --- The Anti-Join ---
  # The `!` operator in a data.table join means "not in".
  # This single, highly optimized command selects all rows from `long_dt`
  # that do NOT have a match in `intersecting_dt`. This is the vectorized
  # equivalent of a set difference.
  disjoint_elements <- long_dt[!intersecting_dt]
  print("-- 7 -- disjoint_elements")
  print(disjoint_elements)

  return(disjoint_elements)
}

vectorized_setdiff_step4 <- function(disjoint_dt, pairs_dt) {

  # --- 4a. Aggregate the disjoint elements into lists ---
  # Group the long table by pair_id and the original variable ('term1' or 'term2')
  # and aggregate the elements into a list for each.
  aggregated_disjoint <- disjoint_dt[, .(disjoint_set = list(element)), by = .(pair_id, variable)]

  # --- 4b. Separate the results for term1 and term2 ---
  disjoint1 <- aggregated_disjoint[variable == "term1", .(pair_id, elements1_disjoint = disjoint_set)]
  disjoint2 <- aggregated_disjoint[variable == "term2", .(pair_id, elements2_disjoint = disjoint_set)]

  # --- 4c. Perform a full join to merge the two sets ---
  # A full join (`all = TRUE`) is crucial. It correctly handles cases where a
  # pair_id exists in one table but not the other (i.e., one of the disjoint sets is empty).
  data.table::setkey(disjoint1, pair_id)
  data.table::setkey(disjoint2, pair_id)
  merged_disjoint <- merge(disjoint1, disjoint2, by = "pair_id", all = TRUE)

  # --- 4d. Join back with the original pairs table ---
  # This ensures we have a row for every original pair, even those with no disjoint elements.
  pairs_with_id <- pairs_dt[, .(pair_id = .I, term1, term2)]
  data.table::setkey(pairs_with_id, pair_id)
  final_output <- merged_disjoint[pairs_with_id, on = "pair_id"]

  # --- 4e. Clean up NA values ---
  # The full join will create NA for list-columns where a disjoint set was empty.
  # We replace these with a properly typed empty list (`list(character(0))`).
  final_output[sapply(elements1_disjoint, is.null), elements1_disjoint := list(list(character(0)))]
  final_output[sapply(elements2_disjoint, is.null), elements2_disjoint := list(list(character(0)))]

  return(final_output[, .(term1, term2, elements1_disjoint, elements2_disjoint)])
}

# Create a small test case
test_pairs <- data.table(
  term1 = c("GO:A", "GO:C"),
  term2 = c("GO:B", "GO:D")
)
bipartite_edges <- data.table(
  term = c("GO:A", "GO:A", "GO:B", "GO:C", "GO:D"),
  element = c("E1", "E2", "E2", "E3", "E4")
)

# Run step1
long_dt <- vectorized_setdiff_step1(test_pairs, bipartite_edges)
# Run step2
intersecting_dt <- vectorized_setdiff_step2(long_dt)
# Run step3
disjoint_dt <- vectorized_setdiff_step3(long_dt, intersecting_dt)
# Run step4
final_wide_dt <- vectorized_setdiff_step4(disjoint_dt, test_pairs)
print("-- 8 -- final_wide_dt")
print(final_wide_dt)
str(final_wide_dt)

get_disjoint_sets_fast <- function(pairs_dt, bipartite_edges) {

  # --- STEP 1: Reshape data into a long format ---
  pairs_with_id <- pairs_dt[, .(pair_id = .I, term1, term2)]

  long_pairs <- melt(pairs_with_id,
                     id.vars = "pair_id",
                     measure.vars = c("term1", "term2"),
                     value.name = "term")

  setkey(long_pairs, term)
  setkey(bipartite_edges, term)
  all_elements_by_pair <- bipartite_edges[long_pairs, on = "term", allow.cartesian = TRUE]


  # --- STEP 2: Find intersecting elements ---
  intersecting_elements <- all_elements_by_pair[, .N, by = .(pair_id, element)][N == 2]
  setkey(intersecting_elements, pair_id, element)


  # --- STEP 3: Find disjoint elements with an anti-join ---
  setkey(all_elements_by_pair, pair_id, element)
  disjoint_elements <- all_elements_by_pair[!intersecting_elements]


  # --- STEP 4: Robustly reshape back to a wide format ---
  # Aggregate the disjoint elements for term1 and term2 separately
  disjoint1 <- disjoint_elements[variable == "term1", .(elements1_disjoint = list(element)), by = pair_id]
  disjoint2 <- disjoint_elements[variable == "term2", .(elements2_disjoint = list(element)), by = pair_id]

  # Join these results back to the original pairs table.
  # This is a more robust method than a full join of list-columns.
  setkey(pairs_with_id, pair_id)
  setkey(disjoint1, pair_id)
  setkey(disjoint2, pair_id)

  final_output <- disjoint2[disjoint1[pairs_with_id]]

  # Replace any NA list-columns (from pairs with no disjoint sets) with a typed empty list
  final_output[is.na(elements1_disjoint), elements1_disjoint := list(list(character(0)))]
  final_output[is.na(elements2_disjoint), elements2_disjoint := list(list(character(0)))]

  # Return the final, clean table
  return(final_output[, .(term1, term2, elements1_disjoint, elements2_disjoint)])
}


# --- THE FINAL BENCHMARK SCRIPT ---

message("--- Loading sample data ---")
data(sample_ecg)
data(sample_ects)

message("---  Getting candidate pairs ---")

# start_time <- Sys.time()
# candidate_pairs_v1 <- get_candidate_pairs_v1()
# elapsed_time_v1 <- as.numeric(Sys.time()-start_time)
# print(paste("v1:", paste(round(elapsed_time_v1, 2), "seconds")))
#
# start_time <- Sys.time()
# candidate_pairs_v2 <- get_candidate_pairs_v2()
# elapsed_time_v2 <- as.numeric(Sys.time()-start_time)
# print(paste("v2:", paste(round(elapsed_time_v2, 2), "seconds")))

start_time <- Sys.time()
candidate_pairs_v3 <- get_candidate_pairs_v3()
elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
message(paste("getting pairs v3:", paste(round(elapsed_time, 2), "seconds")))
#
# # data.table::setorder(candidate_pairs_v1, term1, term2)
# # data.table::setorder(candidate_pairs_v2, term1, term2)
# data.table::setorder(candidate_pairs_v3, term1, term2)

# comparison_result <- all.equal(candidate_pairs_v1, candidate_pairs_v3)
#
# if (isTRUE(comparison_result)) {
#   message("SUCCESS: The outputs of v1 and v3 are identical.")
# } else {
#   message("FAILURE: The outputs are different. Details below:")
#   print(comparison_result)
# }

benchmark_file <- "data-raw/scored_pairs_v1.rds"

if (!file.exists(benchmark_file)) {
  message("--- Scoring pairs with v1 and saving to file ---")
  start_time_v1 <- Sys.time()
  scored_pairs_v1 <- score_pairs_v1(candidate_pairs_v3)
  end_time_v1 <- Sys.time()
  ellapsed_time <- as.numeric(end_time_v1 - start_time_v1, units = "secs")
  message(paste("scoring pairs v1:", round(ellapsed_time, 2), "seconds"))
  saveRDS(scored_pairs_v1, benchmark_file)
} else {
  message("--- Loading pre-computed v1 benchmark results from file... ---")
  scored_pairs_v1 <- readRDS(benchmark_file)
}

# --- Running the Comparison ---

# Take a subset of the candidate pairs for a quick but meaningful test
test_pairs <- head(candidate_pairs_v3, 1000)

# message("\n--- Benchmarking get_disjoint_sets ---")
#
# message("Running v1 (slow but safe loop)...")
# start_time_v1 <- Sys.time()
# disjoint_v1 <- get_disjoint_sets_v1_slow(test_pairs, bipartite_edges)
#
#
# message("Running v3 (fast S4 method)...")
# start_time_v3 <- Sys.time()
# disjoint_v3 <- get_disjoint_sets_v3_fast(test_pairs, bipartite_edges)
#
#
# message("\n--- Comparing v1 and v3 outputs for identity ---")
# # Sort both results to ensure a fair comparison
# setorder(disjoint_v1, term1, term2)
# setorder(disjoint_v3, term1, term2)
#
# comparison_result <- all.equal(disjoint_v1, disjoint_v3, check.attributes = FALSE)
#
# if (isTRUE(comparison_result)) {
#   message("SUCCESS: The outputs of the slow and fast methods are identical.")
# } else {
#   message("FAILURE: The outputs are different. Details below:")
#   print(comparison_result)
# }


message("\n--- Benchmarking get_disjoint_sets ---")

bipartite_edges <- as.data.table(to_dataframe(sample_ects))
test_pairs <- head(candidate_pairs_v3, 100000)


message("Running v1 (slow but safe loop)...")
start_time_v1 <- Sys.time()
disjoint_v1 <- get_disjoint_sets_v1_slow(test_pairs, bipartite_edges)
end_time_v1 <- Sys.time()
time_diff_v1 <- as.numeric(end_time_v1 - start_time_v1, units = "secs")
print(paste("v1 time:", round(time_diff_v1, 4), "seconds"))


message("Running the final fast vectorized method...")
start_time_fast <- Sys.time()
disjoint_fast <- get_disjoint_sets_fast(test_pairs, bipartite_edges)
end_time_fast <- Sys.time()
time_diff_fast <- as.numeric(end_time_fast - start_time_fast, units = "secs")
print(paste("v3 time:", round(time_diff_fast, 4), "seconds"))


message("\n--- Comparing v1 and fast method outputs for identity ---")
setorder(disjoint_v1, term1, term2)
setorder(disjoint_fast, term1, term2)

# Custom comparison for list-columns
are_identical <- all(sapply(1:nrow(disjoint_v1), function(i) {
  setequal(disjoint_v1$elements1_disjoint[[i]], disjoint_fast$elements1_disjoint[[i]]) &&
    setequal(disjoint_v1$elements2_disjoint[[i]], disjoint_fast$elements2_disjoint[[i]])
}))

if (isTRUE(are_identical)) {
  message("SUCCESS: The fast method is correct and produces identical output.")
} else {
  message("FAILURE: The outputs are different.")
}
