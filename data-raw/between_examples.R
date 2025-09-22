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
    if (i %% 1000 == 0) {
      message(paste("Scoring pair", i, "of", nrow(candidate_pairs_dt)))
    }

    # Get the terms for the current row
    term1 <- candidate_pairs_dt$term1[i]
    term2 <- candidate_pairs_dt$term2[i]

    # Get the element sets using the fast lookup list
    elements1 <- term_to_elements_list[[term1]]
    elements2 <- term_to_elements_list[[term2]]

    # Calculate the observed edge count
    observed <- get_edge_count_between(sample_ecg, elements1, elements2)

    # Get the full statistics
    stats <- edge_count_statistics(ecp, elements1, elements2, observed)

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

score_pairs_v2 <- function(candidate_pairs_dt) {

  # --- Setup: This part is identical to v1 ---
  ecp <- ECProb(sample_ecg)
  bipartite_edges <- as.data.table(to_dataframe(sample_ects))
  bipartite_edges[, term := as.character(term)]
  term_to_elements_list <- split(bipartite_edges$element, bipartite_edges$term)

  # --- The "Fast" Vectorized Operation ---
  # The `by = 1:nrow(...)` trick tells data.table to run the code
  # in the `j` expression for each row, which is much faster than an R loop.
  results_dt <- candidate_pairs_dt[, {

    # Get the element sets using the fast lookup list
    elements1 <- term_to_elements_list[[term1]]
    elements2 <- term_to_elements_list[[term2]]

    # Calculate the observed edge count
    observed <- get_edge_count_between(sample_ecg, elements1, elements2)

    # Get the full statistics
    stats <- edge_count_statistics(ecp, elements1, elements2, observed)

    # Return a list of the results for this row.
    list(
      observed_edges = as.integer(observed),
      p_value = stats$p_value,
      lambda = stats$lambda,
      log2_anscombe = stats$log2_Anscombe_ratio
    )
  }, by = 1:nrow(candidate_pairs_dt)] # This tells data.table to iterate row-by-row

  # The 'by' operation creates a temporary grouping column that we need to remove.
  results_dt[, `1:nrow(candidate_pairs_dt)` := NULL]

  # Combine the results with the original pairs
  final_dt <- cbind(candidate_pairs_dt, results_dt)

  return(final_dt)
}


message("--- Loading sample data ---")
data(sample_ecg)
data(sample_ects)

message("---  Getting candidate pairs ---")
# start_time <- Sys.time()
# candidate_pairs_v1 <- get_candidate_pairs_v1()
# elapsed_time_v1 <- Sys.time()-start_time
# print(paste("v1:", elapsed_time_v1))
#
# start_time <- Sys.time()
# candidate_pairs_v2 <- get_candidate_pairs_v2()
# elapsed_time_v2 <- Sys.time()-start_time
# print(paste("v2:", elapsed_time_v2))

start_time <- Sys.time()
candidate_pairs_v3 <- get_candidate_pairs_v3()
elapsed_time_v3 <- Sys.time() - start_time
print(paste("v3:", elapsed_time_v3))

# data.table::setorder(candidate_pairs_v1, term1, term2)
# data.table::setorder(candidate_pairs_v2, term1, term2)
data.table::setorder(candidate_pairs_v3, term1, term2)

# comparison_result <- all.equal(candidate_pairs_v1, candidate_pairs_v3)
#
# if (isTRUE(comparison_result)) {
#   message("SUCCESS: The outputs of v1 and v3 are identical.")
# } else {
#   message("FAILURE: The outputs are different. Details below:")
#   print(comparison_result)
# }

scored_pairs_dt <- score_pairs_v1(candidate_pairs_v3)

# Load the lookup tables
data(sample_term_lookup)
data(sample_gene_symbol_lookup) # Not used here, but for annotating element lists

# Add new columns with the term names using the fast lookup
scored_pairs_dt[, term1_name := sample_term_lookup[term1]]
scored_pairs_dt[, term2_name := sample_term_lookup[term2]]

# Save the annotated data frame to a file
final_scored_pairs_dt <- scored_pairs_dt[p_value < 0.001]
data.table::fwrite(final_scored_pairs_dt, "data-raw/scored_pairs_annotated.tsv", sep = "\t")

# View the result
print(head(scored_pairs_dt))


