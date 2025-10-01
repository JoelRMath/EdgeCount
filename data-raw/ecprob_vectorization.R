library(EdgeCount)
library(data.table)

# ---
# FUNCTION DEFINITIONS
# ---

#' Get all connected term pairs and their observed edge counts (Fast Method)
get_pairs_with_counts <- function() {
  net_edges <- data.table(to_dataframe(sample_ecg))
  setnames(net_edges, c("from", "to"), c("element1", "element2"))

  bip_edges <- as.data.table(to_dataframe(sample_ects))
  bip_edges[, `:=`(term = as.character(term), element = as.character(element))]

  merged1 <- net_edges[bip_edges, on = .(element1 = element), nomatch = 0, allow.cartesian = TRUE]
  setnames(merged1, "term", "term1")
  merged2 <- merged1[bip_edges, on = .(element2 = element), nomatch = 0, allow.cartesian = TRUE]
  setnames(merged2, "term", "term2")

  pairs_with_counts <- merged2[term1 != term2,
                               .(
                                 term_canon_1 = pmin(term1, term2), term_canon_2 = pmax(term1, term2),
                                 edge_canon_1 = pmin(element1, element2), edge_canon_2 = pmax(element1, element2)
                               )
  ][,
    unique(.SD, by = c("term_canon_1", "term_canon_2", "edge_canon_1", "edge_canon_2"))
  ][,
    .(observed_edges = .N), by = .(term1 = term_canon_1, term2 = term_canon_2)
  ]
  return(pairs_with_counts)
}

#' Gold Standard: A slow but safe for-loop implementation
get_disjoint_sets_slow <- function(pairs_dt, bipartite_edges) {
  bipartite_edges[, term := as.character(term)]
  term_to_elements_list <- split(bipartite_edges$element, bipartite_edges$term)
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
  final_dt <- copy(pairs_dt)
  final_dt[, `:=`(elements1_disjoint = results1, elements2_disjoint = results2)]
  return(final_dt)
}

#' High-performance, fully vectorized implementation
get_disjoint_sets_fast_v0 <- function(pairs_dt, bipartite_edges) {
  pairs_with_id <- pairs_dt[, .(pair_id = .I, term1, term2)]

  long_pairs <- melt(pairs_with_id,
                     id.vars = "pair_id",
                     measure.vars = c("term1", "term2"),
                     value.name = "term")

  setkey(long_pairs, term)
  setkey(bipartite_edges, term)
  all_elements_by_pair <- bipartite_edges[long_pairs, on = "term", allow.cartesian = TRUE]

  intersecting_elements <- all_elements_by_pair[, .N, by = .(pair_id, element)][N == 2]
  setkey(intersecting_elements, pair_id, element)

  setkey(all_elements_by_pair, pair_id, element)
  disjoint_elements <- all_elements_by_pair[!intersecting_elements]

  disjoint1 <- disjoint_elements[variable == "term1", .(elements1_disjoint = list(element)), by = pair_id]
  disjoint2 <- disjoint_elements[variable == "term2", .(elements2_disjoint = list(element)), by = pair_id]

  setkey(disjoint1, pair_id)
  setkey(disjoint2, pair_id)
  merged_disjoint <- merge(disjoint1, disjoint2, by = "pair_id", all = TRUE)

  setkey(pairs_with_id, pair_id)
  final_output <- merged_disjoint[pairs_with_id, on = "pair_id"]

  final_output[is.na(elements1_disjoint), elements1_disjoint := list(list(character(0)))]
  final_output[is.na(elements2_disjoint), elements2_disjoint := list(list(character(0)))]

  return(final_output[, .(term1, term2, elements1_disjoint, elements2_disjoint)])
}


# ---
# MAIN BENCHMARK SCRIPT
# ---

message("--- Loading sample data ---")
data(sample_ecg)
data(sample_ects)

message("--- Getting candidate pairs ---")
candidate_pairs <- get_pairs_with_counts()
# Use a large subset for a meaningful benchmark
print(candidate_pairs)
test_pairs <- head(candidate_pairs, 100000)
bipartite_edges <- as.data.table(to_dataframe(sample_ects))

message("\n--- Benchmarking get_disjoint_sets ---")

message("Running slow but safe loop...")
start_time_slow <- Sys.time()
disjoint_slow <- get_disjoint_sets_slow(test_pairs, bipartite_edges)
end_time_slow <- Sys.time()
time_diff_slow <- as.numeric(end_time_slow - start_time_slow, units = "secs")
print(paste("Slow method time:", round(time_diff_slow, 4), "seconds"))


message("Running fast vectorized method...")
start_time_fast <- Sys.time()
disjoint_fast <- get_disjoint_sets_fast_v0(test_pairs, bipartite_edges)
end_time_fast <- Sys.time()
time_diff_fast <- as.numeric(end_time_fast - start_time_fast, units = "secs")
print(paste("Fast method time:", round(time_diff_fast, 4), "seconds"))


message("\n--- Comparing outputs for identity ---")
setorder(disjoint_slow, term1, term2)
setorder(disjoint_fast, term1, term2)

# Use a custom comparison function for list-columns
are_identical <- all(sapply(1:nrow(disjoint_slow), function(i) {
  setequal(disjoint_slow$elements1_disjoint[[i]], disjoint_fast$elements1_disjoint[[i]]) &&
    setequal(disjoint_slow$elements2_disjoint[[i]], disjoint_fast$elements2_disjoint[[i]])
}))

if (isTRUE(are_identical)) {
  message("SUCCESS: The outputs of the slow and fast methods are identical.")
} else {
  message("FAILURE: The outputs of the slow and fast methods are different.")
}

################## LAMBDA ###############

# --- 1. Setup: Create sample data ---

# A list where each element represents a term and contains the degrees of its members.
# This mimics the structure you described.
term_degree_list <- list(
  TermA = c(10, 5, 2),      # size = 3
  TermB = c(8, 1),          # size = 2
  TermC = c(20, 15, 1, 4)   # size = 4
)

# The single, long vector of all degrees, sorted by term (your vector K)
K <- unlist(term_degree_list, use.names = FALSE)
# K is now: 10, 5, 2, 8, 1, 20, 15, 1, 4

# The main cumulative sum of degrees (your sK)
sK <- cumsum(K)
# sK is now: 10, 15, 17, 25, 26, 46, 61, 62, 66


# --- 2. The High-Performance Vectorized Calculation ---

# Get the size of each term block (fast)
term_sizes <- lengths(term_degree_list) # c(3, 2, 4)

# Get the end index of each block in the long vector (your sL, also fast)
end_indices <- cumsum(term_sizes) # c(3, 5, 9)

# Get the cumulative sum values at the end of each block (fast)
sK_at_ends <- sK[end_indices] # c(17, 26, 66)

# Calculate the sum for each block using a lagged difference (fast)
# The first sum is just the first value. The subsequent sums are the differences.
target_sums <- c(sK_at_ends[1], diff(sK_at_ends))


# --- 3. Verification ---
# Let's check against a slower but more obvious method
gold_standard_sums <- sapply(term_degree_list, sum)

# Compare the results
print("Results from the fast vectorized method:")
print(target_sums)
#> [1] 17  9 40

print("Results from the 'gold standard' sapply method:")
print(gold_standard_sums)
#> TermA TermB TermC
#>    17     9    40

# They are identical.
print(all.equal(unname(gold_standard_sums), target_sums))
#> [1] TRUE
