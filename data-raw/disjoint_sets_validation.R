library(EdgeCount)
library(data.table)

# --- 1. FUNCTION DEFINITIONS (The Code to Test) ---

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

get_disjoint_connected_pairs_optimized <- function(ecp, ects) {

  # 1. Prepare Edge Lists
  network_edges <- data.table(to_dataframe(ecp))
  setnames(network_edges, c("from", "to"), c("element1", "element2"))

  bipartite_edges <- as.data.table(to_dataframe(ects))
  bipartite_edges[, `:=`(term = as.character(term), element = as.character(element))]

  # 2. Map element edges to term pairs (Create the 'merged2' table)
  merged1 <- network_edges[bipartite_edges, on = .(element1 = element), nomatch = 0, allow.cartesian = TRUE]
  setnames(merged1, "term", "term1")

  merged2 <- merged1[bipartite_edges, on = .(element2 = element), nomatch = 0, allow.cartesian = TRUE]
  setnames(merged2, "term", "term2")

  # Filter out self-term loops
  merged2 <- merged2[term1 != term2]

  # --- THE OPTIMIZATION: Filter for Disjointness using Anti-Joins ---

  # To use data.table joins efficiently, we need keys
  setkey(bipartite_edges, element, term)

  # Check 1: Is element1 also in term2?
  ids_shared1 <- merged2[bipartite_edges, on = .(element1 = element, term2 = term),
                         nomatch = 0, which = TRUE]

  # Check 2: Is element2 also in term1?
  ids_shared2 <- merged2[bipartite_edges, on = .(element2 = element, term1 = term),
                         nomatch = 0, which = TRUE]

  # Combine indices of all "bad" (non-disjoint) edges
  bad_indices <- unique(c(ids_shared1, ids_shared2))

  # 3. Create the clean list of disjoint edges
  if (length(bad_indices) > 0) {
    valid_edges <- merged2[-bad_indices]
  } else {
    valid_edges <- merged2
  }

  # 4. Canonicalize pairs and count observed edges
  final_results <- valid_edges[, .(
    set1 = pmin(term1, term2),
    set2 = pmax(term1, term2)
  )]

  # Count the number of disjoint edges for each pair
  final_results <- final_results[, .(observed_edges = .N), by = .(set1, set2)]

  return(final_results)
}


# --- 2. GOLD STANDARD (Slow, Loop-Based Reference) ---

get_disjoint_connected_pairs_reference <- function(ecp, ects) {

  # Get candidates first (we assume this part is correct/safe as it is a simple join)
  candidate_pairs <- get_candidate_pairs(ecp, ects)

  # Prepare lookups
  bipartite_edges <- as.data.table(to_dataframe(ects))
  term_members <- split(bipartite_edges$element, bipartite_edges$term)

  results_list <- vector("list", nrow(candidate_pairs))

  # Iterate through every candidate pair and manually verify
  ct <- 0
  for (i in 1:nrow(candidate_pairs)) {
    ct <- ct +1
    if (ct == 100000){
      print(paste0(i,"/",nrow(candidate_pairs)))
      ct <- 0
    }
    t1 <- candidate_pairs$set1[i]
    t2 <- candidate_pairs$set2[i]

    el1 <- term_members[[t1]]
    el2 <- term_members[[t2]]

    # Explicitly calculate disjoint sets
    el1_d <- setdiff(el1, el2)
    el2_d <- setdiff(el2, el1)

    if (length(el1_d) == 0 || length(el2_d) == 0) {
      next # No disjoint elements, skip
    }

    # Use trusted package function to count edges
    obs <- get_edge_count_between(ecp, el1_d, el2_d)

    if (obs > 0) {
      results_list[[i]] <- data.table(set1 = t1, set2 = t2, observed_edges = obs)
    }
  }

  final_dt <- rbindlist(results_list)
  if (nrow(final_dt) > 0) {
    # Canonicalize set1/set2 just to be safe for comparison
    final_dt[, `:=`(set1_canon = pmin(set1, set2), set2_canon = pmax(set1, set2))]
    final_dt[, `:=`(set1 = set1_canon, set2 = set2_canon, set1_canon=NULL, set2_canon=NULL)]
  }

  return(final_dt)
}


# --- 3. EXECUTION AND VERIFICATION ---

message("--- Loading sample data ---")
data(sample_ecg)
data(sample_ects)
ecp <- ECProb(sample_ecg)

message("\n--- Running Benchmarks on Full Data ---")

message("1. Running Optimized Method...")
t1 <- system.time(optimized_res <- get_disjoint_connected_pairs_optimized(ecp, sample_ects))
print(t1)

message("2. Running Reference Method (this might take a moment)...")
# For speed in testing, let's just check consistency on the output of optimized
# vs reference.
t2 <- system.time(reference_res <- get_disjoint_connected_pairs_reference(ecp, sample_ects))
print(t2)

message("\n--- Comparison ---")
# Harmonize for comparison
setorder(optimized_res, set1, set2)
setorder(reference_res, set1, set2)

# Ensure columns are in same order
col_order <- c("set1", "set2", "observed_edges")
optimized_res <- optimized_res[, ..col_order]
reference_res <- reference_res[, ..col_order]

if (all.equal(optimized_res, reference_res)) {
  message("SUCCESS: Optimized method matches Reference method exactly on sample_ects.")
} else {
  message("FAILURE: Mismatch found on sample_ects.")
  print(all.equal(optimized_res, reference_res))
}


# --- 4. Edge Case Test (Intersection Only) ---
message("\n--- Running Edge Case Test (Intersection Only) ---")
# Create a scenario where two terms share an element 'B'.
# 'B' is connected to 'A' (in Term1) and 'C' (in Term2).
# BUT 'A' is NOT connected to 'C'.
# There are NO edges between the disjoint parts ({A} and {C}).
# The only path is A-B-C, which goes through the intersection.
# This pair should NOT appear in the result.

edge_df_toy <- data.frame(p1 = c("A", "B"), p2 = c("B", "C"), stringsAsFactors = FALSE)
ecp_toy <- ECProb(ECGraph(edge_df_toy))

term_df_toy <- data.frame(
  term = c("Term1", "Term1", "Term2", "Term2"),
  element = c("A", "B", "B", "C"), # B is the intersection
  stringsAsFactors = FALSE
)
ects_toy <- ECTermScoring(term_df_toy)

toy_res <- get_disjoint_connected_pairs_optimized(ecp_toy, ects_toy)

if (nrow(toy_res) == 0) {
  message("SUCCESS: Edge case handled correctly. Pair connected only via intersection was filtered out.")
} else {
  message("FAILURE: Edge case failed. Returned rows:")
  print(toy_res)
}
