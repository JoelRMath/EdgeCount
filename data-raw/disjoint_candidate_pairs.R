library(EdgeCount)
library(data.table)

#' @title Get Candidate Pairs (Original)
#' @description Finds all pairs of terms that are connected by at least one
#' edge in the interaction graph. This is a fast, approximate "superset" query
#' and may include pairs connected only via their intersection.
#' @param ecg An ECGraph/ECProb object (the interaction network).
#' @param ects An ECTermScoring object (the bipartite graph).
#' @return A data.table with "set1" and "set2" columns.
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


#' @title Get Disjoint Observed Edges (NEW Lightweight Helper)
#' @description Calculates *only* the observed edge count between the
#' disjoint parts of set pairs.
#' @param object An ECProb object.
#' @param pairs_dt A data table of pairs ("set1", "set2").
#' @param set_membership_dt A data table mapping "set_id" to "element".
#' @return A data.table with "set1", "set2", and "observed_edges".
get_disjoint_observed_edges <- function(object, pairs_dt, set_membership_dt) {

  # --- Step 1: Get Disjoint Sets ---
  # This helper is defined in the ECProb.R file
  disjoint_sets_dt <- get_disjoint_sets(pairs_dt, set_membership_dt)
  disjoint_sets_dt[, pair_id := .I]

  # --- Step 2: Get Graph Edges ---
  ecg_edges <- data.table(to_dataframe(object))
  setnames(ecg_edges, c("from", "to"), c("e1", "e2"))
  ecg_edges[, `:=`(canon1 = pmin(e1, e2), canon2 = pmax(e1, e2))]
  ecg_edges <- unique(ecg_edges, by = c("canon1", "canon2"))
  setkey(ecg_edges, canon1, canon2)

  # --- Step 3: Find all *possible* edges between disjoint sets ---
  long_disjoint1 <- disjoint_sets_dt[, .(element = unlist(elements1_disjoint)), by = pair_id]
  long_disjoint2 <- disjoint_sets_dt[, .(element = unlist(elements2_disjoint)), by = pair_id]

  setkey(long_disjoint1, pair_id)
  setkey(long_disjoint2, pair_id)

  possible_edges <- long_disjoint1[long_disjoint2, on = "pair_id", allow.cartesian = TRUE, nomatch=0]

  if (nrow(possible_edges) == 0) {
    # No possible edges, return an empty table
    return(data.table(set1=character(), set2=character(), observed_edges=integer()))
  }

  setnames(possible_edges, c("element", "i.element"), c("element1", "element2"))
  possible_edges[, `:=`(canon1 = pmin(element1, element2), canon2 = pmax(element1, element2))]

  # --- Step 4: Join with *actual* graph edges ---
  observed_edges_long <- ecg_edges[possible_edges, on = .(canon1, canon2), nomatch = 0]

  # --- Step 5: Count observed edges per pair ---
  observed_edges_dt <- observed_edges_long[, .(observed_edges = .N), by = pair_id]

  # --- Step 6: Join back to original pair names ---
  pairs_with_id <- pairs_dt[, .(pair_id = .I, set1, set2)]
  setkey(pairs_with_id, pair_id)
  setkey(observed_edges_dt, pair_id)

  final_results <- observed_edges_dt[pairs_with_id, on = "pair_id", nomatch = 0]

  return(final_results[, .(set1, set2, observed_edges)])
}


#' @title Get Disjoint Connected Pairs (Optimized Version)
#' @description Finds all pairs of terms that have at least one edge
#' connecting their *disjoint* parts.
#' @param ecp An ECProb object (the interaction network).
#' @param ects An ECTermScoring object (the bipartite graph).
#' @return A data.table containing the pairs ("set1", "set2") and their
#'   disjoint `observed_edges` count.
get_disjoint_connected_pairs <- function(ecp, ects) {

  # 1. Get all *potential* candidate pairs using the fast method.
  message("Step 1: Finding all potential candidate pairs...")
  candidate_pairs <- get_candidate_pairs(ecp, ects)

  # 2. Get the full term-element membership data
  set_membership_dt <- as.data.table(to_dataframe(ects))
  setnames(set_membership_dt, c("term", "element"), c("set_id", "element"))

  # 3. Run the new, *lightweight* function to get observed edge counts
  message("Step 2: Calculating observed edges for disjoint sets...")
  all_pair_stats <- get_disjoint_observed_edges(
    object = ecp,
    pairs_dt = candidate_pairs,
    set_membership_dt = set_membership_dt
  )

  # 4. Filter for pairs where a disjoint connection was actually found.
  #    This step is now implicitly handled by get_disjoint_observed_edges
  #    (it only returns pairs with at least one observed edge).
  message("Step 3: Returning pairs with observed_edges > 0.")

  return(all_pair_stats)
}

# ---
# EXAMPLE USAGE
# ---
message("--- Loading sample data ---")
data(sample_ecg)
data(sample_ects)
ecp <- ECProb(sample_ecg)

message("\n--- Running the original, fast 'superset' method ---")
# This will include pairs connected by their intersection
candidates <- get_candidate_pairs(ecp, sample_ects)
message(paste("Found", nrow(candidates), "total candidate pairs."))

message("\n--- Running the new, rigorous 'disjoint' method ---")
# This will be a smaller, more correct set
disjoint_pairs <- get_disjoint_connected_pairs(ecp, sample_ects)
message(paste("Found", nrow(disjoint_pairs), "pairs with true disjoint connections."))

print(head(disjoint_pairs))
