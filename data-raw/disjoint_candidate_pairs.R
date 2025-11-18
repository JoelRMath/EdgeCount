library(EdgeCount)
library(data.table)

# --- Helper Function (Included for standalone compatibility) ---
get_disjoint_sets <- function(pairs_dt, set_membership_dt) {
  pairs_with_id <- pairs_dt[, .(pair_id = .I, set1, set2)]
  long_pairs <- melt(pairs_with_id, id.vars = "pair_id", measure.vars = c("set1", "set2"), value.name = "set_id")
  setkey(long_pairs, set_id)
  setkey(set_membership_dt, set_id)
  all_elements_by_pair <- set_membership_dt[long_pairs, on = "set_id", allow.cartesian = TRUE]
  intersecting_elements <- all_elements_by_pair[, .N, by = .(pair_id, element)][N == 2]
  setkey(intersecting_elements, pair_id, element)
  setkey(all_elements_by_pair, pair_id, element)
  disjoint_elements <- all_elements_by_pair[!intersecting_elements]
  disjoint1 <- disjoint_elements[variable == "set1", .(elements1_disjoint = list(element)), by = pair_id]
  disjoint2 <- disjoint_elements[variable == "set2", .(elements2_disjoint = list(element)), by = pair_id]
  setkey(disjoint1, pair_id)
  setkey(disjoint2, pair_id)
  merged_disjoint <- merge(disjoint1, disjoint2, by = "pair_id", all = TRUE)
  setkey(pairs_with_id, pair_id)
  final_output <- merged_disjoint[pairs_with_id, on = "pair_id"]
  final_output[sapply(elements1_disjoint, is.null), elements1_disjoint := list(list(character(0)))]
  final_output[sapply(elements2_disjoint, is.null), elements2_disjoint := list(list(character(0)))]
  return(final_output[, .(set1, set2, elements1_disjoint, elements2_disjoint)])
}


#' @title Get Candidate Pairs (Original)
#' @description Finds all pairs of terms that are connected by at least one
#' edge in the interaction graph. This is a fast, approximate "superset" query.
#' @param ecg An ECGraph/ECProb object.
#' @param ects An ECTermScoring object.
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


get_disjoint_connected_pairs_proto <- function(ecp, ects) {

  message("Step 1: Mapping graph edges to term pairs...")
  # 1. Prepare Edge Lists
  network_edges <- data.table(to_dataframe(ecp))
  setnames(network_edges, c("from", "to"), c("element1", "element2"))

  bipartite_edges <- as.data.table(to_dataframe(ects))
  bipartite_edges[, `:=`(term = as.character(term), element = as.character(element))]

  # 2. Map element edges to term pairs (Create the 'merged2' table)
  # This table contains every instance of an edge connecting Term1 and Term2
  merged1 <- network_edges[bipartite_edges, on = .(element1 = element), nomatch = 0, allow.cartesian = TRUE]
  setnames(merged1, "term", "term1")

  merged2 <- merged1[bipartite_edges, on = .(element2 = element), nomatch = 0, allow.cartesian = TRUE]
  setnames(merged2, "term", "term2")

  # Filter out self-term loops immediately to reduce size
  merged2 <- merged2[term1 != term2]

  # --- THE OPTIMIZATION: Filter for Disjointness ---
  message("Step 2: Filtering for disjoint connections using anti-joins...")

  # We have rows: element1 -- element2, where e1 in term1, e2 in term2.
  # A connection is disjoint if: e1 NOT in term2  AND  e2 NOT in term1.

  # To use data.table joins efficiently, we need keys
  setkey(bipartite_edges, element, term)

  # Check 1: Is element1 also in term2?
  # We perform a join. If a match is found, it's an intersection (bad).
  # We flag these rows.
  # Note: merged2 has (element1, term2). Bipartite has (element, term).
  ids_shared1 <- merged2[bipartite_edges, on = .(element1 = element, term2 = term),
                         nomatch = 0, which = TRUE] # 'which=TRUE' returns indices

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

  message("Step 3: Aggregating unique pairs...")

  # 4. Canonicalize pairs and count observed edges
  # We use pmin/pmax to ensure (A, B) and (B, A) are treated as the same pair.
  final_results <- valid_edges[, .(
    set1 = pmin(term1, term2),
    set2 = pmax(term1, term2)
  )]

  # Count the number of disjoint edges for each pair
  final_results <- final_results[, .(observed_edges = .N), by = .(set1, set2)]

  return(final_results)
}

# ---
# EXAMPLE USAGE
# ---
message("--- Loading sample data ---")
data(sample_ecg)
data(sample_ects)
ecp <- ECProb(sample_ecg)

message("\n--- Running the original, fast 'superset' method ---")
system.time(candidates <- get_candidate_pairs(ecp, sample_ects))
message(paste("Found", nrow(candidates), "total candidate pairs."))

message("\n--- Running the new, rigorous 'disjoint' method ---")
system.time(disjoint_pairs <- get_disjoint_connected_pairs_proto(ecp, sample_ects))
message(paste("Found", nrow(disjoint_pairs), "pairs with true disjoint connections."))

print(head(disjoint_pairs))
