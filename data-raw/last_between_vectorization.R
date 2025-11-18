#' @title Get Disjoint Sets for Multiple Pairs (Vectorized)
#'
#' @description Calculates the disjoint element sets for a list of set pairs.
#' This is a high-performance, vectorized helper function useful for
#' between-set analysis. It takes pairs of sets and determines which elements
#' are unique to each set in the pair (removing the intersection).
#'
#' @param pairs_dt A `data.table` with two columns ("set1", "set2") defining
#'   the pairs to analyze.
#' @param set_membership_dt A `data.table` mapping set IDs to their elements.
#'   Must have "set_id" and "element" columns.
#'
#' @return A `data.table` with columns: "set1", "set2", "elements1_disjoint" (list),
#'   and "elements2_disjoint" (list).
#' @export
get_disjoint_sets <- function(pairs_dt, set_membership_dt) {

  pairs_with_id <- pairs_dt[, .(pair_id = .I, set1, set2)]

  long_pairs <- melt(pairs_with_id,
                     id.vars = "pair_id",
                     measure.vars = c("set1", "set2"),
                     value.name = "set_id")

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


#' @title Vectorized calculation of fast between statistics
#'
#' @description Calculates fast between statistics for a collection of set pairs
#' using a high-performance, vectorized algorithm.
#'
#' @details This function is designed for efficiency when processing many term
#' pairs at once. It uses a series of `data.table` joins and aggregations to
#' avoid slow R-level loops. Note that this method uses the "fast" lambda
#' approximation.
#'
#' @param object An ECProb object.
#' @param pairs_dt A data table of pairs of element sets, with two columns
#'   "set1" and "set2" giving the set IDs.
#' @param set_membership_dt A data table mapping set IDs to their elements.
#'   Must have "set_id" and "element" columns.
#'
#' @return A final `data.table` with the original pairs and columns for all
#'   calculated statistics (`set1`, `size1`, `set2`, `size2`, `observed_edges`,
#'   `lambda`, `p_value`, `log2_Anscombe_ratio`).
#' @export
#' @examples
#' # TBD
setGeneric("calculate_between_stats_fast_vectorized",
           function(object, pairs_dt, set_membership_dt) standardGeneric("calculate_between_stats_fast_vectorized"))

#' @describeIn calculate_between_stats_fast_vectorized Method for ECProb objects
setMethod("calculate_between_stats_fast_vectorized",
          "ECProb",
          function(object, pairs_dt, set_membership_dt) {

            # --- Step 1: Get Disjoint Sets and Intermediate Values ---
            # Call the standalone function (no object needed)
            disjoint_sets_dt <- get_disjoint_sets(pairs_dt, set_membership_dt)
            disjoint_sets_dt[, pair_id := .I]

            all_element_degrees_dt <- data.table(
              element = object@names,
              degree = unlist(object@degrees)
            )
            data.table::setkey(all_element_degrees_dt, element)

            # Unnest the lists to get long-format element tables for joining
            long_disjoint1 <- disjoint_sets_dt[, .(element = unlist(elements1_disjoint)), by = pair_id]
            data.table::setkey(long_disjoint1, element)

            # Pre-calculate degree sums for Set 1
            sums1_dt <- all_element_degrees_dt[long_disjoint1, on = "element", nomatch = 0][,
                                                                                            .(sum_degrees1 = sum(degree, na.rm = TRUE)), by = pair_id]

            long_disjoint2 <- disjoint_sets_dt[, .(element = unlist(elements2_disjoint)), by = pair_id]
            data.table::setkey(long_disjoint2, element)

            # Pre-calculate degree sums for Set 2
            sums2_dt <- all_element_degrees_dt[long_disjoint2, on = "element", nomatch = 0][,
                                                                                            .(sum_degrees2 = sum(degree, na.rm = TRUE)), by = pair_id]

            # --- THE OPTIMIZATION: Replace Cartesian Product with Edge List Join ---

            # 1. Prepare a bidirectional edge list (u -> v AND v -> u)
            #    This allows us to find connections regardless of direction in the original data.
            ecg_edges <- data.table(to_dataframe(object))
            setnames(ecg_edges, c("from", "to"), c("u", "v"))
            # Double the edges to make them bidirectional for the join
            ecg_edges_dbl <- rbind(ecg_edges, ecg_edges[, .(u=v, v=u)])
            setkey(ecg_edges_dbl, u)

            # 2. Join Set 1 elements to the Edge List
            #    This finds all neighbors of elements in Set 1.
            #    We reuse 'long_disjoint1' (cols: pair_id, element)
            #    Join condition: element == u
            set1_neighbors <- ecg_edges_dbl[long_disjoint1, on = .(u = element), nomatch = 0, allow.cartesian = TRUE]
            # Result cols: v (neighbor), pair_id

            # 3. Join Neighbors to Set 2 elements
            #    This filters the neighbors to keep ONLY those that are in Set 2 *for the same pair*.
            #    We reuse 'long_disjoint2' (cols: pair_id, element)
            #    Join condition: v == element AND pair_id == pair_id
            setkey(long_disjoint2, pair_id, element)
            # Rename 'v' to 'element' for the join
            setnames(set1_neighbors, "v", "element")

            observed_edges_long <- long_disjoint2[set1_neighbors, on = .(pair_id, element), nomatch = 0]

            # 4. Aggregate to count observed edges
            observed_edges_dt <- observed_edges_long[, .(observed_edges = .N), by = pair_id]

            # --- End Optimization ---

            all_pairs_ids <- data.table(pair_id = 1:nrow(disjoint_sets_dt))
            setkey(all_pairs_ids, pair_id)
            setkey(observed_edges_dt, pair_id)
            observed_edges_dt <- observed_edges_dt[all_pairs_ids]
            observed_edges_dt[is.na(observed_edges), observed_edges := 0L]

            # --- Step 2: Join all intermediate results to create the final table ---
            data.table::setkey(sums1_dt, pair_id)
            data.table::setkey(sums2_dt, pair_id)
            degree_sums_dt <- merge(sums1_dt, sums2_dt, by = "pair_id", all = TRUE)
            degree_sums_dt[is.na(sum_degrees1), sum_degrees1 := 0]
            degree_sums_dt[is.na(sum_degrees2), sum_degrees2 := 0]

            data.table::setkey(disjoint_sets_dt, pair_id)
            final_dt <- degree_sums_dt[disjoint_sets_dt, on = "pair_id"]

            data.table::setkey(observed_edges_dt, pair_id)
            final_dt <- observed_edges_dt[final_dt, on = "pair_id"]

            # --- Step 3: Perform the final calculations on this complete table ---
            final_dt[, `:=`(
              size1 = lengths(elements1_disjoint),
              size2 = lengths(elements2_disjoint)
            )]
            final_dt[, max_possible_edges := as.numeric(size1) * as.numeric(size2)]
            final_dt[, lambda := (sum_degrees1 * sum_degrees2) / (2 * object@graph_size)]

            final_dt[, p_value := calculate_p_value(object, observed_edges, max_possible_edges, lambda)]
            final_dt[, log2_Anscombe_ratio := 0.5 * (log2(observed_edges + 3/8) - log2(lambda + 3/8))]


            # --- Return the final columns ---
            return(final_dt[, .(set1, size1, set2, size2, observed_edges, lambda, p_value, log2_Anscombe_ratio)])
          })
