#' @title ECProb S4 Class and Constructor
#'
#' @description Extends the ECGraph class to include pre-calculated properties
#' necessary for edge-count probability calculations.
#'
#' @name ECProb
#' @aliases ECProb-class
#'
#' @slot graph_size numeric. The total number of edges in the graph (M).
#' @slot graph_order numeric. The total number of vertices in the graph (N).
#' @slot average_degree numeric. The average degree of the vertices in the graph.
#' @slot adj list. (Inherited from ECGraph) The adjacency list.
#' @slot degrees list. (Inherited from ECGraph) The degree of each vertex.
#' @slot names character. (Inherited from ECGraph) The names of all vertices.
#'
#' @param ecgraph An ECGraph object.
#'
#' @return An object of class ECProb.
#' @seealso
#' The base class for graph representation: \code{\link{ECGraph}}
#'
#' @exportClass ECProb
#' @export ECProb
#' @examples
#' edge_df <- data.frame(p1 = c("A", "B", "C"), p2 = c("B", "C", "D"))
#' ecg <- ECGraph(edge_df)
#' ecp <- ECProb(ecg)
#' ecp@graph_size # 3
#'
setClass("ECProb",
         slots = list(
           graph_size = "numeric",
           graph_order = "numeric",
           average_degree = "numeric"
         ),
         contains = "ECGraph"
)

ECProb <- function(ecgraph) {

  graph_size <- sum(unlist(ecgraph@degrees))/2
  graph_order <- length(ecgraph@names)
  average_degree <- if (graph_order > 0) 2*graph_size/graph_order else 0

  new("ECProb", ecgraph, graph_size = graph_size, graph_order = graph_order,
      average_degree = average_degree)
}

# --------------------------------------------------------------------------- #
# Non-Vectorized (Single Set/Pair) Methods
# --------------------------------------------------------------------------- #

#' @title Calculate Lambda Between Two Vertex Sets (Naive Method)
#' @description Calculates the lambda parameter (expected edge count) between two
#' vertex sets using a direct, unoptimized summation.
#' @details The naive method has a time complexity of O(n*m), where n and m are
#' the sizes of the two sets. It is accurate but can be very slow for large sets.
#' Its primary purpose is for validation and speed-testing against optimized methods.
#' @param object An ECProb object.
#' @param set1 A character vector of vertex names.
#' @param set2 A character vector of vertex names.
#' @return The numeric lambda parameter.
#' @export
setGeneric(
  "calculate_lambda_between_naive",
  function(object, set1, set2) standardGeneric("calculate_lambda_between_naive")
)

#' @describeIn calculate_lambda_between_naive Method for ECProb objects.
setMethod(
  "calculate_lambda_between_naive",
  "ECProb",
  function(object, set1, set2) {
    # Ensure inputs are unique
    set1 <- unique(set1)
    set2 <- unique(set2)

    valid_set1 <- set1[set1 %in% object@names]
    valid_set2 <- set2[set2 %in% object@names]

    set1_d <- setdiff(valid_set1, valid_set2)
    set2_d <- setdiff(valid_set2, valid_set1)

    if (length(set1_d) == 0 | length(set2_d) == 0){
      return(0)
    }
    degrees_set1 <- unlist(object@degrees[set1_d])
    degrees_set2 <- unlist(object@degrees[set2_d])

    lambda <- 0
    for (i in 1:length(set1_d)) {
      for (j in 1:length(set2_d)) {
        lambda <- lambda +  min(1, (degrees_set1[i] * degrees_set2[j]) / (2 * object@graph_size))
      }
    }
    return(as.numeric(lambda))
  }
)


#' @title Calculate Lambda Between Two Vertex Sets (Optimized Method)
#' @description Calculates the lambda parameter between two vertex sets using a
#' speed-optimized and accurate algorithm. This is the recommended method for accurate calculation.
#' @details This method has a time complexity of O(n*log(n) + m*log(m)).
#' @param object An ECProb object.
#' @param set1 A character vector of vertex names.
#' @param set2 A character vector of vertex names.
#' @return The numeric lambda parameter.
#' @export
setGeneric(
  "calculate_lambda_between",
  function(object, set1, set2) standardGeneric("calculate_lambda_between")
)

#' @describeIn calculate_lambda_between Method for ECProb objects.
setMethod(
  "calculate_lambda_between",
  "ECProb",
  function(object, set1, set2) {
    set1 <- unique(set1)
    set2 <- unique(set2)

    valid_set1 <- set1[set1 %in% object@names]
    valid_set2 <- set2[set2 %in% object@names]

    set1_d <- setdiff(valid_set1, valid_set2)
    set2_d <- setdiff(valid_set2, valid_set1)

    if (length(set1_d) == 0 | length(set2_d) == 0){
      return(0)
    }
    k1 <- sort(unlist(object@degrees[set1_d]), decreasing = TRUE)
    k2 <- sort(unlist(object@degrees[set2_d]), decreasing = TRUE)
    M <- object@graph_size

    find_k1t <- function(i, k1_vec, k2_vec, M_val){
      for (j in 1:length(k2_vec)) {
        if (k1_vec[i] * k2_vec[j] <= 2 * M_val) {
          return(j)
        }
      }
      return(length(k2_vec) + 1)
    }
    vfind_k1t <- Vectorize(find_k1t, vectorize.args = "i")
    k1t <- vfind_k1t(1:length(k1), k1, k2, M)

    k2c <- c(rev(cumsum(rev(k2))), 0)

    lambda <- sum((k1t - 1) + (k1 * k2c[k1t]) / (2 * M))

    return(as.numeric(lambda))
  }
)

#' @title Calculate Lambda Between Two Vertex Sets (Fast Method)
#' @description Calculates an approximation of the lambda parameter between two
#' vertex sets.
#' @details This method is the fastest, with a complexity of O(n+m).
#' @param object An ECProb object.
#' @param set1 A character vector of vertex names.
#' @param set2 A character vector of vertex names.
#' @return The approximate numeric lambda parameter.
#' @export
setGeneric(
  "calculate_lambda_between_fast",
  function(object, set1, set2) standardGeneric("calculate_lambda_between_fast")
)

#' @describeIn calculate_lambda_between_fast Method for ECProb objects.
setMethod(
  "calculate_lambda_between_fast",
  "ECProb",
  function(object, set1, set2) {
    set1 <- unique(set1)
    set2 <- unique(set2)

    valid_set1 <- set1[set1 %in% object@names]
    valid_set2 <- set2[set2 %in% object@names]

    set1_d <- setdiff(valid_set1, valid_set2)
    set2_d <- setdiff(valid_set2, valid_set1)

    if (length(set1_d) == 0 | length(set2_d) == 0){
      return(0)
    }
    k1 <- unlist(object@degrees[set1_d])
    k2 <- unlist(object@degrees[set2_d])
    M <- object@graph_size
    lambda <- sum(k1)*sum(k2)/(2*M)
    return(as.numeric(lambda))
  }
)

#' @title Calculate Lambda Within a Vertex Set (Fast Method)
#' @description Calculates an approximation of the lambda parameter for edges
#' within a single set of vertices.
#' @details This method is the fastest, with a complexity of O(n).
#' @param object An ECProb object.
#' @param set A character vector of vertex names.
#' @return The approximate numeric lambda parameter.
#' @export
setGeneric(
  "calculate_lambda_in_fast",
  function(object, set) standardGeneric("calculate_lambda_in_fast")
)

#' @describeIn calculate_lambda_in_fast Method for ECProb objects.
setMethod(
  "calculate_lambda_in_fast",
  "ECProb",
  function(object, set) {
    valid_set <- unique(set[set %in% object@names])
    if (length(valid_set) <= 1){
      return(0)
    }
    M <- object@graph_size
    k <- unlist(object@degrees[valid_set])
    sum_k <- sum(k)
    lambda <- (sum_k^2 - sum(k^2)) / (4 * M)

    return(as.numeric(lambda))
  }
)

#' @title Calculate Lambda Within a Vertex Set (Naive Method)
#' @description Calculates the lambda parameter for edges within a single set
#' of vertices using a direct, unoptimized summation.
#' @details The naive method has a time complexity of O(n^2).
#' @param object An ECProb object.
#' @param set A character vector of vertex names.
#' @return The numeric lambda parameter.
#' @export
setGeneric(
  "calculate_lambda_in_naive",
  function(object, set) standardGeneric("calculate_lambda_in_naive")
)

#' @describeIn calculate_lambda_in_naive Method for ECProb objects.
setMethod(
  "calculate_lambda_in_naive",
  "ECProb",
  function(object, set) {
    valid_set <- unique(set[set %in% object@names])
    if (length(valid_set) <= 1){
      return(0)
    }
    M <- object@graph_size
    k <- unlist(object@degrees[valid_set])
    lambda <- 0
    m <- length(k)
    if (m > 1) {
      for (i in 1:(m-1)){
        for (j in (i+1):m){
          lambda <- lambda + min(1, k[i]*k[j]/(2*M))
        }
      }
    }
    return(as.numeric(lambda))
  }
)

#' @title Calculate Lambda Within a Vertex Set (Optimized Method)
#' @description Calculates the lambda parameter for edges within a single set
#' of vertices using a speed-optimized and accurate algorithm.
#' @details This method has a time complexity of O(n*log(n)).
#' @param object An ECProb object.
#' @param set A character vector of vertex names.
#' @return The numeric lambda parameter.
#' @export
setGeneric(
  "calculate_lambda_in",
  function(object, set) standardGeneric("calculate_lambda_in")
)

#' @describeIn calculate_lambda_in Method for ECProb objects.
setMethod(
  "calculate_lambda_in",
  "ECProb",
  function(object, set) {

    valid_set <- unique(set[set %in% object@names])
    if (length(valid_set) <= 1){
      return(0)
    }

    find_kt_i <- function(i, k, m, M) {
      for (j in (i + 1):m) {
        if (k[i] * k[j] <= 2 * M) {
          return(j)
        }
      }
      return(m+1)
    }
    vfind_kt_i <- Vectorize(find_kt_i, vectorize.args = "i")

    M <- object@graph_size
    k <- unlist(object@degrees[valid_set])
    m <- length(k)
    k <- sort(k, decreasing = TRUE)

    kc <- c(rev(cumsum(rev(k))), 0)

    if (m <= 1) return(0)

    kt <- vfind_kt_i(c(1:(m-1)), k, m, M)
    kt <- c(kt, m)

    lambda <- sum((kt[1:(m - 1)] - (1:(m - 1)) - 1) + (k[1:(m - 1)] * kc[kt[1:(m - 1)]]) / (2 * M))

    return(as.numeric(lambda))
  }
)


# --------------------------------------------------------------------------- #
# Vectorized Methods
# --------------------------------------------------------------------------- #

#' @title Get Disjoint Sets for Multiple Pairs (Vectorized)
#' @description An internal function that calculates the
#' disjoint element sets for a list of set pairs.
#' @param pairs_dt A `data.table` with two columns ("set1", "set2").
#' @param set_membership_dt A `data.table` with "set_id" and "element" columns.
#' @return A `data.table` with "set1", "set2", "elements1_disjoint" (list),
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

#' @title Get Disjoint Connected Pairs
#'
#' @description Finds all pairs of sets (terms) that have at least one edge
#' connecting their *disjoint* parts in the interaction network (\code{ECProb}).
#' This method uses fast anti-joins on the observed edge list to filter out connections
#' that only exist due to overlapping elements (intersections), returning only pairs
#' with true disjoint connectivity.
#'
#' @param object An \code{ECProb} object (the interaction network).
#' @param set_membership A \code{data.frame} or \code{data.table} representing the
#'   bipartite graph (membership). It must have two columns corresponding to
#'   the set/term ID and the element ID. Columns are assumed to be "term" and "element"
#'   if named, or the first two columns otherwise.
#' @param sets A character vector of set IDs to focus on. If \code{NULL} (default),
#'   all sets in \code{set_membership} are considered.
#'
#' @return A \code{data.table} containing the pairs ("set1", "set2") and their
#'   disjoint \code{observed_edges} count.
#' @export
setGeneric("get_disjoint_connected_pairs", function(object, set_membership, sets = NULL) standardGeneric("get_disjoint_connected_pairs"))

#' @describeIn get_disjoint_connected_pairs Method for ECProb and data.frame/data.table membership.
setMethod("get_disjoint_connected_pairs", signature(object = "ECProb", set_membership = "data.frame"),
          function(object, set_membership, sets = NULL) {

            # 1. Prepare Edge Lists
            network_edges <- data.table(to_dataframe(object))
            setnames(network_edges, c("from", "to"), c("element1", "element2"))

            bipartite_edges <- as.data.table(set_membership)

            # Standardize column names (Term ID, Element ID)
            if (!all(c("term", "element") %in% names(bipartite_edges))) {
              setnames(bipartite_edges, c(1, 2), c("term", "element"))
            }
            bipartite_edges[, `:=`(term = as.character(term), element = as.character(element))]

            # Filter for specific sets if requested
            if (!is.null(sets)) {
              bipartite_edges <- bipartite_edges[term %in% sets]
            }

            if (nrow(bipartite_edges) == 0) {
              return(data.table(set1=character(), set2=character(), observed_edges=integer()))
            }

            # 2. Map element edges to term pairs (Create the 'merged2' table)
            # This table contains every instance of an edge connecting Term1 and Term2
            merged1 <- network_edges[bipartite_edges, on = .(element1 = element), nomatch = 0, allow.cartesian = TRUE]
            setnames(merged1, "term", "term1")

            merged2 <- merged1[bipartite_edges, on = .(element2 = element), nomatch = 0, allow.cartesian = TRUE]
            setnames(merged2, "term", "term2")

            # Filter out self-term loops immediately
            merged2 <- merged2[term1 != term2]

            if (nrow(merged2) == 0) {
              return(data.table(set1=character(), set2=character(), observed_edges=integer()))
            }

            # --- 3. The Anti-Join Optimization: Filter for Disjointness ---

            # To use data.table joins efficiently, we need keys
            setkey(bipartite_edges, element, term)

            # Check 1: Is element1 also in term2? (Bad connection)
            ids_shared1 <- merged2[bipartite_edges, on = .(element1 = element, term2 = term),
                                   nomatch = 0, which = TRUE]

            # Check 2: Is element2 also in term1? (Bad connection)
            ids_shared2 <- merged2[bipartite_edges, on = .(element2 = element, term1 = term),
                                   nomatch = 0, which = TRUE]

            # Combine indices of all "bad" (non-disjoint) edges
            bad_indices <- unique(c(ids_shared1, ids_shared2))

            # 4. Create the clean list of disjoint edges
            if (length(bad_indices) > 0) {
              valid_edges <- merged2[-bad_indices]
            } else {
              valid_edges <- merged2
            }

            # 5. Canonicalize pairs and count observed edges
            final_results <- valid_edges[, .(
              set1 = pmin(term1, term2),
              set2 = pmax(term1, term2)
            )]

            # Count the number of disjoint edges for each pair
            final_results <- final_results[, .(observed_edges = .N), by = .(set1, set2)]

            return(final_results)
          })

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

#' @title Vectorized calculation of fast within-set statistics
#'
#' @description Calculates fast within-set statistics for a collection of sets
#' using a high-performance, vectorized algorithm.
#'
#' @param object An ECProb object.
#' @param sets_dt A data table with a single column ("set_id")
#' @param set_membership_dt A data table with "set_id" and "element" columns.
#'
#' @return A data table with edge-count stattistics (`observed_edge_count`, `lambda`, etc.).
#' @export
#' @examples
#' # Create sample data
#' ecg <- ECGraph(data.frame(p1=c("E1","E3"), p2=c("E2","E4")))
#' ecp <- ECProb(ecg)
#' sets <- data.table(set_id = c("SA", "SB"))
#' memberships <- data.table(
#'   set_id = c("SA","SA","SA","SB","SB"),
#'   element = c("E1","E2","E3","E3","E4")
#' )
#' # Run vectorized calculation
#' calculate_in_stats_fast_vectorized(ecp, sets, memberships)
#'
setGeneric("calculate_in_stats_fast_vectorized",
           function(object, sets_dt, set_membership_dt) standardGeneric("calculate_in_stats_fast_vectorized"))
#' @describeIn calculate_in_stats_fast_vectorized Method for ECProb objects
setMethod("calculate_in_stats_fast_vectorized",
          "ECProb",
          function(object, sets_dt, set_membership_dt) {

            # Ensure input is unique
            set_membership_dt <- unique(set_membership_dt, by = c("set_id", "element"))

            # --- Step 1: Calculate Observed Edge Counts ---
            ecg_edges <- data.table(to_dataframe(object))
            setnames(ecg_edges, c("from", "to"), c("e1", "e2"))
            ecg_edges[, `:=`(canon1 = pmin(e1, e2), canon2 = pmax(e1, e2))]
            ecg_edges <- unique(ecg_edges, by = c("canon1", "canon2"))
            setkey(ecg_edges, canon1, canon2)

            setkey(set_membership_dt, set_id)
            possible_edges <- set_membership_dt[set_membership_dt, on = "set_id", allow.cartesian = TRUE]
            possible_edges <- possible_edges[element < i.element]

            if (nrow(possible_edges) > 0) {
              possible_edges[, `:=`(canon1 = pmin(element, i.element), canon2 = pmax(element, i.element))]
              observed_edges_long <- ecg_edges[possible_edges, on = .(canon1, canon2), nomatch = 0]
              observed_edges_dt <- observed_edges_long[, .(observed_edges = .N), by = set_id]
            } else {
              observed_edges_dt <- data.table(set_id = character(), observed_edges = integer())
            }

            # --- Step 2: Calculate Lambda Components ---
            all_element_degrees_dt <- data.table(element = object@names, degree = unlist(object@degrees))
            setkey(all_element_degrees_dt, element)

            setkey(set_membership_dt, element)
            sets_with_degrees <- all_element_degrees_dt[set_membership_dt, on = "element", nomatch = 0]

            term_summary <- sets_with_degrees[, .(
              sum_of_degrees = sum(degree, na.rm = TRUE),
              sum_of_sq_degrees = sum(degree^2, na.rm = TRUE),
              set_size = .N
            ), by = set_id]

            # --- Step 3: Join all results and perform final calculations ---
            final_dt <- copy(sets_dt)

            final_dt[observed_edges_dt, on = "set_id", observed_edge_count := i.observed_edges]
            final_dt[term_summary, on = "set_id", `:=`(
              sum_of_degrees = i.sum_of_degrees,
              sum_of_sq_degrees = i.sum_of_sq_degrees,
              set_size = i.set_size
            )]

            # Clean NAs from joins for all columns
            cols_to_clean <- c("observed_edge_count", "sum_of_degrees", "sum_of_sq_degrees", "set_size")
            for (col in cols_to_clean) {
              # final_dt[is.na(get(col)), (col) := 0]
              if (col %in% names(final_dt)) {
                if (is.integer(final_dt[[col]])) {
                  final_dt[is.na(get(col)), (col) := 0L]
                } else {
                  final_dt[is.na(get(col)), (col) := 0.0]
                }
              }
            }

            # Calculate final statistics
            final_dt[, lambda := (sum_of_degrees^2 - sum_of_sq_degrees) / (4 * object@graph_size)]
            final_dt[, max_possible_edges := set_size * (set_size - 1) / 2]
            final_dt[, p_value := calculate_p_value(object, observed_edge_count, max_possible_edges, lambda)]
            final_dt[, log2_Anscombe_ratio := 0.5 * (log2(observed_edge_count + 3/8) - log2(lambda + 3/8))]

            return(final_dt[, .(set_id, observed_edge_count, lambda, p_value, log2_Anscombe_ratio, set_size, max_possible_edges)])
          })

# --------------------------------------------------------------------------- #
# General Utility Functions
# --------------------------------------------------------------------------- #

#' @title Calculate Full Edge Count Statistics
#' @description A wrapper that calculates a suite of statistics
#' for an observed edge count.
#' @param object An ECProb object.
#' @param set1 A character vector of vertex names.
#' @param set2 (Optional) A second character vector for between-set analysis.
#' @param observed_edge_count The observed number of edges.
#' @param lambda_method The method for lambda calculation ("accurate", "optimized", "fast").
#' @return A list containing the p-value, lambda, log2_Anscombe_ratio and log2_relative_change.
#' @export
setGeneric("edge_count_statistics",
           function(object, set1, set2 = NULL, observed_edge_count, lambda_method = "optimized")
             standardGeneric("edge_count_statistics"))

#' @describeIn edge_count_statistics Method for ECProb objects.
setMethod(
  "edge_count_statistics",
  "ECProb",
  function(object, set1, set2 = NULL, observed_edge_count, lambda_method = "optimized") {

    valid_set1 <- unique(set1[set1 %in% object@names])
    max_possible_edges <- 0
    size1 <- 0
    size2 <- 0

    if (is.null(set2)) {
      lambda <- switch(lambda_method,
                       accurate = calculate_lambda_in_naive(object, valid_set1),
                       optimized = calculate_lambda_in(object, valid_set1),
                       fast = calculate_lambda_in_fast(object, valid_set1),
                       stop("Invalid method for within-set analysis."))
      max_possible_edges <- length(valid_set1) * (length(valid_set1) - 1) / 2
      size1 <- length(valid_set1)
    } else {
      valid_set2 <- unique(set2[set2 %in% object@names])
      lambda <- switch(lambda_method,
                       accurate = calculate_lambda_between_naive(object, valid_set1, valid_set2),
                       optimized = calculate_lambda_between(object, valid_set1, valid_set2),
                       fast = calculate_lambda_between_fast(object, valid_set1, valid_set2),
                       stop("Invalid method for between-set analysis."))
      max_possible_edges <- length(valid_set1) * length(valid_set2)
      size1 <- length(valid_set1)
      size2 <- length(valid_set2)
    }

    p_value <- calculate_p_value(object, observed_edge_count, max_possible_edges, lambda)

    log2_Anscombe_ratio <- NA_real_
    if (!is.na(lambda) && (lambda + 3/8) > 0 && (observed_edge_count + 3/8) > 0) {
      log2_Anscombe_ratio <- 0.5 * (log2(observed_edge_count + 3/8) - log2(lambda + 3/8))
    }

    return(list(p_value = p_value,
                observed_edge_count = observed_edge_count,
                log2_Anscombe_ratio = log2_Anscombe_ratio,
                lambda = lambda,
                max_possible_edges = max_possible_edges,
                size1 = size1,
                size2 = size2))
  })

#' @title Calculate P-value for an Observed Edge Count
#'
#' @description Calculates the p-value for an observed edge count based on a
#' truncated Poisson distribution.
#'
#' @param object An ECProb object (for S4 dispatch).
#' @param z The observed number of edges.
#' @param m The maximum possible number of edges.
#' @param lambda The Poisson parameter (expected edge count).
#'
#' @return The calculated p-value.
#' @export
setGeneric(
  "calculate_p_value",
  function(object, z, m, lambda)
    standardGeneric("calculate_p_value")
)

#' @describeIn calculate_p_value Method for ECProb objects.
setMethod(
  "calculate_p_value",
  "ECProb",
  function(object, z, m, lambda) {
    # This version is fully vectorized to handle entire columns of data at once.
    p_values <- ifelse(is.na(lambda) | lambda < 0, NA_real_, {
      alpha <- stats::ppois(m, lambda, lower.tail = TRUE)
      ifelse(alpha == 0, 1.0, {
        (stats::ppois(z - 1, lambda, lower.tail = FALSE) - stats::ppois(m, lambda, lower.tail = FALSE)) / alpha
      })
    })
    return(as.numeric(p_values))
  }
)


#' @title Summarize Graph's Suitability for Fast Lambda Approximation
#'
#' @description Provides statistics on the pairwise Bernoulli parameters (p_ij)
#' to help assess if fast lambda approximation methods are suitable for the graph.
#' This method also reports the proportion of vertices that are problematic for fast
#' lambda approximation in the sense that their degree is greater than sqrt(2*M), M
#' being the graph size.
#'
#' @param object An ECProb object.
#'
#' @return A list of summary statistics: 'p_ij_over_1' = proportion of vertex pairs such
#' that p_ij > 1 and 'prop_problematic_vertices' = proportion of vertices with degree > sqrt(2*M)
#' @export
setGeneric("summarize_suitability_fast",
           function(object) standardGeneric("summarize_suitability_fast"))

#' @describeIn summarize_suitability_fast Method for ECProb objects
setMethod("summarize_suitability_fast",
          "ECProb",
          function(object) {

            # proportion of problematic pairs
            degrees <- unlist(object@degrees)
            N <- length(degrees)
            if (N == 0) {
              return(list(
                pij_over_1 = 0,
                prop_problematic_vertices = 0,
                summary_pij = summary(numeric(0))
              ))
            }

            M <- object@graph_size
            two_M <- 2 * M

            degree_distribution <- table(degrees)
            k <- as.numeric(names(degree_distribution))
            p <- as.numeric(degree_distribution) / N
            prod_matrix <- outer(k, k, "*")/(2*M)
            p_matrix <- outer(p, p)
            diag(prod_matrix) <- 0
            diag(p_matrix) <- 0
            prod_vec <- as.vector(prod_matrix)
            p_vec <- as.vector(p_matrix)
            total_p <- sum(p_vec)
            prop_bad <- sum(p_vec[prod_vec >= 1])
            pij_over_1 <- if(total_p == 0) 0 else (prop_bad / total_p)

            # proportion of problematic vertices
            threshold <- sqrt(2 * M)
            problematic_vertices <- degrees[degrees > threshold]
            prop_problematic_vertices <- length(problematic_vertices) / N

            # -- proportions of unique distinct degree pairs (undirected) ---
            l_k <- length(k)
            n_k <- as.numeric(degree_distribution)
            count_matrix <- outer(n_k, n_k)
            diag(count_matrix) <- n_k * (n_k - 1) / 2
            undirected_matrix <- matrix(0.5, l_k, l_k) + diag(0.5, l_k)
            undirected_vec <- as.vector(undirected_matrix)
            count_vec <- as.vector(count_matrix)
            count_vec <- count_vec*undirected_vec
            prop_vec <- count_vec / sum(count_vec)

            # --- pij distribution for undirected edges ---
            pij_matrix <- outer(k, k) / two_M
            pij_vec <- as.vector(pij_matrix)
            dt <- data.table(pij = pij_vec,
                             count = count_vec,
                             prop = prop_vec)

            dt <- dt[, .(
              count = as.integer(sum(count)),
              prop = sum(prop)
            ), by = .(pij)]

            pij_over_1 = sum(dt[pij >= 1, prop])

            # --- proportion of problematic vertices ---
            threshold <- sqrt(two_M)
            problematic_vertices <- degrees[degrees > threshold]
            prop_problematic_vertices <- length(problematic_vertices) / N

            # -- results including summary of pij distribution
            return(list(
              pij_over_1 = pij_over_1,
              prop_problematic_vertices = prop_problematic_vertices,
              summary_pij = summary(rep(dt$pij, dt$count)),
              pij_distribution = dt[count > 0]
            ))
          })
