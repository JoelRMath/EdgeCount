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
#' @description An internal, high-performance function that calculates the
#' disjoint element sets for a list of set pairs.
#' @param pairs_dt A `data.table` with two columns ("set1", "set2").
#' @param set_membership_dt A `data.table` with "set_id" and "element" columns.
#' @return A `data.table` with "set1", "set2", "elements1_disjoint" (list),
#'   and "elements2_disjoint" (list).
#' @keywords internal
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

#' @title Vectorized calculation of fast between-set statistics
#'
#' @description Calculates fast between-set statistics for a collection of set pairs
#' using a high-performance, vectorized algorithm.
#'
#' @param object An ECProb object.
#' @param pairs_dt A data table of pairs of element sets.
#' @param set_membership_dt A data table mapping set IDs to their elements.
#'
#' @return A final `data.table` with stats (`observed_edges`, `lambda`, etc.).
#' @export
#' @examples
#' # Create sample data
#' ecg <- ECGraph(data.frame(p1=c("E1","E3"), p2=c("E4","E5")))
#' ecp <- ECProb(ecg)
#' pairs <- data.table(set1 = c("SA", "SC"), set2 = c("SB", "SD"))
#' memberships <- data.table(
#'   set_id = c("SA","SA","SB","SC","SC","SD"),
#'   element = c("E1","E2","E2","E3","E4","E5")
#' )
#' # Run vectorized calculation
#' calculate_between_stats_fast_vectorized(ecp, pairs, memberships)
#'
setGeneric("calculate_between_stats_fast_vectorized",
           function(object, pairs_dt, set_membership_dt) standardGeneric("calculate_between_stats_fast_vectorized"))

#' @describeIn calculate_between_stats_fast_vectorized Method for ECProb objects
setMethod("calculate_between_stats_fast_vectorized",
          "ECProb",
          function(object, pairs_dt, set_membership_dt) {

            disjoint_sets_dt <- get_disjoint_sets(pairs_dt, set_membership_dt)
            disjoint_sets_dt[, pair_id := .I]

            all_element_degrees_dt <- data.table(
              element = object@names,
              degree = unlist(object@degrees)
            )
            data.table::setkey(all_element_degrees_dt, element)

            long_disjoint1 <- disjoint_sets_dt[, .(element = unlist(elements1_disjoint)), by = pair_id]
            data.table::setkey(long_disjoint1, element)
            sums1_dt <- all_element_degrees_dt[long_disjoint1, on = "element", nomatch = 0][,
                                                                                            .(sum_degrees1 = sum(degree, na.rm = TRUE)), by = pair_id]

            long_disjoint2 <- disjoint_sets_dt[, .(element = unlist(elements2_disjoint)), by = pair_id]
            data.table::setkey(long_disjoint2, element)
            sums2_dt <- all_element_degrees_dt[long_disjoint2, on = "element", nomatch = 0][,
                                                                                            .(sum_degrees2 = sum(degree, na.rm = TRUE)), by = pair_id]

            ecg_edges <- data.table(to_dataframe(object))
            setnames(ecg_edges, c("from", "to"), c("e1", "e2"))
            ecg_edges[, `:=`(canon1 = pmin(e1, e2), canon2 = pmax(e1, e2))]
            ecg_edges <- unique(ecg_edges, by = c("canon1", "canon2"))
            setkey(ecg_edges, canon1, canon2)

            setkey(long_disjoint1, pair_id)
            setkey(long_disjoint2, pair_id)
            possible_edges <- long_disjoint1[long_disjoint2, on = "pair_id", allow.cartesian = TRUE, nomatch=0]
            if (nrow(possible_edges) > 0) {
              setnames(possible_edges, c("element", "i.element"), c("element1", "element2"))
              possible_edges[, `:=`(canon1 = pmin(element1, element2), canon2 = pmax(element1, element2))]
              observed_edges_long <- ecg_edges[possible_edges, on = .(canon1, canon2), nomatch = 0]
              observed_edges_dt <- observed_edges_long[, .(observed_edges = .N), by = pair_id]
            } else {
              observed_edges_dt <- data.table(pair_id = integer(), observed_edges = integer())
            }

            all_pairs_ids <- data.table(pair_id = 1:nrow(disjoint_sets_dt))
            setkey(all_pairs_ids, pair_id)
            setkey(observed_edges_dt, pair_id)
            observed_edges_dt <- observed_edges_dt[all_pairs_ids]
            observed_edges_dt[is.na(observed_edges), observed_edges := 0L]

            final_dt <- copy(disjoint_sets_dt)

            final_dt[sums1_dt, on = "pair_id", sum_degrees1 := i.sum_degrees1]
            final_dt[sums2_dt, on = "pair_id", sum_degrees2 := i.sum_degrees2]

            final_dt[is.na(sum_degrees1), sum_degrees1 := 0]
            final_dt[is.na(sum_degrees2), sum_degrees2 := 0]

            final_dt[observed_edges_dt, on = "pair_id", observed_edges := i.observed_edges]

            final_dt[, `:=`(
              size1 = lengths(elements1_disjoint),
              size2 = lengths(elements2_disjoint)
            )]
            final_dt[, max_possible_edges := as.numeric(size1) * as.numeric(size2)]
            final_dt[, lambda := (sum_degrees1 * sum_degrees2) / (2 * object@graph_size)]
            final_dt[, p_value := calculate_p_value(object, observed_edges, max_possible_edges, lambda)]
            final_dt[, log2_Anscombe_ratio := 0.5 * (log2(observed_edges + 3/8) - log2(lambda + 3/8))]

            return(final_dt[, .(set1, set2, observed_edges, lambda, p_value, log2_Anscombe_ratio)])
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
#' @return A final `data.table` with stats (`observed_edges`, `lambda`, etc.).
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

            set_membership_dt <- unique(set_membership_dt, by = c("set_id", "element"))

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

            all_element_degrees_dt <- data.table(
              element = object@names,
              degree = unlist(object@degrees)
            )
            setkey(all_element_degrees_dt, element)

            setkey(set_membership_dt, element)
            sets_with_degrees <- all_element_degrees_dt[set_membership_dt, on = "element", nomatch = 0]

            term_summary <- sets_with_degrees[, .(
              sum_of_degrees = sum(degree, na.rm = TRUE),
              sum_of_sq_degrees = sum(degree^2, na.rm = TRUE),
              set_size = .N
            ), by = set_id]

            final_dt <- copy(sets_dt)

            final_dt[observed_edges_dt, on = "set_id", observed_edges := i.observed_edges]
            final_dt[term_summary, on = "set_id", `:=`(
              sum_of_degrees = i.sum_of_degrees,
              sum_of_sq_degrees = i.sum_of_sq_degrees,
              set_size = i.set_size
            )]

            final_dt[is.na(observed_edges), observed_edges := 0L]
            final_dt[is.na(sum_of_degrees), sum_of_degrees := 0]
            final_dt[is.na(sum_of_sq_degrees), sum_of_sq_degrees := 0]
            final_dt[is.na(set_size), set_size := 0]

            final_dt[, lambda := (sum_of_degrees^2 - sum_of_sq_degrees) / (4 * object@graph_size)]
            final_dt[, max_possible_edges := set_size * (set_size - 1) / 2]
            final_dt[, p_value := calculate_p_value(object, observed_edges, max_possible_edges, lambda)]
            final_dt[, log2_Anscombe_ratio := 0.5 * (log2(observed_edges + 3/8) - log2(lambda + 3/8))]

            return(final_dt[, .(set_id, observed_edges, lambda, p_value, log2_Anscombe_ratio)])
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

    if (is.null(set2)) {
      lambda <- switch(lambda_method,
                       accurate = calculate_lambda_in_naive(object, valid_set1),
                       optimized = calculate_lambda_in(object, valid_set1),
                       fast = calculate_lambda_in_fast(object, valid_set1),
                       stop("Invalid method for within-set analysis."))
      max_possible_edges <- length(valid_set1) * (length(valid_set1) - 1) / 2
    } else {
      valid_set2 <- unique(set2[set2 %in% object@names])
      lambda <- switch(lambda_method,
                       accurate = calculate_lambda_between_naive(object, valid_set1, valid_set2),
                       optimized = calculate_lambda_between(object, valid_set1, valid_set2),
                       fast = calculate_lambda_between_fast(object, valid_set1, valid_set2),
                       stop("Invalid method for between-set analysis."))
      max_possible_edges <- length(valid_set1) * length(valid_set2)
    }

    p_value <- calculate_p_value(object, observed_edge_count, max_possible_edges, lambda)

    log2_Anscombe_ratio <- NA_real_
    if (!is.na(lambda) && (lambda + 3/8) > 0 && (observed_edge_count + 3/8) > 0) {
      log2_Anscombe_ratio <- 0.5 * (log2(observed_edge_count + 3/8) - log2(lambda + 3/8))
    }

    log2_relative_change <- NA_real_
    if (!is.na(lambda) && lambda > 0 && observed_edge_count > 0){
      log2_relative_change <- log2(observed_edge_count) - log2(lambda)
    }

    return(list(p_value = p_value,
                observed_edge_count = observed_edge_count,
                log2_Anscombe_ratio = log2_Anscombe_ratio,
                log2_relative_change = log2_relative_change,
                lambda = lambda))
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
#'
#' @param object An ECProb object.
#'
#' @return A list of summary statistics for the p_ij distribution.
#' @export
setGeneric("summarize_suitability_fast",
           function(object) standardGeneric("summarize_suitability_fast"))

#' @describeIn summarize_suitability_fast Method for ECProb objects
setMethod("summarize_suitability_fast",
          "ECProb",
          function(object) {
            degrees <- unlist(object@degrees)
            N <- length(degrees)
            M <- object@graph_size
            degree_distribution <- table(degrees)

            k <- as.numeric(names(degree_distribution))
            p <- as.numeric(degree_distribution) / N
            prod <- as.vector(outer(k, k, "*"))/(2*M)
            q <- as.vector(outer(p, p))
            capped <- pmin(1, prod)
            df <- data.frame(prod = prod, q = q, capped = capped)
            prop <- sum(df$q[df$capped == 1])

            return(list(
              pij_over_1 = prop,
              summary_pij = summary(prod),
              summary_capped_pij = summary(capped)
            ))
          })

