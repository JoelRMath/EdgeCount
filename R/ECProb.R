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
#' Lambda calculation methods:
#' \itemize{
#'   \item Lambda between two sets
#'   \itemize{
#'     \item \code{\link{calculate_lambda_between}} (optimized and accurate)
#'     \item \code{\link{calculate_lambda_between_naive}} (slow, for validation only)
#'     \item \code{\link{calculate_lambda_between_fast}} (fast approximation)
#'    }
#'    \item Lamda within a set
#'    \itemize{
#'     \item \code{\link{calculate_lambda_in}} (optimized and accurate)
#'     \item \code{\link{calculate_lambda_in_naive}} (slow, for validation only)
#'     \item \code{\link{calculate_lambda_in_fast}} (fast approximation)
#'    }
#' }
#'
#' Probability and statistics functions:
#' \itemize{
#'   \item \code{\link{calculate_p_value}}
#'   \item \code{\link{edge_count_statistics}}
#'   \item \code{\link{summarize_suitability_fast}}
#' }
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
  graph_order <- length(unlist(ecgraph@degrees))
  average_degree <- 2*graph_size/graph_order
  
  new("ECProb", ecgraph, graph_size = graph_size, graph_order = graph_order,
      average_degree = average_degree)
}

#' @title Calculate Lambda Between Two Vertex Sets (Naive Method)
#'
#' @description Calculates the lambda parameter (expected edge count) between two
#' vertex sets using a direct, unoptimized summation.
#' @details The naive method has a time complexity of O(n*m), where n and m are
#' the sizes of the two sets. It is accurate but can be very slow for large sets.
#' Its primary purpose is for validation and speed-testing against optimized methods.
#'
#' @param object An ECProb object.
#' @param set1 A character vector of vertex names.
#' @param set2 A character vector of vertex names.
#'
#' @return The numeric lambda parameter.
#' @export
#' @examples
#' ecg <- ECGraph(data.frame(p1=c("A","B","C"), p2=c("X","Y","Z")))
#' ecp <- ECProb(ecg)
#' calculate_lambda_between_naive(ecp, c("A", "B"), c("X", "Y"))
#'
setGeneric(
  "calculate_lambda_between_naive",
  function(object, set1, set2) standardGeneric("calculate_lambda_between_naive")
)

#' @describeIn calculate_lambda_between_naive Method for ECProb objects.
setMethod(
  "calculate_lambda_between_naive",
  "ECProb",
  function(object, set1, set2) {
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
#'
#' @description Calculates the lambda parameter between two vertex sets using a
#' speed-optimized and accurate algorithm. This is the recommended method for accurate calculation.
#' @details This method has a time complexity of O(n*log(n) + m*log(m)) due to the
#' sorting of degree sequences, which is significantly faster than the O(n*m)
#' complexity of the naive method for large sets. It produces the exact same
#' result as the naive method.
#'
#' @param object An ECProb object.
#' @param set1 A character vector of vertex names.
#' @param set2 A character vector of vertex names.
#'
#' @return The numeric lambda parameter.
#' @export
#' @examples
#' ecg <- ECGraph(data.frame(p1=c("A","B","C"), p2=c("X","Y","Z")))
#' ecp <- ECProb(ecg)
#' calculate_lambda_between(ecp, c("A", "B"), c("X", "Y"))
#'
setGeneric(
  "calculate_lambda_between",
  function(object, set1, set2) standardGeneric("calculate_lambda_between")
)

#' @describeIn calculate_lambda_between Method for ECProb objects.
setMethod(
  "calculate_lambda_between",
  "ECProb",
  function(object, set1, set2) {
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
#'
#' @description Calculates an approximation of the lambda parameter between two
#' vertex sets.
#' @details This method is the fastest, with a complexity of O(n+m). It is an
#' approximation based on the assumption that the Bernoulli probability of an edge
#' between any two vertices i and j, p_ij = (k_i * k_j) / (2M), is always less than 1.
#' The accuracy of this method decreases as the proportion of vertex pairs for which
#' p_ij > 1 (let's call this proportion alpha) increases. You can estimate alpha
#' for your graph using \code{\link{summarize_suitability_fast}}.
#'
#' @param object An ECProb object.
#' @param set1 A character vector of vertex names.
#' @param set2 A character vector of vertex names.
#'
#' @return The approximate numeric lambda parameter.
#' @export
#' @examples
#' ecg <- ECGraph(data.frame(p1=c("A","B","C"), p2=c("X","Y","Z")))
#' ecp <- ECProb(ecg)
#' calculate_lambda_between_fast(ecp, c("A", "B"), c("X", "Y"))
#'
setGeneric(
  "calculate_lambda_between_fast",
  function(object, set1, set2) standardGeneric("calculate_lambda_between_fast")
)

#' @describeIn calculate_lambda_between_fast Method for ECProb objects.
setMethod(
  "calculate_lambda_between_fast",
  "ECProb",
  function(object, set1, set2) {
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
#'
#' @description Calculates an approximation of the lambda parameter for edges
#' within a single set of vertices.
#' @details This method is the fastest, with a complexity of O(n). It is an
#' approximation based on the assumption that the Bernoulli probability of an edge
#' between any two vertices i and j, p_ij = (k_i * k_j) / (2M), is always less than 1.
#' The accuracy of this method decreases as the proportion of vertex pairs for which
#' p_ij > 1 increases. You can estimate this proportion for your graph using
#' \code{\link{summarize_suitability_fast}}.
#'
#' @param object An ECProb object.
#' @param set A character vector of vertex names.
#'
#' @return The approximate numeric lambda parameter.
#' @export
#' @examples
#' ecg <- ECGraph(data.frame(p1=c("A","B","C"), p2=c("B","C","A")))
#' ecp <- ECProb(ecg)
#' calculate_lambda_in_fast(ecp, c("A", "B"))
#'
setGeneric(
  "calculate_lambda_in_fast",
  function(object, set) standardGeneric("calculate_lambda_in_fast")
)

#' @describeIn calculate_lambda_in_fast Method for ECProb objects.
setMethod(
  "calculate_lambda_in_fast",
  "ECProb",
  function(object, set) {
    valid_set <- set[set %in% object@names]
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
#'
#' @description Calculates the lambda parameter for edges within a single set
#' of vertices using a direct, unoptimized summation.
#' @details The naive method has a time complexity of O(n^2), where n is the
#' size of the set. It is accurate but can be very slow for large sets.
#' Its primary purpose is for validation and speed-testing purposes.
#'
#' @param object An ECProb object.
#' @param set A character vector of vertex names.
#'
#' @return The numeric lambda parameter.
#' @export
#' @examples
#' ecg <- ECGraph(data.frame(p1=c("A","B","C"), p2=c("B","C","A")))
#' ecp <- ECProb(ecg)
#' calculate_lambda_in_naive(ecp, c("A", "B"))
#'
setGeneric(
  "calculate_lambda_in_naive",
  function(object, set) standardGeneric("calculate_lambda_in_naive")
)

#' @describeIn calculate_lambda_in_naive Method for ECProb objects.
setMethod(
  "calculate_lambda_in_naive",
  "ECProb",
  function(object, set) {
    valid_set <- set[set %in% object@names]
    if (length(valid_set) <= 1){
      return(0)
    }
    M <- object@graph_size
    k <- unlist(object@degrees[valid_set])
    lambda <- 0
    m <- length(k)
    for (i in 1:(m-1)){
      for (j in (i+1):m){
        lambda <- lambda + min(1, k[i]*k[j]/(2*M))
      }
    }
    return(as.numeric(lambda))
  }
)

#' @title Calculate Lambda Within a Vertex Set (Optimized Method)
#'
#' @description Calculates the lambda parameter for edges within a single set
#' of vertices using a speed-optimized and accurate algorithm.
#' @details This method has a time complexity of O(n*log(n)) due to the sorting
#' of the degree sequence, which is significantly faster than the O(n^2)
#' complexity of the naive method for large sets.
#'
#' @param object An ECProb object.
#' @param set A character vector of vertex names.
#'
#' @return The numeric lambda parameter.
#' @export
#' @examples
#' ecg <- ECGraph(data.frame(p1=c("A","B","C"), p2=c("B","C","A")))
#' ecp <- ECProb(ecg)
#' calculate_lambda_in(ecp, c("A", "B"))
#'
setGeneric(
  "calculate_lambda_in",
  function(object, set) standardGeneric("calculate_lambda_in")
)

#' @describeIn calculate_lambda_in Method for ECProb objects.
setMethod(
  "calculate_lambda_in",
  "ECProb",
  function(object, set) {
    
    find_kt_i <- function(i, k, m, M) {
      for (j in (i + 1):m) {
        if (k[i] * k[j] <= 2 * M) {
          return(j)
        }
      }
      return(m+1)
    }
    vfind_kt_i <- Vectorize(find_kt_i, vectorize.args = "i")
    
    valid_set <- set[set %in% object@names]
    if (length(valid_set) <= 1){
      return(0)
    }
    M <- object@graph_size
    k <- unlist(object@degrees[valid_set])
    m <- length(k)
    k <- sort(k, decreasing = TRUE)
    
    kc <- c(rev(cumsum(rev(k))), 0)
    
    kt <- vfind_kt_i(c(1:(m-1)), k, m, M)
    kt <- c(kt, m)
    
    lambda <- sum((kt[1:(m - 1)] - (1:(m - 1)) - 1) + (k[1:(m - 1)] * kc[kt[1:(m - 1)]]) / (2 * M))
    
    return(as.numeric(lambda))
  }
)

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
#' @examples
#' ecg <- ECGraph(data.frame(p1="A", p2="B")) # Dummy graph
#' ecp <- ECProb(ecg)
#' # P-value for observing 5 edges when max is 10 and lambda is 2
#' calculate_p_value(ecp, z = 5, m = 10, lambda = 2)
#'
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
    if(is.na(lambda) || lambda < 0) return(NA_real_)
    alpha <- stats::ppois(m, lambda, lower.tail = TRUE)
    if (alpha == 0) return(1) # Avoid division by zero if m is very small compared to lambda
    p <- (stats::ppois(z-1, lambda, lower.tail = FALSE) - stats::ppois(m, lambda, lower.tail = FALSE)) / alpha
    return(as.numeric(p))
  }
)

#' @title Calculate Full Edge Count Statistics
#'
#' @description A wrapper that calculates a suite of statistics (p-value, z-score, etc.)
#' for an observed edge count. It automatically determines the correct lambda and
#' maximum edge count based on the input parameters.
#'
#' @param object An ECProb object.
#' @param set1 A character vector of vertex names.
#' @param set2 (Optional) A second character vector for between-set analysis.
#' @param observed_edge_count The observed number of edges.
#' @param lambda_method The method for lambda calculation ("accurate", "optimized", "fast").
#'
#' @return A list containing the p-value, z-score, and other statistics.
#' @export
#' @examples
#' ecg <- ECGraph(data.frame(p1=c("A","B","C"), p2=c("B","C","A")))
#' ecp <- ECProb(ecg)
#' # Calculate stats for observing 3 edges within the set {A, B, C}
#' edge_count_statistics(ecp, set1 = c("A", "B", "C"), observed_edge_count = 3)
#'
setGeneric("edge_count_statistics",
           function(object, set1, set2 = NULL, observed_edge_count, lambda_method = "optimized")
             standardGeneric("edge_count_statistics"))

#' @describeIn edge_count_statistics Method for ECProb objects.
setMethod(
  "edge_count_statistics",
  "ECProb",
  function(object, set1, set2 = NULL, observed_edge_count, lambda_method = "optimized") {
    
    valid_set1 <- set1[set1 %in% object@names]
    
    if (is.null(set2)) {
      lambda <- switch(lambda_method,
                       accurate = calculate_lambda_in_naive(object, valid_set1),
                       optimized = calculate_lambda_in(object, valid_set1),
                       fast = calculate_lambda_in_fast(object, valid_set1),
                       stop("Invalid method for within-set analysis."))
      max_possible_edges <- length(valid_set1) * (length(valid_set1) - 1) / 2
    } else {
      valid_set2 <- set2[set2 %in% object@names]
      lambda <- switch(lambda_method,
                       accurate = calculate_lambda_between_naive(object, valid_set1, valid_set2),
                       optimized = calculate_lambda_between(object, valid_set1, valid_set2),
                       fast = calculate_lambda_between_fast(object, valid_set1, valid_set2),
                       stop("Invalid method for between-set analysis."))
      max_possible_edges <- length(valid_set1) * length(valid_set2)
    }
    
    p_value <- calculate_p_value(object, observed_edge_count, max_possible_edges, lambda)
    
    z_score <- NA_real_
    if (!is.na(lambda) && lambda > 0){
      z_score <- (observed_edge_count - lambda) / sqrt(lambda)
    }
    
    log2_Anscombe_ratio <- NA_real_
    if (!is.na(lambda) && (lambda + 3/8) > 0 && (observed_edge_count + 3/8) > 0) {
      log2_Anscombe_ratio <- 0.5 * (log2(observed_edge_count + 3/8) - log2(lambda + 3/8))
    }
    
    log2_relative_change <- NA_real_
    if (!is.na(lambda) && lambda > 0 && observed_edge_count > 0){
      log2_relative_change <- log2(observed_edge_count) - log2(lambda)
    }
    
    return(list(p_value = p_value,
                z_score = z_score,
                log2_Anscombe_ratio = log2_Anscombe_ratio,
                log2_relative_change = log2_relative_change))
  })

#' @title Summarize Graph's Suitability for Fast Lambda Approximation
#'
#' @description Provides statistics on the pairwise Bernoulli parameters (p_ij)
#' to help assess if fast lambda approximation methods are suitable for the graph.
#'
#' @param object An ECProb object.
#'
#' @return A list of summary statistics for the p_ij distribution.
#' @export
#' @examples
#' # A dense graph might be less suitable for fast approximation
#' ecg_dense <- ECGraph(as.data.frame(as.matrix(igraph::as_adjacency_matrix(
#'   igraph::make_full_graph(10)
#' ))))
#' ecp_dense <- ECProb(ecg_dense)
#' summarize_suitability_fast(ecp_dense)
#'
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
