#' EdgeCount: A package for Edge-Count based analyses
#'
#' The EdgeCount package provides S4 classes and methods for performing
#' two types of analyses: - 1) connectedness within or between sets of vertices based on
#' a null model of random graph with given expected degrees (RGGED) - 2) scoring of ranked
#' lists of elements based on a bipartite graph of element-term memberships and the RGGED.
#'
#' @section Key Functions and Classes:
#' The main components of the package are:
#' \itemize{
#'   \item \code{\link{ECGraph}}: The base class for representing graphs.
#'   \item \code{\link{ECProb}}: Extends ECGraph with probability scoring methods based
#'   on a model of Random Graph with Given Expected Degrees (RGGED).
#'   \item \code{\link{ECTermScoring}}: The main class for term-enrichment analyses based
#'   on term-element bipartite graphs and probabilistic methods (RGGED).
#'   \itemize{
#'   \item \code{\link{terms_ecset_statistics}}
#'    \item \code{\link{terms_ecranks_statistics}}
#'    \item \code{\link{run_vsea}}
#'   }
#' }
#'@references
#' 1. Chung & Lu (2002). The average distances in random graphs with given
#' expected degrees. PNAS, 99 (25) 15879-15882.
#' \url{https://www.pnas.org/doi/10.1073/pnas.252631999}
#'
#' 2. Pradines et al. (2005). Analyzing Protein Lists with Large Networks:
#' Edge-Count Probabilities in Random Graphs with Given Expected Degrees.
#' J. Comp. Biol. 12(2):113-28.
#' \url{https://www.liebertpub.com/doi/10.1089/cmb.2005.12.113}
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/joelRMath/EdgeCount}
#'   \item Report bugs at \url{https://github.com/joelRMath/EdgeCount/issues}
#' }
#'
#' @keywords internal
"_PACKAGE"

#' @import methods
#' @importFrom stats median na.omit setNames
#' @importFrom utils head stack tail read.table
#' @importFrom parallel detectCores mclapply
utils::globalVariables(c(
  ".", "term_id", "set_id", "element", "term", "set1", "set2", "pair_id",
  "elements1_disjoint", "elements2_disjoint", "sum_degrees1", "sum_degrees2",
  "observed_edges", "lambda", "p_value", "log2_Anscombe_ratio", "max_possible_edges",
  "size1", "size2", "e1", "e2", "canon1", "canon2", "sum_of_degrees",
  "sum_of_sq_degrees", "set_size", "observed_edge_count", "term1", "term2",
  "t_idx", "e_idx", "perm_id", "rank", "e_deg", "univ_cumsum", "t_deg",
  "observed_ec", "max_ec", "score", "min_score", "max_score", "median_score",
  "rank_at_max", "rank_at_min", "mean_pos", "mean_neg", "nes", "null_nes",
  "fdr_q_value", "E_FalsePos", "Rank", "FDR_raw", "tm", "sc", "rk", "ec",
  "degree", "cumsum_degrees", "global_rank", "rank_in_term", "..final_cols",
  "i.term_degree", "i.term_size", "i.sum_degrees_set", "i.set_size", "i.sum_of_sq_degrees",
  "i.sum_of_degrees", "i.observed_edges", "i.mean_pos", "i.mean_neg", "i.t_deg_calc",
  "input_set_id", "term_degree", "sum_degrees_set", "mean_null_pos", "mean_null_neg", "N",
  "variable", "term_size", "v", "u", "i.element", "count", "prop", "pij", "element_id"
))
NULL
