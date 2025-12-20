#' @include ECProb.R
#' @title ECTermScoring S4 Class
#'
#' @description The main class for EdgeCount analysis with bipartite graphs, representing a bipartite
#' network of term-element memberships. It extends the ECProb class.
#'
#' @slot ecprob An object of class ECProb, representing the underlying
#'    undirected version of the bipartite graph.
#' @slot elements A character vector of all unique vertex IDs identified as elements.
#' @slot terms A character vector of all unique vertex IDs identified as terms.
#'
#' @name ECTermScoring-class
#' @rdname ECTermScoring-class
#' @aliases ECTermScoring-class
#' @exportClass ECTermScoring
#' @seealso \code{\link{ECTermScoring}} for the constructor.
setClass("ECTermScoring",
         slots = list(
           ecprob    = "ECProb",
           elements  = "character",
           terms     = "character"
         )
)

#' @title Constructor for ECTermScoring Objects
#'
#' @description Creates an \code{\linkS4class{ECTermScoring}} object from a data frame or file.
#'
#' @param term_element_edges A data frame or path to a file. Needs two columns:
#'    one for term IDs and one for element IDs, representing membership.
#'    If a file path, it's assumed to be tab-separated with a header by default.
#' @param col_term The name or index of the column containing term IDs
#'    (in `term_element_edges`). Defaults to 1.
#' @param col_element The name or index of the column containing element IDs
#'    (in `term_element_edges`). Defaults to 2.
#' @param ... Additional arguments passed to `utils::read.table` if
#'    `term_element_edges` is a file path (e.g., `sep`, `header`).
#'
#' @return An object of class \code{\linkS4class{ECTermScoring}}.
#' @seealso
#' Primary analysis functions:
#' \itemize{
#'    \item \code{\link{terms_ecset_statistics}}: Score terms against a single set of elements.
#'    \item \code{\link{terms_ecset_statistics_fdr}}: Wrapper for single-set scoring with Empirical FDR correction.
#'    \item \code{\link{run_vsea}}: Perform Vertex Set Enrichment Analysis (VSEA) on a ranked list of elements.
#'    \item \code{\link{terms_ecranks_statistics}}: Calculate raw enrichment profiles for ranked lists.
#' }
#'
#' FDR and Vectorized functions:
#' \itemize{
#'    \item \code{\link{calculate_ecset_fdr}}: Calculate Empirical FDR for single-set enrichment results.
#'    \item \code{\link{terms_ecset_statistics_vectorized}}: Efficiently score terms against multiple element sets.
#' }
#'
#' Graph management and Utilities:
#' \itemize{
#'    \item \code{\link{trim_bipartite_terms}}: Trim high-degree terms to improve fast approximation suitability.
#'    \item \code{\link{remove_isolated_elements}}: Remove elements with no term connections.
#'    \item \code{\link{remove_empty_terms}}: Remove terms with no element connections.
#'    \item \code{\link{reduce_universe_by_elements}}: Filter the graph to a specific set of elements.
#'    \item \code{\link{reduce_universe_by_terms}}: Filter the graph to a specific set of terms.
#'    \item \code{\link{summarize_suitability_bipartite}}: Check graph suitability for fast approximation.
#'    \item \code{\link{to_dataframe}}: Convert the bipartite graph to a data frame edge list.
#'    \item \code{\link{show}}: Summary display of the object.
#' }
#'
#' Underlying graph classes: \code{\link{ECGraph}}, \code{\link{ECProb}}
#' @export
#' @examples
#' # Create a sample term-element membership data frame
#' te_df <- data.frame(
#'    term = c("TermA", "TermA", "TermB", "TermC", "TermC"),
#'    element = c("Elem1", "Elem2", "Elem2", "Elem3", "Elem4")
#' )
#'
#' # Create an ECTermScoring object
#' ects <- ECTermScoring(te_df)
#' print(ects@terms)
#' print(ects@elements)
#'
ECTermScoring <- function(term_element_edges, col_term = 1, col_element = 2, ...) {

  if (is.character(term_element_edges) && length(term_element_edges) == 1) {
    # File path
    if (!file.exists(term_element_edges)) {
      stop("File not found: ", term_element_edges)
    }
    default_read_args <- list(
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = ""
    )
    user_read_args <- list(...)
    final_read_args <- utils::modifyList(default_read_args, user_read_args)
    final_read_args$file <- term_element_edges

    tryCatch({
      edge_df <- do.call(utils::read.table, final_read_args)
    }, error = function(e) {
      stop("Error reading file '", term_element_edges, "': ", e$message,
           "\nArguments used for read.table: ", paste(names(final_read_args), final_read_args, sep="=", collapse=", "))
    })
  } else if (is.data.frame(term_element_edges)) {
    edge_df <- term_element_edges
  } else {
    stop("'term_element_edges' must be a data frame or a file path.")
  }
  if (ncol(edge_df) < 2) {
    stop("Input 'term_element_edges' must have at least two columns.")
  }

  validate_col_spec <- function(col_spec, df_names, df_ncols, col_name_for_error) {
    if (is.numeric(col_spec)) {
      if (col_spec < 1 || col_spec > df_ncols) {
        stop(col_name_for_error, " index (", col_spec, ") is out of bounds for the data frame with ", df_ncols, " columns.")
      }
    } else if (is.character(col_spec)) {
      if (!(col_spec %in% df_names)) {
        stop(col_name_for_error, " '", col_spec, "' not found in data frame column names: ", paste(df_names, collapse=", "))
      }
    } else {
      stop(col_name_for_error, " must be a numeric index or a character column name.")
    }
  }
  validate_col_spec(col_term, names(edge_df), ncol(edge_df), "col_term")
  validate_col_spec(col_element, names(edge_df), ncol(edge_df), "col_element")

  # Vertex IDs are character
  term_ids_char <- as.character(edge_df[[col_term]])
  element_ids_char <- as.character(edge_df[[col_element]])
  bipartite_edge_list_for_ecgraph <- data.frame(
    vertex1 = term_ids_char,
    vertex2 = element_ids_char,
    stringsAsFactors = FALSE
  )

  # ECGraph and ECProb objects for the bipartite graph
  ecgraph_bipartite <- ECGraph(edges_input = bipartite_edge_list_for_ecgraph, col1="vertex1", col2="vertex2")
  ecprob_bipartite <- ECProb(ecgraph_bipartite)

  # Identify unique element and term vertices from the original input columns
  potential_element_vertices <- unique(element_ids_char)
  potential_term_vertices <- unique(term_ids_char)

  # Filter to those actually present in the graph
  actual_graph_vertices <- ecprob_bipartite@names
  final_element_vertices <- potential_element_vertices[potential_element_vertices %in% actual_graph_vertices]
  final_term_vertices    <- potential_term_vertices[potential_term_vertices %in% actual_graph_vertices]

  # A vertex cannot be both an element and a term
  if (any(final_term_vertices %in% final_element_vertices)) {
    overlapping_vertices <- intersect(final_term_vertices, final_element_vertices)
    stop("Some vertex IDs are identified as both elements and terms within the constructed graph: ",
         paste(utils::head(overlapping_vertices, 5), collapse=", "),
         ifelse(length(overlapping_vertices) > 5, "...", ""))
  }

  new("ECTermScoring",
      ecprob    = ecprob_bipartite,
      elements  = final_element_vertices,
      terms     = final_term_vertices
  )
}

#' @describeIn ECTermScoring Show method for ECTermScoring objects.
#' @param object An ECTermScoring object
setMethod("show", "ECTermScoring", function(object) {
  n_terms <- length(object@terms)
  n_elems <- length(object@elements)
  # In a bipartite graph context, M is the total number of connections
  n_edges <- object@ecprob@graph_size

  cat("An object of class \"ECTermScoring\" (Bipartite)\n")
  cat(paste0("------------------------------------------------\n"))
  cat(paste0("Statistics:\n"))
  cat(paste0("  Terms:    ", format(n_terms, big.mark=","), "\n"))
  cat(paste0("  Elements: ", format(n_elems, big.mark=","), "\n"))
  cat(paste0("  Edges:    ", format(n_edges, big.mark=","), "\n"))

  if (n_terms > 0) {
    cat(paste0("  Terms ex: ", paste(head(object@terms, 4), collapse=", "),
               if(n_terms > 4) "..." else "", "\n"))
  }

  cat(paste0("\nSlot Summary:\n"))
  cat(paste0("  @ecprob   : <ECProb> object (The underlying graph engine)\n"))
  cat(paste0("  @terms    : Character vector [", format(n_terms, big.mark=","), "]\n"))
  cat(paste0("  @elements : Character vector [", format(n_elems, big.mark=","), "]\n"))
  cat(paste0("------------------------------------------------\n"))
})

# --------------------------------------------------------------------------- #
# Internal Calculation Engine (Not Exported)
# --------------------------------------------------------------------------- #

# Internal helper: calculates stats but keeps result as a flat data.table
# to_dataframe/edges extraction happens OUTSIDE this function to optimize repeated calls.
calc_ecset_stats_core <- function(object, input_sets_dt, bipartite_edges_dt) {

  # 1.
  input_sets_dt_unique <- unique(input_sets_dt, by = c("set_id", "element"))
  data.table::setkey(input_sets_dt_unique, element)
  data.table::setkey(bipartite_edges_dt, element)

  # 2.
  all_connections <- bipartite_edges_dt[input_sets_dt_unique, on = "element", nomatch = 0, allow.cartesian = TRUE]

  # 3.
  observed_edges_dt <- all_connections[, .(observed_edge_count = .N), by = .(input_set_id = set_id, term_id)]

  # 4.
  all_element_degrees <- unlist(object@ecprob@degrees)

  term_degrees <- all_element_degrees[bipartite_edges_dt[, unique(term_id)]]
  term_summary <- data.table::data.table(term_id = names(term_degrees), term_degree = term_degrees)

  # 5.
  valid_input_elements <- all_connections[, .(input_set_id = set_id, element)] |> unique()
  input_set_summary <- valid_input_elements[,
                                            .(sum_degrees_set = sum(all_element_degrees[element], na.rm = TRUE),
                                              set_size = .N),
                                            by = input_set_id
  ]

  # 6.
  final_dt <- data.table::copy(observed_edges_dt)
  final_dt[term_summary, on = "term_id", term_degree := i.term_degree]
  final_dt[input_set_summary, on = "input_set_id", `:=`(sum_degrees_set = i.sum_degrees_set, max_possible_edges = i.set_size)]

  final_dt[, lambda := (term_degree * sum_degrees_set) / (2 * object@ecprob@graph_size)]

  final_dt[, p_value := calculate_p_value(object@ecprob, observed_edge_count, max_possible_edges, lambda)]
  final_dt[, log2_Anscombe_ratio := 0.5 * (log2(observed_edge_count + 3/8) - log2(lambda + 3/8))]

  data.table::setnames(final_dt, "input_set_id", "set1")
  data.table::setnames(final_dt, "term_id", "set2")

  return(final_dt)
}

#' @title Score Terms Against an Element Set
#'
#' @description Calculates enrichment statistics (p-value, z-score, etc.) for
#' all terms connected to a given set of elements.
#'
#' @param object An ECTermScoring object.
#' @param element_set A character vector of element names.
#' @param lambda_method The method for lambda calculation ("accurate", "optimized", "fast").
#'
#' @return A data frame where each row corresponds to a connected term,
#'   containing its enrichment statistics. Returns NULL if no valid elements or
#'   connected terms are found.
#' @export
#' @examples
#' te_df <- data.frame(term=c("T1","T1","T2"), element=c("E1","E2","E2"))
#' ects <- ECTermScoring(te_df)
#' # Score terms connected to the set {E1, E2}
#' terms_ecset_statistics(ects, element_set = c("E1", "E2"))
#'
setGeneric("terms_ecset_statistics",
           function(object, element_set, lambda_method = "fast")
             standardGeneric("terms_ecset_statistics"))
#' @describeIn terms_ecset_statistics Method for ECTermScoring objects.
setMethod(
  "terms_ecset_statistics",
  "ECTermScoring",
  function(object, element_set, lambda_method = "fast") {

    valid_element_set <- unique(element_set[element_set %in% object@elements])
    if (length(valid_element_set) < 1){
      warning("No valid elements from the input set found in the ECTermScoring object.")
      return(NULL)
    }

    connected_terms <- get_neighbors(object@ecprob, valid_element_set)
    relevant_terms <- intersect(connected_terms, object@terms)
    if (length(relevant_terms) < 1){
      warning("No terms are connected to the provided element set.")
      return(NULL)
    }

    ect_stats_single_term <- function(ecprob_obj, single_term_id, current_element_set, current_lambda_method) {
      max_possible_edges <- length(current_element_set)
      lambda <- switch(current_lambda_method,
                       accurate = calculate_lambda_between_naive(ecprob_obj, c(single_term_id), current_element_set),
                       optimized = calculate_lambda_between(ecprob_obj, c(single_term_id), current_element_set),
                       fast = calculate_lambda_between_fast(ecprob_obj, c(single_term_id), current_element_set),
                       stop("Invalid lambda_method specified in ect_stats_single_term.")
      )
      observed_edges <- get_edge_count_between(ecprob_obj, c(single_term_id), current_element_set)
      p_value <- calculate_p_value(ecprob_obj, observed_edges, max_possible_edges, lambda)
      observed_edge_count <- observed_edges
      log2_Anscombe_ratio <- NA_real_
      if (!is.na(lambda) && (lambda + 3/8) > 0 && (observed_edges + 3/8) > 0) {
        log2_Anscombe_ratio <- 0.5 * (log2(observed_edges + 3/8) - log2(lambda + 3/8))
      }
      return(list(p_value = p_value,
                  observed_edge_count = observed_edge_count,
                  lambda = lambda,
                  log2_Anscombe_ratio = log2_Anscombe_ratio
                  ))
    }

    vectorized_stats_calculator <- Vectorize(
      ect_stats_single_term,
      vectorize.args = "single_term_id",
      SIMPLIFY = FALSE
    )

    all_term_scores_list <- vectorized_stats_calculator(
      ecprob_obj = object@ecprob,
      single_term_id = relevant_terms,
      current_element_set = valid_element_set,
      current_lambda_method = lambda_method
    )

    names(all_term_scores_list) <- relevant_terms
    if (length(all_term_scores_list) > 0) {
      results_df <- do.call(rbind, lapply(names(all_term_scores_list), function(term_name) {
        data.frame(
          term_id = term_name,
          p_value = all_term_scores_list[[term_name]]$p_value,
          lambda = all_term_scores_list[[term_name]]$lambda,
          observed_edge_count = all_term_scores_list[[term_name]]$observed_edge_count,
          log2_Anscombe_ratio = all_term_scores_list[[term_name]]$log2_Anscombe_ratio,
          set_size = length(valid_element_set),
          stringsAsFactors = FALSE
        )
      }))
      return(results_df)
    } else {
      return(NULL)
    }
  })

#' @title Vectorized Scoring of Terms Against Multiple Element Sets
#'
#' @description A vectorized function that calculates enrichment
#' statistics for all terms against a collection of input element sets.
#'
#' @details This function is designed for efficiency when testing many input sets
#' against the term database at once.
#'
#' @param object An ECTermScoring object.
#' @param input_sets An object defining the input element sets. This can be either:
#'   \itemize{
#'     \item A `data.table` with two columns: `set_id` and `element`.
#'     \item A named `list` where names are the set IDs and values are character
#'       vectors of elements.
#'   }
#' @param lambda_method The method for lambda calculation. Only "fast"
#'   is supported for this vectorized function.
#'
#' @return A named list of `data.table`s. Each name corresponds to an
#'   `input_set_id`, and each `data.table` contains the connectedness statistics
#'   for all connected terms (see terms_ecset_statistics()).
#' @export
#' @examples
#' # Load sample data included with the package
#' library(data.table)
#' data(sample_ects)
#'
#' # --- Example 1: Using a data.table as input ---
#' input_dt <- data.table(
#'   set_id = c("Set_A", "Set_A", "Set_B"),
#'   element = c(sample_ects@elements[1], sample_ects@elements[2], sample_ects@elements[3])
#' )
#' results_dt <- terms_ecset_statistics_vectorized(sample_ects, input_dt)
#' print(results_dt$Set_A)
#'
#' # --- Example 2: Using a named list as input ---
#' input_list <- list(
#'   Set_A = c(sample_ects@elements[1], sample_ects@elements[2]),
#'   Set_B = c(sample_ects@elements[3])
#' )
#' results_list <- terms_ecset_statistics_vectorized(sample_ects, input_list)
#' print(results_list$Set_A)
#'
setGeneric("terms_ecset_statistics_vectorized",
           function(object, input_sets, lambda_method = "fast")
             standardGeneric("terms_ecset_statistics_vectorized"))
#' @describeIn terms_ecset_statistics_vectorized Method for ECTermScoring objects.
setMethod("terms_ecset_statistics_vectorized",
          "ECTermScoring",
          function(object, input_sets, lambda_method = "fast") {

             if (is.list(input_sets) && !is.data.frame(input_sets)) {
              input_sets_dt <- data.table::as.data.table(utils::stack(input_sets))
              data.table::setnames(input_sets_dt, c("values", "ind"), c("element", "set_id"))
            } else {
              input_sets_dt <- data.table::as.data.table(input_sets)
            }

            bipartite_edges <- data.table::as.data.table(to_dataframe(object))
            data.table::setnames(bipartite_edges, "term", "term_id")
            bipartite_edges[, term_id := as.character(term_id)]

            final_dt <- calc_ecset_stats_core(object, input_sets_dt, bipartite_edges)

            return(split(final_dt, by = "set1"))
          })

# --------------------------------------------------------------------------- #
# Empirical FDR Method
# --------------------------------------------------------------------------- #

#' @title Calculate Empirical FDR for Set Enrichment
#'
#' @description Calculates Normalized Enrichment Scores (NES) and False Discovery
#' Rates (FDR) for term enrichment by comparing observed scores against a null
#' distribution generated from random sets.
#'
#' @details This method addresses the non-uniform distribution of p-values in
#' sparse graphs under the RGGED null model.
#'
#' @param object An \code{\link{ECTermScoring}} object.
#' @param real_results_dt A \code{data.table} or \code{data.frame} containing the
#' scoring results for a single set (e.g., the output of \code{\link{terms_ecset_statistics}}).
#' Must contain columns `observed_edge_count` and `log2_Anscombe_ratio`.
#' @param n_permutations Integer. The number of random sets to generate for the
#' null distribution. Defaults to 1000.
#' @param seed Integer (Optional). A random seed for reproducibility.
#'
#' @return A \code{data.table} containing the original results with two additional columns:
#' \itemize{
#'   \item \code{nes}: Normalized Enrichment Score.
#'   \item \code{fdr_q_value}: The empirical False Discovery Rate q-value.
#' }
#' The table is sorted by \code{fdr_q_value} (ascending) and \code{nes} (descending).
#'
#' @seealso \code{\link{terms_ecset_statistics}}, \code{\link{run_vsea}}
#' @export
#' @import data.table
calculate_ecset_fdr <- function(object, real_results_dt, n_permutations = 1000, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  real_results_dt <- data.table::as.data.table(real_results_dt)
  real_results_dt <- data.table::copy(real_results_dt)

  if ("term_id" %in% names(real_results_dt) && !"set2" %in% names(real_results_dt)) {
    data.table::setnames(real_results_dt, "term_id", "set2")
  }

  if ("set_size" %in% names(real_results_dt)) {
    set_size <- max(real_results_dt$set_size)
  } else {
    stop("Input results must contain a 'set_size' column to define the null model.")
  }

  if (set_size < 1) {
    warning("Set size is 0 after filtering. Cannot perform permutations.")
    return(real_results_dt)
  }

  bipartite_edges <- data.table::as.data.table(to_dataframe(object))
  data.table::setnames(bipartite_edges, "term", "term_id")
  bipartite_edges[, term_id := as.character(term_id)]

  random_elements <- unlist(replicate(n_permutations, sample(object@elements, set_size), simplify = FALSE))

  null_inputs_dt <- data.table::data.table(
    set_id = rep(paste0("R", 1:n_permutations), each = set_size),
    element = random_elements
  )

  null_stats_dt <- calc_ecset_stats_core(object, null_inputs_dt, bipartite_edges)

  score_col <- "log2_Anscombe_ratio"

  mean_nulls <- null_stats_dt[, .(
    mean_pos = mean(get(score_col)[get(score_col) >= 0], na.rm = TRUE),
    mean_neg = mean(abs(get(score_col)[get(score_col) < 0]), na.rm = TRUE)
  ), by = .(term_id = set2)]

  mean_nulls[is.nan(mean_pos), mean_pos := 1]
  mean_nulls[is.nan(mean_neg), mean_neg := 1]

  real_results_dt[mean_nulls, on = .(set2 = term_id), `:=`(mean_null_pos = i.mean_pos, mean_null_neg = i.mean_neg)]

  real_results_dt[is.na(mean_null_pos), mean_null_pos := 1]
  real_results_dt[is.na(mean_null_neg), mean_null_neg := 1]

  real_results_dt[, nes := ifelse(get(score_col) >= 0,
                                  get(score_col) / mean_null_pos,
                                  get(score_col) / mean_null_neg)]

  null_stats_dt[mean_nulls, on = .(set2 = term_id), `:=`(mean_null_pos = i.mean_pos, mean_null_neg = i.mean_neg)]
  null_stats_dt[is.na(mean_null_pos), mean_null_pos := 1]
  null_stats_dt[is.na(mean_null_neg), mean_null_neg := 1]

  null_stats_dt[, null_nes := ifelse(get(score_col) >= 0,
                                     get(score_col) / mean_null_pos,
                                     get(score_col) / mean_null_neg)]

  all_null_nes <- na.omit(null_stats_dt$null_nes)
  null_nes_pos <- all_null_nes[all_null_nes >= 0]
  null_nes_neg <- all_null_nes[all_null_nes < 0]

  nes_real_vec <- real_results_dt$nes

  if(any(is.na(nes_real_vec))) {
    warning("NA values found in NES. Replacing with 0 for FDR calculation.")
    nes_real_vec[is.na(nes_real_vec)] <- 0
  }

  fdr_values <- vapply(nes_real_vec, function(score) {
    if (score >= 0) {
      if (length(null_nes_pos) == 0) return(1)
      frac_null <- (sum(null_nes_pos >= score, na.rm = TRUE) + 1) / (length(null_nes_pos) + 1)
      frac_real <- sum(nes_real_vec >= 0 & nes_real_vec >= score, na.rm = TRUE) / sum(nes_real_vec >= 0, na.rm = TRUE)

      return(if (frac_real == 0) 1 else frac_null / frac_real)
    } else {
      if (length(null_nes_neg) == 0) return(1)
      frac_null <- (sum(null_nes_neg <= score, na.rm = TRUE) + 1) / (length(null_nes_neg) + 1)
      frac_real <- sum(nes_real_vec < 0 & nes_real_vec <= score, na.rm = TRUE) / sum(nes_real_vec < 0, na.rm = TRUE)

      return(if (frac_real == 0) 1 else frac_null / frac_real)
    }
  }, FUN.VALUE = numeric(1))

  fdr_values[fdr_values > 1] <- 1
  real_results_dt[, fdr_q_value := fdr_values]

  real_results_dt[, c("mean_null_pos", "mean_null_neg") := NULL]
  data.table::setorder(real_results_dt, fdr_q_value, -nes)

  if ("set2" %in% names(real_results_dt)) {
    data.table::setnames(real_results_dt, "set2", "term_id")
  }
  setorder(real_results_dt, -log2_Anscombe_ratio)
  return(real_results_dt)
}

#' @title Score Terms Against an Element Set with FDR Correction
#'
#' @description A convenience wrapper that performs enrichment analysis
#' for a single set of elements. It calculates raw enrichment statistics AND
#' computes Normalized Enrichment Scores (NES) and Empirical False
#' Discovery Rates (FDR).
#'
#' @details This function enforces the "fast" lambda approximation method for both
#' the real set and the random simulations.
#'
#' @param object An \code{\link{ECTermScoring}} object.
#' @param element_set A character vector of element names.
#' @param n_permutations Integer. The number of random sets to generate for the
#' null distribution. Defaults to 1000.
#' @param seed Integer (Optional). A random seed for reproducibility.
#'
#' @return A \code{data.table} containing the enrichment results, including:
#' \itemize{
#'   \item \code{term_id}: The term identifier.
#'   \item \code{observed_edge_count}: Number of edges between the set and the term.
#'   \item \code{log2_Anscombe_ratio}: The effect size.
#'   \item \code{p_value}: The raw RGGED p-value (fast approximation).
#'   \item \code{nes}: Normalized Enrichment Score.
#'   \item \code{fdr_q_value}: The empirical False Discovery Rate q-value.
#' }
#' Returns \code{NULL} if no valid connected terms are found.
#'
#' @seealso \code{\link{terms_ecset_statistics}}, \code{\link{calculate_ecset_fdr}}
#' @export
#' @examples
#' data(sample_ects)
#' target_term <- "GO:0005787"
#' set_1 <- sample_ects@ecprob@adj[[target_term]]
#'
#' # Run full analysis in one step
#' results <- terms_ecset_statistics_fdr(sample_ects, set_1, n_permutations = 100, seed = 123)
#' print(head(results))
#'
setGeneric("terms_ecset_statistics_fdr",
           function(object, element_set, n_permutations = 1000, seed = NULL)
             standardGeneric("terms_ecset_statistics_fdr"))

#' @describeIn terms_ecset_statistics_fdr Method for ECTermScoring objects.
setMethod("terms_ecset_statistics_fdr",
          "ECTermScoring",
          function(object, element_set, n_permutations = 1000, seed = NULL) {

            raw_results <- terms_ecset_statistics(object, element_set, lambda_method = "fast")

            if (is.null(raw_results)) {
              return(NULL)
            }

            final_results <- calculate_ecset_fdr(object,
                                                     raw_results,
                                                     n_permutations = n_permutations,
                                                     seed = seed)

            return(final_results)
          })

#' @title Score Terms Against a Ranked List of Elements
#'
#' @description Calculates a running enrichment score for terms based on a ranked list of elements.
#' @details This method implements a fast, vectorized algorithm. For statistical
#' rigor, the analysis is conditioned on the universe of elements provided in the
#' `element_ranks` argument.
#'
#' @param object An ECTermScoring object.
#' @param element_ranks A named list or vector of element ranks.
#'
#' @return A named list where each element is a `data.table` of statistics for a term.
#' @export
#' @examples
#' # Load sample data included with the package
#' data(sample_ects)
#'
#' # Create a ranked list using a subset of the elements
#' set.seed(1)
#' elements_subset <- sample(sample_ects@elements, 2000)
#' ranked_list <- setNames(seq_along(elements_subset), elements_subset)
#'
#' # Run the analysis
#' ranked_scores <- terms_ecranks_statistics(sample_ects, ranked_list)
#'
#' # View the results for a specific term
#' terms <- get_neighbors(sample_ects@ecprob, elements_subset[1])
#' term_sizes <- lengths(sample_ects@ecprob@adj[terms])
#' term <- terms[which.max(term_sizes)]
#' print(term)
#' print(ranked_scores[[term]])
#'
setGeneric("terms_ecranks_statistics",
           function(object, element_ranks)
             standardGeneric("terms_ecranks_statistics"))
#' @describeIn terms_ecranks_statistics Method for ECTermScoring objects.
setMethod(
  "terms_ecranks_statistics",
  "ECTermScoring",
  function(object, element_ranks) {

    valid_elements <- intersect(names(element_ranks), object@elements)
    if (length(valid_elements) < 1) {
      stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
    }
    object <- reduce_universe_by_elements(object, valid_elements)
    element_ranks <- rank(element_ranks[valid_elements])

    all_element_degrees <- unlist(object@ecprob@degrees)
    ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
    setorder(ranks_dt, global_rank)
    ranks_dt[, degree := all_element_degrees[element_id]]
    ranks_dt[, cumsum_degrees := cumsum(degree)]
    bipartite_edges <- as.data.table(to_dataframe(object))
    setnames(bipartite_edges, c("term", "element"), c("term_id", "element_id"))
    bipartite_edges[, term_id := as.character(term_id)]

    setkey(bipartite_edges, element_id)
    setkey(ranks_dt, element_id)
    final_dt <- ranks_dt[bipartite_edges, on = "element_id"]
    setorder(final_dt, term_id, global_rank)
    final_dt[, rank_in_term := 1:.N, by = term_id]
    term_sizes <- lengths(object@ecprob@adj[object@terms])
    term_summary <- data.table(
      term_id = names(term_sizes),
      term_size = term_sizes,
      term_degree = term_sizes
    )
    final_dt[term_summary, on = "term_id", `:=`(term_size = i.term_size, term_degree = i.term_degree)]

    final_dt[, `:=`(
      observed_ec = rank_in_term,
      max_ec = pmin(term_size, global_rank),
      lambda = (term_degree / (2 * object@ecprob@graph_size)) * cumsum_degrees
    )]
    final_dt[, `:=`(
      p_value = calculate_p_value(object@ecprob, observed_ec, max_ec, lambda),
      log2_Anscombe_ratio = 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))
    )]

    results_list <- split(final_dt, by = "term_id")
    final_cols <- c("element_id", "global_rank", "rank_in_term",
                    "observed_ec", "max_ec", "term_size",
                    "lambda", "p_value", "log2_Anscombe_ratio")
    results_list_final <- lapply(results_list, function(dt) {
      dt[, ..final_cols]
    })
    return(results_list_final)
  })

#' @title Full Summary of Rank Statistics (
#'
#' @description A function that takes the output of
#' `terms_ecranks_statistics` and calculates a summary (min, max, median
#' scores and their context) for each term.
#'
#' @param term_scores_list A named list of data.tables from `terms_ecranks_statistics`.
#' @param scoring_statistic The column name of the score to summarize.
#'
#' @return A `data.table` with one row per term, containing the rich summary.
summarize_ranks_full <- function(term_scores_list, scoring_statistic = "log2_Anscombe_ratio") {

  long_dt <- rbindlist(term_scores_list, idcol = "term_id")

  summary_dt <- long_dt[, {
    idx_min <- which.min(get(scoring_statistic))
    idx_max <- which.max(get(scoring_statistic))
    idx_median <- which.min(abs(get(scoring_statistic) - median(get(scoring_statistic), na.rm = TRUE)))
    .(
      term_size = term_size[1],
      min_score = get(scoring_statistic)[idx_min],
      element_at_min = element[idx_min],
      rank_at_min = global_rank[idx_min],
      median_score = get(scoring_statistic)[idx_median],
      element_at_median = element[idx_median],
      rank_at_median = global_rank[idx_median],
      max_score = get(scoring_statistic)[idx_max],
      element_at_max = element[idx_max],
      rank_at_max = global_rank[idx_max]
    )
  }, by = term_id]
  return(summary_dt)
}

#' @title Perform Vertex Set Enrichment Analysis
#'
#' @description Runs a ranked-list enrichment analysis (Vertex Set Enrichment Analysis),
#' calculating a Normalized Enrichment Score (NES) and a False Discovery Rate (FDR) q-value
#' for each term using a permutation-based null distribution.
#'
#' @details This version focuses solely on identifying positive enrichment,
#' i.e., the maximum score for each term
#'
#' @param object An ECTermScoring object.
#' @param element_ranks A named list or vector of element ranks.
#' @param scoring_statistic The score to use for enrichment (e.g., "log2_Anscombe_ratio").
#' @param n_permutations The number of permutations to generate the null
#'   distribution, 1000 by default.
#' @param seed An optional random seed for reproducibility.
#' @param n_cores Number of cores to use for parallel processing. Defaults to
#'   detectCores() - 1. Set to 1 for serial execution (or on Windows).
#'
#' @return A list containing a single data frame: \code{max_score_summary}.
#'   The data frame is sorted by \code{fdr_q_value} and contains:
#'   \itemize{
#'     \item \code{term_id}: The term identifier.
#'     \item \code{max_score}: The maximum running enrichment score observed.
#'     \item \code{nes}: Normalized Enrichment Score.
#'     \item \code{fdr_q_value}: Empirical False Discovery Rate.
#'     \item \code{term_size}: Number of elements in the term.
#'     \item \code{rank_at_max}: The element rank at which the score peaked.
#'     \item \code{observed_ec}: The number of observed edges at the peak.
#'   }
#'
#' @export
#' @examples
#' library(EdgeCount)
#' library(data.table)
#'
#' # Load data
#' data(sample_ects)
#'
#' # Random uniform ranking
#' element_ranks <- setNames(
#'   sample(seq_along(sample_ects@elements)),
#'   sample_ects@elements
#' )
#'
#' # VSEA analysis
#' vsea_results <- run_vsea(
#'   sample_ects,
#'   element_ranks,
#'   n_permutations = 10, # Low for speed
#'   n_cores = 2,
#'   seed = 1
#' )
#'
#' print(head(vsea_results$max_score_summary))
setGeneric("run_vsea",
           function(object, element_ranks, scoring_statistic = "log2_Anscombe_ratio", n_permutations = 1000, n_cores = NULL, seed = NULL)
             standardGeneric("run_vsea"))
#' @describeIn run_vsea Method for ECTermScoring objects.
setMethod("run_vsea",
          "ECTermScoring",
          function(object,
                   element_ranks,
                   scoring_statistic = "log2_Anscombe_ratio",
                   n_permutations = 1000,
                   n_cores = NULL,
                   seed = NULL) {

            old_threads <- data.table::getDTthreads()
            data.table::setDTthreads(1)
            on.exit(data.table::setDTthreads(old_threads))
            if (!is.null(seed)) set.seed(seed)
            if (is.null(n_cores)){
              n_cores <- parallel::detectCores() - 1
              if (is.na(n_cores) || n_cores < 1) n_cores <- 1
            }
            if (n_cores < 1) n_cores <- 1
            if (.Platform$OS.type == "windows") n_cores <- 1

            if (!is.null(names(element_ranks))) {
              element_ranking <- names(sort(element_ranks))
            } else {
               element_ranking <- as.character(element_ranks)
            }

            elements_char <- as.character(object@elements)
            valid_elements <- intersect(element_ranking, elements_char)
            if (length(valid_elements) < 1) {
              stop("No valid elements found in the ECTermScoring object.")
            }

            # reduce universe, if necessary
            if (length(valid_elements) < length(elements_char)) {
              object <- reduce_universe_by_elements(object, valid_elements)
            }
            elements <- object@elements
            terms    <- object@terms
            Ne <- length(elements)
            Nt <- length(terms)
            degrees <- unlist(object@ecprob@degrees[elements], use.names = FALSE)

            element_name_to_idx <- setNames(seq_len(Ne), elements)
            term_map <- setNames(seq_len(Nt), terms)
            aligned_adj <- object@ecprob@adj[elements]
            adj_int <- lapply(aligned_adj, function(x) {
              if (is.null(x) || length(x) == 0) return(integer(0))
              unname(term_map[as.character(x)])
            })

            M <- length(unlist(adj_int))
            u_terms <- unlist(adj_int, use.names = FALSE)
            global_kt <- tabulate(u_terms, nbins = Nt)

            # simulations
            total_len <- as.numeric(Ne) * n_permutations
            sim_vector <- integer(total_len)
            for (i in 1:n_permutations) {
              start_idx <- (i - 1) * Ne + 1
              end_idx   <- start_idx + Ne - 1
              sim_vector[start_idx:end_idx] <- dqrng::dqsample.int(Ne)
            }

            # cumulative degree matrix
            sim_degrees <- degrees[sim_vector]
            dim(sim_degrees) <- c(Ne, n_permutations)
            sim_cum_degrees <- matrixStats::colCumsums(sim_degrees)
            rm(sim_degrees); gc()

            # scoring function
            score_simulation <- function(sim_idx) {

              start_idx <- (sim_idx - 1) * Ne + 1
              end_idx   <- start_idx + Ne - 1
              perm_elements <- sim_vector[start_idx:end_idx]

              terms_in_rank_order <- adj_int[perm_elements]
              edge_terms <- unlist(terms_in_rank_order, use.names = FALSE)
              edge_ranks <- rep(seq_len(Ne), lengths(terms_in_rank_order))

              ord <- order(edge_terms, method = "radix")
              sorted_terms <- edge_terms[ord]
              sorted_ranks <- edge_ranks[ord]

              run_lens <- rle(sorted_terms)$lengths
              x_vals <- sequence(run_lens)
              term_kt_subset <- global_kt[sorted_terms]

              lambdas <- (sim_cum_degrees[sorted_ranks, sim_idx] * term_kt_subset) / (2 * M)
              scores <- 0.5 * (log2(x_vals + 0.375) - log2(lambdas + 0.375))

              dt <- data.table::data.table(tm = sorted_terms, sc = scores)
              return(dt[, max(sc), by = tm]$V1)
            }

            # score simulations
            results_list <- parallel::mclapply(1:n_permutations, score_simulation, mc.cores = n_cores)
            null_scores_pool <- unlist(results_list, use.names = FALSE)
            null_dist_sorted <- sort(null_scores_pool)

            obs_indices <- element_name_to_idx[element_ranking]
            obs_degrees <- degrees[obs_indices]
            real_cum_degrees <- cumsum(obs_degrees)

            terms_in_rank_order <- adj_int[obs_indices]
            edge_terms <- unlist(terms_in_rank_order, use.names = FALSE)
            edge_ranks <- rep(seq_len(Ne), lengths(terms_in_rank_order))

            ord <- order(edge_terms, method = "radix")
            sorted_terms <- edge_terms[ord]
            sorted_ranks <- edge_ranks[ord]

            run_lens <- rle(sorted_terms)$lengths
            x_vals <- sequence(run_lens)
            term_kt_subset <- global_kt[sorted_terms]

            lambdas <- (real_cum_degrees[sorted_ranks] * term_kt_subset) / (2 * M)
            scores <- 0.5 * (log2(x_vals + 0.375) - log2(lambdas + 0.375))

            dt_obs <- data.table::data.table(tm = sorted_terms, sc = scores, ec = x_vals, rk = sorted_ranks)

            obs_agg <- dt_obs[, .(
              max_score = max(sc),
              observed_ec = ec[which.max(sc)],
              rank_at_max = rk[which.max(sc)]
            ), by = tm]

            data.table::setorder(obs_agg, -max_score)

            # fdr
            n_null_total <- length(null_dist_sorted)
            n_null_ge <- n_null_total - findInterval(obs_agg$max_score - 1e-10, null_dist_sorted)

            obs_agg[, E_FalsePos := n_null_ge / n_permutations]
            obs_agg[, Rank := .I]
            obs_agg[, FDR_raw := E_FalsePos / Rank]
            obs_agg[FDR_raw > 1, FDR_raw := 1]
            obs_agg[, fdr_q_value := rev(cummin(rev(FDR_raw)))]
            obs_agg[, term_id := terms[tm]]

            null_pos_mean <- mean(null_dist_sorted[null_dist_sorted > 0])
            if (is.na(null_pos_mean)) null_pos_mean <- 1
            obs_agg[, nes := max_score / null_pos_mean]

            obs_agg[, term_size := global_kt[tm]]

            final_cols <- c("term_id", "max_score", "nes", "fdr_q_value", "term_size", "rank_at_max", "observed_ec")

            return(list(
              max_score_summary = obs_agg[, ..final_cols]
            ))
          })

#' @title Trim High-Degree Vertices from a Bipartite Graph
#'
#' @description Iteratively removes the highest-degree terms from an
#' ECTermScoring object until the condition for fast lambda approximation is met.
#' Elements are not directly removed, but may be excluded from the final object
#' if they become disconnected from all remaining terms.
#'
#' @details This function provides a data-driven way to prepare a bipartite graph
#' for analyses that use fast lambda calculation methods. The trimming condition is
#' `(max_ke * max_kt) < threshold * 2M`, where `max_ke` is the current maximum
#' element degree, `max_kt` is the current maximum term degree, and `M` is the
#' current total number of edges. In each step, it removes the term(s) with the
#' highest degree and dynamically updates the degrees of the affected elements.
#'
#' @param object An ECTermScoring object.
#' @param threshold A numeric threshold for the trimming condition. The default
#'   value of 1 corresponds to the theoretical condition. Lowering this value will
#'   result in more aggressive trimming.
#'
#' @return A list containing three elements:
#' \itemize{
#'   \item `trimmed_object`: A new, smaller ECTermScoring object.
#'   \item `removed_terms`: A character vector of the term names that were actively removed.
#'   \item `removed_elements`: A character vector of the element names that were
#'     removed as a consequence of becoming disconnected from all terms.
#' }
#' @export
setGeneric("trim_bipartite_terms", function(object, threshold = 1.0) standardGeneric("trim_bipartite_terms"))
#' @describeIn trim_bipartite_terms Method for ECTermScoring objects.
setMethod("trim_bipartite_terms", "ECTermScoring", function(object, threshold = 1.0) {

  all_degrees <- unlist(object@ecprob@degrees)
  current_term_degrees <- all_degrees[object@terms]
  current_element_degrees <- all_degrees[object@elements]
  two_m_current <- object@ecprob@graph_size * 2
  removed_terms <- character(0)

  while (length(current_term_degrees) > 0 && length(current_element_degrees) > 0) {
    max_kt <- max(current_term_degrees)
    max_ke <- max(current_element_degrees)
    if (max_ke * max_kt < threshold * two_m_current) {
      break
    }
    terms_to_remove <- names(current_term_degrees[current_term_degrees == max_kt])
    removed_terms <- c(removed_terms, terms_to_remove)
    neighbors_to_update <- unlist(object@ecprob@adj[terms_to_remove])
    if(length(neighbors_to_update) > 0) {
      neighbor_counts <- table(neighbors_to_update)
      affected_elements <- names(neighbor_counts)
      current_element_degrees[affected_elements] <- current_element_degrees[affected_elements] - neighbor_counts
    }
    two_m_current <- two_m_current - (length(terms_to_remove) * max_kt * 2) # Each edge removal reduces sum of degrees by 2
    current_term_degrees <- current_term_degrees[!names(current_term_degrees) %in% terms_to_remove]
    current_element_degrees <- current_element_degrees[current_element_degrees > 0] # Remove elements that are now disconnected
  }

  kept_terms <- setdiff(object@terms, removed_terms)
  if (length(kept_terms) == 0) {
    empty_df <- data.frame(term=character(), element=character())
    final_removed_elements <- object@elements
    return(list(
      trimmed_object = ECTermScoring(empty_df),
      removed_terms = removed_terms,
      removed_elements = final_removed_elements
    ))
  }

  term_neighbors_list <- object@ecprob@adj[kept_terms]
  new_edge_df <- utils::stack(term_neighbors_list)
  names(new_edge_df) <- c("element", "term")

  new_ects <- ECTermScoring(new_edge_df[, c("term", "element")])
  final_removed_elements <- setdiff(object@elements, new_ects@elements)
  return(list(
    trimmed_object = new_ects,
    removed_terms = removed_terms,
    removed_elements = final_removed_elements
  ))
})

#' title Remove Isolated Elements from an ECTermScoring Object
#'
#' @description Creates a new ECTermScoring object that excludes all elements with
#' zero degree (no element). This method only removes isolated *elements*;
#' isolated *terms* (terms with no associated elements) are kept.
#'
#' @param object An ECTermScoring object.
#'
#' @return A new, smaller ECTermScoring object.
#' @export
#' @examples
#' # Create an object where E3 is an isolated element
#' te_df <- data.frame(
#'   term = c("T1", "T1"),
#'   element = c("E1", "E2")
#' )
#' # Manually create a universe that includes an extra element
#' full_universe_df <- data.frame(
#'    v1 = c("T1", "T1", "T2"), # T2 is an isolated term
#'    v2 = c("E1", "E2", "E3")  # E3 is an isolated element
#' )
#' ecg <- ECGraph(full_universe_df)
#' ects_base <- ECTermScoring(te_df) # E3 is not in this edge list
#'
setGeneric("remove_isolated_elements", function(object) standardGeneric("remove_isolated_elements"))
#' @describeIn remove_isolated_elements Method for ECTermScoring objects.
setMethod("remove_isolated_elements", "ECTermScoring", function(object) {

  element_degrees <- unlist(object@ecprob@degrees[object@elements])
  kept_elements <- names(element_degrees[element_degrees > 0])
  if (length(kept_elements) == length(object@elements)) {
    return(object)
  }
  term_neighbors_list <- object@ecprob@adj[object@terms]
  term_neighbors_clean <- lapply(term_neighbors_list, function(neighbors) {
    neighbors[neighbors %in% kept_elements]
  })
  new_edge_df <- utils::stack(term_neighbors_clean)
  if(nrow(new_edge_df) == 0) {
    final_df <- data.frame(term=character(), element=character())
  } else {
    names(new_edge_df) <- c("element", "term")
    final_df <- new_edge_df[, c("term", "element")]
  }
  ECTermScoring(final_df)
})

#' @title Remove Empty Terms from an ECTermScoring Object
#'
#' @description Creates a new ECTermScoring object that excludes all terms with a
#' degree of zero (i.e., terms with no associated elements).
#'
#' @param object An ECTermScoring object.
#'
#' @return A new, smaller ECTermScoring object with empty terms removed.
#' @export
#' @examples
#' # Create an object where T2 is an empty term
#' te_df <- data.frame(
#'   term = c("T1", "T1", "T2"),
#'   element = c("E1", "E2", "E1") # T2 has an edge, but let's make it empty
#' )
#' # To properly create an empty term, we build a graph with it
#' full_universe_df <- data.frame(
#'    v1 = c("T1", "T1", "T2"),
#'    v2 = c("E1", "E2", "E1")
#' )
#' ecg <- ECGraph(full_universe_df) # ECGraph constructor will find T2
#'
#' # Create the ECTermScoring object from an edge list that omits T2's edges
#' ects_with_empty_term <- ECTermScoring(data.frame(term="T1", element=c("E1","E2")))
#'
#' # For this example, let's assume `ects_with_empty_term` was created
#' # in a way that it knows about "T2" but has no elements for it.
#' # ects_trimmed <- remove_empty_terms(ects_with_empty_term)
#' # print(ects_trimmed@terms) # Would not contain "T2"
#'
setGeneric("remove_empty_terms", function(object) standardGeneric("remove_empty_terms"))
#' @describeIn remove_empty_terms Method for ECTermScoring objects.
setMethod("remove_empty_terms", "ECTermScoring", function(object) {

  term_degrees <- unlist(object@ecprob@degrees[object@terms])
  kept_terms <- names(term_degrees[term_degrees > 0])
  if (length(kept_terms) == length(object@terms)) {
    return(object)
  }
  term_neighbors_list <- object@ecprob@adj[kept_terms]
  new_edge_df <- utils::stack(term_neighbors_list)
  if(nrow(new_edge_df) == 0) {
    final_df <- data.frame(term=character(), element=character())
  } else {
    names(new_edge_df) <- c("element", "term")
    final_df <- new_edge_df[, c("term", "element")]
  }
  ECTermScoring(final_df)
})

#' @title Convert ECTermScoring to Data Frame
#'
#' @description Creates a data frame of edges from an ECTermScoring object.
#' The returned data frame has two columns representing term-element edges.
#' Isolated elements (no associated terms) and empty terms (no associated elements)
#' are not included in the output.
#'
#' @param object An ECTermScoring object.
#'
#' @return A data frame with two columns "term" and "element" representing
#' the edge set of the bipartite graph.
#' @export
#' @examples
#' # The to_dataframe generic is defined with the ECGraph class.
#' # This method provides the implementation for ECTermScoring objects.
#' te_df <- data.frame(
#'   term = c("T1", "T1", "T2"),
#'   element = c("E1", "E2", "E2")
#' )
#' ects <- ECTermScoring(te_df)
#' edge_list <- to_dataframe(ects)
#' print(edge_list)
#'
#' @describeIn to_dataframe Method for ECTermScoring objects.
setMethod("to_dataframe",
          "ECTermScoring",
          function(object) {
            edge_df <- utils::stack(object@ecprob@adj[object@terms])
            if(nrow(edge_df) == 0) {
              return(data.frame(term = character(), element = character()))
            }
            names(edge_df) <- c("element", "term")
            edge_df[, c("term", "element")]
          })
#' @title Reduce the Universe of an ECTermScoring Object to a Set of Elements
#'
#' @description Creates a new, smaller ECTermScoring object by restricting it
#' to a specified subset of elements.
#'
#' @details This is a helper for ensuring statistical tests are correctly
#' conditioned. It takes an ECTermScoring object and a vector of elements, finds
#' the intersection, and returns a new object where all properties (degrees, graph
#' size, etc.) are recalculated based on this smaller universe. This method is an
#' internal helper for rank scoring.
#'
#' @param object An ECTermScoring object.
#' @param elements_to_keep A character vector of element IDs to keep.
#'
#' @return A new, smaller ECTermScoring object.
#' @export
setGeneric("reduce_universe_by_elements", function(object, elements_to_keep) standardGeneric("reduce_universe_by_elements"))
#' @describeIn reduce_universe_by_elements Method for ECTermScoring objects.
setMethod("reduce_universe_by_elements", "ECTermScoring", function(object, elements_to_keep) {

  bipartite_edges <- as.data.table(to_dataframe(object))
  reduced_edges <- bipartite_edges[element %in% elements_to_keep]
  return(ECTermScoring(reduced_edges))
})

#' @title Reduce the Universe of an ECTermScoring Object to a Set of Terms
#'
#' @description Creates a new, smaller ECTermScoring object by restricting it
#' to a specified subset of terms.
#'
#' @details This is a convenience method which also ensures that statistical tests are correctly
#' conditioned. It takes an ECTermScoring object and a vector of terms, finds
#' the intersection, and returns a new object where all properties (degrees, graph
#' size, etc.) are recalculated based on the input constraint.
#'
#' @param object An ECTermScoring object.
#' @param terms_to_keep A character vector of term IDs to keep.
#'
#' @return A new, smaller ECTermScoring object.
#' @export
setGeneric("reduce_universe_by_terms", function(object, terms_to_keep) standardGeneric("reduce_universe_by_terms"))
#' @describeIn reduce_universe_by_terms Method for ECTermScoring objects.
setMethod("reduce_universe_by_terms", "ECTermScoring", function(object, terms_to_keep) {
  bipartite_edges <- as.data.table(to_dataframe(object))
  reduced_edges <- bipartite_edges[term %in% terms_to_keep]
  return(ECTermScoring(reduced_edges))
})
#' @title Summarize Bipartite Graph Suitability for Fast Lambda Approximation
#'
#' @description Provides statistics on the pairwise Bernoulli parameters (p_ij)
#' specifically for Term-Element pairs. This helps assess if fast lambda
#' approximation methods are suitable for the bipartite graph structure.
#'
#' @details Unlike \code{summarize_suitability_fast} which checks all vertex pairs,
#' this method only checks pairs where one vertex is a Term and the other is an Element.
#' This is the relevant check for bipartite enrichment analyses.
#'
#' @param object An ECTermScoring object.
#'
#' @return A list of summary statistics:
#' \itemize{
#'   \item `pij_over_1`: The proportion of (Term, Element) pairs such that p_ij >= 1.
#'   \item `prop_problematic_terms`: Proportion of terms with degree > sqrt(2*M).
#'   \item `prop_problematic_elements`: Proportion of elements with degree > sqrt(2*M).
#' }
#' @export
setGeneric("summarize_suitability_bipartite",
           function(object) standardGeneric("summarize_suitability_bipartite"))
#' @describeIn summarize_suitability_bipartite Method for ECTermScoring objects.
setMethod("summarize_suitability_bipartite",
          "ECTermScoring",
          function(object) {

            all_degrees <- unlist(object@ecprob@degrees)
            term_degrees <- all_degrees[object@terms]
            element_degrees <- all_degrees[object@elements]
            Nt <- length(term_degrees)
            Ne <- length(element_degrees)
            if (Nt == 0 || Ne == 0) {
              return(list(
                pij_over_1 = 0,
                prop_problematic_terms = 0,
                prop_problematic_elements = 0
              ))
            }
            M <- object@ecprob@graph_size
            two_M <- 2 * M
            t_dist <- table(term_degrees)
            k_t <- as.numeric(names(t_dist))
            n_t <- as.numeric(t_dist)
            e_dist <- table(element_degrees)
            k_e <- as.numeric(names(e_dist))
            n_e <- as.numeric(e_dist)
            p_t <- n_t / Nt
            p_e <- n_e / Ne
            prod_matrix <- outer(k_t, k_e, "*") / two_M
            prob_matrix <- outer(p_t, p_e, "*")
            count_matrix <- outer(n_t, n_e)
            count_vec <- as.vector(count_matrix)
            prod_vec <- as.vector(prod_matrix)
            prob_vec <- as.vector(prob_matrix)
            dt <- data.table(pij = prod_vec,
                             count = count_vec,
                             prop = prob_vec)
            dt <- dt[, .(
              count = as.integer(sum(count)),
              prop = sum(prop)
            ), by = .(pij)]
            pij_over_1 <- sum(prob_vec[prod_vec >= 1])
            threshold <- sqrt(two_M)
            prop_prob_terms <- mean(term_degrees > threshold)
            prop_prob_elems <- mean(element_degrees > threshold)
            return(list(
              pij_over_1 = pij_over_1,
              prop_problematic_terms = prop_prob_terms,
              prop_problematic_elements = prop_prob_elems,
              summary_pij = summary(rep(dt$pij, dt$count)),
              pij_distribution = dt[count > 0]
            ))
          })
