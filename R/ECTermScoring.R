#' @title ECTermScoring S4 Class and Constructor
#'
#' @description The main class for EdgeCount analysis, representing a bipartite
#' graph of term-element memberships. It extends the ECProb class and provides
#' methods for scoring terms against sets or ranked lists of elements.
#'
#' @name ECTermScoring
#' @aliases ECTermScoring-class
#'
#' @slot ecprob An object of class ECProb, representing the underlying
#'   undirected version of the bipartite graph.
#' @slot elements A character vector of all unique vertex IDs identified as elements.
#' @slot terms A character vector of all unique vertex IDs identified as terms.
#'
#' @param term_element_edges A data frame or path to a file. Needs two columns:
#'   one for term IDs and one for element IDs, representing membership.
#'   If a file path, it's assumed to be tab-separated with a header by default.
#' @param col_term The name or index of the column containing term IDs
#'   (in `term_element_edges`). Defaults to 1.
#' @param col_element The name or index of the column containing element IDs
#'   (in `term_element_edges`). Defaults to 2.
#' @param ... Additional arguments passed to `utils::read.table` if
#'   `term_element_edges` is a file path (e.g., `sep`, `header`).
#'
#' @return An object of class ECTermScoring.
#' @seealso
#' Primary analysis functions:
#' \itemize{
#'   \item \code{\link{terms_ecset_statistics}}: Score terms against a single set of elements.
#'   \item \code{\link{terms_ecranks_statistics}}: Score terms against a ranked list of elements.
#'   \item \code{\link{table_terms_ecranks_statistics}}: Summarize and rank results from ranked list scoring.
#' }
#' Underlying graph classes: \code{\link{ECGraph}}, \code{\link{ECProb}}
#'
#' @exportClass ECTermScoring
#' @export ECTermScoring
#' @examples
#' # Create a sample term-element membership data frame
#' te_df <- data.frame(
#'   term = c("TermA", "TermA", "TermB", "TermC", "TermC"),
#'   element = c("Elem1", "Elem2", "Elem2", "Elem3", "Elem4")
#' )
#'
#' # Create an ECTermScoring object
#' ects <- ECTermScoring(te_df)
#' print(ects@terms)
#' print(ects@elements)
#'
setClass("ECTermScoring",
         slots = list(
           ecprob    = "ECProb",
           elements  = "character",
           terms     = "character"
         )
)

ECTermScoring <- function(term_element_edges, col_term = 1, col_element = 2, ...) {

  if (is.character(term_element_edges) && length(term_element_edges) == 1) {
    # File path
    if (!file.exists(term_element_edges)) {
      stop("File not found: ", term_element_edges)
    }

    # Default arguments for read.table
    default_read_args <- list(
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = ""
    )

    # If user-provided arguments from ...
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
    vertex1 = term_ids_char,    # Terms will be 'vertex1'
    vertex2 = element_ids_char, # Elements will be 'vertex2'
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
           function(object, element_set, lambda_method = "optimized")
             standardGeneric("terms_ecset_statistics"))

#' @describeIn terms_ecset_statistics Method for ECTermScoring objects.
setMethod(
  "terms_ecset_statistics",
  "ECTermScoring",
  function(object, element_set, lambda_method = "optimized") {

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

      log2_relative_change <- NA_real_
      if (!is.na(lambda) && lambda > 0 && observed_edges > 0) {
        log2_relative_change <- log2(observed_edges) - log2(lambda)
      }

      return(list(p_value = p_value,
                  observed_edge_count = observed_edge_count,
                  lambda = lambda,
                  log2_Anscombe_ratio = log2_Anscombe_ratio,
                  log2_relative_change = log2_relative_change))
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
          log2_Anscombe = all_term_scores_list[[term_name]]$log2_Anscombe_ratio,
          log2_RelChange = all_term_scores_list[[term_name]]$log2_relative_change,
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
#' @description A high-performance, vectorized function that calculates enrichment
#' statistics for all terms against a collection of input element sets.
#'
#' @details This function is designed for efficiency when testing many input sets
#' (e.g., multiple experimental signatures) against the term database at once.
#' It uses a join-based approach to first identify all relevant `(input_set, term)`
#' pairs and then calculates statistics in a fully vectorized manner.
#'
#' @param object An ECTermScoring object.
#' @param input_sets An object defining the input element sets. This can be either:
#'   \itemize{
#'     \item A `data.table` with two columns: `set_id` and `element`.
#'     \item A named `list` where names are the set IDs and values are character
#'       vectors of elements.
#'   }
#' @param lambda_method The method for lambda calculation. Currently, only "fast"
#'   is supported for this vectorized function.
#'
#' @return A named list of `data.table`s. Each name corresponds to an
#'   `input_set_id`, and each `data.table` contains the enrichment statistics
#'   for all connected terms.
#' @export
#' @examples
#' # Load sample data included with the package
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
              input_sets_dt <- as.data.table(utils::stack(input_sets))
              setnames(input_sets_dt, c("values", "ind"), c("element", "set_id"))
            } else {
              input_sets_dt <- as.data.table(input_sets)
            }

            bipartite_edges <- as.data.table(to_dataframe(object))
            setnames(bipartite_edges, "term", "term_id")
            bipartite_edges[, term_id := as.character(term_id)]

            input_sets_dt_unique <- unique(input_sets_dt, by = c("set_id", "element"))

            setkey(input_sets_dt_unique, element)
            setkey(bipartite_edges, element)
            all_connections <- bipartite_edges[input_sets_dt_unique, on = "element", nomatch = 0, allow.cartesian = TRUE]

            observed_edges_dt <- all_connections[, .(observed_edges = .N), by = .(input_set_id = set_id, term_id)]

            all_element_degrees <- unlist(object@ecprob@degrees)

            term_degrees <- all_element_degrees[bipartite_edges[, unique(term_id)]]
            term_summary <- data.table(term_id = names(term_degrees), term_degree = term_degrees)

            valid_input_elements <- all_connections[, .(input_set_id = set_id, element)] |> unique()
            input_set_summary <- valid_input_elements[,
                                                      .(sum_degrees_set = sum(all_element_degrees[element], na.rm = TRUE)),
                                                      by = input_set_id
            ]

            final_dt <- copy(observed_edges_dt)

            final_dt[term_summary, on = "term_id", term_degree := i.term_degree]
            final_dt[input_set_summary, on = "input_set_id", sum_degrees_set := i.sum_degrees_set]

            input_set_sizes <- valid_input_elements[, .(set_size = .N), by = input_set_id]
            final_dt[input_set_sizes, on = "input_set_id", max_possible_edges := i.set_size]

            final_dt[, lambda := (term_degree * sum_degrees_set) / (2 * object@ecprob@graph_size)]
            final_dt[, p_value := calculate_p_value(object@ecprob, observed_edges, max_possible_edges, lambda)]
            final_dt[, log2_Anscombe_ratio := 0.5 * (log2(observed_edges + 3/8) - log2(lambda + 3/8))]

            setnames(final_dt, "input_set_id", "set1")
            setnames(final_dt, "term_id", "set2")

            results_list <- split(final_dt, by = "set1")

            return(results_list)
          })


#' @title Score Terms Against a Ranked List of Elements
#'
#' @description Calculates a running enrichment score for terms based on a ranked list of elements.
#' @details This method implements a fast algorithm for ranked list analysis, similar in
#' principle to GSEA. For each term, it calculates a profile of statistics (Anscombe-ratio, lambda, p-value, etc.)
#' at each rank position occupied by an element from that term. The time complexity is
#' approximately O(N + sum_k k*log(k)), where N is the total number of ranked elements and k is the
#' size of a term, as it uses an efficient cumulative sum approach for lambda calculation.
#'
#' @param object An ECTermScoring object.
#' @param element_ranks A named list or vector where names are element IDs and values
#'   are their numeric ranks (1 being the highest rank).
#' @param terms A character vector of term IDs to score. If NULL (default),
#'   all terms in the object are processed.
#'
#' @return A named list where each name is a term ID. Each element of the list is a
#'   data frame containing the running enrichment statistics for that term.
#' @export
#' @examples
#' # Create a sample ECTermScoring object
#' te_df <- data.frame(
#'   term = c("TermA", "TermA", "TermB"),
#'   element = c("Elem1", "Elem3", "Elem2")
#' )
#' ects <- ECTermScoring(te_df)
#'
#' # Create a sample ranked list
#' ranks <- c(Elem1 = 1, Elem2 = 2, Elem3 = 3, Elem4 = 4)
#'
#' # Get the ranked scoring profile for TermA
#' ranked_scores <- terms_ecranks_statistics(ects, ranks, terms = "TermA")
#' print(ranked_scores)
#'
setGeneric("terms_ecranks_statistics",
           function(object, element_ranks, terms = NULL)
             standardGeneric("terms_ecranks_statistics"))

#' @describeIn terms_ecranks_statistics Method for ECTermScoring objects.
setMethod(
  "terms_ecranks_statistics",
  "ECTermScoring",
  function(object, element_ranks, terms = NULL) {

    input_elements <- names(element_ranks)
    d_elements <- setdiff(input_elements, object@elements)
    if (length(d_elements) > 0){
      warning(paste("some input elements are not in the ECTErmScoring element universe:",
                    toString(d_elements)))
    }
    input_elements <- input_elements[input_elements %in% object@elements]
    if (length(input_elements) < 1){
      stop("no input element in the ECTErmScoring element universe")
    }

    input_terms <- NULL
    if (is.null(terms)){
      input_terms <- object@terms
    } else {
      d_terms <- setdiff(terms, object@terms)
      if (length(d_terms) > 0){
        warning(paste("some input terms are not in the ECTErmScoring term universe:",
                      toString(d_terms)))
      }
      input_terms <- terms[terms %in% object@terms]
      if (length(input_terms) < 1){
        stop("no input term in the ECTErmScoring term universe")
      }
    }

    ecprob <- object@ecprob
    M_g <- ecprob@graph_size
    N_e <- length(object@elements)

    df <- data.frame(elements = names(element_ranks), ranks = as.numeric(unlist(element_ranks)))
    df <- df[order(df$ranks),]
    ranked_elements <- unlist(df$elements)
    K <- unlist(ecprob@degrees[ranked_elements])
    cumul_sum_K <- cumsum(K)

    score_one_term <- function(obj, term, element_to_ranks, cumul_K, M, N){

      elements_term <- get_neighbors(obj, term)
      sz <- length(elements_term)
      size_term <- rep(sz, sz)
      ranks_term <- unlist(element_to_ranks[elements_term])
      df_term <- data.frame(ranks = ranks_term, elements = elements_term)
      df_term <- df_term[order(df_term$ranks),]
      elements_term <- df_term$elements
      ranks_term <- df_term$ranks
      K_term <- rep(obj@degrees[[term]], length(elements_term))
      one_over_2M <- rep((1/(2*M)), length(elements_term))
      lambda_term <- K_term * one_over_2M * cumul_sum_K[ranks_term]
      observed_ec_term <- 1:length(elements_term)
      max_ec_term <- ranks_term
      log2_anscombe_ratio_term <- 0.5 * (log2(observed_ec_term + 3/8) - log2(lambda_term + 3/8))
      log2_relative_change_term <- log2(observed_ec_term) - log2(lambda_term)
      p_value_term <- mapply(calculate_p_value,
                             z = observed_ec_term,
                             m = max_ec_term,
                             lambda = lambda_term,
                             MoreArgs = list(object = obj),
                             SIMPLIFY = TRUE)
      df <- data.frame(element = elements_term,
                       element_relative_rank = ranks_term/N,
                       lambda = lambda_term,
                       observed_edge_count = observed_ec_term,
                       max_ec = max_ec_term,
                       term_size = size_term,
                       log2_Anscombe_ratio = log2_anscombe_ratio_term,
                       log2_relative_change = log2_relative_change_term,
                       p_value = p_value_term)
      return(df)
    }

    all_results_list <- lapply(input_terms, function(term_id_iter) {
      score_one_term(
        obj = ecprob,
        term = term_id_iter,
        element_to_ranks = element_ranks,
        cumul_K = cumul_sum_K,
        M = M_g,
        N = N_e
      )
    })
    names(all_results_list) <- input_terms

    return(all_results_list)
  })

#' @title Summarize and Rank Term Scoring Profiles
#'
#' @description Processes the list output from \code{\link{terms_ecranks_statistics}}
#' to produce a single, ranked data frame. For each term, it calculates summary
#' statistics (min, median, max) for a chosen score and identifies the elements
#' where these scores occurred. \strong{WARNING: The scores this function produces are
#' not corrected for multiple testing. For a statistically rigorous analysis that includes a
#' Normalized Enrichment Score (NES) and an FDR q-value, use the main
#' analysis function}
#' @param term_scores_list A list of data frames, the output of \code{\link{terms_ecranks_statistics}}.
#' @param scoring_statistic The name of the column to use for calculating summary statistics.
#' @param rank_by A character string specifying which summary statistic to use for
#'   sorting the final data frame. Must be one of "min", "median", or "max".
#'
#' @return A data frame with one row per term, containing a rich summary including
#'   the min, median, and max scores, and the elements/ranks where they were found.
#' @export
#' @examples
#' # Create a sample ECTermScoring object and ranks
#' te_df <- data.frame(
#'   term = c("TermA", "TermA", "TermB"),
#'   element = c("Elem1", "Elem3", "Elem2")
#' )
#' ects <- ECTermScoring(te_df)
#' ranks <- c(Elem1 = 1, Elem2 = 2, Elem3 = 3, Elem4 = 4)
#'
#' # 1. Get the ranked scoring profiles
#' ranked_scores <- terms_ecranks_statistics(ects, ranks)
#'
#' # 2. Summarize the results, ranking by the maximum log2_Anscombe_ratio
#' summary_table <- table_terms_ecranks_statistics(
#'   ranked_scores,
#'   scoring_statistic = "log2_Anscombe_ratio",
#'   rank_by = "max"
#' )
#' print(summary_table)
#'
table_terms_ecranks_statistics <- function(term_scores_list,
                                           scoring_statistic = "log2_Anscombe_ratio",
                                           rank_by = "median") {

  if (!is.list(term_scores_list)){
    stop("term_scores_list must be a list")
  }

  if (length(term_scores_list) == 0) {
    warning("Input 'term_scores_list' is empty.")
    return(data.frame())
  }

  allowed_stats = c("lambda", "observed_edge_count", "log2_Anscombe_ratio", "log2_relative_change", "p_value", "element_relative_rank")
  if (!scoring_statistic %in% allowed_stats){
    warning(paste0("Invalid 'scoring_statistic': ", scoring_statistic, "Allowed choices are: ", toString(allowed_stats)))
    return(data.frame())
  }

  required_cols <- c("element", "element_relative_rank", "observed_edge_count",
                     "lambda", "p_value", "log2_Anscombe_ratio",
                     "log2_relative_change")

  first_valid_df <- NULL
  for(item in term_scores_list){
    if(is.data.frame(item) && nrow(item) > 0){
      first_valid_df <- item
      break
    }
  }

  missing_cols <- setdiff(required_cols, names(first_valid_df))
  if (length(missing_cols) > 0) {
    warning(paste0("Data frame in 'term_scores_list' is missing required columns: ",
                   toString(missing_cols)))
  }

  if(is.null(first_valid_df) || !(scoring_statistic %in% names(first_valid_df))){
    stop(paste0("Specified 'scoring_statistic': '", scoring_statistic,
                "' not found in the result data frames. Valid columns include: ",
                toString(intersect(allowed_stats, names(first_valid_df)))))
  }

  get_term_row <- function(term_id){

    term_df <- term_scores_list[[term_id]]

    na_row <- data.frame(
      term_id = term_id,
      term_size = NA_integer_,
      min_score = NA_real_, element_at_min = NA_character_, rank_at_min = NA_real_,
      median_score = NA_real_, element_at_median = NA_character_, rank_at_median = NA_real_,
      max_score = NA_real_, element_at_max = NA_character_, rank_at_max = NA_real_,
      stringsAsFactors = FALSE
    )

    if (!is.data.frame(term_df) || nrow(term_df) == 0) {
      return(na_row)
    }

    if (!(scoring_statistic %in% names(term_df))) {
      warning(paste("Column '", scoring_statistic, "' not found for term:", term_id, ". Skipping term."))
      return(na_row)
    }

    scores_vec <- as.numeric(term_df[[scoring_statistic]])
    if (all(is.na(scores_vec))) {
      return(na_row)
    }

    is_valid_score <- is.finite(scores_vec)
    finite_scores <- scores_vec[is_valid_score]

    min_val <- min(finite_scores, na.rm = TRUE)
    median_val <- stats::median(finite_scores, na.rm = TRUE)
    max_val <- max(finite_scores, na.rm = TRUE)

    min_idx <- which(scores_vec == min_val)[1]
    max_idx <- which(scores_vec == max_val)[1]
    median_idx <- which.min(abs(scores_vec - median_val))[1]

    summary_row <- data.frame(
      term_id = term_id,
      term_size = term_df$term_size[1],
      min_score = min_val,
      element_at_min = term_df$element[min_idx],
      rank_at_min = term_df$element_relative_rank[min_idx],
      median_score = median_val,
      element_at_median = term_df$element[median_idx],
      rank_at_median = term_df$element_relative_rank[median_idx],
      max_score = max_val,
      element_at_max = term_df$element[max_idx],
      rank_at_max = term_df$element_relative_rank[max_idx],
      stringsAsFactors = FALSE
    )

    return(summary_row)
  }

  summary_list <- lapply(names(term_scores_list), get_term_row)
  summary_df <- dplyr::bind_rows(summary_list)

  if (nrow(summary_df) > 0) {
    sort_col <- switch(rank_by,
                       "max" = "max_score",
                       "min" = "min_score",
                       "median" = "median_score",
                       "max_score") # Default to max_score if invalid input

    decreasing_order <- TRUE
    if (rank_by == "min" || (scoring_statistic == "p_value" && rank_by != "max")) {
      decreasing_order <- FALSE
    }

    if (any(!is.na(summary_df[[sort_col]]))) {
      summary_df <- summary_df[order(summary_df[[sort_col]], decreasing = decreasing_order), ]
    }

    rownames(summary_df) <- NULL
  }

  return(summary_df)
}

#' @title Perform Full VSEA-style Permutation Analysis for Multiple Metrics
#'
#' @description Runs a complete ranked-list enrichment analysis, calculating a
#' Normalized Enrichment Score (NES) and a False Discovery Rate (FDR) q-value
#' for each term using a permutation-based null distribution. This function
#' efficiently calculates results for the **max**, **min**, and **median** summary
#' scores from a single set of permutations.
#'
#' @param ects_object An ECTermScoring object.
#' @param element_ranks A named list or vector of element ranks.
#' @param scoring_statistic The score to use for enrichment (e.g., "log2_Anscombe_ratio").
#' @param n_permutations The number of permutations to generate the null
#'   distribution. 1000 is a common choice for publication.
#' @param seed An optional random seed for reproducibility.
#'
#' @return A named list of three data frames: `max_score_summary`,
#'   `min_score_summary`, and `median_score_summary`. Each data frame is sorted
#'   by its respective NES and contains a rich summary including the raw score,
#'   NES, and FDR q-value.
#'
#' @references
#' 1. Subramanian, A., et al. (2005). Gene set enrichment analysis: a knowledge-based
#' approach for interpreting genome-wide expression profiles. PNAS.
#' \url{https://www.pnas.org/doi/10.1073/pnas.0506580102}
#'
#' #' 2. Pradines et al. (2005). Analyzing Protein Lists with Large Networks:
#' Edge-Count Probabilities in Random Graphs with Given Expected Degrees.
#' J. Comp. Biol. 12(2):113-28.
#' \url{https://www.liebertpub.com/doi/10.1089/cmb.2005.12.113}
#'
#' @export
#' @examples
#' \dontrun{
#' # This is a computationally intensive function.
#' te_df <- data.frame(
#'   term = c("TermA", "TermA", "TermB"),
#'   element = c("Elem1", "Elem3", "Elem2")
#' )
#' ects <- ECTermScoring(te_df)
#' ranks <- c(Elem1 = 1, Elem2 = 2, Elem3 = 3, Elem4 = 4)
#'
#' # Run the full analysis for all three summary metrics
#' vsea_results <- run_vsea_analysis(
#'   ects,
#'   ranks,
#'   scoring_statistic = "log2_Anscombe_ratio",
#'   n_permutations = 100 # Use a smaller number for quick tests
#' )
#'
#' # View the top results for enrichment at the top of the list
#' print(head(vsea_results$max_score_summary))
#'
#' # View the top results for enrichment at the bottom of the list
#' print(head(vsea_results$min_score_summary))
#' }
#'
run_vsea_analysis <- function(ects_object,
                              element_ranks,
                              scoring_statistic = "log2_Anscombe_ratio",
                              n_permutations = 1000,
                              seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Step 1: Calculate the REAL scores for the actual ranked list ---
  message("Calculating real enrichment scores...")
  real_scores_list <- terms_ecranks_statistics(ects_object, element_ranks)

  # Get all three summary tables up front
  real_summary_max    <- table_terms_ecranks_statistics(real_scores_list, scoring_statistic, "max")
  real_summary_min    <- table_terms_ecranks_statistics(real_scores_list, scoring_statistic, "min")
  real_summary_median <- table_terms_ecranks_statistics(real_scores_list, scoring_statistic, "median")

  # Extract the real scores into named vectors
  real_scores_max    <- stats::setNames(real_summary_max$max_score, real_summary_max$term_id)
  real_scores_min    <- stats::setNames(real_summary_min$min_score, real_summary_min$term_id)
  real_scores_median <- stats::setNames(real_summary_median$median_score, real_summary_median$term_id)

  # --- Step 2: Build the NULL distribution via permutation ---
  message(paste("Running", n_permutations, "permutations..."))

  # Helper function
  get_all_summary_scores <- function(term_df, stat) {
    if (is.null(term_df) || nrow(term_df) == 0) {
      return(c(max_score = NA_real_, min_score = NA_real_, median_score = NA_real_))
    }
    scores <- term_df[[stat]]
    scores <- scores[is.finite(scores)]
    if (length(scores) == 0) {
      return(c(max_score = NA_real_, min_score = NA_real_, median_score = NA_real_))
    }
    return(c(
      max_score = max(scores, na.rm = TRUE),
      min_score = min(scores, na.rm = TRUE),
      median_score = stats::median(scores, na.rm = TRUE)
    ))
  }

  # List where each element is a matrix from a single permutation
  # (rows = terms, cols = max/min/median)
  perm_results_list <- replicate(n_permutations, {
    permuted_ranks <- sample(element_ranks)
    names(permuted_ranks) <- names(element_ranks)
    permuted_scores_list <- terms_ecranks_statistics(ects_object, permuted_ranks)
    t(sapply(permuted_scores_list, get_all_summary_scores, stat = scoring_statistic))
  }, simplify = FALSE)

  # Combine the list of matrices into three separate null score matrices
  term_names <- names(real_scores_list)
  null_scores_max    <- do.call(cbind, lapply(perm_results_list, function(m) m[term_names, "max_score"]))
  null_scores_min    <- do.call(cbind, lapply(perm_results_list, function(m) m[term_names, "min_score"]))
  null_scores_median <- do.call(cbind, lapply(perm_results_list, function(m) m[term_names, "median_score"]))


  # --- Step 3 & 4: Calculate NES and FDR for each metric ---
  message("Calculating NES and FDR for each summary metric...")

  calculate_nes_fdr <- function(real_scores, null_scores_matrix) {
    # Normalize real scores
    mean_null_pos <- apply(null_scores_matrix, 1, function(x) mean(x[x > 0], na.rm = TRUE))
    mean_null_neg <- apply(null_scores_matrix, 1, function(x) mean(abs(x[x < 0]), na.rm = TRUE))
    mean_null_pos[is.nan(mean_null_pos)] <- 1
    mean_null_neg[is.nan(mean_null_neg)] <- 1

    nes_real <- ifelse(
      real_scores > 0,
      real_scores / mean_null_pos[names(real_scores)],
      real_scores / mean_null_neg[names(real_scores)]
    )

    # Normalize null scores
    null_nes_matrix <- ifelse(
      null_scores_matrix > 0,
      null_scores_matrix / mean_null_pos,
      null_scores_matrix / mean_null_neg
    )
    all_null_nes <- as.vector(null_nes_matrix)
    all_null_nes <- all_null_nes[is.finite(all_null_nes)]
    null_nes_pos <- all_null_nes[all_null_nes > 0]
    null_nes_neg <- all_null_nes[all_null_nes < 0]

    # Calculate FDR
    fdr_pos <- vapply(nes_real[nes_real > 0], function(score) {
      if (length(null_nes_pos) == 0) return(1)
      frac_null <- sum(null_nes_pos >= score) / length(null_nes_pos)
      frac_real <- sum(nes_real > 0 & nes_real >= score) / sum(nes_real > 0)
      if (frac_real == 0) return(1)
      return(frac_null / frac_real)
    }, FUN.VALUE = numeric(1))

    fdr_neg <- vapply(nes_real[nes_real < 0], function(score) {
      if (length(null_nes_neg) == 0) return(1)
      frac_null <- sum(null_nes_neg <= score) / length(null_nes_neg)
      frac_real <- sum(nes_real < 0 & nes_real <= score) / sum(nes_real < 0)
      if (frac_real == 0) return(1)
      return(frac_null / frac_real)
    }, FUN.VALUE = numeric(1))

    fdr_q_values <- c(fdr_pos, fdr_neg)
    fdr_q_values[fdr_q_values > 1] <- 1

    return(list(nes = nes_real, fdr_q_value = fdr_q_values))
  }

  results_max    <- calculate_nes_fdr(real_scores_max, null_scores_max)
  results_min    <- calculate_nes_fdr(real_scores_min, null_scores_min)
  results_median <- calculate_nes_fdr(real_scores_median, null_scores_median)

  # --- Step 5: Combine, sort, and return ---
  add_results_and_sort <- function(summary_df, results) {
    summary_df$nes <- results$nes[summary_df$term_id]
    summary_df$fdr_q_value <- results$fdr_q_value[summary_df$term_id]
    summary_df <- summary_df[order(summary_df$fdr_q_value, decreasing = FALSE), ]
    rownames(summary_df) <- NULL
    return(summary_df)
  }

  final_max    <- add_results_and_sort(real_summary_max, results_max)
  final_min    <- add_results_and_sort(real_summary_min, results_min)
  final_median <- add_results_and_sort(real_summary_median, results_median)

  return(list(
    max_score_summary = final_max,
    min_score_summary = final_min,
    median_score_summary = final_median
  ))
}

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

  # --- STEP 1: Initialize mutable degree vectors for terms and elements ---
  all_degrees <- unlist(object@ecprob@degrees)
  current_term_degrees <- all_degrees[object@terms]
  current_element_degrees <- all_degrees[object@elements]

  two_m_current <- object@ecprob@graph_size * 2
  removed_terms <- character(0)

  # --- STEP 2: Iteratively remove terms and update element degrees ---
  while (length(current_term_degrees) > 0 && length(current_element_degrees) > 0) {
    max_kt <- max(current_term_degrees)
    max_ke <- max(current_element_degrees)
    # print(max_ke * max_kt / two_m_current)

    # If the condition is met, we're done trimming.
    if (max_ke * max_kt < threshold * two_m_current) {
      break
    }

    # Identify the term(s) with the current max degree
    terms_to_remove <- names(current_term_degrees[current_term_degrees == max_kt])
    removed_terms <- c(removed_terms, terms_to_remove)

    # Find neighbors of the removed terms to update their degrees
    neighbors_to_update <- unlist(object@ecprob@adj[terms_to_remove])

    # Efficiently decrement the degrees of affected elements
    if(length(neighbors_to_update) > 0) {
      neighbor_counts <- table(neighbors_to_update)
      affected_elements <- names(neighbor_counts)
      current_element_degrees[affected_elements] <- current_element_degrees[affected_elements] - neighbor_counts
    }

    # Update state for the next iteration
    two_m_current <- two_m_current - (length(terms_to_remove) * max_kt * 2) # Each edge removal reduces sum of degrees by 2
    current_term_degrees <- current_term_degrees[!names(current_term_degrees) %in% terms_to_remove]
    current_element_degrees <- current_element_degrees[current_element_degrees > 0] # Remove elements that are now disconnected
  }

  # --- STEP 3: Rebuild the new, trimmed ECTermScoring object ---
  kept_terms <- setdiff(object@terms, removed_terms)

  if (length(kept_terms) == 0) {
    # Return an empty object if no terms are left
    empty_df <- data.frame(term=character(), element=character())
    final_removed_elements <- object@elements
    return(list(
      trimmed_object = ECTermScoring(empty_df),
      removed_terms = removed_terms,
      removed_elements = final_removed_elements
    ))
  }

  # Reconstruct the term->element edge list from the original object
  term_neighbors_list <- object@ecprob@adj[kept_terms]
  new_edge_df <- utils::stack(term_neighbors_list)
  names(new_edge_df) <- c("element", "term")

  # Create the new ECTermScoring object using its constructor.
  # The constructor will automatically determine the final set of elements.
  new_ects <- ECTermScoring(new_edge_df[, c("term", "element")])

  # Determine which elements were removed as a consequence
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
    # No empty terms to remove, return original object
    return(object)
  }

  # Reconstruct the term->element edge list from the original object,
  # but only for the terms we are keeping.
  term_neighbors_list <- object@ecprob@adj[kept_terms]

  new_edge_df <- utils::stack(term_neighbors_list)

  if(nrow(new_edge_df) == 0) {
    # Handle case where removing terms leaves no edges
    final_df <- data.frame(term=character(), element=character())
  } else {
    names(new_edge_df) <- c("element", "term")
    final_df <- new_edge_df[, c("term", "element")]
  }

  # Create a new ECTermScoring object using its constructor.
  # The constructor will correctly build the new graph and element list.
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
            # stack() correctly creates the unique term -> element edge list
            edge_df <- utils::stack(object@ecprob@adj[object@terms])

            if(nrow(edge_df) == 0) {
              return(data.frame(term = character(), element = character()))
            }

            names(edge_df) <- c("element", "term")

            edge_df[, c("term", "element")]
          })
