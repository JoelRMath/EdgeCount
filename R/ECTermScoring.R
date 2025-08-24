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
    
    # Ensure element_set contains valid elements known to the ECTermScoring object
    valid_element_set <- unique(element_set[element_set %in% object@elements])
    if (length(valid_element_set) < 1){
      warning("No valid elements from the input set found in the ECTermScoring object.")
      return(NULL)
    }
    
    # Get terms connected to any element in the valid_element_set
    connected_terms <- get_neighbors(object@ecprob, valid_element_set)
    relevant_terms <- intersect(connected_terms, object@terms)
    
    if (length(relevant_terms) < 1){
      warning("No terms are connected to the provided element set.")
      return(NULL)
    }
    
    # Calculate stats for a single term
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
      
      z_score <- NA_real_
      if (!is.na(lambda) && lambda > 0) {
        z_score <- (observed_edges - lambda) / sqrt(lambda)
      }
      
      log2_Anscombe_ratio <- NA_real_
      if (!is.na(lambda) && (lambda + 3/8) > 0 && (observed_edges + 3/8) > 0) {
        log2_Anscombe_ratio <- 0.5 * (log2(observed_edges + 3/8) - log2(lambda + 3/8))
      }
      
      log2_relative_change <- NA_real_
      if (!is.na(lambda) && lambda > 0 && observed_edges > 0) {
        log2_relative_change <- log2(observed_edges) - log2(lambda)
      }
      
      return(list(p_value = p_value,
                  z_score = z_score,
                  lambda = lambda,
                  observed_edges = observed_edges,
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
          z_score = all_term_scores_list[[term_name]]$z_score,
          lambda = all_term_scores_list[[term_name]]$lambda,
          observed_edges = all_term_scores_list[[term_name]]$observed_edges,
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

#' @title Score Terms Against a Ranked List of Elements
#'
#' @description Calculates a running enrichment score for terms based on a ranked list of elements.
#' @details This method implements a fast algorithm for ranked list analysis, similar in
#' principle to GSEA. For each term, it calculates a profile of statistics (z-score, p-value, etc.)
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
      z_score_term <- (observed_ec_term - lambda_term)/sqrt(lambda_term)
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
                       observed_ec = observed_ec_term,
                       max_ec = max_ec_term,
                       term_size = size_term,
                       z_score = z_score_term,
                       log2_anscombe_ratio = log2_anscombe_ratio_term,
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
#' where these scores occurred.
#'
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
#' # 2. Summarize the results, ranking by the maximum z-score
#' summary_table <- table_terms_ecranks_statistics(
#'   ranked_scores,
#'   scoring_statistic = "z_score",
#'   rank_by = "max"
#' )
#' print(summary_table)
#'
table_terms_ecranks_statistics <- function(term_scores_list,
                                           scoring_statistic = "z_score",
                                           rank_by = "max") {
  
  if (!is.list(term_scores_list)){
    stop("term_scores_list must be a list")
  }
  
  if (length(term_scores_list) == 0) {
    warning("Input 'term_scores_list' is empty.")
    return(data.frame())
  }
  
  allowed_stats = c("lambda", "z_score", "log2_anscombe_ratio", "log2_relative_change", "p_value", "element_relative_rank")
  if (!scoring_statistic %in% allowed_stats){
    warning(paste0("Invalid 'scoring_statistic': ", scoring_statistic, "Allowed choices are: ", toString(allowed_stats)))
    return(data.frame())
  }
  
  required_cols <- c("element", "element_relative_rank", "observed_ec",
                     "lambda", "p_value", "z_score", "log2_anscombe_ratio",
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
