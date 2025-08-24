#' @title ECGraph S4 Class and Constructor
#'
#' @description Represents an undirected graph and provides a constructor to create it.
#' The ECGraph object stores the graph's adjacency list, vertex degrees, and names.
#' The constructor accepts a file path, a data frame, or an igraph object.
#'
#' @name ECGraph
#' @aliases ECGraph-class
#'
#' @slot adj list. An adjacency list where names are vertex IDs and values are
#'   character vectors of neighboring vertex IDs.
#' @slot degrees list. A named list where names are vertex IDs and values are
#'   the integer degrees of vertices.
#' @slot names character. A character vector of all unique vertex names in the graph.
#'
#' @param edges_input A character string specifying the path to a file,
#'   a data frame with two columns representing edges, or an \code{igraph} graph object.
#'   Vertex IDs will be coerced to character.
#' @param col1 The name or index of the first column if `edges_input` is a
#'   data frame or file. Defaults to 1.
#' @param col2 The name or index of the second column if `edges_input` is a
#'   data frame or file. Defaults to 2.
#' @param header A logical value for reading from a file. Defaults to TRUE.
#' @param sep The field separator character for reading from a file. Defaults to "\t".
#'
#' @return An object of class ECGraph.
#' @importFrom igraph is_igraph as_data_frame as.undirected
#' @seealso
#' Functions that operate on ECGraph objects:
#' \itemize{
#'   \item \code{\link{to_dataframe}}: Convert an ECGraph to a data frame edge list.
#'   \item \code{\link{get_neighbors}}: Retrieve all neighbors for a set of vertices.
#'   \item \code{\link{get_edge_count_in}}: Count edges within a set of vertices.
#'   \item \code{\link{get_edge_count_between}}: Count edges between two sets of vertices.
#'   \item \code{\link{get_edge_count_in_max}}: Maximum possible edge count within a set (complete graph).
#'   \item \code{\link{get_edge_count_between_max}}: Maximum possible edge count between two sets (complete graph).
#'   \item \code{\link{get_edge_count_in_max_fds}}: Maximum possible edge count within a set (fixed degree sequence).
#'   \item \code{\link{get_edge_count_between_max_fds}}: Maximum possible edge count between two sets (fixed degree sequence).
#' }
#'
#' @exportClass ECGraph
#' @export ECGraph
#' @examples
#' # --- Data ---
#' edge_df <- data.frame(
#'   vertex1 = c("A", "A", "B", "C", "E"), # E is an isolated vertex
#'   vertex2 = c("B", "C", "C", "D", "E"), # E is connected only to itself (self-loop)
#'   stringsAsFactors = FALSE
#' )
#' temp_file <- tempfile(fileext = ".tsv")
#' write.table(edge_df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' # --- ECGraph from a data frame ---
#' graph_from_df <- ECGraph(edge_df)
#' print("Adjacency list from data frame:")
#' print(graph_from_df@adj)
#'
#' # --- ECGraph from a file path ---
#' graph_from_file <- ECGraph(temp_file)
#'
#' # --- ECGraph from an igraph object ---
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   ig <- igraph::graph_from_data_frame(edge_df, directed = FALSE)
#'   graph_from_ig <- ECGraph(ig)
#' }
#'
#' unlink(temp_file)
#'
setClass("ECGraph",
         slots = list(
           adj = "list",
           degrees = "list",
           names = "character"
         )
)

ECGraph <- function(edges_input,
                    col1 = 1,
                    col2 = 2,
                    header = TRUE,
                    sep = "\t") {
  
  if (is.character(edges_input) && length(edges_input) == 1) {
    # file path
    if (!file.exists(edges_input)) {
      stop("File not found: ", edges_input)
    }
    tryCatch({
      df <- utils::read.table(edges_input,
                              header = header,
                              sep = sep,
                              stringsAsFactors = FALSE,
                              quote = "")
    }, error = function(e) {
      stop("Error reading file '", edges_input, "': ", e$message)
    })
    if (ncol(df) < 2) {
      stop("Input file '", edges_input, "' must have at least two columns.")
    }
  } else if (is.data.frame(edges_input)) {
    # data frame
    df <- edges_input
  } else if (requireNamespace("igraph", quietly = TRUE) && igraph::is_igraph(edges_input)) {
    # igraph object
    df <- igraph::as_data_frame(igraph::as.undirected(edges_input), what = "edges")
    col1 <- "from"
    col2 <- "to"
  } else {
    stop("'edges_input' must be a file path, a data frame, or an igraph object.")
  }
  
  # input checks
  if (ncol(df) < 2) {
    stop("Input data frame must have at least two columns.")
  }
  if (is.numeric(col1) && (col1 < 1 || col1 > ncol(df))) stop("col1 index is out of bounds for the data frame.")
  if (is.numeric(col2) && (col2 < 1 || col2 > ncol(df))) stop("col2 index is out of bounds for the data frame.")
  if (is.character(col1) && !(col1 %in% names(df))) stop("col1 '", col1, "' not found in data frame names.")
  if (is.character(col2) && !(col2 %in% names(df))) stop("col2 '", col2, "' not found in data frame names.")
  
  node_a_char <- as.character(df[[col1]])
  node_b_char <- as.character(df[[col2]])
  
  edges_symmetric <- data.frame(
    w1 = c(node_a_char, node_b_char),
    w2 = c(node_b_char, node_a_char),
    stringsAsFactors = FALSE
  )
  
  # no self-loops
  edges_symmetric <- edges_symmetric[edges_symmetric$w1 != edges_symmetric$w2, ]
  
  if (nrow(edges_symmetric) == 0) {
    adj_list <- list()
    degrees_list <- list()
    all_vertex_names <- character(0)
  } else {
    # adjacency list
    adj_list <- split(edges_symmetric$w2, factor(edges_symmetric$w1))
    
    # neighbors are unique for each vertex and sorted
    adj_list <- lapply(adj_list, function(x) sort(unique(x)))
    
    # degrees
    all_vertex_names <- names(adj_list)
    degrees_list <- lapply(adj_list, length)
    names(degrees_list) <- all_vertex_names
    
    all_mentioned_nodes <- unique(c(node_a_char, node_b_char))
    all_mentioned_nodes <- all_mentioned_nodes[all_mentioned_nodes != ""]
    
    # initializes degrees
    temp_degrees <- stats::setNames(rep(0, length(all_mentioned_nodes)), all_mentioned_nodes)
    
    # updates with actual degrees
    if(length(degrees_list) > 0) {
      temp_degrees[names(degrees_list)] <- unlist(degrees_list)
    }
    degrees_list <- as.list(temp_degrees)
    
    # updates adj_list to include empty neighbor lists for isolated nodes
    for (node_name in all_mentioned_nodes) {
      if (is.null(adj_list[[node_name]])) {
        adj_list[[node_name]] <- character(0)
      }
    }
    # re-sorts adj_list by name
    adj_list <- adj_list[sort(names(adj_list))]
    degrees_list <- degrees_list[sort(names(degrees_list))]
    all_vertex_names <- sort(names(degrees_list))
  }
  
  new("ECGraph", adj = adj_list, degrees = degrees_list, names = all_vertex_names)
}

#' @title Convert ECGraph to Data Frame
#'
#' @description Creates a data frame of edges from an ECGraph object.
#' The returned data frame has two columns representing the undirected
#' edges contained in the ECGraph object. Isolated vertices are not included.
#'
#' @param object An ECGraph object.
#'
#' @return A data frame with two columns "from" and "to" representing the edge set.
#' @export
#' @examples
#' edge_df <- data.frame(vertex1 = c("A", "B"), vertex2 = c("B", "C"))
#' ecg <- ECGraph(edge_df)
#' edge_list_df <- to_dataframe(ecg)
#' print(edge_list_df)
#'
setGeneric("to_dataframe",
           function(object) standardGeneric("to_dataframe"))

#' @describeIn to_dataframe Method for ECGraph objects.
setMethod("to_dataframe",
          "ECGraph",
          function(object) {
            edge_df <- stack(object@adj)
            names(edge_df) <- c("to", "from")
            # Return a consistent column order
            edge_df[, c("from", "to")]
          })

#' @title Get Neighbors of Vertices
#'
#' @description Retrieves a vector of all unique neighbors for a given
#' set of vertices.
#'
#' @param object An ECGraph object.
#' @param vertices A character vector of vertex names.
#'
#' @return A character vector of unique neighbor names.
#' @export
#' @examples
#' edge_df <- data.frame(v1 = c("A", "B", "C"), v2 = c("B", "C", "D"))
#' ecg <- ECGraph(edge_df)
#' # Get neighbors of vertices A and C
#' get_neighbors(ecg, c("A", "C")) # Returns "B" and "D"
#'
setGeneric("get_neighbors",
           function(object, vertices) standardGeneric("get_neighbors"))

#' @describeIn get_neighbors Method for ECGraph objects.
setMethod("get_neighbors",
          "ECGraph",
          function(object, vertices) {
            # Filter for vertices that exist in the graph
            valid_vertices <- vertices[vertices %in% object@names]
            if (length(valid_vertices) == 0) return(character(0))
            
            neighbors_list <- object@adj[valid_vertices]
            all_neighbors <- unlist(neighbors_list)
            unique(all_neighbors)
          })

#' @title Get Edge Count Within a Set of Vertices
#'
#' @description Calculates the number of edges connecting vertices *within* a given set.
#'
#' @param object An ECGraph object.
#' @param vertices A character vector of vertex names.
#'
#' @return An integer count of the edges within the set.
#' @export
#' @examples
#' edge_df <- data.frame(p1 = c("A", "B", "C", "D"), p2 = c("B", "C", "A", "E"))
#' ecg <- ECGraph(edge_df)
#' # Count edges within the set {A, B, C} (A-B, B-C, C-A)
#' get_edge_count_in(ecg, c("A", "B", "C")) # Returns 3
#'
setGeneric("get_edge_count_in",
           function(object, vertices) standardGeneric("get_edge_count_in"))

#' @describeIn get_edge_count_in Method for ECGraph objects.
setMethod("get_edge_count_in",
          "ECGraph",
          function(object, vertices) {
            # Filter for vertices that exist in the graph
            valid_vertices <- vertices[vertices %in% object@names]
            if (length(valid_vertices) == 0) return(0)
            
            adj_subset <- object@adj[valid_vertices]
            edge_count <- sum(unlist(lapply(adj_subset,
                                            function(neighbors) {
                                              sum(neighbors %in% valid_vertices)
                                            })))
            # Each edge is counted twice (once from each direction), so divide by 2
            edge_count / 2
          })


#' @title Get Edge Count Between Two Sets of Vertices
#'
#' @description Calculates the number of edges connecting two specified sets of vertices.
#' Edges within the intersection of the sets are not counted.
#'
#' @param object An ECGraph object.
#' @param set1 A character vector of vertex names for the first set.
#' @param set2 A character vector of vertex names for the second set.
#'
#' @return An integer count of the edges between the two sets.
#' @export
#' @examples
#' edge_df <- data.frame(p1 = c("A", "B", "C"), p2 = c("X", "Y", "A"))
#' ecg <- ECGraph(edge_df)
#' get_edge_count_between(ecg, c("A", "B"), c("X", "Y")) # Returns 2 (A-X, B-Y)
#' get_edge_count_between(ecg, c("A", "B"), c("C", "X")) # Returns 2 (A-X, A-C)
#'
setGeneric("get_edge_count_between",
           function(object, set1, set2) standardGeneric("get_edge_count_between"))

#' @describeIn get_edge_count_between Method for ECGraph objects.
setMethod("get_edge_count_between",
          "ECGraph",
          function(object, set1, set2) {
            # Filter for vertices that exist in the graph
            valid_set1 <- set1[set1 %in% object@names]
            valid_set2 <- set2[set2 %in% object@names]
            
            set1_d <- setdiff(valid_set1, valid_set2)
            set2_d <- setdiff(valid_set2, valid_set1)
            
            if (length(set1_d) == 0 || length(set2_d) == 0) return(0)
            
            adj_subset <- object@adj[set1_d]
            edge_count <- sum(unlist(lapply(adj_subset,
                                            function(neighbors) {
                                              sum(neighbors %in% set2_d)
                                            })))
            edge_count
          })

#' @title Get Maximum Possible Edge Count Within a Set
#'
#' @description Calculates the maximum possible number of edges in a simple graph
#' formed by a given set of vertices (i.e., a complete graph). Only vertices
#' present in the ECGraph object are considered.
#'
#' @param object An ECGraph object.
#' @param vertices A character vector of vertex names.
#'
#' @return The maximum possible number of edges.
#' @export
#' @examples
#' ecg <- ECGraph(data.frame(p1=c("A","B","C","D"), p2=c("B","C","D","A")))
#' # The function ignores vertices not in the graph, like "X"
#' # It calculates max edges for the valid set {A, B, C}, which is 3.
#' get_edge_count_in_max(ecg, c("A", "B", "C", "X"))
#'
setGeneric("get_edge_count_in_max",
           function(object, vertices) standardGeneric("get_edge_count_in_max"))

#' @describeIn get_edge_count_in_max Method for ECGraph objects.
setMethod("get_edge_count_in_max",
          "ECGraph",
          function(object, vertices) {
            # FIX: Only consider vertices that are actually in the graph
            valid_vertices <- vertices[vertices %in% object@names]
            n <- length(valid_vertices)
            if (n <= 1) {
              return(0)
            } else {
              return(n * (n - 1) / 2)
            }
          })

#' @title Get Maximum Possible Edge Count Between Two Sets
#'
#' @description Calculates the maximum possible number of edges between two disjoint
#' sets of vertices (i.e., a complete bipartite graph). Only vertices
#' present in the ECGraph object are considered.
#'
#' @param object An ECGraph object.
#' @param set1 A character vector of vertex names.
#' @param set2 A character vector of vertex names.
#'
#' @return The maximum possible number of edges between the sets.
#' @export
#' @examples
#' ecg <- ECGraph(data.frame(p1=c("A","B","X","Y"), p2=c("B","A","Y","X")))
#' # "Z" is ignored as it's not in the graph.
#' # The result is length(c("A","B")) * length(c("X","Y")) = 2 * 2 = 4.
#' get_edge_count_between_max(ecg, c("A", "B"), c("X", "Y", "Z"))
#'
setGeneric("get_edge_count_between_max",
           function(object, set1, set2) standardGeneric("get_edge_count_between_max"))

#' @describeIn get_edge_count_between_max Method for ECGraph objects.
setMethod("get_edge_count_between_max",
          "ECGraph",
          function(object, set1, set2) {
            # FIX: Only consider vertices that are actually in the graph
            valid_set1 <- set1[set1 %in% object@names]
            valid_set2 <- set2[set2 %in% object@names]
            length(valid_set1) * length(valid_set2)
          })

#' @title Get Max Edge Count Within a Set (Fixed Degree Sequence)
#'
#' @description Calculates the maximum possible number of edges within a set of
#' vertices, given their prescribed degrees (Fixed Degree Sequence model).
#'
#' @param object An ECGraph object.
#' @param vertices A character vector of vertex names.
#'
#' @return The maximum possible number of edges, while preserving vertex degrees.
#' @export
#' @examples
#' # Example from a known graphic sequence
#' ecg <- ECGraph(data.frame(p1=c("A","B","C"), p2=c("B","C","A"))) # All degrees are 2
#' get_edge_count_in_max_fds(ecg, c("A", "B", "C")) # Returns 3
#'
setGeneric("get_edge_count_in_max_fds",
           function(object, vertices) standardGeneric("get_edge_count_in_max_fds"))

#' @describeIn get_edge_count_in_max_fds Method for ECGraph objects.
setMethod("get_edge_count_in_max_fds",
          "ECGraph",
          function(object, vertices) {
            # Filter for valid vertices and get their degrees
            valid_vertices <- vertices[vertices %in% object@names]
            if (length(valid_vertices) < 2) return(0)
            
            R <- unlist(object@degrees[valid_vertices])
            w <- 0
            while (TRUE) {
              R <- sort(R, decreasing = TRUE)
              R <- R[R >= 1]
              m <- length(R)
              R <- pmin(R, m - 1)
              if (length(R) > 0) {
                w <- w + R[1]
              } else {
                return(as.numeric(w))
              }
              r1 <- R[1]
              R <- c(head(R[-1], r1) - 1, tail(R[-1], max(0, length(R) - 1 - r1)))
              R <- R[R > 0]
              if (length(R) < 2) {
                return(as.numeric(w))
              }
            }
            return(as.numeric(w))
          })

#' @title Get Max Edge Count Between Sets (Fixed Degree Sequence)
#'
#' @description Calculates the maximum possible number of edges between two disjoint
#' sets of vertices, given their prescribed degrees (Fixed Degree Sequence model).
#'
#' @param object An ECGraph object.
#' @param set1 A character vector of vertex names.
#' @param set2 A character vector of vertex names.
#'
#' @return The maximum possible number of edges, while preserving vertex degrees.
#' @export
#' @examples
#' # Bipartite graph example
#' edge_df <- data.frame(p1=c("A","A","B"), p2=c("X","Y","Y"))
#' ecg <- ECGraph(edge_df)
#' get_edge_count_between_max_fds(ecg, c("A","B"), c("X","Y")) # Returns 3
#'
setGeneric("get_edge_count_between_max_fds",
           function(object, set1, set2) standardGeneric("get_edge_count_between_max_fds"))

#' @describeIn get_edge_count_between_max_fds Method for ECGraph objects.
setMethod("get_edge_count_between_max_fds",
          "ECGraph",
          function(object, set1, set2) {
            # Filter for valid vertices
            valid_set1 <- set1[set1 %in% object@names]
            valid_set2 <- set2[set2 %in% object@names]
            
            set1_disjoint <- setdiff(valid_set1, valid_set2)
            set2_disjoint <- setdiff(valid_set2, valid_set1)
            
            if (length(set1_disjoint) == 0 || length(set2_disjoint) == 0) {
              return(0)
            }
            
            R1_iter <- unlist(object@degrees[set1_disjoint])
            R2_iter <- unlist(object@degrees[set2_disjoint])
            w <- 0
            while (TRUE) {
              R1_iter <- sort(R1_iter, decreasing = TRUE)
              R2_iter <- sort(R2_iter, decreasing = TRUE)
              R1_iter <- R1_iter[R1_iter >= 1]
              R2_iter <- R2_iter[R2_iter >= 1]
              m1_iter <- length(R1_iter)
              m2_iter <- length(R2_iter)
              if (m1_iter == 0 || m2_iter == 0) {
                return(as.numeric(w))
              }
              R1_for_selection <- pmin(R1_iter, m2_iter)
              R2_for_selection <- pmin(R2_iter, m1_iter)
              R1_for_selection <- R1_for_selection[R1_for_selection >= 1]
              R2_for_selection <- R2_for_selection[R2_for_selection >= 1]
              if (length(R1_for_selection) == 0 && length(R2_for_selection) == 0) {
                return(as.numeric(w))
              }
              s_indicates_R1_is_Rs <- TRUE
              if (length(R1_for_selection) == 0) {
                if (length(R2_for_selection) > 0) s_indicates_R1_is_Rs <- FALSE else return(as.numeric(w))
              } else if (length(R2_for_selection) > 0) {
                if (R1_for_selection[1] < R2_for_selection[1]) {
                  s_indicates_R1_is_Rs <- FALSE
                }
              }
              r_s_1_connections <- 0
              if (s_indicates_R1_is_Rs) {
                r_s_1_connections <- R1_for_selection[1]
              } else {
                r_s_1_connections <- R2_for_selection[1]
              }
              if (r_s_1_connections < 1) {
                return(as.numeric(w))
              }
              w <- w + r_s_1_connections
              if (s_indicates_R1_is_Rs) {
                R1_iter <- R1_iter[-1]
                if (length(R2_iter) > 0 && r_s_1_connections > 0) {
                  num_to_decrement <- min(r_s_1_connections, length(R2_iter))
                  decremented_part <- head(R2_iter, num_to_decrement) - 1
                  tail_part <- tail(R2_iter, max(0, length(R2_iter) - num_to_decrement))
                  R2_iter <- c(decremented_part, tail_part)
                }
              } else {
                R2_iter <- R2_iter[-1]
                if (length(R1_iter) > 0 && r_s_1_connections > 0) {
                  num_to_decrement <- min(r_s_1_connections, length(R1_iter))
                  decremented_part <- head(R1_iter, num_to_decrement) - 1
                  tail_part <- tail(R1_iter, max(0, length(R1_iter) - num_to_decrement))
                  R1_iter <- c(decremented_part, tail_part)
                }
              }
            }
          })
