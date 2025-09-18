library(data.table)
library(org.Hs.eg.db)
library(GO.db)


# Creates a consistent pair (same element universe)
# of sample datasets for the EdgeCount package.
# Returns a list: $sample_ecg, $sample_ects

create_example_data <- function(network_threshold = 0.01,
                                bipartite_threshold = 0.01,
                                min_term_size = 2) {

  # Trims the Network Graph
  network_dt_full <- data.table::fread("data-raw/raw_network.txt")
  network_dt_subset <- network_dt_full[best_combined_score >= 900, .(gene1, gene2)]
  ecg_candidate <- ECGraph(network_dt_subset)
  ecg_candidate <- trim_ecgraph(ecg_candidate, network_threshold)[["trimmed_graph"]]
  ecg_candidate <- remove_isolated_vertices(ecg_candidate)

  # Gets the initial candidate edge list for the network
  net_edges <- data.table(to_dataframe(ecg_candidate))

  # Trims the Bipartite Graph
  bipartite_dt_full <- data.table::fread("data-raw/raw_bipartite.txt")
  ects_candidate <- ECTermScoring(bipartite_dt_full)
  ects_candidate <- trim_bipartite_terms(ects_candidate, bipartite_threshold)[["trimmed_object"]]
  ects_candidate <- remove_empty_terms(ects_candidate)
  ects_candidate <- remove_isolated_elements(ects_candidate)

  # Gets the initial candidate edge list for the bipartite graph
  bip_edges <- data.table(to_dataframe(ects_candidate))

  term_sizes <- bip_edges[, .N, by = term]
  kept_terms <- term_sizes[N >= min_term_size, term]
  bip_edges <- bip_edges[term %in% kept_terms]

  last_universe_size <- 0

  while(TRUE) {
    net_universe <- unique(c(net_edges$from, net_edges$to))
    bip_universe <- unique(bip_edges$element)

    current_universe_size <- length(intersect(net_universe, bip_universe))
    message(paste("Current common universe size:", current_universe_size))

    # If stable.
    if (current_universe_size == last_universe_size) {
      break
    }

    last_universe_size <- current_universe_size
    common_universe <- intersect(net_universe, bip_universe)

    net_edges <- net_edges[from %in% common_universe & to %in% common_universe]
    bip_edges <- bip_edges[element %in% common_universe]

    term_counts <- bip_edges[, .N, by = term]
    bip_edges <- bip_edges[term %in% term_counts[N > 0, term]]
  }

  final_common_universe <- common_universe

  final_ecg <- ECGraph(net_edges)
  final_ects <- ECTermScoring(bip_edges)

  message("\n--- Verification ---")
  final_check_universe <- intersect(final_ecg@names, final_ects@elements)
  message(paste("Final common universe size is", length(final_check_universe)))
  message(paste("Final network has", length(final_ecg@names), "elements."))
  message(paste("... and", sum(unlist(ecg@degrees)/2)), " edges")
  message(paste("Final bipartite graph has", length(final_ects@elements), "elements."))
  message(paste("... and", final_ects@ecprob@graph_size, " edges"))

  return(list(
    sample_ecg = final_ecg,
    sample_ects = final_ects
  ))
}

# Creates annotation data for sample datasets

create_sample_annotations <- function(sample_ects) {

  # GO term names
  unique_go_terms <- sample_ects@terms
  go_term_names <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = unique_go_terms,
    keytype = "GOID",
    columns = c("TERM", "ONTOLOGY")
  )
  names(go_term_names) <- c("term", "term_name", "ontology")


  # Gene symbols
  unique_gene_ids <- sample_ects@elements
  gene_symbols <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique_gene_ids,
    keytype = "ENTREZID",
    columns = "SYMBOL"
  )
  names(gene_symbols) <- c("element", "symbol")

  return(list(
    sample_term_names = go_term_names,
    sample_gene_symbols = gene_symbols
  ))
}

# Creating sample data and annotation

sample_data <- create_example_data()
sample_ecg <- sample_data$sample_ecg
sample_ects <- sample_data$sample_ects
annotation_data <- create_sample_annotations(sample_ects)
sample_term_names <- annotation_data$sample_term_names
sample_gene_symbols <- annotation_data$sample_gene_symbols

usethis::use_data(
  sample_ecg,
  sample_ects,
  sample_term_names,
  sample_gene_symbols,
  overwrite = TRUE
)
