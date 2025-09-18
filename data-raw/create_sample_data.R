library(data.table)
library(org.Hs.eg.db)

parse_string <- function() {

  message("Reading alias file...")
  alias_file <- "data-raw/9606.protein.aliases.v12.0.txt.gz"
  header_line <- readLines(alias_file, n = 1)
  clean_header <- sub("^#", "", header_line)
  header_names <- strsplit(clean_header, "\\s+")[[1]]
  alias_dt <- data.table::fread(alias_file, skip = 1, header = FALSE)
  data.table::setnames(alias_dt, old = names(alias_dt), new = header_names)

  message("Preparing mapping table...")
  entrez_source_name <- "Ensembl_HGNC_entrez_id" # Confirmed source name
  mapping_dt <- alias_dt[source == entrez_source_name, .(
    protein_id = string_protein_id,
    gene_id = alias
  )]
  mapping_dt[, protein_id := sub("9606.", "", protein_id, fixed = TRUE)]
  mapping_dt <- unique(mapping_dt)
  data.table::setkey(mapping_dt, protein_id)

  message("Reading and processing protein network file...")
  network_dt <- data.table::fread("data-raw/9606.protein.links.v12.0.txt")

  network_dt[, protein1 := sub("9606.", "", protein1, fixed = TRUE)]
  network_dt[, protein2 := sub("9606.", "", protein2, fixed = TRUE)]

  network_dt[, c("p_canon_1", "p_canon_2") := .(pmin(protein1, protein2), pmax(protein1, protein2))]

  protein_network_dt <- network_dt[, .(best_score = max(combined_score)), by = .(p_canon_1, p_canon_2)]


  message("Mapping protein IDs to gene IDs...")
  protein_network_dt[mapping_dt, on = .(p_canon_1 = protein_id), gene1 := i.gene_id]
  protein_network_dt[mapping_dt, on = .(p_canon_2 = protein_id), gene2 := i.gene_id]


  message("Cleaning final gene network...")

  gene_network_with_scores <- na.omit(protein_network_dt, cols = c("gene1", "gene2"))

  gene_network_with_scores <- gene_network_with_scores[gene1 != gene2]

  gene_network_with_scores[, c("g_canon_1", "g_canon_2") := .(pmin(gene1, gene2), pmax(gene1, gene2))]

  final_gene_network <- gene_network_with_scores[,
                                                 .(best_combined_score = max(best_score)),
                                                 by = .(g_canon_1, g_canon_2)
  ]

  data.table::setnames(final_gene_network, c("g_canon_1", "g_canon_2"), c("gene1", "gene2"))
  message("Done.")
  return(final_gene_network)
}


calculate_network_properties_by_threshold <- function(dt) {

  if (!is.data.table(dt) || nrow(dt) == 0) {
    return(data.table(threshold = numeric(), order = integer(), size = integer()))
  }

 summary_dt <- dt[, .(
    size_of_block = .N,
    nodes_in_block = list(unique(c(gene1, gene2)))
  ), by = best_combined_score]

  data.table::setorder(summary_dt, -best_combined_score)

  summary_dt[, cumulative_size := cumsum(size_of_block)]

  cumulative_node_sets <- Reduce(union, summary_dt$nodes_in_block, accumulate = TRUE)

  summary_dt[, cumulative_order := lengths(cumulative_node_sets)]

  final_summary <- summary_dt[, .(
    threshold = best_combined_score,
    order = cumulative_order,
    size = cumulative_size
  )]

  data.table::setorder(final_summary, threshold)

  return(final_summary)
}

dt <- parse_string()
data.table::setorder(dt, best_combined_score)
data.table::fwrite(dt, "data-raw/raw_network.txt", sep = "\t")
network_summary <- calculate_network_properties_by_threshold(dt)
data.table::fwrite(network_summary, "data-raw/network_summary.txt", sep = "\t")

network <- dt[best_combined_score >= 900, .(gene1, gene2)]

xx <- as.list(org.Hs.egGO2EG)
xx <- lapply(xx, unique)
bipartite_dt <- data.table::as.data.table(utils::stack(xx))
data.table::setnames(bipartite_dt,
                     old = c("values", "ind"),
                     new = c("element", "term"))
bipartite_dt <- bipartite_dt[, .(term, element)]
bipartite_dt[, term := as.character(term)]
data.table::fwrite(bipartite_dt, "data-raw/raw_bipartite.txt", sep = "\t")

