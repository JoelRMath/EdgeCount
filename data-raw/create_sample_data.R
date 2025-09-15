network_dt <- data.table::fread("data-raw/raw_network.txt")
bipartite_dt <- data.table::fread("data-raw/raw_bipartite.txt")

network_universe <- unique(c(network_dt[[1]], network_dt[[2]]))
bipartite_universe <- unique(bipartite_dt$element)

common_universe <- intersect(network_universe, bipartite_universe)

network_subset_dt <- network_dt[element1 %in% common_universe & element2 %in% common_universe]
bipartite_subset_dt <- bipartite_dt[element %in% common_universe]


ecg <- ECGraph(network_subset_dt)
lst <- trim_ecgraph(ecg, 0.5)
ecg <- remove_isolated_vertices(lst[["trimmed_graph"]])

ects <- ECTermScoring(bipartite_subset_dt)
lst <- trim_bipartite_terms(ects, 0.2)
ects <- remove_isolated_elements(lst[["trimmed_object"]])
ects <- remove_empty_terms(ects)
