network_dt <- data.table::fread("data-raw/raw_network.txt")
bipartite_dt <- data.table::fread("data-raw/raw_bipartite.txt")

network_universe <- unique(c(network_dt[[1]], network_dt[[2]]))
bipartite_universe <- unique(bipartite_dt$element)

common_universe <- intersect(network_universe, bipartite_universe)

network_subset_dt <- network_dt[element1 %in% common_universe & element2 %in% common_universe]
bipartite_subset_dt <- bipartite_dt[element %in% common_universe]


ecg <- ECGraph(network_subset_dt)

trim_ecgraph <- function(ecg, threshold){

  df <- data.frame(degrees = as.numeric(unlist(ecg@degrees)),
                   id = as.character(names(ecg@degrees)))
  df <- df[order(df$degrees, decreasing = TRUE),]
  two_m <- sum(df$degrees)
  max_k <- df[1,1]
  criterion <- max_k^2/two_m
  print(criterion)
  current_row <- 1
  removed_vertices <- NULL
  while(criterion > threshold){
    vertex_to_remove <- df[current_row, 2]
    removed_vertices <- c(removed_vertices, vertex_to_remove)
    hood_to_remove <- get_neighbors(ecg, vertex_to_remove)
    for (i in 1:length(hood_to_remove)){
      v <- hood_to_remove[i]
      ecg@adj$v <- setdiff(ecg@adj$v, vertex_to_remove)
      df[df$id == v,1] <- df[df$id == v,1] - 1
    }
    df[1,1] <- 0
    df <- df[order(df$degrees, decreasing = TRUE),]
    two_m <- sum(df$degrees)
    max_k <- df[1,1]
    criterion <- max_k^2/two_m
    print(paste(max_k,criterion))
  }
  remaining_vertices <- setdiff(ecg@names, removed_vertices)
  ecg@adj <-  ecg@adj[remaining_vertices]
  # next, create the new ecg
}

trim_ecgraph(ecg, 0.5)
