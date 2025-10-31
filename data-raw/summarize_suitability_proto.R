library(EdgeCount)
library(data.table)

summarize_suitability_fast_proto <- function(object) {

  degrees <- unlist(object@degrees)
  N <- length(degrees)
  if (N == 0) {
    return(list(
      pij_over_1 = 0,
      prop_problematic_vertices = 0,
      summary_pij = summary(numeric(0))
    ))
  }

  M <- object@graph_size
  two_M <- 2 * M

  degree_distribution <- table(degrees)
  k <- as.numeric(names(degree_distribution))
  n_i <- as.numeric(degree_distribution) # n_i is the count for k_i
  p <- n_i / N

  # --- Create probability and product matrices ---
  prod_matrix <- outer(k, k, "*") / two_M
  p_matrix <- outer(p, p)

  # --- Correct the diagonal ---

  # 1. Calculate the probability of all self-pairs for each degree class
  # P(self-pair of degree k_i) = n_i * (1/N)^2 = (n_i/N) * (1/N) = p_i / N
  diag_p_self_pair <- p / N

  # 2. Subtract this self-pair probability from the diagonal of p_matrix
  # This correctly leaves only the probability of non-self-pairs
  # where both nodes happen to have the same degree.
  diag(p_matrix) <- diag(p_matrix) - diag_p_self_pair
  print(diag_p_self_pair)

  # 3. Set the product matrix diagonal to 0 (self-pairs are not counted)
  diag(prod_matrix) <- 0

  # --- Calculate the final proportions ---
  prod_vec <- as.vector(prod_matrix)
  p_vec <- as.vector(p_matrix)

  # This is the total probability of picking any non-self-pair
  total_p_nonself <- sum(p_vec)

  # This is the probability of picking a "bad" non-self-pair
  prop_bad_nonself <- sum(p_vec[prod_vec >= 1])

  # The final metric is the conditional probability
  pij_over_1 <- if(total_p_nonself == 0) 0 else (prop_bad_nonself / total_p_nonself)

  # --- Proportion of problematic vertices ---
  threshold <- sqrt(two_M)
  problematic_vertices <- degrees[degrees > threshold]
  prop_problematic_vertices <- length(problematic_vertices) / N

  return(list(
    pij_over_1 = pij_over_1,
    prop_problematic_vertices = prop_problematic_vertices,
    summary_pij = summary(prod_vec[prod_vec > 0]) # Summary of non-self pairs
  ))
}

summarize_suitability_fast_proto_2 <- function(object) {

  degrees <- unlist(object@degrees)
  N <- length(degrees)
  if (N == 0) {
    return(list(
      pij_over_1 = 0,
      prop_problematic_vertices = 0,
      summary_pij = summary(numeric(0))
    ))
  }

  M <- object@graph_size
  two_M <- 2 * M

  degree_distribution <- table(degrees)
  k <- as.numeric(names(degree_distribution))

  # -- proportions of unique distinct degree pairs (undirected) ---
  l_k <- length(k)
  n_k <- as.numeric(degree_distribution)
  count_matrix <- outer(n_k, n_k)
  diag(count_matrix) <- n_k * (n_k - 1) / 2
  undirected_matrix <- matrix(0.5, l_k, l_k) + diag(0.5, l_k)
  undirected_vec <- as.vector(undirected_matrix)
  count_vec <- as.vector(count_matrix)
  count_vec <- count_vec*undirected_vec
  prop_vec <- count_vec / sum(count_vec)

  # --- pij distribution for undirected edges ---
  pij_matrix <- outer(k, k) / two_M
  pij_vec <- as.vector(pij_matrix)
  dt <- data.table(pij = pij_vec,
                   count = count_vec,
                   prop = prop_vec)

  dt <- dt[, .(
    count = as.integer(sum(count)),
    prop = sum(prop)
  ), by = .(pij)]

  pij_over_1 = sum(dt[pij >= 1, prop])

  # --- proportion of problematic vertices ---
  threshold <- sqrt(two_M)
  problematic_vertices <- degrees[degrees > threshold]
  prop_problematic_vertices <- length(problematic_vertices) / N

  # -- results including summary of pij distribution
  return(list(
    pij_over_1 = pij_over_1,
    prop_problematic_vertices = prop_problematic_vertices,
    summary_pij = summary(rep(dt$pij, dt$count)),
    pij_distribution = dt

  ))
}

summarize_suitability_fast_proto_3 <- function(object) {

  degrees <- unlist(object@degrees)
  N <- length(degrees)
  if (N <= 1) {
    return(list(
      pij_over_1 = 0,
      prop_problematic_vertices = 0,
      summary_pij_nonself = summary(numeric(0))
    ))
  }

  M <- object@graph_size
  two_M <- 2 * M

  degree_distribution <- table(degrees)
  k <- as.numeric(names(degree_distribution))
  n_k <- as.numeric(degree_distribution)

  # 1. Create the k x k matrix of p_ij values
  prod_matrix <- outer(k, k) / two_M

  # 2. Create the k x k matrix of non-self-pair *COUNTS*

  # Start with all pairs (i, j)
  count_matrix <- outer(n_k, n_k)

  # Correct the diagonal: replace n_i*n_i with n_i*(n_i-1)/2
  diag(count_matrix) <- n_k * (n_k - 1) / 2

  # 3. Get the distribution of p_ij for non-self-pairs

  # Extract the p_ij values from the upper triangle (to avoid double counting)
  pij_vec <- prod_matrix[upper.tri(prod_matrix, diag = TRUE)]

  # Extract the corresponding *counts* from the upper triangle
  count_vec <- count_matrix[upper.tri(count_matrix, diag = TRUE)]

  # The total number of non-self-pairs is just sum(count_vec)
  # This should equal N * (N - 1) / 2
  total_pairs <- sum(count_vec)

  # Calculate the proportion for each p_ij bin
  prop_vec <- count_vec / total_pairs

  # 4. Create the final distribution table
  df <- data.table(pij = pij_vec, count = count_vec, prop = prop_vec)

  # 5. The final conditional probability
  pij_over_1 <- sum(df[pij >= 1, prop])

  # --- Proportion of problematic vertices ---
  threshold <- sqrt(two_M)
  problematic_vertices <- degrees[degrees > threshold]
  prop_problematic_vertices <- length(problematic_vertices) / N

  return(list(
    pij_over_1 = pij_over_1,
    prop_problematic_vertices = prop_problematic_vertices,
    summary_pij_nonself = summary(rep(df$pij, df$count)), # Weighted summary
    pij_distribution = df[count > 0] # Show the full distribution
  ))
}
# ---
# SCRIPT TO RUN THE FUNCTION
# ---
edge_df <- data.frame(
  p1 = c("A", "A", "A", "A", "A", "B", "B", "B", "B", "C", "C", "C"),
  p2 = c("B", "C", "D", "E", "F", "C", "G", "H", "I", "J", "K", "L"),
  stringsAsFactors = FALSE
)

ecp <- ECProb(ECGraph(edge_df))

suitability_stats <- summarize_suitability_fast_proto_2(ecp)
print(suitability_stats)


