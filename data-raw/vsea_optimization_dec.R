library(dqrng)
library(matrixStats)
library(Matrix)
library(EdgeCount)
library(parallel)
library(data.table)

# Prevent data.table from spawning threads inside mclapply workers
setDTthreads(1)

#' Run VSEA (Vectorized Set Enrichment Analysis)
#'
#' @param ects An object of class 'ECTermScoring' containing the bipartite graph.
#' @param element_ranking A character vector of element names (e.g., genes) in the order
#'        to be tested (e.g., sorted by differential expression).
#' @param n_sim Number of simulations for the null distribution.
#' @param n_cores Number of cores for parallel processing. Defaults to detectCores() - 1.
#'
#' @return A data.table containing the enriched terms, statistics, and FDR.
run_vsea <- function(ects, element_ranking, n_sim = 1000, n_cores = NULL) {

  # --- 1. Setup & Data Validation ---
  elements <- ects@elements
  terms    <- ects@terms

  Ne <- length(elements)
  Nt <- length(terms)

  message(paste("Setup: Ne =", Ne, "| Nt =", Nt, "| Simulations =", n_sim))

  # Ensure strict character conversion for reliable matching
  element_ranking <- as.character(element_ranking)
  elements_char   <- as.character(elements)

  # Check coverage
  if (!all(element_ranking %in% elements_char)) {
    stop("Some elements in 'element_ranking' are not found in the ECTS object.")
  }
  if (length(element_ranking) != Ne) {
    warning("Length of 'element_ranking' differs from total elements in graph. Using intersection.")
  }

  # Map Internal Graph Indices (1..Ne) based on the order in the OBJECT
  degrees <- unlist(ects@ecprob@degrees[elements], use.names = FALSE)

  # Create a Map: Name -> Internal Index
  # Force elements to be names (character) to allow name-based lookup later
  element_name_to_idx <- setNames(seq_len(Ne), elements_char)

  # --- 2. Prepare Graph Structures (Internal Indices) ---
  message("--- Preparing Graph Structures ---")

  # Term Map: Name -> Index
  term_map <- setNames(seq_len(Nt), terms)

  # Adjacency: Element Index -> Vector of Term Indices
  aligned_adj <- ects@ecprob@adj[elements]

  adj_int <- lapply(aligned_adj, function(x) {
    if (is.null(x) || length(x) == 0) return(integer(0))
    # CRITICAL: Force character lookup to avoid positional indexing if 'x' is numeric
    unname(term_map[as.character(x)])
  })

  M <- length(unlist(adj_int)) # Total Edges in graph

  # Global Term Degrees (k_t)
  u_terms    <- unlist(adj_int, use.names = FALSE)
  global_kt  <- tabulate(u_terms, nbins = Nt)

  # --- 3. Build Null Model (Simulations) ---
  message("--- Generating Null Distribution ---")

  total_len <- as.numeric(Ne) * n_sim
  sim_vector <- integer(total_len)

  # A. Generate Random Permutations
  for (i in 1:n_sim) {
    start_idx <- (i - 1) * Ne + 1
    end_idx   <- start_idx + Ne - 1
    sim_vector[start_idx:end_idx] <- dqrng::dqsample.int(Ne)
  }

  # B. Build Cumulative Degree Matrix (Null)
  sim_degrees <- degrees[sim_vector]
  dim(sim_degrees) <- c(Ne, n_sim)
  sim_cum_degrees <- matrixStats::colCumsums(sim_degrees)
  rm(sim_degrees); gc()

  # C. Simulation-Wise Scoring Function
  score_simulation <- function(sim_idx) {
    start_idx <- (sim_idx - 1) * Ne + 1
    end_idx   <- start_idx + Ne - 1

    perm_elements <- sim_vector[start_idx:end_idx]

    # Expand to Edges (Natural Sort)
    terms_in_rank_order <- adj_int[perm_elements]

    edge_terms <- unlist(terms_in_rank_order, use.names = FALSE)
    edge_ranks <- rep(seq_len(Ne), lengths(terms_in_rank_order))

    # Group by Term (Stable Sort)
    ord <- order(edge_terms, method = "radix")
    sorted_terms <- edge_terms[ord]
    sorted_ranks <- edge_ranks[ord]

    # Vectorized Score
    run_lens <- rle(sorted_terms)$lengths
    x_vals <- sequence(run_lens)
    term_kt_subset <- global_kt[sorted_terms]

    lambdas <- (sim_cum_degrees[sorted_ranks, sim_idx] * term_kt_subset) / (2 * M)
    scores <- 0.5 * (log2(x_vals + 0.375) - log2(lambdas + 0.375))

    # Aggregate Max
    dt <- data.table::data.table(tm = sorted_terms, sc = scores)
    max_scores <- dt[, max(sc), by = tm]$V1
    return(max_scores)
  }

  # D. Run Simulations (Parallel)
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1
    if(n_cores < 1) n_cores <- 1
  }

  results_list <- parallel::mclapply(1:n_sim, score_simulation, mc.cores = n_cores)

  # E. Pool Null Distribution
  null_scores_pool <- unlist(results_list, use.names = FALSE)
  null_dist_sorted <- sort(null_scores_pool)

  message(paste("Null Distribution Built. Size:", length(null_dist_sorted)))

  # --- 4. Calculate Observed Scores (Using User Ranking) ---
  message("--- Calculating Observed Scores & FDR ---")

  # A. Map User Ranking to Internal Indices
  # CRITICAL: Force character lookup.
  obs_indices <- element_name_to_idx[as.character(element_ranking)]

  # B. Prepare Data in Observed Order
  obs_degrees <- degrees[obs_indices]
  real_cum_degrees <- cumsum(obs_degrees)

  terms_in_rank_order <- adj_int[obs_indices]

  edge_terms <- unlist(terms_in_rank_order, use.names = FALSE)
  edge_ranks <- rep(seq_len(Ne), lengths(terms_in_rank_order))

  # C. Scoring Logic (Real Data)
  ord <- order(edge_terms, method = "radix")
  sorted_terms <- edge_terms[ord]
  sorted_ranks <- edge_ranks[ord]

  run_lens <- rle(sorted_terms)$lengths
  x_vals <- sequence(run_lens)
  term_kt_subset <- global_kt[sorted_terms]

  lambdas <- (real_cum_degrees[sorted_ranks] * term_kt_subset) / (2 * M)
  scores <- 0.5 * (log2(x_vals + 0.375) - log2(lambdas + 0.375))

  # D. Aggregate & FDR
  # Added ec = x_vals to track the observed edge count (1, 2, ... k) corresponding to each score
  dt_obs <- data.table::data.table(tm = sorted_terms, sc = scores, ec = x_vals)

  # Extract max score AND the observed edge count where that max occurred
  obs_agg <- dt_obs[, .(max_score = max(sc), observed_ec = ec[which.max(sc)]), by = tm]

  setorder(obs_agg, -max_score)

  n_null_total <- length(null_dist_sorted)
  n_null_ge <- n_null_total - findInterval(obs_agg$max_score - 1e-10, null_dist_sorted)

  obs_agg[, E_FalsePos := n_null_ge / n_sim]
  obs_agg[, Rank := .I]
  obs_agg[, FDR_raw := E_FalsePos / Rank]
  obs_agg[FDR_raw > 1, FDR_raw := 1]
  obs_agg[, FDR := rev(cummin(rev(FDR_raw)))]

  obs_agg[, TermName := terms[tm]]

  return(obs_agg)
}

# --- Usage Example: Positive Control ---

positive_control <- function(){

  data("sample_ects")

  # 1. Identify a valid target term to "spike"
  all_edge_terms <- unlist(sample_ects@ecprob@adj, use.names = FALSE)
  term_counts <- table(all_edge_terms)

  # Filter strictly for terms that EXIST in the object's term slot
  valid_terms <- intersect(names(term_counts), sample_ects@terms)
  term_counts <- term_counts[valid_terms]

  # Pick the first term with ~20 elements
  candidate_terms <- names(term_counts)[term_counts >= 5 & term_counts <= 15]

  if(length(candidate_terms) == 0) stop("No valid candidate terms found for control.")
  target_term <- candidate_terms[1]

  message(paste("POSITIVE CONTROL: Spiking term", target_term, "with", term_counts[target_term], "elements."))

  # 2. Find the elements belonging to this term
  all_elements <- sample_ects@elements
  adj <- sample_ects@ecprob@adj[all_elements]

  # Slow lookup for demo purposes only (Function uses optimized lookup)
  is_associated <- sapply(adj, function(x) target_term %in% x)
  term_elements <- all_elements[is_associated]

  # 3. Construct the Spiked Ranking
  # Place the term's elements at the VERY TOP (Ranks 1..k)
  background_elements <- setdiff(all_elements, term_elements)
  spiked_ranking <- c(term_elements, sample(background_elements))

  # 4. Run Analysis
  # We expect 'target_term' to be at Rank 1 with FDR 0
  res <- run_vsea(sample_ects, element_ranking = spiked_ranking, n_sim = 1000)

  # 5. Check Results
  print("--- Top 5 Enriched Terms ---")
  print(head(res, 10))

  print("--- Target Term Stats ---")
  print(res[TermName == target_term])
}

negative_control <- function(){

  data("sample_ects")
  element_ranks <- sample(sample_ects@elements)
  res <- run_vsea(sample_ects, element_ranking = element_ranks, n_sim = 1000)
  print(head(res, 10))
}

print("********* NEGATIVE CONTROL *****************")
print(system.time(
  negative_control()
))
print("************ POSITIVE CONTROL *********************")
print(system.time(
  positive_control()
))

