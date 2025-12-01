library(dqrng)
library(matrixStats)
library(Matrix)
library(EdgeCount)
library(parallel)

# --- 1. Data Loading & Constants ---
data("sample_ects")

elements <- sample_ects@elements
terms    <- sample_ects@terms

Ne <- length(elements)
Nt <- length(terms)
degrees <- unlist(sample_ects@ecprob@degrees[elements], use.names = FALSE)

n_sim <- 1000
total_len <- as.numeric(Ne) * n_sim

message(paste("Setup: Ne =", Ne, "| Nt =", Nt, "| Simulations =", n_sim))

# --- 2. Generate Random Permutations ---
message("--- Generating Random Permutations ---")
sim_vector <- integer(total_len)
for (i in 1:n_sim) {
  start_idx <- (i - 1) * Ne + 1
  end_idx   <- start_idx + Ne - 1
  sim_vector[start_idx:end_idx] <- dqrng::dqsample.int(Ne)
}

# --- 3. Build Cumulative Degree Matrix ---
message("--- Building Cumulative Degree Matrix ---")
sim_degrees <- degrees[sim_vector]
dim(sim_degrees) <- c(Ne, n_sim)
sim_cum_degrees <- matrixStats::colCumsums(sim_degrees)
rm(sim_degrees); gc() # Clean up

# --- 4. Prepare Integer Adjacency (Element -> Terms) ---
# We need a fast lookup: "Which terms does Element X connect to?"
message("--- Building Integer Adjacency Structure ---")
term_map <- setNames(seq_len(Nt), terms)
aligned_adj <- sample_ects@ecprob@adj[elements]

# Convert names to integers
# This list allows us to expand Elements into (Term, Rank) pairs instantly
adj_int <- lapply(aligned_adj, function(x) {
  if (is.null(x) || length(x) == 0) return(integer(0))
  unname(term_map[x])
})

M <- length(unlist(adj_int)) # Total Edges

# --- 5. Simulation-Wise Scoring Function ---
# Instead of looping Terms, we process ONE simulation completely.
# This exploits the fact that 'sim_vector' is already sorted by Rank (1..Ne).

score_simulation <- function(sim_idx) {

  # A. Retrieve Data for this Simulation
  start_idx <- (sim_idx - 1) * Ne + 1
  end_idx   <- start_idx + Ne - 1

  # The permutation IS the elements ordered by rank (1, 2, ..., Ne)
  perm_elements <- sim_vector[start_idx:end_idx]

  # B. Expand to Edges (The "Natural Sort" Trick)
  # We iterate ranks 1..Ne. We look up terms for each element.
  # The resulting Ranks are GUARANTEED to be sorted ascendingly for every term.

  # Get terms for the elements in rank order
  # This creates a list of terms corresponding to Rank 1, Rank 2...
  terms_in_rank_order <- adj_int[perm_elements]

  # Flatten to vectors
  # edge_terms: The term ID for every edge found
  # edge_ranks: The rank (1..Ne) where that edge was found
  edge_terms <- unlist(terms_in_rank_order, use.names = FALSE)
  edge_ranks <- rep(seq_len(Ne), lengths(terms_in_rank_order))

  # C. Group by Term (Stable Sort)
  # We need to group data by Term to calculate scores.
  # CRITICAL: We use "radix" sort which is STABLE.
  # This preserves the relative order of 'edge_ranks', keeping them sorted!
  ord <- order(edge_terms, method = "radix")

  sorted_terms <- edge_terms[ord]
  sorted_ranks <- edge_ranks[ord]

  # D. Vectorized Score Calculation
  # We need 'x' (1, 2, 3...) for each term group.
  # Since sorted_terms is grouped, we can generate sequences efficiently.

  # Find run lengths (how many edges each term has in this sim)
  # rle is very fast on sorted data
  run_lens <- rle(sorted_terms)$lengths

  # Generate x = 1, 2, 3, 1, 2...
  # sequence() is a fast C primitive
  x_vals <- sequence(run_lens)

  # Lookup Lambda
  # We use the specific column for this simulation
  # sim_cum_degrees[rank, sim_idx]
  lambdas <- (sim_cum_degrees[sorted_ranks, sim_idx] * length(edge_ranks)) / (2 * M)
  # Note: The formula uses k_t (degree of term).
  # But wait, k_t is constant per term. We can get it from run_lens?
  # NO. run_lens is the observed count in this sim (which is <= k_t).
  # We need the ACTUAL global k_t for the term.

  # Fix for k_t:
  # We need a global lookup for k_t.
  # Pre-calculated global k_t map (1..Nt)
  # global_kt <- lengths(term_to_elements) # We need to pass this in or calc it globally
  # For now, let's use a quick lookup table
  # (In production, pass 'global_kt' to this function to avoid recalculation)

  # Let's assume global_kt vector exists in environment (calculated below)
  term_kt <- global_kt[sorted_terms]
  lambdas <- (sim_cum_degrees[sorted_ranks, sim_idx] * term_kt) / (2 * M)

  # Calculate Scores
  scores <- 0.5 * (log2(x_vals + 0.375) - log2(lambdas + 0.375))

  # E. Aggregate Max Score per Term
  # We have 'scores' and 'sorted_terms'. We want max per term.
  # Since it's sorted, tapply is reasonably fast, but 'data.table' would be instant.
  # Using base R 'tapply' is the safest dependency-free approach.

  # This returns a named vector (Name=TermID, Value=MaxScore)
  max_scores <- tapply(scores, sorted_terms, max)

  return(max_scores)
}

# --- 6. Pre-calculate Global Constants ---
# Reconstruct the Term->Element map just to get accurate k_t counts
# (Or we can derive it from adj_int more slowly, this is fast enough)
u_elements <- rep(seq_len(Ne), lengths(aligned_adj))
u_terms    <- unlist(adj_int, use.names = FALSE)
global_kt  <- tabulate(u_terms, nbins = Nt)

# --- 7. Execution (Simulation-Wise) ---
message("--- Computing Null Scores (Sim-Wise Optimization) ---")

n_cores <- parallel::detectCores() - 1
if(n_cores < 1) n_cores <- 1

print(system.time({
  # Iterate over SIMULATIONS, not Terms.
  # Result is a List of Named Vectors
  results_list <- parallel::mclapply(1:n_sim, score_simulation, mc.cores = n_cores)
}))

message("--- Pooling and Sorting Global Null ---")
print(system.time({
  # Flatten the list of results
  # Note: This pool contains only terms that had at least 1 edge in a simulation.
  # Terms with 0 edges in a sim have score -Inf (implicitly), so excluding them is fine/correct
  # for the "Max Effect Size" distribution (max(-Inf) is irrelevant).

  null_scores_pool <- unlist(results_list, use.names = FALSE)
  null_dist_sorted <- sort(null_scores_pool)
}))

message(paste("Null Distribution Built. Size:", length(null_dist_sorted)))
