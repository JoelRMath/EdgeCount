library(data.table)
library(EdgeCount)

# Function to prepare the graph data and mappings
vsea_gj_step1_setup <- function(object, element_ranks) {

  # 1. Validation & Reduction
  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) stop("No valid elements found.")

  # Crucial: This ensures we only work with the relevant universe
  analysis_object <- reduce_universe_by_elements(object, valid_elements)

  # 2. Extract Edges
  dt_edges <- as.data.table(to_dataframe(analysis_object))
  setnames(dt_edges, c("term", "element"), c("term_id", "element_id"))

  # 3. Create Integer Maps (The "Dictionaries")
  u_terms <- unique(dt_edges$term_id)
  u_elems <- unique(dt_edges$element_id) # These are string IDs

  # Map String -> Integer (1..N)
  term_map <- setNames(seq_along(u_terms), u_terms)
  elem_map <- setNames(seq_along(u_elems), u_elems)

  # 4. Convert Edges to Integer Indices
  dt_edges[, `:=`(
    t_idx = term_map[term_id],
    e_idx = elem_map[element_id]
  )]

  # 5. Pre-calculate Degrees using Integer Indexing
  # We need these for the cumulative sum calculation later
  # deg_vec[i] corresponds to the degree of the element with e_idx == i
  deg_vec <- unlist(analysis_object@ecprob@degrees[u_elems])

  # Return everything we need for Step 2
  return(list(
    analysis_object = analysis_object,
    valid_elements = valid_elements,
    dt_edges = dt_edges, # The topology (Terms <-> Elements)
    deg_vec = deg_vec,   # Element degrees (for Lambda)
    elem_map = elem_map, # To map ranks correctly
    n_elem = length(u_elems)
  ))
}

vsea_gj_step2_perms <- function(step1_out, element_ranks, n_permutations = 2) {

  n_elem <- step1_out$n_elem

  # 1. Get Real Ranks aligned to e_idx (1..N)
  # elem_map names are the original IDs ordered 1..N
  ordered_ids <- names(step1_out$elem_map)

  # Extract ranks in the correct order (1..N)
  # rank() handles ties/NAs if any, though input should be clean
  rank_vec_real <- rank(element_ranks[ordered_ids])

  # 2. Build Real Data Table (perm_id = 0)
  dt_real <- data.table(
    perm_id = 0L,
    e_idx = 1:n_elem,
    rank = as.vector(rank_vec_real) # Ensure simple vector
  )

  # 3. Build Null Data Table (perm_id = 1..N)
  # Shuffle the values of rank_vec_real
  perm_matrix <- replicate(n_permutations, sample(rank_vec_real), simplify = "matrix")

  dt_nulls <- data.table(
    perm_id = rep(1:n_permutations, each = n_elem),
    e_idx = rep(1:n_elem, n_permutations),
    rank = as.vector(perm_matrix)
  )

  # Combine
  dt_perms <- rbindlist(list(dt_real, dt_nulls))

  # 4. Attach Degrees
  # deg_vec is already aligned to e_idx 1..N
  dt_perms[, e_deg := step1_out$deg_vec[e_idx]]

  # 5. Calculate Global Cumulative Sum (The Fix)
  # We must sort by Permutation -> Rank to integrate over the ranked list
  setorder(dt_perms, perm_id, rank)

  # Calculate running sum of degrees
  dt_perms[, univ_cumsum := cumsum(e_deg), by = perm_id]

  return(dt_perms)
}

vsea_gj_step3_calc <- function(step1_out, step2_out) {

  dt_edges <- step1_out$dt_edges
  dt_perms <- step2_out
  setkey(dt_edges, e_idx)
  setkey(dt_perms, e_idx)

  # 1. Giant Join
  dt_joined <- dt_perms[dt_edges, on = "e_idx", allow.cartesian = TRUE]
  setorder(dt_joined, perm_id, t_idx, rank)

  # 2. FIX: Robust Term Degree Calculation
  # Count edges per term index directly from the edge table
  term_deg_dt <- step1_out$dt_edges[, .(t_deg_calc = .N), by = t_idx]

  # Add t_deg to the joined table via join (safest method)
  dt_joined[term_deg_dt, on = "t_idx", t_deg := i.t_deg_calc]

  # 3. Calculate Scores
  # Denominator is 2M (Sum of Element Degrees + Sum of Term Degrees)
  # In bipartite, Sum(Elem) == Sum(Term), so 2 * Sum(Elem) is correct.
  global_M2 <- 2 * sum(step1_out$deg_vec)

  dt_joined[, `:=`(
    observed_ec = 1:.N
  ), by = .(perm_id, t_idx)]

  dt_joined[, lambda := (t_deg / global_M2) * univ_cumsum]
  dt_joined[, max_ec := pmin(t_deg, rank)]
  dt_joined[, score := 0.5 * (log2(observed_ec + 0.375) - log2(lambda + 0.375))]

  return(dt_joined)
}

# --- CORRECTED FULL DEBUG FUNCTION ---
vsea_gj_debug_full <- function(object, element_ranks, n_permutations = 10, seed = 123) {

  if (!is.null(seed)) set.seed(seed)

  valid_elements <- intersect(names(element_ranks), object@elements)
  analysis_object <- reduce_universe_by_elements(object, valid_elements)

  dt_edges <- as.data.table(to_dataframe(analysis_object))
  setnames(dt_edges, c("term", "element"), c("term_id", "element_id"))

  u_terms <- unique(dt_edges$term_id)
  u_elems <- unique(dt_edges$element_id)

  term_map <- setNames(seq_along(u_terms), u_terms)
  elem_map <- setNames(seq_along(u_elems), u_elems)

  dt_edges[, `:=`(t_idx = term_map[term_id], e_idx = elem_map[element_id])]

  deg_vec <- unlist(analysis_object@ecprob@degrees[u_elems])

  # FIX: Robust Term Degree Calculation
  term_deg_dt <- dt_edges[, .(t_deg_calc = .N), by = t_idx]

  rank_vec_real <- rank(element_ranks[u_elems])
  n_elem <- length(u_elems)

  perm_matrix <- replicate(n_permutations, sample(rank_vec_real), simplify = "matrix")

  dt_perms <- data.table(
    perm_id = rep(1:n_permutations, each = n_elem),
    e_idx = rep(1:n_elem, n_permutations),
    rank = as.vector(perm_matrix)
  )

  dt_real <- data.table(perm_id = 0L, e_idx = 1:n_elem, rank = as.integer(rank_vec_real))
  dt_perms <- rbindlist(list(dt_real, dt_perms))

  dt_perms[, e_deg := deg_vec[e_idx]]
  setorder(dt_perms, perm_id, rank)
  dt_perms[, univ_cumsum := cumsum(e_deg), by = perm_id]

  setkey(dt_edges, e_idx)
  setkey(dt_perms, e_idx)
  dt_joined <- dt_perms[dt_edges, on = "e_idx", allow.cartesian = TRUE]
  setorder(dt_joined, perm_id, t_idx, rank)

  # Add Term Degrees
  dt_joined[term_deg_dt, on = "t_idx", t_deg := i.t_deg_calc]

  global_M2 <- 2 * analysis_object@ecprob@graph_size

  dt_joined[, `:=`(
    observed_ec = 1:.N
  ), by = .(perm_id, t_idx)]

  dt_joined[, lambda := (t_deg / global_M2) * univ_cumsum]
  dt_joined[, max_ec := pmin(t_deg, rank)]
  dt_joined[, score := 0.5 * (log2(observed_ec + 0.375) - log2(lambda + 0.375))]

  dt_joined[, term_id := u_terms[t_idx]]

  return(dt_joined)
}

###########################

term_element_df <- data.frame(
  Term = c("T1", "T1", "T2", "T2", "T3"),
  Element = c("11", "12", "12", "13", "14"),
  stringsAsFactors = FALSE
)
toy_ects <- ECTermScoring(term_element_df)
set.seed(123)
# Create a dummy ranking
all_elements <- toy_ects@elements
element_ranks <- setNames(sample(seq_along(all_elements)), all_elements)

# Run Step 1
step1_out <- vsea_gj_step1_setup(toy_ects, element_ranks)

# Verification
print(paste("Elements in universe:", step1_out$n_elem))
print(paste("Length of Degree Vector:", length(step1_out$deg_vec)))
print("Edge Table:")
#print(step1_out$dt_edges)

# Step 2

set.seed(123)
step2_out <- vsea_gj_step2_perms(step1_out, element_ranks, n_permutations = 2)
#print(step2_out)

# Step 3

step3_out <- vsea_gj_step3_calc(step1_out, step2_out)
# print(step3_out)
# print("--- Comparison: Term 1 across Permutations ---")
# print(step3_out[t_idx == 1, .(perm_id, rank, score)])

# full version

full_out <- vsea_gj_debug_full(toy_ects, element_ranks, n_permutations = 2, seed = 123)

# compare full and step 3

setcolorder(step3_out, sort(names(step3_out)))
setcolorder(full_out, sort(names(full_out)))

# 2. Check Equality
# check.attributes = FALSE allows ignoring row.names or slight metadata differences
is_identical <- all.equal(step3_out, full_out, check.attributes = FALSE)

if (isTRUE(is_identical)) {
  print("SUCCESS: The Full Debug Function matches the Step-by-Step logic exactly.")
} else {
  print("FAILURE: Mismatch detected.")
  print(is_identical)
}

# 3. Visual Check of Scores
# We want to see if Perm 0 (Real) scores differ from Perm 1 (Null)
print("--- Scores for T1 across Permutations ---")
print(full_out[t_idx == 1, .(perm_id, rank, score)])


