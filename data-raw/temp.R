library(data.table)
library(EdgeCount)

# 1. Setup Data
data("sample_ects")
set.seed(123)
# Create a dummy ranking
all_elements <- sample_ects@elements
scores <- setNames(runif(length(all_elements)), all_elements)
element_ranks <- rank(-scores)

# 2. Internal Prep (Mimicking run_vsea_analysis)
valid_elements <- intersect(names(element_ranks), sample_ects@elements)
analysis_object <- reduce_universe_by_elements(sample_ects, valid_elements)

# Map IDs
dt_edges <- as.data.table(to_dataframe(analysis_object))
setnames(dt_edges, c("term", "element"), c("term_id", "element_id"))
u_elems <- unique(dt_edges$element_id)
elem_map <- setNames(seq_along(u_elems), u_elems)
dt_edges[, e_idx := elem_map[element_id]]

# 3. Generate Ranks (Real vs Null)
rank_vec_real <- rank(element_ranks[u_elems])
n_elem <- length(u_elems)
n_permutations <- 5 # Small test

# Generate Permutations
perm_matrix <- replicate(n_permutations, sample(rank_vec_real), simplify = "matrix")

# Build Table
dt_perms <- data.table(
  perm_id = rep(1:n_permutations, each = n_elem),
  e_idx = rep(1:n_elem, n_permutations),
  rank = as.vector(perm_matrix)
)

dt_real <- data.table(perm_id = 0L, e_idx = 1:n_elem, rank = as.integer(rank_vec_real))
dt_all <- rbindlist(list(dt_real, dt_perms))

# 4. THE CHECK
# Pick one element (e.g., e_idx = 1) and check its rank across permutations
# If they are all identical, the shuffle failed.
check_ranks <- dt_all[e_idx == 1, .(perm_id, rank)]
print("--- Ranks for Element 1 across Permutations (Should vary) ---")
print(check_ranks)

# 5. Check Join Logic
setkey(dt_edges, e_idx)
setkey(dt_all, e_idx)
dt_joined <- dt_all[dt_edges, on = "e_idx", allow.cartesian = TRUE]

# Check if Perm IDs expanded correctly for a single Term
term_1_idx <- dt_edges$term_id[1] # Pick first term
subset_join <- dt_joined[e_idx == 1] # Look at element 1 rows
print("--- Joined Rows for Element 1 (Should have Perms 0..5) ---")
print(head(subset_join, 10))
