library(EdgeCount)
library(data.table)

# --- A small, self-contained toy example for step-by-step validation ---

# 1. The main graph
edges_dt <- data.table(p1 = c("E1", "E3", "E6"), p2 = c("E2", "E4", "E7"))
ecg <- ECGraph(edges_dt)
ecp <- ECProb(ecg)

# 2. The list of sets we want to test
sets_dt <- data.table(set_id = c("SA", "SB", "SC"))

# 3. The membership of elements in those sets
set_membership_dt <- data.table(
  set_id = c("SA", "SA", "SA", "SB", "SB", "SB", "SC"),
  element = c("E1", "E2", "E3", "E3", "E4", "E5", "E6") # E5 is not in the graph
)

message("--- STEP 1: Calculate Observed Edge Counts ---")

# 1a. Get the canonical edge list from the main graph
ecg_edges <- data.table(to_dataframe(ecp))
setnames(ecg_edges, c("from", "to"), c("e1", "e2"))
ecg_edges[, `:=`(canon1 = pmin(e1, e2), canon2 = pmax(e1, e2))]
ecg_edges <- unique(ecg_edges, by = c("canon1", "canon2"))
setkey(ecg_edges, canon1, canon2)

message("\n-- Canonical Edges in Graph --")
print(ecg_edges)

# 1b. Self-join the set memberships to get all possible intra-set element pairs
set_membership_dt_unique <- unique(set_membership_dt, by = c("set_id", "element"))
setkey(set_membership_dt_unique, set_id)
possible_edges <- set_membership_dt_unique[set_membership_dt_unique, on = "set_id", allow.cartesian = TRUE]
possible_edges <- possible_edges[element < i.element] # Keep unique pairs only

message("\n-- All Possible Intra-Set Edges --")
print(possible_edges)

# 1c. Create canonical representation and join with real edges
possible_edges[, `:=`(canon1 = pmin(element, i.element), canon2 = pmax(element, i.element))]
observed_edges_long <- ecg_edges[possible_edges, on = .(canon1, canon2), nomatch = 0]

message("\n-- Observed Intra-Set Edges (the intersection) --")
print(observed_edges_long)

# 1d. Aggregate to get the final count for each set
observed_edges_dt <- observed_edges_long[, .(observed_edges = .N), by = set_id]

message("\n-- Final Observed Edge Counts per Set --")
print(observed_edges_dt)


message("\n\n--- STEP 2: Calculate Lambda Components ---")

# 2a. Create a lookup table for element degrees
all_element_degrees_dt <- data.table(
  element = ecp@names,
  degree = unlist(ecp@degrees)
)
setkey(all_element_degrees_dt, element)

# 2b. Join to get the degree for each element in each set (filtering out unknown elements)
setkey(set_membership_dt_unique, element)
sets_with_degrees <- all_element_degrees_dt[set_membership_dt_unique, on = "element", nomatch = 0]

message("\n-- Set Memberships with Degrees (elements not in graph are dropped) --")
print(sets_with_degrees)

# 2c. Aggregate to get sum of degrees and sum of squared degrees
term_summary <- sets_with_degrees[, .(
  sum_of_degrees = sum(degree, na.rm = TRUE),
  sum_of_sq_degrees = sum(degree^2, na.rm = TRUE),
  set_size = .N
), by = set_id]

message("\n-- Final Summary of Degrees per Set --")
print(term_summary)


message("\n\n--- STEP 3: Final Combination ---")
# This demonstrates the final join and calculation logic from the main function

final_dt <- copy(sets_dt)
final_dt[observed_edges_dt, on = "set_id", observed_edges := i.observed_edges]
final_dt[term_summary, on = "set_id", `:=`(
  sum_of_degrees = i.sum_of_degrees,
  sum_of_sq_degrees = i.sum_of_sq_degrees,
  set_size = i.set_size
)]
final_dt[is.na(observed_edges), observed_edges := 0L]
# ... (rest of the NA cleaning and calculations) ...

message("\n-- Final Combined Table (before stats calculation) --")
print(final_dt)
