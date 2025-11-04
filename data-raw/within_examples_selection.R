# This script is for selecting compelling "opposite pairs" of GO terms
# for the tutorial. The goal is to find pairs of terms that have the
# same size and same observed edge count, but widely different
# expected edge counts (lambda), which highlights the importance of
# vertex degrees.

library(data.table)

# 1. Load the data
file_path <- "data-raw/res/go_term_in_stats.tsv"
if (!file.exists(file_path)) {
  stop(paste("File not found. Please run the 'data-raw/in_examples.R' script first to generate:", file_path))
}
all_stats <- fread(file_path)

# 2. Initial Filtering
# We only want to compare terms where there is at least some connectivity
filtered_stats <- all_stats[observed_edges >= 3]

# 3. Group by (set_size, observed_edges) and find min/max lambda stats
#    For each group, we find the terms with the min and max lambda,
#    and the count of terms in that group.
summary_dt <- filtered_stats[, {
  # Find the indices (row numbers) of the min and max lambda
  idx_min <- which.min(lambda)
  idx_max <- which.max(lambda)

  # Use the indices to get all related stats
  .(
    count_in_group = .N,

    min_lambda = lambda[idx_min],
    set_id_at_min = set_id[idx_min], # <-- ADDED
    term_at_min = term_name[idx_min],
    pval_at_min = p_value[idx_min],
    log2AR_at_min = log2_Anscombe_ratio[idx_min],

    max_lambda = lambda[idx_max],
    set_id_at_max = set_id[idx_max], # <-- ADDED
    term_at_max = term_name[idx_max],
    pval_at_max = p_value[idx_max],
    log2AR_at_max = log2_Anscombe_ratio[idx_max]
  )
}, by = .(set_size, observed_edges)]


# 4. Filter for groups that have at least two terms to compare
summary_dt <- summary_dt[count_in_group > 1]

# 5. Calculate the lambda ratio to find the most "opposite" pairs
#    We handle the (unlikely) case of min_lambda being 0 to avoid Inf.
summary_dt[, lambda_ratio := ifelse(
  min_lambda > 0,
  max_lambda / min_lambda,
  Inf
)]

# 6. Sort by the most extreme ratios first
setorder(summary_dt, -lambda_ratio)

# 7. Save the full results table
output_dir <- "data-raw/res"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
output_file <- file.path(output_dir, "within_examples_selection.tsv")
fwrite(summary_dt, output_file, quote = FALSE, sep = "\t")

# 8. Print the top 20 most illustrative pairs
message("--- Top 20 Most 'Opposite' Term Pairs ---")
message("(Pairs with same size and observed edges, but max difference in lambda)")
message(paste("Full results saved to:", output_file))
print(head(summary_dt, 20))
