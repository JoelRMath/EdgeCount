library(EdgeCount)
library(data.table)

terms_ecranks_statistics_for_complexity <-  function(object, element_ranks) {

  # --- Step 1: Condition the analysis on the provided ranked list ---
  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) {
    stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
  }
  object <- reduce_universe_by_elements(object, valid_elements)
  element_ranks <- rank(element_ranks[valid_elements])

  # start of core computation #
  start_time <- Sys.time()

  # --- Step 2: Pre-computation of core data structures ---
  all_element_degrees <- unlist(object@ecprob@degrees)
  ranks_dt <- data.table(element_id = names(element_ranks), global_rank = element_ranks)
  setorder(ranks_dt, global_rank)
  ranks_dt[, degree := all_element_degrees[element_id]]
  ranks_dt[, cumsum_degrees := cumsum(degree)]

  bipartite_edges <- as.data.table(to_dataframe(object))
  setnames(bipartite_edges, c("term", "element"), c("term_id", "element_id"))
  bipartite_edges[, term_id := as.character(term_id)]

  # --- Step 3: Join all information into a single flat table ---
  setkey(bipartite_edges, element_id)
  setkey(ranks_dt, element_id)
  final_dt <- ranks_dt[bipartite_edges, on = "element_id"]

  setorder(final_dt, term_id, global_rank)
  final_dt[, rank_in_term := 1:.N, by = term_id]

  term_sizes <- lengths(object@ecprob@adj[object@terms])
  term_summary <- data.table(
    term_id = names(term_sizes),
    term_size = term_sizes,
    term_degree = term_sizes
  )
  final_dt[term_summary, on = "term_id", `:=`(term_size = i.term_size, term_degree = i.term_degree)]

  # --- Step 4: Final, fully vectorized statistics calculation ---
  final_dt[, `:=`(
    observed_ec = rank_in_term,
    max_ec = pmin(term_size, global_rank),
    lambda = (term_degree / (2 * object@ecprob@graph_size)) * cumsum_degrees
  )]

  final_dt[, `:=`(
    p_value = calculate_p_value(object@ecprob, observed_ec, max_ec, lambda),
    log2_Anscombe_ratio = 0.5 * (log2(observed_ec + 3/8) - log2(lambda + 3/8))
  )]

  # end of core computation #
  end_time <- Sys.time()
  run_time <- as.numeric(end_time - start_time, units = "secs")

  # --- STEP 5: Reshape output to the required named list format ---
  results_list <- split(final_dt, by = "term_id")

  # Define and reorder columns for final output
  final_cols <- c("element_id", "global_rank", "rank_in_term",
                  "observed_ec", "max_ec", "term_size",
                  "lambda", "p_value", "log2_Anscombe_ratio")

  results_list_final <- lapply(results_list, function(dt) {
    # It is safer to create the final table with a fresh subset
    dt[, ..final_cols]
  })

  results_list_final[["run_time"]] <- run_time

  return(results_list_final)
}

get_complexity_metric <- function(ects){

  ne <- length(ects@elements)
  ne_log_ne <- ne*log(ne)
  kt <- lengths(ects@ecprob@adj[ects@terms])
  kt_log_kt <- kt*log(kt)
  sum_kt_log_kt <- sum(kt_log_kt)
  return(list(
    ne_log_ne = ne_log_ne,
    sum_kt_log_kt = sum_kt_log_kt
  ))
}

data("sample_ects")
all_terms <- sample_ects@terms
term_dt <- data.table(term = all_terms,
                      term_degree = unlist(sample_ects@ecprob@degrees[all_terms]))


term_size_min <- c(5:14)
n_repeats <- 20
ne_log_ne <- NULL
sum_kt_log_kt <- NULL
run_time <- NULL
n_elements <- seq(from = 500, to = 5000, by = 500)
for (n in n_elements){
  message(n)
  elements <- sample(sample_ects@elements, n)
  ects <- reduce_universe_by_elements(sample_ects, elements)
  metrics <- get_complexity_metric(ects)
  mean_run_time <- 0.
  for (i in 1:n_repeats){
    element_ranks <- setNames(sample(1:length(elements)), elements)
    result <- terms_ecranks_statistics_for_complexity(ects, element_ranks)
    mean_run_time <- mean_run_time + result[["run_time"]]
  }
  mean_run_time <- mean_run_time / n_repeats
  ne_log_ne <- c(ne_log_ne, metrics[["ne_log_ne"]])
  sum_kt_log_kt <- c(sum_kt_log_kt, metrics[["sum_kt_log_kt"]])
  run_time <- c(run_time, mean_run_time)
}
df <- data.frame(run_time = run_time,
                 ne_log_ne = ne_log_ne,
                 sum_kt_log_kt = sum_kt_log_kt)
model_multiple <- lm(
  run_time ~ sum_kt_log_kt + ne_log_ne,
  data = df
)
print(summary(model_multiple))
plot(df$sum_kt_log_kt, df$run_time)
