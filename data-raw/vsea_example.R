library(EdgeCount)
library(data.table)

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

get_run_time <- function(ects, element_ranks, n_replicates, n_permutations){

  run_times <- replicate(n_replicates, {
    system.time({
      vsea_results <- run_vsea_analysis(ects, scoring_statistic = "log2_Anscombe_ratio", element_ranks = element_ranks, n_permutations = n_permutations)
    })[["elapsed"]]
  })
  mean(run_times)
}

data("sample_ects")
term_dt <- data.table(term = sample_ects@terms,
                                term_degree = unlist(sample_ects@ecprob@degrees[sample_ects@terms]))
nelogne <- NULL
sktlogkt <- NULL
run_time <- NULL
term_size_min <- c(2:14)
for (term_size in term_size_min){
  message(term_size)
  term_selection_dt <- term_dt[term_degree >= term_size]
  selected_terms <- unlist(term_selection_dt$term)
  ects <- reduce_universe_by_terms(sample_ects, selected_terms)
  lst <- get_complexity_metric(ects)
  nelogne <- c(nelogne, lst$ne_log_ne)
  sktlogkt <- c(sktlogkt, lst$sum_kt_log_kt)
  elements <- ects@elements
  element_ranks <- setNames(sample(1:length(elements)), elements)
  time_to_run <- get_run_time(ects, element_ranks = element_ranks, n_permutations = 20, n_replicates = 10)
  run_time <- c(run_time, time_to_run)
}

df <- data.frame(nelogne = nelogne,
                 sktlogkt = sktlogkt,
                 run_time = run_time)
df <- df[order(df$nelogne),]
complexity <- nelogne + sktlogkt
plot(sktlogkt, run_time, xlim = c(0, max(c(nelogne, sktlogkt))))
points(nelogne, run_time)
model_multiple <- lm(
  run_time ~ sktlogkt + nelogne,
  data = df
)
print(summary(model_multiple))
# plot(model_multiple)
print(df)
