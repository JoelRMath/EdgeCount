library(EdgeCount)
library(data.table)

terms_ecranks_statistics_gs_loop <- function(object, element_ranks) {

  valid_elements <- intersect(names(element_ranks), object@elements)
  if (length(valid_elements) < 1) {
    stop("None of the elements in `element_ranks` are in the ECTermScoring object.")
  }
  object <- reduce_universe_by_elements(object, valid_elements)
  element_ranks <- rank(element_ranks[valid_elements]) # Re-rank from 1 to N_valid

  ecprob <- object@ecprob
  M_g <- ecprob@graph_size
  N_e <- length(object@elements)
  all_terms <- object@terms

  df <- data.frame(element = names(element_ranks),
                   rank = as.numeric(element_ranks),
                   stringsAsFactors = FALSE)
  df <- df[order(df$rank),]
  ranked_elements_sorted <- df$element

  K <- unlist(ecprob@degrees[ranked_elements_sorted])
  cumul_sum_K <- cumsum(K)

  score_one_term_safe <- function(term_id) {

    elements_in_term <- get_neighbors(ecprob, term_id)

    elements_in_term <- intersect(elements_in_term, names(element_ranks))

    if (length(elements_in_term) == 0) return(NULL)

    ranks_for_term <- element_ranks[elements_in_term]

    ranks_for_term_sorted <- sort(ranks_for_term)
    elements_term_sorted <- names(ranks_for_term_sorted)

    sz <- length(elements_term_sorted)

    observed_ec_term <- 1:sz
    term_size <- rep(sz, sz)
    max_ec_term <- pmin(term_size, ranks_for_term_sorted)
    term_degree <- rep(ecprob@degrees[[term_id]], sz)

    lambda_term <- (term_degree / (2 * M_g)) * cumul_sum_K[ranks_for_term_sorted]

    p_value_term <- calculate_p_value(ecprob, observed_ec_term, max_ec_term, lambda_term)
    log2_anscombe_term <- 0.5 * (log2(observed_ec_term + 3/8) - log2(lambda_term + 3/8))
    log2_rel_change_term <- log2(observed_ec_term) - log2(lambda_term)

    df <- data.frame(
      element_id = elements_term_sorted,
      global_rank = element_ranks[elements_term_sorted],
      rank_in_term = rank(ranks_for_term_sorted),
      observed_ec = observed_ec_term,
      max_ec = max_ec_term,
      term_size = term_size,
      lambda = lambda_term,
      p_value = p_value_term,
      log2_Anscombe_ratio = log2_anscombe_term,
      stringsAsFactors = FALSE
    )
    row.names(df) <- NULL
    return(as.data.table(df))
  }

  all_results_list <- lapply(all_terms, score_one_term_safe)
  names(all_results_list) <- all_terms

  all_results_list <- all_results_list[!sapply(all_results_list, is.null)]

  return(all_results_list)
}

set.seed(3)
data("sample_ects")
ects <- sample_ects
elements <- ects@elements
elements <- sample(elements)
element_ranks <- setNames(sample(1:length(elements)), elements)

success <- TRUE
n_tries <- 10

for (i in 1:n_tries){

  element_ranks <- setNames(sample(1:length(elements)), elements)

  message("-- running gs ",i,"/",n_tries," --")
  result_gs <- terms_ecranks_statistics_gs_loop(ects, element_ranks)
  message("-- running s4 ",i,"/",n_tries," --")
  result_s4 <- terms_ecranks_statistics(ects, element_ranks)

  sorted_slow <- result_gs[sort(names(result_gs))]
  sorted_fast <- result_s4[sort(names(result_s4))]

  if (!isTRUE(all.equal(sorted_slow, sorted_fast))){
    success <- FALSE
  }
}
message("success = ", success)

