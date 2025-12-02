#' Sample Gene Symbols Lookup
#'
#' A named character vector mapping Entrez IDs to Gene Symbols.
#' @format A named character vector.
"sample_gene_symbols"

#' Sample Term Names Lookup
#'
#' A named character vector mapping Term IDs (GO IDs) to descriptive names.
#' @format A named character vector.
"sample_term_names"

#' Sample Protein Interaction Network
#'
#' A curated subset of the human STRING network (version 12.0) with interaction
#' scores of 900 or higher. The network has been trimmed to create a smaller,
#' more manageable example for the package's vignettes and tests.
#'
#' @format An object of class \code{\link{ECGraph}} with 7244 nodes (elements)
#'   and 286,341 edges.
#' @source \url{https://string-db.org/}
#' Data from STRING is licensed under CC BY 4.0:
#' \url{https://creativecommons.org/licenses/by/4.0/}
"sample_ecg"

#' Sample GO Term to Gene Mappings
#'
#' A curated subset of Gene Ontology (GO) term-to-gene mappings. The data has
#' been trimmed and filtered to share a common element universe with the
#' `sample_ecg` dataset.
#'
#' @format An object of class \code{\link{ECTermScoring}} with 7244 elements
#'   and 9796 terms.
#' @source Gene Ontology Consortium (\url{http://geneontology.org/}) and org.Hs.eg.db.
#' GO data is licensed under CC BY 4.0:
#' \url{https://creativecommons.org/licenses/by/4.0/}
"sample_ects"

#' GO Term Name Lookup Table
#'
#' A named character vector that maps GO term IDs (e.g., "GO:0005575") to their
#' corresponding descriptive names (e.g., "cellular component"). This is a
#' helper dataset for annotating results.
#'
#' @format A named character vector. The names are GO IDs, and the values are
#'   the term names.
#' @source GO.db Bioconductor package
"sample_term_lookup"

#' Entrez Gene Symbol Lookup Table
#'
#' A named character vector that maps Entrez Gene IDs (e.g., "7157") to their
#' corresponding official gene symbols (e.g., "TP53"). This is a helper dataset
#' for annotating results.
#'
#' @format A named character vector. The names are Entrez IDs, and the values are
#'   the gene symbols.
#' @source org.Hs.eg.db Bioconductor package
"sample_gene_symbol_lookup"
