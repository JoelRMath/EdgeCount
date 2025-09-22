library(data.table)
library(GO.db)

# --- Define some high-level categories to look for ---
# We use grep to find them, so these are patterns.
category_map <- list(
  "Metabolism" = "metabolic process",
  "Signaling" = "signaling",
  "Development" = "developmental process",
  "Immune System" = "immune system process",
  "Complex/Subunit" = "complex|subunit",
  "Localization" = "localization",
  "Transport" = "transport",
  "Regulation" = "regulation of"
)

# --- Function to find all ancestors of a GO term ---
get_ancestors <- function(goid) {
  # This uses the specific parent maps from GO.db
  ancestors <- c(
    if(exists(goid, GOBPANCESTOR)) as.list(GOBPANCESTOR[goid])[[1]],
    if(exists(goid, GOMFANCESTOR)) as.list(GOMFANCESTOR[goid])[[1]],
    if(exists(goid, GOCCANCESTOR)) as.list(GOCCANCESTOR[goid])[[1]]
  )
  return(unique(ancestors[ancestors != "all"]))
}

# --- Load your pre-computed summary file ---
summary_dt <- fread("data-raw/term_pair_connectivity_summary.tsv")

# --- Find a category for each pair ---
categories <- sapply(1:nrow(summary_dt), function(i) {
  term1 <- summary_dt$term1[i]
  term2 <- summary_dt$term2[i]

  # Get all ancestors for both terms
  ancestors1 <- get_ancestors(term1)
  ancestors2 <- get_ancestors(term2)

  # Find the common ancestors
  common_ancestors <- intersect(ancestors1, ancestors2)

  if (length(common_ancestors) > 0) {
    # Get the names of the common ancestors
    ancestor_names <- select(GO.db, keys = common_ancestors, columns = "TERM", keytype = "GOID")$TERM

    # Check if any ancestor name matches our predefined categories
    for (cat_name in names(category_map)) {
      if (any(grepl(category_map[[cat_name]], ancestor_names, ignore.case = TRUE))) {
        return(cat_name) # Return the first matching category
      }
    }
  }

  return("Other") # Default category if no match is found
})

# Add the new category column to your data
summary_dt$category <- categories

# Save the new, categorized file
fwrite(summary_dt, "data-raw/term_pairs_categorized.tsv", sep = "\t")

message("Categorization complete. New file saved.")
