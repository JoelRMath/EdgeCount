# EdgeCount

[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-Tutorial-blue)](https://joelrmath.github.io/EdgeCount/articles/EdgeCount-Tutorial.html)

`EdgeCount` provides a statistical framework for enrichment analysis based on graph connectedness. It is specifically optimized for large, sparse graphs, such as protein interaction networks, social networks, or term-element membership graphs (bipartite graphs).

## Methodology

The package combines the **Random Graph with Given Expected Degrees (RGGED)** model with analytical expressions for connectedness effect size and p-values. 

A key feature is the introduction of **Vertex Set Enrichment Analysis (VSEA)**, a permutation-based method inspired by GSEA. Unlike standard enrichment methods, VSEA uniquely accounts for the "promiscuity" of elements by adjusting for the degree of each element (i.e., the number of terms an element belongs to).

## Key Functionality

The primary functions allow users to:

* **Model Validation:** Evaluate a graph's suitability for specific `EdgeCount` statistical models.
* **Internal Connectedness:** Test for significant edge density within one or more vertex sets.
* **External Connectedness:** Test for connectedness between distinct pairs of sets.
* **Scoring:** Score a collection of terms against a single set or a data table of element sets.
* **VSEA:** Perform Vertex Set Enrichment Analysis on a ranked list of elements.
* **Results Summarization:** Calculate Normalized Enrichment Scores (NES) and FDR q-values for VSEA results.

## Installation

Install the development version from GitHub:

```r
# if (!require("devtools")) install.packages("devtools")
devtools::install_github("JoelRMath/EdgeCount")
```

Documentation

For a comprehensive guide and worked examples, please refer to the official documentation:

🚀 EdgeCount Tutorial: A step-by-step walkthrough with real-world examples.

📖 Theory Vignette (PDF): Detailed mathematical background and model specifications.

Primary Publication

If you use EdgeCount in your research, please cite:

J. Pradines, V. Farutin, S. Rowley, and V. Dancik. "Analyzing protein lists with large networks: edge-count probabilities in random graphs with given expected degrees." Journal of Computational Biology, 12(2):113–128, 2005. DOI: 10.1089/cmb.2005.12.113

Maintained by Joel R. Pradines
