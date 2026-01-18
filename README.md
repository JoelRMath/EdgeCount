# EdgeCount

[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-Tutorial-blue)](https://joelrmath.github.io/EdgeCount/articles/EdgeCount-Tutorial.html)

`EdgeCount` provides a statistical framework for enrichment analysis based on graph connectedness. It is specifically optimized for large, sparse graphs, such as protein interaction networks, social networks, or term-element membership graphs (bipartite graphs).

## 📖 High-Level Overview

**EdgeCount** is an R package designed to uncover statistically significant patterns within large, complex systems.

While the mathematical framework was originally published to solve problems in Computational Biology (analyzing protein interaction networks), the logic is **domain-agnostic**. It solves a universal problem in data science: **distinguishing true signals from noise in highly interconnected data.**

### The Problem: Popularity Bias
In many large datasets—whether social networks, retail transactions, or biological systems—some elements are naturally "louder" than others simply because they are popular (they have a high "node degree").
* In **Social Networks**, an influencer connects to everyone.
* In **Retail Analytics**, common items (like water or batteries) appear in almost every basket.
* In **Biology**, "promiscuous" proteins bind to many partners.

Standard analysis often falsely flags these popular items as "significant" connections.

### The Solution
**EdgeCount** uses a **Random Graph with Given Expected Degrees (RGGED)** model to rigorously correct for this bias. It allows a Data Scientist to ask: *"Is the connection between Group A and Group B actually surprising, or is it just a result of the underlying popularity of the items involved?"*

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
