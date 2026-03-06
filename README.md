# comparpq

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN status](https://www.r-pkg.org/badges/version/comparpq)](https://CRAN.R-project.org/package=comparpq)
<!-- badges: end -->

**comparpq** is an extension to the [{MiscMetabar}](https://adrientaudiere.github.io/MiscMetabar/) and [{phyloseq}](https://joey711.github.io/phyloseq/) packages that provides tools for comparing phyloseq objects. The package focuses on taxonomic comparison, accuracy metrics computation, and interactive visualizations for microbiome data analysis.

## Installation

You can install the development version of comparpq from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("adrientaudiere/comparpq")
```

## Main Features

### Compare phyloseq objects

- **`list_phyloseq`**: S7 class for storing and comparing multiple phyloseq objects, with automatic detection of comparison type (REPRODUCIBILITY, ROBUSTNESS, NESTED_ROBUSTNESS, REPLICABILITY, EXPLORATION, SEPARATE_ANALYSIS)
- **`compare_refseq()`**: Compare reference sequences between two phyloseq objects, identifying shared/unique ASVs by name and by DNA sequence content
- **`filter_common_lpq()`**: Filter phyloseq objects to keep only shared samples and/or taxa
- **`apply_to_lpq()`**: Apply any function to all phyloseq objects in a `list_phyloseq`

### Statistical Analysis

- **`adonis_lpq()`**: PERMANOVA analysis across multiple phyloseq objects
- **`estim_diff_pq()` & `estim_diff_lpq()`**: Estimation statistics for categorical comparisons using bootstrap confidence intervals (via [dabestr](https://cran.r-project.org/package=dabestr))
- **`estim_cor_pq()` & `estim_cor_lpq()`**: Estimation statistics for continuous correlations with bootstrap confidence intervals
- **`ancombc_lpq()`**, **`aldex_lpq()`**, **`multipatt_lpq()`**, **`maaslin3_pq()`**: Differential abundance analysis methods

### Taxonomic Comparison and Accuracy Metrics
- **`tc_metrics_mock()`**: Compute accuracy metrics comparing taxonomic assignations against mock communities
- **`tc_metrics_mock_vec()`**: Vectorized version for efficient metric computation
- **`tc_points_matrix()`**: Create matrices for taxonomic comparison visualization

### Visualization Tools
- **`bubbles_pq()`**: Interactive bubble plots using Observable HQ for phyloseq objects
- **`tc_bar()` & `tc_circle()`**: Bar charts and circular plots for taxonomic comparisons
- **`tc_sankey()`**: Sankey diagrams for taxonomic rank relationships
- **`tc_linked_trees()`**: Side-by-side linked taxonomy trees for comparing two assignments
- **`tc_congruence_metrics()`**: Congruence metrics between two taxonomic assignments
- **`rainplot_taxo_na()`**: Rain plots for visualizing NA patterns in taxonomic data
- **`formattable_lpq()` & `formattable_lpq_full()`**: Formatted summary tables with color bars for `list_phyloseq`
- **`upset_lpq()`**: UpSet plots for visualizing sample/taxa overlap across phyloseq objects
- **`gg_aldex_plot()` & `gg_maaslin3_plot()`**: ggplot2-based visualizations for differential abundance results

### Data Manipulation
- **Mock community preparation**: `add_external_seq_pq()`, `add_shuffle_seq_pq()`, `multiply_counts_pq()`, `permute_da_pq()`, `midasim_pq()` for creating and manipulating fake taxa to benchmark taxonomic assignment or differential abundance methods
- **Taxonomic table utilities**: `rename_ranks_pq()`, `select_ranks_pq()`, `resolve_taxo_conflict()`, `taxtab_replace_pattern_by_NA()`
- **Taxonomy to tree**: `taxo2tree()` to build a phylogenetic tree from taxonomy table

## Quick Start

```r
library(comparpq)

# Load example data
data("Glom_otu")

# Create an interactive bubble plot
bubbles_pq(Glom_otu, 
          rank_color = "Family", 
          min_nb_seq = 1000)

# For taxonomic comparison analysis with mock communities:
# 1. First prepare your data with fake taxa
physeq_with_fake <- add_shuffle_seq_pq(Glom_otu, nb_seq_fake = 50)

# 2. Perform taxonomic assignment (using your preferred method)
# 3. Compare results against known taxonomy
# tc_metrics_mock(physeq_with_fake, ranks_df, true_values_df)
```

## Documentation

- **Package website**: https://adrientaudiere.github.io/comparpq/
- **Getting started vignette**: See the package vignettes for comprehensive examples
- **Function documentation**: All functions have detailed help pages accessible via `?function_name`

## Related Packages

This package extends the functionality of:
- [{MiscMetabar}](https://adrientaudiere.github.io/MiscMetabar/): Miscellaneous functions for metabarcoding analysis
- [{phyloseq}](https://joey711.github.io/phyloseq/): Handling and analysis of high-throughput microbiome census data

## Citation

If you use comparpq in your research, please cite:

```r
citation("comparpq")
```

## License

GPL (>= 3)

## Issues and Contributions

Please report bugs and feature requests at https://github.com/adrientaudiere/comparpq/issues

## Author

**Adrien Taudière** (aut, cre, cph)  
ORCID: [0000-0003-1088-1182](https://orcid.org/0000-0003-1088-1182)  
Email: adrien.taudiere@zaclys.net

