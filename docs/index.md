# comparpq

**comparpq** is an extension to the
[{MiscMetabar}](https://adrientaudiere.github.io/MiscMetabar/) and
[{phyloseq}](https://joey711.github.io/phyloseq/) packages that provides
tools for comparing phyloseq objects. The package focuses on taxonomic
comparison, accuracy metrics computation, and interactive visualizations
for microbiome data analysis.

## Installation

You can install the development version of comparpq from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("adrientaudiere/comparpq")
```

## Main Features

### Taxonomic Comparison and Accuracy Metrics

- **[`tc_metrics_mock()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock.md)**:
  Compute accuracy metrics comparing taxonomic assignations against mock
  communities
- **[`tc_metrics_mock_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock_vec.md)**:
  Vectorized version for efficient metric computation
- **[`tc_points_matrix()`](https://adrientaudiere.github.io/comparpq/reference/tc_points_matrix.md)**:
  Create matrices for taxonomic comparison visualization

### Visualization Tools

- **[`bubbles_pq()`](https://adrientaudiere.github.io/comparpq/reference/bubbles_pq.md)**:
  Interactive bubble plots using Observable HQ for phyloseq objects
- **[`tc_bar()`](https://adrientaudiere.github.io/comparpq/reference/tc_bar.md)
  &
  [`tc_circle()`](https://adrientaudiere.github.io/comparpq/reference/tc_circle.md)**:
  Bar charts and circular plots for taxonomic comparisons
- **[`rainplot_taxo_na()`](https://adrientaudiere.github.io/comparpq/reference/rainplot_taxo_na.md)**:
  Rain plots for visualizing NA patterns in taxonomic data

### Data Manipulation

- **Fake taxa creation**:
  [`add_external_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_external_seq_pq.md),
  [`add_shuffle_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_shuffle_seq_pq.md)
  for mock community preparation
- **Taxonomic table utilities**:
  [`rename_ranks_pq()`](https://adrientaudiere.github.io/comparpq/reference/rename_ranks_pq.md),
  [`select_ranks_pq()`](https://adrientaudiere.github.io/comparpq/reference/select_ranks_pq.md),
  [`resolve_taxo_conflict()`](https://adrientaudiere.github.io/comparpq/reference/resolve_taxo_conflict.md)
- **Pattern replacement**:
  [`taxtab_replace_pattern_by_NA()`](https://adrientaudiere.github.io/comparpq/reference/taxtab_replace_pattern_by_NA.md)
  for data cleaning

## Quick Start

``` r
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

- **Package website**: <https://adrientaudiere.github.io/comparpq/>
- **Getting started vignette**: See the package vignettes for
  comprehensive examples
- **Function documentation**: All functions have detailed help pages
  accessible via `?function_name`

## Related Packages

This package extends the functionality of: -
[{MiscMetabar}](https://adrientaudiere.github.io/MiscMetabar/):
Miscellaneous functions for metabarcoding analysis -
[{phyloseq}](https://joey711.github.io/phyloseq/): Handling and analysis
of high-throughput microbiome census data

## Citation

If you use comparpq in your research, please cite:

``` r
citation("comparpq")
```

## License

GPL (\>= 3)

## Issues and Contributions

Please report bugs and feature requests at
<https://github.com/adrientaudiere/comparpq/issues>

## Author

**Adrien Taudière** (aut, cre, cph)  
ORCID: [0000-0003-1088-1182](https://orcid.org/0000-0003-1088-1182)  
Email: <adrien.taudiere@zaclys.net>
