# Plot MaAsLin3 results

Creates ggplot2 visualizations of MaAsLin3 differential abundance
results. Supports summary plots (default), volcano plots, and forest
plots.

## Usage

``` r
gg_maaslin3_plot(
  res,
  type = c("summary", "volcano", "forest"),
  model = c("abundance", "prevalence", "both"),
  metadata_filter = NULL,
  signif_threshold = 0.1,
  use_qval = TRUE,
  top_n = 25,
  show_labels = FALSE,
  point_size = 2,
  colors = NULL,
  show_annotations = TRUE,
  reference_label = NULL,
  coef_plot_vars = NULL,
  heatmap_vars = NULL,
  normalization = "TSS",
  transform = "LOG"
)
```

## Arguments

- res:

  (list, required) The full result object from
  [`maaslin3_pq()`](https://adrientaudiere.github.io/comparpq/reference/maaslin3_pq.md).
  For types other than "summary", can also be directly the `$results`
  data frame.

- type:

  (character, default "summary") Type of plot. One of:

  - `"summary"`: Default MaAsLin3 summary plot with coefficients and
    heatmap (uses
    [`maaslin3::maaslin_plot_results()`](https://rdrr.io/pkg/maaslin3/man/maaslin_plot_results.html)
    internally)

  - `"volcano"`: Effect size (coef) vs -log10(p-value)

  - `"forest"`: Effect sizes with confidence intervals

- model:

  (character, default "abundance") Which model results to plot. One of
  `"abundance"`, `"prevalence"`, or `"both"`. Used for volcano and
  forest plots.

- metadata_filter:

  (character, default NULL) Filter results to a specific metadata
  variable. If NULL, uses the first non-intercept variable. Used for
  volcano and forest plots.

- signif_threshold:

  (numeric, default 0.1) Significance threshold for q-value
  (FDR-corrected). For summary plot, this is passed to
  `max_significance`. For other plots, used to color significant
  features.

- use_qval:

  (logical, default TRUE) If TRUE, use q-values (FDR-corrected) for
  significance. If FALSE, use raw p-values. Used for volcano and forest.

- top_n:

  (integer, default 25) For summary plot, number of top features to
  display (`summary_plot_first_n`). For forest plot, show only the top N
  features by absolute effect size.

- show_labels:

  (logical, default FALSE) For volcano plots, show feature labels for
  significant points.

- point_size:

  (numeric, default 2) Size of points in volcano plot.

- colors:

  (character vector, default NULL) Colors for significant/
  non-significant points. If NULL, uses c("grey60", "firebrick3").

- show_annotations:

  (logical, default TRUE) If TRUE, adds annotations indicating which
  group has higher values on each side of the plot. For simple formulas
  like `~ Height`, shows "More in High" vs "More in Low".

- reference_label:

  (character, default NULL) Label for the reference group. If NULL,
  attempts to extract from the metadata column name.

- coef_plot_vars:

  (character vector, default NULL) For summary plot only. Variables to
  include in coefficient plot. If NULL, uses all non-intercept
  variables.

- heatmap_vars:

  (character vector, default NULL) For summary plot only. Variables to
  include in heatmap. If NULL, determined automatically.

- normalization:

  (character, default "TSS") Normalization method used. Only needed for
  summary plot. One of "TSS", "CLR", or "NONE".

- transform:

  (character, default "LOG") Transform method used. Only needed for
  summary plot. One of "LOG", "PLOG", or "NONE".

## Value

A ggplot2 object.

## Details

\#TODO VERY experimental
[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

**Summary plot** (default): Uses
[`maaslin3::maaslin_plot_results()`](https://rdrr.io/pkg/maaslin3/man/maaslin_plot_results.html)
to create the standard MaAsLin3 visualization with sorted per-feature
coefficients plotted with standard errors, and a heatmap for additional
variables. This requires the full maaslin3 result object (not just the
results data frame).

**Volcano plot**: Shows the relationship between effect size
(coefficient) and statistical significance. Points above the horizontal
dashed line are significant at the specified threshold. Points are
colored by significance.

**Forest plot**: Shows effect sizes with 95 percent confidence intervals
for the top features. Features are ordered by effect size. Significant
associations are highlighted.

The `coef` in MaAsLin3 abundance models represents log2 fold change: a
one-unit change in the variable corresponds to a 2^coef fold change in
relative abundance.

## See also

[`maaslin3_pq()`](https://adrientaudiere.github.io/comparpq/reference/maaslin3_pq.md),
[`gg_aldex_plot()`](https://adrientaudiere.github.io/comparpq/reference/gg_aldex_plot.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# Run MaAsLin3
res <- maaslin3_pq(
  data_fungi,
  formula = "~ Height",
  reference = list(Height = "Low")
)

# Summary plot (default) - uses maaslin3's native visualization
gg_maaslin3_plot(res)

# Volcano plot
gg_maaslin3_plot(res, type = "volcano")

# Forest plot of top 15 features
gg_maaslin3_plot(res, type = "forest", top_n = 15)

# Plot prevalence model results
gg_maaslin3_plot(res, type = "volcano", model = "prevalence")

# Customize significance threshold
gg_maaslin3_plot(res, type = "volcano", signif_threshold = 0.05, show_labels = TRUE)

 # Complete Example with HMP2 Dataset

 library(maaslin3)
 # ============================================================
 # Convert maaslin3 default HMP2 dataset to phyloseq
 # ============================================================

 # Read the HMP2 default dataset from maaslin3 package
 taxa_table_name <- system.file("extdata", "HMP2_taxonomy.tsv", package = "maaslin3")
 taxa_table <- read.csv(taxa_table_name, sep = "\t", row.names = 1)

 metadata_name <- system.file("extdata", "HMP2_metadata.tsv", package = "maaslin3")
 metadata <- read.csv(metadata_name, sep = "\t", row.names = 1)

 # Set factor levels
 metadata$antibiotics <- factor(metadata$antibiotics, levels = c("No", "Yes"))

 # Create phyloseq object (HMP2 data has samples as rows, taxa as columns)
 otu <- otu_table(as.matrix(taxa_table), taxa_are_rows = FALSE)
 sam <- sample_data(metadata)
 species_names <- colnames(taxa_table)
 tax_df <- data.frame(
   Species = species_names,
   Genus = sapply(strsplit(species_names, "_"), \(x) x[1]),
   row.names = species_names
 )
 tax <- tax_table(as.matrix(tax_df))
 physeq_hmp2 <- phyloseq(otu, sam, tax)

 # ============================================================
 # Run MaAsLin3 analysis
 # ============================================================

 res <- maaslin3_pq(
   physeq_hmp2,
   formula = "~ antibiotics",
   reference = list(antibiotics = "No"),
   output = "output/maaslin3_hmp2",
   correction_for_sample_size = FALSE
 )

 # ============================================================
 # Compare plots
 # ============================================================

 # Summary plot (NEW DEFAULT) - uses maaslin3's native visualization
 gg_maaslin3_plot(res)

 # Alternative plot types
 gg_maaslin3_plot(res, type = "volcano")
 gg_maaslin3_plot(res, type = "forest", top_n = 15)
} # }
```
