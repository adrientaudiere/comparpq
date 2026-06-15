# Community sharing plot: modalities as pie nodes with multi-metric links

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Draws a figure with one node per modality of `fact` (2 to 4 supported),
positioned on a regular polygon. Each node is a pie chart showing the
taxonomic composition at rank `pie_taxrank`. Between each pair of nodes,
one curved link per metric is drawn; linewidth is rescaled within each
metric to `linewidth_range`. A legend below the figure gives each
metric's min / mean / max over all pairs.

## Usage

``` r
community_sharing_pq(
  physeq,
  fact,
  metrics = default_sharing_metrics(),
  pie_taxrank = "Class",
  pie_r = 0.28,
  label_offset = 0.18,
  curvature_range = c(-0.35, 0.35),
  linewidth_range = c(0.8, 4.5),
  n_perm = 0,
  sig_threshold = 0.05,
  nonsig_alpha = 0.15,
  seed = NULL,
  palette = "Set3",
  max_taxa = 12,
  other_color = "grey70",
  show_na_modality = FALSE,
  show_na = TRUE,
  na_color = "grey40",
  pie_border_color = "black",
  pie_border_width = 0.6,
  title = NULL,
  base_size = 12
)
```

## Arguments

- physeq:

  A `phyloseq` object.

- fact:

  Character. Name of a `sample_data(physeq)` column used to group
  samples into modalities. Must have 2 to 4 unique values.

- metrics:

  Named list of metric definitions. See
  [`default_sharing_metrics()`](https://adrientaudiere.github.io/comparpq/reference/default_sharing_metrics.md)
  and
  [`make_sharing_metric()`](https://adrientaudiere.github.io/comparpq/reference/make_sharing_metric.md)
  for the format. Default uses 4 metrics: shared species, Bray-Curtis
  similarity, Jaccard similarity, and proportion of common genera.

- pie_taxrank:

  Character. Taxonomic rank for the pie charts. Default `"Class"`.

- pie_r:

  Numeric. Pie radius in data units. Default `0.28`.

- label_offset:

  Numeric. Distance between pie edge and node label. Default `0.18`.

- curvature_range:

  Numeric of length 2. Range of `geom_curve` curvatures used to fan out
  metrics between a pair of nodes. Default `c(-0.35, 0.35)`. The sign is
  flipped per pair so that the first metric always fans toward the plot
  centre (and the last toward the border), regardless of which side of
  the plot the pair is on.

- linewidth_range:

  Numeric of length 2. `linewidth_range[1]` (minimum linewidth)
  corresponds to the weakest similarity or sharing; `linewidth_range[2]`
  (maximum linewidth) corresponds to the strongest. All built-in metrics
  are expressed as similarities (higher value = more similar communities
  = thicker link). For bounded metrics (`bounds = c(0, 1)`), the mapping
  is global across runs; for unbounded metrics it is scaled within the
  observed range. Default `c(0.8, 4.5)`.

- n_perm:

  Integer. Number of label-permutation iterations for significance
  testing. Default `0` (no test). See Details.

- sig_threshold:

  Numeric in `(0, 1)`. p-value threshold below which a link is
  considered significant. Default `0.05`. Only used when `n_perm > 0`.

- nonsig_alpha:

  Numeric in `[0, 1]`. Alpha (transparency) applied to non-significant
  links. Default `0.15`. Significant links keep alpha `0.75`. Only used
  when `n_perm > 0`.

- seed:

  Integer or `NULL`. Passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the
  permutation loop for reproducibility. Default `NULL` (no seed).

- palette:

  Either the name of an RColorBrewer qualitative palette (e.g. `"Set3"`)
  or a character vector of fill colors used for the top taxa.

- max_taxa:

  Integer. Keep colors for the `max_taxa` most abundant taxa (summed
  across all modalities) and collapse the rest into an `"Other"`
  category filled with `other_color`. Default `12`.

- other_color:

  Fill color for the `"Other"` category. Default `"grey70"`.

- show_na_modality:

  Logical. If `TRUE`, samples whose `fact` value is `NA` are grouped
  into an extra `"NA"` modality and displayed as an additional node. The
  total number of modalities (including `"NA"`) must still be between 2
  and 4. Default `FALSE`.

- show_na:

  Logical. If `TRUE` (default), taxa with an `NA`/empty value at
  `pie_taxrank` are shown as a dedicated `"NA"` category filled with
  `na_color`. If `FALSE`, they are dropped.

- na_color:

  Fill color for the `"NA"` category. Default `"grey40"`.

- pie_border_color:

  Color of the circle drawn around each pie. Default `"black"`.

- pie_border_width:

  Linewidth of the pie border. Default `0.6`.

- title:

  Plot title (optional).

- base_size:

  Base font size. Default `12`.

## Value

A `ggplot` object.

## Details

**Permutation null model.** When `n_perm > 0`, significance is assessed
by label permutation: the `fact` column in `sample_data(physeq)` is
shuffled uniformly at random among samples (preserving group sizes),
then all `prep()` functions and metric computations are re-run on the
permuted data. This is repeated `n_perm` times. The empirical p-value
for each (pair, metric) combination is the proportion of permuted values
greater than or equal to the observed value (one-sided upper-tail test).
This null model is appropriate for all built-in metrics because they
compare aggregated community profiles between groups, and the
aggregation depends on which samples are assigned to each group.
Non-significant links are drawn with alpha `nonsig_alpha` (faded);
significant links use alpha `0.75`.

**Performance.** Each permutation re-runs all `prep()` functions. For
metrics that call
[`phyloseq::tax_glom()`](https://rdrr.io/pkg/phyloseq/man/tax_glom.html)
(e.g. `genus_prop`), this can be slow on large datasets. Recommended
range: `n_perm = 99` to `n_perm = 199`.

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: 4 metrics, pie charts at Class rank
community_sharing_pq(data_fungi, fact = "Height")

# Custom metrics: add robust Aitchison similarity (vegan >= 2.6)
mets <- c(
  default_sharing_metrics()[1],
  list(
    rob_aitch = make_sharing_metric(
      label = "Robust Aitchison similarity",
      color = "#FF7F00",
      prep  = function(physeq, fact, modalities) {
        otu_sp <- .agg_by_mod(physeq, fact, modalities)
        as.matrix(1 - vegan::vegdist(t(otu_sp), method = "robust.aitchison"))
      },
      fn    = \(a, b, otu_sp, cache) 1 - cache[a, b]
    )
  )
)
community_sharing_pq(data_fungi, fact = "Height", metrics = mets)

# Only show a single metric (Jaccard)
community_sharing_pq(
  data_fungi,
  fact    = "Height",
  metrics = default_sharing_metrics()["jac_sim"]
)

# Use Phylum rank for the pie charts and a different palette
community_sharing_pq(
  data_fungi,
  fact        = "Height",
  pie_taxrank = "Phylum",
  palette     = "Set2",
  max_taxa    = 8
)

# Hide NA taxa in pies and add a title
community_sharing_pq(
  data_fungi,
  fact    = "Height",
  show_na = FALSE,
  title   = "Community sharing across heights"
)

# Include samples with NA height as a fourth modality
community_sharing_pq(
  data_fungi,
  fact             = "Height",
  show_na_modality = TRUE
)

# Permutation significance test (99 permutations, faded non-sig links)
community_sharing_pq(
  data_fungi,
  fact   = "Height",
  n_perm = 99,
  seed   = 42
)
} # }
```
