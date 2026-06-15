# Companion bar chart for `community_sharing_pq()`

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Computes the same pairwise metrics as
[`community_sharing_pq()`](https://adrientaudiere.github.io/comparpq/reference/community_sharing_pq.md)
and displays them as grouped bars. Useful for precise numerical
comparison alongside the network figure.

## Usage

``` r
community_sharing_barplot_pq(
  physeq,
  fact,
  metrics = default_sharing_metrics(),
  facet_by = c("metric", "pair"),
  show_na_modality = FALSE,
  base_size = 12,
  title = NULL
)
```

## Arguments

- physeq:

  A `phyloseq` object.

- fact:

  Character. Name of a `sample_data(physeq)` column. Must have 2 to 4
  unique values.

- metrics:

  Named list of metric definitions from
  [`make_sharing_metric()`](https://adrientaudiere.github.io/comparpq/reference/make_sharing_metric.md).
  Default:
  [`default_sharing_metrics()`](https://adrientaudiere.github.io/comparpq/reference/default_sharing_metrics.md).

- facet_by:

  Character. One of `"metric"` (default: one panel per metric, x-axis =
  pair label) or `"pair"` (one panel per pair, x-axis = metric label).
  Bar fill colors come from each metric's `color` field.

- show_na_modality:

  Logical. Same meaning as in
  [`community_sharing_pq()`](https://adrientaudiere.github.io/comparpq/reference/community_sharing_pq.md).
  Default `FALSE`.

- base_size:

  Base font size. Default `12`.

- title:

  Plot title. Default `NULL`.

## Value

A `ggplot` object.

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# One panel per metric
community_sharing_barplot_pq(data_fungi, fact = "Height")

# One panel per pair
community_sharing_barplot_pq(data_fungi, fact = "Height", facet_by = "pair")
} # }
```
