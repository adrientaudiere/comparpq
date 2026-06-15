# Default metrics for `community_sharing_pq()`

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Returns the 4 default metrics: shared species count, Bray-Curtis
similarity, Jaccard binary similarity, and proportion of shared genera.

`bray_sim`, `jac_sim`, and `genus_prop` use `bounds = c(0, 1)` for
globally consistent linewidth scaling. `shared_sp` uses `bounds = NULL`
(count metric without a fixed upper bound; rescaled within the observed
range).

## Usage

``` r
default_sharing_metrics()
```

## Value

A named list of metric definitions.

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# Default set of 4 metrics
mets <- default_sharing_metrics()
community_sharing_pq(data_fungi, fact = "Height", metrics = mets)

# Extend defaults with a robust Aitchison similarity metric
mets_plus <- c(
  default_sharing_metrics(),
  list(
    rob_aitch = make_sharing_metric(
      label = "Robust Aitchison similarity",
      color = "#FF7F00",
      prep  = function(physeq, fact, modalities) {
        otu_sp <- .agg_by_mod(physeq, fact, modalities)
        as.matrix(vegan::vegdist(t(otu_sp), method = "robust.aitchison"))
      },
      fn    = \(a, b, otu_sp, cache) 1 - cache[a, b]
    )
  )
)
community_sharing_pq(data_fungi, fact = "Height", metrics = mets_plus)

# Keep only a subset of the defaults
community_sharing_pq(
  data_fungi,
  fact    = "Height",
  metrics = default_sharing_metrics()[c("shared_sp", "bray_sim")]
)
} # }
```
