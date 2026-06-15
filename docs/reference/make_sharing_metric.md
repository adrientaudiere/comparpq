# Build a single metric definition for `community_sharing_pq()`

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Usage

``` r
make_sharing_metric(
  label,
  color,
  fn,
  fmt = "%.2f",
  prep = NULL,
  bounds = NULL
)
```

## Arguments

- label:

  Character. Human-readable label used in the legend.

- color:

  Character. Color of the curve connecting modalities.

- fn:

  Function with signature `function(a, b, otu_sp, cache)` returning a
  single numeric value. `a`, `b` are modality names; `otu_sp` is the
  species-level OTU table aggregated by modality; `cache` is the object
  returned by `prep` (or `NULL`).

- fmt:

  `sprintf` format string for legend stats (min / mean / max). Default
  `"%.2f"`; use `"%.0f"` for integer counts.

- prep:

  Optional. Function `function(physeq, fact, modalities)` called once to
  precompute helper data shared across all pairs (e.g. genus-level
  aggregation). Returns an object passed as `cache` to `fn`.

- bounds:

  Either `NULL` (default) or a numeric vector `c(lower, upper)` giving
  the theoretical range of the metric. When non-`NULL`, linewidth is
  scaled to `linewidth_range` using these fixed bounds, so the same
  metric value always maps to the same linewidth regardless of the
  observed data. When `NULL`, linewidth is rescaled within the observed
  range of values. Use `bounds = c(0, 1)` for metrics naturally
  expressed as proportions or similarities (e.g. Jaccard, Bray-Curtis).
  All metrics should be expressed as similarities: higher value = more
  similar communities = thicker link.

## Value

A list with the metric definition.

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# Custom metric: robust Aitchison similarity (vegan >= 2.6).
# The method uses matrix completion and needs the full matrix, so we
# compute the full distance once in `prep` and look up each pair in `fn`.
rob_aitch <- make_sharing_metric(
  label = "Robust Aitchison similarity",
  color = "#FF7F00",
  fmt   = "%.2f",
  prep  = function(physeq, fact, modalities) {
    otu_sp <- .agg_by_mod(physeq, fact, modalities)
    as.matrix(vegan::vegdist(t(otu_sp), method = "robust.aitchison"))
  },
  fn    = \(a, b, otu_sp, cache) 1 - cache[a, b]
)

# Custom metric: number of unique (non-shared) species
unique_sp <- make_sharing_metric(
  label = "Unique species count",
  color = "#A65628",
  fmt   = "%.0f",
  fn    = \(a, b, otu_sp, cache) {
    sum(xor(otu_sp[, a] > 0, otu_sp[, b] > 0))
  }
)

# Custom metric with prep step: Sorensen similarity at Family rank
sor_family <- make_sharing_metric(
  label = "Sorensen at Family rank",
  color = "#F781BF",
  fmt   = "%.2f",
  prep  = function(physeq, fact, modalities) {
    d_fam <- phyloseq::tax_glom(physeq, "Family", NArm = FALSE)
    .agg_by_mod(d_fam, fact, modalities)
  },
  fn    = \(a, b, otu_sp, cache) {
    1 - as.numeric(vegan::vegdist(
      t(cache[, c(a, b)]),
      method = "bray", binary = TRUE
    ))
  }
)
} # }
```
