# Compare two taxonomic assignments of a phyloseq in one call

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Convenience wrapper that assembles, for a single pair of `tax_table`
columns (typically the same rank assigned by two databases), the
congruence metrics and the four main comparison views of comparpq. It
returns a named list bundling
[`tc_congruence_metrics()`](https://adrientaudiere.github.io/comparpq/reference/tc_congruence_metrics.md),
a filtered contingency table, and the
[`tc_bar()`](https://adrientaudiere.github.io/comparpq/reference/tc_bar.md),
[`tc_sankey()`](https://adrientaudiere.github.io/comparpq/reference/tc_sankey.md),
[`tc_heatmap()`](https://adrientaudiere.github.io/comparpq/reference/tc_heatmap.md)
plots plus a closure that draws the
[`tc_circle()`](https://adrientaudiere.github.io/comparpq/reference/tc_circle.md)
chord diagram on demand (the latter uses base graphics rather than
ggplot2).

## Usage

``` r
compare_taxo_db(physeq, rank_1, rank_2, rank_3 = NULL, min_n = 5)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq::phyloseq-class()`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object.

- rank_1:

  (character, required) `tax_table` column of the first database (e.g.
  `"Genus_SILVA"`).

- rank_2:

  (character, required) `tax_table` column of the second database (e.g.
  `"Genus_KSGP"`).

- rank_3:

  (character) Column used to color
  [`tc_bar()`](https://adrientaudiere.github.io/comparpq/reference/tc_bar.md).
  Defaults to the parent rank of `rank_1` in the same database when it
  can be resolved from the column naming (`"<Rank>_<Db>"`), otherwise
  `rank_1`.

- min_n:

  (integer, default 5) Minimum number of taxa for a contingency cell to
  be kept (cells with `n > min_n` are returned).

## Value

A named list with:

- congruence:

  result of
  [`tc_congruence_metrics()`](https://adrientaudiere.github.io/comparpq/reference/tc_congruence_metrics.md).

- contingency:

  data.frame of cells with more than `min_n` taxa, sorted by decreasing
  count.

- bar:

  a
  [`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  from
  [`tc_bar()`](https://adrientaudiere.github.io/comparpq/reference/tc_bar.md).

- sankey:

  a
  [`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  from
  [`tc_sankey()`](https://adrientaudiere.github.io/comparpq/reference/tc_sankey.md).

- circle:

  a function of no argument that draws
  [`tc_circle()`](https://adrientaudiere.github.io/comparpq/reference/tc_circle.md).

- heatmap:

  a
  [`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  from
  [`tc_heatmap()`](https://adrientaudiere.github.io/comparpq/reference/tc_heatmap.md).

## See also

[`count_taxo_congruence()`](https://adrientaudiere.github.io/comparpq/reference/count_taxo_congruence.md),
[`plot_congruence_counts()`](https://adrientaudiere.github.io/comparpq/reference/plot_congruence_counts.md),
[`tc_congruence_metrics()`](https://adrientaudiere.github.io/comparpq/reference/tc_congruence_metrics.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
pq <- MiscMetabar::subset_taxa_pq(Glom_otu, phyloseq::taxa_sums(Glom_otu) > 5000)
#> Cleaning suppress 0 taxa (  ) and 1 sample(s) ( samp_Blanc-PCR-racines ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1147
#> Number of filtered-out ASV 955
#> Number of kept ASV 192
#> Number of kept samples 443
res <- compare_taxo_db(pq, "Genus", "Genus__eukaryome_Glomero")
#> Warning: The `fun.y` argument of `stat_summary()` is deprecated as of ggplot2 3.3.0.
#> ℹ Please use the `fun` argument instead.
#> ℹ The deprecated feature was likely used in the comparpq package.
#>   Please report the issue at
#>   <https://github.com/adrientaudiere/comparpq/issues>.
res$contingency
#>              Genus Genus__eukaryome_Glomero  n
#> 39            <NA>                     <NA> 95
#> 38          Glomus                     <NA> 52
#> 12            <NA>            Entrophospora 10
#> 3             <NA>             Archaeospora  6
#> 10 Claroideoglomus            Entrophospora  6
res$bar

# }
```
