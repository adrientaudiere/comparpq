# Stacked barplot of taxonomic-assignment congruence across ranks

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Draws, for several taxonomic ranks at once, a stacked barplot of the
percentage of taxa (or sequences) falling into each congruence category
when comparing two databases (see
[`count_taxo_congruence()`](https://adrientaudiere.github.io/comparpq/reference/count_taxo_congruence.md)).
The `only_<rank>_<db>` categories are normalized to `only_<db>` so that
the color of a database is consistent across ranks (e.g.
`only_Order_EUK` and `only_Family_EUK` share the `only_EUK` color).

Columns can be supplied in two ways: explicitly through `ranks_1` /
`ranks_2` (paired column names), or automatically through `suffix_1` /
`suffix_2` (database suffixes such as `"_EUK"`), in which case columns
are paired by taxonomic rank.

## Usage

``` r
plot_congruence_counts(
  physeq,
  ranks_1 = NULL,
  ranks_2 = NULL,
  suffix_1 = NULL,
  suffix_2 = NULL,
  ranks = NULL,
  n_seq = FALSE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq::phyloseq-class()`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object.

- ranks_1:

  (character vector) `tax_table` column names for the first database,
  one per comparison. Ignored when `suffix_1` is provided.

- ranks_2:

  (character vector) `tax_table` column names for the second database,
  same length as `ranks_1`. Ignored when `suffix_2` is provided.

- suffix_1:

  (character) Column suffix (e.g. `"_EUK"`) selecting the columns of the
  first database automatically. When provided together with `suffix_2`,
  `ranks_1` / `ranks_2` are ignored and columns are paired by rank.

- suffix_2:

  (character) Column suffix (e.g. `"_Unite"`) for the second database.

- ranks:

  (character vector) Rank labels to plot, in the desired y-axis order
  (broadest to finest). Accepts standard ranks (`"Species"` -\>
  `Species_<db>`) and custom prefixes (`"genusSpeciesEpithet"`). A pair
  is kept only when both `<prefix><suffix_1>` and `<prefix><suffix_2>`
  exist; others are skipped with a warning. `NULL` (default)
  auto-selects all standard ranks present in both databases, in
  taxonomic order.

- n_seq:

  (logical, default `FALSE`) `FALSE`: bars show the percentage of taxa
  (`pct_asv`); `TRUE`: bars show the percentage of sequences
  (`pct_seq`).

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object (stacked bars, one bar per rank, filled by normalized congruence
category).

## See also

[`count_taxo_congruence()`](https://adrientaudiere.github.io/comparpq/reference/count_taxo_congruence.md),
[`compare_taxo_db()`](https://adrientaudiere.github.io/comparpq/reference/compare_taxo_db.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
plot_congruence_counts(
  Glom_otu,
  ranks_1 = c("Order", "Family", "Genus"),
  ranks_2 = c(
    "Order__eukaryome_Glomero",
    "Family__eukaryome_Glomero",
    "Genus__eukaryome_Glomero"
  )
)

# }
```
