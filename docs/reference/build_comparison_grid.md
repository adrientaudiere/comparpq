# Build the grid of pairwise database comparisons across phyloseq objects

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Scans the `tax_table` of each phyloseq object for columns named
`"<Rank>_<Db>"` (e.g. `"Genus_SILVA"`, `"Genus_KSGP"`) and enumerates,
for each object and each rank, every pairwise combination of databases.
The resulting long-format grid has one row per comparison and can be
iterated over to run
[`compare_taxo_db()`](https://adrientaudiere.github.io/comparpq/reference/compare_taxo_db.md)
or
[`tc_congruence_metrics()`](https://adrientaudiere.github.io/comparpq/reference/tc_congruence_metrics.md)
systematically across objects, ranks and database pairs.

## Usage

``` r
build_comparison_grid(pq_list, ranks = c("Order", "Family", "Genus"))
```

## Arguments

- pq_list:

  (required) A named list of
  [`phyloseq::phyloseq-class()`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  objects, or a
  [`list_phyloseq()`](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)
  object. Names are used in the `physeq` column of the grid; unnamed
  lists are auto-named `pq1`, `pq2`, ...

- ranks:

  (character, default `c("Order", "Family", "Genus")`) Taxonomic ranks
  to consider. Ranks with fewer than two database columns in a given
  phyloseq object are skipped for that object.

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with one row per pairwise comparison and the columns `physeq`, `rank`,
`db_1`, `col_1`, `db_2` and `col_2`, where `col_1`/`col_2` are the
`tax_table` column names to pass to
[`compare_taxo_db()`](https://adrientaudiere.github.io/comparpq/reference/compare_taxo_db.md).
Returns an empty tibble with the same columns when no rank carries two
databases.

## See also

[`compare_taxo_db()`](https://adrientaudiere.github.io/comparpq/reference/compare_taxo_db.md),
[`count_taxo_congruence()`](https://adrientaudiere.github.io/comparpq/reference/count_taxo_congruence.md),
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
mat <- as(pq@tax_table, "matrix")
pq@tax_table <- phyloseq::tax_table(cbind(
  mat,
  Genus_SILVA = mat[, "Genus"],
  Genus_KSGP = mat[, "Genus"]
))
build_comparison_grid(list(glom = pq))
#> # A tibble: 3 × 6
#>   physeq rank  db_1  col_1       db_2               col_2                   
#>   <chr>  <chr> <chr> <chr>       <chr>              <chr>                   
#> 1 glom   Genus KSGP  Genus_KSGP  SILVA              Genus_SILVA             
#> 2 glom   Genus KSGP  Genus_KSGP  _eukaryome_Glomero Genus__eukaryome_Glomero
#> 3 glom   Genus SILVA Genus_SILVA _eukaryome_Glomero Genus__eukaryome_Glomero
# }
```
