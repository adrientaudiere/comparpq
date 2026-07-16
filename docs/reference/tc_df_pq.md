# Contingency table of two taxonomic ranks

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Creates a cross-tabulation (contingency table) comparing two taxonomic
ranks from a phyloseq object. Useful for comparing taxonomic assignments
from different databases, algorithms, or taxonomic levels.

## Usage

``` r
tc_df_pq(physeq, rank_1 = "Family", rank_2 = "Class", ...)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- rank_1:

  (character, default "Family") The name of the first taxonomic rank
  (column in tax_table slot) for the cross-tabulation rows.

- rank_2:

  (character, default "Class") The name of the second taxonomic rank
  (column in tax_table slot) for the cross-tabulation columns.

- ...:

  Additional arguments passed to
  [`gtsummary::tbl_cross()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_cross.html).

## Value

A gtsummary tbl_cross object displaying the cross-tabulation of the two
taxonomic ranks.

## Author

Adrien Taudière

## Examples

``` r
tc_df_pq(data_fungi_mini)


  

```
