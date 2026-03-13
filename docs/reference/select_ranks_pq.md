# Select taxonomic ranks in a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful in tidyverse-like pipeline.

## Usage

``` r
select_ranks_pq(physeq, ...)
```

## Arguments

- physeq:

  (required) A phyloseq object.

- ...:

  One or more unquoted expressions separated by commas. Variable names
  can be used as if they were positions in the data frame, so
  expressions like x:y can be used to select a range of variables. See
  ?dplyr::select.

## Value

A phyloseq object

## Author

Adrien Taudière

## Examples

``` r
select_ranks_pq(data_fungi, Order, Family)@tax_table |>
  dim()
#> [1] 1420    2
select_ranks_pq(data_fungi, !Order)@tax_table |>
  colnames()
#>  [1] "Domain"             "Phylum"             "Class"             
#>  [4] "Family"             "Genus"              "Species"           
#>  [7] "Trophic.Mode"       "Guild"              "Trait"             
#> [10] "Confidence.Ranking" "Genus_species"     
#
```
