# Circle of correspondence between two taxonomic levels

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful to compare taxonomy from two different source/db/algo
side-by-side

## Usage

``` r
tc_circle(
  physeq,
  rank_1 = NULL,
  rank_2 = NULL,
  suffix_1 = "_1",
  suffix_2 = "_2"
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- rank_1:

  (character or integer) Define the first taxonomic rank as the number
  or the name of the column in tax_table slot.

- rank_2:

  (character or integer) Define the second taxonomic rank as the number
  or the name of the column in tax_table slot.

- suffix_1:

  (character, default "\_1") A suffix to add to rank_1 values.

- suffix_2:

  (character, default "\_2") A suffix to add to rank_2 values.

## Value

A ggplot2 object

## Author

Adrien Taudière

## Examples

``` r
tc_circle(
  Glom_otu,
  "Genus__eukaryome_Glomero",
  "Genus",
  suffix_1 = "_Euk",
  suffix_2 = "_Marjaam"
)


tc_circle(
  Glom_otu,
  "Family__eukaryome_Glomero",
  "Family",
  suffix_1 = "_Euk",
  suffix_2 = "_Marjaam"
)
```
