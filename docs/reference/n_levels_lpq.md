# Count unique taxonomic levels across phyloseq objects

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Creates a summary table showing the number of unique taxonomic values
(levels) for each taxonomic rank across all phyloseq objects in a
list_phyloseq. This is useful for comparing taxonomic resolution and
diversity across different datasets or classification methods.

## Usage

``` r
n_levels_lpq(x, taxonomic_ranks, na.rm = TRUE)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- taxonomic_ranks:

  (character vector, required) Names of taxonomic ranks to count. Must
  be present in the tax_table of ALL phyloseq objects in the list.

- na.rm:

  (logical, default TRUE) If TRUE, NA values are excluded when counting
  unique levels.

## Value

A data frame with:

- Rows:

  One row per phyloseq object (named by the phyloseq name)

- Columns:

  One column per taxonomic rank, containing the count of unique values
  for that rank in that phyloseq object

## See also

[`upset_lpq()`](https://adrientaudiere.github.io/comparpq/reference/upset_lpq.md)

## Author

Adrien Taudière

## Examples

``` r
lpq <- list_phyloseq(list(fungi = data_fungi, fungi_mini = data_fungi_mini))
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

n_levels_lpq(lpq, c("Phylum", "Class", "Order", "Family", "Genus"))
#>            Phylum Class Order Family Genus
#> fungi           6    25    73    163   252
#> fungi_mini      1     3     9     20    22
```
