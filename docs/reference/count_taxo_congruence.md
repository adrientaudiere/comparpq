# Count taxa and sequences by congruence of two taxonomic assignments

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Classifies each taxon of a phyloseq object into a congruence category by
comparing two `tax_table` columns (typically the same taxonomic rank
assigned by two different databases, algorithms or reference sets), and
counts both the number of taxa (ASV/OTU) and the number of sequences
(reads) falling in each category. The five categories are: `both_equal`
(identical non-missing assignment), `both_na` (missing in both),
`only_<rank_1>` (assigned only by the first column), `only_<rank_2>`
(assigned only by the second column), and `different` (both assigned but
disagreeing). Empty strings and the literal values `"NA"` / `"NA_NA"`
are treated as missing.

## Usage

``` r
count_taxo_congruence(physeq, rank_1, rank_2)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq::phyloseq-class()`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object.

- rank_1:

  (character, required) Name of the `tax_table` column for the first
  assignment (e.g. `"Genus"` or `"Genus_SILVA"`).

- rank_2:

  (character, required) Name of the `tax_table` column for the second
  assignment (e.g. `"Genus__eukaryome_Glomero"`).

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with one row per category and the columns `category`, `n_asv` (number of
taxa), `n_seq` (summed sequence counts), `pct_asv` and `pct_seq`
(percentages, rounded to two decimals). The `only_<rank_1>` and
`only_<rank_2>` category names embed the actual column names passed as
arguments.

## See also

[`plot_congruence_counts()`](https://adrientaudiere.github.io/comparpq/reference/plot_congruence_counts.md),
[`compare_taxo_db()`](https://adrientaudiere.github.io/comparpq/reference/compare_taxo_db.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
count_taxo_congruence(Glom_otu, "Genus", "Genus__eukaryome_Glomero")
#> # A tibble: 5 × 5
#>   category                      n_asv   n_seq pct_asv pct_seq
#>   <chr>                         <dbl>   <dbl>   <dbl>   <dbl>
#> 1 both_equal                        4    4143    0.35    0.02
#> 2 both_na                         664 6401309   57.9    31.7 
#> 3 only_Genus                      318 5922501   27.7    29.4 
#> 4 only_Genus__eukaryome_Glomero    73 2676722    6.36   13.3 
#> 5 different                        88 5163760    7.67   25.6 
# }
```
