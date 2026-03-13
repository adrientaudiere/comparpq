# Replace taxonomic value with a given pattern by NA

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use [`base::gsub()`](https://rdrr.io/r/base/grep.html) for substitution.

## Usage

``` r
taxtab_replace_pattern_by_NA(
  physeq,
  patterns = c(".*_incertae_sedis", "unclassified.*"),
  taxonomic_ranks = NULL,
  progress_bar = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- patterns:

  (character vector, default c(".*\_incertae_sedis", "unclassified.*"))
  A vector of patterns to select taxonomic values to convert to NA.

- taxonomic_ranks:

  (character vector, default NULL) A vector of taxonomic ranks where the
  substitution should occur. If left to NULL, all taxonomic ranks are
  modified.

- progress_bar:

  (logical, default FALSE) If TRUE, print progress during the
  calculation.

- ...:

  Additional arguments passed to
  [`base::gsub()`](https://rdrr.io/r/base/grep.html).

## Value

A phyloseq object

## Author

Adrien Taudière

## Examples

``` r
data_fungi@tax_table["ASV85", "Family"]
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>       Family                              
#> ASV85 "Hymenochaetales_fam_Incertae_sedis"
data_fungi2 <- taxtab_replace_pattern_by_NA(data_fungi, "fam_Incertae_sedis", taxonomic_ranks = "Family")
data_fungi2@tax_table["ASV85", "Family"]
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>       Family
#> ASV85 NA    

# By default patterns ".*_incertae_sedis" and "unclassified.*" are replaced by NA
data_fungi3 <-
  taxtab_replace_pattern_by_NA(data_fungi, ignore.case = TRUE, progress_bar = TRUE)
#>   |                                                          |                                                  |   0%  |                                                          |==                                                |   4%  |                                                          |====                                              |   8%  |                                                          |======                                            |  12%  |                                                          |========                                          |  17%  |                                                          |==========                                        |  21%  |                                                          |============                                      |  25%  |                                                          |===============                                   |  29%  |                                                          |=================                                 |  33%  |                                                          |===================                               |  38%  |                                                          |=====================                             |  42%  |                                                          |=======================                           |  46%  |                                                          |=========================                         |  50%  |                                                          |===========================                       |  54%  |                                                          |=============================                     |  58%  |                                                          |===============================                   |  62%  |                                                          |=================================                 |  67%  |                                                          |===================================               |  71%  |                                                          |======================================            |  75%  |                                                          |========================================          |  79%  |                                                          |==========================================        |  83%  |                                                          |============================================      |  88%  |                                                          |==============================================    |  92%  |                                                          |================================================  |  96%  |                                                          |==================================================| 100%
data_fungi3@tax_table["ASV85", "Family"]
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>       Family
#> ASV85 NA    
```
