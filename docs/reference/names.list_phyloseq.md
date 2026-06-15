# Names of phyloseq objects in a `list_phyloseq`

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Usage

``` r
# S3 method for class 'list_phyloseq'
names(x)
```

## Arguments

- x:

  A `list_phyloseq` object.

## Value

A character vector of names.

## Examples

``` r
lpq <- list_phyloseq(list(a = data_fungi, b = data_fungi))
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: REPRODUCIBILITY
#> ℹ 185 common samples, 1420 common taxa
#> ✔ list_phyloseq created (REPRODUCIBILITY)
names(lpq) # c("a", "b")
#> [1] "a" "b"
```
