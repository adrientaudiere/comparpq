# Remove a phyloseq object from a list_phyloseq

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Usage

``` r
remove_phyloseq(x, name, verbose = TRUE)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- name:

  (character or integer) Name or index of the phyloseq object to remove.

## Value

A new list_phyloseq object without the removed phyloseq

## Examples

``` r
lpq <- list_phyloseq(list(run1 = data_fungi, run2 = data_fungi_mini),
  verbose = FALSE
)
lpq2 <- remove_phyloseq(lpq, "run2", verbose = FALSE)
length(lpq2)
#> [1] 1
```
