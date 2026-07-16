# Add a phyloseq object to a list_phyloseq

Append a phyloseq object to an existing `list_phyloseq`, returning a new
object that includes the added element (and its recomputed comparison
characteristics).
[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Usage

``` r
add_phyloseq(x, physeq, name = NULL, verbose = TRUE)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- physeq:

  (required) A phyloseq object to add.

- name:

  (character, default NULL) Optional name for the new phyloseq object.
  If NULL, a name is generated automatically.

- verbose:

  (logical, default TRUE) If TRUE, print progress messages.

## Value

A new list_phyloseq object with the added phyloseq

## Examples

``` r
lpq <- list_phyloseq(list(run1 = data_fungi), verbose = FALSE)
lpq2 <- add_phyloseq(lpq, data_fungi_mini, name = "run2", verbose = FALSE)
length(lpq2)
#> [1] 2
```
