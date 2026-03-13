# Filter phyloseq objects to keep only shared samples and/or taxa

Filters each phyloseq object in a list_phyloseq to retain only the
samples and/or taxa that are common across all objects. This is useful
for making direct comparisons on a common basis.

## Usage

``` r
filter_common_lpq(
  x,
  filter_samples = TRUE,
  filter_taxa = FALSE,
  clean = TRUE,
  verbose = TRUE
)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- filter_samples:

  (logical, default TRUE) If TRUE, filter to keep only samples present
  in all phyloseq objects.

- filter_taxa:

  (logical, default FALSE) If TRUE, filter to keep only taxa present in
  all phyloseq objects.

- clean:

  (logical, default TRUE) If TRUE, apply
  [`MiscMetabar::clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html)
  after filtering to remove empty samples/taxa.

- verbose:

  (logical, default TRUE) If TRUE, print information about the filtering
  process.

## Value

A new list_phyloseq object with filtered phyloseq objects

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function is particularly useful for:

- **NESTED_ROBUSTNESS** comparisons: filter to common samples when
  comparing original vs rarefied data

- **EXPLORATION** comparisons: filter to common samples/taxa when
  comparing different sample groups

- Any comparison where you need to ensure all phyloseq objects contain
  the same samples or taxa

## Examples

``` r
lpq <- list_phyloseq(list(
  original = data_fungi,
  mini = data_fungi_mini,
  rarefied = rarefy_even_depth(data_fungi)
))
#> You set `rngseed` to FALSE. Make sure you've set & recorded
#>  the random seed of your session for reproducibility.
#> See `?set.seed`
#> ...
#> 1007OTUs were removed because they are no longer 
#> present in any sample after random subsampling
#> ...
#> ℹ Building summary table for 3 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 42 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

# Filter to keep only common samples (useful for nested comparisons)
lpq_filtered <- filter_common_lpq(lpq, filter_samples = TRUE, verbose = FALSE)
#> Cleaning suppress 10 taxa and 0 samples.
#> Cleaning suppress 61 taxa and 0 samples.

# Filter to keep only common taxa
lpq_filtered <- filter_common_lpq(lpq, filter_samples = FALSE, filter_taxa = TRUE)
#> Filtering to 42 common taxa.
#> Cleaning suppress 0 taxa and 48 samples.
#>   original: 185 -> 137 samples, 1420 -> 42 taxa
#>   mini: 137 -> 137 samples, 45 -> 42 taxa
#> Cleaning suppress 0 taxa and 115 samples.
#>   rarefied: 185 -> 70 samples, 413 -> 42 taxa
#> ℹ Building summary table for 3 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 70 common samples, 42 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

# Filter both samples and taxa
lpq_filtered <- filter_common_lpq(lpq, filter_samples = TRUE, filter_taxa = TRUE)
#> Filtering to 137 common samples.
#> Filtering to 42 common taxa.
#>   original: 185 -> 137 samples, 1420 -> 42 taxa
#>   mini: 137 -> 137 samples, 45 -> 42 taxa
#> Cleaning suppress 0 taxa and 67 samples.
#>   rarefied: 185 -> 70 samples, 413 -> 42 taxa
#> ℹ Building summary table for 3 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 70 common samples, 42 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)
```
