# Apply a function to all phyloseq objects in a list_phyloseq

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Applies a function to each phyloseq object in a list_phyloseq and
returns a new list_phyloseq with the transformed objects. This is useful
for applying cleaning, filtering, or transformation functions uniformly
across all phyloseq objects.

## Usage

``` r
apply_to_lpq(x, .f, ..., compute_dist = NULL, verbose = TRUE)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- .f:

  (function, required) A function to apply to each phyloseq object. The
  function must take a phyloseq object as its first argument and return
  a phyloseq object.

- ...:

  Additional arguments passed to `.f`.

- compute_dist:

  (logical, default NULL) Whether to compute pairwise k-mer distances in
  the rebuilt list_phyloseq. When NULL (default), automatically inherits
  the setting from `x`: distances are recomputed only if they were
  computed in the original object (i.e. `x@comparison$refseq_comparison`
  is not NULL).

- verbose:

  (logical, default TRUE) If TRUE, print information about the
  transformation applied to each phyloseq object.

## Value

A new list_phyloseq object with the transformed phyloseq objects. The
comparison parameters (`same_primer_seq_tech`, `same_bioinfo_pipeline`)
are preserved from the original object.

## Details

Common functions to apply include:

- [`MiscMetabar::clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html):
  Remove empty samples and taxa

- `phyloseq::taxa_as_rows()`: Ensure taxa are rows in otu_table

- [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html):
  Rarefy to even depth

- [`phyloseq::transform_sample_counts()`](https://rdrr.io/pkg/phyloseq/man/transformcounts.html):
  Transform counts (e.g., relative abundance)

- [`phyloseq::subset_taxa()`](https://rdrr.io/pkg/phyloseq/man/subset_taxa-methods.html):
  Filter taxa based on criteria

- [`phyloseq::subset_samples()`](https://rdrr.io/pkg/phyloseq/man/subset_samples-methods.html):
  Filter samples based on criteria

## See also

[`filter_common_lpq()`](https://adrientaudiere.github.io/comparpq/reference/filter_common_lpq.md),
[`update_list_phyloseq()`](https://adrientaudiere.github.io/comparpq/reference/update_list_phyloseq.md)

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

# Apply clean_pq to all phyloseq objects
lpq_clean <- apply_to_lpq(lpq, MiscMetabar::clean_pq)
#>   fungi: 185 -> 185 samples, 1420 -> 1420 taxa
#>   fungi_mini: 137 -> 137 samples, 45 -> 45 taxa
#> Applied `MiscMetabar::clean_pq` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

# Apply taxa_as_rows
lpq_rows <- apply_to_lpq(lpq, MiscMetabar::taxa_as_rows)
#> Taxa are now in rows.
#>   fungi: 185 -> 185 samples, 1420 -> 1420 taxa
#>   fungi_mini: 137 -> 137 samples, 45 -> 45 taxa
#> Applied `MiscMetabar::taxa_as_rows` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

# Apply rarefy_even_depth with a specific rngseed
lpq_rar <- apply_to_lpq(lpq, rarefy_even_depth, rngseed = 21)
#> `set.seed(21)` was used to initialize repeatable random subsampling.
#> Please record this for your records so others can reproduce.
#> Try `set.seed(21); .Random.seed` for the full vector
#> ...
#> 1000OTUs were removed because they are no longer 
#> present in any sample after random subsampling
#> ...
#>   fungi: 185 -> 185 samples, 1420 -> 420 taxa
#> `set.seed(21)` was used to initialize repeatable random subsampling.
#> Please record this for your records so others can reproduce.
#> Try `set.seed(21); .Random.seed` for the full vector
#> ...
#> 7OTUs were removed because they are no longer 
#> present in any sample after random subsampling
#> ...
#>   fungi_mini: 137 -> 137 samples, 45 -> 38 taxa
#> Applied `rarefy_even_depth` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 35 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

lpq_rar
#> list_phyloseq object with 2 phyloseq objects
#> 
#> --- Summary ---
#> # A tibble: 2 × 17
#>   name   n_samples n_taxa n_sequences n_occurence mean_seq_length min_seq_length
#>   <chr>      <int>  <int>       <dbl>       <int>           <dbl>          <int>
#> 1 fungi        185    420        1110         695            323.            255
#> 2 fungi…       137     38         137         137            344.            303
#> # ℹ 10 more variables: mean_seq_per_sample <dbl>, sd_seq_per_sample <dbl>,
#> #   min_seq_per_sample <dbl>, max_seq_per_sample <dbl>,
#> #   mean_seq_per_taxon <dbl>, sd_seq_per_taxon <dbl>, has_sam_data <lgl>,
#> #   has_tax_table <lgl>, has_refseq <lgl>, has_phy_tree <lgl>
#> 
#> --- Comparison characteristics ---
#> Type of comparison: NESTED_ROBUSTNESS 
#> Nested samples detected (one phyloseq derived from another).
#>   Nesting: fungi_mini (137 samples) is nested in fungi (185 samples)
#>   Useful to test robustness to data processing (e.g., rarefaction).
#>   Comparisons should focus on the common (nested) samples. 
#> 
#> Same primer/seq tech: TRUE 
#> Same bioinfo pipeline: TRUE 
#> Same sample_data structure: TRUE 
#> Same samples: FALSE 
#> Nested samples: TRUE 
#> Same taxa: FALSE 
#> Common samples: 137 
#> Common taxa: 35 
#> 
#> --- Reference sequence comparison ---
#> fungi_vs_fungi_mini: 35 shared seqs, 385 unique in fungi, 3 unique in fungi_mini

# Transform to relative abundance
lpq_rel <- apply_to_lpq(
  lpq,
  phyloseq::transform_sample_counts,
  function(x) x / sum(x)
)
#>   fungi: 185 -> 185 samples, 1420 -> 1420 taxa
#>   fungi_mini: 137 -> 137 samples, 45 -> 45 taxa
#> Applied `phyloseq::transform_sample_counts` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

# Chain multiple transformations
lpq_processed <- lpq |>
  apply_to_lpq(MiscMetabar::clean_pq) |>
  apply_to_lpq(MiscMetabar::taxa_as_rows)
#>   fungi: 185 -> 185 samples, 1420 -> 1420 taxa
#>   fungi_mini: 137 -> 137 samples, 45 -> 45 taxa
#> Applied `MiscMetabar::clean_pq` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)
#> Taxa are now in rows.
#>   fungi: 185 -> 185 samples, 1420 -> 1420 taxa
#>   fungi_mini: 137 -> 137 samples, 45 -> 45 taxa
#> Applied `MiscMetabar::taxa_as_rows` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)
```
