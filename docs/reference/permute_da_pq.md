# Simulate differential abundance by redistributing OTU counts

Creates a fake differential abundance signal by redistributing counts
within matched samples: selected "DA" taxa receive counts taken from
non-selected taxa. Library sizes (row sums) are preserved, but taxa
totals (column sums) are allowed to change to create a detectable DA
signal.

This approach is suitable for testing DA methods (ALDEx2, ANCOM-BC)
because it creates a true compositional shift while maintaining the
library size normalization these methods rely on.

## Usage

``` r
permute_da_pq(
  physeq,
  fact,
  conditions,
  effect_size = 3,
  prop_taxa = 0.05,
  min_prevalence = 0.1,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- fact:

  (character, required) Column name in `sample_data` for grouping.

- conditions:

  (character vector, required) Levels of `fact` where DA should be
  simulated.

- effect_size:

  (numeric, default 3) The fold-change in relative abundance for
  selected taxa in matched samples. Values \> 1 increase abundance.

- prop_taxa:

  (numeric, default 0.05) Proportion of taxa to be selected as DA taxa.

- min_prevalence:

  (numeric, default 0.1) Minimum prevalence in matched samples for a
  taxon to be eligible.

- seed:

  (integer, default NULL) Random seed for reproducibility.

- verbose:

  (logical, default FALSE) Print progress info.

## Value

A phyloseq object with modified OTU counts creating a DA signal. The
selected DA taxa are stored in `attr(result, "da_taxa")`.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The algorithm works as follows:

1.  Select a subset of prevalent taxa as "DA" taxa

2.  For matched samples only:

    - Multiply selected taxa counts by `effect_size`

    - Scale down non-selected taxa to preserve library size

3.  Round to integers

Unlike
[`multiply_counts_pq()`](https://adrientaudiere.github.io/comparpq/reference/multiply_counts_pq.md),
this function guarantees that library sizes are exactly preserved
(except for rounding), creating a pure compositional shift that DA
methods can detect.

## See also

[`multiply_counts_pq()`](https://adrientaudiere.github.io/comparpq/reference/multiply_counts_pq.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
d_da <- permute_da_pq(
  data_fungi,
  fact = "Height",
  conditions = "High",
  effect_size = 3,
  prop_taxa = 0.05,
  seed = 123
)

# Library sizes are preserved
all.equal(sample_sums(data_fungi), sample_sums(d_da))

# Test with ALDEx2
res <- MiscMetabar::aldex_pq(d_da,
  bifactor = "Height",
  modalities = c("Low", "High")
)
gg_aldex_plot(res, type = "volcano")

# Check which taxa were made DA
attr(d_da, "da_taxa")
} # }
```
