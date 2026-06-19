# Compare refseq slots between two phyloseq objects

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Performs a detailed comparison of the reference sequences (`refseq`
slot) between two phyloseq objects. Identifies shared and unique
ASVs/OTUs both by taxa name and by actual DNA sequence content.

This is useful to detect subtle differences between two phyloseq objects
that may share the same samples but differ in their ASV/OTU composition,
e.g. after different bioinformatics pipelines or filtering steps.

When unique (non-shared) sequences exist, the mean nearest-neighbor
k-mer distance from each unique sequence to the other set is computed.
This gives a sense of how different the unmatched sequences are from
their closest counterpart in the other object.

## Usage

``` r
compare_refseq(
  physeq1,
  physeq2 = NULL,
  name1 = NULL,
  name2 = NULL,
  k = 5,
  max_seqs = 500,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- physeq1:

  (phyloseq or list_phyloseq, required) First phyloseq object, which
  must have a `refseq` slot. Alternatively, a
  [list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)
  object; in that case the first two phyloseq objects are used (with
  their names) and `physeq2` is ignored.

- physeq2:

  (phyloseq, default NULL) Second phyloseq object. Must have a `refseq`
  slot. Ignored when `physeq1` is a
  [list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md).

- name1:

  (character, default NULL) Label for the first phyloseq object in the
  output. If NULL, inferred from the
  [list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)
  names or defaults to `"physeq1"`.

- name2:

  (character, default NULL) Label for the second phyloseq object in the
  output. If NULL, inferred from the
  [list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)
  names or defaults to `"physeq2"`.

- k:

  (integer, default 5) k-mer size for nearest-neighbor distance
  computation on unique sequences.

- max_seqs:

  (integer, default 500) Maximum number of query and target sequences to
  use for k-mer distance computation. If there are more, a random sample
  is drawn. Set to `Inf` to use all (may be slow).

- seed:

  (integer, default NULL) Random seed for reproducible sampling when
  `max_seqs` is exceeded.

- verbose:

  (logical, default TRUE) If TRUE, print a summary of the comparison.

## Value

A list of class `"compare_refseq"` with components:

- `n_taxa`:

  Named integer vector with the number of taxa in each object.

- `shared_names`:

  Character vector of taxa names present in both objects.

- `unique_names_1`:

  Taxa names only in `physeq1`.

- `unique_names_2`:

  Taxa names only in `physeq2`.

- `shared_seqs`:

  Character vector of unique DNA sequences found in both objects
  (regardless of taxa name).

- `unique_seqs_1`:

  DNA sequences only in `physeq1`.

- `unique_seqs_2`:

  DNA sequences only in `physeq2`.

- `same_name_diff_seq`:

  A data frame of taxa with the same name but different sequences across
  objects. Columns: `taxa_name`, `seq_physeq1`, `seq_physeq2`.

- `same_seq_diff_name`:

  A data frame of sequences shared between objects under different
  names. Columns: `sequence`, `name_physeq1`, `name_physeq2`.

- `mean_nn_kmer_dist_1`:

  Mean nearest-neighbor k-mer distance from sequences unique to
  `physeq1` to the closest sequence in `physeq2`. `NA` if no unique
  sequences.

- `mean_nn_kmer_dist_2`:

  Mean nearest-neighbor k-mer distance from sequences unique to
  `physeq2` to the closest sequence in `physeq1`. `NA` if no unique
  sequences.

## Author

Adrien Taudiere

## Examples

``` r
# Compare a phyloseq object with itself (all shared)
res <- compare_refseq(data_fungi_mini, data_fungi_mini)
#> == Reference Sequence Comparison ==
#> physeq1 : 45 taxa
#> physeq2 : 45 taxa
#> 
#> -- By taxa name --
#>   Shared : 45 
#>   Only in physeq1 : 0 
#>   Only in physeq2 : 0 
#> 
#> -- By DNA sequence --
#>   Shared : 45 
#>   Only in physeq1 : 0 
#>   Only in physeq2 : 0 
#> 
res
#> == Reference Sequence Comparison ==
#> physeq1 : 45 taxa
#> physeq2 : 45 taxa
#> 
#> -- By taxa name --
#>   Shared : 45 
#>   Only in physeq1 : 0 
#>   Only in physeq2 : 0 
#> 
#> -- By DNA sequence --
#>   Shared : 45 
#>   Only in physeq1 : 0 
#>   Only in physeq2 : 0 
#> 

# Compare with a subset
sub <- prune_taxa(taxa_names(data_fungi_mini)[1:20], data_fungi_mini)
res2 <- compare_refseq(data_fungi_mini, sub,
  name1 = "full", name2 = "subset"
)
#> == Reference Sequence Comparison ==
#> full : 45 taxa
#> subset : 20 taxa
#> 
#> -- By taxa name --
#>   Shared : 20 
#>   Only in full : 25 
#>   Only in subset : 0 
#> 
#> -- By DNA sequence --
#>   Shared : 20 
#>   Only in full : 25 
#>   Only in subset : 0 
#> 
#> 
#> -- Mean nearest-neighbor k-mer distance (k=5) --
#>   Unique in full -> subset : 0.121 
res2
#> == Reference Sequence Comparison ==
#> full : 45 taxa
#> subset : 20 taxa
#> 
#> -- By taxa name --
#>   Shared : 20 
#>   Only in full : 25 
#>   Only in subset : 0 
#> 
#> -- By DNA sequence --
#>   Shared : 20 
#>   Only in full : 25 
#>   Only in subset : 0 
#> 
#> 
#> -- Mean nearest-neighbor k-mer distance (k=5) --
#>   Unique in full -> subset : 0.121 

# From a list_phyloseq (uses the first two objects and their names)
lpq <- list_phyloseq(list(full = data_fungi_mini, subset = sub))
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: REPRODUCIBILITY
#> ℹ 137 common samples, 20 common taxa
#> ✔ list_phyloseq created (REPRODUCIBILITY)
res3 <- compare_refseq(lpq)
#> == Reference Sequence Comparison ==
#> full : 45 taxa
#> subset : 20 taxa
#> 
#> -- By taxa name --
#>   Shared : 20 
#>   Only in full : 25 
#>   Only in subset : 0 
#> 
#> -- By DNA sequence --
#>   Shared : 20 
#>   Only in full : 25 
#>   Only in subset : 0 
#> 
#> 
#> -- Mean nearest-neighbor k-mer distance (k=5) --
#>   Unique in full -> subset : 0.121 
res3
#> == Reference Sequence Comparison ==
#> full : 45 taxa
#> subset : 20 taxa
#> 
#> -- By taxa name --
#>   Shared : 20 
#>   Only in full : 25 
#>   Only in subset : 0 
#> 
#> -- By DNA sequence --
#>   Shared : 20 
#>   Only in full : 25 
#>   Only in subset : 0 
#> 
#> 
#> -- Mean nearest-neighbor k-mer distance (k=5) --
#>   Unique in full -> subset : 0.121 
```
