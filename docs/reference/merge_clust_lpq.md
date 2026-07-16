# Merge a list_phyloseq into one phyloseq by clustering reference sequences

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Merges all phyloseq objects from a
[list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)
into a **single** phyloseq object while **keeping every sample
separate**. Taxa are unified *across* objects by clustering their
reference sequences (`refseq` slot) with a clustering algorithm (vsearch
by default), so that sequences that are similar enough (controlled by
`id`) become a single merged taxon whose counts are summed within each
original sample.

This differs from
[`merge_lpq()`](https://adrientaudiere.github.io/comparpq/reference/merge_lpq.md),
which collapses each phyloseq object into a single sample and matches
taxa by **exact** sequence identity. Use `merge_clust_lpq()` when you
want to pool the raw samples of several independently built datasets
into one object and reconcile their ASVs/OTUs through fuzzy
(identity-threshold) clustering.

Samples keep their original names; when the same sample name occurs in
more than one object, the parent object name is appended as a suffix
(`<sample>_<object>`) to keep sample names unique. A column
(`source_name` by default) recording the parent object is added to the
`sample_data`.

## Usage

``` r
merge_clust_lpq(
  x,
  method = "vsearch",
  id = 0.97,
  source_col = "source_name",
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  (list_phyloseq or list, required) A
  [list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)
  object, or a (preferably named) list of phyloseq objects, all of which
  must have a `refseq` slot.

- method:

  (character, default `"vsearch"`) Clustering method passed to
  [`MiscMetabar::postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.html).
  One of `"clusterize"`, `"vsearch"`, `"swarm"` or `"mmseqs2"`.

- id:

  (numeric, default 0.97) Sequence identity threshold for the clustering
  (used by the `vsearch`/`mmseqs2` methods). Ignored by `swarm`.

- source_col:

  (character, default `"source_name"`) Name of the `sample_data` column
  that will store the parent object name for each sample.

- verbose:

  (logical, default TRUE) Print information about the merge.

- ...:

  Further arguments passed to
  [`MiscMetabar::postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.html)
  (e.g. `tax_adjust`, `rank_propagation`, `nproc`).

## Value

A phyloseq object with:

- `otu_table`:

  One column per original sample (across all objects, suffixed on name
  collision), one row per clustered taxon.

- `sample_data`:

  One row per original sample, with a `source_col` column giving the
  parent object name.

- `tax_table`:

  Taxonomy carried from the clustering representative of each merged
  taxon.

- `refseq`:

  Representative sequence of each clustered taxon.

## See also

[`merge_lpq()`](https://adrientaudiere.github.io/comparpq/reference/merge_lpq.md),
[list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md),
[`MiscMetabar::postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.html)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires vsearch to be installed by default (MiscMetabar::install_vsearch())
library(MiscMetabar)
pq1 <- postcluster_pq(data_fungi_mini, method = "vsearch", id = 0.97)
pq2 <- clean_pq(prune_samples(sample_names(pq1)[1:4], pq1))

lpq <- list_phyloseq(list(run_a = pq2, run_b = pq1))

merged <- merge_clust_lpq(lpq, id = 0.97)
merged
table(sample_data(merged)$source_name)
} # }
```
