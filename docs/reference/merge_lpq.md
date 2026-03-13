# Merge a list_phyloseq into a single phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Merges all phyloseq objects from a
[list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)
into a single phyloseq object where each original phyloseq becomes one
sample. Abundances are summed across samples within each phyloseq
object.

Taxa are matched across phyloseq objects either by reference sequences
(`refseq` slot, default) or by taxa names. Matching by refseq is
preferred because taxa names are often inconsistent across independently
built phyloseq objects (e.g., ASV_1 in one object is not the same taxon
as ASV_1 in another).

## Usage

``` r
merge_lpq(
  x,
  match_by = c("refseq", "names"),
  tax_priority = 1L,
  verbose = TRUE
)
```

## Arguments

- x:

  (list_phyloseq, required) A
  [list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)
  object.

- match_by:

  (character, default `"refseq"`) How to match taxa across phyloseq
  objects. One of:

  - `"refseq"`: match by DNA sequence in the `refseq` slot
    (recommended). All phyloseq objects must have a `refseq` slot.

  - `"names"`: match by taxa names. Use only when taxa names are
    consistent across objects (e.g., same pipeline, same database).

- tax_priority:

  (character or integer, default 1L) Which phyloseq object's taxonomy to
  use when taxa are matched. Either a name from the list or an integer
  index. When a taxon appears in multiple objects, the taxonomy from the
  priority object is used; if absent there, the first available taxonomy
  is used.

- verbose:

  (logical, default TRUE) Print information about the merge.

## Value

A phyloseq object with:

- `otu_table`:

  One column per original phyloseq object (summed across its samples),
  one row per unique taxon.

- `sample_data`:

  One row per original phyloseq, with a column `source_name` containing
  the list_phyloseq names.

- `tax_table`:

  Taxonomy from the priority object (or first available).

- `refseq`:

  Present when `match_by = "refseq"`.

## See also

[list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md),
[`simple_venn_pq()`](https://adrientaudiere.github.io/comparpq/reference/simple_venn_pq.md)

## Author

Adrien Taudière

## Examples

``` r
lpq <- list_phyloseq(list(
  fungi = data_fungi_mini,
  fungi2 = data_fungi_mini
))
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: REPRODUCIBILITY
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (REPRODUCIBILITY)

# Merge by refseq (default)
merged <- merge_lpq(lpq)
#> Merging 2 phyloseq objects by refseq: 45 + 45 taxa -> 45 unique sequences.
merged
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 2 samples ]
#> sample_data() Sample Data:       [ 2 samples by 2 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 45 reference sequences ]

# Merge by taxa names
merged_names <- merge_lpq(lpq, match_by = "names")
#> Merging 2 phyloseq objects by names: 45 + 45 taxa -> 45 unique names.
```
