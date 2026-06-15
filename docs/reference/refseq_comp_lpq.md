# Compare refseq sequences across all objects in a list_phyloseq

For each pair of phyloseq objects in `x`, computes pairwise sequence
similarity between `@refseq` slots using k-mer Jaccard similarity, then
groups near-identical sequences across all objects via connected
components (union-find). Returns a list of results per threshold, and
one Venn diagram per threshold (requires ggVennDiagram).

Similarity is approximated as `p ≈ (2j / (1 + j))^(1/k)` where `j` is
the Jaccard index of shared k-mer presence and `k` is the k-mer size.
The approximation is accurate within ~2 percentage points for identities
in the 90–100% range.

## Usage

``` r
refseq_comp_lpq(x, thresholds = c(0, 1, 2, 3, 5, 10), k = 7, verbose = TRUE)
```

## Arguments

- x:

  (required) A
  [list_phyloseq](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)
  object. All phyloseq objects must have a populated `@refseq` slot.
  Supports 2–7 objects (Venn diagrams are suppressed for \> 7).

- thresholds:

  (numeric vector, default `c(0, 1, 2, 3, 5, 10)`) Maximum dissimilarity
  thresholds in percent. `0` means exact string matches only; `10` means
  sequences within 10% divergence are grouped together.

- k:

  (integer, default `7`) K-mer size for building presence/absence
  profiles. Larger `k` increases specificity; smaller `k` increases
  sensitivity for divergent sequences.

- verbose:

  (logical, default `TRUE`) If `TRUE`, print progress messages.

## Value

A list of class `"refseq_comp_lpq_result"` with:

- `venn_plots`:

  Named list of ggVennDiagram objects (one per threshold, names follow
  `"<threshold>_pct"`), or `NULL` if ggVennDiagram is not installed or N
  \> 7.

- `results`:

  List of per-threshold results, each containing `threshold`,
  `venn_data` (named list of component IDs per set), and `n_shared`
  (named list of shared cluster counts per pair).

- `labels`:

  Character vector of phyloseq object names.

- `thresholds`:

  The thresholds used.

- `k`:

  The k-mer size used.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Author

Adrien Taudière

## Examples

``` r
lpq <- list_phyloseq(list(
  a = data_fungi_mini,
  b = prune_taxa(taxa_names(data_fungi_mini)[1:20], data_fungi_mini)
))
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: REPRODUCIBILITY
#> ℹ 137 common samples, 20 common taxa
#> ✔ list_phyloseq created (REPRODUCIBILITY)
res <- refseq_comp_lpq(lpq)
#> Building k-mer profiles (k = 7) for 2 objects...
#> Computing 1 pairwise similarity matrices...
#>   a vs b
```
