# Find taxa whose reference sequences match primer sequences

Searches every sequence in `@refseq` for occurrences of the supplied
primers (forward **and** reverse complement) using IUPAC-aware matching
via
[`Biostrings::vcountPattern()`](https://rdrr.io/pkg/Biostrings/man/matchPattern.html).
Returns a data frame of taxa that match at least one primer, which can
be passed directly to
[`tidypq::filter_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_taxa_pq.html)
to prune them from the phyloseq object.

## Usage

``` r
find_primers_pq(physeq, primers, verbose = TRUE)
```

## Arguments

- physeq:

  (required) A
  [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
  object with a populated `@refseq` slot.

- primers:

  (required) A named character vector of primer sequences. IUPAC
  ambiguity codes (M, R, Y, S, W, K, B, D, H, V, N) are supported.

- verbose:

  (logical, default `TRUE`) If `TRUE`, print a summary message.

## Value

A `data.frame` (or `NULL` if no matches) with columns:

- `taxon`:

  Character. Taxa name as in `taxa_names(physeq)`.

- `matched_primers`:

  Character. Comma-separated names of matching primers.

- `n_reads`:

  Numeric. Total read count across all samples
  (`phyloseq::taxa_sums(physeq)`).

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Author

Adrien Taudière

## Examples

``` r
primers <- c(
  mcrA_fwd = "GGTGGTGTMGGDTTCACMCARTA",
  mcrA_rev = "CGTTCATBGCGTAGTTVGGRTAGT"
)
bad <- find_primers_pq(data_fungi_mini, primers)
#> 0 taxa matched at least one primer out of 45.
bad
#> NULL

# Prune contaminated taxa (requires tidypq):
# if (!is.null(bad)) {
#   tidypq::filter_taxa_pq(
#     data_fungi_mini,
#     !taxa_names(data_fungi_mini) %in% bad$taxon,
#     clean_phyloseq_object = TRUE
#   )
# }
```
