# Add external sequences to a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful to compute true negative values with functions
[`tc_metrics_mock()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock.md)
and
[`tc_metrics_mock_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock_vec.md).

Note the that the tax_table for additional sequences is full of NA and
that the corresponding otu_table is full of 0.

## Usage

``` r
add_external_seq_pq(physeq, ext_seqs, prefix = "external_")
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- ext_seqs:

  (DNAStringSet, required) A DNAStringSet object containing external
  sequences to add.

- prefix:

  (character, default "external\_") A prefix to add to the taxa name.

## Value

A phyloseq object

## Author

Adrien Taudière

## Examples

``` r
add_external_seq_pq(
  data_fungi_mini,
  Biostrings::readDNAStringSet(
    system.file("extdata/ex_little.fasta", package = "MiscMetabar")
  )
)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 47 taxa and 137 samples ]
#> sample_data() Sample Data:       [ 137 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 47 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 47 reference sequences ]
```
