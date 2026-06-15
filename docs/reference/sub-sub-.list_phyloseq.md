# Extract a single phyloseq object from a `list_phyloseq`

Extract a single phyloseq object from a `list_phyloseq`

## Usage

``` r
# S3 method for class 'list_phyloseq'
x[[i]]
```

## Arguments

- x:

  A `list_phyloseq` object.

- i:

  Index (integer or character) selecting one element.

## Value

The selected `phyloseq` object.

## Examples

``` r
lpq <- list_phyloseq(list(a = data_fungi, b = data_fungi))
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: REPRODUCIBILITY
#> ℹ 185 common samples, 1420 common taxa
#> ✔ list_phyloseq created (REPRODUCIBILITY)
lpq[[1]]
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
