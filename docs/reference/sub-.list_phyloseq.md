# Subset a `list_phyloseq` object

Subset a `list_phyloseq` object

## Usage

``` r
# S3 method for class 'list_phyloseq'
x[i]
```

## Arguments

- x:

  A `list_phyloseq` object.

- i:

  Index (integer or character) selecting elements.

## Value

A new `list_phyloseq` built from the selected phyloseq objects.

## Examples

``` r
lpq <- list_phyloseq(list(a = data_fungi, b = data_fungi))
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: REPRODUCIBILITY
#> ℹ 185 common samples, 1420 common taxa
#> ✔ list_phyloseq created (REPRODUCIBILITY)
lpq[1]
#> list_phyloseq object with 1 phyloseq objects
#> 
#> --- Summary ---
#> # A tibble: 1 × 17
#>   name  n_samples n_taxa n_sequences n_occurence mean_seq_length min_seq_length
#>   <chr>     <int>  <int>       <dbl>       <int>           <dbl>          <int>
#> 1 a           185   1420     1839124       12499            318.            251
#> # ℹ 10 more variables: mean_seq_per_sample <dbl>, sd_seq_per_sample <dbl>,
#> #   min_seq_per_sample <dbl>, max_seq_per_sample <dbl>,
#> #   mean_seq_per_taxon <dbl>, sd_seq_per_taxon <dbl>, has_sam_data <lgl>,
#> #   has_tax_table <lgl>, has_refseq <lgl>, has_phy_tree <lgl>
#> 
#> --- Comparison characteristics ---
#> Type of comparison: REPRODUCIBILITY 
#> Same samples, same primer/technology, same bioinformatics pipeline.
#>   Used to test reproducibility of results when running
#>   the exact same analysis multiple times. 
#> 
#> Same primer/seq tech: TRUE 
#> Same bioinfo pipeline: TRUE 
#> Same sample_data structure: NA 
#> Same samples: TRUE 
#> Nested samples: FALSE 
#> Same taxa: TRUE 
#> Common samples: 185 
#> Common taxa: 1420 
```
