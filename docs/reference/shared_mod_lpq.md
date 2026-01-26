# Display shared sample_data modalities

Create a table showing which sample_data variable modalities are shared
across phyloseq objects in a list_phyloseq.

Create a table showing which sample_data variable modalities are shared
across phyloseq objects in a list_phyloseq.

## Usage

``` r
shared_mod_lpq(x, max_modalities = NULL)

shared_mod_lpq(x, max_modalities = NULL)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- max_modalities:

  (integer, default 10) Maximum number of modalities to display per
  variable.

## Value

A tibble or NULL if no shared modalities exist

A tibble or NULL if no shared modalities exist

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Examples

``` r
lpq <- list_phyloseq(list(data1 = data_fungi, data2 = data_fungi_mini))

shared_mod_lpq(lpq)
#> # A tibble: 5 × 3
#>   Variable     N_shared Shared_modalities                                       
#>   <chr>           <int> <chr>                                                   
#> 1 X                 137 A10-005-B_S188_MERGED.fastq.gz, A10-005-H_S189_MERGED.f…
#> 2 Sample_names      137 A10-005-B_S188, A10-005-H_S189, A10-005-M_S190, A12-007…
#> 3 Tree_name         100 A10-005, A12-007, A15-004, A8-005, AB29-abm-X, AC27-013…
#> 4 Height              4 Low, High, Middle, NA                                   
#> 5 Diameter           70 52, 28,4, 30,7, 32,8, 33,3, 99, 32, 55,4, 115,5, -, 10,…
shared_mod_lpq(lpq, 10)
#> # A tibble: 5 × 3
#>   Variable     N_shared Shared_modalities                                       
#>   <chr>           <int> <chr>                                                   
#> 1 X                 137 A10-005-B_S188_MERGED.fastq.gz, A10-005-H_S189_MERGED.f…
#> 2 Sample_names      137 A10-005-B_S188, A10-005-H_S189, A10-005-M_S190, A12-007…
#> 3 Tree_name         100 A10-005, A12-007, A15-004, A8-005, AB29-abm-X, AC27-013…
#> 4 Height              4 Low, High, Middle, NA                                   
#> 5 Diameter           70 52, 28,4, 30,7, 32,8, 33,3, 99, 32, 55,4, 115,5, -, ...…
lpq <- list_phyloseq(list(data1 = data_fungi, data2 = data_fungi_mini))

shared_mod_lpq(lpq)
#> # A tibble: 5 × 3
#>   Variable     N_shared Shared_modalities                                       
#>   <chr>           <int> <chr>                                                   
#> 1 X                 137 A10-005-B_S188_MERGED.fastq.gz, A10-005-H_S189_MERGED.f…
#> 2 Sample_names      137 A10-005-B_S188, A10-005-H_S189, A10-005-M_S190, A12-007…
#> 3 Tree_name         100 A10-005, A12-007, A15-004, A8-005, AB29-abm-X, AC27-013…
#> 4 Height              4 Low, High, Middle, NA                                   
#> 5 Diameter           70 52, 28,4, 30,7, 32,8, 33,3, 99, 32, 55,4, 115,5, -, 10,…
shared_mod_lpq(lpq, 10)
#> # A tibble: 5 × 3
#>   Variable     N_shared Shared_modalities                                       
#>   <chr>           <int> <chr>                                                   
#> 1 X                 137 A10-005-B_S188_MERGED.fastq.gz, A10-005-H_S189_MERGED.f…
#> 2 Sample_names      137 A10-005-B_S188, A10-005-H_S189, A10-005-M_S190, A12-007…
#> 3 Tree_name         100 A10-005, A12-007, A15-004, A8-005, AB29-abm-X, AC27-013…
#> 4 Height              4 Low, High, Middle, NA                                   
#> 5 Diameter           70 52, 28,4, 30,7, 32,8, 33,3, 99, 32, 55,4, 115,5, -, ...…
```
