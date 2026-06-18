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
#> <list_phyloseq>
#>  @ phyloseq_list:List of 1
#>  .. $ a:Formal class 'phyloseq' [package "phyloseq"] with 5 slots
#>  ..  .. ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
#>  ..  .. .. .. ..@ .Data        : num [1:185, 1:1420] 1 23 1 281 260 0 5 1 5 5 ...
#>  ..  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#>  ..  .. .. .. .. .. ..$ : chr [1:185] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  .. .. .. .. .. ..$ : chr [1:1420] "ASV2" "ASV6" "ASV7" "ASV8" ...
#>  ..  .. .. .. ..@ taxa_are_rows: logi FALSE
#>  ..  .. .. .. ..$ dim     : int [1:2] 185 1420
#>  ..  .. .. .. ..$ dimnames:List of 2
#>  ..  .. .. .. .. ..$ : chr [1:185] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  .. .. .. .. ..$ : chr [1:1420] "ASV2" "ASV6" "ASV7" "ASV8" ...
#>  ..  .. ..@ tax_table:Formal class 'taxonomyTable' [package "phyloseq"] with 1 slot
#>  ..  .. .. .. ..@ .Data: chr [1:1420, 1:12] "Fungi" "Fungi" "Fungi" "Fungi" ...
#>  ..  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#>  ..  .. .. .. .. .. ..$ : chr [1:1420] "ASV2" "ASV6" "ASV7" "ASV8" ...
#>  ..  .. .. .. .. .. ..$ : chr [1:12] "Domain" "Phylum" "Class" "Order" ...
#>  ..  .. .. .. ..$ dim     : int [1:2] 1420 12
#>  ..  .. .. .. ..$ dimnames:List of 2
#>  ..  .. .. .. .. ..$ : chr [1:1420] "ASV2" "ASV6" "ASV7" "ASV8" ...
#>  ..  .. .. .. .. ..$ : chr [1:12] "Domain" "Phylum" "Class" "Order" ...
#>  ..  .. ..@ sam_data :'data.frame':  185 obs. of  7 variables:
#> Formal class 'sample_data' [package "phyloseq"] with 4 slots
#>  ..  .. .. .. ..@ .Data    :List of 7
#>  ..  .. .. .. .. ..$ : chr [1:185] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  .. .. .. .. ..$ : chr [1:185] "A10-005-B_S188" "A10-005-H_S189" "A10-005-M_S190" "A12-007_S191" ...
#>  ..  .. .. .. .. ..$ : chr [1:185] "A10-005" "A10-005" "A10-005" "A12-007" ...
#>  ..  .. .. .. .. ..$ : int [1:185] 188 189 190 191 2 3 4 5 6 7 ...
#>  ..  .. .. .. .. ..$ : chr [1:185] "Low" "High" "Middle" NA ...
#>  ..  .. .. .. .. ..$ : chr [1:185] "52" "52" "52" "28,4" ...
#>  ..  .. .. .. .. ..$ : int [1:185] 15 15 15 0 0 0 0 0 5 15 ...
#>  ..  .. .. .. ..@ names    : chr [1:7] "X" "Sample_names" "Tree_name" "Sample_id" ...
#>  ..  .. .. .. ..@ row.names: chr [1:185] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  .. .. .. ..@ .S3Class : chr "data.frame"
#>  ..  .. ..@ phy_tree : NULL
#>  ..  .. ..@ refseq   :Formal class 'DNAStringSet' [package "Biostrings"] with 5 slots
#>  ..  .. .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
#>  ..  .. .. .. .. .. ..@ xp_list                    :List of 1
#>  ..  .. .. .. .. .. .. ..$ :<pointer: (nil)> 
#>  ..  .. .. .. .. .. ..@ .link_to_cached_object_list:List of 1
#>  ..  .. .. .. .. .. .. ..$ :<environment: 0x62f57ff61378> 
#>  ..  .. .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
#>  ..  .. .. .. .. .. ..@ group          : int [1:1420] 1 1 1 1 1 1 1 1 1 1 ...
#>  ..  .. .. .. .. .. ..@ start          : int [1:1420] 1 313 614 963 1320 1620 1950 2241 2599 2939 ...
#>  ..  .. .. .. .. .. ..@ width          : int [1:1420] 312 301 349 357 300 330 291 358 340 305 ...
#>  ..  .. .. .. .. .. ..@ NAMES          : chr [1:1420] "ASV2" "ASV6" "ASV7" "ASV8" ...
#>  ..  .. .. .. .. .. ..@ elementType    : chr "ANY"
#>  ..  .. .. .. .. .. ..@ elementMetadata: NULL
#>  ..  .. .. .. .. .. ..@ metadata       : list()
#>  ..  .. .. .. ..@ elementType    : chr "DNAString"
#>  ..  .. .. .. ..@ elementMetadata: NULL
#>  ..  .. .. .. ..@ metadata       : list()
#>  @ summary_table: tibble [1 × 17] (S3: tbl_df/tbl/data.frame)
#>  $ name               : chr "a"
#>  $ n_samples          : int 185
#>  $ n_taxa             : int 1420
#>  $ n_sequences        : num 1839124
#>  $ n_occurence        : int 12499
#>  $ mean_seq_length    : num 318
#>  $ min_seq_length     : int 251
#>  $ mean_seq_per_sample: num 9941
#>  $ sd_seq_per_sample  : num 9678
#>  $ min_seq_per_sample : num 6
#>  $ max_seq_per_sample : num 39593
#>  $ mean_seq_per_taxon : num 1295
#>  $ sd_seq_per_taxon   : num 4581
#>  $ has_sam_data       : logi TRUE
#>  $ has_tax_table      : logi TRUE
#>  $ has_refseq         : logi TRUE
#>  $ has_phy_tree       : logi FALSE
#>  @ comparison   :List of 19
#>  .. $ type_of_comparison        : chr "REPRODUCIBILITY"
#>  .. $ type_description          : chr "Same samples, same primer/technology, same bioinformatics pipeline.\n  Used to test reproducibility of results "| __truncated__
#>  .. $ same_primer_seq_tech      : logi TRUE
#>  .. $ same_bioinfo_pipeline     : logi TRUE
#>  .. $ same_sam_data_structure   : logi NA
#>  .. $ same_samples              : logi TRUE
#>  .. $ nested_samples            : logi FALSE
#>  .. $ nesting_structure         : NULL
#>  .. $ same_taxa                 : logi TRUE
#>  .. $ n_common_samples          : int 185
#>  .. $ common_samples            : chr [1:185] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  .. $ n_common_taxa             : int 1420
#>  .. $ common_taxa               : chr [1:1420] "ASV2" "ASV6" "ASV7" "ASV8" ...
#>  .. $ shared_sam_data_modalities: list()
#>  .. $ all_have_sam_data         : logi TRUE
#>  .. $ all_have_tax_table        : logi TRUE
#>  .. $ all_have_refseq           : logi TRUE
#>  .. $ all_have_phy_tree         : logi FALSE
#>  .. $ refseq_comparison         : NULL
```
