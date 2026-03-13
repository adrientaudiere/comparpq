# Apply a function to all phyloseq objects in a list_phyloseq

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Applies a function to each phyloseq object in a list_phyloseq and
returns a new list_phyloseq with the transformed objects. This is useful
for applying cleaning, filtering, or transformation functions uniformly
across all phyloseq objects.

## Usage

``` r
apply_to_lpq(x, .f, ..., verbose = TRUE)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- .f:

  (function, required) A function to apply to each phyloseq object. The
  function must take a phyloseq object as its first argument and return
  a phyloseq object.

- ...:

  Additional arguments passed to `.f`.

- verbose:

  (logical, default TRUE) If TRUE, print information about the
  transformation applied to each phyloseq object.

## Value

A new list_phyloseq object with the transformed phyloseq objects. The
comparison parameters (`same_primer_seq_tech`, `same_bioinfo_pipeline`)
are preserved from the original object.

## Details

Common functions to apply include:

- [`MiscMetabar::clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.html):
  Remove empty samples and taxa

- `phyloseq::taxa_as_rows()`: Ensure taxa are rows in otu_table

- [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html):
  Rarefy to even depth

- [`phyloseq::transform_sample_counts()`](https://rdrr.io/pkg/phyloseq/man/transformcounts.html):
  Transform counts (e.g., relative abundance)

- [`phyloseq::subset_taxa()`](https://rdrr.io/pkg/phyloseq/man/subset_taxa-methods.html):
  Filter taxa based on criteria

- [`phyloseq::subset_samples()`](https://rdrr.io/pkg/phyloseq/man/subset_samples-methods.html):
  Filter samples based on criteria

## See also

[`filter_common_lpq()`](https://adrientaudiere.github.io/comparpq/reference/filter_common_lpq.md),
[`update_list_phyloseq()`](https://adrientaudiere.github.io/comparpq/reference/update_list_phyloseq.md)

## Author

Adrien Taudière

## Examples

``` r
lpq <- list_phyloseq(list(fungi = data_fungi, fungi_mini = data_fungi_mini))
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

# Apply clean_pq to all phyloseq objects
lpq_clean <- apply_to_lpq(lpq, MiscMetabar::clean_pq)
#>   fungi: 185 -> 185 samples, 1420 -> 1420 taxa
#>   fungi_mini: 137 -> 137 samples, 45 -> 45 taxa
#> Applied `MiscMetabar::clean_pq` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

# Apply taxa_as_rows
lpq_rows <- apply_to_lpq(lpq, MiscMetabar::taxa_as_rows)
#> Taxa are now in rows.
#>   fungi: 185 -> 185 samples, 1420 -> 1420 taxa
#>   fungi_mini: 137 -> 137 samples, 45 -> 45 taxa
#> Applied `MiscMetabar::taxa_as_rows` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

# Apply rarefy_even_depth with a specific rngseed
lpq_rar <- apply_to_lpq(lpq, rarefy_even_depth, rngseed = 21)
#> `set.seed(21)` was used to initialize repeatable random subsampling.
#> Please record this for your records so others can reproduce.
#> Try `set.seed(21); .Random.seed` for the full vector
#> ...
#> 1000OTUs were removed because they are no longer 
#> present in any sample after random subsampling
#> ...
#>   fungi: 185 -> 185 samples, 1420 -> 420 taxa
#> `set.seed(21)` was used to initialize repeatable random subsampling.
#> Please record this for your records so others can reproduce.
#> Try `set.seed(21); .Random.seed` for the full vector
#> ...
#> 7OTUs were removed because they are no longer 
#> present in any sample after random subsampling
#> ...
#>   fungi_mini: 137 -> 137 samples, 45 -> 38 taxa
#> Applied `rarefy_even_depth` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 35 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

lpq_rar
#> <comparpq::list_phyloseq>
#>  @ phyloseq_list:List of 2
#>  .. $ fungi     :Formal class 'phyloseq' [package "phyloseq"] with 5 slots
#>  ..  .. ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
#>  ..  .. .. .. ..@ .Data        : num [1:185, 1:420] 0 0 0 0 1 0 0 0 0 0 ...
#>  ..  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#>  ..  .. .. .. .. .. ..$ : chr [1:185] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  .. .. .. .. .. ..$ : chr [1:420] "ASV2" "ASV6" "ASV7" "ASV8" ...
#>  ..  .. .. .. ..@ taxa_are_rows: logi FALSE
#>  ..  .. .. .. ..$ dim     : int [1:2] 185 420
#>  ..  .. .. .. ..$ dimnames:List of 2
#>  ..  .. .. .. .. ..$ : chr [1:185] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  .. .. .. .. ..$ : chr [1:420] "ASV2" "ASV6" "ASV7" "ASV8" ...
#>  ..  .. ..@ tax_table:Formal class 'taxonomyTable' [package "phyloseq"] with 1 slot
#>  ..  .. .. .. ..@ .Data: chr [1:420, 1:12] "Fungi" "Fungi" "Fungi" "Fungi" ...
#>  ..  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#>  ..  .. .. .. .. .. ..$ : chr [1:420] "ASV2" "ASV6" "ASV7" "ASV8" ...
#>  ..  .. .. .. .. .. ..$ : chr [1:12] "Domain" "Phylum" "Class" "Order" ...
#>  ..  .. .. .. ..$ dim     : int [1:2] 420 12
#>  ..  .. .. .. ..$ dimnames:List of 2
#>  ..  .. .. .. .. ..$ : chr [1:420] "ASV2" "ASV6" "ASV7" "ASV8" ...
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
#>  ..  .. .. .. .. .. .. ..$ :<externalptr> 
#>  ..  .. .. .. .. .. ..@ .link_to_cached_object_list:List of 1
#>  ..  .. .. .. .. .. .. ..$ :<environment: 0x556d79d199d8> 
#>  ..  .. .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
#>  ..  .. .. .. .. .. ..@ group          : int [1:420] 1 1 1 1 1 1 1 1 1 1 ...
#>  ..  .. .. .. .. .. ..@ start          : int [1:420] 1 313 614 963 1320 1620 1950 2241 2599 2939 ...
#>  ..  .. .. .. .. .. ..@ width          : int [1:420] 312 301 349 357 300 330 291 358 340 305 ...
#>  ..  .. .. .. .. .. ..@ NAMES          : chr [1:420] "ASV2" "ASV6" "ASV7" "ASV8" ...
#>  ..  .. .. .. .. .. ..@ elementType    : chr "ANY"
#>  ..  .. .. .. .. .. ..@ elementMetadata: NULL
#>  ..  .. .. .. .. .. ..@ metadata       : list()
#>  ..  .. .. .. ..@ elementType    : chr "DNAString"
#>  ..  .. .. .. ..@ elementMetadata: NULL
#>  ..  .. .. .. ..@ metadata       : list()
#>  .. $ fungi_mini:Formal class 'phyloseq' [package "phyloseq"] with 5 slots
#>  ..  .. ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
#>  ..  .. .. .. ..@ .Data        : num [1:38, 1:137] 0 0 0 0 0 1 0 0 0 0 ...
#>  ..  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#>  ..  .. .. .. .. .. ..$ : chr [1:38] "ASV7" "ASV8" "ASV12" "ASV18" ...
#>  ..  .. .. .. .. .. ..$ : chr [1:137] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  .. .. .. ..@ taxa_are_rows: logi TRUE
#>  ..  .. .. .. ..$ dim     : int [1:2] 38 137
#>  ..  .. .. .. ..$ dimnames:List of 2
#>  ..  .. .. .. .. ..$ : chr [1:38] "ASV7" "ASV8" "ASV12" "ASV18" ...
#>  ..  .. .. .. .. ..$ : chr [1:137] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  .. ..@ tax_table:Formal class 'taxonomyTable' [package "phyloseq"] with 1 slot
#>  ..  .. .. .. ..@ .Data: chr [1:38, 1:12] "Fungi" "Fungi" "Fungi" "Fungi" ...
#>  ..  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#>  ..  .. .. .. .. .. ..$ : chr [1:38] "ASV7" "ASV8" "ASV12" "ASV18" ...
#>  ..  .. .. .. .. .. ..$ : chr [1:12] "Domain" "Phylum" "Class" "Order" ...
#>  ..  .. .. .. ..$ dim     : int [1:2] 38 12
#>  ..  .. .. .. ..$ dimnames:List of 2
#>  ..  .. .. .. .. ..$ : chr [1:38] "ASV7" "ASV8" "ASV12" "ASV18" ...
#>  ..  .. .. .. .. ..$ : chr [1:12] "Domain" "Phylum" "Class" "Order" ...
#>  ..  .. ..@ sam_data :'data.frame':  137 obs. of  7 variables:
#> Formal class 'sample_data' [package "phyloseq"] with 4 slots
#>  ..  .. .. .. ..@ .Data    :List of 7
#>  ..  .. .. .. .. ..$ : chr [1:137] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  .. .. .. .. ..$ : chr [1:137] "A10-005-B_S188" "A10-005-H_S189" "A10-005-M_S190" "A12-007_S191" ...
#>  ..  .. .. .. .. ..$ : chr [1:137] "A10-005" "A10-005" "A10-005" "A12-007" ...
#>  ..  .. .. .. .. ..$ : int [1:137] 188 189 190 191 2 3 4 6 7 8 ...
#>  ..  .. .. .. .. ..$ : chr [1:137] "Low" "High" "Middle" NA ...
#>  ..  .. .. .. .. ..$ : chr [1:137] "52" "52" "52" "28,4" ...
#>  ..  .. .. .. .. ..$ : int [1:137] 15 15 15 0 0 0 0 5 15 15 ...
#>  ..  .. .. .. ..@ names    : chr [1:7] "X" "Sample_names" "Tree_name" "Sample_id" ...
#>  ..  .. .. .. ..@ row.names: chr [1:137] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  .. .. .. ..@ .S3Class : chr "data.frame"
#>  ..  .. ..@ phy_tree : NULL
#>  ..  .. ..@ refseq   :Formal class 'DNAStringSet' [package "Biostrings"] with 5 slots
#>  ..  .. .. .. ..@ pool           :Formal class 'SharedRaw_Pool' [package "XVector"] with 2 slots
#>  ..  .. .. .. .. .. ..@ xp_list                    :List of 1
#>  ..  .. .. .. .. .. .. ..$ :<externalptr> 
#>  ..  .. .. .. .. .. ..@ .link_to_cached_object_list:List of 1
#>  ..  .. .. .. .. .. .. ..$ :<environment: 0x556d80c990f8> 
#>  ..  .. .. .. ..@ ranges         :Formal class 'GroupedIRanges' [package "XVector"] with 7 slots
#>  ..  .. .. .. .. .. ..@ group          : int [1:38] 1 1 1 1 1 1 1 1 1 1 ...
#>  ..  .. .. .. .. .. ..@ start          : int [1:38] 614 963 1620 2241 3898 4248 4595 5297 5903 6589 ...
#>  ..  .. .. .. .. .. ..@ width          : int [1:38] 349 357 330 358 350 347 345 303 383 358 ...
#>  ..  .. .. .. .. .. ..@ NAMES          : chr [1:38] "ASV7" "ASV8" "ASV12" "ASV18" ...
#>  ..  .. .. .. .. .. ..@ elementType    : chr "ANY"
#>  ..  .. .. .. .. .. ..@ elementMetadata: NULL
#>  ..  .. .. .. .. .. ..@ metadata       : list()
#>  ..  .. .. .. ..@ elementType    : chr "DNAString"
#>  ..  .. .. .. ..@ elementMetadata: NULL
#>  ..  .. .. .. ..@ metadata       : list()
#>  @ summary_table: tibble [2 × 17] (S3: tbl_df/tbl/data.frame)
#>  $ name               : chr [1:2] "fungi" "fungi_mini"
#>  $ n_samples          : int [1:2] 185 137
#>  $ n_taxa             : int [1:2] 420 38
#>  $ n_sequences        : num [1:2] 1110 137
#>  $ n_occurence        : int [1:2] 695 137
#>  $ mean_seq_length    : num [1:2] 323 344
#>  $ min_seq_length     : int [1:2] 255 303
#>  $ mean_seq_per_sample: num [1:2] 6 1
#>  $ sd_seq_per_sample  : num [1:2] 0 0
#>  $ min_seq_per_sample : num [1:2] 6 1
#>  $ max_seq_per_sample : num [1:2] 6 1
#>  $ mean_seq_per_taxon : num [1:2] 2.64 3.61
#>  $ sd_seq_per_taxon   : num [1:2] 3.62 4.02
#>  $ has_sam_data       : logi [1:2] TRUE TRUE
#>  $ has_tax_table      : logi [1:2] TRUE TRUE
#>  $ has_refseq         : logi [1:2] TRUE TRUE
#>  $ has_phy_tree       : logi [1:2] FALSE FALSE
#>  @ comparison   :List of 19
#>  .. $ type_of_comparison        : chr "NESTED_ROBUSTNESS"
#>  .. $ type_description          : chr "Nested samples detected (one phyloseq derived from another).\n  Nesting: fungi_mini (137 samples) is nested in "| __truncated__
#>  .. $ same_primer_seq_tech      : logi TRUE
#>  .. $ same_bioinfo_pipeline     : logi TRUE
#>  .. $ same_sam_data_structure   : logi TRUE
#>  .. $ same_samples              : logi FALSE
#>  .. $ nested_samples            : logi TRUE
#>  .. $ nesting_structure         : chr "fungi_mini (137 samples) is nested in fungi (185 samples)"
#>  .. $ same_taxa                 : logi FALSE
#>  .. $ n_common_samples          : int 137
#>  .. $ common_samples            : chr [1:137] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  .. $ n_common_taxa             : int 35
#>  .. $ common_taxa               : chr [1:35] "ASV7" "ASV8" "ASV12" "ASV18" ...
#>  .. $ shared_sam_data_modalities:List of 5
#>  ..  ..$ X           : chr [1:137] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  ..  ..$ Sample_names: chr [1:137] "A10-005-B_S188" "A10-005-H_S189" "A10-005-M_S190" "A12-007_S191" ...
#>  ..  ..$ Tree_name   : chr [1:100] "A10-005" "A12-007" "A15-004" "A8-005" ...
#>  ..  ..$ Height      : chr [1:4] "Low" "High" "Middle" NA
#>  ..  ..$ Diameter    : chr [1:70] "52" "28,4" "30,7" "32,8" ...
#>  .. $ all_have_sam_data         : logi TRUE
#>  .. $ all_have_tax_table        : logi TRUE
#>  .. $ all_have_refseq           : logi TRUE
#>  .. $ all_have_phy_tree         : logi FALSE
#>  .. $ refseq_comparison         :List of 1
#>  ..  ..$ fungi_vs_fungi_mini:List of 14
#>  ..  .. ..$ name1              : chr "fungi"
#>  ..  .. ..$ name2              : chr "fungi_mini"
#>  ..  .. ..$ k                  : num 5
#>  ..  .. ..$ n_taxa             : Named int [1:2] 420 38
#>  ..  .. .. ..- attr(*, "names")= chr [1:2] "fungi" "fungi_mini"
#>  ..  .. ..$ shared_names       : chr [1:35] "ASV7" "ASV8" "ASV12" "ASV18" ...
#>  ..  .. ..$ unique_names_1     : chr [1:385] "ASV2" "ASV6" "ASV10" "ASV13" ...
#>  ..  .. ..$ unique_names_2     : chr [1:3] "ASV68" "ASV83" "ASV93"
#>  ..  .. ..$ shared_seqs        : chr [1:35] "AAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACCTTGCGCCCTTTGGTATTCCGAAGGGCACACCTGTTTGAGTGTCGTGAAATT"| __truncated__ "AAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACCTTGCGCCCTTTGGTATTCCGAAGGGCACACCTGTTTGAGTGTCGTGAAATT"| __truncated__ "AAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACCTTGCGCTCCTTGGTATTCCGAGGAGCATGCCTGTTTGAGTGTCGTGAAATT"| __truncated__ "AAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACCTTGCGCCCTTTGGTATTCCGAAGGGCACACCTGTTTGAGTGTCGTGAAATT"| __truncated__ ...
#>  ..  .. ..$ unique_seqs_1      : chr [1:385] "AAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCATTAGTATTCTAGTGGGCATGCCTGTTCGAGCGTCATTTTCAA"| __truncated__ "AAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCATTAGTATTCTAGTGGGCATGCCTGTTCGAGCGTCATTTCAAC"| __truncated__ "AAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCTAGTACTCTAGCGGGCATGCCTGTTCGAGCGTCATTTCAAC"| __truncated__ "AAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCTCTGGTATTCCGGGGGGCATGCCTGTTCGAGCGTCATTACAAC"| __truncated__ ...
#>  ..  .. ..$ unique_seqs_2      : chr [1:3] "AAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACCTTGCGCTCCTTGGTATTCCGAGGAGCATGCCTGTTTGAGTGTCATTAAATT"| __truncated__ "AAATGCGATAAGTAATGTGAATTGCAGTATTCAGTGAATCATCGAATCTTTGAACGCATCTTGCGCCCTTTGGTATTCCGAAGGGCATGCCTGTTTGAGTGTCATGAAAAT"| __truncated__ "AAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCGTCGAATCTTTGAACGCACCTTGCGCCCTTTGGTGTTCCGAAGGGCACACCTGTTTGAGTGTCGTGAAATT"| __truncated__
#>  ..  .. ..$ same_name_diff_seq :'data.frame':    0 obs. of  3 variables:
#>  ..  .. .. ..$ taxa_name  : chr(0) 
#>  ..  .. .. ..$ seq_physeq1: chr(0) 
#>  ..  .. .. ..$ seq_physeq2: chr(0) 
#>  ..  .. ..$ same_seq_diff_name :'data.frame':    0 obs. of  3 variables:
#>  ..  .. .. ..$ sequence    : chr(0) 
#>  ..  .. .. ..$ name_physeq1: chr(0) 
#>  ..  .. .. ..$ name_physeq2: chr(0) 
#>  ..  .. ..$ mean_nn_kmer_dist_1: num 0.241
#>  ..  .. ..$ mean_nn_kmer_dist_2: num 0.00391
#>  ..  .. ..- attr(*, "class")= chr "compare_refseq"

# Transform to relative abundance
lpq_rel <- apply_to_lpq(
  lpq,
  phyloseq::transform_sample_counts,
  function(x) x / sum(x)
)
#>   fungi: 185 -> 185 samples, 1420 -> 1420 taxa
#>   fungi_mini: 137 -> 137 samples, 45 -> 45 taxa
#> Applied `phyloseq::transform_sample_counts` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)

# Chain multiple transformations
lpq_processed <- lpq |>
  apply_to_lpq(MiscMetabar::clean_pq) |>
  apply_to_lpq(MiscMetabar::taxa_as_rows)
#>   fungi: 185 -> 185 samples, 1420 -> 1420 taxa
#>   fungi_mini: 137 -> 137 samples, 45 -> 45 taxa
#> Applied `MiscMetabar::clean_pq` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)
#> Taxa are now in rows.
#>   fungi: 185 -> 185 samples, 1420 -> 1420 taxa
#>   fungi_mini: 137 -> 137 samples, 45 -> 45 taxa
#> Applied `MiscMetabar::taxa_as_rows` to 2 phyloseq objects.
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: NESTED_ROBUSTNESS
#> ℹ 137 common samples, 45 common taxa
#> ✔ list_phyloseq created (NESTED_ROBUSTNESS)
```
