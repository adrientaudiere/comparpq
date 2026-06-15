# Add fake sequences by shuffling existing ones in a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful to compute true negative values with functions
[`tc_metrics_mock()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock.md)
and
[`tc_metrics_mock_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock_vec.md).

Note the that the tax_table for additional sequences is full of NA and
that the corresponding otu_table is full of 0.

## Usage

``` r
add_shuffle_seq_pq(physeq, n_fake = NULL, prop_fake = NULL, prefix = "fake_")
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- n_fake:

  (integer, default NULL) A number of fake taxa to add. Must specify
  either `n_fake` or `prop_fake`, not both.

- prop_fake:

  (numeric, default NULL) A proportion of fake taxa to add. Must specify
  either `n_fake` or `prop_fake`, not both.

- prefix:

  (character, default "fake\_") A prefix to add to the taxa name.

## Value

A phyloseq object

## Author

Adrien Taudière

## Examples

``` r
d_fake_F <- data_fungi_mini |>
  add_shuffle_seq_pq(prop_fake = 0.1)
#> Taxa are now in columns.
#> Cleaning suppress 0 taxa (  ) and 65 sample(s) ( A10-005-B_S188_MERGED.fastq.gz / A12-007-B_S2_MERGED.fastq.gz / A15-004_S3_MERGED.fastq.gz / A8-005_S4_MERGED.fastq.gz / AC27-013_S7_MERGED.fastq.gz / AD26-005-H_S10_MERGED.fastq.gz / AD30-ABMX-M_S12_MERGED.fastq.gz / B18-006-B_S19_MERGED.fastq.gz / BB19-006-H_S22_MERGED.fastq.gz / BB6-019-M_S25_MERGED.fastq.gz / BH9-021_S33_MERGED.fastq.gz / BL7-006-H_S37_MERGED.fastq.gz / BO8-005_S42_MERGED.fastq.gz / BP11-001-M_S45_MERGED.fastq.gz / BP12-025-B_S46_MERGED.fastq.gz / BQ3-019_S48_MERGED.fastq.gz / BQ4-018-B_S49_MERGED.fastq.gz / BQ4-018-H_S50_MERGED.fastq.gz / BR8-005_S53_MERGED.fastq.gz / BT-006-M_S55_MERGED.fastq.gz / BT7-006_S56_MERGED.fastq.gz / C1-001_S61_MERGED.fastq.gz / C9-005_S65_MERGED.fastq.gz / CB8-019-B_S69_MERGED.fastq.gz / CB8-019-H_S70_MERGED.fastq.gz / CB8-019-M_S71_MERGED.fastq.gz / D18-003-B_S78_MERGED.fastq.gz / D61-010-B_S82_MERGED.fastq.gz / D9-027-H_S84_MERGED.fastq.gz / D9-027-M_S85_MERGED.fastq.gz / DP4-ABM001_S90_MERGED.fastq.gz / DU3-045-B_S94_MERGED.fastq.gz / DW4-007_S95_MERGED.fastq.gz / DY5-004-B_S96_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / F7-015-M_S106_MERGED.fastq.gz / FOMES19-H_S108_MERGED.fastq.gz / FOMES19-M_S109_MERGED.fastq.gz / H10-018-M_S110_MERGED.fastq.gz / H24-NVABM1-H_S111_MERGED.fastq.gz / J18-004-B_S114_MERGED.fastq.gz / J18-004-M_S116_MERGED.fastq.gz / L19X-B_S119_MERGED.fastq.gz / L19X-H_S120_MERGED.fastq.gz / L19X-M_S121_MERGED.fastq.gz / L23-002-M_S124_MERGED.fastq.gz / N19X-M_S128_MERGED.fastq.gz / NVABM-0397_S138_MERGED.fastq.gz / O20-X-M_S141_MERGED.fastq.gz / O24-003-B_S145_MERGED.fastq.gz / O27-012_S151_MERGED.fastq.gz / O9-005-B_S152_MERGED.fastq.gz / P27-015-M_S154_MERGED.fastq.gz / Q27-ABM003-B_S156_MERGED.fastq.gz / T28-011_S161_MERGED.fastq.gz / W26-001-B_S165_MERGED.fastq.gz / W9-025-M_S169_MERGED.fastq.gz / X24-009-B_S170_MERGED.fastq.gz / X24-009-H_S171_MERGED.fastq.gz / X24-009-M_S172_MERGED.fastq.gz / X24-010_S173_MERGED.fastq.gz / X29-004-B_S174_MERGED.fastq.gz / Y28-002-B_S178_MERGED.fastq.gz / Z29-001-H_S185_MERGED.fastq.gz / Z30-ABM560-M_S187_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 45
#> Number of filtered-out ASV 41
#> Number of kept ASV 4
#> Number of kept samples 72
#> Taxa are now in rows.
ntaxa(d_fake_F) - ntaxa(data_fungi_mini)
#> [1] 4
```
