# Estimation statistics for categorical comparisons on a phyloseq object

Computes diversity metrics (Hill numbers by default) per sample and
compares them across groups defined by a categorical variable using
estimation statistics (effect sizes + bootstrap confidence intervals)
via the dabestr package.

This approach replaces traditional p-value-based hypothesis testing with
Gardner-Altman or Cumming estimation plots, following the estimation
statistics framework (Ho et al. 2019).

## Usage

``` r
estim_diff_pq(
  physeq,
  fact,
  hill_scales = c(0, 1, 2),
  custom_fn = NULL,
  effect_type = "cohens_d",
  idx = NULL,
  ci = 95,
  resamples = 5000,
  na_remove = TRUE,
  ...
)
```

## Arguments

- physeq:

  (phyloseq, required) A phyloseq object.

- fact:

  (character, required) The name of a categorical column in
  `sample_data` to use as the grouping factor.

- hill_scales:

  (numeric vector, default `c(0, 1, 2)`) The q values for Hill number
  computation: 0 = richness, 1 = Shannon exponential, 2 = inverse
  Simpson.

- custom_fn:

  (function, default NULL) An optional custom diversity function. Must
  take a phyloseq object and return a named numeric vector (names =
  sample names) or a data.frame with one row per sample. If provided,
  `hill_scales` is ignored.

- effect_type:

  (character, default `"mean_diff"`) The type of effect size to compute.
  One of: `"mean_diff"`, `"median_diff"`, `"cohens_d"`, `"hedges_g"`,
  `"cliffs_delta"`.

- idx:

  (list or character vector, default NULL) The group ordering for
  comparisons. If NULL, uses factor levels with first level as control.
  For 2 groups: `c("Control", "Treatment")`. For 3+ groups:
  `list(c("Ctrl", "T1", "T2"))`.

- ci:

  (numeric, default 95) Confidence interval level (0-100).

- resamples:

  (integer, default 5000) Number of bootstrap resamples.

- na_remove:

  (logical, default TRUE) If TRUE, samples with NA in `fact` are removed
  before analysis.

- ...:

  Additional arguments passed to dabestr plotting functions.

## Value

A list of class `"estim_diff_pq_result"` with components:

- data:

  The diversity data.frame used for analysis

- dabest_objects:

  A named list of dabestr objects (one per metric)

- plots:

  A named list of dabestr plots (one per metric)

- summary:

  A tibble summarizing all effect sizes and CIs with columns: `metric`,
  `comparison`, `effect_size`, `ci_lower`, `ci_upper`,
  `pvalue_permtest`, `pvalue_welch`, `pvalue_mann_whitney`

- effect_type:

  The effect size type used

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The function uses the dabestr package to produce estimation plots. For
two groups, Gardner-Altman plots are produced (`float_contrast = TRUE`).
For three or more groups, Cumming plots are used
(`float_contrast = FALSE`).

## References

Ho, J., Tumkaya, T., Aryal, S., Choi, H., & Claridge-Chang, A. (2019).
Moving beyond P values: data analysis with estimation graphics. *Nature
Methods*, 16(7), 565-566.

## See also

[`estim_cor_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_cor_pq.md),
[`estim_diff_lpq()`](https://adrientaudiere.github.io/comparpq/reference/estim_diff_lpq.md),
[`adonis_lpq()`](https://adrientaudiere.github.io/comparpq/reference/adonis_lpq.md)

## Examples

``` r
library(phyloseq)
data("data_fungi", package = "MiscMetabar")

pq <- subset_samples(data_fungi, !is.na(Height))
pq <- clean_pq(pq)
#> Cleaning suppress 144 taxa and 0 samples.

res <- estim_diff_pq(pq, fact = "Height")
res
#> $data
#>                                   Hill_0    Hill_1    Hill_2
#> A10-005-B_S188_MERGED.fastq.gz        67  7.967879  3.177791
#> A10-005-H_S189_MERGED.fastq.gz        70  2.701628  1.908759
#> A10-005-M_S190_MERGED.fastq.gz        54 18.735415 13.172882
#> A12-007-B_S2_MERGED.fastq.gz          47  8.918381  5.224328
#> AB29-ABMX-H_S6_MERGED.fastq.gz        70 10.225479  5.504729
#> AD26-005-B_S9_MERGED.fastq.gz        150  8.934307  3.358402
#> AD26-005-H_S10_MERGED.fastq.gz       110  4.931153  2.385448
#> AD26-005-M_S11_MERGED.fastq.gz       125 27.640250 16.450947
#> AD30-ABMX-M_S12_MERGED.fastq.gz       61  8.148864  4.476903
#> AD32-007-M_S13_MERGED.fastq.gz       101  8.391872  3.323816
#> ADABM30X-B_S14_MERGED.fastq.gz        90 15.230481 10.761130
#> ADABM30X-H_S15_MERGED.fastq.gz       132 17.895618 10.730501
#> ADABM30X-M_S16_MERGED.fastq.gz        63 11.336162  6.831980
#> B18-006-B_S19_MERGED.fastq.gz         30  4.579201  2.462089
#> BA17-050-B_S21_MERGED.fastq.gz       162 14.968563  7.461988
#> BB19-006-H_S22_MERGED.fastq.gz       131  8.238795  4.533925
#> BB6-019-B_S23_MERGED.fastq.gz         88  6.009367  3.163064
#> BB6-019-H_S24_MERGED.fastq.gz        119  7.119430  2.622819
#> BB6-019-M_S25_MERGED.fastq.gz         88 25.769213 15.059862
#> BE9-006-B_S27_MERGED.fastq.gz         17  5.795489  3.344332
#> BE9-006-H_S28_MERGED.fastq.gz         11  5.501697  3.711273
#> BE9-006-M_S29_MERGED.fastq.gz         22 11.412889  8.606762
#> BG7-010-B_S30_MERGED.fastq.gz         39 15.258149  6.577062
#> BG7-010-H_S31_MERGED.fastq.gz         86 11.098284  4.417107
#> BG7-010-M_S32_MERGED.fastq.gz         52 15.365859  6.605331
#> BJ17-007-M_S34_MERGED.fastq.gz       158 10.857994  5.787965
#> BL7-006-B_S36_MERGED.fastq.gz         63 18.708521  6.979186
#> BL7-006-H_S37_MERGED.fastq.gz         50 22.572352 10.972762
#> BL7-006-M_S38_MERGED.fastq.gz         36 16.925041  9.666779
#> BP11-001-B_S43_MERGED.fastq.gz        84  6.626587  3.173056
#> BP11-001-H_S44_MERGED.fastq.gz       119  8.063167  2.946433
#> BP11-001-M_S45_MERGED.fastq.gz        73  3.914861  1.957519
#> BP12-025-B_S46_MERGED.fastq.gz        56  7.776651  4.052291
#> BQ4-018-B_S49_MERGED.fastq.gz        163  7.997291  4.608135
#> BQ4-018-H_S50_MERGED.fastq.gz         72  9.181735  4.824880
#> BQ4-018-M_S51_MERGED.fastq.gz        156 16.097605  8.275389
#> BT-006-M_S55_MERGED.fastq.gz          61  9.672663  5.999040
#> BV11-002-B_S57_MERGED.fastq.gz        98  6.944072  4.873753
#> BV11-002-H_S58_MERGED.fastq.gz       179 13.800674  5.592337
#> BV11-002-M_S59_MERGED.fastq.gz        76 11.799149  7.391348
#> C21-NV1-B_S62_MERGED.fastq.gz         35 11.286637  5.932252
#> C21-NV1-H_S63_MERGED.fastq.gz         73  6.527443  2.620203
#> C21-NV1-M_S64_MERGED.fastq.gz          8  2.866495  1.772771
#> CB8-019-B_S69_MERGED.fastq.gz         16 11.706178  8.812950
#> CB8-019-H_S70_MERGED.fastq.gz         33  1.397644  1.109174
#> CB8-019-M_S71_MERGED.fastq.gz         20  3.894227  2.051775
#> D18-003-B_S78_MERGED.fastq.gz         66  9.631651  5.428374
#> D18-003-H_S79_MERGED.fastq.gz         87  2.659561  1.478776
#> D18-003-M_S80_MERGED.fastq.gz         88  9.906082  4.506529
#> D61-010-B_S82_MERGED.fastq.gz         21  1.524569  1.214513
#> D9-027-B_S83_MERGED.fastq.gz          99 19.433579  7.589108
#> D9-027-H_S84_MERGED.fastq.gz          84 12.646257  6.675085
#> D9-027-M_S85_MERGED.fastq.gz         111 24.089806  9.609966
#> DJ2-008-B_S87_MERGED.fastq.gz         29  7.607661  4.461234
#> DJ2-008-H_S88_MERGED.fastq.gz         34 16.691044 11.587030
#> DJ2-008-M_S89_MERGED.fastq.gz         50  5.454623  2.292853
#> DS1-ABM002-B_S91_MERGED.fastq.gz      77  7.417789  3.083001
#> DS1-ABM002-H_S92_MERGED.fastq.gz      43  3.466204  2.269978
#> DS1-ABM002-M_S93_MERGED.fastq.gz     132 18.528548 10.705628
#> DU3-045-B_S94_MERGED.fastq.gz        134 15.427321  7.151758
#> DY5-004-B_S96_MERGED.fastq.gz         61  7.040679  3.558646
#> DY5-004-H_S97_MERGED.fastq.gz         27  2.638123  1.663783
#> DY5-004-M_S98_MERGED.fastq.gz         19  4.699036  2.188861
#> E9-009-B_S100_MERGED.fastq.gz         21 17.723121 15.144144
#> E9-009-H_S101_MERGED.fastq.gz         20 14.639441 12.096774
#> E9-009-M_S102_MERGED.fastq.gz         25  2.707789  1.585136
#> EC2-013-B_S104_MERGED.fastq.gz        89  2.935296  1.535002
#> F7-015-M_S106_MERGED.fastq.gz         52  6.853456  3.521697
#> FOMES19-H_S108_MERGED.fastq.gz        23  1.558617  1.197679
#> FOMES19-M_S109_MERGED.fastq.gz        21  1.480968  1.160913
#> H10-018-M_S110_MERGED.fastq.gz         4  1.016264  1.003876
#> H24-NVABM1-H_S111_MERGED.fastq.gz     62 14.250086  6.975467
#> J18-004-B_S114_MERGED.fastq.gz        21 13.693906 10.157100
#> J18-004-H_S115_MERGED.fastq.gz        47  5.330861  3.285425
#> J18-004-M_S116_MERGED.fastq.gz        34 14.097513  8.876414
#> K18-002-H_S117_MERGED.fastq.gz        74  1.635048  1.171398
#> L19X-B_S119_MERGED.fastq.gz           46 19.456407 12.677086
#> L19X-H_S120_MERGED.fastq.gz           93 13.777911  7.445490
#> L19X-M_S121_MERGED.fastq.gz           84  9.310051  5.264056
#> L23-002-B_S122_MERGED.fastq.gz       147 17.929015  6.462027
#> L23-002-H_S123_MERGED.fastq.gz        55 10.557775  6.035842
#> L23-002-M_S124_MERGED.fastq.gz        27  3.268428  2.469693
#> N19X-B_S126_MERGED.fastq.gz           34 14.626454  7.827476
#> N19X-H_S127_MERGED.fastq.gz          194 18.968399  7.731261
#> N19X-M_S128_MERGED.fastq.gz          113 22.360256 11.184175
#> N22-001-B_S129_MERGED.fastq.gz         4  2.790604  2.314286
#> N23-002-B_S130_MERGED.fastq.gz       135  1.667228  1.165840
#> N23-002-H_S131_MERGED.fastq.gz        74 13.579320  7.933936
#> N23-002-M_S132_MERGED.fastq.gz         3  2.021246  1.750455
#> NVABM-0163-H_S135_MERGED.fastq.gz     79  9.891670  3.980326
#>                                                                   X
#> A10-005-B_S188_MERGED.fastq.gz       A10-005-B_S188_MERGED.fastq.gz
#> A10-005-H_S189_MERGED.fastq.gz       A10-005-H_S189_MERGED.fastq.gz
#> A10-005-M_S190_MERGED.fastq.gz       A10-005-M_S190_MERGED.fastq.gz
#> A12-007-B_S2_MERGED.fastq.gz           A12-007-B_S2_MERGED.fastq.gz
#> AB29-ABMX-H_S6_MERGED.fastq.gz       AB29-ABMX-H_S6_MERGED.fastq.gz
#> AD26-005-B_S9_MERGED.fastq.gz         AD26-005-B_S9_MERGED.fastq.gz
#> AD26-005-H_S10_MERGED.fastq.gz       AD26-005-H_S10_MERGED.fastq.gz
#> AD26-005-M_S11_MERGED.fastq.gz       AD26-005-M_S11_MERGED.fastq.gz
#> AD30-ABMX-M_S12_MERGED.fastq.gz     AD30-ABMX-M_S12_MERGED.fastq.gz
#> AD32-007-M_S13_MERGED.fastq.gz       AD32-007-M_S13_MERGED.fastq.gz
#> ADABM30X-B_S14_MERGED.fastq.gz       ADABM30X-B_S14_MERGED.fastq.gz
#> ADABM30X-H_S15_MERGED.fastq.gz       ADABM30X-H_S15_MERGED.fastq.gz
#> ADABM30X-M_S16_MERGED.fastq.gz       ADABM30X-M_S16_MERGED.fastq.gz
#> B18-006-B_S19_MERGED.fastq.gz         B18-006-B_S19_MERGED.fastq.gz
#> BA17-050-B_S21_MERGED.fastq.gz       BA17-050-B_S21_MERGED.fastq.gz
#> BB19-006-H_S22_MERGED.fastq.gz       BB19-006-H_S22_MERGED.fastq.gz
#> BB6-019-B_S23_MERGED.fastq.gz         BB6-019-B_S23_MERGED.fastq.gz
#> BB6-019-H_S24_MERGED.fastq.gz         BB6-019-H_S24_MERGED.fastq.gz
#> BB6-019-M_S25_MERGED.fastq.gz         BB6-019-M_S25_MERGED.fastq.gz
#> BE9-006-B_S27_MERGED.fastq.gz         BE9-006-B_S27_MERGED.fastq.gz
#> BE9-006-H_S28_MERGED.fastq.gz         BE9-006-H_S28_MERGED.fastq.gz
#> BE9-006-M_S29_MERGED.fastq.gz         BE9-006-M_S29_MERGED.fastq.gz
#> BG7-010-B_S30_MERGED.fastq.gz         BG7-010-B_S30_MERGED.fastq.gz
#> BG7-010-H_S31_MERGED.fastq.gz         BG7-010-H_S31_MERGED.fastq.gz
#> BG7-010-M_S32_MERGED.fastq.gz         BG7-010-M_S32_MERGED.fastq.gz
#> BJ17-007-M_S34_MERGED.fastq.gz       BJ17-007-M_S34_MERGED.fastq.gz
#> BL7-006-B_S36_MERGED.fastq.gz         BL7-006-B_S36_MERGED.fastq.gz
#> BL7-006-H_S37_MERGED.fastq.gz         BL7-006-H_S37_MERGED.fastq.gz
#> BL7-006-M_S38_MERGED.fastq.gz         BL7-006-M_S38_MERGED.fastq.gz
#> BP11-001-B_S43_MERGED.fastq.gz       BP11-001-B_S43_MERGED.fastq.gz
#> BP11-001-H_S44_MERGED.fastq.gz       BP11-001-H_S44_MERGED.fastq.gz
#> BP11-001-M_S45_MERGED.fastq.gz       BP11-001-M_S45_MERGED.fastq.gz
#> BP12-025-B_S46_MERGED.fastq.gz       BP12-025-B_S46_MERGED.fastq.gz
#> BQ4-018-B_S49_MERGED.fastq.gz         BQ4-018-B_S49_MERGED.fastq.gz
#> BQ4-018-H_S50_MERGED.fastq.gz         BQ4-018-H_S50_MERGED.fastq.gz
#> BQ4-018-M_S51_MERGED.fastq.gz         BQ4-018-M_S51_MERGED.fastq.gz
#> BT-006-M_S55_MERGED.fastq.gz           BT-006-M_S55_MERGED.fastq.gz
#> BV11-002-B_S57_MERGED.fastq.gz       BV11-002-B_S57_MERGED.fastq.gz
#> BV11-002-H_S58_MERGED.fastq.gz       BV11-002-H_S58_MERGED.fastq.gz
#> BV11-002-M_S59_MERGED.fastq.gz       BV11-002-M_S59_MERGED.fastq.gz
#> C21-NV1-B_S62_MERGED.fastq.gz         C21-NV1-B_S62_MERGED.fastq.gz
#> C21-NV1-H_S63_MERGED.fastq.gz         C21-NV1-H_S63_MERGED.fastq.gz
#> C21-NV1-M_S64_MERGED.fastq.gz         C21-NV1-M_S64_MERGED.fastq.gz
#> CB8-019-B_S69_MERGED.fastq.gz         CB8-019-B_S69_MERGED.fastq.gz
#> CB8-019-H_S70_MERGED.fastq.gz         CB8-019-H_S70_MERGED.fastq.gz
#> CB8-019-M_S71_MERGED.fastq.gz         CB8-019-M_S71_MERGED.fastq.gz
#> D18-003-B_S78_MERGED.fastq.gz         D18-003-B_S78_MERGED.fastq.gz
#> D18-003-H_S79_MERGED.fastq.gz         D18-003-H_S79_MERGED.fastq.gz
#> D18-003-M_S80_MERGED.fastq.gz         D18-003-M_S80_MERGED.fastq.gz
#> D61-010-B_S82_MERGED.fastq.gz         D61-010-B_S82_MERGED.fastq.gz
#> D9-027-B_S83_MERGED.fastq.gz           D9-027-B_S83_MERGED.fastq.gz
#> D9-027-H_S84_MERGED.fastq.gz           D9-027-H_S84_MERGED.fastq.gz
#> D9-027-M_S85_MERGED.fastq.gz           D9-027-M_S85_MERGED.fastq.gz
#> DJ2-008-B_S87_MERGED.fastq.gz         DJ2-008-B_S87_MERGED.fastq.gz
#> DJ2-008-H_S88_MERGED.fastq.gz         DJ2-008-H_S88_MERGED.fastq.gz
#> DJ2-008-M_S89_MERGED.fastq.gz         DJ2-008-M_S89_MERGED.fastq.gz
#> DS1-ABM002-B_S91_MERGED.fastq.gz   DS1-ABM002-B_S91_MERGED.fastq.gz
#> DS1-ABM002-H_S92_MERGED.fastq.gz   DS1-ABM002-H_S92_MERGED.fastq.gz
#> DS1-ABM002-M_S93_MERGED.fastq.gz   DS1-ABM002-M_S93_MERGED.fastq.gz
#> DU3-045-B_S94_MERGED.fastq.gz         DU3-045-B_S94_MERGED.fastq.gz
#> DY5-004-B_S96_MERGED.fastq.gz         DY5-004-B_S96_MERGED.fastq.gz
#> DY5-004-H_S97_MERGED.fastq.gz         DY5-004-H_S97_MERGED.fastq.gz
#> DY5-004-M_S98_MERGED.fastq.gz         DY5-004-M_S98_MERGED.fastq.gz
#> E9-009-B_S100_MERGED.fastq.gz         E9-009-B_S100_MERGED.fastq.gz
#> E9-009-H_S101_MERGED.fastq.gz         E9-009-H_S101_MERGED.fastq.gz
#> E9-009-M_S102_MERGED.fastq.gz         E9-009-M_S102_MERGED.fastq.gz
#> EC2-013-B_S104_MERGED.fastq.gz       EC2-013-B_S104_MERGED.fastq.gz
#> F7-015-M_S106_MERGED.fastq.gz         F7-015-M_S106_MERGED.fastq.gz
#> FOMES19-H_S108_MERGED.fastq.gz       FOMES19-H_S108_MERGED.fastq.gz
#> FOMES19-M_S109_MERGED.fastq.gz       FOMES19-M_S109_MERGED.fastq.gz
#> H10-018-M_S110_MERGED.fastq.gz       H10-018-M_S110_MERGED.fastq.gz
#> H24-NVABM1-H_S111_MERGED.fastq.gz H24-NVABM1-H_S111_MERGED.fastq.gz
#> J18-004-B_S114_MERGED.fastq.gz       J18-004-B_S114_MERGED.fastq.gz
#> J18-004-H_S115_MERGED.fastq.gz       J18-004-H_S115_MERGED.fastq.gz
#> J18-004-M_S116_MERGED.fastq.gz       J18-004-M_S116_MERGED.fastq.gz
#> K18-002-H_S117_MERGED.fastq.gz       K18-002-H_S117_MERGED.fastq.gz
#> L19X-B_S119_MERGED.fastq.gz             L19X-B_S119_MERGED.fastq.gz
#> L19X-H_S120_MERGED.fastq.gz             L19X-H_S120_MERGED.fastq.gz
#> L19X-M_S121_MERGED.fastq.gz             L19X-M_S121_MERGED.fastq.gz
#> L23-002-B_S122_MERGED.fastq.gz       L23-002-B_S122_MERGED.fastq.gz
#> L23-002-H_S123_MERGED.fastq.gz       L23-002-H_S123_MERGED.fastq.gz
#> L23-002-M_S124_MERGED.fastq.gz       L23-002-M_S124_MERGED.fastq.gz
#> N19X-B_S126_MERGED.fastq.gz             N19X-B_S126_MERGED.fastq.gz
#> N19X-H_S127_MERGED.fastq.gz             N19X-H_S127_MERGED.fastq.gz
#> N19X-M_S128_MERGED.fastq.gz             N19X-M_S128_MERGED.fastq.gz
#> N22-001-B_S129_MERGED.fastq.gz       N22-001-B_S129_MERGED.fastq.gz
#> N23-002-B_S130_MERGED.fastq.gz       N23-002-B_S130_MERGED.fastq.gz
#> N23-002-H_S131_MERGED.fastq.gz       N23-002-H_S131_MERGED.fastq.gz
#> N23-002-M_S132_MERGED.fastq.gz       N23-002-M_S132_MERGED.fastq.gz
#> NVABM-0163-H_S135_MERGED.fastq.gz NVABM-0163-H_S135_MERGED.fastq.gz
#>                                        Sample_names   Tree_name Sample_id
#> A10-005-B_S188_MERGED.fastq.gz       A10-005-B_S188     A10-005       188
#> A10-005-H_S189_MERGED.fastq.gz       A10-005-H_S189     A10-005       189
#> A10-005-M_S190_MERGED.fastq.gz       A10-005-M_S190     A10-005       190
#> A12-007-B_S2_MERGED.fastq.gz           A12-007-B_S2     A12-007         2
#> AB29-ABMX-H_S6_MERGED.fastq.gz       AB29-ABMX-H_S6  AB29-abm-X         6
#> AD26-005-B_S9_MERGED.fastq.gz         AD26-005-B_S9    AD26-005         9
#> AD26-005-H_S10_MERGED.fastq.gz       AD26-005-H_S10    AD26-005        10
#> AD26-005-M_S11_MERGED.fastq.gz       AD26-005-M_S11    AD26-005        11
#> AD30-ABMX-M_S12_MERGED.fastq.gz     AD30-ABMX-M_S12  AD30-abm-X        12
#> AD32-007-M_S13_MERGED.fastq.gz       AD32-007-M_S13    AD32-007        13
#> ADABM30X-B_S14_MERGED.fastq.gz       ADABM30X-B_S14  AD30-abm-X        14
#> ADABM30X-H_S15_MERGED.fastq.gz       ADABM30X-H_S15  AD30-abm-X        15
#> ADABM30X-M_S16_MERGED.fastq.gz       ADABM30X-M_S16  AD30-abm-X        16
#> B18-006-B_S19_MERGED.fastq.gz         B18-006-B_S19     B18-006        19
#> BA17-050-B_S21_MERGED.fastq.gz       BA17-050-B_S21    BA17-050        21
#> BB19-006-H_S22_MERGED.fastq.gz       BB19-006-H_S22    BB19-006        22
#> BB6-019-B_S23_MERGED.fastq.gz         BB6-019-B_S23     BB6-019        23
#> BB6-019-H_S24_MERGED.fastq.gz         BB6-019-H_S24     BB6-019        24
#> BB6-019-M_S25_MERGED.fastq.gz         BB6-019-M_S25     BB6-019        25
#> BE9-006-B_S27_MERGED.fastq.gz         BE9-006-B_S27     BE9-006        27
#> BE9-006-H_S28_MERGED.fastq.gz         BE9-006-H_S28     BE9-006        28
#> BE9-006-M_S29_MERGED.fastq.gz         BE9-006-M_S29     BE9-006        29
#> BG7-010-B_S30_MERGED.fastq.gz         BG7-010-B_S30     BG7-010        30
#> BG7-010-H_S31_MERGED.fastq.gz         BG7-010-H_S31     BG7-010        31
#> BG7-010-M_S32_MERGED.fastq.gz         BG7-010-M_S32     BG7-010        32
#> BJ17-007-M_S34_MERGED.fastq.gz       BJ17-007-M_S34    BJ17-007        34
#> BL7-006-B_S36_MERGED.fastq.gz         BL7-006-B_S36     BL7-006        36
#> BL7-006-H_S37_MERGED.fastq.gz         BL7-006-H_S37     BL7-006        37
#> BL7-006-M_S38_MERGED.fastq.gz         BL7-006-M_S38     BL7-006        38
#> BP11-001-B_S43_MERGED.fastq.gz       BP11-001-B_S43    BP11-001        43
#> BP11-001-H_S44_MERGED.fastq.gz       BP11-001-H_S44    BP11-001        44
#> BP11-001-M_S45_MERGED.fastq.gz       BP11-001-M_S45    BP11-001        45
#> BP12-025-B_S46_MERGED.fastq.gz       BP12-025-B_S46    BP12-025        46
#> BQ4-018-B_S49_MERGED.fastq.gz         BQ4-018-B_S49     BQ4-018        49
#> BQ4-018-H_S50_MERGED.fastq.gz         BQ4-018-H_S50     BQ4-018        50
#> BQ4-018-M_S51_MERGED.fastq.gz         BQ4-018-M_S51     BQ4-018        51
#> BT-006-M_S55_MERGED.fastq.gz           BT-006-M_S55     BT7-006        55
#> BV11-002-B_S57_MERGED.fastq.gz       BV11-002-B_S57    BU11-002        57
#> BV11-002-H_S58_MERGED.fastq.gz       BV11-002-H_S58    BU11-002        58
#> BV11-002-M_S59_MERGED.fastq.gz       BV11-002-M_S59    BU11-002        59
#> C21-NV1-B_S62_MERGED.fastq.gz         C21-NV1-B_S62    C21-nv-1        62
#> C21-NV1-H_S63_MERGED.fastq.gz         C21-NV1-H_S63    C21-nv-1        63
#> C21-NV1-M_S64_MERGED.fastq.gz         C21-NV1-M_S64    C21-nv-1        64
#> CB8-019-B_S69_MERGED.fastq.gz         CB8-019-B_S69     CB8-019        69
#> CB8-019-H_S70_MERGED.fastq.gz         CB8-019-H_S70     CB8-019        70
#> CB8-019-M_S71_MERGED.fastq.gz         CB8-019-M_S71     CB8-019        71
#> D18-003-B_S78_MERGED.fastq.gz         D18-003-B_S78     D18-003        78
#> D18-003-H_S79_MERGED.fastq.gz         D18-003-H_S79     D18-003        79
#> D18-003-M_S80_MERGED.fastq.gz         D18-003-M_S80     D18-003        80
#> D61-010-B_S82_MERGED.fastq.gz         D61-010-B_S82     DG1-010        82
#> D9-027-B_S83_MERGED.fastq.gz           D9-027-B_S83      D9-027        83
#> D9-027-H_S84_MERGED.fastq.gz           D9-027-H_S84      D9-027        84
#> D9-027-M_S85_MERGED.fastq.gz           D9-027-M_S85      D9-027        85
#> DJ2-008-B_S87_MERGED.fastq.gz         DJ2-008-B_S87     DJ2-008        87
#> DJ2-008-H_S88_MERGED.fastq.gz         DJ2-008-H_S88     DJ2-008        88
#> DJ2-008-M_S89_MERGED.fastq.gz         DJ2-008-M_S89     DJ2-008        89
#> DS1-ABM002-B_S91_MERGED.fastq.gz   DS1-ABM002-B_S91   DS1abm002        91
#> DS1-ABM002-H_S92_MERGED.fastq.gz   DS1-ABM002-H_S92   DS1abm002        92
#> DS1-ABM002-M_S93_MERGED.fastq.gz   DS1-ABM002-M_S93   DS1abm002        93
#> DU3-045-B_S94_MERGED.fastq.gz         DU3-045-B_S94     DU3-045        94
#> DY5-004-B_S96_MERGED.fastq.gz         DY5-004-B_S96     DY5-004        96
#> DY5-004-H_S97_MERGED.fastq.gz         DY5-004-H_S97     DY5-004        97
#> DY5-004-M_S98_MERGED.fastq.gz         DY5-004-M_S98     DY5-004        98
#> E9-009-B_S100_MERGED.fastq.gz         E9-009-B_S100      E9-009       100
#> E9-009-H_S101_MERGED.fastq.gz         E9-009-H_S101      E9-009       101
#> E9-009-M_S102_MERGED.fastq.gz         E9-009-M_S102      E9-009       102
#> EC2-013-B_S104_MERGED.fastq.gz       EC2-013-B_S104     EC2-013       104
#> F7-015-M_S106_MERGED.fastq.gz         F7-015-M_S106      F7-015       106
#> FOMES19-H_S108_MERGED.fastq.gz       FOMES19-H_S108     FOMES19       108
#> FOMES19-M_S109_MERGED.fastq.gz       FOMES19-M_S109     FOMES19       109
#> H10-018-M_S110_MERGED.fastq.gz       H10-018-M_S110     H10-018       110
#> H24-NVABM1-H_S111_MERGED.fastq.gz H24-NVABM1-H_S111 H24-nvabm-1       111
#> J18-004-B_S114_MERGED.fastq.gz       J18-004-B_S114     J18-004       114
#> J18-004-H_S115_MERGED.fastq.gz       J18-004-H_S115     J18-004       115
#> J18-004-M_S116_MERGED.fastq.gz       J18-004-M_S116     J18-004       116
#> K18-002-H_S117_MERGED.fastq.gz       K18-002-H_S117     K18-002       117
#> L19X-B_S119_MERGED.fastq.gz             L19X-B_S119       L19-X       119
#> L19X-H_S120_MERGED.fastq.gz             L19X-H_S120       L19-X       120
#> L19X-M_S121_MERGED.fastq.gz             L19X-M_S121       L19-X       121
#> L23-002-B_S122_MERGED.fastq.gz       L23-002-B_S122     L23-002       122
#> L23-002-H_S123_MERGED.fastq.gz       L23-002-H_S123     L23-002       123
#> L23-002-M_S124_MERGED.fastq.gz       L23-002-M_S124     L23-002       124
#> N19X-B_S126_MERGED.fastq.gz             N19X-B_S126       N19-X       126
#> N19X-H_S127_MERGED.fastq.gz             N19X-H_S127       N19-X       127
#> N19X-M_S128_MERGED.fastq.gz             N19X-M_S128       N19-X       128
#> N22-001-B_S129_MERGED.fastq.gz       N22-001-B_S129     M22-001       129
#> N23-002-B_S130_MERGED.fastq.gz       N23-002-B_S130     N23-002       130
#> N23-002-H_S131_MERGED.fastq.gz       N23-002-H_S131     N23-002       131
#> N23-002-M_S132_MERGED.fastq.gz       N23-002-M_S132     N23-002       132
#> NVABM-0163-H_S135_MERGED.fastq.gz NVABM-0163-H_S135   nvabm0163       135
#>                                   Height Diameter Time
#> A10-005-B_S188_MERGED.fastq.gz       Low       52   15
#> A10-005-H_S189_MERGED.fastq.gz      High       52   15
#> A10-005-M_S190_MERGED.fastq.gz    Middle       52   15
#> A12-007-B_S2_MERGED.fastq.gz         Low     28,4    0
#> AB29-ABMX-H_S6_MERGED.fastq.gz      High       99    5
#> AD26-005-B_S9_MERGED.fastq.gz        Low    115,5   15
#> AD26-005-H_S10_MERGED.fastq.gz      High    115,5   15
#> AD26-005-M_S11_MERGED.fastq.gz    Middle    115,5   15
#> AD30-ABMX-M_S12_MERGED.fastq.gz   Middle        -    5
#> AD32-007-M_S13_MERGED.fastq.gz    Middle       52   NA
#> ADABM30X-B_S14_MERGED.fastq.gz       Low        -    5
#> ADABM30X-H_S15_MERGED.fastq.gz      High        -    5
#> ADABM30X-M_S16_MERGED.fastq.gz    Middle        -    5
#> B18-006-B_S19_MERGED.fastq.gz        Low     32,5   10
#> BA17-050-B_S21_MERGED.fastq.gz       Low       30   10
#> BB19-006-H_S22_MERGED.fastq.gz      High     37,5   15
#> BB6-019-B_S23_MERGED.fastq.gz        Low     33,8    5
#> BB6-019-H_S24_MERGED.fastq.gz       High     33,8    5
#> BB6-019-M_S25_MERGED.fastq.gz     Middle     33,8    5
#> BE9-006-B_S27_MERGED.fastq.gz        Low     23,1    0
#> BE9-006-H_S28_MERGED.fastq.gz       High     23,1    0
#> BE9-006-M_S29_MERGED.fastq.gz     Middle     23,1    0
#> BG7-010-B_S30_MERGED.fastq.gz        Low       32    0
#> BG7-010-H_S31_MERGED.fastq.gz       High       32    0
#> BG7-010-M_S32_MERGED.fastq.gz     Middle       32    0
#> BJ17-007-M_S34_MERGED.fastq.gz    Middle       32    5
#> BL7-006-B_S36_MERGED.fastq.gz        Low     35,6    0
#> BL7-006-H_S37_MERGED.fastq.gz       High     35,6    0
#> BL7-006-M_S38_MERGED.fastq.gz     Middle     35,6    0
#> BP11-001-B_S43_MERGED.fastq.gz       Low     47,5   15
#> BP11-001-H_S44_MERGED.fastq.gz      High     47,5   15
#> BP11-001-M_S45_MERGED.fastq.gz    Middle     47,5   15
#> BP12-025-B_S46_MERGED.fastq.gz       Low     64,5   10
#> BQ4-018-B_S49_MERGED.fastq.gz        Low       34    0
#> BQ4-018-H_S50_MERGED.fastq.gz       High       34    0
#> BQ4-018-M_S51_MERGED.fastq.gz     Middle       34    0
#> BT-006-M_S55_MERGED.fastq.gz      Middle       41    0
#> BV11-002-B_S57_MERGED.fastq.gz       Low       33    5
#> BV11-002-H_S58_MERGED.fastq.gz      High       33    5
#> BV11-002-M_S59_MERGED.fastq.gz    Middle       33    5
#> C21-NV1-B_S62_MERGED.fastq.gz        Low       30    0
#> C21-NV1-H_S63_MERGED.fastq.gz       High       30    0
#> C21-NV1-M_S64_MERGED.fastq.gz     Middle       30    0
#> CB8-019-B_S69_MERGED.fastq.gz        Low     33,3    0
#> CB8-019-H_S70_MERGED.fastq.gz       High     33,3    0
#> CB8-019-M_S71_MERGED.fastq.gz     Middle     33,3    0
#> D18-003-B_S78_MERGED.fastq.gz        Low     45,5    5
#> D18-003-H_S79_MERGED.fastq.gz       High     45,5    5
#> D18-003-M_S80_MERGED.fastq.gz     Middle     45,5    5
#> D61-010-B_S82_MERGED.fastq.gz        Low     36,3   15
#> D9-027-B_S83_MERGED.fastq.gz         Low     18,9   NA
#> D9-027-H_S84_MERGED.fastq.gz        High     18,9   NA
#> D9-027-M_S85_MERGED.fastq.gz      Middle     18,9   NA
#> DJ2-008-B_S87_MERGED.fastq.gz        Low     76,4    0
#> DJ2-008-H_S88_MERGED.fastq.gz       High     76,4    0
#> DJ2-008-M_S89_MERGED.fastq.gz     Middle     76,4    0
#> DS1-ABM002-B_S91_MERGED.fastq.gz     Low       25   15
#> DS1-ABM002-H_S92_MERGED.fastq.gz    High       25   15
#> DS1-ABM002-M_S93_MERGED.fastq.gz  Middle       25   15
#> DU3-045-B_S94_MERGED.fastq.gz        Low       38   15
#> DY5-004-B_S96_MERGED.fastq.gz        Low     60,4    0
#> DY5-004-H_S97_MERGED.fastq.gz       High     60,4    0
#> DY5-004-M_S98_MERGED.fastq.gz     Middle     60,4    0
#> E9-009-B_S100_MERGED.fastq.gz        Low     29,5    0
#> E9-009-H_S101_MERGED.fastq.gz       High     29,5    0
#> E9-009-M_S102_MERGED.fastq.gz     Middle     29,5    0
#> EC2-013-B_S104_MERGED.fastq.gz       Low       52   15
#> F7-015-M_S106_MERGED.fastq.gz     Middle       50   10
#> FOMES19-H_S108_MERGED.fastq.gz      High     <NA>   NA
#> FOMES19-M_S109_MERGED.fastq.gz    Middle     <NA>   NA
#> H10-018-M_S110_MERGED.fastq.gz    Middle     88,5   15
#> H24-NVABM1-H_S111_MERGED.fastq.gz   High       30    5
#> J18-004-B_S114_MERGED.fastq.gz       Low     27,6    0
#> J18-004-H_S115_MERGED.fastq.gz      High     27,6    0
#> J18-004-M_S116_MERGED.fastq.gz    Middle     27,6    0
#> K18-002-H_S117_MERGED.fastq.gz      High     33,5   10
#> L19X-B_S119_MERGED.fastq.gz          Low       30   NA
#> L19X-H_S120_MERGED.fastq.gz         High       30   NA
#> L19X-M_S121_MERGED.fastq.gz       Middle       30   NA
#> L23-002-B_S122_MERGED.fastq.gz       Low     40,5   10
#> L23-002-H_S123_MERGED.fastq.gz      High     40,5   10
#> L23-002-M_S124_MERGED.fastq.gz    Middle     40,5   10
#> N19X-B_S126_MERGED.fastq.gz          Low     33,5   NA
#> N19X-H_S127_MERGED.fastq.gz         High     33,5   NA
#> N19X-M_S128_MERGED.fastq.gz       Middle     33,5   NA
#> N22-001-B_S129_MERGED.fastq.gz       Low       34    0
#> N23-002-B_S130_MERGED.fastq.gz       Low     21,5   15
#> N23-002-H_S131_MERGED.fastq.gz      High     21,5   15
#> N23-002-M_S132_MERGED.fastq.gz    Middle     21,5   15
#> NVABM-0163-H_S135_MERGED.fastq.gz   High       40   10
#>                                                          .sample_id
#> A10-005-B_S188_MERGED.fastq.gz       A10-005-B_S188_MERGED.fastq.gz
#> A10-005-H_S189_MERGED.fastq.gz       A10-005-H_S189_MERGED.fastq.gz
#> A10-005-M_S190_MERGED.fastq.gz       A10-005-M_S190_MERGED.fastq.gz
#> A12-007-B_S2_MERGED.fastq.gz           A12-007-B_S2_MERGED.fastq.gz
#> AB29-ABMX-H_S6_MERGED.fastq.gz       AB29-ABMX-H_S6_MERGED.fastq.gz
#> AD26-005-B_S9_MERGED.fastq.gz         AD26-005-B_S9_MERGED.fastq.gz
#> AD26-005-H_S10_MERGED.fastq.gz       AD26-005-H_S10_MERGED.fastq.gz
#> AD26-005-M_S11_MERGED.fastq.gz       AD26-005-M_S11_MERGED.fastq.gz
#> AD30-ABMX-M_S12_MERGED.fastq.gz     AD30-ABMX-M_S12_MERGED.fastq.gz
#> AD32-007-M_S13_MERGED.fastq.gz       AD32-007-M_S13_MERGED.fastq.gz
#> ADABM30X-B_S14_MERGED.fastq.gz       ADABM30X-B_S14_MERGED.fastq.gz
#> ADABM30X-H_S15_MERGED.fastq.gz       ADABM30X-H_S15_MERGED.fastq.gz
#> ADABM30X-M_S16_MERGED.fastq.gz       ADABM30X-M_S16_MERGED.fastq.gz
#> B18-006-B_S19_MERGED.fastq.gz         B18-006-B_S19_MERGED.fastq.gz
#> BA17-050-B_S21_MERGED.fastq.gz       BA17-050-B_S21_MERGED.fastq.gz
#> BB19-006-H_S22_MERGED.fastq.gz       BB19-006-H_S22_MERGED.fastq.gz
#> BB6-019-B_S23_MERGED.fastq.gz         BB6-019-B_S23_MERGED.fastq.gz
#> BB6-019-H_S24_MERGED.fastq.gz         BB6-019-H_S24_MERGED.fastq.gz
#> BB6-019-M_S25_MERGED.fastq.gz         BB6-019-M_S25_MERGED.fastq.gz
#> BE9-006-B_S27_MERGED.fastq.gz         BE9-006-B_S27_MERGED.fastq.gz
#> BE9-006-H_S28_MERGED.fastq.gz         BE9-006-H_S28_MERGED.fastq.gz
#> BE9-006-M_S29_MERGED.fastq.gz         BE9-006-M_S29_MERGED.fastq.gz
#> BG7-010-B_S30_MERGED.fastq.gz         BG7-010-B_S30_MERGED.fastq.gz
#> BG7-010-H_S31_MERGED.fastq.gz         BG7-010-H_S31_MERGED.fastq.gz
#> BG7-010-M_S32_MERGED.fastq.gz         BG7-010-M_S32_MERGED.fastq.gz
#> BJ17-007-M_S34_MERGED.fastq.gz       BJ17-007-M_S34_MERGED.fastq.gz
#> BL7-006-B_S36_MERGED.fastq.gz         BL7-006-B_S36_MERGED.fastq.gz
#> BL7-006-H_S37_MERGED.fastq.gz         BL7-006-H_S37_MERGED.fastq.gz
#> BL7-006-M_S38_MERGED.fastq.gz         BL7-006-M_S38_MERGED.fastq.gz
#> BP11-001-B_S43_MERGED.fastq.gz       BP11-001-B_S43_MERGED.fastq.gz
#> BP11-001-H_S44_MERGED.fastq.gz       BP11-001-H_S44_MERGED.fastq.gz
#> BP11-001-M_S45_MERGED.fastq.gz       BP11-001-M_S45_MERGED.fastq.gz
#> BP12-025-B_S46_MERGED.fastq.gz       BP12-025-B_S46_MERGED.fastq.gz
#> BQ4-018-B_S49_MERGED.fastq.gz         BQ4-018-B_S49_MERGED.fastq.gz
#> BQ4-018-H_S50_MERGED.fastq.gz         BQ4-018-H_S50_MERGED.fastq.gz
#> BQ4-018-M_S51_MERGED.fastq.gz         BQ4-018-M_S51_MERGED.fastq.gz
#> BT-006-M_S55_MERGED.fastq.gz           BT-006-M_S55_MERGED.fastq.gz
#> BV11-002-B_S57_MERGED.fastq.gz       BV11-002-B_S57_MERGED.fastq.gz
#> BV11-002-H_S58_MERGED.fastq.gz       BV11-002-H_S58_MERGED.fastq.gz
#> BV11-002-M_S59_MERGED.fastq.gz       BV11-002-M_S59_MERGED.fastq.gz
#> C21-NV1-B_S62_MERGED.fastq.gz         C21-NV1-B_S62_MERGED.fastq.gz
#> C21-NV1-H_S63_MERGED.fastq.gz         C21-NV1-H_S63_MERGED.fastq.gz
#> C21-NV1-M_S64_MERGED.fastq.gz         C21-NV1-M_S64_MERGED.fastq.gz
#> CB8-019-B_S69_MERGED.fastq.gz         CB8-019-B_S69_MERGED.fastq.gz
#> CB8-019-H_S70_MERGED.fastq.gz         CB8-019-H_S70_MERGED.fastq.gz
#> CB8-019-M_S71_MERGED.fastq.gz         CB8-019-M_S71_MERGED.fastq.gz
#> D18-003-B_S78_MERGED.fastq.gz         D18-003-B_S78_MERGED.fastq.gz
#> D18-003-H_S79_MERGED.fastq.gz         D18-003-H_S79_MERGED.fastq.gz
#> D18-003-M_S80_MERGED.fastq.gz         D18-003-M_S80_MERGED.fastq.gz
#> D61-010-B_S82_MERGED.fastq.gz         D61-010-B_S82_MERGED.fastq.gz
#> D9-027-B_S83_MERGED.fastq.gz           D9-027-B_S83_MERGED.fastq.gz
#> D9-027-H_S84_MERGED.fastq.gz           D9-027-H_S84_MERGED.fastq.gz
#> D9-027-M_S85_MERGED.fastq.gz           D9-027-M_S85_MERGED.fastq.gz
#> DJ2-008-B_S87_MERGED.fastq.gz         DJ2-008-B_S87_MERGED.fastq.gz
#> DJ2-008-H_S88_MERGED.fastq.gz         DJ2-008-H_S88_MERGED.fastq.gz
#> DJ2-008-M_S89_MERGED.fastq.gz         DJ2-008-M_S89_MERGED.fastq.gz
#> DS1-ABM002-B_S91_MERGED.fastq.gz   DS1-ABM002-B_S91_MERGED.fastq.gz
#> DS1-ABM002-H_S92_MERGED.fastq.gz   DS1-ABM002-H_S92_MERGED.fastq.gz
#> DS1-ABM002-M_S93_MERGED.fastq.gz   DS1-ABM002-M_S93_MERGED.fastq.gz
#> DU3-045-B_S94_MERGED.fastq.gz         DU3-045-B_S94_MERGED.fastq.gz
#> DY5-004-B_S96_MERGED.fastq.gz         DY5-004-B_S96_MERGED.fastq.gz
#> DY5-004-H_S97_MERGED.fastq.gz         DY5-004-H_S97_MERGED.fastq.gz
#> DY5-004-M_S98_MERGED.fastq.gz         DY5-004-M_S98_MERGED.fastq.gz
#> E9-009-B_S100_MERGED.fastq.gz         E9-009-B_S100_MERGED.fastq.gz
#> E9-009-H_S101_MERGED.fastq.gz         E9-009-H_S101_MERGED.fastq.gz
#> E9-009-M_S102_MERGED.fastq.gz         E9-009-M_S102_MERGED.fastq.gz
#> EC2-013-B_S104_MERGED.fastq.gz       EC2-013-B_S104_MERGED.fastq.gz
#> F7-015-M_S106_MERGED.fastq.gz         F7-015-M_S106_MERGED.fastq.gz
#> FOMES19-H_S108_MERGED.fastq.gz       FOMES19-H_S108_MERGED.fastq.gz
#> FOMES19-M_S109_MERGED.fastq.gz       FOMES19-M_S109_MERGED.fastq.gz
#> H10-018-M_S110_MERGED.fastq.gz       H10-018-M_S110_MERGED.fastq.gz
#> H24-NVABM1-H_S111_MERGED.fastq.gz H24-NVABM1-H_S111_MERGED.fastq.gz
#> J18-004-B_S114_MERGED.fastq.gz       J18-004-B_S114_MERGED.fastq.gz
#> J18-004-H_S115_MERGED.fastq.gz       J18-004-H_S115_MERGED.fastq.gz
#> J18-004-M_S116_MERGED.fastq.gz       J18-004-M_S116_MERGED.fastq.gz
#> K18-002-H_S117_MERGED.fastq.gz       K18-002-H_S117_MERGED.fastq.gz
#> L19X-B_S119_MERGED.fastq.gz             L19X-B_S119_MERGED.fastq.gz
#> L19X-H_S120_MERGED.fastq.gz             L19X-H_S120_MERGED.fastq.gz
#> L19X-M_S121_MERGED.fastq.gz             L19X-M_S121_MERGED.fastq.gz
#> L23-002-B_S122_MERGED.fastq.gz       L23-002-B_S122_MERGED.fastq.gz
#> L23-002-H_S123_MERGED.fastq.gz       L23-002-H_S123_MERGED.fastq.gz
#> L23-002-M_S124_MERGED.fastq.gz       L23-002-M_S124_MERGED.fastq.gz
#> N19X-B_S126_MERGED.fastq.gz             N19X-B_S126_MERGED.fastq.gz
#> N19X-H_S127_MERGED.fastq.gz             N19X-H_S127_MERGED.fastq.gz
#> N19X-M_S128_MERGED.fastq.gz             N19X-M_S128_MERGED.fastq.gz
#> N22-001-B_S129_MERGED.fastq.gz       N22-001-B_S129_MERGED.fastq.gz
#> N23-002-B_S130_MERGED.fastq.gz       N23-002-B_S130_MERGED.fastq.gz
#> N23-002-H_S131_MERGED.fastq.gz       N23-002-H_S131_MERGED.fastq.gz
#> N23-002-M_S132_MERGED.fastq.gz       N23-002-M_S132_MERGED.fastq.gz
#> NVABM-0163-H_S135_MERGED.fastq.gz NVABM-0163-H_S135_MERGED.fastq.gz
#>  [ reached 'max' / getOption("max.print") -- omitted 41 rows ]
#> 
#> $dabest_objects
#> $dabest_objects$Hill_0
#> DABESTR v2025.3.15
#> ==================
#> 
#> Good afternoon!
#> The current time is 14:36 PM on Friday March 13, 2026.
#> 
#> The character(0) Cohen's d between Low and High is -0.12 [95%CI -0.56, 0.309].
#> The p-value of the two-sided permutation t-test is 0.5813, calculated for legacy purposes only.
#> 
#> The character(0) Cohen's d between Middle and High is -0.187 [95%CI -0.617, 0.232].
#> The p-value of the two-sided permutation t-test is 0.3911, calculated for legacy purposes only.
#> 
#> 5000 bootstrap samples were taken; the confidence interval is bias-corrected and accelerated.
#> Any p-value reported is the probability of observing the effect size (or greater),
#> assuming the null hypothesis of zero difference is true.
#> For each p-value, 5000 reshuffles of the control and test labels were performed.
#> 
#> 
#> $dabest_objects$Hill_1
#> DABESTR v2025.3.15
#> ==================
#> 
#> Good afternoon!
#> The current time is 14:36 PM on Friday March 13, 2026.
#> 
#> The character(0) Cohen's d between Low and High is -0.016 [95%CI -0.434, 0.418].
#> The p-value of the two-sided permutation t-test is 0.9417, calculated for legacy purposes only.
#> 
#> The character(0) Cohen's d between Middle and High is 0.163 [95%CI -0.254, 0.585].
#> The p-value of the two-sided permutation t-test is 0.4475, calculated for legacy purposes only.
#> 
#> 5000 bootstrap samples were taken; the confidence interval is bias-corrected and accelerated.
#> Any p-value reported is the probability of observing the effect size (or greater),
#> assuming the null hypothesis of zero difference is true.
#> For each p-value, 5000 reshuffles of the control and test labels were performed.
#> 
#> 
#> $dabest_objects$Hill_2
#> DABESTR v2025.3.15
#> ==================
#> 
#> Good afternoon!
#> The current time is 14:36 PM on Friday March 13, 2026.
#> 
#> The character(0) Cohen's d between Low and High is -0.039 [95%CI -0.47, 0.387].
#> The p-value of the two-sided permutation t-test is 0.8581, calculated for legacy purposes only.
#> 
#> The character(0) Cohen's d between Middle and High is 0.187 [95%CI -0.23, 0.592].
#> The p-value of the two-sided permutation t-test is 0.3828, calculated for legacy purposes only.
#> 
#> 5000 bootstrap samples were taken; the confidence interval is bias-corrected and accelerated.
#> Any p-value reported is the probability of observing the effect size (or greater),
#> assuming the null hypothesis of zero difference is true.
#> For each p-value, 5000 reshuffles of the control and test labels were performed.
#> 
#> 
#> 
#> $plots
#> $plots$Hill_0

#> 
#> $plots$Hill_1

#> 
#> $plots$Hill_2

#> 
#> 
#> $summary
#> # A tibble: 6 × 8
#>   metric comparison   effect_size ci_lower ci_upper pvalue_permtest pvalue_welch
#>   <chr>  <chr>              <dbl>    <dbl>    <dbl>           <dbl>        <dbl>
#> 1 Hill_0 High vs Low      -0.120    -0.560    0.309           0.577        0.581
#> 2 Hill_0 High vs Mid…     -0.187    -0.617    0.232           0.386        0.391
#> 3 Hill_1 High vs Low      -0.0158   -0.434    0.418           0.942        0.942
#> 4 Hill_1 High vs Mid…      0.163    -0.254    0.585           0.454        0.448
#> 5 Hill_2 High vs Low      -0.0387   -0.470    0.387           0.869        0.858
#> 6 Hill_2 High vs Mid…      0.187    -0.230    0.592           0.396        0.383
#> # ℹ 1 more variable: pvalue_mann_whitney <dbl>
#> 
#> $effect_type
#> [1] "cohens_d"
#> 
#> attr(,"class")
#> [1] "estim_diff_pq_result"
res$plots$Hill_0

res$summary
#> # A tibble: 6 × 8
#>   metric comparison   effect_size ci_lower ci_upper pvalue_permtest pvalue_welch
#>   <chr>  <chr>              <dbl>    <dbl>    <dbl>           <dbl>        <dbl>
#> 1 Hill_0 High vs Low      -0.120    -0.560    0.309           0.577        0.581
#> 2 Hill_0 High vs Mid…     -0.187    -0.617    0.232           0.386        0.391
#> 3 Hill_1 High vs Low      -0.0158   -0.434    0.418           0.942        0.942
#> 4 Hill_1 High vs Mid…      0.163    -0.254    0.585           0.454        0.448
#> 5 Hill_2 High vs Low      -0.0387   -0.470    0.387           0.869        0.858
#> 6 Hill_2 High vs Mid…      0.187    -0.230    0.592           0.396        0.383
#> # ℹ 1 more variable: pvalue_mann_whitney <dbl>
```
