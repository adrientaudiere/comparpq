# Multiply OTU counts conditionally based on sample metadata

Multiplies OTU table counts by different values depending on the level
of a factor column in `sample_data`. Samples whose factor value does not
match any condition (including `NA` samples) are left unchanged unless
explicitly targeted.

## Usage

``` r
multiply_counts_pq(
  physeq,
  fact,
  conditions,
  multipliers,
  prop_taxa = 0.5,
  seed = NULL,
  compensate = FALSE,
  min_prevalence = 0,
  round = TRUE,
  verbose = TRUE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (character, required) The name of a column in `sample_data(physeq)`.

- conditions:

  (character vector, required) Values of `fact` to match. Use `NA` to
  target samples with missing values.

- multipliers:

  (numeric vector, required) Multiplying factors corresponding to each
  element of `conditions`. Must be the same length as `conditions`.

- prop_taxa:

  (numeric, default 0.5) Proportion of taxa (between 0 and 1) to which
  the multiplication is applied. Taxa are randomly sampled. Set to 1 to
  apply to all taxa.

- seed:

  (integer, default NULL) Random seed for reproducible taxa sampling
  when `prop_taxa < 1`.

- compensate:

  (logical, default FALSE) If TRUE, scale down the non-selected taxa so
  that total library size per sample is preserved. This creates a pure
  compositional shift: selected taxa gain relative abundance while the
  rest lose it. This mode is essential for testing differential
  abundance methods (ALDEx2, ANCOM-BC) which normalize for library size.

- min_prevalence:

  (numeric, default 0) Minimum prevalence (proportion of matched samples
  with non-zero counts) for a taxon to be eligible for selection.
  Setting this \> 0 (e.g. 0.5) ensures that only taxa actually present
  in the target samples are inflated, producing a stronger and more
  realistic DA signal.

- round:

  (logical, default TRUE) If TRUE, round the resulting counts to
  integers (since OTU tables typically contain integer counts).

- verbose:

  (logical, default TRUE) If TRUE, print a message with the number of
  modified taxa.

## Value

A phyloseq object with modified OTU counts. The selected taxa are stored
as an attribute accessible via `attr(result, "taxa_modified")`.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Each sample is checked against `conditions` in order. When the sample's
value in `fact` matches a condition, its counts are multiplied by the
corresponding multiplier. Samples that do not match any condition are
left unchanged.

When `prop_taxa < 1`, only a random subset of taxa is affected by the
multiplication. This is useful to simulate differential abundance where
only some taxa respond to a condition.

**Simulating differential abundance for DA methods:** DA methods like
ALDEx2 and ANCOM-BC work on compositions (relative abundances). Simply
multiplying counts changes library size, which these methods normalize
away. Use `compensate = TRUE` to create a real compositional shift:
selected taxa are inflated, then non- selected taxa are scaled down so
total counts per sample stay constant. Combine with `min_prevalence > 0`
and small `prop_taxa` (e.g. 0.05) for a realistic DA signal on few,
prevalent taxa.

**Compensation limits:** When the multiplier is very high and selected
taxa already have high counts, the non-selected taxa may not have enough
counts to compensate. In such cases, non-selected taxa are zeroed and
selected taxa are scaled down to preserve library size. A warning is
issued when this occurs. To avoid this, use moderate multipliers (2-3)
or select taxa with lower initial abundance.

To target `NA` values in the factor column, include `NA` in `conditions`
(e.g., `conditions = c("High", NA)`).

## Author

Adrien Taudière

## Examples

``` r
# Multiply counts by 1.2 for "High" samples in half of taxa
d <- multiply_counts_pq(data_fungi,
  fact = "Height",
  conditions = "High", multipliers = 2
)
#> Modified 710 taxa in 41 matched samples

# Multiple conditions with different multipliers
d2 <- multiply_counts_pq(data_fungi,
  fact = "Height",
  conditions = c("High", "Low"), multipliers = c(1.5, 0.8)
)
#> Modified 710 taxa in 86 matched samples

# Target NA samples
d3 <- multiply_counts_pq(data_fungi,
  fact = "Height",
  conditions = c("High", NA), multipliers = c(1.2, 0.5)
)
#> Modified 710 taxa in 95 matched samples

sum(data_fungi@otu_table)
#> [1] 1839124
sum(d@otu_table)
#> [1] 2116884
sum(d2@otu_table)
#> [1] 1911558
sum(d3@otu_table)
#> [1] 1723185

# Apply only to 30% of taxa with a reproducible seed set
d4 <- multiply_counts_pq(data_fungi,
  fact = "Height",
  conditions = "High", multipliers = 2, prop_taxa = 0.3, seed = 42
)
#> Modified 426 taxa in 41 matched samples

# Simulate DA for compositional methods (ALDEx2, ANCOM-BC):
# Use compensate=TRUE to create a real compositional shift
# Use min_prevalence to select only prevalent taxa
# Use small prop_taxa to affect few taxa strongly
d_da <- multiply_counts_pq(
  data_fungi,
  fact = "Height",
  conditions = "High",
  multipliers = 5,
  prop_taxa = 0.05,
  min_prevalence = 0.5,
  compensate = TRUE,
  seed = 123
)
#> Warning: 4 sample(s) could not be fully compensated (multiplier too high relative to non-selected taxa abundance). Non-selected taxa were zeroed and selected taxa were scaled down to preserve library size.
#> Modified 1 taxa in 41 matched samples

# Library sizes are preserved
sum(sample_sums(data_fungi))
#> [1] 1839124
sum(sample_sums(d_da))
#> [1] 1839184

# Test with maaslin3
res <- maaslin3_pq(d_da,
  formula = "~Height"
)
#> 2026-06-17 13:44:57.80 INFO::Writing function arguments to log file
#> 2026-06-17 13:44:57.82 INFO::Verifying options selected are valid
#> 2026-06-17 13:44:57.82 INFO::Determining format of input files
#> 2026-06-17 13:44:57.82 INFO::Input format is data samples as rows and metadata samples as rows
#> 2026-06-17 13:44:57.83 INFO::Running selected normalization method: TSS
#> 2026-06-17 13:44:57.87 INFO::Writing normalized data to file res_maaslin3/features/data_norm.tsv
#> 2026-06-17 13:44:57.91 INFO::Filter data based on min abundance, min prevalence, and max prevalence
#> 2026-06-17 13:44:57.91 INFO::Total samples in data: 185
#> 2026-06-17 13:44:57.91 INFO::Min samples required with min abundance for a feature not to be filtered: 0.000000
#> 2026-06-17 13:44:57.91 INFO::Max samples allowed with min abundance for a feature not to be filtered: 186.850000
#> 2026-06-17 13:44:58.15 INFO::Total filtered features: 0
#> 2026-06-17 13:44:58.15 INFO::Filtered feature names from abundance, min prevalence, and max prevalence filtering:
#> 2026-06-17 13:44:58.17 INFO::Total features filtered by non-zero variance filtering: 208
#> 2026-06-17 13:44:58.17 INFO::Filtered feature names from variance filtering: ASV54, ASV108, ASV122, ASV152, ASV168, ASV207, ASV230, ASV242, ASV253, ASV280, ASV284, ASV298, ASV305, ASV325, ASV332, ASV336, ASV379, ASV381, ASV384, ASV386, ASV393, ASV394, ASV397, ASV406, ASV407, ASV408, ASV430, ASV432, ASV480, ASV484, ASV486, ASV487, ASV495, ASV497, ASV524, ASV528, ASV529, ASV532, ASV537, ASV547, ASV552, ASV557, ASV565, ASV568, ASV588, ASV603, ASV619, ASV620, ASV622, ASV630, ASV636, ASV651, ASV656, ASV658, ASV681, ASV682, ASV689, ASV699, ASV701, ASV716, ASV721, ASV739, ASV745, ASV751, ASV761, ASV763, ASV765, ASV778, ASV804, ASV806, ASV820, ASV826, ASV827, ASV833, ASV835, ASV849, ASV850, ASV851, ASV856, ASV862, ASV866, ASV867, ASV868, ASV869, ASV900, ASV908, ASV924, ASV925, ASV929, ASV931, ASV944, ASV946, ASV952, ASV956, ASV957, ASV960, ASV968, ASV978, ASV985, ASV1000, ASV1013, ASV1033, ASV1045, ASV1055, ASV1056, ASV1061, ASV1062, ASV1064, ASV1081, ASV1089, ASV1092, ASV1098, ASV1102, ASV1104, ASV1117, ASV1118, ASV1124, ASV1157, ASV1166, ASV1174, ASV1178, ASV1181, ASV1188, ASV1197, ASV1207, ASV1215, ASV1226, ASV1235, ASV1255, ASV1266, ASV1280, ASV1281, ASV1291, ASV1295, ASV1298, ASV1322, ASV1324, ASV1325, ASV1331, ASV1333, ASV1344, ASV1349, ASV1354, ASV1358, ASV1364, ASV1377, ASV1385, ASV1394, ASV1395, ASV1402, ASV1405, ASV1407, ASV1416, ASV1417, ASV1426, ASV1432, ASV1444, ASV1447, ASV1448, ASV1451, ASV1453, ASV1456, ASV1457, ASV1465, ASV1470, ASV1471, ASV1476, ASV1480, ASV1481, ASV1482, ASV1485, ASV1489, ASV1524, ASV1529, ASV1535, ASV1538, ASV1539, ASV1541, ASV1546, ASV1559, ASV1561, ASV1563, ASV1582, ASV1583, ASV1587, ASV1591, ASV1595, ASV1597, ASV1600, ASV1610, ASV1611, ASV1618, ASV1627, ASV1638, ASV1642, ASV1648, ASV1662, ASV1665, ASV1672, ASV1680, ASV1705, ASV1708, ASV1709, ASV1710, ASV1714, ASV1716, ASV1719, ASV1726
#> 2026-06-17 13:44:58.17 INFO::Writing filtered data to file res_maaslin3/features/filtered_data.tsv
#> 2026-06-17 13:44:58.21 INFO::Running selected transform method: LOG
#> 2026-06-17 13:44:58.22 INFO::Writing normalized, filtered, transformed data to file res_maaslin3/features/data_transformed.tsv
#> 2026-06-17 13:44:58.26 INFO::Factor detected for categorial metadata 'Height'. Using as-is.
#> 2026-06-17 13:44:58.26 INFO::Applying z-score to standardize continuous metadata
#> 2026-06-17 13:44:58.26 INFO::Running the linear model component
#> 2026-06-17 13:44:58.28 INFO::Fitting model to feature number 1, ASV2
#> 2026-06-17 13:44:58.29 INFO::Fitting model to feature number 2, ASV6
#> 2026-06-17 13:44:58.29 INFO::Fitting model to feature number 3, ASV7
#> 2026-06-17 13:44:58.29 INFO::Fitting model to feature number 4, ASV8
#> 2026-06-17 13:44:58.30 INFO::Fitting model to feature number 5, ASV10
#> 2026-06-17 13:44:58.30 INFO::Fitting model to feature number 6, ASV12
#> 2026-06-17 13:44:58.31 INFO::Fitting model to feature number 7, ASV13
#> 2026-06-17 13:44:58.31 INFO::Fitting model to feature number 8, ASV18
#> 2026-06-17 13:44:58.31 INFO::Fitting model to feature number 9, ASV19
#> 2026-06-17 13:44:58.32 INFO::Fitting model to feature number 10, ASV22
#> 2026-06-17 13:44:58.32 INFO::Fitting model to feature number 11, ASV23
#> 2026-06-17 13:44:58.32 WARNING::Fitting problem for feature 11 returning NA
#> 2026-06-17 13:44:58.33 INFO::Fitting model to feature number 12, ASV24
#> 2026-06-17 13:44:58.33 INFO::Fitting model to feature number 13, ASV25
#> 2026-06-17 13:44:58.33 INFO::Fitting model to feature number 14, ASV26
#> 2026-06-17 13:44:58.34 INFO::Fitting model to feature number 15, ASV27
#> 2026-06-17 13:44:58.34 INFO::Fitting model to feature number 16, ASV28
#> 2026-06-17 13:44:58.35 INFO::Fitting model to feature number 17, ASV29
#> 2026-06-17 13:44:58.35 INFO::Fitting model to feature number 18, ASV31
#> 2026-06-17 13:44:58.35 WARNING::Fitting problem for feature 18 returning NA
#> 2026-06-17 13:44:58.35 INFO::Fitting model to feature number 19, ASV32
#> 2026-06-17 13:44:58.36 INFO::Fitting model to feature number 20, ASV33
#> 2026-06-17 13:44:58.36 INFO::Fitting model to feature number 21, ASV34
#> 2026-06-17 13:44:58.37 INFO::Fitting model to feature number 22, ASV35
#> 2026-06-17 13:44:58.37 INFO::Fitting model to feature number 23, ASV38
#> 2026-06-17 13:44:58.37 INFO::Fitting model to feature number 24, ASV41
#> 2026-06-17 13:44:58.38 INFO::Fitting model to feature number 25, ASV42
#> 2026-06-17 13:44:58.38 INFO::Fitting model to feature number 26, ASV43
#> 2026-06-17 13:44:58.38 INFO::Fitting model to feature number 27, ASV45
#> 2026-06-17 13:44:58.39 INFO::Fitting model to feature number 28, ASV46
#> 2026-06-17 13:44:58.39 INFO::Fitting model to feature number 29, ASV47
#> 2026-06-17 13:44:58.39 INFO::Fitting model to feature number 30, ASV48
#> 2026-06-17 13:44:58.39 WARNING::Fitting problem for feature 30 returning NA
#> 2026-06-17 13:44:58.39 INFO::Fitting model to feature number 31, ASV49
#> 2026-06-17 13:44:58.40 INFO::Fitting model to feature number 32, ASV50
#> 2026-06-17 13:44:58.40 WARNING::Fitting problem for feature 32 returning NA
#> 2026-06-17 13:44:58.40 INFO::Fitting model to feature number 33, ASV51
#> 2026-06-17 13:44:58.41 INFO::Fitting model to feature number 34, ASV52
#> 2026-06-17 13:44:58.41 INFO::Fitting model to feature number 35, ASV53
#> 2026-06-17 13:44:58.41 INFO::Fitting model to feature number 36, ASV55
#> 2026-06-17 13:44:58.42 INFO::Fitting model to feature number 37, ASV56
#> 2026-06-17 13:44:58.42 INFO::Fitting model to feature number 38, ASV57
#> 2026-06-17 13:44:58.42 INFO::Fitting model to feature number 39, ASV58
#> 2026-06-17 13:44:58.43 INFO::Fitting model to feature number 40, ASV59
#> 2026-06-17 13:44:58.43 INFO::Fitting model to feature number 41, ASV60
#> 2026-06-17 13:44:58.43 INFO::Fitting model to feature number 42, ASV61
#> 2026-06-17 13:44:58.44 INFO::Fitting model to feature number 43, ASV62
#> 2026-06-17 13:44:58.44 INFO::Fitting model to feature number 44, ASV63
#> 2026-06-17 13:44:58.45 INFO::Fitting model to feature number 45, ASV64
#> 2026-06-17 13:44:58.45 INFO::Fitting model to feature number 46, ASV65
#> 2026-06-17 13:44:58.45 INFO::Fitting model to feature number 47, ASV67
#> 2026-06-17 13:44:58.45 INFO::Fitting model to feature number 48, ASV68
#> 2026-06-17 13:44:58.45 INFO::Fitting model to feature number 49, ASV69
#> 2026-06-17 13:44:58.46 INFO::Fitting model to feature number 50, ASV70
#> 2026-06-17 13:44:58.46 INFO::Fitting model to feature number 51, ASV71
#> 2026-06-17 13:44:58.47 INFO::Fitting model to feature number 52, ASV72
#> 2026-06-17 13:44:58.47 INFO::Fitting model to feature number 53, ASV75
#> 2026-06-17 13:44:58.47 INFO::Fitting model to feature number 54, ASV77
#> 2026-06-17 13:44:58.48 WARNING::Fitting problem for feature 54 returning NA
#> 2026-06-17 13:44:58.48 INFO::Fitting model to feature number 55, ASV78
#> 2026-06-17 13:44:58.48 INFO::Fitting model to feature number 56, ASV79
#> 2026-06-17 13:44:58.49 INFO::Fitting model to feature number 57, ASV80
#> 2026-06-17 13:44:58.49 INFO::Fitting model to feature number 58, ASV81
#> 2026-06-17 13:44:58.49 INFO::Fitting model to feature number 59, ASV82
#> 2026-06-17 13:44:58.50 INFO::Fitting model to feature number 60, ASV83
#> 2026-06-17 13:44:58.50 INFO::Fitting model to feature number 61, ASV84
#> 2026-06-17 13:44:58.50 INFO::Fitting model to feature number 62, ASV85
#> 2026-06-17 13:44:58.51 INFO::Fitting model to feature number 63, ASV87
#> 2026-06-17 13:44:58.51 INFO::Fitting model to feature number 64, ASV89
#> 2026-06-17 13:44:58.51 INFO::Fitting model to feature number 65, ASV90
#> 2026-06-17 13:44:58.51 INFO::Fitting model to feature number 66, ASV91
#> 2026-06-17 13:44:58.52 INFO::Fitting model to feature number 67, ASV92
#> 2026-06-17 13:44:58.52 INFO::Fitting model to feature number 68, ASV93
#> 2026-06-17 13:44:58.52 WARNING::Fitting problem for feature 68 returning NA
#> 2026-06-17 13:44:58.52 INFO::Fitting model to feature number 69, ASV94
#> 2026-06-17 13:44:58.53 INFO::Fitting model to feature number 70, ASV95
#> 2026-06-17 13:44:58.53 INFO::Fitting model to feature number 71, ASV96
#> 2026-06-17 13:44:58.54 INFO::Fitting model to feature number 72, ASV98
#> 2026-06-17 13:44:58.54 INFO::Fitting model to feature number 73, ASV99
#> 2026-06-17 13:44:58.54 INFO::Fitting model to feature number 74, ASV100
#> 2026-06-17 13:44:58.55 INFO::Fitting model to feature number 75, ASV101
#> 2026-06-17 13:44:58.55 INFO::Fitting model to feature number 76, ASV103
#> 2026-06-17 13:44:58.55 INFO::Fitting model to feature number 77, ASV104
#> 2026-06-17 13:44:58.56 INFO::Fitting model to feature number 78, ASV105
#> 2026-06-17 13:44:58.56 INFO::Fitting model to feature number 79, ASV106
#> 2026-06-17 13:44:58.56 INFO::Fitting model to feature number 80, ASV107
#> 2026-06-17 13:44:58.57 INFO::Fitting model to feature number 81, ASV109
#> 2026-06-17 13:44:58.57 INFO::Fitting model to feature number 82, ASV111
#> 2026-06-17 13:44:58.58 INFO::Fitting model to feature number 83, ASV112
#> 2026-06-17 13:44:58.58 INFO::Fitting model to feature number 84, ASV113
#> 2026-06-17 13:44:58.58 INFO::Fitting model to feature number 85, ASV115
#> 2026-06-17 13:44:58.59 INFO::Fitting model to feature number 86, ASV116
#> 2026-06-17 13:44:58.59 INFO::Fitting model to feature number 87, ASV118
#> 2026-06-17 13:44:58.59 INFO::Fitting model to feature number 88, ASV119
#> 2026-06-17 13:44:58.60 INFO::Fitting model to feature number 89, ASV120
#> 2026-06-17 13:44:58.60 INFO::Fitting model to feature number 90, ASV121
#> 2026-06-17 13:44:58.60 INFO::Fitting model to feature number 91, ASV124
#> 2026-06-17 13:44:58.61 INFO::Fitting model to feature number 92, ASV126
#> 2026-06-17 13:44:58.61 INFO::Fitting model to feature number 93, ASV127
#> 2026-06-17 13:44:58.61 WARNING::Fitting problem for feature 93 returning NA
#> 2026-06-17 13:44:58.62 INFO::Fitting model to feature number 94, ASV128
#> 2026-06-17 13:44:58.62 INFO::Fitting model to feature number 95, ASV129
#> 2026-06-17 13:44:58.62 INFO::Fitting model to feature number 96, ASV130
#> 2026-06-17 13:44:58.63 INFO::Fitting model to feature number 97, ASV131
#> 2026-06-17 13:44:58.63 INFO::Fitting model to feature number 98, ASV132
#> 2026-06-17 13:44:58.63 INFO::Fitting model to feature number 99, ASV133
#> 2026-06-17 13:44:58.64 INFO::Fitting model to feature number 100, ASV135
#> 2026-06-17 13:44:58.64 INFO::Fitting model to feature number 101, ASV136
#> 2026-06-17 13:44:58.64 INFO::Fitting model to feature number 102, ASV137
#> 2026-06-17 13:44:58.65 INFO::Fitting model to feature number 103, ASV138
#> 2026-06-17 13:44:58.65 INFO::Fitting model to feature number 104, ASV139
#> 2026-06-17 13:44:58.65 INFO::Fitting model to feature number 105, ASV142
#> 2026-06-17 13:44:58.66 INFO::Fitting model to feature number 106, ASV143
#> 2026-06-17 13:44:58.66 INFO::Fitting model to feature number 107, ASV144
#> 2026-06-17 13:44:58.66 INFO::Fitting model to feature number 108, ASV145
#> 2026-06-17 13:44:58.67 INFO::Fitting model to feature number 109, ASV146
#> 2026-06-17 13:44:58.67 INFO::Fitting model to feature number 110, ASV147
#> 2026-06-17 13:44:58.67 INFO::Fitting model to feature number 111, ASV148
#> 2026-06-17 13:44:58.68 INFO::Fitting model to feature number 112, ASV149
#> 2026-06-17 13:44:58.68 INFO::Fitting model to feature number 113, ASV150
#> 2026-06-17 13:44:58.68 INFO::Fitting model to feature number 114, ASV151
#> 2026-06-17 13:44:58.69 INFO::Fitting model to feature number 115, ASV153
#> 2026-06-17 13:44:58.69 INFO::Fitting model to feature number 116, ASV154
#> 2026-06-17 13:44:58.69 INFO::Fitting model to feature number 117, ASV157
#> 2026-06-17 13:44:58.69 INFO::Fitting model to feature number 118, ASV158
#> 2026-06-17 13:44:58.70 INFO::Fitting model to feature number 119, ASV159
#> 2026-06-17 13:44:58.70 INFO::Fitting model to feature number 120, ASV162
#> 2026-06-17 13:44:58.71 INFO::Fitting model to feature number 121, ASV163
#> 2026-06-17 13:44:58.71 INFO::Fitting model to feature number 122, ASV164
#> 2026-06-17 13:44:58.71 INFO::Fitting model to feature number 123, ASV166
#> 2026-06-17 13:44:58.72 INFO::Fitting model to feature number 124, ASV167
#> 2026-06-17 13:44:58.72 INFO::Fitting model to feature number 125, ASV170
#> 2026-06-17 13:44:58.72 INFO::Fitting model to feature number 126, ASV171
#> 2026-06-17 13:44:58.73 INFO::Fitting model to feature number 127, ASV172
#> 2026-06-17 13:44:58.73 INFO::Fitting model to feature number 128, ASV173
#> 2026-06-17 13:44:58.74 INFO::Fitting model to feature number 129, ASV175
#> 2026-06-17 13:44:58.74 INFO::Fitting model to feature number 130, ASV176
#> 2026-06-17 13:44:58.74 INFO::Fitting model to feature number 131, ASV178
#> 2026-06-17 13:44:58.75 INFO::Fitting model to feature number 132, ASV179
#> 2026-06-17 13:44:58.75 INFO::Fitting model to feature number 133, ASV181
#> 2026-06-17 13:44:58.75 INFO::Fitting model to feature number 134, ASV182
#> 2026-06-17 13:44:58.76 INFO::Fitting model to feature number 135, ASV183
#> 2026-06-17 13:44:58.76 INFO::Fitting model to feature number 136, ASV186
#> 2026-06-17 13:44:58.76 INFO::Fitting model to feature number 137, ASV187
#> 2026-06-17 13:44:58.77 INFO::Fitting model to feature number 138, ASV189
#> 2026-06-17 13:44:58.77 INFO::Fitting model to feature number 139, ASV190
#> 2026-06-17 13:44:58.77 INFO::Fitting model to feature number 140, ASV192
#> 2026-06-17 13:44:58.78 INFO::Fitting model to feature number 141, ASV193
#> 2026-06-17 13:44:58.78 INFO::Fitting model to feature number 142, ASV194
#> 2026-06-17 13:44:58.78 INFO::Fitting model to feature number 143, ASV195
#> 2026-06-17 13:44:58.78 INFO::Fitting model to feature number 144, ASV196
#> 2026-06-17 13:44:58.79 INFO::Fitting model to feature number 145, ASV197
#> 2026-06-17 13:44:58.79 INFO::Fitting model to feature number 146, ASV198
#> 2026-06-17 13:44:58.79 WARNING::Fitting problem for feature 146 returning NA
#> 2026-06-17 13:44:58.80 INFO::Fitting model to feature number 147, ASV199
#> 2026-06-17 13:44:58.80 INFO::Fitting model to feature number 148, ASV201
#> 2026-06-17 13:44:58.80 INFO::Fitting model to feature number 149, ASV202
#> 2026-06-17 13:44:58.80 INFO::Fitting model to feature number 150, ASV203
#> 2026-06-17 13:44:58.80 INFO::Fitting model to feature number 151, ASV205
#> 2026-06-17 13:44:58.81 INFO::Fitting model to feature number 152, ASV208
#> 2026-06-17 13:44:58.81 INFO::Fitting model to feature number 153, ASV209
#> 2026-06-17 13:44:58.81 INFO::Fitting model to feature number 154, ASV210
#> 2026-06-17 13:44:58.82 INFO::Fitting model to feature number 155, ASV211
#> 2026-06-17 13:44:58.82 INFO::Fitting model to feature number 156, ASV212
#> 2026-06-17 13:44:58.82 INFO::Fitting model to feature number 157, ASV214
#> 2026-06-17 13:44:58.82 INFO::Fitting model to feature number 158, ASV215
#> 2026-06-17 13:44:58.83 INFO::Fitting model to feature number 159, ASV216
#> 2026-06-17 13:44:58.83 INFO::Fitting model to feature number 160, ASV217
#> 2026-06-17 13:44:58.83 INFO::Fitting model to feature number 161, ASV219
#> 2026-06-17 13:44:58.84 INFO::Fitting model to feature number 162, ASV221
#> 2026-06-17 13:44:58.84 INFO::Fitting model to feature number 163, ASV222
#> 2026-06-17 13:44:58.84 WARNING::Fitting problem for feature 163 returning NA
#> 2026-06-17 13:44:58.85 INFO::Fitting model to feature number 164, ASV223
#> 2026-06-17 13:44:58.85 INFO::Fitting model to feature number 165, ASV224
#> 2026-06-17 13:44:58.85 INFO::Fitting model to feature number 166, ASV226
#> 2026-06-17 13:44:58.85 INFO::Fitting model to feature number 167, ASV227
#> 2026-06-17 13:44:58.86 INFO::Fitting model to feature number 168, ASV228
#> 2026-06-17 13:44:58.86 INFO::Fitting model to feature number 169, ASV229
#> 2026-06-17 13:44:58.86 INFO::Fitting model to feature number 170, ASV231
#> 2026-06-17 13:44:58.86 WARNING::Fitting problem for feature 170 returning NA
#> 2026-06-17 13:44:58.87 INFO::Fitting model to feature number 171, ASV233
#> 2026-06-17 13:44:58.87 INFO::Fitting model to feature number 172, ASV234
#> 2026-06-17 13:44:58.87 INFO::Fitting model to feature number 173, ASV235
#> 2026-06-17 13:44:58.88 INFO::Fitting model to feature number 174, ASV237
#> 2026-06-17 13:44:58.88 INFO::Fitting model to feature number 175, ASV238
#> 2026-06-17 13:44:58.89 INFO::Fitting model to feature number 176, ASV239
#> 2026-06-17 13:44:58.89 INFO::Fitting model to feature number 177, ASV240
#> 2026-06-17 13:44:58.89 INFO::Fitting model to feature number 178, ASV243
#> 2026-06-17 13:44:58.90 INFO::Fitting model to feature number 179, ASV244
#> 2026-06-17 13:44:58.90 INFO::Fitting model to feature number 180, ASV245
#> 2026-06-17 13:44:58.90 INFO::Fitting model to feature number 181, ASV246
#> 2026-06-17 13:44:58.90 INFO::Fitting model to feature number 182, ASV247
#> 2026-06-17 13:44:58.91 INFO::Fitting model to feature number 183, ASV248
#> 2026-06-17 13:44:58.91 INFO::Fitting model to feature number 184, ASV249
#> 2026-06-17 13:44:58.91 WARNING::Fitting problem for feature 184 returning NA
#> 2026-06-17 13:44:58.91 INFO::Fitting model to feature number 185, ASV251
#> 2026-06-17 13:44:58.92 INFO::Fitting model to feature number 186, ASV254
#> 2026-06-17 13:44:58.92 INFO::Fitting model to feature number 187, ASV255
#> 2026-06-17 13:44:58.92 INFO::Fitting model to feature number 188, ASV256
#> 2026-06-17 13:44:58.93 INFO::Fitting model to feature number 189, ASV257
#> 2026-06-17 13:44:58.93 INFO::Fitting model to feature number 190, ASV258
#> 2026-06-17 13:44:58.93 INFO::Fitting model to feature number 191, ASV260
#> 2026-06-17 13:44:58.93 INFO::Fitting model to feature number 192, ASV261
#> 2026-06-17 13:44:58.94 INFO::Fitting model to feature number 193, ASV262
#> 2026-06-17 13:44:58.94 INFO::Fitting model to feature number 194, ASV263
#> 2026-06-17 13:44:58.95 INFO::Fitting model to feature number 195, ASV264
#> 2026-06-17 13:44:58.95 INFO::Fitting model to feature number 196, ASV265
#> 2026-06-17 13:44:58.95 INFO::Fitting model to feature number 197, ASV266
#> 2026-06-17 13:44:58.96 INFO::Fitting model to feature number 198, ASV267
#> 2026-06-17 13:44:58.96 INFO::Fitting model to feature number 199, ASV270
#> 2026-06-17 13:44:58.96 INFO::Fitting model to feature number 200, ASV271
#> 2026-06-17 13:44:58.97 INFO::Fitting model to feature number 201, ASV272
#> 2026-06-17 13:44:58.97 INFO::Fitting model to feature number 202, ASV273
#> 2026-06-17 13:44:58.97 INFO::Fitting model to feature number 203, ASV274
#> 2026-06-17 13:44:58.98 INFO::Fitting model to feature number 204, ASV277
#> 2026-06-17 13:44:58.98 INFO::Fitting model to feature number 205, ASV278
#> 2026-06-17 13:44:58.98 INFO::Fitting model to feature number 206, ASV279
#> 2026-06-17 13:44:58.99 INFO::Fitting model to feature number 207, ASV282
#> 2026-06-17 13:44:58.99 INFO::Fitting model to feature number 208, ASV283
#> 2026-06-17 13:44:59.00 INFO::Fitting model to feature number 209, ASV285
#> 2026-06-17 13:44:59.00 INFO::Fitting model to feature number 210, ASV286
#> 2026-06-17 13:44:59.00 INFO::Fitting model to feature number 211, ASV287
#> 2026-06-17 13:44:59.00 INFO::Fitting model to feature number 212, ASV288
#> 2026-06-17 13:44:59.01 INFO::Fitting model to feature number 213, ASV292
#> 2026-06-17 13:44:59.01 INFO::Fitting model to feature number 214, ASV293
#> 2026-06-17 13:44:59.01 INFO::Fitting model to feature number 215, ASV294
#> 2026-06-17 13:44:59.01 INFO::Fitting model to feature number 216, ASV295
#> 2026-06-17 13:44:59.02 INFO::Fitting model to feature number 217, ASV297
#> 2026-06-17 13:44:59.02 INFO::Fitting model to feature number 218, ASV300
#> 2026-06-17 13:44:59.02 INFO::Fitting model to feature number 219, ASV302
#> 2026-06-17 13:44:59.03 INFO::Fitting model to feature number 220, ASV307
#> 2026-06-17 13:44:59.03 INFO::Fitting model to feature number 221, ASV309
#> 2026-06-17 13:44:59.03 INFO::Fitting model to feature number 222, ASV310
#> 2026-06-17 13:44:59.04 INFO::Fitting model to feature number 223, ASV311
#> 2026-06-17 13:44:59.04 INFO::Fitting model to feature number 224, ASV312
#> 2026-06-17 13:44:59.04 INFO::Fitting model to feature number 225, ASV313
#> 2026-06-17 13:44:59.05 INFO::Fitting model to feature number 226, ASV314
#> 2026-06-17 13:44:59.05 INFO::Fitting model to feature number 227, ASV315
#> 2026-06-17 13:44:59.05 INFO::Fitting model to feature number 228, ASV316
#> 2026-06-17 13:44:59.06 INFO::Fitting model to feature number 229, ASV318
#> 2026-06-17 13:44:59.06 INFO::Fitting model to feature number 230, ASV320
#> 2026-06-17 13:44:59.06 INFO::Fitting model to feature number 231, ASV322
#> 2026-06-17 13:44:59.06 INFO::Fitting model to feature number 232, ASV323
#> 2026-06-17 13:44:59.07 INFO::Fitting model to feature number 233, ASV327
#> 2026-06-17 13:44:59.07 INFO::Fitting model to feature number 234, ASV328
#> 2026-06-17 13:44:59.08 INFO::Fitting model to feature number 235, ASV329
#> 2026-06-17 13:44:59.08 INFO::Fitting model to feature number 236, ASV330
#> 2026-06-17 13:44:59.08 INFO::Fitting model to feature number 237, ASV333
#> 2026-06-17 13:44:59.09 INFO::Fitting model to feature number 238, ASV334
#> 2026-06-17 13:44:59.09 INFO::Fitting model to feature number 239, ASV337
#> 2026-06-17 13:44:59.09 INFO::Fitting model to feature number 240, ASV338
#> 2026-06-17 13:44:59.10 INFO::Fitting model to feature number 241, ASV339
#> 2026-06-17 13:44:59.10 INFO::Fitting model to feature number 242, ASV340
#> 2026-06-17 13:44:59.10 INFO::Fitting model to feature number 243, ASV341
#> 2026-06-17 13:44:59.11 INFO::Fitting model to feature number 244, ASV342
#> 2026-06-17 13:44:59.11 INFO::Fitting model to feature number 245, ASV344
#> 2026-06-17 13:44:59.12 INFO::Fitting model to feature number 246, ASV345
#> 2026-06-17 13:44:59.12 INFO::Fitting model to feature number 247, ASV346
#> 2026-06-17 13:44:59.12 INFO::Fitting model to feature number 248, ASV347
#> 2026-06-17 13:44:59.13 INFO::Fitting model to feature number 249, ASV348
#> 2026-06-17 13:44:59.13 INFO::Fitting model to feature number 250, ASV349
#> 2026-06-17 13:44:59.14 INFO::Fitting model to feature number 251, ASV350
#> 2026-06-17 13:44:59.14 INFO::Fitting model to feature number 252, ASV351
#> 2026-06-17 13:44:59.14 INFO::Fitting model to feature number 253, ASV353
#> 2026-06-17 13:44:59.15 INFO::Fitting model to feature number 254, ASV355
#> 2026-06-17 13:44:59.15 INFO::Fitting model to feature number 255, ASV356
#> 2026-06-17 13:44:59.16 INFO::Fitting model to feature number 256, ASV357
#> 2026-06-17 13:44:59.16 INFO::Fitting model to feature number 257, ASV358
#> 2026-06-17 13:44:59.16 INFO::Fitting model to feature number 258, ASV359
#> 2026-06-17 13:44:59.17 INFO::Fitting model to feature number 259, ASV360
#> 2026-06-17 13:44:59.17 INFO::Fitting model to feature number 260, ASV362
#> 2026-06-17 13:44:59.17 INFO::Fitting model to feature number 261, ASV364
#> 2026-06-17 13:44:59.18 INFO::Fitting model to feature number 262, ASV365
#> 2026-06-17 13:44:59.18 INFO::Fitting model to feature number 263, ASV366
#> 2026-06-17 13:44:59.19 INFO::Fitting model to feature number 264, ASV367
#> 2026-06-17 13:44:59.19 INFO::Fitting model to feature number 265, ASV368
#> 2026-06-17 13:44:59.19 INFO::Fitting model to feature number 266, ASV369
#> 2026-06-17 13:44:59.20 INFO::Fitting model to feature number 267, ASV371
#> 2026-06-17 13:44:59.20 INFO::Fitting model to feature number 268, ASV372
#> 2026-06-17 13:44:59.20 INFO::Fitting model to feature number 269, ASV373
#> 2026-06-17 13:44:59.21 WARNING::Fitting problem for feature 269 returning NA
#> 2026-06-17 13:44:59.21 INFO::Fitting model to feature number 270, ASV375
#> 2026-06-17 13:44:59.21 INFO::Fitting model to feature number 271, ASV376
#> 2026-06-17 13:44:59.22 INFO::Fitting model to feature number 272, ASV377
#> 2026-06-17 13:44:59.22 INFO::Fitting model to feature number 273, ASV378
#> 2026-06-17 13:44:59.23 INFO::Fitting model to feature number 274, ASV380
#> 2026-06-17 13:44:59.23 INFO::Fitting model to feature number 275, ASV382
#> 2026-06-17 13:44:59.23 INFO::Fitting model to feature number 276, ASV383
#> 2026-06-17 13:44:59.24 INFO::Fitting model to feature number 277, ASV385
#> 2026-06-17 13:44:59.24 INFO::Fitting model to feature number 278, ASV388
#> 2026-06-17 13:44:59.24 INFO::Fitting model to feature number 279, ASV390
#> 2026-06-17 13:44:59.24 INFO::Fitting model to feature number 280, ASV395
#> 2026-06-17 13:44:59.25 INFO::Fitting model to feature number 281, ASV399
#> 2026-06-17 13:44:59.25 INFO::Fitting model to feature number 282, ASV401
#> 2026-06-17 13:44:59.25 INFO::Fitting model to feature number 283, ASV402
#> 2026-06-17 13:44:59.25 INFO::Fitting model to feature number 284, ASV404
#> 2026-06-17 13:44:59.26 INFO::Fitting model to feature number 285, ASV405
#> 2026-06-17 13:44:59.26 INFO::Fitting model to feature number 286, ASV409
#> 2026-06-17 13:44:59.27 INFO::Fitting model to feature number 287, ASV411
#> 2026-06-17 13:44:59.27 INFO::Fitting model to feature number 288, ASV412
#> 2026-06-17 13:44:59.27 INFO::Fitting model to feature number 289, ASV413
#> 2026-06-17 13:44:59.28 INFO::Fitting model to feature number 290, ASV414
#> 2026-06-17 13:44:59.28 WARNING::Fitting problem for feature 290 returning NA
#> 2026-06-17 13:44:59.28 INFO::Fitting model to feature number 291, ASV415
#> 2026-06-17 13:44:59.28 INFO::Fitting model to feature number 292, ASV416
#> 2026-06-17 13:44:59.29 INFO::Fitting model to feature number 293, ASV417
#> 2026-06-17 13:44:59.29 INFO::Fitting model to feature number 294, ASV418
#> 2026-06-17 13:44:59.29 INFO::Fitting model to feature number 295, ASV419
#> 2026-06-17 13:44:59.30 INFO::Fitting model to feature number 296, ASV420
#> 2026-06-17 13:44:59.30 INFO::Fitting model to feature number 297, ASV422
#> 2026-06-17 13:44:59.30 INFO::Fitting model to feature number 298, ASV423
#> 2026-06-17 13:44:59.31 INFO::Fitting model to feature number 299, ASV424
#> 2026-06-17 13:44:59.31 INFO::Fitting model to feature number 300, ASV426
#> 2026-06-17 13:44:59.31 INFO::Fitting model to feature number 301, ASV427
#> 2026-06-17 13:44:59.32 INFO::Fitting model to feature number 302, ASV428
#> 2026-06-17 13:44:59.32 INFO::Fitting model to feature number 303, ASV429
#> 2026-06-17 13:44:59.32 INFO::Fitting model to feature number 304, ASV431
#> 2026-06-17 13:44:59.32 INFO::Fitting model to feature number 305, ASV433
#> 2026-06-17 13:44:59.33 INFO::Fitting model to feature number 306, ASV434
#> 2026-06-17 13:44:59.33 INFO::Fitting model to feature number 307, ASV438
#> 2026-06-17 13:44:59.33 INFO::Fitting model to feature number 308, ASV439
#> 2026-06-17 13:44:59.34 INFO::Fitting model to feature number 309, ASV440
#> 2026-06-17 13:44:59.34 INFO::Fitting model to feature number 310, ASV442
#> 2026-06-17 13:44:59.34 INFO::Fitting model to feature number 311, ASV443
#> 2026-06-17 13:44:59.35 INFO::Fitting model to feature number 312, ASV444
#> 2026-06-17 13:44:59.35 INFO::Fitting model to feature number 313, ASV445
#> 2026-06-17 13:44:59.35 INFO::Fitting model to feature number 314, ASV446
#> 2026-06-17 13:44:59.36 INFO::Fitting model to feature number 315, ASV448
#> 2026-06-17 13:44:59.36 INFO::Fitting model to feature number 316, ASV449
#> 2026-06-17 13:44:59.36 INFO::Fitting model to feature number 317, ASV451
#> 2026-06-17 13:44:59.37 INFO::Fitting model to feature number 318, ASV452
#> 2026-06-17 13:44:59.37 INFO::Fitting model to feature number 319, ASV453
#> 2026-06-17 13:44:59.37 INFO::Fitting model to feature number 320, ASV454
#> 2026-06-17 13:44:59.38 INFO::Fitting model to feature number 321, ASV456
#> 2026-06-17 13:44:59.38 INFO::Fitting model to feature number 322, ASV457
#> 2026-06-17 13:44:59.38 INFO::Fitting model to feature number 323, ASV458
#> 2026-06-17 13:44:59.38 INFO::Fitting model to feature number 324, ASV459
#> 2026-06-17 13:44:59.39 INFO::Fitting model to feature number 325, ASV460
#> 2026-06-17 13:44:59.39 INFO::Fitting model to feature number 326, ASV461
#> 2026-06-17 13:44:59.39 INFO::Fitting model to feature number 327, ASV462
#> 2026-06-17 13:44:59.40 INFO::Fitting model to feature number 328, ASV463
#> 2026-06-17 13:44:59.40 INFO::Fitting model to feature number 329, ASV464
#> 2026-06-17 13:44:59.40 INFO::Fitting model to feature number 330, ASV465
#> 2026-06-17 13:44:59.41 WARNING::Fitting problem for feature 330 returning NA
#> 2026-06-17 13:44:59.41 INFO::Fitting model to feature number 331, ASV466
#> 2026-06-17 13:44:59.41 INFO::Fitting model to feature number 332, ASV467
#> 2026-06-17 13:44:59.42 INFO::Fitting model to feature number 333, ASV468
#> 2026-06-17 13:44:59.42 INFO::Fitting model to feature number 334, ASV469
#> 2026-06-17 13:44:59.42 INFO::Fitting model to feature number 335, ASV470
#> 2026-06-17 13:44:59.43 INFO::Fitting model to feature number 336, ASV472
#> 2026-06-17 13:44:59.43 INFO::Fitting model to feature number 337, ASV473
#> 2026-06-17 13:44:59.43 INFO::Fitting model to feature number 338, ASV474
#> 2026-06-17 13:44:59.44 INFO::Fitting model to feature number 339, ASV475
#> 2026-06-17 13:44:59.44 INFO::Fitting model to feature number 340, ASV476
#> 2026-06-17 13:44:59.44 INFO::Fitting model to feature number 341, ASV477
#> 2026-06-17 13:44:59.45 INFO::Fitting model to feature number 342, ASV478
#> 2026-06-17 13:44:59.54 INFO::Fitting model to feature number 343, ASV479
#> 2026-06-17 13:44:59.55 INFO::Fitting model to feature number 344, ASV481
#> 2026-06-17 13:44:59.56 INFO::Fitting model to feature number 345, ASV483
#> 2026-06-17 13:44:59.56 INFO::Fitting model to feature number 346, ASV488
#> 2026-06-17 13:44:59.56 INFO::Fitting model to feature number 347, ASV489
#> 2026-06-17 13:44:59.56 INFO::Fitting model to feature number 348, ASV490
#> 2026-06-17 13:44:59.57 INFO::Fitting model to feature number 349, ASV491
#> 2026-06-17 13:44:59.57 INFO::Fitting model to feature number 350, ASV492
#> 2026-06-17 13:44:59.58 INFO::Fitting model to feature number 351, ASV493
#> 2026-06-17 13:44:59.58 INFO::Fitting model to feature number 352, ASV494
#> 2026-06-17 13:44:59.58 INFO::Fitting model to feature number 353, ASV496
#> 2026-06-17 13:44:59.59 INFO::Fitting model to feature number 354, ASV498
#> 2026-06-17 13:44:59.59 INFO::Fitting model to feature number 355, ASV499
#> 2026-06-17 13:44:59.59 INFO::Fitting model to feature number 356, ASV500
#> 2026-06-17 13:44:59.60 INFO::Fitting model to feature number 357, ASV501
#> 2026-06-17 13:44:59.60 INFO::Fitting model to feature number 358, ASV502
#> 2026-06-17 13:44:59.61 INFO::Fitting model to feature number 359, ASV504
#> 2026-06-17 13:44:59.61 INFO::Fitting model to feature number 360, ASV505
#> 2026-06-17 13:44:59.61 INFO::Fitting model to feature number 361, ASV507
#> 2026-06-17 13:44:59.62 INFO::Fitting model to feature number 362, ASV508
#> 2026-06-17 13:44:59.62 INFO::Fitting model to feature number 363, ASV509
#> 2026-06-17 13:44:59.63 INFO::Fitting model to feature number 364, ASV511
#> 2026-06-17 13:44:59.63 INFO::Fitting model to feature number 365, ASV512
#> 2026-06-17 13:44:59.63 INFO::Fitting model to feature number 366, ASV514
#> 2026-06-17 13:44:59.64 INFO::Fitting model to feature number 367, ASV515
#> 2026-06-17 13:44:59.64 INFO::Fitting model to feature number 368, ASV516
#> 2026-06-17 13:44:59.64 INFO::Fitting model to feature number 369, ASV517
#> 2026-06-17 13:44:59.65 INFO::Fitting model to feature number 370, ASV519
#> 2026-06-17 13:44:59.65 INFO::Fitting model to feature number 371, ASV520
#> 2026-06-17 13:44:59.65 INFO::Fitting model to feature number 372, ASV521
#> 2026-06-17 13:44:59.66 INFO::Fitting model to feature number 373, ASV522
#> 2026-06-17 13:44:59.66 INFO::Fitting model to feature number 374, ASV523
#> 2026-06-17 13:44:59.66 INFO::Fitting model to feature number 375, ASV526
#> 2026-06-17 13:44:59.67 INFO::Fitting model to feature number 376, ASV527
#> 2026-06-17 13:44:59.67 INFO::Fitting model to feature number 377, ASV530
#> 2026-06-17 13:44:59.67 INFO::Fitting model to feature number 378, ASV531
#> 2026-06-17 13:44:59.67 INFO::Fitting model to feature number 379, ASV533
#> 2026-06-17 13:44:59.68 INFO::Fitting model to feature number 380, ASV534
#> 2026-06-17 13:44:59.68 INFO::Fitting model to feature number 381, ASV535
#> 2026-06-17 13:44:59.69 INFO::Fitting model to feature number 382, ASV536
#> 2026-06-17 13:44:59.69 INFO::Fitting model to feature number 383, ASV538
#> 2026-06-17 13:44:59.69 INFO::Fitting model to feature number 384, ASV539
#> 2026-06-17 13:44:59.70 INFO::Fitting model to feature number 385, ASV540
#> 2026-06-17 13:44:59.70 INFO::Fitting model to feature number 386, ASV541
#> 2026-06-17 13:44:59.70 INFO::Fitting model to feature number 387, ASV542
#> 2026-06-17 13:44:59.71 INFO::Fitting model to feature number 388, ASV543
#> 2026-06-17 13:44:59.71 INFO::Fitting model to feature number 389, ASV544
#> 2026-06-17 13:44:59.71 INFO::Fitting model to feature number 390, ASV545
#> 2026-06-17 13:44:59.72 INFO::Fitting model to feature number 391, ASV546
#> 2026-06-17 13:44:59.72 INFO::Fitting model to feature number 392, ASV548
#> 2026-06-17 13:44:59.72 INFO::Fitting model to feature number 393, ASV549
#> 2026-06-17 13:44:59.72 INFO::Fitting model to feature number 394, ASV550
#> 2026-06-17 13:44:59.73 INFO::Fitting model to feature number 395, ASV551
#> 2026-06-17 13:44:59.73 INFO::Fitting model to feature number 396, ASV553
#> 2026-06-17 13:44:59.74 INFO::Fitting model to feature number 397, ASV554
#> 2026-06-17 13:44:59.74 INFO::Fitting model to feature number 398, ASV555
#> 2026-06-17 13:44:59.74 INFO::Fitting model to feature number 399, ASV556
#> 2026-06-17 13:44:59.74 INFO::Fitting model to feature number 400, ASV558
#> 2026-06-17 13:44:59.75 INFO::Fitting model to feature number 401, ASV559
#> 2026-06-17 13:44:59.75 INFO::Fitting model to feature number 402, ASV560
#> 2026-06-17 13:44:59.75 INFO::Fitting model to feature number 403, ASV561
#> 2026-06-17 13:44:59.76 INFO::Fitting model to feature number 404, ASV562
#> 2026-06-17 13:44:59.76 INFO::Fitting model to feature number 405, ASV563
#> 2026-06-17 13:44:59.76 INFO::Fitting model to feature number 406, ASV564
#> 2026-06-17 13:44:59.77 INFO::Fitting model to feature number 407, ASV566
#> 2026-06-17 13:44:59.77 INFO::Fitting model to feature number 408, ASV567
#> 2026-06-17 13:44:59.77 INFO::Fitting model to feature number 409, ASV569
#> 2026-06-17 13:44:59.78 INFO::Fitting model to feature number 410, ASV570
#> 2026-06-17 13:44:59.78 INFO::Fitting model to feature number 411, ASV571
#> 2026-06-17 13:44:59.78 INFO::Fitting model to feature number 412, ASV572
#> 2026-06-17 13:44:59.79 INFO::Fitting model to feature number 413, ASV573
#> 2026-06-17 13:44:59.79 INFO::Fitting model to feature number 414, ASV574
#> 2026-06-17 13:44:59.80 INFO::Fitting model to feature number 415, ASV576
#> 2026-06-17 13:44:59.80 INFO::Fitting model to feature number 416, ASV577
#> 2026-06-17 13:44:59.80 INFO::Fitting model to feature number 417, ASV578
#> 2026-06-17 13:44:59.81 WARNING::Fitting problem for feature 417 returning NA
#> 2026-06-17 13:44:59.81 INFO::Fitting model to feature number 418, ASV579
#> 2026-06-17 13:44:59.81 INFO::Fitting model to feature number 419, ASV580
#> 2026-06-17 13:44:59.81 INFO::Fitting model to feature number 420, ASV581
#> 2026-06-17 13:44:59.81 INFO::Fitting model to feature number 421, ASV582
#> 2026-06-17 13:44:59.82 WARNING::Fitting problem for feature 421 returning NA
#> 2026-06-17 13:44:59.82 INFO::Fitting model to feature number 422, ASV584
#> 2026-06-17 13:44:59.82 INFO::Fitting model to feature number 423, ASV585
#> 2026-06-17 13:44:59.83 INFO::Fitting model to feature number 424, ASV586
#> 2026-06-17 13:44:59.83 INFO::Fitting model to feature number 425, ASV587
#> 2026-06-17 13:44:59.83 INFO::Fitting model to feature number 426, ASV589
#> 2026-06-17 13:44:59.84 INFO::Fitting model to feature number 427, ASV590
#> 2026-06-17 13:44:59.84 INFO::Fitting model to feature number 428, ASV591
#> 2026-06-17 13:44:59.84 INFO::Fitting model to feature number 429, ASV592
#> 2026-06-17 13:44:59.85 INFO::Fitting model to feature number 430, ASV593
#> 2026-06-17 13:44:59.85 INFO::Fitting model to feature number 431, ASV594
#> 2026-06-17 13:44:59.85 INFO::Fitting model to feature number 432, ASV595
#> 2026-06-17 13:44:59.85 INFO::Fitting model to feature number 433, ASV596
#> 2026-06-17 13:44:59.86 INFO::Fitting model to feature number 434, ASV597
#> 2026-06-17 13:44:59.86 INFO::Fitting model to feature number 435, ASV598
#> 2026-06-17 13:44:59.87 INFO::Fitting model to feature number 436, ASV599
#> 2026-06-17 13:44:59.87 INFO::Fitting model to feature number 437, ASV600
#> 2026-06-17 13:44:59.87 INFO::Fitting model to feature number 438, ASV602
#> 2026-06-17 13:44:59.88 INFO::Fitting model to feature number 439, ASV604
#> 2026-06-17 13:44:59.88 INFO::Fitting model to feature number 440, ASV605
#> 2026-06-17 13:44:59.88 INFO::Fitting model to feature number 441, ASV607
#> 2026-06-17 13:44:59.89 INFO::Fitting model to feature number 442, ASV608
#> 2026-06-17 13:44:59.89 INFO::Fitting model to feature number 443, ASV610
#> 2026-06-17 13:44:59.89 INFO::Fitting model to feature number 444, ASV612
#> 2026-06-17 13:44:59.90 INFO::Fitting model to feature number 445, ASV614
#> 2026-06-17 13:44:59.90 INFO::Fitting model to feature number 446, ASV615
#> 2026-06-17 13:44:59.91 INFO::Fitting model to feature number 447, ASV616
#> 2026-06-17 13:44:59.91 INFO::Fitting model to feature number 448, ASV617
#> 2026-06-17 13:44:59.91 INFO::Fitting model to feature number 449, ASV618
#> 2026-06-17 13:44:59.92 INFO::Fitting model to feature number 450, ASV621
#> 2026-06-17 13:44:59.92 INFO::Fitting model to feature number 451, ASV624
#> 2026-06-17 13:44:59.93 INFO::Fitting model to feature number 452, ASV625
#> 2026-06-17 13:44:59.93 INFO::Fitting model to feature number 453, ASV626
#> 2026-06-17 13:44:59.93 INFO::Fitting model to feature number 454, ASV627
#> 2026-06-17 13:44:59.94 INFO::Fitting model to feature number 455, ASV628
#> 2026-06-17 13:44:59.94 INFO::Fitting model to feature number 456, ASV629
#> 2026-06-17 13:44:59.94 INFO::Fitting model to feature number 457, ASV631
#> 2026-06-17 13:44:59.95 INFO::Fitting model to feature number 458, ASV632
#> 2026-06-17 13:44:59.95 INFO::Fitting model to feature number 459, ASV633
#> 2026-06-17 13:44:59.95 INFO::Fitting model to feature number 460, ASV634
#> 2026-06-17 13:44:59.96 INFO::Fitting model to feature number 461, ASV635
#> 2026-06-17 13:44:59.96 INFO::Fitting model to feature number 462, ASV637
#> 2026-06-17 13:44:59.96 INFO::Fitting model to feature number 463, ASV639
#> 2026-06-17 13:44:59.97 INFO::Fitting model to feature number 464, ASV640
#> 2026-06-17 13:44:59.97 INFO::Fitting model to feature number 465, ASV641
#> 2026-06-17 13:44:59.97 INFO::Fitting model to feature number 466, ASV643
#> 2026-06-17 13:44:59.98 INFO::Fitting model to feature number 467, ASV644
#> 2026-06-17 13:44:59.98 INFO::Fitting model to feature number 468, ASV645
#> 2026-06-17 13:44:59.98 INFO::Fitting model to feature number 469, ASV647
#> 2026-06-17 13:44:59.99 INFO::Fitting model to feature number 470, ASV649
#> 2026-06-17 13:44:59.99 INFO::Fitting model to feature number 471, ASV650
#> 2026-06-17 13:44:59.99 INFO::Fitting model to feature number 472, ASV652
#> 2026-06-17 13:44:59.99 INFO::Fitting model to feature number 473, ASV653
#> 2026-06-17 13:45:00.00 INFO::Fitting model to feature number 474, ASV654
#> 2026-06-17 13:45:00.00 INFO::Fitting model to feature number 475, ASV657
#> 2026-06-17 13:45:00.00 INFO::Fitting model to feature number 476, ASV659
#> 2026-06-17 13:45:00.00 INFO::Fitting model to feature number 477, ASV660
#> 2026-06-17 13:45:00.01 INFO::Fitting model to feature number 478, ASV661
#> 2026-06-17 13:45:00.01 INFO::Fitting model to feature number 479, ASV662
#> 2026-06-17 13:45:00.01 INFO::Fitting model to feature number 480, ASV664
#> 2026-06-17 13:45:00.02 INFO::Fitting model to feature number 481, ASV665
#> 2026-06-17 13:45:00.02 INFO::Fitting model to feature number 482, ASV666
#> 2026-06-17 13:45:00.02 INFO::Fitting model to feature number 483, ASV667
#> 2026-06-17 13:45:00.03 INFO::Fitting model to feature number 484, ASV668
#> 2026-06-17 13:45:00.03 INFO::Fitting model to feature number 485, ASV669
#> 2026-06-17 13:45:00.03 INFO::Fitting model to feature number 486, ASV670
#> 2026-06-17 13:45:00.04 INFO::Fitting model to feature number 487, ASV671
#> 2026-06-17 13:45:00.04 INFO::Fitting model to feature number 488, ASV672
#> 2026-06-17 13:45:00.04 INFO::Fitting model to feature number 489, ASV673
#> 2026-06-17 13:45:00.05 INFO::Fitting model to feature number 490, ASV674
#> 2026-06-17 13:45:00.05 INFO::Fitting model to feature number 491, ASV675
#> 2026-06-17 13:45:00.05 INFO::Fitting model to feature number 492, ASV676
#> 2026-06-17 13:45:00.06 INFO::Fitting model to feature number 493, ASV677
#> 2026-06-17 13:45:00.06 INFO::Fitting model to feature number 494, ASV678
#> 2026-06-17 13:45:00.06 INFO::Fitting model to feature number 495, ASV679
#> 2026-06-17 13:45:00.06 INFO::Fitting model to feature number 496, ASV680
#> 2026-06-17 13:45:00.07 INFO::Fitting model to feature number 497, ASV683
#> 2026-06-17 13:45:00.07 INFO::Fitting model to feature number 498, ASV684
#> 2026-06-17 13:45:00.07 INFO::Fitting model to feature number 499, ASV685
#> 2026-06-17 13:45:00.07 INFO::Fitting model to feature number 500, ASV686
#> 2026-06-17 13:45:00.08 INFO::Fitting model to feature number 501, ASV687
#> 2026-06-17 13:45:00.08 INFO::Fitting model to feature number 502, ASV688
#> 2026-06-17 13:45:00.08 INFO::Fitting model to feature number 503, ASV691
#> 2026-06-17 13:45:00.09 INFO::Fitting model to feature number 504, ASV692
#> 2026-06-17 13:45:00.09 INFO::Fitting model to feature number 505, ASV693
#> 2026-06-17 13:45:00.09 INFO::Fitting model to feature number 506, ASV694
#> 2026-06-17 13:45:00.09 INFO::Fitting model to feature number 507, ASV695
#> 2026-06-17 13:45:00.10 INFO::Fitting model to feature number 508, ASV696
#> 2026-06-17 13:45:00.10 INFO::Fitting model to feature number 509, ASV697
#> 2026-06-17 13:45:00.10 INFO::Fitting model to feature number 510, ASV698
#> 2026-06-17 13:45:00.11 INFO::Fitting model to feature number 511, ASV700
#> 2026-06-17 13:45:00.11 INFO::Fitting model to feature number 512, ASV702
#> 2026-06-17 13:45:00.12 INFO::Fitting model to feature number 513, ASV703
#> 2026-06-17 13:45:00.12 INFO::Fitting model to feature number 514, ASV704
#> 2026-06-17 13:45:00.12 INFO::Fitting model to feature number 515, ASV705
#> 2026-06-17 13:45:00.13 INFO::Fitting model to feature number 516, ASV707
#> 2026-06-17 13:45:00.13 INFO::Fitting model to feature number 517, ASV710
#> 2026-06-17 13:45:00.13 INFO::Fitting model to feature number 518, ASV711
#> 2026-06-17 13:45:00.13 INFO::Fitting model to feature number 519, ASV712
#> 2026-06-17 13:45:00.14 INFO::Fitting model to feature number 520, ASV713
#> 2026-06-17 13:45:00.14 INFO::Fitting model to feature number 521, ASV714
#> 2026-06-17 13:45:00.14 INFO::Fitting model to feature number 522, ASV715
#> 2026-06-17 13:45:00.15 INFO::Fitting model to feature number 523, ASV717
#> 2026-06-17 13:45:00.15 INFO::Fitting model to feature number 524, ASV718
#> 2026-06-17 13:45:00.15 INFO::Fitting model to feature number 525, ASV720
#> 2026-06-17 13:45:00.16 INFO::Fitting model to feature number 526, ASV722
#> 2026-06-17 13:45:00.16 INFO::Fitting model to feature number 527, ASV723
#> 2026-06-17 13:45:00.16 INFO::Fitting model to feature number 528, ASV724
#> 2026-06-17 13:45:00.17 INFO::Fitting model to feature number 529, ASV726
#> 2026-06-17 13:45:00.17 INFO::Fitting model to feature number 530, ASV727
#> 2026-06-17 13:45:00.18 INFO::Fitting model to feature number 531, ASV728
#> 2026-06-17 13:45:00.18 WARNING::Fitting problem for feature 531 returning NA
#> 2026-06-17 13:45:00.18 INFO::Fitting model to feature number 532, ASV729
#> 2026-06-17 13:45:00.18 INFO::Fitting model to feature number 533, ASV730
#> 2026-06-17 13:45:00.19 INFO::Fitting model to feature number 534, ASV731
#> 2026-06-17 13:45:00.19 INFO::Fitting model to feature number 535, ASV732
#> 2026-06-17 13:45:00.19 INFO::Fitting model to feature number 536, ASV733
#> 2026-06-17 13:45:00.20 INFO::Fitting model to feature number 537, ASV734
#> 2026-06-17 13:45:00.20 INFO::Fitting model to feature number 538, ASV736
#> 2026-06-17 13:45:00.20 INFO::Fitting model to feature number 539, ASV737
#> 2026-06-17 13:45:00.20 INFO::Fitting model to feature number 540, ASV738
#> 2026-06-17 13:45:00.21 INFO::Fitting model to feature number 541, ASV740
#> 2026-06-17 13:45:00.21 INFO::Fitting model to feature number 542, ASV741
#> 2026-06-17 13:45:00.22 INFO::Fitting model to feature number 543, ASV742
#> 2026-06-17 13:45:00.22 INFO::Fitting model to feature number 544, ASV743
#> 2026-06-17 13:45:00.22 INFO::Fitting model to feature number 545, ASV744
#> 2026-06-17 13:45:00.22 INFO::Fitting model to feature number 546, ASV746
#> 2026-06-17 13:45:00.23 INFO::Fitting model to feature number 547, ASV747
#> 2026-06-17 13:45:00.23 INFO::Fitting model to feature number 548, ASV748
#> 2026-06-17 13:45:00.23 INFO::Fitting model to feature number 549, ASV749
#> 2026-06-17 13:45:00.24 INFO::Fitting model to feature number 550, ASV752
#> 2026-06-17 13:45:00.24 INFO::Fitting model to feature number 551, ASV753
#> 2026-06-17 13:45:00.24 INFO::Fitting model to feature number 552, ASV754
#> 2026-06-17 13:45:00.25 INFO::Fitting model to feature number 553, ASV755
#> 2026-06-17 13:45:00.25 INFO::Fitting model to feature number 554, ASV756
#> 2026-06-17 13:45:00.25 INFO::Fitting model to feature number 555, ASV757
#> 2026-06-17 13:45:00.26 INFO::Fitting model to feature number 556, ASV758
#> 2026-06-17 13:45:00.26 INFO::Fitting model to feature number 557, ASV760
#> 2026-06-17 13:45:00.26 INFO::Fitting model to feature number 558, ASV762
#> 2026-06-17 13:45:00.27 INFO::Fitting model to feature number 559, ASV764
#> 2026-06-17 13:45:00.27 INFO::Fitting model to feature number 560, ASV766
#> 2026-06-17 13:45:00.28 INFO::Fitting model to feature number 561, ASV767
#> 2026-06-17 13:45:00.28 INFO::Fitting model to feature number 562, ASV768
#> 2026-06-17 13:45:00.28 INFO::Fitting model to feature number 563, ASV769
#> 2026-06-17 13:45:00.29 INFO::Fitting model to feature number 564, ASV770
#> 2026-06-17 13:45:00.29 INFO::Fitting model to feature number 565, ASV771
#> 2026-06-17 13:45:00.29 INFO::Fitting model to feature number 566, ASV772
#> 2026-06-17 13:45:00.30 INFO::Fitting model to feature number 567, ASV773
#> 2026-06-17 13:45:00.30 INFO::Fitting model to feature number 568, ASV774
#> 2026-06-17 13:45:00.31 INFO::Fitting model to feature number 569, ASV775
#> 2026-06-17 13:45:00.31 INFO::Fitting model to feature number 570, ASV776
#> 2026-06-17 13:45:00.31 INFO::Fitting model to feature number 571, ASV777
#> 2026-06-17 13:45:00.31 INFO::Fitting model to feature number 572, ASV779
#> 2026-06-17 13:45:00.32 INFO::Fitting model to feature number 573, ASV780
#> 2026-06-17 13:45:00.32 INFO::Fitting model to feature number 574, ASV781
#> 2026-06-17 13:45:00.33 INFO::Fitting model to feature number 575, ASV782
#> 2026-06-17 13:45:00.33 INFO::Fitting model to feature number 576, ASV783
#> 2026-06-17 13:45:00.34 INFO::Fitting model to feature number 577, ASV784
#> 2026-06-17 13:45:00.34 INFO::Fitting model to feature number 578, ASV785
#> 2026-06-17 13:45:00.34 INFO::Fitting model to feature number 579, ASV786
#> 2026-06-17 13:45:00.35 WARNING::Fitting problem for feature 579 returning NA
#> 2026-06-17 13:45:00.35 INFO::Fitting model to feature number 580, ASV787
#> 2026-06-17 13:45:00.35 INFO::Fitting model to feature number 581, ASV788
#> 2026-06-17 13:45:00.36 INFO::Fitting model to feature number 582, ASV790
#> 2026-06-17 13:45:00.36 INFO::Fitting model to feature number 583, ASV792
#> 2026-06-17 13:45:00.36 INFO::Fitting model to feature number 584, ASV795
#> 2026-06-17 13:45:00.36 INFO::Fitting model to feature number 585, ASV796
#> 2026-06-17 13:45:00.37 INFO::Fitting model to feature number 586, ASV797
#> 2026-06-17 13:45:00.37 INFO::Fitting model to feature number 587, ASV798
#> 2026-06-17 13:45:00.37 INFO::Fitting model to feature number 588, ASV799
#> 2026-06-17 13:45:00.38 INFO::Fitting model to feature number 589, ASV801
#> 2026-06-17 13:45:00.38 INFO::Fitting model to feature number 590, ASV802
#> 2026-06-17 13:45:00.38 INFO::Fitting model to feature number 591, ASV803
#> 2026-06-17 13:45:00.39 INFO::Fitting model to feature number 592, ASV805
#> 2026-06-17 13:45:00.39 INFO::Fitting model to feature number 593, ASV807
#> 2026-06-17 13:45:00.39 INFO::Fitting model to feature number 594, ASV808
#> 2026-06-17 13:45:00.40 INFO::Fitting model to feature number 595, ASV810
#> 2026-06-17 13:45:00.40 INFO::Fitting model to feature number 596, ASV811
#> 2026-06-17 13:45:00.40 INFO::Fitting model to feature number 597, ASV814
#> 2026-06-17 13:45:00.40 INFO::Fitting model to feature number 598, ASV815
#> 2026-06-17 13:45:00.41 INFO::Fitting model to feature number 599, ASV816
#> 2026-06-17 13:45:00.41 INFO::Fitting model to feature number 600, ASV817
#> 2026-06-17 13:45:00.41 INFO::Fitting model to feature number 601, ASV819
#> 2026-06-17 13:45:00.42 INFO::Fitting model to feature number 602, ASV821
#> 2026-06-17 13:45:00.42 INFO::Fitting model to feature number 603, ASV822
#> 2026-06-17 13:45:00.43 INFO::Fitting model to feature number 604, ASV823
#> 2026-06-17 13:45:00.43 INFO::Fitting model to feature number 605, ASV824
#> 2026-06-17 13:45:00.43 INFO::Fitting model to feature number 606, ASV828
#> 2026-06-17 13:45:00.44 INFO::Fitting model to feature number 607, ASV829
#> 2026-06-17 13:45:00.44 INFO::Fitting model to feature number 608, ASV830
#> 2026-06-17 13:45:00.44 INFO::Fitting model to feature number 609, ASV831
#> 2026-06-17 13:45:00.44 INFO::Fitting model to feature number 610, ASV832
#> 2026-06-17 13:45:00.45 INFO::Fitting model to feature number 611, ASV834
#> 2026-06-17 13:45:00.45 INFO::Fitting model to feature number 612, ASV836
#> 2026-06-17 13:45:00.46 INFO::Fitting model to feature number 613, ASV837
#> 2026-06-17 13:45:00.46 INFO::Fitting model to feature number 614, ASV838
#> 2026-06-17 13:45:00.46 INFO::Fitting model to feature number 615, ASV839
#> 2026-06-17 13:45:00.47 INFO::Fitting model to feature number 616, ASV840
#> 2026-06-17 13:45:00.47 INFO::Fitting model to feature number 617, ASV841
#> 2026-06-17 13:45:00.47 INFO::Fitting model to feature number 618, ASV842
#> 2026-06-17 13:45:00.48 INFO::Fitting model to feature number 619, ASV843
#> 2026-06-17 13:45:00.48 INFO::Fitting model to feature number 620, ASV844
#> 2026-06-17 13:45:00.48 INFO::Fitting model to feature number 621, ASV845
#> 2026-06-17 13:45:00.49 INFO::Fitting model to feature number 622, ASV847
#> 2026-06-17 13:45:00.49 WARNING::Fitting problem for feature 622 returning NA
#> 2026-06-17 13:45:00.49 INFO::Fitting model to feature number 623, ASV848
#> 2026-06-17 13:45:00.49 INFO::Fitting model to feature number 624, ASV852
#> 2026-06-17 13:45:00.50 INFO::Fitting model to feature number 625, ASV853
#> 2026-06-17 13:45:00.50 INFO::Fitting model to feature number 626, ASV854
#> 2026-06-17 13:45:00.50 INFO::Fitting model to feature number 627, ASV855
#> 2026-06-17 13:45:00.51 INFO::Fitting model to feature number 628, ASV857
#> 2026-06-17 13:45:00.51 INFO::Fitting model to feature number 629, ASV858
#> 2026-06-17 13:45:00.52 INFO::Fitting model to feature number 630, ASV859
#> 2026-06-17 13:45:00.52 INFO::Fitting model to feature number 631, ASV860
#> 2026-06-17 13:45:00.52 INFO::Fitting model to feature number 632, ASV861
#> 2026-06-17 13:45:00.53 WARNING::Fitting problem for feature 632 returning NA
#> 2026-06-17 13:45:00.53 INFO::Fitting model to feature number 633, ASV863
#> 2026-06-17 13:45:00.53 INFO::Fitting model to feature number 634, ASV865
#> 2026-06-17 13:45:00.54 INFO::Fitting model to feature number 635, ASV870
#> 2026-06-17 13:45:00.54 INFO::Fitting model to feature number 636, ASV873
#> 2026-06-17 13:45:00.54 INFO::Fitting model to feature number 637, ASV874
#> 2026-06-17 13:45:00.55 INFO::Fitting model to feature number 638, ASV875
#> 2026-06-17 13:45:00.55 INFO::Fitting model to feature number 639, ASV876
#> 2026-06-17 13:45:00.55 INFO::Fitting model to feature number 640, ASV877
#> 2026-06-17 13:45:00.56 INFO::Fitting model to feature number 641, ASV878
#> 2026-06-17 13:45:00.56 INFO::Fitting model to feature number 642, ASV879
#> 2026-06-17 13:45:00.56 INFO::Fitting model to feature number 643, ASV880
#> 2026-06-17 13:45:00.56 INFO::Fitting model to feature number 644, ASV883
#> 2026-06-17 13:45:00.57 INFO::Fitting model to feature number 645, ASV884
#> 2026-06-17 13:45:00.57 INFO::Fitting model to feature number 646, ASV887
#> 2026-06-17 13:45:00.58 INFO::Fitting model to feature number 647, ASV888
#> 2026-06-17 13:45:00.58 INFO::Fitting model to feature number 648, ASV890
#> 2026-06-17 13:45:00.58 INFO::Fitting model to feature number 649, ASV891
#> 2026-06-17 13:45:00.58 INFO::Fitting model to feature number 650, ASV892
#> 2026-06-17 13:45:00.59 INFO::Fitting model to feature number 651, ASV895
#> 2026-06-17 13:45:00.59 INFO::Fitting model to feature number 652, ASV897
#> 2026-06-17 13:45:00.60 INFO::Fitting model to feature number 653, ASV898
#> 2026-06-17 13:45:00.60 INFO::Fitting model to feature number 654, ASV899
#> 2026-06-17 13:45:00.60 INFO::Fitting model to feature number 655, ASV901
#> 2026-06-17 13:45:00.61 INFO::Fitting model to feature number 656, ASV902
#> 2026-06-17 13:45:00.61 INFO::Fitting model to feature number 657, ASV903
#> 2026-06-17 13:45:00.61 WARNING::Fitting problem for feature 657 returning NA
#> 2026-06-17 13:45:00.61 INFO::Fitting model to feature number 658, ASV904
#> 2026-06-17 13:45:00.62 INFO::Fitting model to feature number 659, ASV905
#> 2026-06-17 13:45:00.62 INFO::Fitting model to feature number 660, ASV906
#> 2026-06-17 13:45:00.62 INFO::Fitting model to feature number 661, ASV907
#> 2026-06-17 13:45:00.63 WARNING::Fitting problem for feature 661 returning NA
#> 2026-06-17 13:45:00.63 INFO::Fitting model to feature number 662, ASV909
#> 2026-06-17 13:45:00.63 INFO::Fitting model to feature number 663, ASV910
#> 2026-06-17 13:45:00.63 INFO::Fitting model to feature number 664, ASV911
#> 2026-06-17 13:45:00.64 INFO::Fitting model to feature number 665, ASV913
#> 2026-06-17 13:45:00.64 INFO::Fitting model to feature number 666, ASV914
#> 2026-06-17 13:45:00.65 INFO::Fitting model to feature number 667, ASV915
#> 2026-06-17 13:45:00.65 INFO::Fitting model to feature number 668, ASV916
#> 2026-06-17 13:45:00.65 INFO::Fitting model to feature number 669, ASV917
#> 2026-06-17 13:45:00.66 INFO::Fitting model to feature number 670, ASV918
#> 2026-06-17 13:45:00.66 INFO::Fitting model to feature number 671, ASV919
#> 2026-06-17 13:45:00.66 INFO::Fitting model to feature number 672, ASV920
#> 2026-06-17 13:45:00.67 INFO::Fitting model to feature number 673, ASV921
#> 2026-06-17 13:45:00.67 INFO::Fitting model to feature number 674, ASV922
#> 2026-06-17 13:45:00.67 INFO::Fitting model to feature number 675, ASV923
#> 2026-06-17 13:45:00.68 INFO::Fitting model to feature number 676, ASV926
#> 2026-06-17 13:45:00.68 INFO::Fitting model to feature number 677, ASV927
#> 2026-06-17 13:45:00.69 INFO::Fitting model to feature number 678, ASV928
#> 2026-06-17 13:45:00.69 INFO::Fitting model to feature number 679, ASV930
#> 2026-06-17 13:45:00.69 INFO::Fitting model to feature number 680, ASV932
#> 2026-06-17 13:45:00.70 INFO::Fitting model to feature number 681, ASV934
#> 2026-06-17 13:45:00.70 INFO::Fitting model to feature number 682, ASV935
#> 2026-06-17 13:45:00.70 INFO::Fitting model to feature number 683, ASV936
#> 2026-06-17 13:45:00.71 INFO::Fitting model to feature number 684, ASV937
#> 2026-06-17 13:45:00.71 INFO::Fitting model to feature number 685, ASV938
#> 2026-06-17 13:45:00.72 INFO::Fitting model to feature number 686, ASV939
#> 2026-06-17 13:45:00.72 INFO::Fitting model to feature number 687, ASV940
#> 2026-06-17 13:45:00.73 INFO::Fitting model to feature number 688, ASV941
#> 2026-06-17 13:45:00.73 INFO::Fitting model to feature number 689, ASV942
#> 2026-06-17 13:45:00.74 INFO::Fitting model to feature number 690, ASV943
#> 2026-06-17 13:45:00.74 INFO::Fitting model to feature number 691, ASV945
#> 2026-06-17 13:45:00.74 INFO::Fitting model to feature number 692, ASV947
#> 2026-06-17 13:45:00.75 INFO::Fitting model to feature number 693, ASV948
#> 2026-06-17 13:45:00.75 INFO::Fitting model to feature number 694, ASV949
#> 2026-06-17 13:45:00.75 INFO::Fitting model to feature number 695, ASV950
#> 2026-06-17 13:45:00.76 INFO::Fitting model to feature number 696, ASV951
#> 2026-06-17 13:45:00.76 INFO::Fitting model to feature number 697, ASV953
#> 2026-06-17 13:45:00.76 INFO::Fitting model to feature number 698, ASV955
#> 2026-06-17 13:45:00.76 INFO::Fitting model to feature number 699, ASV958
#> 2026-06-17 13:45:00.77 INFO::Fitting model to feature number 700, ASV959
#> 2026-06-17 13:45:00.77 INFO::Fitting model to feature number 701, ASV961
#> 2026-06-17 13:45:00.77 INFO::Fitting model to feature number 702, ASV962
#> 2026-06-17 13:45:00.78 INFO::Fitting model to feature number 703, ASV963
#> 2026-06-17 13:45:00.78 INFO::Fitting model to feature number 704, ASV964
#> 2026-06-17 13:45:00.78 INFO::Fitting model to feature number 705, ASV966
#> 2026-06-17 13:45:00.79 INFO::Fitting model to feature number 706, ASV967
#> 2026-06-17 13:45:00.79 INFO::Fitting model to feature number 707, ASV969
#> 2026-06-17 13:45:00.80 INFO::Fitting model to feature number 708, ASV970
#> 2026-06-17 13:45:00.80 INFO::Fitting model to feature number 709, ASV971
#> 2026-06-17 13:45:00.80 INFO::Fitting model to feature number 710, ASV972
#> 2026-06-17 13:45:00.80 INFO::Fitting model to feature number 711, ASV973
#> 2026-06-17 13:45:00.80 INFO::Fitting model to feature number 712, ASV974
#> 2026-06-17 13:45:00.81 INFO::Fitting model to feature number 713, ASV975
#> 2026-06-17 13:45:00.81 INFO::Fitting model to feature number 714, ASV976
#> 2026-06-17 13:45:00.82 INFO::Fitting model to feature number 715, ASV977
#> 2026-06-17 13:45:00.82 INFO::Fitting model to feature number 716, ASV979
#> 2026-06-17 13:45:00.82 INFO::Fitting model to feature number 717, ASV980
#> 2026-06-17 13:45:00.82 INFO::Fitting model to feature number 718, ASV981
#> 2026-06-17 13:45:00.83 INFO::Fitting model to feature number 719, ASV983
#> 2026-06-17 13:45:00.83 INFO::Fitting model to feature number 720, ASV984
#> 2026-06-17 13:45:00.84 INFO::Fitting model to feature number 721, ASV986
#> 2026-06-17 13:45:00.84 INFO::Fitting model to feature number 722, ASV987
#> 2026-06-17 13:45:00.84 INFO::Fitting model to feature number 723, ASV988
#> 2026-06-17 13:45:00.85 INFO::Fitting model to feature number 724, ASV989
#> 2026-06-17 13:45:00.85 INFO::Fitting model to feature number 725, ASV990
#> 2026-06-17 13:45:00.85 INFO::Fitting model to feature number 726, ASV992
#> 2026-06-17 13:45:00.86 WARNING::Fitting problem for feature 726 returning NA
#> 2026-06-17 13:45:00.86 INFO::Fitting model to feature number 727, ASV993
#> 2026-06-17 13:45:00.86 INFO::Fitting model to feature number 728, ASV994
#> 2026-06-17 13:45:00.86 INFO::Fitting model to feature number 729, ASV995
#> 2026-06-17 13:45:00.87 INFO::Fitting model to feature number 730, ASV996
#> 2026-06-17 13:45:00.87 INFO::Fitting model to feature number 731, ASV997
#> 2026-06-17 13:45:00.87 INFO::Fitting model to feature number 732, ASV998
#> 2026-06-17 13:45:00.87 INFO::Fitting model to feature number 733, ASV999
#> 2026-06-17 13:45:00.88 INFO::Fitting model to feature number 734, ASV1001
#> 2026-06-17 13:45:00.88 INFO::Fitting model to feature number 735, ASV1002
#> 2026-06-17 13:45:00.89 INFO::Fitting model to feature number 736, ASV1003
#> 2026-06-17 13:45:00.89 INFO::Fitting model to feature number 737, ASV1004
#> 2026-06-17 13:45:00.89 INFO::Fitting model to feature number 738, ASV1005
#> 2026-06-17 13:45:00.90 INFO::Fitting model to feature number 739, ASV1006
#> 2026-06-17 13:45:00.90 INFO::Fitting model to feature number 740, ASV1007
#> 2026-06-17 13:45:00.91 INFO::Fitting model to feature number 741, ASV1008
#> 2026-06-17 13:45:00.91 INFO::Fitting model to feature number 742, ASV1010
#> 2026-06-17 13:45:00.91 INFO::Fitting model to feature number 743, ASV1011
#> 2026-06-17 13:45:00.97 INFO::Fitting model to feature number 744, ASV1014
#> 2026-06-17 13:45:00.97 INFO::Fitting model to feature number 745, ASV1015
#> 2026-06-17 13:45:00.97 WARNING::Fitting problem for feature 745 returning NA
#> 2026-06-17 13:45:00.98 INFO::Fitting model to feature number 746, ASV1016
#> 2026-06-17 13:45:00.98 INFO::Fitting model to feature number 747, ASV1017
#> 2026-06-17 13:45:00.98 INFO::Fitting model to feature number 748, ASV1018
#> 2026-06-17 13:45:00.99 INFO::Fitting model to feature number 749, ASV1020
#> 2026-06-17 13:45:00.99 INFO::Fitting model to feature number 750, ASV1021
#> 2026-06-17 13:45:00.99 INFO::Fitting model to feature number 751, ASV1022
#> 2026-06-17 13:45:01.00 INFO::Fitting model to feature number 752, ASV1023
#> 2026-06-17 13:45:01.00 INFO::Fitting model to feature number 753, ASV1024
#> 2026-06-17 13:45:01.00 INFO::Fitting model to feature number 754, ASV1025
#> 2026-06-17 13:45:01.00 INFO::Fitting model to feature number 755, ASV1026
#> 2026-06-17 13:45:01.00 INFO::Fitting model to feature number 756, ASV1027
#> 2026-06-17 13:45:01.01 INFO::Fitting model to feature number 757, ASV1028
#> 2026-06-17 13:45:01.01 INFO::Fitting model to feature number 758, ASV1029
#> 2026-06-17 13:45:01.01 INFO::Fitting model to feature number 759, ASV1030
#> 2026-06-17 13:45:01.02 INFO::Fitting model to feature number 760, ASV1031
#> 2026-06-17 13:45:01.02 INFO::Fitting model to feature number 761, ASV1032
#> 2026-06-17 13:45:01.03 INFO::Fitting model to feature number 762, ASV1035
#> 2026-06-17 13:45:01.03 INFO::Fitting model to feature number 763, ASV1036
#> 2026-06-17 13:45:01.03 INFO::Fitting model to feature number 764, ASV1037
#> 2026-06-17 13:45:01.04 INFO::Fitting model to feature number 765, ASV1038
#> 2026-06-17 13:45:01.04 INFO::Fitting model to feature number 766, ASV1039
#> 2026-06-17 13:45:01.04 INFO::Fitting model to feature number 767, ASV1040
#> 2026-06-17 13:45:01.04 INFO::Fitting model to feature number 768, ASV1041
#> 2026-06-17 13:45:01.05 INFO::Fitting model to feature number 769, ASV1042
#> 2026-06-17 13:45:01.05 INFO::Fitting model to feature number 770, ASV1044
#> 2026-06-17 13:45:01.05 INFO::Fitting model to feature number 771, ASV1046
#> 2026-06-17 13:45:01.06 INFO::Fitting model to feature number 772, ASV1047
#> 2026-06-17 13:45:01.06 INFO::Fitting model to feature number 773, ASV1048
#> 2026-06-17 13:45:01.06 INFO::Fitting model to feature number 774, ASV1049
#> 2026-06-17 13:45:01.07 INFO::Fitting model to feature number 775, ASV1052
#> 2026-06-17 13:45:01.07 INFO::Fitting model to feature number 776, ASV1053
#> 2026-06-17 13:45:01.08 INFO::Fitting model to feature number 777, ASV1057
#> 2026-06-17 13:45:01.08 INFO::Fitting model to feature number 778, ASV1058
#> 2026-06-17 13:45:01.08 INFO::Fitting model to feature number 779, ASV1059
#> 2026-06-17 13:45:01.09 INFO::Fitting model to feature number 780, ASV1063
#> 2026-06-17 13:45:01.09 INFO::Fitting model to feature number 781, ASV1065
#> 2026-06-17 13:45:01.10 INFO::Fitting model to feature number 782, ASV1066
#> 2026-06-17 13:45:01.10 INFO::Fitting model to feature number 783, ASV1067
#> 2026-06-17 13:45:01.10 INFO::Fitting model to feature number 784, ASV1068
#> 2026-06-17 13:45:01.10 INFO::Fitting model to feature number 785, ASV1069
#> 2026-06-17 13:45:01.11 INFO::Fitting model to feature number 786, ASV1070
#> 2026-06-17 13:45:01.11 INFO::Fitting model to feature number 787, ASV1071
#> 2026-06-17 13:45:01.11 INFO::Fitting model to feature number 788, ASV1072
#> 2026-06-17 13:45:01.11 INFO::Fitting model to feature number 789, ASV1073
#> 2026-06-17 13:45:01.12 INFO::Fitting model to feature number 790, ASV1074
#> 2026-06-17 13:45:01.12 INFO::Fitting model to feature number 791, ASV1076
#> 2026-06-17 13:45:01.13 INFO::Fitting model to feature number 792, ASV1078
#> 2026-06-17 13:45:01.13 INFO::Fitting model to feature number 793, ASV1079
#> 2026-06-17 13:45:01.13 INFO::Fitting model to feature number 794, ASV1082
#> 2026-06-17 13:45:01.13 INFO::Fitting model to feature number 795, ASV1084
#> 2026-06-17 13:45:01.14 INFO::Fitting model to feature number 796, ASV1085
#> 2026-06-17 13:45:01.14 INFO::Fitting model to feature number 797, ASV1086
#> 2026-06-17 13:45:01.14 INFO::Fitting model to feature number 798, ASV1087
#> 2026-06-17 13:45:01.15 INFO::Fitting model to feature number 799, ASV1090
#> 2026-06-17 13:45:01.15 INFO::Fitting model to feature number 800, ASV1091
#> 2026-06-17 13:45:01.16 WARNING::Fitting problem for feature 800 returning NA
#> 2026-06-17 13:45:01.16 INFO::Fitting model to feature number 801, ASV1093
#> 2026-06-17 13:45:01.16 INFO::Fitting model to feature number 802, ASV1095
#> 2026-06-17 13:45:01.17 INFO::Fitting model to feature number 803, ASV1096
#> 2026-06-17 13:45:01.17 INFO::Fitting model to feature number 804, ASV1097
#> 2026-06-17 13:45:01.17 INFO::Fitting model to feature number 805, ASV1099
#> 2026-06-17 13:45:01.17 INFO::Fitting model to feature number 806, ASV1100
#> 2026-06-17 13:45:01.18 INFO::Fitting model to feature number 807, ASV1101
#> 2026-06-17 13:45:01.18 INFO::Fitting model to feature number 808, ASV1103
#> 2026-06-17 13:45:01.19 INFO::Fitting model to feature number 809, ASV1105
#> 2026-06-17 13:45:01.19 INFO::Fitting model to feature number 810, ASV1107
#> 2026-06-17 13:45:01.19 INFO::Fitting model to feature number 811, ASV1108
#> 2026-06-17 13:45:01.19 INFO::Fitting model to feature number 812, ASV1109
#> 2026-06-17 13:45:01.20 INFO::Fitting model to feature number 813, ASV1110
#> 2026-06-17 13:45:01.20 INFO::Fitting model to feature number 814, ASV1111
#> 2026-06-17 13:45:01.21 INFO::Fitting model to feature number 815, ASV1112
#> 2026-06-17 13:45:01.21 INFO::Fitting model to feature number 816, ASV1114
#> 2026-06-17 13:45:01.21 INFO::Fitting model to feature number 817, ASV1115
#> 2026-06-17 13:45:01.22 INFO::Fitting model to feature number 818, ASV1116
#> 2026-06-17 13:45:01.22 INFO::Fitting model to feature number 819, ASV1120
#> 2026-06-17 13:45:01.22 INFO::Fitting model to feature number 820, ASV1121
#> 2026-06-17 13:45:01.23 INFO::Fitting model to feature number 821, ASV1122
#> 2026-06-17 13:45:01.23 INFO::Fitting model to feature number 822, ASV1125
#> 2026-06-17 13:45:01.24 INFO::Fitting model to feature number 823, ASV1126
#> 2026-06-17 13:45:01.24 INFO::Fitting model to feature number 824, ASV1127
#> 2026-06-17 13:45:01.24 INFO::Fitting model to feature number 825, ASV1128
#> 2026-06-17 13:45:01.25 INFO::Fitting model to feature number 826, ASV1129
#> 2026-06-17 13:45:01.25 INFO::Fitting model to feature number 827, ASV1131
#> 2026-06-17 13:45:01.25 INFO::Fitting model to feature number 828, ASV1132
#> 2026-06-17 13:45:01.25 WARNING::Fitting problem for feature 828 returning NA
#> 2026-06-17 13:45:01.25 INFO::Fitting model to feature number 829, ASV1133
#> 2026-06-17 13:45:01.26 INFO::Fitting model to feature number 830, ASV1134
#> 2026-06-17 13:45:01.26 WARNING::Fitting problem for feature 830 returning NA
#> 2026-06-17 13:45:01.26 INFO::Fitting model to feature number 831, ASV1135
#> 2026-06-17 13:45:01.26 INFO::Fitting model to feature number 832, ASV1137
#> 2026-06-17 13:45:01.27 INFO::Fitting model to feature number 833, ASV1138
#> 2026-06-17 13:45:01.27 INFO::Fitting model to feature number 834, ASV1139
#> 2026-06-17 13:45:01.28 INFO::Fitting model to feature number 835, ASV1141
#> 2026-06-17 13:45:01.28 INFO::Fitting model to feature number 836, ASV1143
#> 2026-06-17 13:45:01.28 INFO::Fitting model to feature number 837, ASV1144
#> 2026-06-17 13:45:01.28 INFO::Fitting model to feature number 838, ASV1146
#> 2026-06-17 13:45:01.29 INFO::Fitting model to feature number 839, ASV1147
#> 2026-06-17 13:45:01.29 INFO::Fitting model to feature number 840, ASV1148
#> 2026-06-17 13:45:01.29 INFO::Fitting model to feature number 841, ASV1150
#> 2026-06-17 13:45:01.30 INFO::Fitting model to feature number 842, ASV1151
#> 2026-06-17 13:45:01.30 INFO::Fitting model to feature number 843, ASV1152
#> 2026-06-17 13:45:01.30 INFO::Fitting model to feature number 844, ASV1154
#> 2026-06-17 13:45:01.31 INFO::Fitting model to feature number 845, ASV1155
#> 2026-06-17 13:45:01.31 INFO::Fitting model to feature number 846, ASV1156
#> 2026-06-17 13:45:01.31 INFO::Fitting model to feature number 847, ASV1158
#> 2026-06-17 13:45:01.32 INFO::Fitting model to feature number 848, ASV1159
#> 2026-06-17 13:45:01.32 WARNING::Fitting problem for feature 848 returning NA
#> 2026-06-17 13:45:01.32 INFO::Fitting model to feature number 849, ASV1160
#> 2026-06-17 13:45:01.32 INFO::Fitting model to feature number 850, ASV1161
#> 2026-06-17 13:45:01.33 INFO::Fitting model to feature number 851, ASV1162
#> 2026-06-17 13:45:01.33 INFO::Fitting model to feature number 852, ASV1163
#> 2026-06-17 13:45:01.34 INFO::Fitting model to feature number 853, ASV1164
#> 2026-06-17 13:45:01.34 INFO::Fitting model to feature number 854, ASV1165
#> 2026-06-17 13:45:01.34 WARNING::Fitting problem for feature 854 returning NA
#> 2026-06-17 13:45:01.34 INFO::Fitting model to feature number 855, ASV1167
#> 2026-06-17 13:45:01.35 INFO::Fitting model to feature number 856, ASV1168
#> 2026-06-17 13:45:01.35 INFO::Fitting model to feature number 857, ASV1169
#> 2026-06-17 13:45:01.35 WARNING::Fitting problem for feature 857 returning NA
#> 2026-06-17 13:45:01.35 INFO::Fitting model to feature number 858, ASV1171
#> 2026-06-17 13:45:01.36 INFO::Fitting model to feature number 859, ASV1172
#> 2026-06-17 13:45:01.36 INFO::Fitting model to feature number 860, ASV1173
#> 2026-06-17 13:45:01.37 INFO::Fitting model to feature number 861, ASV1175
#> 2026-06-17 13:45:01.37 INFO::Fitting model to feature number 862, ASV1176
#> 2026-06-17 13:45:01.37 INFO::Fitting model to feature number 863, ASV1177
#> 2026-06-17 13:45:01.38 INFO::Fitting model to feature number 864, ASV1179
#> 2026-06-17 13:45:01.38 WARNING::Fitting problem for feature 864 returning NA
#> 2026-06-17 13:45:01.38 INFO::Fitting model to feature number 865, ASV1180
#> 2026-06-17 13:45:01.39 INFO::Fitting model to feature number 866, ASV1182
#> 2026-06-17 13:45:01.39 INFO::Fitting model to feature number 867, ASV1184
#> 2026-06-17 13:45:01.39 INFO::Fitting model to feature number 868, ASV1185
#> 2026-06-17 13:45:01.39 INFO::Fitting model to feature number 869, ASV1186
#> 2026-06-17 13:45:01.40 INFO::Fitting model to feature number 870, ASV1187
#> 2026-06-17 13:45:01.40 INFO::Fitting model to feature number 871, ASV1189
#> 2026-06-17 13:45:01.41 INFO::Fitting model to feature number 872, ASV1190
#> 2026-06-17 13:45:01.41 INFO::Fitting model to feature number 873, ASV1192
#> 2026-06-17 13:45:01.41 INFO::Fitting model to feature number 874, ASV1193
#> 2026-06-17 13:45:01.42 INFO::Fitting model to feature number 875, ASV1194
#> 2026-06-17 13:45:01.42 INFO::Fitting model to feature number 876, ASV1195
#> 2026-06-17 13:45:01.42 INFO::Fitting model to feature number 877, ASV1198
#> 2026-06-17 13:45:01.43 INFO::Fitting model to feature number 878, ASV1199
#> 2026-06-17 13:45:01.43 INFO::Fitting model to feature number 879, ASV1200
#> 2026-06-17 13:45:01.43 INFO::Fitting model to feature number 880, ASV1203
#> 2026-06-17 13:45:01.44 INFO::Fitting model to feature number 881, ASV1204
#> 2026-06-17 13:45:01.44 INFO::Fitting model to feature number 882, ASV1205
#> 2026-06-17 13:45:01.44 INFO::Fitting model to feature number 883, ASV1206
#> 2026-06-17 13:45:01.45 INFO::Fitting model to feature number 884, ASV1208
#> 2026-06-17 13:45:01.45 INFO::Fitting model to feature number 885, ASV1209
#> 2026-06-17 13:45:01.45 INFO::Fitting model to feature number 886, ASV1210
#> 2026-06-17 13:45:01.46 INFO::Fitting model to feature number 887, ASV1211
#> 2026-06-17 13:45:01.46 INFO::Fitting model to feature number 888, ASV1212
#> 2026-06-17 13:45:01.46 INFO::Fitting model to feature number 889, ASV1213
#> 2026-06-17 13:45:01.47 INFO::Fitting model to feature number 890, ASV1214
#> 2026-06-17 13:45:01.47 INFO::Fitting model to feature number 891, ASV1216
#> 2026-06-17 13:45:01.47 INFO::Fitting model to feature number 892, ASV1217
#> 2026-06-17 13:45:01.48 INFO::Fitting model to feature number 893, ASV1218
#> 2026-06-17 13:45:01.48 INFO::Fitting model to feature number 894, ASV1219
#> 2026-06-17 13:45:01.48 INFO::Fitting model to feature number 895, ASV1221
#> 2026-06-17 13:45:01.49 WARNING::Fitting problem for feature 895 returning NA
#> 2026-06-17 13:45:01.49 INFO::Fitting model to feature number 896, ASV1223
#> 2026-06-17 13:45:01.49 INFO::Fitting model to feature number 897, ASV1224
#> 2026-06-17 13:45:01.49 INFO::Fitting model to feature number 898, ASV1225
#> 2026-06-17 13:45:01.50 INFO::Fitting model to feature number 899, ASV1227
#> 2026-06-17 13:45:01.50 INFO::Fitting model to feature number 900, ASV1228
#> 2026-06-17 13:45:01.50 INFO::Fitting model to feature number 901, ASV1229
#> 2026-06-17 13:45:01.50 WARNING::Fitting problem for feature 901 returning NA
#> 2026-06-17 13:45:01.50 INFO::Fitting model to feature number 902, ASV1230
#> 2026-06-17 13:45:01.51 INFO::Fitting model to feature number 903, ASV1231
#> 2026-06-17 13:45:01.51 INFO::Fitting model to feature number 904, ASV1232
#> 2026-06-17 13:45:01.52 INFO::Fitting model to feature number 905, ASV1233
#> 2026-06-17 13:45:01.52 INFO::Fitting model to feature number 906, ASV1234
#> 2026-06-17 13:45:01.52 WARNING::Fitting problem for feature 906 returning NA
#> 2026-06-17 13:45:01.52 INFO::Fitting model to feature number 907, ASV1236
#> 2026-06-17 13:45:01.53 INFO::Fitting model to feature number 908, ASV1238
#> 2026-06-17 13:45:01.53 INFO::Fitting model to feature number 909, ASV1239
#> 2026-06-17 13:45:01.54 INFO::Fitting model to feature number 910, ASV1241
#> 2026-06-17 13:45:01.54 INFO::Fitting model to feature number 911, ASV1242
#> 2026-06-17 13:45:01.54 INFO::Fitting model to feature number 912, ASV1243
#> 2026-06-17 13:45:01.55 INFO::Fitting model to feature number 913, ASV1245
#> 2026-06-17 13:45:01.55 INFO::Fitting model to feature number 914, ASV1246
#> 2026-06-17 13:45:01.55 INFO::Fitting model to feature number 915, ASV1247
#> 2026-06-17 13:45:01.55 INFO::Fitting model to feature number 916, ASV1251
#> 2026-06-17 13:45:01.56 INFO::Fitting model to feature number 917, ASV1252
#> 2026-06-17 13:45:01.56 INFO::Fitting model to feature number 918, ASV1253
#> 2026-06-17 13:45:01.57 INFO::Fitting model to feature number 919, ASV1254
#> 2026-06-17 13:45:01.57 INFO::Fitting model to feature number 920, ASV1257
#> 2026-06-17 13:45:01.57 INFO::Fitting model to feature number 921, ASV1258
#> 2026-06-17 13:45:01.58 INFO::Fitting model to feature number 922, ASV1259
#> 2026-06-17 13:45:01.58 INFO::Fitting model to feature number 923, ASV1260
#> 2026-06-17 13:45:01.58 INFO::Fitting model to feature number 924, ASV1261
#> 2026-06-17 13:45:01.59 INFO::Fitting model to feature number 925, ASV1262
#> 2026-06-17 13:45:01.59 INFO::Fitting model to feature number 926, ASV1263
#> 2026-06-17 13:45:01.59 INFO::Fitting model to feature number 927, ASV1264
#> 2026-06-17 13:45:01.60 INFO::Fitting model to feature number 928, ASV1265
#> 2026-06-17 13:45:01.60 INFO::Fitting model to feature number 929, ASV1267
#> 2026-06-17 13:45:01.60 INFO::Fitting model to feature number 930, ASV1268
#> 2026-06-17 13:45:01.60 INFO::Fitting model to feature number 931, ASV1269
#> 2026-06-17 13:45:01.60 INFO::Fitting model to feature number 932, ASV1270
#> 2026-06-17 13:45:01.61 INFO::Fitting model to feature number 933, ASV1271
#> 2026-06-17 13:45:01.61 INFO::Fitting model to feature number 934, ASV1272
#> 2026-06-17 13:45:01.61 INFO::Fitting model to feature number 935, ASV1273
#> 2026-06-17 13:45:01.62 INFO::Fitting model to feature number 936, ASV1274
#> 2026-06-17 13:45:01.62 INFO::Fitting model to feature number 937, ASV1275
#> 2026-06-17 13:45:01.62 INFO::Fitting model to feature number 938, ASV1276
#> 2026-06-17 13:45:01.62 INFO::Fitting model to feature number 939, ASV1278
#> 2026-06-17 13:45:01.63 INFO::Fitting model to feature number 940, ASV1279
#> 2026-06-17 13:45:01.63 INFO::Fitting model to feature number 941, ASV1282
#> 2026-06-17 13:45:01.63 INFO::Fitting model to feature number 942, ASV1283
#> 2026-06-17 13:45:01.64 INFO::Fitting model to feature number 943, ASV1284
#> 2026-06-17 13:45:01.64 INFO::Fitting model to feature number 944, ASV1285
#> 2026-06-17 13:45:01.64 INFO::Fitting model to feature number 945, ASV1286
#> 2026-06-17 13:45:01.64 INFO::Fitting model to feature number 946, ASV1287
#> 2026-06-17 13:45:01.65 INFO::Fitting model to feature number 947, ASV1288
#> 2026-06-17 13:45:01.65 INFO::Fitting model to feature number 948, ASV1289
#> 2026-06-17 13:45:01.66 INFO::Fitting model to feature number 949, ASV1290
#> 2026-06-17 13:45:01.66 INFO::Fitting model to feature number 950, ASV1293
#> 2026-06-17 13:45:01.66 INFO::Fitting model to feature number 951, ASV1294
#> 2026-06-17 13:45:01.66 INFO::Fitting model to feature number 952, ASV1296
#> 2026-06-17 13:45:01.67 INFO::Fitting model to feature number 953, ASV1297
#> 2026-06-17 13:45:01.67 INFO::Fitting model to feature number 954, ASV1300
#> 2026-06-17 13:45:01.68 INFO::Fitting model to feature number 955, ASV1301
#> 2026-06-17 13:45:01.68 INFO::Fitting model to feature number 956, ASV1302
#> 2026-06-17 13:45:01.68 INFO::Fitting model to feature number 957, ASV1303
#> 2026-06-17 13:45:01.69 INFO::Fitting model to feature number 958, ASV1304
#> 2026-06-17 13:45:01.69 INFO::Fitting model to feature number 959, ASV1305
#> 2026-06-17 13:45:01.69 INFO::Fitting model to feature number 960, ASV1307
#> 2026-06-17 13:45:01.70 INFO::Fitting model to feature number 961, ASV1310
#> 2026-06-17 13:45:01.70 INFO::Fitting model to feature number 962, ASV1311
#> 2026-06-17 13:45:01.70 INFO::Fitting model to feature number 963, ASV1312
#> 2026-06-17 13:45:01.71 INFO::Fitting model to feature number 964, ASV1313
#> 2026-06-17 13:45:01.71 INFO::Fitting model to feature number 965, ASV1314
#> 2026-06-17 13:45:01.71 INFO::Fitting model to feature number 966, ASV1315
#> 2026-06-17 13:45:01.72 INFO::Fitting model to feature number 967, ASV1316
#> 2026-06-17 13:45:01.72 INFO::Fitting model to feature number 968, ASV1317
#> 2026-06-17 13:45:01.73 INFO::Fitting model to feature number 969, ASV1319
#> 2026-06-17 13:45:01.73 INFO::Fitting model to feature number 970, ASV1320
#> 2026-06-17 13:45:01.73 INFO::Fitting model to feature number 971, ASV1321
#> 2026-06-17 13:45:01.74 INFO::Fitting model to feature number 972, ASV1323
#> 2026-06-17 13:45:01.74 INFO::Fitting model to feature number 973, ASV1326
#> 2026-06-17 13:45:01.74 INFO::Fitting model to feature number 974, ASV1327
#> 2026-06-17 13:45:01.75 INFO::Fitting model to feature number 975, ASV1328
#> 2026-06-17 13:45:01.75 INFO::Fitting model to feature number 976, ASV1330
#> 2026-06-17 13:45:01.75 INFO::Fitting model to feature number 977, ASV1332
#> 2026-06-17 13:45:01.76 INFO::Fitting model to feature number 978, ASV1334
#> 2026-06-17 13:45:01.76 INFO::Fitting model to feature number 979, ASV1335
#> 2026-06-17 13:45:01.76 INFO::Fitting model to feature number 980, ASV1336
#> 2026-06-17 13:45:01.77 INFO::Fitting model to feature number 981, ASV1337
#> 2026-06-17 13:45:01.77 INFO::Fitting model to feature number 982, ASV1338
#> 2026-06-17 13:45:01.77 INFO::Fitting model to feature number 983, ASV1340
#> 2026-06-17 13:45:01.78 INFO::Fitting model to feature number 984, ASV1341
#> 2026-06-17 13:45:01.78 INFO::Fitting model to feature number 985, ASV1342
#> 2026-06-17 13:45:01.78 INFO::Fitting model to feature number 986, ASV1345
#> 2026-06-17 13:45:01.79 INFO::Fitting model to feature number 987, ASV1350
#> 2026-06-17 13:45:01.79 INFO::Fitting model to feature number 988, ASV1351
#> 2026-06-17 13:45:01.80 INFO::Fitting model to feature number 989, ASV1352
#> 2026-06-17 13:45:01.80 INFO::Fitting model to feature number 990, ASV1353
#> 2026-06-17 13:45:01.80 INFO::Fitting model to feature number 991, ASV1355
#> 2026-06-17 13:45:01.81 INFO::Fitting model to feature number 992, ASV1356
#> 2026-06-17 13:45:01.81 INFO::Fitting model to feature number 993, ASV1359
#> 2026-06-17 13:45:01.81 INFO::Fitting model to feature number 994, ASV1360
#> 2026-06-17 13:45:01.82 INFO::Fitting model to feature number 995, ASV1361
#> 2026-06-17 13:45:01.82 INFO::Fitting model to feature number 996, ASV1362
#> 2026-06-17 13:45:01.83 INFO::Fitting model to feature number 997, ASV1363
#> 2026-06-17 13:45:01.83 INFO::Fitting model to feature number 998, ASV1365
#> 2026-06-17 13:45:01.83 INFO::Fitting model to feature number 999, ASV1366
#> 2026-06-17 13:45:01.84 INFO::Fitting model to feature number 1000, ASV1367
#> 2026-06-17 13:45:01.84 INFO::Fitting model to feature number 1001, ASV1368
#> 2026-06-17 13:45:01.84 INFO::Fitting model to feature number 1002, ASV1369
#> 2026-06-17 13:45:01.85 INFO::Fitting model to feature number 1003, ASV1370
#> 2026-06-17 13:45:01.85 INFO::Fitting model to feature number 1004, ASV1371
#> 2026-06-17 13:45:01.85 WARNING::Fitting problem for feature 1004 returning NA
#> 2026-06-17 13:45:01.85 INFO::Fitting model to feature number 1005, ASV1372
#> 2026-06-17 13:45:01.86 INFO::Fitting model to feature number 1006, ASV1373
#> 2026-06-17 13:45:01.86 INFO::Fitting model to feature number 1007, ASV1374
#> 2026-06-17 13:45:01.86 INFO::Fitting model to feature number 1008, ASV1375
#> 2026-06-17 13:45:01.86 INFO::Fitting model to feature number 1009, ASV1376
#> 2026-06-17 13:45:01.87 INFO::Fitting model to feature number 1010, ASV1378
#> 2026-06-17 13:45:01.87 INFO::Fitting model to feature number 1011, ASV1379
#> 2026-06-17 13:45:01.87 INFO::Fitting model to feature number 1012, ASV1380
#> 2026-06-17 13:45:01.88 INFO::Fitting model to feature number 1013, ASV1381
#> 2026-06-17 13:45:01.88 INFO::Fitting model to feature number 1014, ASV1384
#> 2026-06-17 13:45:01.88 INFO::Fitting model to feature number 1015, ASV1386
#> 2026-06-17 13:45:01.88 INFO::Fitting model to feature number 1016, ASV1387
#> 2026-06-17 13:45:01.88 INFO::Fitting model to feature number 1017, ASV1388
#> 2026-06-17 13:45:01.89 WARNING::Fitting problem for feature 1017 returning NA
#> 2026-06-17 13:45:01.89 INFO::Fitting model to feature number 1018, ASV1389
#> 2026-06-17 13:45:01.89 INFO::Fitting model to feature number 1019, ASV1390
#> 2026-06-17 13:45:01.89 INFO::Fitting model to feature number 1020, ASV1392
#> 2026-06-17 13:45:01.90 INFO::Fitting model to feature number 1021, ASV1393
#> 2026-06-17 13:45:01.90 INFO::Fitting model to feature number 1022, ASV1396
#> 2026-06-17 13:45:01.90 INFO::Fitting model to feature number 1023, ASV1397
#> 2026-06-17 13:45:01.90 INFO::Fitting model to feature number 1024, ASV1398
#> 2026-06-17 13:45:01.91 INFO::Fitting model to feature number 1025, ASV1399
#> 2026-06-17 13:45:01.91 WARNING::Fitting problem for feature 1025 returning NA
#> 2026-06-17 13:45:01.91 INFO::Fitting model to feature number 1026, ASV1400
#> 2026-06-17 13:45:01.91 INFO::Fitting model to feature number 1027, ASV1401
#> 2026-06-17 13:45:01.92 INFO::Fitting model to feature number 1028, ASV1403
#> 2026-06-17 13:45:01.92 INFO::Fitting model to feature number 1029, ASV1404
#> 2026-06-17 13:45:01.92 INFO::Fitting model to feature number 1030, ASV1406
#> 2026-06-17 13:45:01.93 INFO::Fitting model to feature number 1031, ASV1408
#> 2026-06-17 13:45:01.93 INFO::Fitting model to feature number 1032, ASV1409
#> 2026-06-17 13:45:01.93 INFO::Fitting model to feature number 1033, ASV1410
#> 2026-06-17 13:45:01.94 INFO::Fitting model to feature number 1034, ASV1412
#> 2026-06-17 13:45:01.94 INFO::Fitting model to feature number 1035, ASV1413
#> 2026-06-17 13:45:01.94 INFO::Fitting model to feature number 1036, ASV1414
#> 2026-06-17 13:45:01.94 INFO::Fitting model to feature number 1037, ASV1415
#> 2026-06-17 13:45:01.95 WARNING::Fitting problem for feature 1037 returning NA
#> 2026-06-17 13:45:01.95 INFO::Fitting model to feature number 1038, ASV1418
#> 2026-06-17 13:45:01.95 INFO::Fitting model to feature number 1039, ASV1419
#> 2026-06-17 13:45:01.95 INFO::Fitting model to feature number 1040, ASV1420
#> 2026-06-17 13:45:01.96 INFO::Fitting model to feature number 1041, ASV1421
#> 2026-06-17 13:45:01.96 INFO::Fitting model to feature number 1042, ASV1422
#> 2026-06-17 13:45:01.96 INFO::Fitting model to feature number 1043, ASV1423
#> 2026-06-17 13:45:01.97 WARNING::Fitting problem for feature 1043 returning NA
#> 2026-06-17 13:45:01.97 INFO::Fitting model to feature number 1044, ASV1424
#> 2026-06-17 13:45:01.97 INFO::Fitting model to feature number 1045, ASV1425
#> 2026-06-17 13:45:01.97 INFO::Fitting model to feature number 1046, ASV1427
#> 2026-06-17 13:45:01.98 INFO::Fitting model to feature number 1047, ASV1428
#> 2026-06-17 13:45:01.98 INFO::Fitting model to feature number 1048, ASV1430
#> 2026-06-17 13:45:01.98 INFO::Fitting model to feature number 1049, ASV1433
#> 2026-06-17 13:45:01.98 INFO::Fitting model to feature number 1050, ASV1434
#> 2026-06-17 13:45:01.99 INFO::Fitting model to feature number 1051, ASV1435
#> 2026-06-17 13:45:01.99 INFO::Fitting model to feature number 1052, ASV1436
#> 2026-06-17 13:45:01.99 WARNING::Fitting problem for feature 1052 returning NA
#> 2026-06-17 13:45:02.00 INFO::Fitting model to feature number 1053, ASV1437
#> 2026-06-17 13:45:02.00 INFO::Fitting model to feature number 1054, ASV1438
#> 2026-06-17 13:45:02.00 INFO::Fitting model to feature number 1055, ASV1442
#> 2026-06-17 13:45:02.00 INFO::Fitting model to feature number 1056, ASV1443
#> 2026-06-17 13:45:02.01 INFO::Fitting model to feature number 1057, ASV1449
#> 2026-06-17 13:45:02.01 INFO::Fitting model to feature number 1058, ASV1450
#> 2026-06-17 13:45:02.01 INFO::Fitting model to feature number 1059, ASV1452
#> 2026-06-17 13:45:02.01 INFO::Fitting model to feature number 1060, ASV1454
#> 2026-06-17 13:45:02.02 INFO::Fitting model to feature number 1061, ASV1455
#> 2026-06-17 13:45:02.02 INFO::Fitting model to feature number 1062, ASV1458
#> 2026-06-17 13:45:02.02 INFO::Fitting model to feature number 1063, ASV1459
#> 2026-06-17 13:45:02.02 INFO::Fitting model to feature number 1064, ASV1460
#> 2026-06-17 13:45:02.03 INFO::Fitting model to feature number 1065, ASV1461
#> 2026-06-17 13:45:02.03 INFO::Fitting model to feature number 1066, ASV1462
#> 2026-06-17 13:45:02.04 INFO::Fitting model to feature number 1067, ASV1463
#> 2026-06-17 13:45:02.04 INFO::Fitting model to feature number 1068, ASV1466
#> 2026-06-17 13:45:02.04 INFO::Fitting model to feature number 1069, ASV1467
#> 2026-06-17 13:45:02.05 INFO::Fitting model to feature number 1070, ASV1468
#> 2026-06-17 13:45:02.05 INFO::Fitting model to feature number 1071, ASV1469
#> 2026-06-17 13:45:02.05 WARNING::Fitting problem for feature 1071 returning NA
#> 2026-06-17 13:45:02.06 INFO::Fitting model to feature number 1072, ASV1472
#> 2026-06-17 13:45:02.06 INFO::Fitting model to feature number 1073, ASV1477
#> 2026-06-17 13:45:02.06 INFO::Fitting model to feature number 1074, ASV1478
#> 2026-06-17 13:45:02.06 INFO::Fitting model to feature number 1075, ASV1479
#> 2026-06-17 13:45:02.07 INFO::Fitting model to feature number 1076, ASV1483
#> 2026-06-17 13:45:02.07 INFO::Fitting model to feature number 1077, ASV1484
#> 2026-06-17 13:45:02.07 INFO::Fitting model to feature number 1078, ASV1486
#> 2026-06-17 13:45:02.08 INFO::Fitting model to feature number 1079, ASV1487
#> 2026-06-17 13:45:02.08 INFO::Fitting model to feature number 1080, ASV1488
#> 2026-06-17 13:45:02.08 INFO::Fitting model to feature number 1081, ASV1490
#> 2026-06-17 13:45:02.09 INFO::Fitting model to feature number 1082, ASV1492
#> 2026-06-17 13:45:02.09 INFO::Fitting model to feature number 1083, ASV1493
#> 2026-06-17 13:45:02.09 INFO::Fitting model to feature number 1084, ASV1494
#> 2026-06-17 13:45:02.10 INFO::Fitting model to feature number 1085, ASV1495
#> 2026-06-17 13:45:02.10 INFO::Fitting model to feature number 1086, ASV1496
#> 2026-06-17 13:45:02.11 INFO::Fitting model to feature number 1087, ASV1497
#> 2026-06-17 13:45:02.11 INFO::Fitting model to feature number 1088, ASV1498
#> 2026-06-17 13:45:02.11 INFO::Fitting model to feature number 1089, ASV1500
#> 2026-06-17 13:45:02.12 INFO::Fitting model to feature number 1090, ASV1501
#> 2026-06-17 13:45:02.12 INFO::Fitting model to feature number 1091, ASV1502
#> 2026-06-17 13:45:02.12 INFO::Fitting model to feature number 1092, ASV1503
#> 2026-06-17 13:45:02.13 INFO::Fitting model to feature number 1093, ASV1504
#> 2026-06-17 13:45:02.13 WARNING::Fitting problem for feature 1093 returning NA
#> 2026-06-17 13:45:02.13 INFO::Fitting model to feature number 1094, ASV1507
#> 2026-06-17 13:45:02.13 INFO::Fitting model to feature number 1095, ASV1509
#> 2026-06-17 13:45:02.14 INFO::Fitting model to feature number 1096, ASV1510
#> 2026-06-17 13:45:02.14 INFO::Fitting model to feature number 1097, ASV1513
#> 2026-06-17 13:45:02.14 INFO::Fitting model to feature number 1098, ASV1515
#> 2026-06-17 13:45:02.14 INFO::Fitting model to feature number 1099, ASV1516
#> 2026-06-17 13:45:02.15 INFO::Fitting model to feature number 1100, ASV1517
#> 2026-06-17 13:45:02.15 INFO::Fitting model to feature number 1101, ASV1518
#> 2026-06-17 13:45:02.15 INFO::Fitting model to feature number 1102, ASV1521
#> 2026-06-17 13:45:02.16 INFO::Fitting model to feature number 1103, ASV1523
#> 2026-06-17 13:45:02.16 INFO::Fitting model to feature number 1104, ASV1525
#> 2026-06-17 13:45:02.16 INFO::Fitting model to feature number 1105, ASV1526
#> 2026-06-17 13:45:02.16 INFO::Fitting model to feature number 1106, ASV1527
#> 2026-06-17 13:45:02.17 INFO::Fitting model to feature number 1107, ASV1528
#> 2026-06-17 13:45:02.17 INFO::Fitting model to feature number 1108, ASV1531
#> 2026-06-17 13:45:02.17 INFO::Fitting model to feature number 1109, ASV1532
#> 2026-06-17 13:45:02.17 INFO::Fitting model to feature number 1110, ASV1533
#> 2026-06-17 13:45:02.18 INFO::Fitting model to feature number 1111, ASV1534
#> 2026-06-17 13:45:02.18 WARNING::Fitting problem for feature 1111 returning NA
#> 2026-06-17 13:45:02.18 INFO::Fitting model to feature number 1112, ASV1536
#> 2026-06-17 13:45:02.19 INFO::Fitting model to feature number 1113, ASV1537
#> 2026-06-17 13:45:02.19 INFO::Fitting model to feature number 1114, ASV1542
#> 2026-06-17 13:45:02.19 INFO::Fitting model to feature number 1115, ASV1543
#> 2026-06-17 13:45:02.19 INFO::Fitting model to feature number 1116, ASV1544
#> 2026-06-17 13:45:02.19 INFO::Fitting model to feature number 1117, ASV1545
#> 2026-06-17 13:45:02.20 INFO::Fitting model to feature number 1118, ASV1548
#> 2026-06-17 13:45:02.20 INFO::Fitting model to feature number 1119, ASV1549
#> 2026-06-17 13:45:02.20 INFO::Fitting model to feature number 1120, ASV1550
#> 2026-06-17 13:45:02.20 INFO::Fitting model to feature number 1121, ASV1551
#> 2026-06-17 13:45:02.20 INFO::Fitting model to feature number 1122, ASV1552
#> 2026-06-17 13:45:02.21 INFO::Fitting model to feature number 1123, ASV1553
#> 2026-06-17 13:45:02.21 INFO::Fitting model to feature number 1124, ASV1554
#> 2026-06-17 13:45:02.21 INFO::Fitting model to feature number 1125, ASV1557
#> 2026-06-17 13:45:02.22 INFO::Fitting model to feature number 1126, ASV1558
#> 2026-06-17 13:45:02.22 INFO::Fitting model to feature number 1127, ASV1560
#> 2026-06-17 13:45:02.22 INFO::Fitting model to feature number 1128, ASV1562
#> 2026-06-17 13:45:02.23 INFO::Fitting model to feature number 1129, ASV1564
#> 2026-06-17 13:45:02.23 INFO::Fitting model to feature number 1130, ASV1565
#> 2026-06-17 13:45:02.23 INFO::Fitting model to feature number 1131, ASV1568
#> 2026-06-17 13:45:02.23 INFO::Fitting model to feature number 1132, ASV1569
#> 2026-06-17 13:45:02.23 INFO::Fitting model to feature number 1133, ASV1570
#> 2026-06-17 13:45:02.24 INFO::Fitting model to feature number 1134, ASV1571
#> 2026-06-17 13:45:02.24 INFO::Fitting model to feature number 1135, ASV1573
#> 2026-06-17 13:45:02.24 INFO::Fitting model to feature number 1136, ASV1574
#> 2026-06-17 13:45:02.24 INFO::Fitting model to feature number 1137, ASV1575
#> 2026-06-17 13:45:02.25 INFO::Fitting model to feature number 1138, ASV1576
#> 2026-06-17 13:45:02.25 INFO::Fitting model to feature number 1139, ASV1577
#> 2026-06-17 13:45:02.25 INFO::Fitting model to feature number 1140, ASV1578
#> 2026-06-17 13:45:02.25 INFO::Fitting model to feature number 1141, ASV1580
#> 2026-06-17 13:45:02.26 INFO::Fitting model to feature number 1142, ASV1581
#> 2026-06-17 13:45:02.26 INFO::Fitting model to feature number 1143, ASV1584
#> 2026-06-17 13:45:02.26 INFO::Fitting model to feature number 1144, ASV1585
#> 2026-06-17 13:45:02.26 INFO::Fitting model to feature number 1145, ASV1586
#> 2026-06-17 13:45:02.27 INFO::Fitting model to feature number 1146, ASV1588
#> 2026-06-17 13:45:02.27 INFO::Fitting model to feature number 1147, ASV1589
#> 2026-06-17 13:45:02.27 INFO::Fitting model to feature number 1148, ASV1592
#> 2026-06-17 13:45:02.27 INFO::Fitting model to feature number 1149, ASV1593
#> 2026-06-17 13:45:02.28 INFO::Fitting model to feature number 1150, ASV1594
#> 2026-06-17 13:45:02.28 INFO::Fitting model to feature number 1151, ASV1601
#> 2026-06-17 13:45:02.28 INFO::Fitting model to feature number 1152, ASV1602
#> 2026-06-17 13:45:02.29 INFO::Fitting model to feature number 1153, ASV1603
#> 2026-06-17 13:45:02.29 INFO::Fitting model to feature number 1154, ASV1604
#> 2026-06-17 13:45:02.29 WARNING::Fitting problem for feature 1154 returning NA
#> 2026-06-17 13:45:02.29 INFO::Fitting model to feature number 1155, ASV1605
#> 2026-06-17 13:45:02.30 INFO::Fitting model to feature number 1156, ASV1606
#> 2026-06-17 13:45:02.30 INFO::Fitting model to feature number 1157, ASV1607
#> 2026-06-17 13:45:02.30 INFO::Fitting model to feature number 1158, ASV1609
#> 2026-06-17 13:45:02.31 WARNING::Fitting problem for feature 1158 returning NA
#> 2026-06-17 13:45:02.31 INFO::Fitting model to feature number 1159, ASV1612
#> 2026-06-17 13:45:02.31 INFO::Fitting model to feature number 1160, ASV1613
#> 2026-06-17 13:45:02.32 INFO::Fitting model to feature number 1161, ASV1614
#> 2026-06-17 13:45:02.32 INFO::Fitting model to feature number 1162, ASV1615
#> 2026-06-17 13:45:02.32 INFO::Fitting model to feature number 1163, ASV1616
#> 2026-06-17 13:45:02.32 INFO::Fitting model to feature number 1164, ASV1620
#> 2026-06-17 13:45:02.32 INFO::Fitting model to feature number 1165, ASV1621
#> 2026-06-17 13:45:02.33 INFO::Fitting model to feature number 1166, ASV1622
#> 2026-06-17 13:45:02.37 WARNING::Fitting problem for feature 1166 returning NA
#> 2026-06-17 13:45:02.37 INFO::Fitting model to feature number 1167, ASV1624
#> 2026-06-17 13:45:02.37 INFO::Fitting model to feature number 1168, ASV1625
#> 2026-06-17 13:45:02.37 INFO::Fitting model to feature number 1169, ASV1628
#> 2026-06-17 13:45:02.38 INFO::Fitting model to feature number 1170, ASV1629
#> 2026-06-17 13:45:02.38 INFO::Fitting model to feature number 1171, ASV1630
#> 2026-06-17 13:45:02.38 INFO::Fitting model to feature number 1172, ASV1631
#> 2026-06-17 13:45:02.39 INFO::Fitting model to feature number 1173, ASV1632
#> 2026-06-17 13:45:02.39 INFO::Fitting model to feature number 1174, ASV1633
#> 2026-06-17 13:45:02.39 INFO::Fitting model to feature number 1175, ASV1634
#> 2026-06-17 13:45:02.39 INFO::Fitting model to feature number 1176, ASV1635
#> 2026-06-17 13:45:02.40 INFO::Fitting model to feature number 1177, ASV1636
#> 2026-06-17 13:45:02.40 INFO::Fitting model to feature number 1178, ASV1639
#> 2026-06-17 13:45:02.40 INFO::Fitting model to feature number 1179, ASV1641
#> 2026-06-17 13:45:02.40 WARNING::Fitting problem for feature 1179 returning NA
#> 2026-06-17 13:45:02.41 INFO::Fitting model to feature number 1180, ASV1644
#> 2026-06-17 13:45:02.41 INFO::Fitting model to feature number 1181, ASV1645
#> 2026-06-17 13:45:02.41 INFO::Fitting model to feature number 1182, ASV1649
#> 2026-06-17 13:45:02.41 INFO::Fitting model to feature number 1183, ASV1650
#> 2026-06-17 13:45:02.42 INFO::Fitting model to feature number 1184, ASV1651
#> 2026-06-17 13:45:02.42 WARNING::Fitting problem for feature 1184 returning NA
#> 2026-06-17 13:45:02.42 INFO::Fitting model to feature number 1185, ASV1653
#> 2026-06-17 13:45:02.42 INFO::Fitting model to feature number 1186, ASV1654
#> 2026-06-17 13:45:02.43 INFO::Fitting model to feature number 1187, ASV1656
#> 2026-06-17 13:45:02.43 INFO::Fitting model to feature number 1188, ASV1657
#> 2026-06-17 13:45:02.43 INFO::Fitting model to feature number 1189, ASV1658
#> 2026-06-17 13:45:02.44 INFO::Fitting model to feature number 1190, ASV1661
#> 2026-06-17 13:45:02.44 INFO::Fitting model to feature number 1191, ASV1663
#> 2026-06-17 13:45:02.44 WARNING::Fitting problem for feature 1191 returning NA
#> 2026-06-17 13:45:02.44 INFO::Fitting model to feature number 1192, ASV1664
#> 2026-06-17 13:45:02.45 INFO::Fitting model to feature number 1193, ASV1666
#> 2026-06-17 13:45:02.45 INFO::Fitting model to feature number 1194, ASV1667
#> 2026-06-17 13:45:02.45 INFO::Fitting model to feature number 1195, ASV1670
#> 2026-06-17 13:45:02.45 INFO::Fitting model to feature number 1196, ASV1671
#> 2026-06-17 13:45:02.45 INFO::Fitting model to feature number 1197, ASV1673
#> 2026-06-17 13:45:02.46 INFO::Fitting model to feature number 1198, ASV1674
#> 2026-06-17 13:45:02.46 INFO::Fitting model to feature number 1199, ASV1677
#> 2026-06-17 13:45:02.46 INFO::Fitting model to feature number 1200, ASV1681
#> 2026-06-17 13:45:02.47 INFO::Fitting model to feature number 1201, ASV1683
#> 2026-06-17 13:45:02.47 INFO::Fitting model to feature number 1202, ASV1684
#> 2026-06-17 13:45:02.47 INFO::Fitting model to feature number 1203, ASV1687
#> 2026-06-17 13:45:02.47 WARNING::Fitting problem for feature 1203 returning NA
#> 2026-06-17 13:45:02.48 INFO::Fitting model to feature number 1204, ASV1689
#> 2026-06-17 13:45:02.48 INFO::Fitting model to feature number 1205, ASV1690
#> 2026-06-17 13:45:02.48 INFO::Fitting model to feature number 1206, ASV1691
#> 2026-06-17 13:45:02.49 INFO::Fitting model to feature number 1207, ASV1694
#> 2026-06-17 13:45:02.49 INFO::Fitting model to feature number 1208, ASV1697
#> 2026-06-17 13:45:02.49 INFO::Fitting model to feature number 1209, ASV1704
#> 2026-06-17 13:45:02.50 INFO::Fitting model to feature number 1210, ASV1706
#> 2026-06-17 13:45:02.50 INFO::Fitting model to feature number 1211, ASV1707
#> 2026-06-17 13:45:02.50 WARNING::Fitting problem for feature 1211 returning NA
#> 2026-06-17 13:45:02.50 INFO::Fitting model to feature number 1212, ASV1712
#> 2026-06-17 13:45:02.51 WARNING::Fitting problem for feature 1212 returning NA
#> 2026-06-17 13:45:02.52 INFO::Performing tests against medians
#> 2026-06-17 13:45:07.20 INFO::Counting total values for each feature
#> 2026-06-17 13:45:07.31 INFO::Running the logistic model component
#> 2026-06-17 13:45:07.34 INFO::Fitting model to feature number 1, ASV2
#> 2026-06-17 13:45:07.35 INFO::Fitting model to feature number 2, ASV6
#> 2026-06-17 13:45:07.35 INFO::Fitting model to feature number 3, ASV7
#> 2026-06-17 13:45:07.36 INFO::Fitting model to feature number 4, ASV8
#> 2026-06-17 13:45:07.40 INFO::Fitting model to feature number 5, ASV10
#> 2026-06-17 13:45:07.40 INFO::Fitting model to feature number 6, ASV12
#> 2026-06-17 13:45:07.41 INFO::Fitting model to feature number 7, ASV13
#> 2026-06-17 13:45:07.41 INFO::Fitting model to feature number 8, ASV18
#> 2026-06-17 13:45:07.42 INFO::Fitting model to feature number 9, ASV19
#> 2026-06-17 13:45:07.43 INFO::Fitting model to feature number 10, ASV22
#> 2026-06-17 13:45:07.43 INFO::Fitting model to feature number 11, ASV23
#> 2026-06-17 13:45:07.44 INFO::Fitting model to feature number 12, ASV24
#> 2026-06-17 13:45:07.44 INFO::Fitting model to feature number 13, ASV25
#> 2026-06-17 13:45:07.45 INFO::Fitting model to feature number 14, ASV26
#> 2026-06-17 13:45:07.45 INFO::Fitting model to feature number 15, ASV27
#> 2026-06-17 13:45:07.46 INFO::Fitting model to feature number 16, ASV28
#> 2026-06-17 13:45:07.47 INFO::Fitting model to feature number 17, ASV29
#> 2026-06-17 13:45:07.47 INFO::Fitting model to feature number 18, ASV31
#> 2026-06-17 13:45:07.48 INFO::Fitting model to feature number 19, ASV32
#> 2026-06-17 13:45:07.48 INFO::Fitting model to feature number 20, ASV33
#> 2026-06-17 13:45:07.49 INFO::Fitting model to feature number 21, ASV34
#> 2026-06-17 13:45:07.50 INFO::Fitting model to feature number 22, ASV35
#> 2026-06-17 13:45:07.50 INFO::Fitting model to feature number 23, ASV38
#> 2026-06-17 13:45:07.51 INFO::Fitting model to feature number 24, ASV41
#> 2026-06-17 13:45:07.51 INFO::Fitting model to feature number 25, ASV42
#> 2026-06-17 13:45:07.52 INFO::Fitting model to feature number 26, ASV43
#> 2026-06-17 13:45:07.52 INFO::Fitting model to feature number 27, ASV45
#> 2026-06-17 13:45:07.53 INFO::Fitting model to feature number 28, ASV46
#> 2026-06-17 13:45:07.54 INFO::Fitting model to feature number 29, ASV47
#> 2026-06-17 13:45:07.54 INFO::Fitting model to feature number 30, ASV48
#> 2026-06-17 13:45:07.55 INFO::Fitting model to feature number 31, ASV49
#> 2026-06-17 13:45:07.55 INFO::Fitting model to feature number 32, ASV50
#> 2026-06-17 13:45:07.56 INFO::Fitting model to feature number 33, ASV51
#> 2026-06-17 13:45:07.57 INFO::Fitting model to feature number 34, ASV52
#> 2026-06-17 13:45:07.57 INFO::Fitting model to feature number 35, ASV53
#> 2026-06-17 13:45:07.58 INFO::Fitting model to feature number 36, ASV55
#> 2026-06-17 13:45:07.58 INFO::Fitting model to feature number 37, ASV56
#> 2026-06-17 13:45:07.59 INFO::Fitting model to feature number 38, ASV57
#> 2026-06-17 13:45:07.59 INFO::Fitting model to feature number 39, ASV58
#> 2026-06-17 13:45:07.60 INFO::Fitting model to feature number 40, ASV59
#> 2026-06-17 13:45:07.61 INFO::Fitting model to feature number 41, ASV60
#> 2026-06-17 13:45:07.61 INFO::Fitting model to feature number 42, ASV61
#> 2026-06-17 13:45:07.62 INFO::Fitting model to feature number 43, ASV62
#> 2026-06-17 13:45:07.62 INFO::Fitting model to feature number 44, ASV63
#> 2026-06-17 13:45:07.63 INFO::Fitting model to feature number 45, ASV64
#> 2026-06-17 13:45:07.63 INFO::Fitting model to feature number 46, ASV65
#> 2026-06-17 13:45:07.64 INFO::Fitting model to feature number 47, ASV67
#> 2026-06-17 13:45:07.65 INFO::Fitting model to feature number 48, ASV68
#> 2026-06-17 13:45:07.65 INFO::Fitting model to feature number 49, ASV69
#> 2026-06-17 13:45:07.66 INFO::Fitting model to feature number 50, ASV70
#> 2026-06-17 13:45:07.66 INFO::Fitting model to feature number 51, ASV71
#> 2026-06-17 13:45:07.67 INFO::Fitting model to feature number 52, ASV72
#> 2026-06-17 13:45:07.67 INFO::Fitting model to feature number 53, ASV75
#> 2026-06-17 13:45:07.68 INFO::Fitting model to feature number 54, ASV77
#> 2026-06-17 13:45:07.69 INFO::Fitting model to feature number 55, ASV78
#> 2026-06-17 13:45:07.69 INFO::Fitting model to feature number 56, ASV79
#> 2026-06-17 13:45:07.70 INFO::Fitting model to feature number 57, ASV80
#> 2026-06-17 13:45:07.70 INFO::Fitting model to feature number 58, ASV81
#> 2026-06-17 13:45:07.71 INFO::Fitting model to feature number 59, ASV82
#> 2026-06-17 13:45:07.71 INFO::Fitting model to feature number 60, ASV83
#> 2026-06-17 13:45:07.72 INFO::Fitting model to feature number 61, ASV84
#> 2026-06-17 13:45:07.73 INFO::Fitting model to feature number 62, ASV85
#> 2026-06-17 13:45:07.73 INFO::Fitting model to feature number 63, ASV87
#> 2026-06-17 13:45:07.74 INFO::Fitting model to feature number 64, ASV89
#> 2026-06-17 13:45:07.74 INFO::Fitting model to feature number 65, ASV90
#> 2026-06-17 13:45:07.75 INFO::Fitting model to feature number 66, ASV91
#> 2026-06-17 13:45:07.76 INFO::Fitting model to feature number 67, ASV92
#> 2026-06-17 13:45:07.76 INFO::Fitting model to feature number 68, ASV93
#> 2026-06-17 13:45:07.77 INFO::Fitting model to feature number 69, ASV94
#> 2026-06-17 13:45:07.77 INFO::Fitting model to feature number 70, ASV95
#> 2026-06-17 13:45:07.78 INFO::Fitting model to feature number 71, ASV96
#> 2026-06-17 13:45:07.79 INFO::Fitting model to feature number 72, ASV98
#> 2026-06-17 13:45:07.79 INFO::Fitting model to feature number 73, ASV99
#> 2026-06-17 13:45:07.80 INFO::Fitting model to feature number 74, ASV100
#> 2026-06-17 13:45:07.80 INFO::Fitting model to feature number 75, ASV101
#> 2026-06-17 13:45:07.81 INFO::Fitting model to feature number 76, ASV103
#> 2026-06-17 13:45:07.81 INFO::Fitting model to feature number 77, ASV104
#> 2026-06-17 13:45:07.82 INFO::Fitting model to feature number 78, ASV105
#> 2026-06-17 13:45:07.83 INFO::Fitting model to feature number 79, ASV106
#> 2026-06-17 13:45:07.83 INFO::Fitting model to feature number 80, ASV107
#> 2026-06-17 13:45:07.84 INFO::Fitting model to feature number 81, ASV109
#> 2026-06-17 13:45:07.84 INFO::Fitting model to feature number 82, ASV111
#> 2026-06-17 13:45:07.85 INFO::Fitting model to feature number 83, ASV112
#> 2026-06-17 13:45:07.85 INFO::Fitting model to feature number 84, ASV113
#> 2026-06-17 13:45:07.86 INFO::Fitting model to feature number 85, ASV115
#> 2026-06-17 13:45:07.87 INFO::Fitting model to feature number 86, ASV116
#> 2026-06-17 13:45:07.87 INFO::Fitting model to feature number 87, ASV118
#> 2026-06-17 13:45:07.88 INFO::Fitting model to feature number 88, ASV119
#> 2026-06-17 13:45:07.88 INFO::Fitting model to feature number 89, ASV120
#> 2026-06-17 13:45:07.89 INFO::Fitting model to feature number 90, ASV121
#> 2026-06-17 13:45:07.89 INFO::Fitting model to feature number 91, ASV124
#> 2026-06-17 13:45:07.90 INFO::Fitting model to feature number 92, ASV126
#> 2026-06-17 13:45:07.91 INFO::Fitting model to feature number 93, ASV127
#> 2026-06-17 13:45:07.91 INFO::Fitting model to feature number 94, ASV128
#> 2026-06-17 13:45:07.92 INFO::Fitting model to feature number 95, ASV129
#> 2026-06-17 13:45:07.92 INFO::Fitting model to feature number 96, ASV130
#> 2026-06-17 13:45:07.93 INFO::Fitting model to feature number 97, ASV131
#> 2026-06-17 13:45:07.93 INFO::Fitting model to feature number 98, ASV132
#> 2026-06-17 13:45:07.94 INFO::Fitting model to feature number 99, ASV133
#> 2026-06-17 13:45:07.94 INFO::Fitting model to feature number 100, ASV135
#> 2026-06-17 13:45:07.95 INFO::Fitting model to feature number 101, ASV136
#> 2026-06-17 13:45:07.96 INFO::Fitting model to feature number 102, ASV137
#> 2026-06-17 13:45:07.96 INFO::Fitting model to feature number 103, ASV138
#> 2026-06-17 13:45:07.97 INFO::Fitting model to feature number 104, ASV139
#> 2026-06-17 13:45:07.97 INFO::Fitting model to feature number 105, ASV142
#> 2026-06-17 13:45:07.98 INFO::Fitting model to feature number 106, ASV143
#> 2026-06-17 13:45:07.99 INFO::Fitting model to feature number 107, ASV144
#> 2026-06-17 13:45:07.99 INFO::Fitting model to feature number 108, ASV145
#> 2026-06-17 13:45:08.00 INFO::Fitting model to feature number 109, ASV146
#> 2026-06-17 13:45:08.00 INFO::Fitting model to feature number 110, ASV147
#> 2026-06-17 13:45:08.01 INFO::Fitting model to feature number 111, ASV148
#> 2026-06-17 13:45:08.01 INFO::Fitting model to feature number 112, ASV149
#> 2026-06-17 13:45:08.02 INFO::Fitting model to feature number 113, ASV150
#> 2026-06-17 13:45:08.03 INFO::Fitting model to feature number 114, ASV151
#> 2026-06-17 13:45:08.03 INFO::Fitting model to feature number 115, ASV153
#> 2026-06-17 13:45:08.04 INFO::Fitting model to feature number 116, ASV154
#> 2026-06-17 13:45:08.04 INFO::Fitting model to feature number 117, ASV157
#> 2026-06-17 13:45:08.05 INFO::Fitting model to feature number 118, ASV158
#> 2026-06-17 13:45:08.05 INFO::Fitting model to feature number 119, ASV159
#> 2026-06-17 13:45:08.06 INFO::Fitting model to feature number 120, ASV162
#> 2026-06-17 13:45:08.07 INFO::Fitting model to feature number 121, ASV163
#> 2026-06-17 13:45:08.07 INFO::Fitting model to feature number 122, ASV164
#> 2026-06-17 13:45:08.08 INFO::Fitting model to feature number 123, ASV166
#> 2026-06-17 13:45:08.08 INFO::Fitting model to feature number 124, ASV167
#> 2026-06-17 13:45:08.09 INFO::Fitting model to feature number 125, ASV170
#> 2026-06-17 13:45:08.09 INFO::Fitting model to feature number 126, ASV171
#> 2026-06-17 13:45:08.10 INFO::Fitting model to feature number 127, ASV172
#> 2026-06-17 13:45:08.10 INFO::Fitting model to feature number 128, ASV173
#> 2026-06-17 13:45:08.11 INFO::Fitting model to feature number 129, ASV175
#> 2026-06-17 13:45:08.12 INFO::Fitting model to feature number 130, ASV176
#> 2026-06-17 13:45:08.12 INFO::Fitting model to feature number 131, ASV178
#> 2026-06-17 13:45:08.13 INFO::Fitting model to feature number 132, ASV179
#> 2026-06-17 13:45:08.13 INFO::Fitting model to feature number 133, ASV181
#> 2026-06-17 13:45:08.14 INFO::Fitting model to feature number 134, ASV182
#> 2026-06-17 13:45:08.14 INFO::Fitting model to feature number 135, ASV183
#> 2026-06-17 13:45:08.15 INFO::Fitting model to feature number 136, ASV186
#> 2026-06-17 13:45:08.16 INFO::Fitting model to feature number 137, ASV187
#> 2026-06-17 13:45:08.16 INFO::Fitting model to feature number 138, ASV189
#> 2026-06-17 13:45:08.17 INFO::Fitting model to feature number 139, ASV190
#> 2026-06-17 13:45:08.17 INFO::Fitting model to feature number 140, ASV192
#> 2026-06-17 13:45:08.18 INFO::Fitting model to feature number 141, ASV193
#> 2026-06-17 13:45:08.19 INFO::Fitting model to feature number 142, ASV194
#> 2026-06-17 13:45:08.19 INFO::Fitting model to feature number 143, ASV195
#> 2026-06-17 13:45:08.20 INFO::Fitting model to feature number 144, ASV196
#> 2026-06-17 13:45:08.20 INFO::Fitting model to feature number 145, ASV197
#> 2026-06-17 13:45:08.21 INFO::Fitting model to feature number 146, ASV198
#> 2026-06-17 13:45:08.21 INFO::Fitting model to feature number 147, ASV199
#> 2026-06-17 13:45:08.22 INFO::Fitting model to feature number 148, ASV201
#> 2026-06-17 13:45:08.23 INFO::Fitting model to feature number 149, ASV202
#> 2026-06-17 13:45:08.23 INFO::Fitting model to feature number 150, ASV203
#> 2026-06-17 13:45:08.24 INFO::Fitting model to feature number 151, ASV205
#> 2026-06-17 13:45:08.24 INFO::Fitting model to feature number 152, ASV208
#> 2026-06-17 13:45:08.25 INFO::Fitting model to feature number 153, ASV209
#> 2026-06-17 13:45:08.25 INFO::Fitting model to feature number 154, ASV210
#> 2026-06-17 13:45:08.26 INFO::Fitting model to feature number 155, ASV211
#> 2026-06-17 13:45:08.27 INFO::Fitting model to feature number 156, ASV212
#> 2026-06-17 13:45:08.27 INFO::Fitting model to feature number 157, ASV214
#> 2026-06-17 13:45:08.28 INFO::Fitting model to feature number 158, ASV215
#> 2026-06-17 13:45:08.28 INFO::Fitting model to feature number 159, ASV216
#> 2026-06-17 13:45:08.29 INFO::Fitting model to feature number 160, ASV217
#> 2026-06-17 13:45:08.30 INFO::Fitting model to feature number 161, ASV219
#> 2026-06-17 13:45:08.30 INFO::Fitting model to feature number 162, ASV221
#> 2026-06-17 13:45:08.31 INFO::Fitting model to feature number 163, ASV222
#> 2026-06-17 13:45:08.31 INFO::Fitting model to feature number 164, ASV223
#> 2026-06-17 13:45:08.32 INFO::Fitting model to feature number 165, ASV224
#> 2026-06-17 13:45:08.32 INFO::Fitting model to feature number 166, ASV226
#> 2026-06-17 13:45:08.33 INFO::Fitting model to feature number 167, ASV227
#> 2026-06-17 13:45:08.34 INFO::Fitting model to feature number 168, ASV228
#> 2026-06-17 13:45:08.34 INFO::Fitting model to feature number 169, ASV229
#> 2026-06-17 13:45:08.35 INFO::Fitting model to feature number 170, ASV231
#> 2026-06-17 13:45:08.35 INFO::Fitting model to feature number 171, ASV233
#> 2026-06-17 13:45:08.36 INFO::Fitting model to feature number 172, ASV234
#> 2026-06-17 13:45:08.37 INFO::Fitting model to feature number 173, ASV235
#> 2026-06-17 13:45:08.37 INFO::Fitting model to feature number 174, ASV237
#> 2026-06-17 13:45:08.38 INFO::Fitting model to feature number 175, ASV238
#> 2026-06-17 13:45:08.42 INFO::Fitting model to feature number 176, ASV239
#> 2026-06-17 13:45:08.43 INFO::Fitting model to feature number 177, ASV240
#> 2026-06-17 13:45:08.44 INFO::Fitting model to feature number 178, ASV243
#> 2026-06-17 13:45:08.44 INFO::Fitting model to feature number 179, ASV244
#> 2026-06-17 13:45:08.45 INFO::Fitting model to feature number 180, ASV245
#> 2026-06-17 13:45:08.45 INFO::Fitting model to feature number 181, ASV246
#> 2026-06-17 13:45:08.46 INFO::Fitting model to feature number 182, ASV247
#> 2026-06-17 13:45:08.46 INFO::Fitting model to feature number 183, ASV248
#> 2026-06-17 13:45:08.47 INFO::Fitting model to feature number 184, ASV249
#> 2026-06-17 13:45:08.48 INFO::Fitting model to feature number 185, ASV251
#> 2026-06-17 13:45:08.48 INFO::Fitting model to feature number 186, ASV254
#> 2026-06-17 13:45:08.49 INFO::Fitting model to feature number 187, ASV255
#> 2026-06-17 13:45:08.49 INFO::Fitting model to feature number 188, ASV256
#> 2026-06-17 13:45:08.50 INFO::Fitting model to feature number 189, ASV257
#> 2026-06-17 13:45:08.51 INFO::Fitting model to feature number 190, ASV258
#> 2026-06-17 13:45:08.51 INFO::Fitting model to feature number 191, ASV260
#> 2026-06-17 13:45:08.52 INFO::Fitting model to feature number 192, ASV261
#> 2026-06-17 13:45:08.52 INFO::Fitting model to feature number 193, ASV262
#> 2026-06-17 13:45:08.53 INFO::Fitting model to feature number 194, ASV263
#> 2026-06-17 13:45:08.53 INFO::Fitting model to feature number 195, ASV264
#> 2026-06-17 13:45:08.54 INFO::Fitting model to feature number 196, ASV265
#> 2026-06-17 13:45:08.55 INFO::Fitting model to feature number 197, ASV266
#> 2026-06-17 13:45:08.55 INFO::Fitting model to feature number 198, ASV267
#> 2026-06-17 13:45:08.56 INFO::Fitting model to feature number 199, ASV270
#> 2026-06-17 13:45:08.56 INFO::Fitting model to feature number 200, ASV271
#> 2026-06-17 13:45:08.57 INFO::Fitting model to feature number 201, ASV272
#> 2026-06-17 13:45:08.58 INFO::Fitting model to feature number 202, ASV273
#> 2026-06-17 13:45:08.58 INFO::Fitting model to feature number 203, ASV274
#> 2026-06-17 13:45:08.59 INFO::Fitting model to feature number 204, ASV277
#> 2026-06-17 13:45:08.59 INFO::Fitting model to feature number 205, ASV278
#> 2026-06-17 13:45:08.60 INFO::Fitting model to feature number 206, ASV279
#> 2026-06-17 13:45:08.60 INFO::Fitting model to feature number 207, ASV282
#> 2026-06-17 13:45:08.61 INFO::Fitting model to feature number 208, ASV283
#> 2026-06-17 13:45:08.62 INFO::Fitting model to feature number 209, ASV285
#> 2026-06-17 13:45:08.62 INFO::Fitting model to feature number 210, ASV286
#> 2026-06-17 13:45:08.63 INFO::Fitting model to feature number 211, ASV287
#> 2026-06-17 13:45:08.63 INFO::Fitting model to feature number 212, ASV288
#> 2026-06-17 13:45:08.64 INFO::Fitting model to feature number 213, ASV292
#> 2026-06-17 13:45:08.64 INFO::Fitting model to feature number 214, ASV293
#> 2026-06-17 13:45:08.65 INFO::Fitting model to feature number 215, ASV294
#> 2026-06-17 13:45:08.66 INFO::Fitting model to feature number 216, ASV295
#> 2026-06-17 13:45:08.66 INFO::Fitting model to feature number 217, ASV297
#> 2026-06-17 13:45:08.67 INFO::Fitting model to feature number 218, ASV300
#> 2026-06-17 13:45:08.67 INFO::Fitting model to feature number 219, ASV302
#> 2026-06-17 13:45:08.68 INFO::Fitting model to feature number 220, ASV307
#> 2026-06-17 13:45:08.68 INFO::Fitting model to feature number 221, ASV309
#> 2026-06-17 13:45:08.69 INFO::Fitting model to feature number 222, ASV310
#> 2026-06-17 13:45:08.70 INFO::Fitting model to feature number 223, ASV311
#> 2026-06-17 13:45:08.70 INFO::Fitting model to feature number 224, ASV312
#> 2026-06-17 13:45:08.71 INFO::Fitting model to feature number 225, ASV313
#> 2026-06-17 13:45:08.71 INFO::Fitting model to feature number 226, ASV314
#> 2026-06-17 13:45:08.72 INFO::Fitting model to feature number 227, ASV315
#> 2026-06-17 13:45:08.72 INFO::Fitting model to feature number 228, ASV316
#> 2026-06-17 13:45:08.73 INFO::Fitting model to feature number 229, ASV318
#> 2026-06-17 13:45:08.74 INFO::Fitting model to feature number 230, ASV320
#> 2026-06-17 13:45:08.74 INFO::Fitting model to feature number 231, ASV322
#> 2026-06-17 13:45:08.75 INFO::Fitting model to feature number 232, ASV323
#> 2026-06-17 13:45:08.75 INFO::Fitting model to feature number 233, ASV327
#> 2026-06-17 13:45:08.76 INFO::Fitting model to feature number 234, ASV328
#> 2026-06-17 13:45:08.77 INFO::Fitting model to feature number 235, ASV329
#> 2026-06-17 13:45:08.77 INFO::Fitting model to feature number 236, ASV330
#> 2026-06-17 13:45:08.78 INFO::Fitting model to feature number 237, ASV333
#> 2026-06-17 13:45:08.78 INFO::Fitting model to feature number 238, ASV334
#> 2026-06-17 13:45:08.79 INFO::Fitting model to feature number 239, ASV337
#> 2026-06-17 13:45:08.80 INFO::Fitting model to feature number 240, ASV338
#> 2026-06-17 13:45:08.80 INFO::Fitting model to feature number 241, ASV339
#> 2026-06-17 13:45:08.81 INFO::Fitting model to feature number 242, ASV340
#> 2026-06-17 13:45:08.81 INFO::Fitting model to feature number 243, ASV341
#> 2026-06-17 13:45:08.82 INFO::Fitting model to feature number 244, ASV342
#> 2026-06-17 13:45:08.82 INFO::Fitting model to feature number 245, ASV344
#> 2026-06-17 13:45:08.83 INFO::Fitting model to feature number 246, ASV345
#> 2026-06-17 13:45:08.83 INFO::Fitting model to feature number 247, ASV346
#> 2026-06-17 13:45:08.84 INFO::Fitting model to feature number 248, ASV347
#> 2026-06-17 13:45:08.85 INFO::Fitting model to feature number 249, ASV348
#> 2026-06-17 13:45:08.85 INFO::Fitting model to feature number 250, ASV349
#> 2026-06-17 13:45:08.86 INFO::Fitting model to feature number 251, ASV350
#> 2026-06-17 13:45:08.86 INFO::Fitting model to feature number 252, ASV351
#> 2026-06-17 13:45:08.87 INFO::Fitting model to feature number 253, ASV353
#> 2026-06-17 13:45:08.87 INFO::Fitting model to feature number 254, ASV355
#> 2026-06-17 13:45:08.88 INFO::Fitting model to feature number 255, ASV356
#> 2026-06-17 13:45:08.89 INFO::Fitting model to feature number 256, ASV357
#> 2026-06-17 13:45:08.89 INFO::Fitting model to feature number 257, ASV358
#> 2026-06-17 13:45:08.90 INFO::Fitting model to feature number 258, ASV359
#> 2026-06-17 13:45:08.90 INFO::Fitting model to feature number 259, ASV360
#> 2026-06-17 13:45:08.91 INFO::Fitting model to feature number 260, ASV362
#> 2026-06-17 13:45:08.91 INFO::Fitting model to feature number 261, ASV364
#> 2026-06-17 13:45:08.92 INFO::Fitting model to feature number 262, ASV365
#> 2026-06-17 13:45:08.93 INFO::Fitting model to feature number 263, ASV366
#> 2026-06-17 13:45:08.93 INFO::Fitting model to feature number 264, ASV367
#> 2026-06-17 13:45:08.94 INFO::Fitting model to feature number 265, ASV368
#> 2026-06-17 13:45:08.94 INFO::Fitting model to feature number 266, ASV369
#> 2026-06-17 13:45:08.95 INFO::Fitting model to feature number 267, ASV371
#> 2026-06-17 13:45:08.95 INFO::Fitting model to feature number 268, ASV372
#> 2026-06-17 13:45:08.96 INFO::Fitting model to feature number 269, ASV373
#> 2026-06-17 13:45:08.97 INFO::Fitting model to feature number 270, ASV375
#> 2026-06-17 13:45:08.97 INFO::Fitting model to feature number 271, ASV376
#> 2026-06-17 13:45:08.98 INFO::Fitting model to feature number 272, ASV377
#> 2026-06-17 13:45:08.98 INFO::Fitting model to feature number 273, ASV378
#> 2026-06-17 13:45:08.99 INFO::Fitting model to feature number 274, ASV380
#> 2026-06-17 13:45:08.99 INFO::Fitting model to feature number 275, ASV382
#> 2026-06-17 13:45:09.00 INFO::Fitting model to feature number 276, ASV383
#> 2026-06-17 13:45:09.01 INFO::Fitting model to feature number 277, ASV385
#> 2026-06-17 13:45:09.01 INFO::Fitting model to feature number 278, ASV388
#> 2026-06-17 13:45:09.02 INFO::Fitting model to feature number 279, ASV390
#> 2026-06-17 13:45:09.02 INFO::Fitting model to feature number 280, ASV395
#> 2026-06-17 13:45:09.03 INFO::Fitting model to feature number 281, ASV399
#> 2026-06-17 13:45:09.03 INFO::Fitting model to feature number 282, ASV401
#> 2026-06-17 13:45:09.04 INFO::Fitting model to feature number 283, ASV402
#> 2026-06-17 13:45:09.05 INFO::Fitting model to feature number 284, ASV404
#> 2026-06-17 13:45:09.05 INFO::Fitting model to feature number 285, ASV405
#> 2026-06-17 13:45:09.06 INFO::Fitting model to feature number 286, ASV409
#> 2026-06-17 13:45:09.06 INFO::Fitting model to feature number 287, ASV411
#> 2026-06-17 13:45:09.07 INFO::Fitting model to feature number 288, ASV412
#> 2026-06-17 13:45:09.07 INFO::Fitting model to feature number 289, ASV413
#> 2026-06-17 13:45:09.08 INFO::Fitting model to feature number 290, ASV414
#> 2026-06-17 13:45:09.09 INFO::Fitting model to feature number 291, ASV415
#> 2026-06-17 13:45:09.09 INFO::Fitting model to feature number 292, ASV416
#> 2026-06-17 13:45:09.10 INFO::Fitting model to feature number 293, ASV417
#> 2026-06-17 13:45:09.10 INFO::Fitting model to feature number 294, ASV418
#> 2026-06-17 13:45:09.11 INFO::Fitting model to feature number 295, ASV419
#> 2026-06-17 13:45:09.11 INFO::Fitting model to feature number 296, ASV420
#> 2026-06-17 13:45:09.12 INFO::Fitting model to feature number 297, ASV422
#> 2026-06-17 13:45:09.13 INFO::Fitting model to feature number 298, ASV423
#> 2026-06-17 13:45:09.13 INFO::Fitting model to feature number 299, ASV424
#> 2026-06-17 13:45:09.14 INFO::Fitting model to feature number 300, ASV426
#> 2026-06-17 13:45:09.14 INFO::Fitting model to feature number 301, ASV427
#> 2026-06-17 13:45:09.15 INFO::Fitting model to feature number 302, ASV428
#> 2026-06-17 13:45:09.16 INFO::Fitting model to feature number 303, ASV429
#> 2026-06-17 13:45:09.16 INFO::Fitting model to feature number 304, ASV431
#> 2026-06-17 13:45:09.17 INFO::Fitting model to feature number 305, ASV433
#> 2026-06-17 13:45:09.17 INFO::Fitting model to feature number 306, ASV434
#> 2026-06-17 13:45:09.18 INFO::Fitting model to feature number 307, ASV438
#> 2026-06-17 13:45:09.18 INFO::Fitting model to feature number 308, ASV439
#> 2026-06-17 13:45:09.19 INFO::Fitting model to feature number 309, ASV440
#> 2026-06-17 13:45:09.20 INFO::Fitting model to feature number 310, ASV442
#> 2026-06-17 13:45:09.20 INFO::Fitting model to feature number 311, ASV443
#> 2026-06-17 13:45:09.21 INFO::Fitting model to feature number 312, ASV444
#> 2026-06-17 13:45:09.21 INFO::Fitting model to feature number 313, ASV445
#> 2026-06-17 13:45:09.22 INFO::Fitting model to feature number 314, ASV446
#> 2026-06-17 13:45:09.22 INFO::Fitting model to feature number 315, ASV448
#> 2026-06-17 13:45:09.23 INFO::Fitting model to feature number 316, ASV449
#> 2026-06-17 13:45:09.24 INFO::Fitting model to feature number 317, ASV451
#> 2026-06-17 13:45:09.24 INFO::Fitting model to feature number 318, ASV452
#> 2026-06-17 13:45:09.25 INFO::Fitting model to feature number 319, ASV453
#> 2026-06-17 13:45:09.25 INFO::Fitting model to feature number 320, ASV454
#> 2026-06-17 13:45:09.26 INFO::Fitting model to feature number 321, ASV456
#> 2026-06-17 13:45:09.26 INFO::Fitting model to feature number 322, ASV457
#> 2026-06-17 13:45:09.27 INFO::Fitting model to feature number 323, ASV458
#> 2026-06-17 13:45:09.28 INFO::Fitting model to feature number 324, ASV459
#> 2026-06-17 13:45:09.28 INFO::Fitting model to feature number 325, ASV460
#> 2026-06-17 13:45:09.29 INFO::Fitting model to feature number 326, ASV461
#> 2026-06-17 13:45:09.29 INFO::Fitting model to feature number 327, ASV462
#> 2026-06-17 13:45:09.30 INFO::Fitting model to feature number 328, ASV463
#> 2026-06-17 13:45:09.31 INFO::Fitting model to feature number 329, ASV464
#> 2026-06-17 13:45:09.35 INFO::Fitting model to feature number 330, ASV465
#> 2026-06-17 13:45:09.36 INFO::Fitting model to feature number 331, ASV466
#> 2026-06-17 13:45:09.36 INFO::Fitting model to feature number 332, ASV467
#> 2026-06-17 13:45:09.37 INFO::Fitting model to feature number 333, ASV468
#> 2026-06-17 13:45:09.37 INFO::Fitting model to feature number 334, ASV469
#> 2026-06-17 13:45:09.38 INFO::Fitting model to feature number 335, ASV470
#> 2026-06-17 13:45:09.39 INFO::Fitting model to feature number 336, ASV472
#> 2026-06-17 13:45:09.39 INFO::Fitting model to feature number 337, ASV473
#> 2026-06-17 13:45:09.40 INFO::Fitting model to feature number 338, ASV474
#> 2026-06-17 13:45:09.40 INFO::Fitting model to feature number 339, ASV475
#> 2026-06-17 13:45:09.41 INFO::Fitting model to feature number 340, ASV476
#> 2026-06-17 13:45:09.42 INFO::Fitting model to feature number 341, ASV477
#> 2026-06-17 13:45:09.42 INFO::Fitting model to feature number 342, ASV478
#> 2026-06-17 13:45:09.43 INFO::Fitting model to feature number 343, ASV479
#> 2026-06-17 13:45:09.44 INFO::Fitting model to feature number 344, ASV481
#> 2026-06-17 13:45:09.44 INFO::Fitting model to feature number 345, ASV483
#> 2026-06-17 13:45:09.45 INFO::Fitting model to feature number 346, ASV488
#> 2026-06-17 13:45:09.45 INFO::Fitting model to feature number 347, ASV489
#> 2026-06-17 13:45:09.46 INFO::Fitting model to feature number 348, ASV490
#> 2026-06-17 13:45:09.46 INFO::Fitting model to feature number 349, ASV491
#> 2026-06-17 13:45:09.47 INFO::Fitting model to feature number 350, ASV492
#> 2026-06-17 13:45:09.48 INFO::Fitting model to feature number 351, ASV493
#> 2026-06-17 13:45:09.48 INFO::Fitting model to feature number 352, ASV494
#> 2026-06-17 13:45:09.49 INFO::Fitting model to feature number 353, ASV496
#> 2026-06-17 13:45:09.49 INFO::Fitting model to feature number 354, ASV498
#> 2026-06-17 13:45:09.50 INFO::Fitting model to feature number 355, ASV499
#> 2026-06-17 13:45:09.51 INFO::Fitting model to feature number 356, ASV500
#> 2026-06-17 13:45:09.51 INFO::Fitting model to feature number 357, ASV501
#> 2026-06-17 13:45:09.52 INFO::Fitting model to feature number 358, ASV502
#> 2026-06-17 13:45:09.52 INFO::Fitting model to feature number 359, ASV504
#> 2026-06-17 13:45:09.53 INFO::Fitting model to feature number 360, ASV505
#> 2026-06-17 13:45:09.54 INFO::Fitting model to feature number 361, ASV507
#> 2026-06-17 13:45:09.54 INFO::Fitting model to feature number 362, ASV508
#> 2026-06-17 13:45:09.55 INFO::Fitting model to feature number 363, ASV509
#> 2026-06-17 13:45:09.55 INFO::Fitting model to feature number 364, ASV511
#> 2026-06-17 13:45:09.56 INFO::Fitting model to feature number 365, ASV512
#> 2026-06-17 13:45:09.57 INFO::Fitting model to feature number 366, ASV514
#> 2026-06-17 13:45:09.57 INFO::Fitting model to feature number 367, ASV515
#> 2026-06-17 13:45:09.58 INFO::Fitting model to feature number 368, ASV516
#> 2026-06-17 13:45:09.58 INFO::Fitting model to feature number 369, ASV517
#> 2026-06-17 13:45:09.59 INFO::Fitting model to feature number 370, ASV519
#> 2026-06-17 13:45:09.60 INFO::Fitting model to feature number 371, ASV520
#> 2026-06-17 13:45:09.60 INFO::Fitting model to feature number 372, ASV521
#> 2026-06-17 13:45:09.61 INFO::Fitting model to feature number 373, ASV522
#> 2026-06-17 13:45:09.61 INFO::Fitting model to feature number 374, ASV523
#> 2026-06-17 13:45:09.62 INFO::Fitting model to feature number 375, ASV526
#> 2026-06-17 13:45:09.63 INFO::Fitting model to feature number 376, ASV527
#> 2026-06-17 13:45:09.63 INFO::Fitting model to feature number 377, ASV530
#> 2026-06-17 13:45:09.64 INFO::Fitting model to feature number 378, ASV531
#> 2026-06-17 13:45:09.64 INFO::Fitting model to feature number 379, ASV533
#> 2026-06-17 13:45:09.65 INFO::Fitting model to feature number 380, ASV534
#> 2026-06-17 13:45:09.65 INFO::Fitting model to feature number 381, ASV535
#> 2026-06-17 13:45:09.66 INFO::Fitting model to feature number 382, ASV536
#> 2026-06-17 13:45:09.67 INFO::Fitting model to feature number 383, ASV538
#> 2026-06-17 13:45:09.67 INFO::Fitting model to feature number 384, ASV539
#> 2026-06-17 13:45:09.68 INFO::Fitting model to feature number 385, ASV540
#> 2026-06-17 13:45:09.68 INFO::Fitting model to feature number 386, ASV541
#> 2026-06-17 13:45:09.69 INFO::Fitting model to feature number 387, ASV542
#> 2026-06-17 13:45:09.70 INFO::Fitting model to feature number 388, ASV543
#> 2026-06-17 13:45:09.70 INFO::Fitting model to feature number 389, ASV544
#> 2026-06-17 13:45:09.71 INFO::Fitting model to feature number 390, ASV545
#> 2026-06-17 13:45:09.71 INFO::Fitting model to feature number 391, ASV546
#> 2026-06-17 13:45:09.72 INFO::Fitting model to feature number 392, ASV548
#> 2026-06-17 13:45:09.72 INFO::Fitting model to feature number 393, ASV549
#> 2026-06-17 13:45:09.73 INFO::Fitting model to feature number 394, ASV550
#> 2026-06-17 13:45:09.74 INFO::Fitting model to feature number 395, ASV551
#> 2026-06-17 13:45:09.74 INFO::Fitting model to feature number 396, ASV553
#> 2026-06-17 13:45:09.75 INFO::Fitting model to feature number 397, ASV554
#> 2026-06-17 13:45:09.75 INFO::Fitting model to feature number 398, ASV555
#> 2026-06-17 13:45:09.76 INFO::Fitting model to feature number 399, ASV556
#> 2026-06-17 13:45:09.77 INFO::Fitting model to feature number 400, ASV558
#> 2026-06-17 13:45:09.77 INFO::Fitting model to feature number 401, ASV559
#> 2026-06-17 13:45:09.78 INFO::Fitting model to feature number 402, ASV560
#> 2026-06-17 13:45:09.78 INFO::Fitting model to feature number 403, ASV561
#> 2026-06-17 13:45:09.79 INFO::Fitting model to feature number 404, ASV562
#> 2026-06-17 13:45:09.79 INFO::Fitting model to feature number 405, ASV563
#> 2026-06-17 13:45:09.80 INFO::Fitting model to feature number 406, ASV564
#> 2026-06-17 13:45:09.81 INFO::Fitting model to feature number 407, ASV566
#> 2026-06-17 13:45:09.81 INFO::Fitting model to feature number 408, ASV567
#> 2026-06-17 13:45:09.82 INFO::Fitting model to feature number 409, ASV569
#> 2026-06-17 13:45:09.82 INFO::Fitting model to feature number 410, ASV570
#> 2026-06-17 13:45:09.83 INFO::Fitting model to feature number 411, ASV571
#> 2026-06-17 13:45:09.84 INFO::Fitting model to feature number 412, ASV572
#> 2026-06-17 13:45:09.84 INFO::Fitting model to feature number 413, ASV573
#> 2026-06-17 13:45:09.85 INFO::Fitting model to feature number 414, ASV574
#> 2026-06-17 13:45:09.85 INFO::Fitting model to feature number 415, ASV576
#> 2026-06-17 13:45:09.86 INFO::Fitting model to feature number 416, ASV577
#> 2026-06-17 13:45:09.86 INFO::Fitting model to feature number 417, ASV578
#> 2026-06-17 13:45:09.87 INFO::Fitting model to feature number 418, ASV579
#> 2026-06-17 13:45:09.88 INFO::Fitting model to feature number 419, ASV580
#> 2026-06-17 13:45:09.88 INFO::Fitting model to feature number 420, ASV581
#> 2026-06-17 13:45:09.89 INFO::Fitting model to feature number 421, ASV582
#> 2026-06-17 13:45:09.89 INFO::Fitting model to feature number 422, ASV584
#> 2026-06-17 13:45:09.90 INFO::Fitting model to feature number 423, ASV585
#> 2026-06-17 13:45:09.91 INFO::Fitting model to feature number 424, ASV586
#> 2026-06-17 13:45:09.91 INFO::Fitting model to feature number 425, ASV587
#> 2026-06-17 13:45:09.92 INFO::Fitting model to feature number 426, ASV589
#> 2026-06-17 13:45:09.92 INFO::Fitting model to feature number 427, ASV590
#> 2026-06-17 13:45:09.93 INFO::Fitting model to feature number 428, ASV591
#> 2026-06-17 13:45:09.93 INFO::Fitting model to feature number 429, ASV592
#> 2026-06-17 13:45:09.94 INFO::Fitting model to feature number 430, ASV593
#> 2026-06-17 13:45:09.95 INFO::Fitting model to feature number 431, ASV594
#> 2026-06-17 13:45:09.95 INFO::Fitting model to feature number 432, ASV595
#> 2026-06-17 13:45:09.96 INFO::Fitting model to feature number 433, ASV596
#> 2026-06-17 13:45:09.96 INFO::Fitting model to feature number 434, ASV597
#> 2026-06-17 13:45:09.97 INFO::Fitting model to feature number 435, ASV598
#> 2026-06-17 13:45:09.98 INFO::Fitting model to feature number 436, ASV599
#> 2026-06-17 13:45:09.98 INFO::Fitting model to feature number 437, ASV600
#> 2026-06-17 13:45:09.99 INFO::Fitting model to feature number 438, ASV602
#> 2026-06-17 13:45:09.99 INFO::Fitting model to feature number 439, ASV604
#> 2026-06-17 13:45:10.00 INFO::Fitting model to feature number 440, ASV605
#> 2026-06-17 13:45:10.00 INFO::Fitting model to feature number 441, ASV607
#> 2026-06-17 13:45:10.01 INFO::Fitting model to feature number 442, ASV608
#> 2026-06-17 13:45:10.02 INFO::Fitting model to feature number 443, ASV610
#> 2026-06-17 13:45:10.02 INFO::Fitting model to feature number 444, ASV612
#> 2026-06-17 13:45:10.03 INFO::Fitting model to feature number 445, ASV614
#> 2026-06-17 13:45:10.03 INFO::Fitting model to feature number 446, ASV615
#> 2026-06-17 13:45:10.04 INFO::Fitting model to feature number 447, ASV616
#> 2026-06-17 13:45:10.04 INFO::Fitting model to feature number 448, ASV617
#> 2026-06-17 13:45:10.05 INFO::Fitting model to feature number 449, ASV618
#> 2026-06-17 13:45:10.06 INFO::Fitting model to feature number 450, ASV621
#> 2026-06-17 13:45:10.06 INFO::Fitting model to feature number 451, ASV624
#> 2026-06-17 13:45:10.07 INFO::Fitting model to feature number 452, ASV625
#> 2026-06-17 13:45:10.07 INFO::Fitting model to feature number 453, ASV626
#> 2026-06-17 13:45:10.08 INFO::Fitting model to feature number 454, ASV627
#> 2026-06-17 13:45:10.08 INFO::Fitting model to feature number 455, ASV628
#> 2026-06-17 13:45:10.09 INFO::Fitting model to feature number 456, ASV629
#> 2026-06-17 13:45:10.10 INFO::Fitting model to feature number 457, ASV631
#> 2026-06-17 13:45:10.10 INFO::Fitting model to feature number 458, ASV632
#> 2026-06-17 13:45:10.11 INFO::Fitting model to feature number 459, ASV633
#> 2026-06-17 13:45:10.11 INFO::Fitting model to feature number 460, ASV634
#> 2026-06-17 13:45:10.12 INFO::Fitting model to feature number 461, ASV635
#> 2026-06-17 13:45:10.12 INFO::Fitting model to feature number 462, ASV637
#> 2026-06-17 13:45:10.13 INFO::Fitting model to feature number 463, ASV639
#> 2026-06-17 13:45:10.14 INFO::Fitting model to feature number 464, ASV640
#> 2026-06-17 13:45:10.14 INFO::Fitting model to feature number 465, ASV641
#> 2026-06-17 13:45:10.15 INFO::Fitting model to feature number 466, ASV643
#> 2026-06-17 13:45:10.15 INFO::Fitting model to feature number 467, ASV644
#> 2026-06-17 13:45:10.16 INFO::Fitting model to feature number 468, ASV645
#> 2026-06-17 13:45:10.17 INFO::Fitting model to feature number 469, ASV647
#> 2026-06-17 13:45:10.17 INFO::Fitting model to feature number 470, ASV649
#> 2026-06-17 13:45:10.18 INFO::Fitting model to feature number 471, ASV650
#> 2026-06-17 13:45:10.18 INFO::Fitting model to feature number 472, ASV652
#> 2026-06-17 13:45:10.23 INFO::Fitting model to feature number 473, ASV653
#> 2026-06-17 13:45:10.23 INFO::Fitting model to feature number 474, ASV654
#> 2026-06-17 13:45:10.24 INFO::Fitting model to feature number 475, ASV657
#> 2026-06-17 13:45:10.24 INFO::Fitting model to feature number 476, ASV659
#> 2026-06-17 13:45:10.25 INFO::Fitting model to feature number 477, ASV660
#> 2026-06-17 13:45:10.26 INFO::Fitting model to feature number 478, ASV661
#> 2026-06-17 13:45:10.26 INFO::Fitting model to feature number 479, ASV662
#> 2026-06-17 13:45:10.27 INFO::Fitting model to feature number 480, ASV664
#> 2026-06-17 13:45:10.27 INFO::Fitting model to feature number 481, ASV665
#> 2026-06-17 13:45:10.28 INFO::Fitting model to feature number 482, ASV666
#> 2026-06-17 13:45:10.29 INFO::Fitting model to feature number 483, ASV667
#> 2026-06-17 13:45:10.29 INFO::Fitting model to feature number 484, ASV668
#> 2026-06-17 13:45:10.30 INFO::Fitting model to feature number 485, ASV669
#> 2026-06-17 13:45:10.30 INFO::Fitting model to feature number 486, ASV670
#> 2026-06-17 13:45:10.31 INFO::Fitting model to feature number 487, ASV671
#> 2026-06-17 13:45:10.32 INFO::Fitting model to feature number 488, ASV672
#> 2026-06-17 13:45:10.32 INFO::Fitting model to feature number 489, ASV673
#> 2026-06-17 13:45:10.33 INFO::Fitting model to feature number 490, ASV674
#> 2026-06-17 13:45:10.33 INFO::Fitting model to feature number 491, ASV675
#> 2026-06-17 13:45:10.34 INFO::Fitting model to feature number 492, ASV676
#> 2026-06-17 13:45:10.34 INFO::Fitting model to feature number 493, ASV677
#> 2026-06-17 13:45:10.35 INFO::Fitting model to feature number 494, ASV678
#> 2026-06-17 13:45:10.36 INFO::Fitting model to feature number 495, ASV679
#> 2026-06-17 13:45:10.36 INFO::Fitting model to feature number 496, ASV680
#> 2026-06-17 13:45:10.37 INFO::Fitting model to feature number 497, ASV683
#> 2026-06-17 13:45:10.37 INFO::Fitting model to feature number 498, ASV684
#> 2026-06-17 13:45:10.38 INFO::Fitting model to feature number 499, ASV685
#> 2026-06-17 13:45:10.39 INFO::Fitting model to feature number 500, ASV686
#> 2026-06-17 13:45:10.39 INFO::Fitting model to feature number 501, ASV687
#> 2026-06-17 13:45:10.40 INFO::Fitting model to feature number 502, ASV688
#> 2026-06-17 13:45:10.40 INFO::Fitting model to feature number 503, ASV691
#> 2026-06-17 13:45:10.41 INFO::Fitting model to feature number 504, ASV692
#> 2026-06-17 13:45:10.42 INFO::Fitting model to feature number 505, ASV693
#> 2026-06-17 13:45:10.42 INFO::Fitting model to feature number 506, ASV694
#> 2026-06-17 13:45:10.43 INFO::Fitting model to feature number 507, ASV695
#> 2026-06-17 13:45:10.43 INFO::Fitting model to feature number 508, ASV696
#> 2026-06-17 13:45:10.44 INFO::Fitting model to feature number 509, ASV697
#> 2026-06-17 13:45:10.44 INFO::Fitting model to feature number 510, ASV698
#> 2026-06-17 13:45:10.45 INFO::Fitting model to feature number 511, ASV700
#> 2026-06-17 13:45:10.46 INFO::Fitting model to feature number 512, ASV702
#> 2026-06-17 13:45:10.46 INFO::Fitting model to feature number 513, ASV703
#> 2026-06-17 13:45:10.47 INFO::Fitting model to feature number 514, ASV704
#> 2026-06-17 13:45:10.47 INFO::Fitting model to feature number 515, ASV705
#> 2026-06-17 13:45:10.48 INFO::Fitting model to feature number 516, ASV707
#> 2026-06-17 13:45:10.49 INFO::Fitting model to feature number 517, ASV710
#> 2026-06-17 13:45:10.49 INFO::Fitting model to feature number 518, ASV711
#> 2026-06-17 13:45:10.50 INFO::Fitting model to feature number 519, ASV712
#> 2026-06-17 13:45:10.50 INFO::Fitting model to feature number 520, ASV713
#> 2026-06-17 13:45:10.51 INFO::Fitting model to feature number 521, ASV714
#> 2026-06-17 13:45:10.51 INFO::Fitting model to feature number 522, ASV715
#> 2026-06-17 13:45:10.52 INFO::Fitting model to feature number 523, ASV717
#> 2026-06-17 13:45:10.53 INFO::Fitting model to feature number 524, ASV718
#> 2026-06-17 13:45:10.53 INFO::Fitting model to feature number 525, ASV720
#> 2026-06-17 13:45:10.54 INFO::Fitting model to feature number 526, ASV722
#> 2026-06-17 13:45:10.54 INFO::Fitting model to feature number 527, ASV723
#> 2026-06-17 13:45:10.55 INFO::Fitting model to feature number 528, ASV724
#> 2026-06-17 13:45:10.55 INFO::Fitting model to feature number 529, ASV726
#> 2026-06-17 13:45:10.56 INFO::Fitting model to feature number 530, ASV727
#> 2026-06-17 13:45:10.57 INFO::Fitting model to feature number 531, ASV728
#> 2026-06-17 13:45:10.57 INFO::Fitting model to feature number 532, ASV729
#> 2026-06-17 13:45:10.58 INFO::Fitting model to feature number 533, ASV730
#> 2026-06-17 13:45:10.58 INFO::Fitting model to feature number 534, ASV731
#> 2026-06-17 13:45:10.59 INFO::Fitting model to feature number 535, ASV732
#> 2026-06-17 13:45:10.60 INFO::Fitting model to feature number 536, ASV733
#> 2026-06-17 13:45:10.60 INFO::Fitting model to feature number 537, ASV734
#> 2026-06-17 13:45:10.61 INFO::Fitting model to feature number 538, ASV736
#> 2026-06-17 13:45:10.61 INFO::Fitting model to feature number 539, ASV737
#> 2026-06-17 13:45:10.62 INFO::Fitting model to feature number 540, ASV738
#> 2026-06-17 13:45:10.62 INFO::Fitting model to feature number 541, ASV740
#> 2026-06-17 13:45:10.63 INFO::Fitting model to feature number 542, ASV741
#> 2026-06-17 13:45:10.64 INFO::Fitting model to feature number 543, ASV742
#> 2026-06-17 13:45:10.64 INFO::Fitting model to feature number 544, ASV743
#> 2026-06-17 13:45:10.65 INFO::Fitting model to feature number 545, ASV744
#> 2026-06-17 13:45:10.65 INFO::Fitting model to feature number 546, ASV746
#> 2026-06-17 13:45:10.66 INFO::Fitting model to feature number 547, ASV747
#> 2026-06-17 13:45:10.66 INFO::Fitting model to feature number 548, ASV748
#> 2026-06-17 13:45:10.67 INFO::Fitting model to feature number 549, ASV749
#> 2026-06-17 13:45:10.67 INFO::Fitting model to feature number 550, ASV752
#> 2026-06-17 13:45:10.68 INFO::Fitting model to feature number 551, ASV753
#> 2026-06-17 13:45:10.69 INFO::Fitting model to feature number 552, ASV754
#> 2026-06-17 13:45:10.69 INFO::Fitting model to feature number 553, ASV755
#> 2026-06-17 13:45:10.70 INFO::Fitting model to feature number 554, ASV756
#> 2026-06-17 13:45:10.70 INFO::Fitting model to feature number 555, ASV757
#> 2026-06-17 13:45:10.71 INFO::Fitting model to feature number 556, ASV758
#> 2026-06-17 13:45:10.71 INFO::Fitting model to feature number 557, ASV760
#> 2026-06-17 13:45:10.72 INFO::Fitting model to feature number 558, ASV762
#> 2026-06-17 13:45:10.73 INFO::Fitting model to feature number 559, ASV764
#> 2026-06-17 13:45:10.73 INFO::Fitting model to feature number 560, ASV766
#> 2026-06-17 13:45:10.74 INFO::Fitting model to feature number 561, ASV767
#> 2026-06-17 13:45:10.75 INFO::Fitting model to feature number 562, ASV768
#> 2026-06-17 13:45:10.76 INFO::Fitting model to feature number 563, ASV769
#> 2026-06-17 13:45:10.76 INFO::Fitting model to feature number 564, ASV770
#> 2026-06-17 13:45:10.77 INFO::Fitting model to feature number 565, ASV771
#> 2026-06-17 13:45:10.77 INFO::Fitting model to feature number 566, ASV772
#> 2026-06-17 13:45:10.78 INFO::Fitting model to feature number 567, ASV773
#> 2026-06-17 13:45:10.79 INFO::Fitting model to feature number 568, ASV774
#> 2026-06-17 13:45:10.79 INFO::Fitting model to feature number 569, ASV775
#> 2026-06-17 13:45:10.80 INFO::Fitting model to feature number 570, ASV776
#> 2026-06-17 13:45:10.80 INFO::Fitting model to feature number 571, ASV777
#> 2026-06-17 13:45:10.81 INFO::Fitting model to feature number 572, ASV779
#> 2026-06-17 13:45:10.81 INFO::Fitting model to feature number 573, ASV780
#> 2026-06-17 13:45:10.82 INFO::Fitting model to feature number 574, ASV781
#> 2026-06-17 13:45:10.83 INFO::Fitting model to feature number 575, ASV782
#> 2026-06-17 13:45:10.83 INFO::Fitting model to feature number 576, ASV783
#> 2026-06-17 13:45:10.84 INFO::Fitting model to feature number 577, ASV784
#> 2026-06-17 13:45:10.84 INFO::Fitting model to feature number 578, ASV785
#> 2026-06-17 13:45:10.85 INFO::Fitting model to feature number 579, ASV786
#> 2026-06-17 13:45:10.85 INFO::Fitting model to feature number 580, ASV787
#> 2026-06-17 13:45:10.86 INFO::Fitting model to feature number 581, ASV788
#> 2026-06-17 13:45:10.87 INFO::Fitting model to feature number 582, ASV790
#> 2026-06-17 13:45:10.87 INFO::Fitting model to feature number 583, ASV792
#> 2026-06-17 13:45:10.88 INFO::Fitting model to feature number 584, ASV795
#> 2026-06-17 13:45:10.88 INFO::Fitting model to feature number 585, ASV796
#> 2026-06-17 13:45:10.89 INFO::Fitting model to feature number 586, ASV797
#> 2026-06-17 13:45:10.90 INFO::Fitting model to feature number 587, ASV798
#> 2026-06-17 13:45:10.90 INFO::Fitting model to feature number 588, ASV799
#> 2026-06-17 13:45:10.91 INFO::Fitting model to feature number 589, ASV801
#> 2026-06-17 13:45:10.91 INFO::Fitting model to feature number 590, ASV802
#> 2026-06-17 13:45:10.92 INFO::Fitting model to feature number 591, ASV803
#> 2026-06-17 13:45:10.93 INFO::Fitting model to feature number 592, ASV805
#> 2026-06-17 13:45:10.93 INFO::Fitting model to feature number 593, ASV807
#> 2026-06-17 13:45:10.94 INFO::Fitting model to feature number 594, ASV808
#> 2026-06-17 13:45:10.94 INFO::Fitting model to feature number 595, ASV810
#> 2026-06-17 13:45:10.95 INFO::Fitting model to feature number 596, ASV811
#> 2026-06-17 13:45:10.96 INFO::Fitting model to feature number 597, ASV814
#> 2026-06-17 13:45:10.96 INFO::Fitting model to feature number 598, ASV815
#> 2026-06-17 13:45:10.97 INFO::Fitting model to feature number 599, ASV816
#> 2026-06-17 13:45:10.97 INFO::Fitting model to feature number 600, ASV817
#> 2026-06-17 13:45:10.98 INFO::Fitting model to feature number 601, ASV819
#> 2026-06-17 13:45:10.98 INFO::Fitting model to feature number 602, ASV821
#> 2026-06-17 13:45:11.03 INFO::Fitting model to feature number 603, ASV822
#> 2026-06-17 13:45:11.03 INFO::Fitting model to feature number 604, ASV823
#> 2026-06-17 13:45:11.04 INFO::Fitting model to feature number 605, ASV824
#> 2026-06-17 13:45:11.04 INFO::Fitting model to feature number 606, ASV828
#> 2026-06-17 13:45:11.05 INFO::Fitting model to feature number 607, ASV829
#> 2026-06-17 13:45:11.06 INFO::Fitting model to feature number 608, ASV830
#> 2026-06-17 13:45:11.06 INFO::Fitting model to feature number 609, ASV831
#> 2026-06-17 13:45:11.07 INFO::Fitting model to feature number 610, ASV832
#> 2026-06-17 13:45:11.07 INFO::Fitting model to feature number 611, ASV834
#> 2026-06-17 13:45:11.08 INFO::Fitting model to feature number 612, ASV836
#> 2026-06-17 13:45:11.08 INFO::Fitting model to feature number 613, ASV837
#> 2026-06-17 13:45:11.09 INFO::Fitting model to feature number 614, ASV838
#> 2026-06-17 13:45:11.10 INFO::Fitting model to feature number 615, ASV839
#> 2026-06-17 13:45:11.10 INFO::Fitting model to feature number 616, ASV840
#> 2026-06-17 13:45:11.11 INFO::Fitting model to feature number 617, ASV841
#> 2026-06-17 13:45:11.11 INFO::Fitting model to feature number 618, ASV842
#> 2026-06-17 13:45:11.12 INFO::Fitting model to feature number 619, ASV843
#> 2026-06-17 13:45:11.12 INFO::Fitting model to feature number 620, ASV844
#> 2026-06-17 13:45:11.13 INFO::Fitting model to feature number 621, ASV845
#> 2026-06-17 13:45:11.14 INFO::Fitting model to feature number 622, ASV847
#> 2026-06-17 13:45:11.14 INFO::Fitting model to feature number 623, ASV848
#> 2026-06-17 13:45:11.15 INFO::Fitting model to feature number 624, ASV852
#> 2026-06-17 13:45:11.15 INFO::Fitting model to feature number 625, ASV853
#> 2026-06-17 13:45:11.16 INFO::Fitting model to feature number 626, ASV854
#> 2026-06-17 13:45:11.17 INFO::Fitting model to feature number 627, ASV855
#> 2026-06-17 13:45:11.17 INFO::Fitting model to feature number 628, ASV857
#> 2026-06-17 13:45:11.18 INFO::Fitting model to feature number 629, ASV858
#> 2026-06-17 13:45:11.18 INFO::Fitting model to feature number 630, ASV859
#> 2026-06-17 13:45:11.19 INFO::Fitting model to feature number 631, ASV860
#> 2026-06-17 13:45:11.19 INFO::Fitting model to feature number 632, ASV861
#> 2026-06-17 13:45:11.20 INFO::Fitting model to feature number 633, ASV863
#> 2026-06-17 13:45:11.21 INFO::Fitting model to feature number 634, ASV865
#> 2026-06-17 13:45:11.21 INFO::Fitting model to feature number 635, ASV870
#> 2026-06-17 13:45:11.22 INFO::Fitting model to feature number 636, ASV873
#> 2026-06-17 13:45:11.22 INFO::Fitting model to feature number 637, ASV874
#> 2026-06-17 13:45:11.23 INFO::Fitting model to feature number 638, ASV875
#> 2026-06-17 13:45:11.23 INFO::Fitting model to feature number 639, ASV876
#> 2026-06-17 13:45:11.24 INFO::Fitting model to feature number 640, ASV877
#> 2026-06-17 13:45:11.25 INFO::Fitting model to feature number 641, ASV878
#> 2026-06-17 13:45:11.25 INFO::Fitting model to feature number 642, ASV879
#> 2026-06-17 13:45:11.26 INFO::Fitting model to feature number 643, ASV880
#> 2026-06-17 13:45:11.26 INFO::Fitting model to feature number 644, ASV883
#> 2026-06-17 13:45:11.27 INFO::Fitting model to feature number 645, ASV884
#> 2026-06-17 13:45:11.28 INFO::Fitting model to feature number 646, ASV887
#> 2026-06-17 13:45:11.28 INFO::Fitting model to feature number 647, ASV888
#> 2026-06-17 13:45:11.29 INFO::Fitting model to feature number 648, ASV890
#> 2026-06-17 13:45:11.29 INFO::Fitting model to feature number 649, ASV891
#> 2026-06-17 13:45:11.30 INFO::Fitting model to feature number 650, ASV892
#> 2026-06-17 13:45:11.30 INFO::Fitting model to feature number 651, ASV895
#> 2026-06-17 13:45:11.31 INFO::Fitting model to feature number 652, ASV897
#> 2026-06-17 13:45:11.32 INFO::Fitting model to feature number 653, ASV898
#> 2026-06-17 13:45:11.32 INFO::Fitting model to feature number 654, ASV899
#> 2026-06-17 13:45:11.33 INFO::Fitting model to feature number 655, ASV901
#> 2026-06-17 13:45:11.33 INFO::Fitting model to feature number 656, ASV902
#> 2026-06-17 13:45:11.34 INFO::Fitting model to feature number 657, ASV903
#> 2026-06-17 13:45:11.35 INFO::Fitting model to feature number 658, ASV904
#> 2026-06-17 13:45:11.35 INFO::Fitting model to feature number 659, ASV905
#> 2026-06-17 13:45:11.36 INFO::Fitting model to feature number 660, ASV906
#> 2026-06-17 13:45:11.36 INFO::Fitting model to feature number 661, ASV907
#> 2026-06-17 13:45:11.37 INFO::Fitting model to feature number 662, ASV909
#> 2026-06-17 13:45:11.37 INFO::Fitting model to feature number 663, ASV910
#> 2026-06-17 13:45:11.38 INFO::Fitting model to feature number 664, ASV911
#> 2026-06-17 13:45:11.39 INFO::Fitting model to feature number 665, ASV913
#> 2026-06-17 13:45:11.39 INFO::Fitting model to feature number 666, ASV914
#> 2026-06-17 13:45:11.40 INFO::Fitting model to feature number 667, ASV915
#> 2026-06-17 13:45:11.40 INFO::Fitting model to feature number 668, ASV916
#> 2026-06-17 13:45:11.41 INFO::Fitting model to feature number 669, ASV917
#> 2026-06-17 13:45:11.42 INFO::Fitting model to feature number 670, ASV918
#> 2026-06-17 13:45:11.42 INFO::Fitting model to feature number 671, ASV919
#> 2026-06-17 13:45:11.43 INFO::Fitting model to feature number 672, ASV920
#> 2026-06-17 13:45:11.43 INFO::Fitting model to feature number 673, ASV921
#> 2026-06-17 13:45:11.44 INFO::Fitting model to feature number 674, ASV922
#> 2026-06-17 13:45:11.44 INFO::Fitting model to feature number 675, ASV923
#> 2026-06-17 13:45:11.45 INFO::Fitting model to feature number 676, ASV926
#> 2026-06-17 13:45:11.46 INFO::Fitting model to feature number 677, ASV927
#> 2026-06-17 13:45:11.46 INFO::Fitting model to feature number 678, ASV928
#> 2026-06-17 13:45:11.47 INFO::Fitting model to feature number 679, ASV930
#> 2026-06-17 13:45:11.47 INFO::Fitting model to feature number 680, ASV932
#> 2026-06-17 13:45:11.48 INFO::Fitting model to feature number 681, ASV934
#> 2026-06-17 13:45:11.49 INFO::Fitting model to feature number 682, ASV935
#> 2026-06-17 13:45:11.49 INFO::Fitting model to feature number 683, ASV936
#> 2026-06-17 13:45:11.50 INFO::Fitting model to feature number 684, ASV937
#> 2026-06-17 13:45:11.50 INFO::Fitting model to feature number 685, ASV938
#> 2026-06-17 13:45:11.51 INFO::Fitting model to feature number 686, ASV939
#> 2026-06-17 13:45:11.52 INFO::Fitting model to feature number 687, ASV940
#> 2026-06-17 13:45:11.52 INFO::Fitting model to feature number 688, ASV941
#> 2026-06-17 13:45:11.53 INFO::Fitting model to feature number 689, ASV942
#> 2026-06-17 13:45:11.54 INFO::Fitting model to feature number 690, ASV943
#> 2026-06-17 13:45:11.54 INFO::Fitting model to feature number 691, ASV945
#> 2026-06-17 13:45:11.55 INFO::Fitting model to feature number 692, ASV947
#> 2026-06-17 13:45:11.56 INFO::Fitting model to feature number 693, ASV948
#> 2026-06-17 13:45:11.56 INFO::Fitting model to feature number 694, ASV949
#> 2026-06-17 13:45:11.57 INFO::Fitting model to feature number 695, ASV950
#> 2026-06-17 13:45:11.58 INFO::Fitting model to feature number 696, ASV951
#> 2026-06-17 13:45:11.58 INFO::Fitting model to feature number 697, ASV953
#> 2026-06-17 13:45:11.59 INFO::Fitting model to feature number 698, ASV955
#> 2026-06-17 13:45:11.60 INFO::Fitting model to feature number 699, ASV958
#> 2026-06-17 13:45:11.60 INFO::Fitting model to feature number 700, ASV959
#> 2026-06-17 13:45:11.61 INFO::Fitting model to feature number 701, ASV961
#> 2026-06-17 13:45:11.61 INFO::Fitting model to feature number 702, ASV962
#> 2026-06-17 13:45:11.62 INFO::Fitting model to feature number 703, ASV963
#> 2026-06-17 13:45:11.62 INFO::Fitting model to feature number 704, ASV964
#> 2026-06-17 13:45:11.63 INFO::Fitting model to feature number 705, ASV966
#> 2026-06-17 13:45:11.64 INFO::Fitting model to feature number 706, ASV967
#> 2026-06-17 13:45:11.64 INFO::Fitting model to feature number 707, ASV969
#> 2026-06-17 13:45:11.65 INFO::Fitting model to feature number 708, ASV970
#> 2026-06-17 13:45:11.65 INFO::Fitting model to feature number 709, ASV971
#> 2026-06-17 13:45:11.66 INFO::Fitting model to feature number 710, ASV972
#> 2026-06-17 13:45:11.67 INFO::Fitting model to feature number 711, ASV973
#> 2026-06-17 13:45:11.67 INFO::Fitting model to feature number 712, ASV974
#> 2026-06-17 13:45:11.68 INFO::Fitting model to feature number 713, ASV975
#> 2026-06-17 13:45:11.68 INFO::Fitting model to feature number 714, ASV976
#> 2026-06-17 13:45:11.69 INFO::Fitting model to feature number 715, ASV977
#> 2026-06-17 13:45:11.69 INFO::Fitting model to feature number 716, ASV979
#> 2026-06-17 13:45:11.70 INFO::Fitting model to feature number 717, ASV980
#> 2026-06-17 13:45:11.71 INFO::Fitting model to feature number 718, ASV981
#> 2026-06-17 13:45:11.71 INFO::Fitting model to feature number 719, ASV983
#> 2026-06-17 13:45:11.72 INFO::Fitting model to feature number 720, ASV984
#> 2026-06-17 13:45:11.72 INFO::Fitting model to feature number 721, ASV986
#> 2026-06-17 13:45:12.47 INFO::Fitting model to feature number 722, ASV987
#> 2026-06-17 13:45:12.48 INFO::Fitting model to feature number 723, ASV988
#> 2026-06-17 13:45:12.48 INFO::Fitting model to feature number 724, ASV989
#> 2026-06-17 13:45:12.49 INFO::Fitting model to feature number 725, ASV990
#> 2026-06-17 13:45:12.49 INFO::Fitting model to feature number 726, ASV992
#> 2026-06-17 13:45:12.50 INFO::Fitting model to feature number 727, ASV993
#> 2026-06-17 13:45:12.51 INFO::Fitting model to feature number 728, ASV994
#> 2026-06-17 13:45:12.51 INFO::Fitting model to feature number 729, ASV995
#> 2026-06-17 13:45:12.52 INFO::Fitting model to feature number 730, ASV996
#> 2026-06-17 13:45:12.52 INFO::Fitting model to feature number 731, ASV997
#> 2026-06-17 13:45:12.53 INFO::Fitting model to feature number 732, ASV998
#> 2026-06-17 13:45:12.54 INFO::Fitting model to feature number 733, ASV999
#> 2026-06-17 13:45:12.54 INFO::Fitting model to feature number 734, ASV1001
#> 2026-06-17 13:45:12.55 INFO::Fitting model to feature number 735, ASV1002
#> 2026-06-17 13:45:12.55 INFO::Fitting model to feature number 736, ASV1003
#> 2026-06-17 13:45:12.56 INFO::Fitting model to feature number 737, ASV1004
#> 2026-06-17 13:45:12.57 INFO::Fitting model to feature number 738, ASV1005
#> 2026-06-17 13:45:12.57 INFO::Fitting model to feature number 739, ASV1006
#> 2026-06-17 13:45:12.58 INFO::Fitting model to feature number 740, ASV1007
#> 2026-06-17 13:45:12.58 INFO::Fitting model to feature number 741, ASV1008
#> 2026-06-17 13:45:12.59 INFO::Fitting model to feature number 742, ASV1010
#> 2026-06-17 13:45:12.60 INFO::Fitting model to feature number 743, ASV1011
#> 2026-06-17 13:45:12.60 INFO::Fitting model to feature number 744, ASV1014
#> 2026-06-17 13:45:12.61 INFO::Fitting model to feature number 745, ASV1015
#> 2026-06-17 13:45:12.61 INFO::Fitting model to feature number 746, ASV1016
#> 2026-06-17 13:45:12.62 INFO::Fitting model to feature number 747, ASV1017
#> 2026-06-17 13:45:12.63 INFO::Fitting model to feature number 748, ASV1018
#> 2026-06-17 13:45:12.63 INFO::Fitting model to feature number 749, ASV1020
#> 2026-06-17 13:45:12.64 INFO::Fitting model to feature number 750, ASV1021
#> 2026-06-17 13:45:12.64 INFO::Fitting model to feature number 751, ASV1022
#> 2026-06-17 13:45:12.65 INFO::Fitting model to feature number 752, ASV1023
#> 2026-06-17 13:45:12.65 INFO::Fitting model to feature number 753, ASV1024
#> 2026-06-17 13:45:12.66 INFO::Fitting model to feature number 754, ASV1025
#> 2026-06-17 13:45:12.67 INFO::Fitting model to feature number 755, ASV1026
#> 2026-06-17 13:45:12.67 INFO::Fitting model to feature number 756, ASV1027
#> 2026-06-17 13:45:12.68 INFO::Fitting model to feature number 757, ASV1028
#> 2026-06-17 13:45:12.68 INFO::Fitting model to feature number 758, ASV1029
#> 2026-06-17 13:45:12.69 INFO::Fitting model to feature number 759, ASV1030
#> 2026-06-17 13:45:12.70 INFO::Fitting model to feature number 760, ASV1031
#> 2026-06-17 13:45:12.70 INFO::Fitting model to feature number 761, ASV1032
#> 2026-06-17 13:45:12.71 INFO::Fitting model to feature number 762, ASV1035
#> 2026-06-17 13:45:12.71 INFO::Fitting model to feature number 763, ASV1036
#> 2026-06-17 13:45:12.72 INFO::Fitting model to feature number 764, ASV1037
#> 2026-06-17 13:45:12.72 INFO::Fitting model to feature number 765, ASV1038
#> 2026-06-17 13:45:12.73 INFO::Fitting model to feature number 766, ASV1039
#> 2026-06-17 13:45:12.74 INFO::Fitting model to feature number 767, ASV1040
#> 2026-06-17 13:45:12.74 INFO::Fitting model to feature number 768, ASV1041
#> 2026-06-17 13:45:12.75 INFO::Fitting model to feature number 769, ASV1042
#> 2026-06-17 13:45:12.75 INFO::Fitting model to feature number 770, ASV1044
#> 2026-06-17 13:45:12.76 INFO::Fitting model to feature number 771, ASV1046
#> 2026-06-17 13:45:12.77 INFO::Fitting model to feature number 772, ASV1047
#> 2026-06-17 13:45:12.77 INFO::Fitting model to feature number 773, ASV1048
#> 2026-06-17 13:45:12.78 INFO::Fitting model to feature number 774, ASV1049
#> 2026-06-17 13:45:12.78 INFO::Fitting model to feature number 775, ASV1052
#> 2026-06-17 13:45:12.79 INFO::Fitting model to feature number 776, ASV1053
#> 2026-06-17 13:45:12.80 INFO::Fitting model to feature number 777, ASV1057
#> 2026-06-17 13:45:12.80 INFO::Fitting model to feature number 778, ASV1058
#> 2026-06-17 13:45:12.81 INFO::Fitting model to feature number 779, ASV1059
#> 2026-06-17 13:45:12.81 INFO::Fitting model to feature number 780, ASV1063
#> 2026-06-17 13:45:12.82 INFO::Fitting model to feature number 781, ASV1065
#> 2026-06-17 13:45:12.82 INFO::Fitting model to feature number 782, ASV1066
#> 2026-06-17 13:45:12.83 INFO::Fitting model to feature number 783, ASV1067
#> 2026-06-17 13:45:12.84 INFO::Fitting model to feature number 784, ASV1068
#> 2026-06-17 13:45:12.84 INFO::Fitting model to feature number 785, ASV1069
#> 2026-06-17 13:45:12.85 INFO::Fitting model to feature number 786, ASV1070
#> 2026-06-17 13:45:12.85 INFO::Fitting model to feature number 787, ASV1071
#> 2026-06-17 13:45:12.86 INFO::Fitting model to feature number 788, ASV1072
#> 2026-06-17 13:45:12.86 INFO::Fitting model to feature number 789, ASV1073
#> 2026-06-17 13:45:12.87 INFO::Fitting model to feature number 790, ASV1074
#> 2026-06-17 13:45:12.88 INFO::Fitting model to feature number 791, ASV1076
#> 2026-06-17 13:45:12.88 INFO::Fitting model to feature number 792, ASV1078
#> 2026-06-17 13:45:12.89 INFO::Fitting model to feature number 793, ASV1079
#> 2026-06-17 13:45:12.89 INFO::Fitting model to feature number 794, ASV1082
#> 2026-06-17 13:45:12.90 INFO::Fitting model to feature number 795, ASV1084
#> 2026-06-17 13:45:12.91 INFO::Fitting model to feature number 796, ASV1085
#> 2026-06-17 13:45:12.91 INFO::Fitting model to feature number 797, ASV1086
#> 2026-06-17 13:45:12.92 INFO::Fitting model to feature number 798, ASV1087
#> 2026-06-17 13:45:12.92 INFO::Fitting model to feature number 799, ASV1090
#> 2026-06-17 13:45:12.93 INFO::Fitting model to feature number 800, ASV1091
#> 2026-06-17 13:45:12.93 INFO::Fitting model to feature number 801, ASV1093
#> 2026-06-17 13:45:12.94 INFO::Fitting model to feature number 802, ASV1095
#> 2026-06-17 13:45:12.95 INFO::Fitting model to feature number 803, ASV1096
#> 2026-06-17 13:45:12.95 INFO::Fitting model to feature number 804, ASV1097
#> 2026-06-17 13:45:12.96 INFO::Fitting model to feature number 805, ASV1099
#> 2026-06-17 13:45:12.96 INFO::Fitting model to feature number 806, ASV1100
#> 2026-06-17 13:45:12.97 INFO::Fitting model to feature number 807, ASV1101
#> 2026-06-17 13:45:12.98 INFO::Fitting model to feature number 808, ASV1103
#> 2026-06-17 13:45:12.98 INFO::Fitting model to feature number 809, ASV1105
#> 2026-06-17 13:45:12.99 INFO::Fitting model to feature number 810, ASV1107
#> 2026-06-17 13:45:12.99 INFO::Fitting model to feature number 811, ASV1108
#> 2026-06-17 13:45:13.00 INFO::Fitting model to feature number 812, ASV1109
#> 2026-06-17 13:45:13.01 INFO::Fitting model to feature number 813, ASV1110
#> 2026-06-17 13:45:13.01 INFO::Fitting model to feature number 814, ASV1111
#> 2026-06-17 13:45:13.02 INFO::Fitting model to feature number 815, ASV1112
#> 2026-06-17 13:45:13.02 INFO::Fitting model to feature number 816, ASV1114
#> 2026-06-17 13:45:13.03 INFO::Fitting model to feature number 817, ASV1115
#> 2026-06-17 13:45:13.03 INFO::Fitting model to feature number 818, ASV1116
#> 2026-06-17 13:45:13.04 INFO::Fitting model to feature number 819, ASV1120
#> 2026-06-17 13:45:13.05 INFO::Fitting model to feature number 820, ASV1121
#> 2026-06-17 13:45:13.05 INFO::Fitting model to feature number 821, ASV1122
#> 2026-06-17 13:45:13.06 INFO::Fitting model to feature number 822, ASV1125
#> 2026-06-17 13:45:13.07 INFO::Fitting model to feature number 823, ASV1126
#> 2026-06-17 13:45:13.07 INFO::Fitting model to feature number 824, ASV1127
#> 2026-06-17 13:45:13.08 INFO::Fitting model to feature number 825, ASV1128
#> 2026-06-17 13:45:13.08 INFO::Fitting model to feature number 826, ASV1129
#> 2026-06-17 13:45:13.09 INFO::Fitting model to feature number 827, ASV1131
#> 2026-06-17 13:45:13.09 INFO::Fitting model to feature number 828, ASV1132
#> 2026-06-17 13:45:13.10 INFO::Fitting model to feature number 829, ASV1133
#> 2026-06-17 13:45:13.11 INFO::Fitting model to feature number 830, ASV1134
#> 2026-06-17 13:45:13.11 INFO::Fitting model to feature number 831, ASV1135
#> 2026-06-17 13:45:13.12 INFO::Fitting model to feature number 832, ASV1137
#> 2026-06-17 13:45:13.12 INFO::Fitting model to feature number 833, ASV1138
#> 2026-06-17 13:45:13.13 INFO::Fitting model to feature number 834, ASV1139
#> 2026-06-17 13:45:13.14 INFO::Fitting model to feature number 835, ASV1141
#> 2026-06-17 13:45:13.14 INFO::Fitting model to feature number 836, ASV1143
#> 2026-06-17 13:45:13.15 INFO::Fitting model to feature number 837, ASV1144
#> 2026-06-17 13:45:13.15 INFO::Fitting model to feature number 838, ASV1146
#> 2026-06-17 13:45:13.16 INFO::Fitting model to feature number 839, ASV1147
#> 2026-06-17 13:45:13.17 INFO::Fitting model to feature number 840, ASV1148
#> 2026-06-17 13:45:13.17 INFO::Fitting model to feature number 841, ASV1150
#> 2026-06-17 13:45:13.18 INFO::Fitting model to feature number 842, ASV1151
#> 2026-06-17 13:45:13.18 INFO::Fitting model to feature number 843, ASV1152
#> 2026-06-17 13:45:13.19 INFO::Fitting model to feature number 844, ASV1154
#> 2026-06-17 13:45:13.20 INFO::Fitting model to feature number 845, ASV1155
#> 2026-06-17 13:45:13.20 INFO::Fitting model to feature number 846, ASV1156
#> 2026-06-17 13:45:13.21 INFO::Fitting model to feature number 847, ASV1158
#> 2026-06-17 13:45:13.21 INFO::Fitting model to feature number 848, ASV1159
#> 2026-06-17 13:45:13.22 INFO::Fitting model to feature number 849, ASV1160
#> 2026-06-17 13:45:13.23 INFO::Fitting model to feature number 850, ASV1161
#> 2026-06-17 13:45:13.23 INFO::Fitting model to feature number 851, ASV1162
#> 2026-06-17 13:45:13.24 INFO::Fitting model to feature number 852, ASV1163
#> 2026-06-17 13:45:13.24 INFO::Fitting model to feature number 853, ASV1164
#> 2026-06-17 13:45:13.25 INFO::Fitting model to feature number 854, ASV1165
#> 2026-06-17 13:45:13.26 INFO::Fitting model to feature number 855, ASV1167
#> 2026-06-17 13:45:13.26 INFO::Fitting model to feature number 856, ASV1168
#> 2026-06-17 13:45:13.27 INFO::Fitting model to feature number 857, ASV1169
#> 2026-06-17 13:45:13.27 INFO::Fitting model to feature number 858, ASV1171
#> 2026-06-17 13:45:13.28 INFO::Fitting model to feature number 859, ASV1172
#> 2026-06-17 13:45:13.29 INFO::Fitting model to feature number 860, ASV1173
#> 2026-06-17 13:45:13.29 INFO::Fitting model to feature number 861, ASV1175
#> 2026-06-17 13:45:13.30 INFO::Fitting model to feature number 862, ASV1176
#> 2026-06-17 13:45:13.30 INFO::Fitting model to feature number 863, ASV1177
#> 2026-06-17 13:45:13.31 INFO::Fitting model to feature number 864, ASV1179
#> 2026-06-17 13:45:13.32 INFO::Fitting model to feature number 865, ASV1180
#> 2026-06-17 13:45:13.32 INFO::Fitting model to feature number 866, ASV1182
#> 2026-06-17 13:45:13.33 INFO::Fitting model to feature number 867, ASV1184
#> 2026-06-17 13:45:13.33 INFO::Fitting model to feature number 868, ASV1185
#> 2026-06-17 13:45:13.34 INFO::Fitting model to feature number 869, ASV1186
#> 2026-06-17 13:45:13.35 INFO::Fitting model to feature number 870, ASV1187
#> 2026-06-17 13:45:13.35 INFO::Fitting model to feature number 871, ASV1189
#> 2026-06-17 13:45:13.36 INFO::Fitting model to feature number 872, ASV1190
#> 2026-06-17 13:45:13.36 INFO::Fitting model to feature number 873, ASV1192
#> 2026-06-17 13:45:13.37 INFO::Fitting model to feature number 874, ASV1193
#> 2026-06-17 13:45:13.38 INFO::Fitting model to feature number 875, ASV1194
#> 2026-06-17 13:45:13.38 INFO::Fitting model to feature number 876, ASV1195
#> 2026-06-17 13:45:13.39 INFO::Fitting model to feature number 877, ASV1198
#> 2026-06-17 13:45:13.39 INFO::Fitting model to feature number 878, ASV1199
#> 2026-06-17 13:45:13.40 INFO::Fitting model to feature number 879, ASV1200
#> 2026-06-17 13:45:13.41 INFO::Fitting model to feature number 880, ASV1203
#> 2026-06-17 13:45:13.41 INFO::Fitting model to feature number 881, ASV1204
#> 2026-06-17 13:45:13.42 INFO::Fitting model to feature number 882, ASV1205
#> 2026-06-17 13:45:13.42 INFO::Fitting model to feature number 883, ASV1206
#> 2026-06-17 13:45:13.43 INFO::Fitting model to feature number 884, ASV1208
#> 2026-06-17 13:45:13.44 INFO::Fitting model to feature number 885, ASV1209
#> 2026-06-17 13:45:13.44 INFO::Fitting model to feature number 886, ASV1210
#> 2026-06-17 13:45:13.45 INFO::Fitting model to feature number 887, ASV1211
#> 2026-06-17 13:45:13.46 INFO::Fitting model to feature number 888, ASV1212
#> 2026-06-17 13:45:13.46 INFO::Fitting model to feature number 889, ASV1213
#> 2026-06-17 13:45:13.47 INFO::Fitting model to feature number 890, ASV1214
#> 2026-06-17 13:45:13.47 INFO::Fitting model to feature number 891, ASV1216
#> 2026-06-17 13:45:13.48 INFO::Fitting model to feature number 892, ASV1217
#> 2026-06-17 13:45:13.49 INFO::Fitting model to feature number 893, ASV1218
#> 2026-06-17 13:45:13.49 INFO::Fitting model to feature number 894, ASV1219
#> 2026-06-17 13:45:13.50 INFO::Fitting model to feature number 895, ASV1221
#> 2026-06-17 13:45:13.51 INFO::Fitting model to feature number 896, ASV1223
#> 2026-06-17 13:45:13.51 INFO::Fitting model to feature number 897, ASV1224
#> 2026-06-17 13:45:13.52 INFO::Fitting model to feature number 898, ASV1225
#> 2026-06-17 13:45:13.52 INFO::Fitting model to feature number 899, ASV1227
#> 2026-06-17 13:45:13.53 INFO::Fitting model to feature number 900, ASV1228
#> 2026-06-17 13:45:13.53 INFO::Fitting model to feature number 901, ASV1229
#> 2026-06-17 13:45:13.54 INFO::Fitting model to feature number 902, ASV1230
#> 2026-06-17 13:45:13.55 INFO::Fitting model to feature number 903, ASV1231
#> 2026-06-17 13:45:13.55 INFO::Fitting model to feature number 904, ASV1232
#> 2026-06-17 13:45:13.56 INFO::Fitting model to feature number 905, ASV1233
#> 2026-06-17 13:45:13.57 INFO::Fitting model to feature number 906, ASV1234
#> 2026-06-17 13:45:13.57 INFO::Fitting model to feature number 907, ASV1236
#> 2026-06-17 13:45:13.58 INFO::Fitting model to feature number 908, ASV1238
#> 2026-06-17 13:45:13.58 INFO::Fitting model to feature number 909, ASV1239
#> 2026-06-17 13:45:13.59 INFO::Fitting model to feature number 910, ASV1241
#> 2026-06-17 13:45:13.60 INFO::Fitting model to feature number 911, ASV1242
#> 2026-06-17 13:45:13.60 INFO::Fitting model to feature number 912, ASV1243
#> 2026-06-17 13:45:13.61 INFO::Fitting model to feature number 913, ASV1245
#> 2026-06-17 13:45:13.61 INFO::Fitting model to feature number 914, ASV1246
#> 2026-06-17 13:45:13.62 INFO::Fitting model to feature number 915, ASV1247
#> 2026-06-17 13:45:13.62 INFO::Fitting model to feature number 916, ASV1251
#> 2026-06-17 13:45:13.63 INFO::Fitting model to feature number 917, ASV1252
#> 2026-06-17 13:45:13.64 INFO::Fitting model to feature number 918, ASV1253
#> 2026-06-17 13:45:13.64 INFO::Fitting model to feature number 919, ASV1254
#> 2026-06-17 13:45:13.65 INFO::Fitting model to feature number 920, ASV1257
#> 2026-06-17 13:45:13.65 INFO::Fitting model to feature number 921, ASV1258
#> 2026-06-17 13:45:13.66 INFO::Fitting model to feature number 922, ASV1259
#> 2026-06-17 13:45:13.66 INFO::Fitting model to feature number 923, ASV1260
#> 2026-06-17 13:45:13.67 INFO::Fitting model to feature number 924, ASV1261
#> 2026-06-17 13:45:13.68 INFO::Fitting model to feature number 925, ASV1262
#> 2026-06-17 13:45:13.68 INFO::Fitting model to feature number 926, ASV1263
#> 2026-06-17 13:45:13.69 INFO::Fitting model to feature number 927, ASV1264
#> 2026-06-17 13:45:13.69 INFO::Fitting model to feature number 928, ASV1265
#> 2026-06-17 13:45:13.70 INFO::Fitting model to feature number 929, ASV1267
#> 2026-06-17 13:45:13.71 INFO::Fitting model to feature number 930, ASV1268
#> 2026-06-17 13:45:13.71 INFO::Fitting model to feature number 931, ASV1269
#> 2026-06-17 13:45:13.72 INFO::Fitting model to feature number 932, ASV1270
#> 2026-06-17 13:45:13.72 INFO::Fitting model to feature number 933, ASV1271
#> 2026-06-17 13:45:13.73 INFO::Fitting model to feature number 934, ASV1272
#> 2026-06-17 13:45:13.73 INFO::Fitting model to feature number 935, ASV1273
#> 2026-06-17 13:45:13.74 INFO::Fitting model to feature number 936, ASV1274
#> 2026-06-17 13:45:13.75 INFO::Fitting model to feature number 937, ASV1275
#> 2026-06-17 13:45:13.75 INFO::Fitting model to feature number 938, ASV1276
#> 2026-06-17 13:45:13.76 INFO::Fitting model to feature number 939, ASV1278
#> 2026-06-17 13:45:13.77 INFO::Fitting model to feature number 940, ASV1279
#> 2026-06-17 13:45:13.77 INFO::Fitting model to feature number 941, ASV1282
#> 2026-06-17 13:45:13.78 INFO::Fitting model to feature number 942, ASV1283
#> 2026-06-17 13:45:13.78 INFO::Fitting model to feature number 943, ASV1284
#> 2026-06-17 13:45:13.79 INFO::Fitting model to feature number 944, ASV1285
#> 2026-06-17 13:45:13.80 INFO::Fitting model to feature number 945, ASV1286
#> 2026-06-17 13:45:13.80 INFO::Fitting model to feature number 946, ASV1287
#> 2026-06-17 13:45:13.81 INFO::Fitting model to feature number 947, ASV1288
#> 2026-06-17 13:45:13.81 INFO::Fitting model to feature number 948, ASV1289
#> 2026-06-17 13:45:13.82 INFO::Fitting model to feature number 949, ASV1290
#> 2026-06-17 13:45:13.82 INFO::Fitting model to feature number 950, ASV1293
#> 2026-06-17 13:45:13.83 INFO::Fitting model to feature number 951, ASV1294
#> 2026-06-17 13:45:13.84 INFO::Fitting model to feature number 952, ASV1296
#> 2026-06-17 13:45:13.84 INFO::Fitting model to feature number 953, ASV1297
#> 2026-06-17 13:45:13.85 INFO::Fitting model to feature number 954, ASV1300
#> 2026-06-17 13:45:13.85 INFO::Fitting model to feature number 955, ASV1301
#> 2026-06-17 13:45:13.86 INFO::Fitting model to feature number 956, ASV1302
#> 2026-06-17 13:45:13.87 INFO::Fitting model to feature number 957, ASV1303
#> 2026-06-17 13:45:13.87 INFO::Fitting model to feature number 958, ASV1304
#> 2026-06-17 13:45:13.88 INFO::Fitting model to feature number 959, ASV1305
#> 2026-06-17 13:45:13.88 INFO::Fitting model to feature number 960, ASV1307
#> 2026-06-17 13:45:13.89 INFO::Fitting model to feature number 961, ASV1310
#> 2026-06-17 13:45:13.89 INFO::Fitting model to feature number 962, ASV1311
#> 2026-06-17 13:45:13.90 INFO::Fitting model to feature number 963, ASV1312
#> 2026-06-17 13:45:13.91 INFO::Fitting model to feature number 964, ASV1313
#> 2026-06-17 13:45:13.91 INFO::Fitting model to feature number 965, ASV1314
#> 2026-06-17 13:45:13.92 INFO::Fitting model to feature number 966, ASV1315
#> 2026-06-17 13:45:13.92 INFO::Fitting model to feature number 967, ASV1316
#> 2026-06-17 13:45:13.93 INFO::Fitting model to feature number 968, ASV1317
#> 2026-06-17 13:45:13.94 INFO::Fitting model to feature number 969, ASV1319
#> 2026-06-17 13:45:13.94 INFO::Fitting model to feature number 970, ASV1320
#> 2026-06-17 13:45:13.95 INFO::Fitting model to feature number 971, ASV1321
#> 2026-06-17 13:45:13.95 INFO::Fitting model to feature number 972, ASV1323
#> 2026-06-17 13:45:13.96 INFO::Fitting model to feature number 973, ASV1326
#> 2026-06-17 13:45:13.97 INFO::Fitting model to feature number 974, ASV1327
#> 2026-06-17 13:45:13.97 INFO::Fitting model to feature number 975, ASV1328
#> 2026-06-17 13:45:13.98 INFO::Fitting model to feature number 976, ASV1330
#> 2026-06-17 13:45:13.98 INFO::Fitting model to feature number 977, ASV1332
#> 2026-06-17 13:45:13.99 INFO::Fitting model to feature number 978, ASV1334
#> 2026-06-17 13:45:13.99 INFO::Fitting model to feature number 979, ASV1335
#> 2026-06-17 13:45:14.00 INFO::Fitting model to feature number 980, ASV1336
#> 2026-06-17 13:45:14.01 INFO::Fitting model to feature number 981, ASV1337
#> 2026-06-17 13:45:14.01 INFO::Fitting model to feature number 982, ASV1338
#> 2026-06-17 13:45:14.02 INFO::Fitting model to feature number 983, ASV1340
#> 2026-06-17 13:45:14.02 INFO::Fitting model to feature number 984, ASV1341
#> 2026-06-17 13:45:14.03 INFO::Fitting model to feature number 985, ASV1342
#> 2026-06-17 13:45:14.03 INFO::Fitting model to feature number 986, ASV1345
#> 2026-06-17 13:45:14.04 INFO::Fitting model to feature number 987, ASV1350
#> 2026-06-17 13:45:14.05 INFO::Fitting model to feature number 988, ASV1351
#> 2026-06-17 13:45:14.05 INFO::Fitting model to feature number 989, ASV1352
#> 2026-06-17 13:45:14.06 INFO::Fitting model to feature number 990, ASV1353
#> 2026-06-17 13:45:14.06 INFO::Fitting model to feature number 991, ASV1355
#> 2026-06-17 13:45:14.07 INFO::Fitting model to feature number 992, ASV1356
#> 2026-06-17 13:45:14.08 INFO::Fitting model to feature number 993, ASV1359
#> 2026-06-17 13:45:14.08 INFO::Fitting model to feature number 994, ASV1360
#> 2026-06-17 13:45:14.09 INFO::Fitting model to feature number 995, ASV1361
#> 2026-06-17 13:45:14.09 INFO::Fitting model to feature number 996, ASV1362
#> 2026-06-17 13:45:14.10 INFO::Fitting model to feature number 997, ASV1363
#> 2026-06-17 13:45:14.10 INFO::Fitting model to feature number 998, ASV1365
#> 2026-06-17 13:45:14.11 INFO::Fitting model to feature number 999, ASV1366
#> 2026-06-17 13:45:14.12 INFO::Fitting model to feature number 1000, ASV1367
#> 2026-06-17 13:45:14.12 INFO::Fitting model to feature number 1001, ASV1368
#> 2026-06-17 13:45:14.13 INFO::Fitting model to feature number 1002, ASV1369
#> 2026-06-17 13:45:14.13 INFO::Fitting model to feature number 1003, ASV1370
#> 2026-06-17 13:45:14.14 INFO::Fitting model to feature number 1004, ASV1371
#> 2026-06-17 13:45:14.14 INFO::Fitting model to feature number 1005, ASV1372
#> 2026-06-17 13:45:14.15 INFO::Fitting model to feature number 1006, ASV1373
#> 2026-06-17 13:45:14.16 INFO::Fitting model to feature number 1007, ASV1374
#> 2026-06-17 13:45:14.16 INFO::Fitting model to feature number 1008, ASV1375
#> 2026-06-17 13:45:14.17 INFO::Fitting model to feature number 1009, ASV1376
#> 2026-06-17 13:45:14.17 INFO::Fitting model to feature number 1010, ASV1378
#> 2026-06-17 13:45:14.18 INFO::Fitting model to feature number 1011, ASV1379
#> 2026-06-17 13:45:14.19 INFO::Fitting model to feature number 1012, ASV1380
#> 2026-06-17 13:45:14.19 INFO::Fitting model to feature number 1013, ASV1381
#> 2026-06-17 13:45:14.26 INFO::Fitting model to feature number 1014, ASV1384
#> 2026-06-17 13:45:14.26 INFO::Fitting model to feature number 1015, ASV1386
#> 2026-06-17 13:45:14.27 INFO::Fitting model to feature number 1016, ASV1387
#> 2026-06-17 13:45:14.28 INFO::Fitting model to feature number 1017, ASV1388
#> 2026-06-17 13:45:14.28 INFO::Fitting model to feature number 1018, ASV1389
#> 2026-06-17 13:45:14.29 INFO::Fitting model to feature number 1019, ASV1390
#> 2026-06-17 13:45:14.29 INFO::Fitting model to feature number 1020, ASV1392
#> 2026-06-17 13:45:14.30 INFO::Fitting model to feature number 1021, ASV1393
#> 2026-06-17 13:45:14.31 INFO::Fitting model to feature number 1022, ASV1396
#> 2026-06-17 13:45:14.31 INFO::Fitting model to feature number 1023, ASV1397
#> 2026-06-17 13:45:14.32 INFO::Fitting model to feature number 1024, ASV1398
#> 2026-06-17 13:45:14.32 INFO::Fitting model to feature number 1025, ASV1399
#> 2026-06-17 13:45:14.33 INFO::Fitting model to feature number 1026, ASV1400
#> 2026-06-17 13:45:14.34 INFO::Fitting model to feature number 1027, ASV1401
#> 2026-06-17 13:45:14.34 INFO::Fitting model to feature number 1028, ASV1403
#> 2026-06-17 13:45:14.35 INFO::Fitting model to feature number 1029, ASV1404
#> 2026-06-17 13:45:14.35 INFO::Fitting model to feature number 1030, ASV1406
#> 2026-06-17 13:45:14.36 INFO::Fitting model to feature number 1031, ASV1408
#> 2026-06-17 13:45:14.37 INFO::Fitting model to feature number 1032, ASV1409
#> 2026-06-17 13:45:14.37 INFO::Fitting model to feature number 1033, ASV1410
#> 2026-06-17 13:45:14.38 INFO::Fitting model to feature number 1034, ASV1412
#> 2026-06-17 13:45:14.38 INFO::Fitting model to feature number 1035, ASV1413
#> 2026-06-17 13:45:14.39 INFO::Fitting model to feature number 1036, ASV1414
#> 2026-06-17 13:45:14.40 INFO::Fitting model to feature number 1037, ASV1415
#> 2026-06-17 13:45:14.40 INFO::Fitting model to feature number 1038, ASV1418
#> 2026-06-17 13:45:14.41 INFO::Fitting model to feature number 1039, ASV1419
#> 2026-06-17 13:45:14.41 INFO::Fitting model to feature number 1040, ASV1420
#> 2026-06-17 13:45:14.42 INFO::Fitting model to feature number 1041, ASV1421
#> 2026-06-17 13:45:14.42 INFO::Fitting model to feature number 1042, ASV1422
#> 2026-06-17 13:45:14.43 INFO::Fitting model to feature number 1043, ASV1423
#> 2026-06-17 13:45:14.44 INFO::Fitting model to feature number 1044, ASV1424
#> 2026-06-17 13:45:14.44 INFO::Fitting model to feature number 1045, ASV1425
#> 2026-06-17 13:45:14.45 INFO::Fitting model to feature number 1046, ASV1427
#> 2026-06-17 13:45:14.45 INFO::Fitting model to feature number 1047, ASV1428
#> 2026-06-17 13:45:14.46 INFO::Fitting model to feature number 1048, ASV1430
#> 2026-06-17 13:45:14.47 INFO::Fitting model to feature number 1049, ASV1433
#> 2026-06-17 13:45:14.47 INFO::Fitting model to feature number 1050, ASV1434
#> 2026-06-17 13:45:14.48 INFO::Fitting model to feature number 1051, ASV1435
#> 2026-06-17 13:45:14.48 INFO::Fitting model to feature number 1052, ASV1436
#> 2026-06-17 13:45:14.49 INFO::Fitting model to feature number 1053, ASV1437
#> 2026-06-17 13:45:14.50 INFO::Fitting model to feature number 1054, ASV1438
#> 2026-06-17 13:45:14.50 INFO::Fitting model to feature number 1055, ASV1442
#> 2026-06-17 13:45:14.51 INFO::Fitting model to feature number 1056, ASV1443
#> 2026-06-17 13:45:14.51 INFO::Fitting model to feature number 1057, ASV1449
#> 2026-06-17 13:45:14.52 INFO::Fitting model to feature number 1058, ASV1450
#> 2026-06-17 13:45:14.52 INFO::Fitting model to feature number 1059, ASV1452
#> 2026-06-17 13:45:14.53 INFO::Fitting model to feature number 1060, ASV1454
#> 2026-06-17 13:45:14.54 INFO::Fitting model to feature number 1061, ASV1455
#> 2026-06-17 13:45:14.54 INFO::Fitting model to feature number 1062, ASV1458
#> 2026-06-17 13:45:14.55 INFO::Fitting model to feature number 1063, ASV1459
#> 2026-06-17 13:45:14.55 INFO::Fitting model to feature number 1064, ASV1460
#> 2026-06-17 13:45:14.56 INFO::Fitting model to feature number 1065, ASV1461
#> 2026-06-17 13:45:14.57 INFO::Fitting model to feature number 1066, ASV1462
#> 2026-06-17 13:45:14.58 INFO::Fitting model to feature number 1067, ASV1463
#> 2026-06-17 13:45:14.58 INFO::Fitting model to feature number 1068, ASV1466
#> 2026-06-17 13:45:14.59 INFO::Fitting model to feature number 1069, ASV1467
#> 2026-06-17 13:45:14.60 INFO::Fitting model to feature number 1070, ASV1468
#> 2026-06-17 13:45:14.60 INFO::Fitting model to feature number 1071, ASV1469
#> 2026-06-17 13:45:14.61 INFO::Fitting model to feature number 1072, ASV1472
#> 2026-06-17 13:45:14.62 INFO::Fitting model to feature number 1073, ASV1477
#> 2026-06-17 13:45:14.63 INFO::Fitting model to feature number 1074, ASV1478
#> 2026-06-17 13:45:14.63 INFO::Fitting model to feature number 1075, ASV1479
#> 2026-06-17 13:45:14.64 INFO::Fitting model to feature number 1076, ASV1483
#> 2026-06-17 13:45:14.65 INFO::Fitting model to feature number 1077, ASV1484
#> 2026-06-17 13:45:14.65 INFO::Fitting model to feature number 1078, ASV1486
#> 2026-06-17 13:45:14.66 INFO::Fitting model to feature number 1079, ASV1487
#> 2026-06-17 13:45:14.66 INFO::Fitting model to feature number 1080, ASV1488
#> 2026-06-17 13:45:14.67 INFO::Fitting model to feature number 1081, ASV1490
#> 2026-06-17 13:45:14.68 INFO::Fitting model to feature number 1082, ASV1492
#> 2026-06-17 13:45:14.68 INFO::Fitting model to feature number 1083, ASV1493
#> 2026-06-17 13:45:14.69 INFO::Fitting model to feature number 1084, ASV1494
#> 2026-06-17 13:45:14.69 INFO::Fitting model to feature number 1085, ASV1495
#> 2026-06-17 13:45:14.70 INFO::Fitting model to feature number 1086, ASV1496
#> 2026-06-17 13:45:14.71 INFO::Fitting model to feature number 1087, ASV1497
#> 2026-06-17 13:45:14.71 INFO::Fitting model to feature number 1088, ASV1498
#> 2026-06-17 13:45:14.72 INFO::Fitting model to feature number 1089, ASV1500
#> 2026-06-17 13:45:14.72 INFO::Fitting model to feature number 1090, ASV1501
#> 2026-06-17 13:45:14.73 INFO::Fitting model to feature number 1091, ASV1502
#> 2026-06-17 13:45:14.74 INFO::Fitting model to feature number 1092, ASV1503
#> 2026-06-17 13:45:14.74 INFO::Fitting model to feature number 1093, ASV1504
#> 2026-06-17 13:45:14.75 INFO::Fitting model to feature number 1094, ASV1507
#> 2026-06-17 13:45:14.75 INFO::Fitting model to feature number 1095, ASV1509
#> 2026-06-17 13:45:14.76 INFO::Fitting model to feature number 1096, ASV1510
#> 2026-06-17 13:45:14.77 INFO::Fitting model to feature number 1097, ASV1513
#> 2026-06-17 13:45:14.77 INFO::Fitting model to feature number 1098, ASV1515
#> 2026-06-17 13:45:14.78 INFO::Fitting model to feature number 1099, ASV1516
#> 2026-06-17 13:45:14.78 INFO::Fitting model to feature number 1100, ASV1517
#> 2026-06-17 13:45:14.79 INFO::Fitting model to feature number 1101, ASV1518
#> 2026-06-17 13:45:14.80 INFO::Fitting model to feature number 1102, ASV1521
#> 2026-06-17 13:45:14.80 INFO::Fitting model to feature number 1103, ASV1523
#> 2026-06-17 13:45:14.81 INFO::Fitting model to feature number 1104, ASV1525
#> 2026-06-17 13:45:14.81 INFO::Fitting model to feature number 1105, ASV1526
#> 2026-06-17 13:45:14.82 INFO::Fitting model to feature number 1106, ASV1527
#> 2026-06-17 13:45:14.83 INFO::Fitting model to feature number 1107, ASV1528
#> 2026-06-17 13:45:14.83 INFO::Fitting model to feature number 1108, ASV1531
#> 2026-06-17 13:45:14.84 INFO::Fitting model to feature number 1109, ASV1532
#> 2026-06-17 13:45:14.84 INFO::Fitting model to feature number 1110, ASV1533
#> 2026-06-17 13:45:14.85 INFO::Fitting model to feature number 1111, ASV1534
#> 2026-06-17 13:45:14.86 INFO::Fitting model to feature number 1112, ASV1536
#> 2026-06-17 13:45:14.86 INFO::Fitting model to feature number 1113, ASV1537
#> 2026-06-17 13:45:14.87 INFO::Fitting model to feature number 1114, ASV1542
#> 2026-06-17 13:45:14.87 INFO::Fitting model to feature number 1115, ASV1543
#> 2026-06-17 13:45:14.88 INFO::Fitting model to feature number 1116, ASV1544
#> 2026-06-17 13:45:14.88 INFO::Fitting model to feature number 1117, ASV1545
#> 2026-06-17 13:45:14.89 INFO::Fitting model to feature number 1118, ASV1548
#> 2026-06-17 13:45:14.90 INFO::Fitting model to feature number 1119, ASV1549
#> 2026-06-17 13:45:14.90 INFO::Fitting model to feature number 1120, ASV1550
#> 2026-06-17 13:45:14.91 INFO::Fitting model to feature number 1121, ASV1551
#> 2026-06-17 13:45:14.92 INFO::Fitting model to feature number 1122, ASV1552
#> 2026-06-17 13:45:14.92 INFO::Fitting model to feature number 1123, ASV1553
#> 2026-06-17 13:45:14.93 INFO::Fitting model to feature number 1124, ASV1554
#> 2026-06-17 13:45:14.93 INFO::Fitting model to feature number 1125, ASV1557
#> 2026-06-17 13:45:14.94 INFO::Fitting model to feature number 1126, ASV1558
#> 2026-06-17 13:45:14.95 INFO::Fitting model to feature number 1127, ASV1560
#> 2026-06-17 13:45:14.95 INFO::Fitting model to feature number 1128, ASV1562
#> 2026-06-17 13:45:14.96 INFO::Fitting model to feature number 1129, ASV1564
#> 2026-06-17 13:45:14.96 INFO::Fitting model to feature number 1130, ASV1565
#> 2026-06-17 13:45:14.97 INFO::Fitting model to feature number 1131, ASV1568
#> 2026-06-17 13:45:14.98 INFO::Fitting model to feature number 1132, ASV1569
#> 2026-06-17 13:45:14.98 INFO::Fitting model to feature number 1133, ASV1570
#> 2026-06-17 13:45:14.99 INFO::Fitting model to feature number 1134, ASV1571
#> 2026-06-17 13:45:14.99 INFO::Fitting model to feature number 1135, ASV1573
#> 2026-06-17 13:45:15.00 INFO::Fitting model to feature number 1136, ASV1574
#> 2026-06-17 13:45:15.01 INFO::Fitting model to feature number 1137, ASV1575
#> 2026-06-17 13:45:15.01 INFO::Fitting model to feature number 1138, ASV1576
#> 2026-06-17 13:45:15.02 INFO::Fitting model to feature number 1139, ASV1577
#> 2026-06-17 13:45:15.02 INFO::Fitting model to feature number 1140, ASV1578
#> 2026-06-17 13:45:15.03 INFO::Fitting model to feature number 1141, ASV1580
#> 2026-06-17 13:45:15.03 INFO::Fitting model to feature number 1142, ASV1581
#> 2026-06-17 13:45:15.04 INFO::Fitting model to feature number 1143, ASV1584
#> 2026-06-17 13:45:15.05 INFO::Fitting model to feature number 1144, ASV1585
#> 2026-06-17 13:45:15.05 INFO::Fitting model to feature number 1145, ASV1586
#> 2026-06-17 13:45:15.06 INFO::Fitting model to feature number 1146, ASV1588
#> 2026-06-17 13:45:15.06 INFO::Fitting model to feature number 1147, ASV1589
#> 2026-06-17 13:45:15.07 INFO::Fitting model to feature number 1148, ASV1592
#> 2026-06-17 13:45:15.07 INFO::Fitting model to feature number 1149, ASV1593
#> 2026-06-17 13:45:15.08 INFO::Fitting model to feature number 1150, ASV1594
#> 2026-06-17 13:45:15.09 INFO::Fitting model to feature number 1151, ASV1601
#> 2026-06-17 13:45:15.09 INFO::Fitting model to feature number 1152, ASV1602
#> 2026-06-17 13:45:15.10 INFO::Fitting model to feature number 1153, ASV1603
#> 2026-06-17 13:45:15.10 INFO::Fitting model to feature number 1154, ASV1604
#> 2026-06-17 13:45:15.11 INFO::Fitting model to feature number 1155, ASV1605
#> 2026-06-17 13:45:15.11 INFO::Fitting model to feature number 1156, ASV1606
#> 2026-06-17 13:45:15.12 INFO::Fitting model to feature number 1157, ASV1607
#> 2026-06-17 13:45:15.13 INFO::Fitting model to feature number 1158, ASV1609
#> 2026-06-17 13:45:15.13 INFO::Fitting model to feature number 1159, ASV1612
#> 2026-06-17 13:45:15.14 INFO::Fitting model to feature number 1160, ASV1613
#> 2026-06-17 13:45:15.14 INFO::Fitting model to feature number 1161, ASV1614
#> 2026-06-17 13:45:15.15 INFO::Fitting model to feature number 1162, ASV1615
#> 2026-06-17 13:45:15.15 INFO::Fitting model to feature number 1163, ASV1616
#> 2026-06-17 13:45:15.16 INFO::Fitting model to feature number 1164, ASV1620
#> 2026-06-17 13:45:15.17 INFO::Fitting model to feature number 1165, ASV1621
#> 2026-06-17 13:45:15.17 INFO::Fitting model to feature number 1166, ASV1622
#> 2026-06-17 13:45:15.18 INFO::Fitting model to feature number 1167, ASV1624
#> 2026-06-17 13:45:15.18 INFO::Fitting model to feature number 1168, ASV1625
#> 2026-06-17 13:45:15.19 INFO::Fitting model to feature number 1169, ASV1628
#> 2026-06-17 13:45:15.20 INFO::Fitting model to feature number 1170, ASV1629
#> 2026-06-17 13:45:15.20 INFO::Fitting model to feature number 1171, ASV1630
#> 2026-06-17 13:45:15.21 INFO::Fitting model to feature number 1172, ASV1631
#> 2026-06-17 13:45:15.21 INFO::Fitting model to feature number 1173, ASV1632
#> 2026-06-17 13:45:15.22 INFO::Fitting model to feature number 1174, ASV1633
#> 2026-06-17 13:45:15.22 INFO::Fitting model to feature number 1175, ASV1634
#> 2026-06-17 13:45:15.23 INFO::Fitting model to feature number 1176, ASV1635
#> 2026-06-17 13:45:15.24 INFO::Fitting model to feature number 1177, ASV1636
#> 2026-06-17 13:45:15.24 INFO::Fitting model to feature number 1178, ASV1639
#> 2026-06-17 13:45:15.25 INFO::Fitting model to feature number 1179, ASV1641
#> 2026-06-17 13:45:15.25 INFO::Fitting model to feature number 1180, ASV1644
#> 2026-06-17 13:45:15.26 INFO::Fitting model to feature number 1181, ASV1645
#> 2026-06-17 13:45:15.27 INFO::Fitting model to feature number 1182, ASV1649
#> 2026-06-17 13:45:15.27 INFO::Fitting model to feature number 1183, ASV1650
#> 2026-06-17 13:45:15.28 INFO::Fitting model to feature number 1184, ASV1651
#> 2026-06-17 13:45:15.28 INFO::Fitting model to feature number 1185, ASV1653
#> 2026-06-17 13:45:15.29 INFO::Fitting model to feature number 1186, ASV1654
#> 2026-06-17 13:45:15.30 INFO::Fitting model to feature number 1187, ASV1656
#> 2026-06-17 13:45:15.30 INFO::Fitting model to feature number 1188, ASV1657
#> 2026-06-17 13:45:15.31 INFO::Fitting model to feature number 1189, ASV1658
#> 2026-06-17 13:45:15.31 INFO::Fitting model to feature number 1190, ASV1661
#> 2026-06-17 13:45:15.32 INFO::Fitting model to feature number 1191, ASV1663
#> 2026-06-17 13:45:15.32 INFO::Fitting model to feature number 1192, ASV1664
#> 2026-06-17 13:45:15.33 INFO::Fitting model to feature number 1193, ASV1666
#> 2026-06-17 13:45:15.34 INFO::Fitting model to feature number 1194, ASV1667
#> 2026-06-17 13:45:15.34 INFO::Fitting model to feature number 1195, ASV1670
#> 2026-06-17 13:45:15.35 INFO::Fitting model to feature number 1196, ASV1671
#> 2026-06-17 13:45:15.35 INFO::Fitting model to feature number 1197, ASV1673
#> 2026-06-17 13:45:15.36 INFO::Fitting model to feature number 1198, ASV1674
#> 2026-06-17 13:45:15.37 INFO::Fitting model to feature number 1199, ASV1677
#> 2026-06-17 13:45:15.37 INFO::Fitting model to feature number 1200, ASV1681
#> 2026-06-17 13:45:15.38 INFO::Fitting model to feature number 1201, ASV1683
#> 2026-06-17 13:45:15.38 INFO::Fitting model to feature number 1202, ASV1684
#> 2026-06-17 13:45:15.39 INFO::Fitting model to feature number 1203, ASV1687
#> 2026-06-17 13:45:15.40 INFO::Fitting model to feature number 1204, ASV1689
#> 2026-06-17 13:45:15.40 INFO::Fitting model to feature number 1205, ASV1690
#> 2026-06-17 13:45:15.41 INFO::Fitting model to feature number 1206, ASV1691
#> 2026-06-17 13:45:15.41 INFO::Fitting model to feature number 1207, ASV1694
#> 2026-06-17 13:45:15.42 INFO::Fitting model to feature number 1208, ASV1697
#> 2026-06-17 13:45:15.42 INFO::Fitting model to feature number 1209, ASV1704
#> 2026-06-17 13:45:15.43 INFO::Fitting model to feature number 1210, ASV1706
#> 2026-06-17 13:45:15.44 INFO::Fitting model to feature number 1211, ASV1707
#> 2026-06-17 13:45:15.44 INFO::Fitting model to feature number 1212, ASV1712
#> 2026-06-17 13:45:15.57 INFO::Counting total values for each feature
#> 2026-06-17 13:45:15.69 INFO::Re-running abundances for warn_prevalence
#> 2026-06-17 13:45:15.69 INFO::Running selected normalization method: TSS
#> 2026-06-17 13:45:15.79 INFO::Running selected transform method: LOG
#> 2026-06-17 13:45:15.82 INFO::Fitting model to feature number 1, ASV2
#> 2026-06-17 13:45:15.82 INFO::Fitting model to feature number 2, ASV6
#> 2026-06-17 13:45:15.83 INFO::Fitting model to feature number 3, ASV7
#> 2026-06-17 13:45:15.83 INFO::Fitting model to feature number 4, ASV8
#> 2026-06-17 13:45:15.84 INFO::Fitting model to feature number 5, ASV10
#> 2026-06-17 13:45:15.84 INFO::Fitting model to feature number 6, ASV12
#> 2026-06-17 13:45:15.84 INFO::Fitting model to feature number 7, ASV13
#> 2026-06-17 13:45:15.85 INFO::Fitting model to feature number 8, ASV18
#> 2026-06-17 13:45:15.85 INFO::Fitting model to feature number 9, ASV19
#> 2026-06-17 13:45:15.86 INFO::Fitting model to feature number 10, ASV22
#> 2026-06-17 13:45:15.86 INFO::Fitting model to feature number 11, ASV23
#> 2026-06-17 13:45:15.86 WARNING::Fitting problem for feature 11 returning NA
#> 2026-06-17 13:45:15.86 INFO::Fitting model to feature number 12, ASV24
#> 2026-06-17 13:45:15.87 INFO::Fitting model to feature number 13, ASV25
#> 2026-06-17 13:45:15.87 INFO::Fitting model to feature number 14, ASV26
#> 2026-06-17 13:45:15.87 INFO::Fitting model to feature number 15, ASV27
#> 2026-06-17 13:45:15.88 INFO::Fitting model to feature number 16, ASV28
#> 2026-06-17 13:45:15.88 INFO::Fitting model to feature number 17, ASV29
#> 2026-06-17 13:45:15.89 INFO::Fitting model to feature number 18, ASV31
#> 2026-06-17 13:45:15.89 WARNING::Fitting problem for feature 18 returning NA
#> 2026-06-17 13:45:15.89 INFO::Fitting model to feature number 19, ASV32
#> 2026-06-17 13:45:15.89 INFO::Fitting model to feature number 20, ASV33
#> 2026-06-17 13:45:15.90 INFO::Fitting model to feature number 21, ASV34
#> 2026-06-17 13:45:15.90 INFO::Fitting model to feature number 22, ASV35
#> 2026-06-17 13:45:15.91 INFO::Fitting model to feature number 23, ASV38
#> 2026-06-17 13:45:15.91 INFO::Fitting model to feature number 24, ASV41
#> 2026-06-17 13:45:15.91 INFO::Fitting model to feature number 25, ASV42
#> 2026-06-17 13:45:15.91 INFO::Fitting model to feature number 26, ASV43
#> 2026-06-17 13:45:15.92 INFO::Fitting model to feature number 27, ASV45
#> 2026-06-17 13:45:15.92 INFO::Fitting model to feature number 28, ASV46
#> 2026-06-17 13:45:15.92 INFO::Fitting model to feature number 29, ASV47
#> 2026-06-17 13:45:15.93 INFO::Fitting model to feature number 30, ASV48
#> 2026-06-17 13:45:15.93 WARNING::Fitting problem for feature 30 returning NA
#> 2026-06-17 13:45:15.93 INFO::Fitting model to feature number 31, ASV49
#> 2026-06-17 13:45:15.93 INFO::Fitting model to feature number 32, ASV50
#> 2026-06-17 13:45:15.94 WARNING::Fitting problem for feature 32 returning NA
#> 2026-06-17 13:45:15.94 INFO::Fitting model to feature number 33, ASV51
#> 2026-06-17 13:45:15.94 INFO::Fitting model to feature number 34, ASV52
#> 2026-06-17 13:45:15.95 INFO::Fitting model to feature number 35, ASV53
#> 2026-06-17 13:45:15.95 INFO::Fitting model to feature number 36, ASV55
#> 2026-06-17 13:45:15.95 INFO::Fitting model to feature number 37, ASV56
#> 2026-06-17 13:45:15.96 INFO::Fitting model to feature number 38, ASV57
#> 2026-06-17 13:45:15.96 INFO::Fitting model to feature number 39, ASV58
#> 2026-06-17 13:45:15.97 INFO::Fitting model to feature number 40, ASV59
#> 2026-06-17 13:45:15.97 INFO::Fitting model to feature number 41, ASV60
#> 2026-06-17 13:45:15.97 INFO::Fitting model to feature number 42, ASV61
#> 2026-06-17 13:45:15.97 INFO::Fitting model to feature number 43, ASV62
#> 2026-06-17 13:45:15.98 INFO::Fitting model to feature number 44, ASV63
#> 2026-06-17 13:45:15.98 INFO::Fitting model to feature number 45, ASV64
#> 2026-06-17 13:45:15.98 INFO::Fitting model to feature number 46, ASV65
#> 2026-06-17 13:45:15.99 INFO::Fitting model to feature number 47, ASV67
#> 2026-06-17 13:45:15.99 INFO::Fitting model to feature number 48, ASV68
#> 2026-06-17 13:45:15.99 INFO::Fitting model to feature number 49, ASV69
#> 2026-06-17 13:45:16.00 INFO::Fitting model to feature number 50, ASV70
#> 2026-06-17 13:45:16.00 INFO::Fitting model to feature number 51, ASV71
#> 2026-06-17 13:45:16.00 INFO::Fitting model to feature number 52, ASV72
#> 2026-06-17 13:45:16.01 INFO::Fitting model to feature number 53, ASV75
#> 2026-06-17 13:45:16.01 INFO::Fitting model to feature number 54, ASV77
#> 2026-06-17 13:45:16.01 WARNING::Fitting problem for feature 54 returning NA
#> 2026-06-17 13:45:16.02 INFO::Fitting model to feature number 55, ASV78
#> 2026-06-17 13:45:16.02 INFO::Fitting model to feature number 56, ASV79
#> 2026-06-17 13:45:16.02 INFO::Fitting model to feature number 57, ASV80
#> 2026-06-17 13:45:16.03 INFO::Fitting model to feature number 58, ASV81
#> 2026-06-17 13:45:16.03 INFO::Fitting model to feature number 59, ASV82
#> 2026-06-17 13:45:16.03 INFO::Fitting model to feature number 60, ASV83
#> 2026-06-17 13:45:16.04 INFO::Fitting model to feature number 61, ASV84
#> 2026-06-17 13:45:16.04 INFO::Fitting model to feature number 62, ASV85
#> 2026-06-17 13:45:16.04 INFO::Fitting model to feature number 63, ASV87
#> 2026-06-17 13:45:16.04 INFO::Fitting model to feature number 64, ASV89
#> 2026-06-17 13:45:16.05 INFO::Fitting model to feature number 65, ASV90
#> 2026-06-17 13:45:16.05 INFO::Fitting model to feature number 66, ASV91
#> 2026-06-17 13:45:16.05 INFO::Fitting model to feature number 67, ASV92
#> 2026-06-17 13:45:16.06 INFO::Fitting model to feature number 68, ASV93
#> 2026-06-17 13:45:16.06 WARNING::Fitting problem for feature 68 returning NA
#> 2026-06-17 13:45:16.06 INFO::Fitting model to feature number 69, ASV94
#> 2026-06-17 13:45:16.06 INFO::Fitting model to feature number 70, ASV95
#> 2026-06-17 13:45:16.07 INFO::Fitting model to feature number 71, ASV96
#> 2026-06-17 13:45:16.07 INFO::Fitting model to feature number 72, ASV98
#> 2026-06-17 13:45:16.08 INFO::Fitting model to feature number 73, ASV99
#> 2026-06-17 13:45:16.08 INFO::Fitting model to feature number 74, ASV100
#> 2026-06-17 13:45:16.08 INFO::Fitting model to feature number 75, ASV101
#> 2026-06-17 13:45:16.09 INFO::Fitting model to feature number 76, ASV103
#> 2026-06-17 13:45:16.09 INFO::Fitting model to feature number 77, ASV104
#> 2026-06-17 13:45:16.09 INFO::Fitting model to feature number 78, ASV105
#> 2026-06-17 13:45:16.09 INFO::Fitting model to feature number 79, ASV106
#> 2026-06-17 13:45:16.10 INFO::Fitting model to feature number 80, ASV107
#> 2026-06-17 13:45:16.10 INFO::Fitting model to feature number 81, ASV109
#> 2026-06-17 13:45:16.11 INFO::Fitting model to feature number 82, ASV111
#> 2026-06-17 13:45:16.11 INFO::Fitting model to feature number 83, ASV112
#> 2026-06-17 13:45:16.11 INFO::Fitting model to feature number 84, ASV113
#> 2026-06-17 13:45:16.12 INFO::Fitting model to feature number 85, ASV115
#> 2026-06-17 13:45:16.12 INFO::Fitting model to feature number 86, ASV116
#> 2026-06-17 13:45:16.13 INFO::Fitting model to feature number 87, ASV118
#> 2026-06-17 13:45:16.13 INFO::Fitting model to feature number 88, ASV119
#> 2026-06-17 13:45:16.13 INFO::Fitting model to feature number 89, ASV120
#> 2026-06-17 13:45:16.13 INFO::Fitting model to feature number 90, ASV121
#> 2026-06-17 13:45:16.14 INFO::Fitting model to feature number 91, ASV124
#> 2026-06-17 13:45:16.14 INFO::Fitting model to feature number 92, ASV126
#> 2026-06-17 13:45:16.15 INFO::Fitting model to feature number 93, ASV127
#> 2026-06-17 13:45:16.15 WARNING::Fitting problem for feature 93 returning NA
#> 2026-06-17 13:45:16.15 INFO::Fitting model to feature number 94, ASV128
#> 2026-06-17 13:45:16.15 INFO::Fitting model to feature number 95, ASV129
#> 2026-06-17 13:45:16.16 INFO::Fitting model to feature number 96, ASV130
#> 2026-06-17 13:45:16.16 INFO::Fitting model to feature number 97, ASV131
#> 2026-06-17 13:45:16.17 INFO::Fitting model to feature number 98, ASV132
#> 2026-06-17 13:45:16.17 INFO::Fitting model to feature number 99, ASV133
#> 2026-06-17 13:45:16.17 INFO::Fitting model to feature number 100, ASV135
#> 2026-06-17 13:45:16.17 INFO::Fitting model to feature number 101, ASV136
#> 2026-06-17 13:45:16.18 INFO::Fitting model to feature number 102, ASV137
#> 2026-06-17 13:45:16.18 INFO::Fitting model to feature number 103, ASV138
#> 2026-06-17 13:45:16.19 INFO::Fitting model to feature number 104, ASV139
#> 2026-06-17 13:45:16.19 INFO::Fitting model to feature number 105, ASV142
#> 2026-06-17 13:45:16.19 INFO::Fitting model to feature number 106, ASV143
#> 2026-06-17 13:45:16.19 INFO::Fitting model to feature number 107, ASV144
#> 2026-06-17 13:45:16.20 INFO::Fitting model to feature number 108, ASV145
#> 2026-06-17 13:45:16.20 INFO::Fitting model to feature number 109, ASV146
#> 2026-06-17 13:45:16.21 INFO::Fitting model to feature number 110, ASV147
#> 2026-06-17 13:45:16.21 INFO::Fitting model to feature number 111, ASV148
#> 2026-06-17 13:45:16.21 INFO::Fitting model to feature number 112, ASV149
#> 2026-06-17 13:45:16.22 INFO::Fitting model to feature number 113, ASV150
#> 2026-06-17 13:45:16.22 INFO::Fitting model to feature number 114, ASV151
#> 2026-06-17 13:45:16.22 INFO::Fitting model to feature number 115, ASV153
#> 2026-06-17 13:45:16.22 INFO::Fitting model to feature number 116, ASV154
#> 2026-06-17 13:45:16.23 INFO::Fitting model to feature number 117, ASV157
#> 2026-06-17 13:45:16.23 INFO::Fitting model to feature number 118, ASV158
#> 2026-06-17 13:45:16.23 INFO::Fitting model to feature number 119, ASV159
#> 2026-06-17 13:45:16.24 INFO::Fitting model to feature number 120, ASV162
#> 2026-06-17 13:45:16.24 INFO::Fitting model to feature number 121, ASV163
#> 2026-06-17 13:45:16.24 INFO::Fitting model to feature number 122, ASV164
#> 2026-06-17 13:45:16.25 INFO::Fitting model to feature number 123, ASV166
#> 2026-06-17 13:45:16.25 INFO::Fitting model to feature number 124, ASV167
#> 2026-06-17 13:45:16.26 INFO::Fitting model to feature number 125, ASV170
#> 2026-06-17 13:45:16.26 INFO::Fitting model to feature number 126, ASV171
#> 2026-06-17 13:45:16.26 INFO::Fitting model to feature number 127, ASV172
#> 2026-06-17 13:45:16.27 INFO::Fitting model to feature number 128, ASV173
#> 2026-06-17 13:45:16.27 INFO::Fitting model to feature number 129, ASV175
#> 2026-06-17 13:45:16.27 INFO::Fitting model to feature number 130, ASV176
#> 2026-06-17 13:45:16.28 INFO::Fitting model to feature number 131, ASV178
#> 2026-06-17 13:45:16.28 INFO::Fitting model to feature number 132, ASV179
#> 2026-06-17 13:45:16.28 INFO::Fitting model to feature number 133, ASV181
#> 2026-06-17 13:45:16.29 INFO::Fitting model to feature number 134, ASV182
#> 2026-06-17 13:45:16.29 INFO::Fitting model to feature number 135, ASV183
#> 2026-06-17 13:45:16.29 INFO::Fitting model to feature number 136, ASV186
#> 2026-06-17 13:45:16.30 INFO::Fitting model to feature number 137, ASV187
#> 2026-06-17 13:45:16.30 INFO::Fitting model to feature number 138, ASV189
#> 2026-06-17 13:45:16.31 INFO::Fitting model to feature number 139, ASV190
#> 2026-06-17 13:45:16.31 INFO::Fitting model to feature number 140, ASV192
#> 2026-06-17 13:45:16.31 INFO::Fitting model to feature number 141, ASV193
#> 2026-06-17 13:45:16.31 INFO::Fitting model to feature number 142, ASV194
#> 2026-06-17 13:45:16.32 INFO::Fitting model to feature number 143, ASV195
#> 2026-06-17 13:45:16.32 INFO::Fitting model to feature number 144, ASV196
#> 2026-06-17 13:45:16.32 INFO::Fitting model to feature number 145, ASV197
#> 2026-06-17 13:45:16.33 INFO::Fitting model to feature number 146, ASV198
#> 2026-06-17 13:45:16.33 WARNING::Fitting problem for feature 146 returning NA
#> 2026-06-17 13:45:16.33 INFO::Fitting model to feature number 147, ASV199
#> 2026-06-17 13:45:16.33 INFO::Fitting model to feature number 148, ASV201
#> 2026-06-17 13:45:16.34 INFO::Fitting model to feature number 149, ASV202
#> 2026-06-17 13:45:16.34 INFO::Fitting model to feature number 150, ASV203
#> 2026-06-17 13:45:16.34 INFO::Fitting model to feature number 151, ASV205
#> 2026-06-17 13:45:16.34 INFO::Fitting model to feature number 152, ASV208
#> 2026-06-17 13:45:16.35 INFO::Fitting model to feature number 153, ASV209
#> 2026-06-17 13:45:16.35 INFO::Fitting model to feature number 154, ASV210
#> 2026-06-17 13:45:16.35 INFO::Fitting model to feature number 155, ASV211
#> 2026-06-17 13:45:16.36 INFO::Fitting model to feature number 156, ASV212
#> 2026-06-17 13:45:16.36 INFO::Fitting model to feature number 157, ASV214
#> 2026-06-17 13:45:16.36 INFO::Fitting model to feature number 158, ASV215
#> 2026-06-17 13:45:16.36 INFO::Fitting model to feature number 159, ASV216
#> 2026-06-17 13:45:16.37 INFO::Fitting model to feature number 160, ASV217
#> 2026-06-17 13:45:16.37 INFO::Fitting model to feature number 161, ASV219
#> 2026-06-17 13:45:16.37 INFO::Fitting model to feature number 162, ASV221
#> 2026-06-17 13:45:16.38 INFO::Fitting model to feature number 163, ASV222
#> 2026-06-17 13:45:16.38 WARNING::Fitting problem for feature 163 returning NA
#> 2026-06-17 13:45:16.38 INFO::Fitting model to feature number 164, ASV223
#> 2026-06-17 13:45:16.39 INFO::Fitting model to feature number 165, ASV224
#> 2026-06-17 13:45:16.39 INFO::Fitting model to feature number 166, ASV226
#> 2026-06-17 13:45:16.39 INFO::Fitting model to feature number 167, ASV227
#> 2026-06-17 13:45:16.40 INFO::Fitting model to feature number 168, ASV228
#> 2026-06-17 13:45:16.40 INFO::Fitting model to feature number 169, ASV229
#> 2026-06-17 13:45:16.40 INFO::Fitting model to feature number 170, ASV231
#> 2026-06-17 13:45:16.40 WARNING::Fitting problem for feature 170 returning NA
#> 2026-06-17 13:45:16.40 INFO::Fitting model to feature number 171, ASV233
#> 2026-06-17 13:45:16.41 INFO::Fitting model to feature number 172, ASV234
#> 2026-06-17 13:45:16.41 INFO::Fitting model to feature number 173, ASV235
#> 2026-06-17 13:45:16.42 INFO::Fitting model to feature number 174, ASV237
#> 2026-06-17 13:45:16.42 INFO::Fitting model to feature number 175, ASV238
#> 2026-06-17 13:45:16.42 INFO::Fitting model to feature number 176, ASV239
#> 2026-06-17 13:45:16.43 INFO::Fitting model to feature number 177, ASV240
#> 2026-06-17 13:45:16.43 INFO::Fitting model to feature number 178, ASV243
#> 2026-06-17 13:45:16.43 INFO::Fitting model to feature number 179, ASV244
#> 2026-06-17 13:45:16.44 INFO::Fitting model to feature number 180, ASV245
#> 2026-06-17 13:45:16.44 INFO::Fitting model to feature number 181, ASV246
#> 2026-06-17 13:45:16.44 INFO::Fitting model to feature number 182, ASV247
#> 2026-06-17 13:45:16.44 INFO::Fitting model to feature number 183, ASV248
#> 2026-06-17 13:45:16.45 INFO::Fitting model to feature number 184, ASV249
#> 2026-06-17 13:45:16.45 WARNING::Fitting problem for feature 184 returning NA
#> 2026-06-17 13:45:16.45 INFO::Fitting model to feature number 185, ASV251
#> 2026-06-17 13:45:16.45 INFO::Fitting model to feature number 186, ASV254
#> 2026-06-17 13:45:16.46 INFO::Fitting model to feature number 187, ASV255
#> 2026-06-17 13:45:16.46 INFO::Fitting model to feature number 188, ASV256
#> 2026-06-17 13:45:16.46 INFO::Fitting model to feature number 189, ASV257
#> 2026-06-17 13:45:16.46 INFO::Fitting model to feature number 190, ASV258
#> 2026-06-17 13:45:16.47 INFO::Fitting model to feature number 191, ASV260
#> 2026-06-17 13:45:16.47 INFO::Fitting model to feature number 192, ASV261
#> 2026-06-17 13:45:16.47 INFO::Fitting model to feature number 193, ASV262
#> 2026-06-17 13:45:16.48 INFO::Fitting model to feature number 194, ASV263
#> 2026-06-17 13:45:16.48 INFO::Fitting model to feature number 195, ASV264
#> 2026-06-17 13:45:16.49 INFO::Fitting model to feature number 196, ASV265
#> 2026-06-17 13:45:16.49 INFO::Fitting model to feature number 197, ASV266
#> 2026-06-17 13:45:16.49 INFO::Fitting model to feature number 198, ASV267
#> 2026-06-17 13:45:16.50 INFO::Fitting model to feature number 199, ASV270
#> 2026-06-17 13:45:16.50 INFO::Fitting model to feature number 200, ASV271
#> 2026-06-17 13:45:16.50 INFO::Fitting model to feature number 201, ASV272
#> 2026-06-17 13:45:16.51 INFO::Fitting model to feature number 202, ASV273
#> 2026-06-17 13:45:16.51 INFO::Fitting model to feature number 203, ASV274
#> 2026-06-17 13:45:16.51 INFO::Fitting model to feature number 204, ASV277
#> 2026-06-17 13:45:16.52 INFO::Fitting model to feature number 205, ASV278
#> 2026-06-17 13:45:16.52 INFO::Fitting model to feature number 206, ASV279
#> 2026-06-17 13:45:16.53 INFO::Fitting model to feature number 207, ASV282
#> 2026-06-17 13:45:16.53 INFO::Fitting model to feature number 208, ASV283
#> 2026-06-17 13:45:16.53 INFO::Fitting model to feature number 209, ASV285
#> 2026-06-17 13:45:16.53 INFO::Fitting model to feature number 210, ASV286
#> 2026-06-17 13:45:16.54 INFO::Fitting model to feature number 211, ASV287
#> 2026-06-17 13:45:16.54 INFO::Fitting model to feature number 212, ASV288
#> 2026-06-17 13:45:16.54 INFO::Fitting model to feature number 213, ASV292
#> 2026-06-17 13:45:16.55 INFO::Fitting model to feature number 214, ASV293
#> 2026-06-17 13:45:16.55 INFO::Fitting model to feature number 215, ASV294
#> 2026-06-17 13:45:16.55 INFO::Fitting model to feature number 216, ASV295
#> 2026-06-17 13:45:16.56 INFO::Fitting model to feature number 217, ASV297
#> 2026-06-17 13:45:16.56 INFO::Fitting model to feature number 218, ASV300
#> 2026-06-17 13:45:16.56 INFO::Fitting model to feature number 219, ASV302
#> 2026-06-17 13:45:16.57 INFO::Fitting model to feature number 220, ASV307
#> 2026-06-17 13:45:16.57 INFO::Fitting model to feature number 221, ASV309
#> 2026-06-17 13:45:16.57 INFO::Fitting model to feature number 222, ASV310
#> 2026-06-17 13:45:16.58 INFO::Fitting model to feature number 223, ASV311
#> 2026-06-17 13:45:16.58 INFO::Fitting model to feature number 224, ASV312
#> 2026-06-17 13:45:16.58 INFO::Fitting model to feature number 225, ASV313
#> 2026-06-17 13:45:16.59 INFO::Fitting model to feature number 226, ASV314
#> 2026-06-17 13:45:16.59 INFO::Fitting model to feature number 227, ASV315
#> 2026-06-17 13:45:16.59 INFO::Fitting model to feature number 228, ASV316
#> 2026-06-17 13:45:16.60 INFO::Fitting model to feature number 229, ASV318
#> 2026-06-17 13:45:16.60 INFO::Fitting model to feature number 230, ASV320
#> 2026-06-17 13:45:16.60 INFO::Fitting model to feature number 231, ASV322
#> 2026-06-17 13:45:16.60 INFO::Fitting model to feature number 232, ASV323
#> 2026-06-17 13:45:16.61 INFO::Fitting model to feature number 233, ASV327
#> 2026-06-17 13:45:16.61 INFO::Fitting model to feature number 234, ASV328
#> 2026-06-17 13:45:16.62 INFO::Fitting model to feature number 235, ASV329
#> 2026-06-17 13:45:16.62 INFO::Fitting model to feature number 236, ASV330
#> 2026-06-17 13:45:16.62 INFO::Fitting model to feature number 237, ASV333
#> 2026-06-17 13:45:16.63 INFO::Fitting model to feature number 238, ASV334
#> 2026-06-17 13:45:16.63 INFO::Fitting model to feature number 239, ASV337
#> 2026-06-17 13:45:16.63 INFO::Fitting model to feature number 240, ASV338
#> 2026-06-17 13:45:16.64 INFO::Fitting model to feature number 241, ASV339
#> 2026-06-17 13:45:16.64 INFO::Fitting model to feature number 242, ASV340
#> 2026-06-17 13:45:16.65 INFO::Fitting model to feature number 243, ASV341
#> 2026-06-17 13:45:16.65 INFO::Fitting model to feature number 244, ASV342
#> 2026-06-17 13:45:16.65 INFO::Fitting model to feature number 245, ASV344
#> 2026-06-17 13:45:16.66 INFO::Fitting model to feature number 246, ASV345
#> 2026-06-17 13:45:16.66 INFO::Fitting model to feature number 247, ASV346
#> 2026-06-17 13:45:16.67 INFO::Fitting model to feature number 248, ASV347
#> 2026-06-17 13:45:16.67 INFO::Fitting model to feature number 249, ASV348
#> 2026-06-17 13:45:16.67 INFO::Fitting model to feature number 250, ASV349
#> 2026-06-17 13:45:16.68 INFO::Fitting model to feature number 251, ASV350
#> 2026-06-17 13:45:16.68 INFO::Fitting model to feature number 252, ASV351
#> 2026-06-17 13:45:16.69 INFO::Fitting model to feature number 253, ASV353
#> 2026-06-17 13:45:16.69 INFO::Fitting model to feature number 254, ASV355
#> 2026-06-17 13:45:16.70 INFO::Fitting model to feature number 255, ASV356
#> 2026-06-17 13:45:16.70 INFO::Fitting model to feature number 256, ASV357
#> 2026-06-17 13:45:16.70 INFO::Fitting model to feature number 257, ASV358
#> 2026-06-17 13:45:16.71 INFO::Fitting model to feature number 258, ASV359
#> 2026-06-17 13:45:16.71 INFO::Fitting model to feature number 259, ASV360
#> 2026-06-17 13:45:16.71 INFO::Fitting model to feature number 260, ASV362
#> 2026-06-17 13:45:16.72 INFO::Fitting model to feature number 261, ASV364
#> 2026-06-17 13:45:16.72 INFO::Fitting model to feature number 262, ASV365
#> 2026-06-17 13:45:16.72 INFO::Fitting model to feature number 263, ASV366
#> 2026-06-17 13:45:16.73 INFO::Fitting model to feature number 264, ASV367
#> 2026-06-17 13:45:16.73 INFO::Fitting model to feature number 265, ASV368
#> 2026-06-17 13:45:16.74 INFO::Fitting model to feature number 266, ASV369
#> 2026-06-17 13:45:16.74 INFO::Fitting model to feature number 267, ASV371
#> 2026-06-17 13:45:16.75 INFO::Fitting model to feature number 268, ASV372
#> 2026-06-17 13:45:16.75 INFO::Fitting model to feature number 269, ASV373
#> 2026-06-17 13:45:16.75 WARNING::Fitting problem for feature 269 returning NA
#> 2026-06-17 13:45:16.75 INFO::Fitting model to feature number 270, ASV375
#> 2026-06-17 13:45:16.76 INFO::Fitting model to feature number 271, ASV376
#> 2026-06-17 13:45:16.76 INFO::Fitting model to feature number 272, ASV377
#> 2026-06-17 13:45:16.77 INFO::Fitting model to feature number 273, ASV378
#> 2026-06-17 13:45:16.77 INFO::Fitting model to feature number 274, ASV380
#> 2026-06-17 13:45:16.78 INFO::Fitting model to feature number 275, ASV382
#> 2026-06-17 13:45:16.78 INFO::Fitting model to feature number 276, ASV383
#> 2026-06-17 13:45:16.78 INFO::Fitting model to feature number 277, ASV385
#> 2026-06-17 13:45:16.78 INFO::Fitting model to feature number 278, ASV388
#> 2026-06-17 13:45:16.78 INFO::Fitting model to feature number 279, ASV390
#> 2026-06-17 13:45:16.79 INFO::Fitting model to feature number 280, ASV395
#> 2026-06-17 13:45:16.79 INFO::Fitting model to feature number 281, ASV399
#> 2026-06-17 13:45:16.79 INFO::Fitting model to feature number 282, ASV401
#> 2026-06-17 13:45:16.80 INFO::Fitting model to feature number 283, ASV402
#> 2026-06-17 13:45:16.80 INFO::Fitting model to feature number 284, ASV404
#> 2026-06-17 13:45:16.80 INFO::Fitting model to feature number 285, ASV405
#> 2026-06-17 13:45:16.81 INFO::Fitting model to feature number 286, ASV409
#> 2026-06-17 13:45:16.81 INFO::Fitting model to feature number 287, ASV411
#> 2026-06-17 13:45:16.82 INFO::Fitting model to feature number 288, ASV412
#> 2026-06-17 13:45:16.82 INFO::Fitting model to feature number 289, ASV413
#> 2026-06-17 13:45:16.82 INFO::Fitting model to feature number 290, ASV414
#> 2026-06-17 13:45:16.83 WARNING::Fitting problem for feature 290 returning NA
#> 2026-06-17 13:45:16.83 INFO::Fitting model to feature number 291, ASV415
#> 2026-06-17 13:45:16.83 INFO::Fitting model to feature number 292, ASV416
#> 2026-06-17 13:45:16.83 INFO::Fitting model to feature number 293, ASV417
#> 2026-06-17 13:45:16.84 INFO::Fitting model to feature number 294, ASV418
#> 2026-06-17 13:45:16.84 INFO::Fitting model to feature number 295, ASV419
#> 2026-06-17 13:45:16.85 INFO::Fitting model to feature number 296, ASV420
#> 2026-06-17 13:45:16.85 INFO::Fitting model to feature number 297, ASV422
#> 2026-06-17 13:45:16.85 INFO::Fitting model to feature number 298, ASV423
#> 2026-06-17 13:45:16.85 INFO::Fitting model to feature number 299, ASV424
#> 2026-06-17 13:45:16.86 INFO::Fitting model to feature number 300, ASV426
#> 2026-06-17 13:45:16.86 INFO::Fitting model to feature number 301, ASV427
#> 2026-06-17 13:45:16.86 INFO::Fitting model to feature number 302, ASV428
#> 2026-06-17 13:45:16.87 INFO::Fitting model to feature number 303, ASV429
#> 2026-06-17 13:45:16.87 INFO::Fitting model to feature number 304, ASV431
#> 2026-06-17 13:45:16.87 INFO::Fitting model to feature number 305, ASV433
#> 2026-06-17 13:45:16.87 INFO::Fitting model to feature number 306, ASV434
#> 2026-06-17 13:45:16.88 INFO::Fitting model to feature number 307, ASV438
#> 2026-06-17 13:45:16.88 INFO::Fitting model to feature number 308, ASV439
#> 2026-06-17 13:45:16.88 INFO::Fitting model to feature number 309, ASV440
#> 2026-06-17 13:45:16.89 INFO::Fitting model to feature number 310, ASV442
#> 2026-06-17 13:45:16.89 INFO::Fitting model to feature number 311, ASV443
#> 2026-06-17 13:45:16.89 INFO::Fitting model to feature number 312, ASV444
#> 2026-06-17 13:45:16.90 INFO::Fitting model to feature number 313, ASV445
#> 2026-06-17 13:45:16.90 INFO::Fitting model to feature number 314, ASV446
#> 2026-06-17 13:45:16.90 INFO::Fitting model to feature number 315, ASV448
#> 2026-06-17 13:45:16.91 INFO::Fitting model to feature number 316, ASV449
#> 2026-06-17 13:45:16.91 INFO::Fitting model to feature number 317, ASV451
#> 2026-06-17 13:45:16.91 INFO::Fitting model to feature number 318, ASV452
#> 2026-06-17 13:45:16.92 INFO::Fitting model to feature number 319, ASV453
#> 2026-06-17 13:45:16.92 INFO::Fitting model to feature number 320, ASV454
#> 2026-06-17 13:45:16.92 INFO::Fitting model to feature number 321, ASV456
#> 2026-06-17 13:45:16.93 INFO::Fitting model to feature number 322, ASV457
#> 2026-06-17 13:45:16.93 INFO::Fitting model to feature number 323, ASV458
#> 2026-06-17 13:45:16.93 INFO::Fitting model to feature number 324, ASV459
#> 2026-06-17 13:45:16.93 INFO::Fitting model to feature number 325, ASV460
#> 2026-06-17 13:45:16.94 INFO::Fitting model to feature number 326, ASV461
#> 2026-06-17 13:45:16.94 INFO::Fitting model to feature number 327, ASV462
#> 2026-06-17 13:45:16.94 INFO::Fitting model to feature number 328, ASV463
#> 2026-06-17 13:45:16.95 INFO::Fitting model to feature number 329, ASV464
#> 2026-06-17 13:45:16.95 INFO::Fitting model to feature number 330, ASV465
#> 2026-06-17 13:45:16.95 WARNING::Fitting problem for feature 330 returning NA
#> 2026-06-17 13:45:16.95 INFO::Fitting model to feature number 331, ASV466
#> 2026-06-17 13:45:16.96 INFO::Fitting model to feature number 332, ASV467
#> 2026-06-17 13:45:16.96 INFO::Fitting model to feature number 333, ASV468
#> 2026-06-17 13:45:16.97 INFO::Fitting model to feature number 334, ASV469
#> 2026-06-17 13:45:16.97 INFO::Fitting model to feature number 335, ASV470
#> 2026-06-17 13:45:16.97 INFO::Fitting model to feature number 336, ASV472
#> 2026-06-17 13:45:16.98 INFO::Fitting model to feature number 337, ASV473
#> 2026-06-17 13:45:16.98 INFO::Fitting model to feature number 338, ASV474
#> 2026-06-17 13:45:16.98 INFO::Fitting model to feature number 339, ASV475
#> 2026-06-17 13:45:16.99 INFO::Fitting model to feature number 340, ASV476
#> 2026-06-17 13:45:16.99 INFO::Fitting model to feature number 341, ASV477
#> 2026-06-17 13:45:16.99 INFO::Fitting model to feature number 342, ASV478
#> 2026-06-17 13:45:17.00 INFO::Fitting model to feature number 343, ASV479
#> 2026-06-17 13:45:17.00 INFO::Fitting model to feature number 344, ASV481
#> 2026-06-17 13:45:17.01 INFO::Fitting model to feature number 345, ASV483
#> 2026-06-17 13:45:17.01 INFO::Fitting model to feature number 346, ASV488
#> 2026-06-17 13:45:17.01 INFO::Fitting model to feature number 347, ASV489
#> 2026-06-17 13:45:17.01 INFO::Fitting model to feature number 348, ASV490
#> 2026-06-17 13:45:17.02 INFO::Fitting model to feature number 349, ASV491
#> 2026-06-17 13:45:17.02 INFO::Fitting model to feature number 350, ASV492
#> 2026-06-17 13:45:17.02 INFO::Fitting model to feature number 351, ASV493
#> 2026-06-17 13:45:17.03 INFO::Fitting model to feature number 352, ASV494
#> 2026-06-17 13:45:17.03 INFO::Fitting model to feature number 353, ASV496
#> 2026-06-17 13:45:17.03 INFO::Fitting model to feature number 354, ASV498
#> 2026-06-17 13:45:17.04 INFO::Fitting model to feature number 355, ASV499
#> 2026-06-17 13:45:17.04 INFO::Fitting model to feature number 356, ASV500
#> 2026-06-17 13:45:17.04 INFO::Fitting model to feature number 357, ASV501
#> 2026-06-17 13:45:17.05 INFO::Fitting model to feature number 358, ASV502
#> 2026-06-17 13:45:17.05 INFO::Fitting model to feature number 359, ASV504
#> 2026-06-17 13:45:17.05 INFO::Fitting model to feature number 360, ASV505
#> 2026-06-17 13:45:17.06 INFO::Fitting model to feature number 361, ASV507
#> 2026-06-17 13:45:17.06 INFO::Fitting model to feature number 362, ASV508
#> 2026-06-17 13:45:17.07 INFO::Fitting model to feature number 363, ASV509
#> 2026-06-17 13:45:17.07 INFO::Fitting model to feature number 364, ASV511
#> 2026-06-17 13:45:17.07 INFO::Fitting model to feature number 365, ASV512
#> 2026-06-17 13:45:17.08 INFO::Fitting model to feature number 366, ASV514
#> 2026-06-17 13:45:17.08 INFO::Fitting model to feature number 367, ASV515
#> 2026-06-17 13:45:17.08 INFO::Fitting model to feature number 368, ASV516
#> 2026-06-17 13:45:17.09 INFO::Fitting model to feature number 369, ASV517
#> 2026-06-17 13:45:17.09 INFO::Fitting model to feature number 370, ASV519
#> 2026-06-17 13:45:17.10 INFO::Fitting model to feature number 371, ASV520
#> 2026-06-17 13:45:17.10 INFO::Fitting model to feature number 372, ASV521
#> 2026-06-17 13:45:17.10 INFO::Fitting model to feature number 373, ASV522
#> 2026-06-17 13:45:17.10 INFO::Fitting model to feature number 374, ASV523
#> 2026-06-17 13:45:17.11 INFO::Fitting model to feature number 375, ASV526
#> 2026-06-17 13:45:17.11 INFO::Fitting model to feature number 376, ASV527
#> 2026-06-17 13:45:17.11 INFO::Fitting model to feature number 377, ASV530
#> 2026-06-17 13:45:17.12 INFO::Fitting model to feature number 378, ASV531
#> 2026-06-17 13:45:17.12 INFO::Fitting model to feature number 379, ASV533
#> 2026-06-17 13:45:17.12 INFO::Fitting model to feature number 380, ASV534
#> 2026-06-17 13:45:17.13 INFO::Fitting model to feature number 381, ASV535
#> 2026-06-17 13:45:17.13 INFO::Fitting model to feature number 382, ASV536
#> 2026-06-17 13:45:17.13 INFO::Fitting model to feature number 383, ASV538
#> 2026-06-17 13:45:17.14 INFO::Fitting model to feature number 384, ASV539
#> 2026-06-17 13:45:17.14 INFO::Fitting model to feature number 385, ASV540
#> 2026-06-17 13:45:17.14 INFO::Fitting model to feature number 386, ASV541
#> 2026-06-17 13:45:17.15 INFO::Fitting model to feature number 387, ASV542
#> 2026-06-17 13:45:17.15 INFO::Fitting model to feature number 388, ASV543
#> 2026-06-17 13:45:17.15 INFO::Fitting model to feature number 389, ASV544
#> 2026-06-17 13:45:17.16 INFO::Fitting model to feature number 390, ASV545
#> 2026-06-17 13:45:17.16 INFO::Fitting model to feature number 391, ASV546
#> 2026-06-17 13:45:17.16 INFO::Fitting model to feature number 392, ASV548
#> 2026-06-17 13:45:17.17 INFO::Fitting model to feature number 393, ASV549
#> 2026-06-17 13:45:17.17 INFO::Fitting model to feature number 394, ASV550
#> 2026-06-17 13:45:17.17 INFO::Fitting model to feature number 395, ASV551
#> 2026-06-17 13:45:17.18 INFO::Fitting model to feature number 396, ASV553
#> 2026-06-17 13:45:17.18 INFO::Fitting model to feature number 397, ASV554
#> 2026-06-17 13:45:17.18 INFO::Fitting model to feature number 398, ASV555
#> 2026-06-17 13:45:17.18 INFO::Fitting model to feature number 399, ASV556
#> 2026-06-17 13:45:17.19 INFO::Fitting model to feature number 400, ASV558
#> 2026-06-17 13:45:17.19 INFO::Fitting model to feature number 401, ASV559
#> 2026-06-17 13:45:17.19 INFO::Fitting model to feature number 402, ASV560
#> 2026-06-17 13:45:17.20 INFO::Fitting model to feature number 403, ASV561
#> 2026-06-17 13:45:17.20 INFO::Fitting model to feature number 404, ASV562
#> 2026-06-17 13:45:17.20 INFO::Fitting model to feature number 405, ASV563
#> 2026-06-17 13:45:17.21 INFO::Fitting model to feature number 406, ASV564
#> 2026-06-17 13:45:17.21 INFO::Fitting model to feature number 407, ASV566
#> 2026-06-17 13:45:17.22 INFO::Fitting model to feature number 408, ASV567
#> 2026-06-17 13:45:17.22 INFO::Fitting model to feature number 409, ASV569
#> 2026-06-17 13:45:17.22 INFO::Fitting model to feature number 410, ASV570
#> 2026-06-17 13:45:17.23 INFO::Fitting model to feature number 411, ASV571
#> 2026-06-17 13:45:17.23 INFO::Fitting model to feature number 412, ASV572
#> 2026-06-17 13:45:17.23 INFO::Fitting model to feature number 413, ASV573
#> 2026-06-17 13:45:17.23 INFO::Fitting model to feature number 414, ASV574
#> 2026-06-17 13:45:17.24 INFO::Fitting model to feature number 415, ASV576
#> 2026-06-17 13:45:17.24 INFO::Fitting model to feature number 416, ASV577
#> 2026-06-17 13:45:17.25 INFO::Fitting model to feature number 417, ASV578
#> 2026-06-17 13:45:17.25 WARNING::Fitting problem for feature 417 returning NA
#> 2026-06-17 13:45:17.25 INFO::Fitting model to feature number 418, ASV579
#> 2026-06-17 13:45:17.25 INFO::Fitting model to feature number 419, ASV580
#> 2026-06-17 13:45:17.25 INFO::Fitting model to feature number 420, ASV581
#> 2026-06-17 13:45:17.26 INFO::Fitting model to feature number 421, ASV582
#> 2026-06-17 13:45:17.26 WARNING::Fitting problem for feature 421 returning NA
#> 2026-06-17 13:45:17.26 INFO::Fitting model to feature number 422, ASV584
#> 2026-06-17 13:45:17.27 INFO::Fitting model to feature number 423, ASV585
#> 2026-06-17 13:45:17.27 INFO::Fitting model to feature number 424, ASV586
#> 2026-06-17 13:45:17.27 INFO::Fitting model to feature number 425, ASV587
#> 2026-06-17 13:45:17.28 INFO::Fitting model to feature number 426, ASV589
#> 2026-06-17 13:45:17.28 INFO::Fitting model to feature number 427, ASV590
#> 2026-06-17 13:45:17.28 INFO::Fitting model to feature number 428, ASV591
#> 2026-06-17 13:45:17.28 INFO::Fitting model to feature number 429, ASV592
#> 2026-06-17 13:45:17.32 INFO::Fitting model to feature number 430, ASV593
#> 2026-06-17 13:45:17.32 INFO::Fitting model to feature number 431, ASV594
#> 2026-06-17 13:45:17.33 INFO::Fitting model to feature number 432, ASV595
#> 2026-06-17 13:45:17.33 INFO::Fitting model to feature number 433, ASV596
#> 2026-06-17 13:45:17.33 INFO::Fitting model to feature number 434, ASV597
#> 2026-06-17 13:45:17.34 INFO::Fitting model to feature number 435, ASV598
#> 2026-06-17 13:45:17.34 INFO::Fitting model to feature number 436, ASV599
#> 2026-06-17 13:45:17.35 INFO::Fitting model to feature number 437, ASV600
#> 2026-06-17 13:45:17.35 INFO::Fitting model to feature number 438, ASV602
#> 2026-06-17 13:45:17.35 INFO::Fitting model to feature number 439, ASV604
#> 2026-06-17 13:45:17.36 INFO::Fitting model to feature number 440, ASV605
#> 2026-06-17 13:45:17.36 INFO::Fitting model to feature number 441, ASV607
#> 2026-06-17 13:45:17.36 INFO::Fitting model to feature number 442, ASV608
#> 2026-06-17 13:45:17.37 INFO::Fitting model to feature number 443, ASV610
#> 2026-06-17 13:45:17.37 INFO::Fitting model to feature number 444, ASV612
#> 2026-06-17 13:45:17.37 INFO::Fitting model to feature number 445, ASV614
#> 2026-06-17 13:45:17.38 INFO::Fitting model to feature number 446, ASV615
#> 2026-06-17 13:45:17.38 INFO::Fitting model to feature number 447, ASV616
#> 2026-06-17 13:45:17.39 INFO::Fitting model to feature number 448, ASV617
#> 2026-06-17 13:45:17.39 INFO::Fitting model to feature number 449, ASV618
#> 2026-06-17 13:45:17.39 INFO::Fitting model to feature number 450, ASV621
#> 2026-06-17 13:45:17.40 INFO::Fitting model to feature number 451, ASV624
#> 2026-06-17 13:45:17.40 INFO::Fitting model to feature number 452, ASV625
#> 2026-06-17 13:45:17.40 INFO::Fitting model to feature number 453, ASV626
#> 2026-06-17 13:45:17.41 INFO::Fitting model to feature number 454, ASV627
#> 2026-06-17 13:45:17.41 INFO::Fitting model to feature number 455, ASV628
#> 2026-06-17 13:45:17.42 INFO::Fitting model to feature number 456, ASV629
#> 2026-06-17 13:45:17.42 INFO::Fitting model to feature number 457, ASV631
#> 2026-06-17 13:45:17.42 INFO::Fitting model to feature number 458, ASV632
#> 2026-06-17 13:45:17.42 INFO::Fitting model to feature number 459, ASV633
#> 2026-06-17 13:45:17.43 INFO::Fitting model to feature number 460, ASV634
#> 2026-06-17 13:45:17.43 INFO::Fitting model to feature number 461, ASV635
#> 2026-06-17 13:45:17.43 INFO::Fitting model to feature number 462, ASV637
#> 2026-06-17 13:45:17.44 INFO::Fitting model to feature number 463, ASV639
#> 2026-06-17 13:45:17.44 INFO::Fitting model to feature number 464, ASV640
#> 2026-06-17 13:45:17.44 INFO::Fitting model to feature number 465, ASV641
#> 2026-06-17 13:45:17.45 INFO::Fitting model to feature number 466, ASV643
#> 2026-06-17 13:45:17.45 INFO::Fitting model to feature number 467, ASV644
#> 2026-06-17 13:45:17.45 INFO::Fitting model to feature number 468, ASV645
#> 2026-06-17 13:45:17.46 INFO::Fitting model to feature number 469, ASV647
#> 2026-06-17 13:45:17.46 INFO::Fitting model to feature number 470, ASV649
#> 2026-06-17 13:45:17.46 INFO::Fitting model to feature number 471, ASV650
#> 2026-06-17 13:45:17.47 INFO::Fitting model to feature number 472, ASV652
#> 2026-06-17 13:45:17.47 INFO::Fitting model to feature number 473, ASV653
#> 2026-06-17 13:45:17.47 INFO::Fitting model to feature number 474, ASV654
#> 2026-06-17 13:45:17.47 INFO::Fitting model to feature number 475, ASV657
#> 2026-06-17 13:45:17.47 INFO::Fitting model to feature number 476, ASV659
#> 2026-06-17 13:45:17.48 INFO::Fitting model to feature number 477, ASV660
#> 2026-06-17 13:45:17.48 INFO::Fitting model to feature number 478, ASV661
#> 2026-06-17 13:45:17.49 INFO::Fitting model to feature number 479, ASV662
#> 2026-06-17 13:45:17.49 INFO::Fitting model to feature number 480, ASV664
#> 2026-06-17 13:45:17.49 INFO::Fitting model to feature number 481, ASV665
#> 2026-06-17 13:45:17.49 INFO::Fitting model to feature number 482, ASV666
#> 2026-06-17 13:45:17.50 INFO::Fitting model to feature number 483, ASV667
#> 2026-06-17 13:45:17.50 INFO::Fitting model to feature number 484, ASV668
#> 2026-06-17 13:45:17.51 INFO::Fitting model to feature number 485, ASV669
#> 2026-06-17 13:45:17.51 INFO::Fitting model to feature number 486, ASV670
#> 2026-06-17 13:45:17.51 INFO::Fitting model to feature number 487, ASV671
#> 2026-06-17 13:45:17.52 INFO::Fitting model to feature number 488, ASV672
#> 2026-06-17 13:45:17.52 INFO::Fitting model to feature number 489, ASV673
#> 2026-06-17 13:45:17.52 INFO::Fitting model to feature number 490, ASV674
#> 2026-06-17 13:45:17.53 INFO::Fitting model to feature number 491, ASV675
#> 2026-06-17 13:45:17.53 INFO::Fitting model to feature number 492, ASV676
#> 2026-06-17 13:45:17.53 INFO::Fitting model to feature number 493, ASV677
#> 2026-06-17 13:45:17.54 INFO::Fitting model to feature number 494, ASV678
#> 2026-06-17 13:45:17.54 INFO::Fitting model to feature number 495, ASV679
#> 2026-06-17 13:45:17.54 INFO::Fitting model to feature number 496, ASV680
#> 2026-06-17 13:45:17.55 INFO::Fitting model to feature number 497, ASV683
#> 2026-06-17 13:45:17.55 INFO::Fitting model to feature number 498, ASV684
#> 2026-06-17 13:45:17.55 INFO::Fitting model to feature number 499, ASV685
#> 2026-06-17 13:45:17.55 INFO::Fitting model to feature number 500, ASV686
#> 2026-06-17 13:45:17.55 INFO::Fitting model to feature number 501, ASV687
#> 2026-06-17 13:45:17.56 INFO::Fitting model to feature number 502, ASV688
#> 2026-06-17 13:45:17.56 INFO::Fitting model to feature number 503, ASV691
#> 2026-06-17 13:45:17.56 INFO::Fitting model to feature number 504, ASV692
#> 2026-06-17 13:45:17.56 INFO::Fitting model to feature number 505, ASV693
#> 2026-06-17 13:45:17.57 INFO::Fitting model to feature number 506, ASV694
#> 2026-06-17 13:45:17.57 INFO::Fitting model to feature number 507, ASV695
#> 2026-06-17 13:45:17.57 INFO::Fitting model to feature number 508, ASV696
#> 2026-06-17 13:45:17.58 INFO::Fitting model to feature number 509, ASV697
#> 2026-06-17 13:45:17.58 INFO::Fitting model to feature number 510, ASV698
#> 2026-06-17 13:45:17.59 INFO::Fitting model to feature number 511, ASV700
#> 2026-06-17 13:45:17.59 INFO::Fitting model to feature number 512, ASV702
#> 2026-06-17 13:45:17.59 INFO::Fitting model to feature number 513, ASV703
#> 2026-06-17 13:45:17.60 INFO::Fitting model to feature number 514, ASV704
#> 2026-06-17 13:45:17.60 INFO::Fitting model to feature number 515, ASV705
#> 2026-06-17 13:45:17.60 INFO::Fitting model to feature number 516, ASV707
#> 2026-06-17 13:45:17.61 INFO::Fitting model to feature number 517, ASV710
#> 2026-06-17 13:45:17.61 INFO::Fitting model to feature number 518, ASV711
#> 2026-06-17 13:45:17.61 INFO::Fitting model to feature number 519, ASV712
#> 2026-06-17 13:45:17.62 INFO::Fitting model to feature number 520, ASV713
#> 2026-06-17 13:45:17.62 INFO::Fitting model to feature number 521, ASV714
#> 2026-06-17 13:45:17.62 INFO::Fitting model to feature number 522, ASV715
#> 2026-06-17 13:45:17.63 INFO::Fitting model to feature number 523, ASV717
#> 2026-06-17 13:45:17.63 INFO::Fitting model to feature number 524, ASV718
#> 2026-06-17 13:45:17.63 INFO::Fitting model to feature number 525, ASV720
#> 2026-06-17 13:45:17.64 INFO::Fitting model to feature number 526, ASV722
#> 2026-06-17 13:45:17.64 INFO::Fitting model to feature number 527, ASV723
#> 2026-06-17 13:45:17.64 INFO::Fitting model to feature number 528, ASV724
#> 2026-06-17 13:45:17.64 INFO::Fitting model to feature number 529, ASV726
#> 2026-06-17 13:45:17.65 INFO::Fitting model to feature number 530, ASV727
#> 2026-06-17 13:45:17.65 INFO::Fitting model to feature number 531, ASV728
#> 2026-06-17 13:45:17.65 WARNING::Fitting problem for feature 531 returning NA
#> 2026-06-17 13:45:17.66 INFO::Fitting model to feature number 532, ASV729
#> 2026-06-17 13:45:17.66 INFO::Fitting model to feature number 533, ASV730
#> 2026-06-17 13:45:17.66 INFO::Fitting model to feature number 534, ASV731
#> 2026-06-17 13:45:17.67 INFO::Fitting model to feature number 535, ASV732
#> 2026-06-17 13:45:17.67 INFO::Fitting model to feature number 536, ASV733
#> 2026-06-17 13:45:17.67 INFO::Fitting model to feature number 537, ASV734
#> 2026-06-17 13:45:17.68 INFO::Fitting model to feature number 538, ASV736
#> 2026-06-17 13:45:17.68 INFO::Fitting model to feature number 539, ASV737
#> 2026-06-17 13:45:17.68 INFO::Fitting model to feature number 540, ASV738
#> 2026-06-17 13:45:17.68 INFO::Fitting model to feature number 541, ASV740
#> 2026-06-17 13:45:17.69 INFO::Fitting model to feature number 542, ASV741
#> 2026-06-17 13:45:17.69 INFO::Fitting model to feature number 543, ASV742
#> 2026-06-17 13:45:17.70 INFO::Fitting model to feature number 544, ASV743
#> 2026-06-17 13:45:17.70 INFO::Fitting model to feature number 545, ASV744
#> 2026-06-17 13:45:17.70 INFO::Fitting model to feature number 546, ASV746
#> 2026-06-17 13:45:17.70 INFO::Fitting model to feature number 547, ASV747
#> 2026-06-17 13:45:17.70 INFO::Fitting model to feature number 548, ASV748
#> 2026-06-17 13:45:17.71 INFO::Fitting model to feature number 549, ASV749
#> 2026-06-17 13:45:17.71 INFO::Fitting model to feature number 550, ASV752
#> 2026-06-17 13:45:17.71 INFO::Fitting model to feature number 551, ASV753
#> 2026-06-17 13:45:17.72 INFO::Fitting model to feature number 552, ASV754
#> 2026-06-17 13:45:17.72 INFO::Fitting model to feature number 553, ASV755
#> 2026-06-17 13:45:17.72 INFO::Fitting model to feature number 554, ASV756
#> 2026-06-17 13:45:17.73 INFO::Fitting model to feature number 555, ASV757
#> 2026-06-17 13:45:17.73 INFO::Fitting model to feature number 556, ASV758
#> 2026-06-17 13:45:17.73 INFO::Fitting model to feature number 557, ASV760
#> 2026-06-17 13:45:17.74 INFO::Fitting model to feature number 558, ASV762
#> 2026-06-17 13:45:17.74 INFO::Fitting model to feature number 559, ASV764
#> 2026-06-17 13:45:17.74 INFO::Fitting model to feature number 560, ASV766
#> 2026-06-17 13:45:17.75 INFO::Fitting model to feature number 561, ASV767
#> 2026-06-17 13:45:17.75 INFO::Fitting model to feature number 562, ASV768
#> 2026-06-17 13:45:17.75 INFO::Fitting model to feature number 563, ASV769
#> 2026-06-17 13:45:17.76 INFO::Fitting model to feature number 564, ASV770
#> 2026-06-17 13:45:17.76 INFO::Fitting model to feature number 565, ASV771
#> 2026-06-17 13:45:17.76 INFO::Fitting model to feature number 566, ASV772
#> 2026-06-17 13:45:17.77 INFO::Fitting model to feature number 567, ASV773
#> 2026-06-17 13:45:17.77 INFO::Fitting model to feature number 568, ASV774
#> 2026-06-17 13:45:17.77 INFO::Fitting model to feature number 569, ASV775
#> 2026-06-17 13:45:17.78 INFO::Fitting model to feature number 570, ASV776
#> 2026-06-17 13:45:17.78 INFO::Fitting model to feature number 571, ASV777
#> 2026-06-17 13:45:17.78 INFO::Fitting model to feature number 572, ASV779
#> 2026-06-17 13:45:17.79 INFO::Fitting model to feature number 573, ASV780
#> 2026-06-17 13:45:17.79 INFO::Fitting model to feature number 574, ASV781
#> 2026-06-17 13:45:17.80 INFO::Fitting model to feature number 575, ASV782
#> 2026-06-17 13:45:17.80 INFO::Fitting model to feature number 576, ASV783
#> 2026-06-17 13:45:17.81 INFO::Fitting model to feature number 577, ASV784
#> 2026-06-17 13:45:17.81 INFO::Fitting model to feature number 578, ASV785
#> 2026-06-17 13:45:17.81 INFO::Fitting model to feature number 579, ASV786
#> 2026-06-17 13:45:17.82 WARNING::Fitting problem for feature 579 returning NA
#> 2026-06-17 13:45:17.82 INFO::Fitting model to feature number 580, ASV787
#> 2026-06-17 13:45:17.82 INFO::Fitting model to feature number 581, ASV788
#> 2026-06-17 13:45:17.83 INFO::Fitting model to feature number 582, ASV790
#> 2026-06-17 13:45:17.83 INFO::Fitting model to feature number 583, ASV792
#> 2026-06-17 13:45:17.83 INFO::Fitting model to feature number 584, ASV795
#> 2026-06-17 13:45:17.84 INFO::Fitting model to feature number 585, ASV796
#> 2026-06-17 13:45:17.84 INFO::Fitting model to feature number 586, ASV797
#> 2026-06-17 13:45:17.84 INFO::Fitting model to feature number 587, ASV798
#> 2026-06-17 13:45:17.85 INFO::Fitting model to feature number 588, ASV799
#> 2026-06-17 13:45:17.85 INFO::Fitting model to feature number 589, ASV801
#> 2026-06-17 13:45:17.85 INFO::Fitting model to feature number 590, ASV802
#> 2026-06-17 13:45:17.85 INFO::Fitting model to feature number 591, ASV803
#> 2026-06-17 13:45:17.86 INFO::Fitting model to feature number 592, ASV805
#> 2026-06-17 13:45:17.86 INFO::Fitting model to feature number 593, ASV807
#> 2026-06-17 13:45:17.86 INFO::Fitting model to feature number 594, ASV808
#> 2026-06-17 13:45:17.87 INFO::Fitting model to feature number 595, ASV810
#> 2026-06-17 13:45:17.87 INFO::Fitting model to feature number 596, ASV811
#> 2026-06-17 13:45:17.87 INFO::Fitting model to feature number 597, ASV814
#> 2026-06-17 13:45:17.87 INFO::Fitting model to feature number 598, ASV815
#> 2026-06-17 13:45:17.88 INFO::Fitting model to feature number 599, ASV816
#> 2026-06-17 13:45:17.88 INFO::Fitting model to feature number 600, ASV817
#> 2026-06-17 13:45:17.88 INFO::Fitting model to feature number 601, ASV819
#> 2026-06-17 13:45:17.89 INFO::Fitting model to feature number 602, ASV821
#> 2026-06-17 13:45:17.89 INFO::Fitting model to feature number 603, ASV822
#> 2026-06-17 13:45:17.90 INFO::Fitting model to feature number 604, ASV823
#> 2026-06-17 13:45:17.90 INFO::Fitting model to feature number 605, ASV824
#> 2026-06-17 13:45:17.90 INFO::Fitting model to feature number 606, ASV828
#> 2026-06-17 13:45:17.91 INFO::Fitting model to feature number 607, ASV829
#> 2026-06-17 13:45:17.91 INFO::Fitting model to feature number 608, ASV830
#> 2026-06-17 13:45:17.91 INFO::Fitting model to feature number 609, ASV831
#> 2026-06-17 13:45:17.92 INFO::Fitting model to feature number 610, ASV832
#> 2026-06-17 13:45:17.92 INFO::Fitting model to feature number 611, ASV834
#> 2026-06-17 13:45:17.92 INFO::Fitting model to feature number 612, ASV836
#> 2026-06-17 13:45:17.93 INFO::Fitting model to feature number 613, ASV837
#> 2026-06-17 13:45:17.93 INFO::Fitting model to feature number 614, ASV838
#> 2026-06-17 13:45:17.93 INFO::Fitting model to feature number 615, ASV839
#> 2026-06-17 13:45:17.94 INFO::Fitting model to feature number 616, ASV840
#> 2026-06-17 13:45:17.94 INFO::Fitting model to feature number 617, ASV841
#> 2026-06-17 13:45:17.94 INFO::Fitting model to feature number 618, ASV842
#> 2026-06-17 13:45:17.95 INFO::Fitting model to feature number 619, ASV843
#> 2026-06-17 13:45:17.95 INFO::Fitting model to feature number 620, ASV844
#> 2026-06-17 13:45:17.95 INFO::Fitting model to feature number 621, ASV845
#> 2026-06-17 13:45:17.96 INFO::Fitting model to feature number 622, ASV847
#> 2026-06-17 13:45:17.96 WARNING::Fitting problem for feature 622 returning NA
#> 2026-06-17 13:45:17.96 INFO::Fitting model to feature number 623, ASV848
#> 2026-06-17 13:45:17.96 INFO::Fitting model to feature number 624, ASV852
#> 2026-06-17 13:45:17.97 INFO::Fitting model to feature number 625, ASV853
#> 2026-06-17 13:45:17.97 INFO::Fitting model to feature number 626, ASV854
#> 2026-06-17 13:45:17.97 INFO::Fitting model to feature number 627, ASV855
#> 2026-06-17 13:45:17.98 INFO::Fitting model to feature number 628, ASV857
#> 2026-06-17 13:45:17.98 INFO::Fitting model to feature number 629, ASV858
#> 2026-06-17 13:45:17.99 INFO::Fitting model to feature number 630, ASV859
#> 2026-06-17 13:45:17.99 INFO::Fitting model to feature number 631, ASV860
#> 2026-06-17 13:45:17.99 INFO::Fitting model to feature number 632, ASV861
#> 2026-06-17 13:45:18.00 WARNING::Fitting problem for feature 632 returning NA
#> 2026-06-17 13:45:18.00 INFO::Fitting model to feature number 633, ASV863
#> 2026-06-17 13:45:18.00 INFO::Fitting model to feature number 634, ASV865
#> 2026-06-17 13:45:18.00 INFO::Fitting model to feature number 635, ASV870
#> 2026-06-17 13:45:18.01 INFO::Fitting model to feature number 636, ASV873
#> 2026-06-17 13:45:18.01 INFO::Fitting model to feature number 637, ASV874
#> 2026-06-17 13:45:18.01 INFO::Fitting model to feature number 638, ASV875
#> 2026-06-17 13:45:18.02 INFO::Fitting model to feature number 639, ASV876
#> 2026-06-17 13:45:18.02 INFO::Fitting model to feature number 640, ASV877
#> 2026-06-17 13:45:18.03 INFO::Fitting model to feature number 641, ASV878
#> 2026-06-17 13:45:18.03 INFO::Fitting model to feature number 642, ASV879
#> 2026-06-17 13:45:18.03 INFO::Fitting model to feature number 643, ASV880
#> 2026-06-17 13:45:18.03 INFO::Fitting model to feature number 644, ASV883
#> 2026-06-17 13:45:18.04 INFO::Fitting model to feature number 645, ASV884
#> 2026-06-17 13:45:18.04 INFO::Fitting model to feature number 646, ASV887
#> 2026-06-17 13:45:18.04 INFO::Fitting model to feature number 647, ASV888
#> 2026-06-17 13:45:18.05 INFO::Fitting model to feature number 648, ASV890
#> 2026-06-17 13:45:18.05 INFO::Fitting model to feature number 649, ASV891
#> 2026-06-17 13:45:18.05 INFO::Fitting model to feature number 650, ASV892
#> 2026-06-17 13:45:18.06 INFO::Fitting model to feature number 651, ASV895
#> 2026-06-17 13:45:18.06 INFO::Fitting model to feature number 652, ASV897
#> 2026-06-17 13:45:18.06 INFO::Fitting model to feature number 653, ASV898
#> 2026-06-17 13:45:18.07 INFO::Fitting model to feature number 654, ASV899
#> 2026-06-17 13:45:18.07 INFO::Fitting model to feature number 655, ASV901
#> 2026-06-17 13:45:18.07 INFO::Fitting model to feature number 656, ASV902
#> 2026-06-17 13:45:18.08 INFO::Fitting model to feature number 657, ASV903
#> 2026-06-17 13:45:18.08 WARNING::Fitting problem for feature 657 returning NA
#> 2026-06-17 13:45:18.08 INFO::Fitting model to feature number 658, ASV904
#> 2026-06-17 13:45:18.08 INFO::Fitting model to feature number 659, ASV905
#> 2026-06-17 13:45:18.09 INFO::Fitting model to feature number 660, ASV906
#> 2026-06-17 13:45:18.09 INFO::Fitting model to feature number 661, ASV907
#> 2026-06-17 13:45:18.09 WARNING::Fitting problem for feature 661 returning NA
#> 2026-06-17 13:45:18.10 INFO::Fitting model to feature number 662, ASV909
#> 2026-06-17 13:45:18.10 INFO::Fitting model to feature number 663, ASV910
#> 2026-06-17 13:45:18.10 INFO::Fitting model to feature number 664, ASV911
#> 2026-06-17 13:45:18.10 INFO::Fitting model to feature number 665, ASV913
#> 2026-06-17 13:45:18.11 INFO::Fitting model to feature number 666, ASV914
#> 2026-06-17 13:45:18.11 INFO::Fitting model to feature number 667, ASV915
#> 2026-06-17 13:45:18.11 INFO::Fitting model to feature number 668, ASV916
#> 2026-06-17 13:45:18.12 INFO::Fitting model to feature number 669, ASV917
#> 2026-06-17 13:45:18.12 INFO::Fitting model to feature number 670, ASV918
#> 2026-06-17 13:45:18.12 INFO::Fitting model to feature number 671, ASV919
#> 2026-06-17 13:45:18.13 INFO::Fitting model to feature number 672, ASV920
#> 2026-06-17 13:45:18.13 INFO::Fitting model to feature number 673, ASV921
#> 2026-06-17 13:45:18.13 INFO::Fitting model to feature number 674, ASV922
#> 2026-06-17 13:45:18.14 INFO::Fitting model to feature number 675, ASV923
#> 2026-06-17 13:45:18.14 INFO::Fitting model to feature number 676, ASV926
#> 2026-06-17 13:45:18.15 INFO::Fitting model to feature number 677, ASV927
#> 2026-06-17 13:45:18.15 INFO::Fitting model to feature number 678, ASV928
#> 2026-06-17 13:45:18.15 INFO::Fitting model to feature number 679, ASV930
#> 2026-06-17 13:45:18.16 INFO::Fitting model to feature number 680, ASV932
#> 2026-06-17 13:45:18.16 INFO::Fitting model to feature number 681, ASV934
#> 2026-06-17 13:45:18.16 INFO::Fitting model to feature number 682, ASV935
#> 2026-06-17 13:45:18.16 INFO::Fitting model to feature number 683, ASV936
#> 2026-06-17 13:45:18.17 INFO::Fitting model to feature number 684, ASV937
#> 2026-06-17 13:45:18.17 INFO::Fitting model to feature number 685, ASV938
#> 2026-06-17 13:45:18.18 INFO::Fitting model to feature number 686, ASV939
#> 2026-06-17 13:45:18.18 INFO::Fitting model to feature number 687, ASV940
#> 2026-06-17 13:45:18.18 INFO::Fitting model to feature number 688, ASV941
#> 2026-06-17 13:45:18.19 INFO::Fitting model to feature number 689, ASV942
#> 2026-06-17 13:45:18.19 INFO::Fitting model to feature number 690, ASV943
#> 2026-06-17 13:45:18.19 INFO::Fitting model to feature number 691, ASV945
#> 2026-06-17 13:45:18.20 INFO::Fitting model to feature number 692, ASV947
#> 2026-06-17 13:45:18.20 INFO::Fitting model to feature number 693, ASV948
#> 2026-06-17 13:45:18.21 INFO::Fitting model to feature number 694, ASV949
#> 2026-06-17 13:45:18.21 INFO::Fitting model to feature number 695, ASV950
#> 2026-06-17 13:45:18.21 INFO::Fitting model to feature number 696, ASV951
#> 2026-06-17 13:45:18.21 INFO::Fitting model to feature number 697, ASV953
#> 2026-06-17 13:45:18.22 INFO::Fitting model to feature number 698, ASV955
#> 2026-06-17 13:45:18.22 INFO::Fitting model to feature number 699, ASV958
#> 2026-06-17 13:45:18.22 INFO::Fitting model to feature number 700, ASV959
#> 2026-06-17 13:45:18.22 INFO::Fitting model to feature number 701, ASV961
#> 2026-06-17 13:45:18.23 INFO::Fitting model to feature number 702, ASV962
#> 2026-06-17 13:45:18.23 INFO::Fitting model to feature number 703, ASV963
#> 2026-06-17 13:45:18.23 INFO::Fitting model to feature number 704, ASV964
#> 2026-06-17 13:45:18.24 INFO::Fitting model to feature number 705, ASV966
#> 2026-06-17 13:45:18.24 INFO::Fitting model to feature number 706, ASV967
#> 2026-06-17 13:45:18.24 INFO::Fitting model to feature number 707, ASV969
#> 2026-06-17 13:45:18.25 INFO::Fitting model to feature number 708, ASV970
#> 2026-06-17 13:45:18.25 INFO::Fitting model to feature number 709, ASV971
#> 2026-06-17 13:45:18.25 INFO::Fitting model to feature number 710, ASV972
#> 2026-06-17 13:45:18.25 INFO::Fitting model to feature number 711, ASV973
#> 2026-06-17 13:45:18.26 INFO::Fitting model to feature number 712, ASV974
#> 2026-06-17 13:45:18.26 INFO::Fitting model to feature number 713, ASV975
#> 2026-06-17 13:45:18.27 INFO::Fitting model to feature number 714, ASV976
#> 2026-06-17 13:45:18.27 INFO::Fitting model to feature number 715, ASV977
#> 2026-06-17 13:45:18.27 INFO::Fitting model to feature number 716, ASV979
#> 2026-06-17 13:45:18.27 INFO::Fitting model to feature number 717, ASV980
#> 2026-06-17 13:45:18.28 INFO::Fitting model to feature number 718, ASV981
#> 2026-06-17 13:45:18.28 INFO::Fitting model to feature number 719, ASV983
#> 2026-06-17 13:45:18.28 INFO::Fitting model to feature number 720, ASV984
#> 2026-06-17 13:45:18.29 INFO::Fitting model to feature number 721, ASV986
#> 2026-06-17 13:45:18.29 INFO::Fitting model to feature number 722, ASV987
#> 2026-06-17 13:45:18.29 INFO::Fitting model to feature number 723, ASV988
#> 2026-06-17 13:45:18.30 INFO::Fitting model to feature number 724, ASV989
#> 2026-06-17 13:45:18.30 INFO::Fitting model to feature number 725, ASV990
#> 2026-06-17 13:45:18.31 INFO::Fitting model to feature number 726, ASV992
#> 2026-06-17 13:45:18.31 WARNING::Fitting problem for feature 726 returning NA
#> 2026-06-17 13:45:18.31 INFO::Fitting model to feature number 727, ASV993
#> 2026-06-17 13:45:18.31 INFO::Fitting model to feature number 728, ASV994
#> 2026-06-17 13:45:18.32 INFO::Fitting model to feature number 729, ASV995
#> 2026-06-17 13:45:18.32 INFO::Fitting model to feature number 730, ASV996
#> 2026-06-17 13:45:18.32 INFO::Fitting model to feature number 731, ASV997
#> 2026-06-17 13:45:18.32 INFO::Fitting model to feature number 732, ASV998
#> 2026-06-17 13:45:18.33 INFO::Fitting model to feature number 733, ASV999
#> 2026-06-17 13:45:18.33 INFO::Fitting model to feature number 734, ASV1001
#> 2026-06-17 13:45:18.33 INFO::Fitting model to feature number 735, ASV1002
#> 2026-06-17 13:45:18.34 INFO::Fitting model to feature number 736, ASV1003
#> 2026-06-17 13:45:18.34 INFO::Fitting model to feature number 737, ASV1004
#> 2026-06-17 13:45:18.35 INFO::Fitting model to feature number 738, ASV1005
#> 2026-06-17 13:45:18.35 INFO::Fitting model to feature number 739, ASV1006
#> 2026-06-17 13:45:18.35 INFO::Fitting model to feature number 740, ASV1007
#> 2026-06-17 13:45:18.36 INFO::Fitting model to feature number 741, ASV1008
#> 2026-06-17 13:45:18.36 INFO::Fitting model to feature number 742, ASV1010
#> 2026-06-17 13:45:18.36 INFO::Fitting model to feature number 743, ASV1011
#> 2026-06-17 13:45:18.37 INFO::Fitting model to feature number 744, ASV1014
#> 2026-06-17 13:45:18.37 INFO::Fitting model to feature number 745, ASV1015
#> 2026-06-17 13:45:18.38 WARNING::Fitting problem for feature 745 returning NA
#> 2026-06-17 13:45:18.38 INFO::Fitting model to feature number 746, ASV1016
#> 2026-06-17 13:45:18.38 INFO::Fitting model to feature number 747, ASV1017
#> 2026-06-17 13:45:18.38 INFO::Fitting model to feature number 748, ASV1018
#> 2026-06-17 13:45:18.39 INFO::Fitting model to feature number 749, ASV1020
#> 2026-06-17 13:45:18.39 INFO::Fitting model to feature number 750, ASV1021
#> 2026-06-17 13:45:18.39 INFO::Fitting model to feature number 751, ASV1022
#> 2026-06-17 13:45:18.40 INFO::Fitting model to feature number 752, ASV1023
#> 2026-06-17 13:45:18.40 INFO::Fitting model to feature number 753, ASV1024
#> 2026-06-17 13:45:18.40 INFO::Fitting model to feature number 754, ASV1025
#> 2026-06-17 13:45:18.40 INFO::Fitting model to feature number 755, ASV1026
#> 2026-06-17 13:45:18.40 INFO::Fitting model to feature number 756, ASV1027
#> 2026-06-17 13:45:18.41 INFO::Fitting model to feature number 757, ASV1028
#> 2026-06-17 13:45:18.41 INFO::Fitting model to feature number 758, ASV1029
#> 2026-06-17 13:45:18.41 INFO::Fitting model to feature number 759, ASV1030
#> 2026-06-17 13:45:18.42 INFO::Fitting model to feature number 760, ASV1031
#> 2026-06-17 13:45:18.42 INFO::Fitting model to feature number 761, ASV1032
#> 2026-06-17 13:45:18.42 INFO::Fitting model to feature number 762, ASV1035
#> 2026-06-17 13:45:18.43 INFO::Fitting model to feature number 763, ASV1036
#> 2026-06-17 13:45:18.43 INFO::Fitting model to feature number 764, ASV1037
#> 2026-06-17 13:45:18.43 INFO::Fitting model to feature number 765, ASV1038
#> 2026-06-17 13:45:18.44 INFO::Fitting model to feature number 766, ASV1039
#> 2026-06-17 13:45:18.44 INFO::Fitting model to feature number 767, ASV1040
#> 2026-06-17 13:45:18.44 INFO::Fitting model to feature number 768, ASV1041
#> 2026-06-17 13:45:18.44 INFO::Fitting model to feature number 769, ASV1042
#> 2026-06-17 13:45:18.45 INFO::Fitting model to feature number 770, ASV1044
#> 2026-06-17 13:45:18.45 INFO::Fitting model to feature number 771, ASV1046
#> 2026-06-17 13:45:18.45 INFO::Fitting model to feature number 772, ASV1047
#> 2026-06-17 13:45:18.46 INFO::Fitting model to feature number 773, ASV1048
#> 2026-06-17 13:45:18.46 INFO::Fitting model to feature number 774, ASV1049
#> 2026-06-17 13:45:18.47 INFO::Fitting model to feature number 775, ASV1052
#> 2026-06-17 13:45:18.47 INFO::Fitting model to feature number 776, ASV1053
#> 2026-06-17 13:45:18.47 INFO::Fitting model to feature number 777, ASV1057
#> 2026-06-17 13:45:18.48 INFO::Fitting model to feature number 778, ASV1058
#> 2026-06-17 13:45:18.48 INFO::Fitting model to feature number 779, ASV1059
#> 2026-06-17 13:45:18.49 INFO::Fitting model to feature number 780, ASV1063
#> 2026-06-17 13:45:18.49 INFO::Fitting model to feature number 781, ASV1065
#> 2026-06-17 13:45:18.49 INFO::Fitting model to feature number 782, ASV1066
#> 2026-06-17 13:45:18.50 INFO::Fitting model to feature number 783, ASV1067
#> 2026-06-17 13:45:18.50 INFO::Fitting model to feature number 784, ASV1068
#> 2026-06-17 13:45:18.50 INFO::Fitting model to feature number 785, ASV1069
#> 2026-06-17 13:45:18.51 INFO::Fitting model to feature number 786, ASV1070
#> 2026-06-17 13:45:18.51 INFO::Fitting model to feature number 787, ASV1071
#> 2026-06-17 13:45:18.51 INFO::Fitting model to feature number 788, ASV1072
#> 2026-06-17 13:45:18.51 INFO::Fitting model to feature number 789, ASV1073
#> 2026-06-17 13:45:18.52 INFO::Fitting model to feature number 790, ASV1074
#> 2026-06-17 13:45:18.52 INFO::Fitting model to feature number 791, ASV1076
#> 2026-06-17 13:45:18.52 INFO::Fitting model to feature number 792, ASV1078
#> 2026-06-17 13:45:18.53 INFO::Fitting model to feature number 793, ASV1079
#> 2026-06-17 13:45:18.53 INFO::Fitting model to feature number 794, ASV1082
#> 2026-06-17 13:45:18.53 INFO::Fitting model to feature number 795, ASV1084
#> 2026-06-17 13:45:18.54 INFO::Fitting model to feature number 796, ASV1085
#> 2026-06-17 13:45:18.54 INFO::Fitting model to feature number 797, ASV1086
#> 2026-06-17 13:45:18.54 INFO::Fitting model to feature number 798, ASV1087
#> 2026-06-17 13:45:18.55 INFO::Fitting model to feature number 799, ASV1090
#> 2026-06-17 13:45:18.55 INFO::Fitting model to feature number 800, ASV1091
#> 2026-06-17 13:45:18.55 WARNING::Fitting problem for feature 800 returning NA
#> 2026-06-17 13:45:18.55 INFO::Fitting model to feature number 801, ASV1093
#> 2026-06-17 13:45:18.56 INFO::Fitting model to feature number 802, ASV1095
#> 2026-06-17 13:45:18.56 INFO::Fitting model to feature number 803, ASV1096
#> 2026-06-17 13:45:18.57 INFO::Fitting model to feature number 804, ASV1097
#> 2026-06-17 13:45:18.57 INFO::Fitting model to feature number 805, ASV1099
#> 2026-06-17 13:45:18.57 INFO::Fitting model to feature number 806, ASV1100
#> 2026-06-17 13:45:18.58 INFO::Fitting model to feature number 807, ASV1101
#> 2026-06-17 13:45:18.58 INFO::Fitting model to feature number 808, ASV1103
#> 2026-06-17 13:45:18.58 INFO::Fitting model to feature number 809, ASV1105
#> 2026-06-17 13:45:18.59 INFO::Fitting model to feature number 810, ASV1107
#> 2026-06-17 13:45:18.59 INFO::Fitting model to feature number 811, ASV1108
#> 2026-06-17 13:45:18.59 INFO::Fitting model to feature number 812, ASV1109
#> 2026-06-17 13:45:18.60 INFO::Fitting model to feature number 813, ASV1110
#> 2026-06-17 13:45:18.60 INFO::Fitting model to feature number 814, ASV1111
#> 2026-06-17 13:45:18.60 INFO::Fitting model to feature number 815, ASV1112
#> 2026-06-17 13:45:18.61 INFO::Fitting model to feature number 816, ASV1114
#> 2026-06-17 13:45:18.61 INFO::Fitting model to feature number 817, ASV1115
#> 2026-06-17 13:45:18.61 INFO::Fitting model to feature number 818, ASV1116
#> 2026-06-17 13:45:18.62 INFO::Fitting model to feature number 819, ASV1120
#> 2026-06-17 13:45:18.62 INFO::Fitting model to feature number 820, ASV1121
#> 2026-06-17 13:45:18.62 INFO::Fitting model to feature number 821, ASV1122
#> 2026-06-17 13:45:18.63 INFO::Fitting model to feature number 822, ASV1125
#> 2026-06-17 13:45:18.63 INFO::Fitting model to feature number 823, ASV1126
#> 2026-06-17 13:45:18.63 INFO::Fitting model to feature number 824, ASV1127
#> 2026-06-17 13:45:18.64 INFO::Fitting model to feature number 825, ASV1128
#> 2026-06-17 13:45:18.64 INFO::Fitting model to feature number 826, ASV1129
#> 2026-06-17 13:45:18.64 INFO::Fitting model to feature number 827, ASV1131
#> 2026-06-17 13:45:18.64 INFO::Fitting model to feature number 828, ASV1132
#> 2026-06-17 13:45:18.65 WARNING::Fitting problem for feature 828 returning NA
#> 2026-06-17 13:45:18.65 INFO::Fitting model to feature number 829, ASV1133
#> 2026-06-17 13:45:18.65 INFO::Fitting model to feature number 830, ASV1134
#> 2026-06-17 13:45:18.65 WARNING::Fitting problem for feature 830 returning NA
#> 2026-06-17 13:45:18.66 INFO::Fitting model to feature number 831, ASV1135
#> 2026-06-17 13:45:18.66 INFO::Fitting model to feature number 832, ASV1137
#> 2026-06-17 13:45:18.66 INFO::Fitting model to feature number 833, ASV1138
#> 2026-06-17 13:45:18.66 INFO::Fitting model to feature number 834, ASV1139
#> 2026-06-17 13:45:18.67 INFO::Fitting model to feature number 835, ASV1141
#> 2026-06-17 13:45:18.67 INFO::Fitting model to feature number 836, ASV1143
#> 2026-06-17 13:45:18.67 INFO::Fitting model to feature number 837, ASV1144
#> 2026-06-17 13:45:18.68 INFO::Fitting model to feature number 838, ASV1146
#> 2026-06-17 13:45:18.68 INFO::Fitting model to feature number 839, ASV1147
#> 2026-06-17 13:45:18.68 INFO::Fitting model to feature number 840, ASV1148
#> 2026-06-17 13:45:18.69 INFO::Fitting model to feature number 841, ASV1150
#> 2026-06-17 13:45:18.69 INFO::Fitting model to feature number 842, ASV1151
#> 2026-06-17 13:45:18.69 INFO::Fitting model to feature number 843, ASV1152
#> 2026-06-17 13:45:18.69 INFO::Fitting model to feature number 844, ASV1154
#> 2026-06-17 13:45:18.70 INFO::Fitting model to feature number 845, ASV1155
#> 2026-06-17 13:45:18.70 INFO::Fitting model to feature number 846, ASV1156
#> 2026-06-17 13:45:18.70 INFO::Fitting model to feature number 847, ASV1158
#> 2026-06-17 13:45:18.71 INFO::Fitting model to feature number 848, ASV1159
#> 2026-06-17 13:45:18.71 WARNING::Fitting problem for feature 848 returning NA
#> 2026-06-17 13:45:18.71 INFO::Fitting model to feature number 849, ASV1160
#> 2026-06-17 13:45:18.71 INFO::Fitting model to feature number 850, ASV1161
#> 2026-06-17 13:45:18.72 INFO::Fitting model to feature number 851, ASV1162
#> 2026-06-17 13:45:18.72 INFO::Fitting model to feature number 852, ASV1163
#> 2026-06-17 13:45:18.73 INFO::Fitting model to feature number 853, ASV1164
#> 2026-06-17 13:45:18.73 INFO::Fitting model to feature number 854, ASV1165
#> 2026-06-17 13:45:18.73 WARNING::Fitting problem for feature 854 returning NA
#> 2026-06-17 13:45:18.73 INFO::Fitting model to feature number 855, ASV1167
#> 2026-06-17 13:45:18.74 INFO::Fitting model to feature number 856, ASV1168
#> 2026-06-17 13:45:18.74 INFO::Fitting model to feature number 857, ASV1169
#> 2026-06-17 13:45:18.74 WARNING::Fitting problem for feature 857 returning NA
#> 2026-06-17 13:45:18.74 INFO::Fitting model to feature number 858, ASV1171
#> 2026-06-17 13:45:18.75 INFO::Fitting model to feature number 859, ASV1172
#> 2026-06-17 13:45:18.75 INFO::Fitting model to feature number 860, ASV1173
#> 2026-06-17 13:45:18.76 INFO::Fitting model to feature number 861, ASV1175
#> 2026-06-17 13:45:18.76 INFO::Fitting model to feature number 862, ASV1176
#> 2026-06-17 13:45:18.76 INFO::Fitting model to feature number 863, ASV1177
#> 2026-06-17 13:45:18.77 INFO::Fitting model to feature number 864, ASV1179
#> 2026-06-17 13:45:18.77 WARNING::Fitting problem for feature 864 returning NA
#> 2026-06-17 13:45:18.77 INFO::Fitting model to feature number 865, ASV1180
#> 2026-06-17 13:45:18.78 INFO::Fitting model to feature number 866, ASV1182
#> 2026-06-17 13:45:18.78 INFO::Fitting model to feature number 867, ASV1184
#> 2026-06-17 13:45:18.78 INFO::Fitting model to feature number 868, ASV1185
#> 2026-06-17 13:45:18.78 INFO::Fitting model to feature number 869, ASV1186
#> 2026-06-17 13:45:18.79 INFO::Fitting model to feature number 870, ASV1187
#> 2026-06-17 13:45:18.79 INFO::Fitting model to feature number 871, ASV1189
#> 2026-06-17 13:45:18.79 INFO::Fitting model to feature number 872, ASV1190
#> 2026-06-17 13:45:18.80 INFO::Fitting model to feature number 873, ASV1192
#> 2026-06-17 13:45:18.80 INFO::Fitting model to feature number 874, ASV1193
#> 2026-06-17 13:45:18.81 INFO::Fitting model to feature number 875, ASV1194
#> 2026-06-17 13:45:18.81 INFO::Fitting model to feature number 876, ASV1195
#> 2026-06-17 13:45:18.82 INFO::Fitting model to feature number 877, ASV1198
#> 2026-06-17 13:45:18.82 INFO::Fitting model to feature number 878, ASV1199
#> 2026-06-17 13:45:18.82 INFO::Fitting model to feature number 879, ASV1200
#> 2026-06-17 13:45:18.82 INFO::Fitting model to feature number 880, ASV1203
#> 2026-06-17 13:45:18.83 INFO::Fitting model to feature number 881, ASV1204
#> 2026-06-17 13:45:18.86 INFO::Fitting model to feature number 882, ASV1205
#> 2026-06-17 13:45:18.86 INFO::Fitting model to feature number 883, ASV1206
#> 2026-06-17 13:45:18.87 INFO::Fitting model to feature number 884, ASV1208
#> 2026-06-17 13:45:18.87 INFO::Fitting model to feature number 885, ASV1209
#> 2026-06-17 13:45:18.87 INFO::Fitting model to feature number 886, ASV1210
#> 2026-06-17 13:45:18.87 INFO::Fitting model to feature number 887, ASV1211
#> 2026-06-17 13:45:18.88 INFO::Fitting model to feature number 888, ASV1212
#> 2026-06-17 13:45:18.88 INFO::Fitting model to feature number 889, ASV1213
#> 2026-06-17 13:45:18.89 INFO::Fitting model to feature number 890, ASV1214
#> 2026-06-17 13:45:18.89 INFO::Fitting model to feature number 891, ASV1216
#> 2026-06-17 13:45:18.89 INFO::Fitting model to feature number 892, ASV1217
#> 2026-06-17 13:45:18.90 INFO::Fitting model to feature number 893, ASV1218
#> 2026-06-17 13:45:18.90 INFO::Fitting model to feature number 894, ASV1219
#> 2026-06-17 13:45:18.90 INFO::Fitting model to feature number 895, ASV1221
#> 2026-06-17 13:45:18.91 WARNING::Fitting problem for feature 895 returning NA
#> 2026-06-17 13:45:18.91 INFO::Fitting model to feature number 896, ASV1223
#> 2026-06-17 13:45:18.91 INFO::Fitting model to feature number 897, ASV1224
#> 2026-06-17 13:45:18.91 INFO::Fitting model to feature number 898, ASV1225
#> 2026-06-17 13:45:18.92 INFO::Fitting model to feature number 899, ASV1227
#> 2026-06-17 13:45:18.92 INFO::Fitting model to feature number 900, ASV1228
#> 2026-06-17 13:45:18.92 INFO::Fitting model to feature number 901, ASV1229
#> 2026-06-17 13:45:18.92 WARNING::Fitting problem for feature 901 returning NA
#> 2026-06-17 13:45:18.92 INFO::Fitting model to feature number 902, ASV1230
#> 2026-06-17 13:45:18.93 INFO::Fitting model to feature number 903, ASV1231
#> 2026-06-17 13:45:18.93 INFO::Fitting model to feature number 904, ASV1232
#> 2026-06-17 13:45:18.94 INFO::Fitting model to feature number 905, ASV1233
#> 2026-06-17 13:45:18.94 INFO::Fitting model to feature number 906, ASV1234
#> 2026-06-17 13:45:18.94 WARNING::Fitting problem for feature 906 returning NA
#> 2026-06-17 13:45:18.94 INFO::Fitting model to feature number 907, ASV1236
#> 2026-06-17 13:45:18.95 INFO::Fitting model to feature number 908, ASV1238
#> 2026-06-17 13:45:18.95 INFO::Fitting model to feature number 909, ASV1239
#> 2026-06-17 13:45:18.96 INFO::Fitting model to feature number 910, ASV1241
#> 2026-06-17 13:45:18.96 INFO::Fitting model to feature number 911, ASV1242
#> 2026-06-17 13:45:18.96 INFO::Fitting model to feature number 912, ASV1243
#> 2026-06-17 13:45:18.97 INFO::Fitting model to feature number 913, ASV1245
#> 2026-06-17 13:45:18.97 INFO::Fitting model to feature number 914, ASV1246
#> 2026-06-17 13:45:18.97 INFO::Fitting model to feature number 915, ASV1247
#> 2026-06-17 13:45:18.97 INFO::Fitting model to feature number 916, ASV1251
#> 2026-06-17 13:45:18.98 INFO::Fitting model to feature number 917, ASV1252
#> 2026-06-17 13:45:18.98 INFO::Fitting model to feature number 918, ASV1253
#> 2026-06-17 13:45:18.99 INFO::Fitting model to feature number 919, ASV1254
#> 2026-06-17 13:45:18.99 INFO::Fitting model to feature number 920, ASV1257
#> 2026-06-17 13:45:18.99 INFO::Fitting model to feature number 921, ASV1258
#> 2026-06-17 13:45:19.00 INFO::Fitting model to feature number 922, ASV1259
#> 2026-06-17 13:45:19.00 INFO::Fitting model to feature number 923, ASV1260
#> 2026-06-17 13:45:19.00 INFO::Fitting model to feature number 924, ASV1261
#> 2026-06-17 13:45:19.00 INFO::Fitting model to feature number 925, ASV1262
#> 2026-06-17 13:45:19.01 INFO::Fitting model to feature number 926, ASV1263
#> 2026-06-17 13:45:19.01 INFO::Fitting model to feature number 927, ASV1264
#> 2026-06-17 13:45:19.01 INFO::Fitting model to feature number 928, ASV1265
#> 2026-06-17 13:45:19.02 INFO::Fitting model to feature number 929, ASV1267
#> 2026-06-17 13:45:19.02 INFO::Fitting model to feature number 930, ASV1268
#> 2026-06-17 13:45:19.02 INFO::Fitting model to feature number 931, ASV1269
#> 2026-06-17 13:45:19.02 INFO::Fitting model to feature number 932, ASV1270
#> 2026-06-17 13:45:19.03 INFO::Fitting model to feature number 933, ASV1271
#> 2026-06-17 13:45:19.03 INFO::Fitting model to feature number 934, ASV1272
#> 2026-06-17 13:45:19.03 INFO::Fitting model to feature number 935, ASV1273
#> 2026-06-17 13:45:19.03 INFO::Fitting model to feature number 936, ASV1274
#> 2026-06-17 13:45:19.04 INFO::Fitting model to feature number 937, ASV1275
#> 2026-06-17 13:45:19.04 INFO::Fitting model to feature number 938, ASV1276
#> 2026-06-17 13:45:19.04 INFO::Fitting model to feature number 939, ASV1278
#> 2026-06-17 13:45:19.05 INFO::Fitting model to feature number 940, ASV1279
#> 2026-06-17 13:45:19.05 INFO::Fitting model to feature number 941, ASV1282
#> 2026-06-17 13:45:19.05 INFO::Fitting model to feature number 942, ASV1283
#> 2026-06-17 13:45:19.05 INFO::Fitting model to feature number 943, ASV1284
#> 2026-06-17 13:45:19.06 INFO::Fitting model to feature number 944, ASV1285
#> 2026-06-17 13:45:19.06 INFO::Fitting model to feature number 945, ASV1286
#> 2026-06-17 13:45:19.06 INFO::Fitting model to feature number 946, ASV1287
#> 2026-06-17 13:45:19.07 INFO::Fitting model to feature number 947, ASV1288
#> 2026-06-17 13:45:19.07 INFO::Fitting model to feature number 948, ASV1289
#> 2026-06-17 13:45:19.07 INFO::Fitting model to feature number 949, ASV1290
#> 2026-06-17 13:45:19.08 INFO::Fitting model to feature number 950, ASV1293
#> 2026-06-17 13:45:19.08 INFO::Fitting model to feature number 951, ASV1294
#> 2026-06-17 13:45:19.08 INFO::Fitting model to feature number 952, ASV1296
#> 2026-06-17 13:45:19.09 INFO::Fitting model to feature number 953, ASV1297
#> 2026-06-17 13:45:19.09 INFO::Fitting model to feature number 954, ASV1300
#> 2026-06-17 13:45:19.09 INFO::Fitting model to feature number 955, ASV1301
#> 2026-06-17 13:45:19.10 INFO::Fitting model to feature number 956, ASV1302
#> 2026-06-17 13:45:19.10 INFO::Fitting model to feature number 957, ASV1303
#> 2026-06-17 13:45:19.11 INFO::Fitting model to feature number 958, ASV1304
#> 2026-06-17 13:45:19.11 INFO::Fitting model to feature number 959, ASV1305
#> 2026-06-17 13:45:19.11 INFO::Fitting model to feature number 960, ASV1307
#> 2026-06-17 13:45:19.11 INFO::Fitting model to feature number 961, ASV1310
#> 2026-06-17 13:45:19.12 INFO::Fitting model to feature number 962, ASV1311
#> 2026-06-17 13:45:19.12 INFO::Fitting model to feature number 963, ASV1312
#> 2026-06-17 13:45:19.13 INFO::Fitting model to feature number 964, ASV1313
#> 2026-06-17 13:45:19.13 INFO::Fitting model to feature number 965, ASV1314
#> 2026-06-17 13:45:19.13 INFO::Fitting model to feature number 966, ASV1315
#> 2026-06-17 13:45:19.14 INFO::Fitting model to feature number 967, ASV1316
#> 2026-06-17 13:45:19.14 INFO::Fitting model to feature number 968, ASV1317
#> 2026-06-17 13:45:19.14 INFO::Fitting model to feature number 969, ASV1319
#> 2026-06-17 13:45:19.15 INFO::Fitting model to feature number 970, ASV1320
#> 2026-06-17 13:45:19.15 INFO::Fitting model to feature number 971, ASV1321
#> 2026-06-17 13:45:19.16 INFO::Fitting model to feature number 972, ASV1323
#> 2026-06-17 13:45:19.16 INFO::Fitting model to feature number 973, ASV1326
#> 2026-06-17 13:45:19.16 INFO::Fitting model to feature number 974, ASV1327
#> 2026-06-17 13:45:19.17 INFO::Fitting model to feature number 975, ASV1328
#> 2026-06-17 13:45:19.17 INFO::Fitting model to feature number 976, ASV1330
#> 2026-06-17 13:45:19.17 INFO::Fitting model to feature number 977, ASV1332
#> 2026-06-17 13:45:19.18 INFO::Fitting model to feature number 978, ASV1334
#> 2026-06-17 13:45:19.18 INFO::Fitting model to feature number 979, ASV1335
#> 2026-06-17 13:45:19.18 INFO::Fitting model to feature number 980, ASV1336
#> 2026-06-17 13:45:19.19 INFO::Fitting model to feature number 981, ASV1337
#> 2026-06-17 13:45:19.19 INFO::Fitting model to feature number 982, ASV1338
#> 2026-06-17 13:45:19.19 INFO::Fitting model to feature number 983, ASV1340
#> 2026-06-17 13:45:19.20 INFO::Fitting model to feature number 984, ASV1341
#> 2026-06-17 13:45:19.20 INFO::Fitting model to feature number 985, ASV1342
#> 2026-06-17 13:45:19.20 INFO::Fitting model to feature number 986, ASV1345
#> 2026-06-17 13:45:19.21 INFO::Fitting model to feature number 987, ASV1350
#> 2026-06-17 13:45:19.21 INFO::Fitting model to feature number 988, ASV1351
#> 2026-06-17 13:45:19.22 INFO::Fitting model to feature number 989, ASV1352
#> 2026-06-17 13:45:19.22 INFO::Fitting model to feature number 990, ASV1353
#> 2026-06-17 13:45:19.22 INFO::Fitting model to feature number 991, ASV1355
#> 2026-06-17 13:45:19.23 INFO::Fitting model to feature number 992, ASV1356
#> 2026-06-17 13:45:19.23 INFO::Fitting model to feature number 993, ASV1359
#> 2026-06-17 13:45:19.23 INFO::Fitting model to feature number 994, ASV1360
#> 2026-06-17 13:45:19.24 INFO::Fitting model to feature number 995, ASV1361
#> 2026-06-17 13:45:19.24 INFO::Fitting model to feature number 996, ASV1362
#> 2026-06-17 13:45:19.25 INFO::Fitting model to feature number 997, ASV1363
#> 2026-06-17 13:45:19.25 INFO::Fitting model to feature number 998, ASV1365
#> 2026-06-17 13:45:19.25 INFO::Fitting model to feature number 999, ASV1366
#> 2026-06-17 13:45:19.26 INFO::Fitting model to feature number 1000, ASV1367
#> 2026-06-17 13:45:19.26 INFO::Fitting model to feature number 1001, ASV1368
#> 2026-06-17 13:45:19.26 INFO::Fitting model to feature number 1002, ASV1369
#> 2026-06-17 13:45:19.26 INFO::Fitting model to feature number 1003, ASV1370
#> 2026-06-17 13:45:19.27 INFO::Fitting model to feature number 1004, ASV1371
#> 2026-06-17 13:45:19.27 WARNING::Fitting problem for feature 1004 returning NA
#> 2026-06-17 13:45:19.27 INFO::Fitting model to feature number 1005, ASV1372
#> 2026-06-17 13:45:19.28 INFO::Fitting model to feature number 1006, ASV1373
#> 2026-06-17 13:45:19.28 INFO::Fitting model to feature number 1007, ASV1374
#> 2026-06-17 13:45:19.28 INFO::Fitting model to feature number 1008, ASV1375
#> 2026-06-17 13:45:19.28 INFO::Fitting model to feature number 1009, ASV1376
#> 2026-06-17 13:45:19.29 INFO::Fitting model to feature number 1010, ASV1378
#> 2026-06-17 13:45:19.29 INFO::Fitting model to feature number 1011, ASV1379
#> 2026-06-17 13:45:19.29 INFO::Fitting model to feature number 1012, ASV1380
#> 2026-06-17 13:45:19.29 INFO::Fitting model to feature number 1013, ASV1381
#> 2026-06-17 13:45:19.30 INFO::Fitting model to feature number 1014, ASV1384
#> 2026-06-17 13:45:19.30 INFO::Fitting model to feature number 1015, ASV1386
#> 2026-06-17 13:45:19.30 INFO::Fitting model to feature number 1016, ASV1387
#> 2026-06-17 13:45:19.30 INFO::Fitting model to feature number 1017, ASV1388
#> 2026-06-17 13:45:19.31 WARNING::Fitting problem for feature 1017 returning NA
#> 2026-06-17 13:45:19.31 INFO::Fitting model to feature number 1018, ASV1389
#> 2026-06-17 13:45:19.31 INFO::Fitting model to feature number 1019, ASV1390
#> 2026-06-17 13:45:19.31 INFO::Fitting model to feature number 1020, ASV1392
#> 2026-06-17 13:45:19.32 INFO::Fitting model to feature number 1021, ASV1393
#> 2026-06-17 13:45:19.32 INFO::Fitting model to feature number 1022, ASV1396
#> 2026-06-17 13:45:19.32 INFO::Fitting model to feature number 1023, ASV1397
#> 2026-06-17 13:45:19.32 INFO::Fitting model to feature number 1024, ASV1398
#> 2026-06-17 13:45:19.33 INFO::Fitting model to feature number 1025, ASV1399
#> 2026-06-17 13:45:19.33 WARNING::Fitting problem for feature 1025 returning NA
#> 2026-06-17 13:45:19.33 INFO::Fitting model to feature number 1026, ASV1400
#> 2026-06-17 13:45:19.33 INFO::Fitting model to feature number 1027, ASV1401
#> 2026-06-17 13:45:19.34 INFO::Fitting model to feature number 1028, ASV1403
#> 2026-06-17 13:45:19.34 INFO::Fitting model to feature number 1029, ASV1404
#> 2026-06-17 13:45:19.34 INFO::Fitting model to feature number 1030, ASV1406
#> 2026-06-17 13:45:19.35 INFO::Fitting model to feature number 1031, ASV1408
#> 2026-06-17 13:45:19.35 INFO::Fitting model to feature number 1032, ASV1409
#> 2026-06-17 13:45:19.35 INFO::Fitting model to feature number 1033, ASV1410
#> 2026-06-17 13:45:19.36 INFO::Fitting model to feature number 1034, ASV1412
#> 2026-06-17 13:45:19.36 INFO::Fitting model to feature number 1035, ASV1413
#> 2026-06-17 13:45:19.36 INFO::Fitting model to feature number 1036, ASV1414
#> 2026-06-17 13:45:19.36 INFO::Fitting model to feature number 1037, ASV1415
#> 2026-06-17 13:45:19.37 WARNING::Fitting problem for feature 1037 returning NA
#> 2026-06-17 13:45:19.37 INFO::Fitting model to feature number 1038, ASV1418
#> 2026-06-17 13:45:19.37 INFO::Fitting model to feature number 1039, ASV1419
#> 2026-06-17 13:45:19.37 INFO::Fitting model to feature number 1040, ASV1420
#> 2026-06-17 13:45:19.38 INFO::Fitting model to feature number 1041, ASV1421
#> 2026-06-17 13:45:19.38 INFO::Fitting model to feature number 1042, ASV1422
#> 2026-06-17 13:45:19.38 INFO::Fitting model to feature number 1043, ASV1423
#> 2026-06-17 13:45:19.38 WARNING::Fitting problem for feature 1043 returning NA
#> 2026-06-17 13:45:19.39 INFO::Fitting model to feature number 1044, ASV1424
#> 2026-06-17 13:45:19.39 INFO::Fitting model to feature number 1045, ASV1425
#> 2026-06-17 13:45:19.39 INFO::Fitting model to feature number 1046, ASV1427
#> 2026-06-17 13:45:19.40 INFO::Fitting model to feature number 1047, ASV1428
#> 2026-06-17 13:45:19.40 INFO::Fitting model to feature number 1048, ASV1430
#> 2026-06-17 13:45:19.40 INFO::Fitting model to feature number 1049, ASV1433
#> 2026-06-17 13:45:19.40 INFO::Fitting model to feature number 1050, ASV1434
#> 2026-06-17 13:45:19.41 INFO::Fitting model to feature number 1051, ASV1435
#> 2026-06-17 13:45:19.41 INFO::Fitting model to feature number 1052, ASV1436
#> 2026-06-17 13:45:19.41 WARNING::Fitting problem for feature 1052 returning NA
#> 2026-06-17 13:45:19.42 INFO::Fitting model to feature number 1053, ASV1437
#> 2026-06-17 13:45:19.42 INFO::Fitting model to feature number 1054, ASV1438
#> 2026-06-17 13:45:19.42 INFO::Fitting model to feature number 1055, ASV1442
#> 2026-06-17 13:45:19.42 INFO::Fitting model to feature number 1056, ASV1443
#> 2026-06-17 13:45:19.42 INFO::Fitting model to feature number 1057, ASV1449
#> 2026-06-17 13:45:19.43 INFO::Fitting model to feature number 1058, ASV1450
#> 2026-06-17 13:45:19.43 INFO::Fitting model to feature number 1059, ASV1452
#> 2026-06-17 13:45:19.43 INFO::Fitting model to feature number 1060, ASV1454
#> 2026-06-17 13:45:19.43 INFO::Fitting model to feature number 1061, ASV1455
#> 2026-06-17 13:45:19.44 INFO::Fitting model to feature number 1062, ASV1458
#> 2026-06-17 13:45:19.44 INFO::Fitting model to feature number 1063, ASV1459
#> 2026-06-17 13:45:19.44 INFO::Fitting model to feature number 1064, ASV1460
#> 2026-06-17 13:45:19.45 INFO::Fitting model to feature number 1065, ASV1461
#> 2026-06-17 13:45:19.45 INFO::Fitting model to feature number 1066, ASV1462
#> 2026-06-17 13:45:19.45 INFO::Fitting model to feature number 1067, ASV1463
#> 2026-06-17 13:45:19.46 INFO::Fitting model to feature number 1068, ASV1466
#> 2026-06-17 13:45:19.46 INFO::Fitting model to feature number 1069, ASV1467
#> 2026-06-17 13:45:19.47 INFO::Fitting model to feature number 1070, ASV1468
#> 2026-06-17 13:45:19.47 INFO::Fitting model to feature number 1071, ASV1469
#> 2026-06-17 13:45:19.47 WARNING::Fitting problem for feature 1071 returning NA
#> 2026-06-17 13:45:19.47 INFO::Fitting model to feature number 1072, ASV1472
#> 2026-06-17 13:45:19.47 INFO::Fitting model to feature number 1073, ASV1477
#> 2026-06-17 13:45:19.48 INFO::Fitting model to feature number 1074, ASV1478
#> 2026-06-17 13:45:19.48 INFO::Fitting model to feature number 1075, ASV1479
#> 2026-06-17 13:45:19.49 INFO::Fitting model to feature number 1076, ASV1483
#> 2026-06-17 13:45:19.49 INFO::Fitting model to feature number 1077, ASV1484
#> 2026-06-17 13:45:19.49 INFO::Fitting model to feature number 1078, ASV1486
#> 2026-06-17 13:45:19.50 INFO::Fitting model to feature number 1079, ASV1487
#> 2026-06-17 13:45:19.50 INFO::Fitting model to feature number 1080, ASV1488
#> 2026-06-17 13:45:19.50 INFO::Fitting model to feature number 1081, ASV1490
#> 2026-06-17 13:45:19.50 INFO::Fitting model to feature number 1082, ASV1492
#> 2026-06-17 13:45:19.51 INFO::Fitting model to feature number 1083, ASV1493
#> 2026-06-17 13:45:19.51 INFO::Fitting model to feature number 1084, ASV1494
#> 2026-06-17 13:45:19.52 INFO::Fitting model to feature number 1085, ASV1495
#> 2026-06-17 13:45:19.52 INFO::Fitting model to feature number 1086, ASV1496
#> 2026-06-17 13:45:19.52 INFO::Fitting model to feature number 1087, ASV1497
#> 2026-06-17 13:45:19.53 INFO::Fitting model to feature number 1088, ASV1498
#> 2026-06-17 13:45:19.53 INFO::Fitting model to feature number 1089, ASV1500
#> 2026-06-17 13:45:19.53 INFO::Fitting model to feature number 1090, ASV1501
#> 2026-06-17 13:45:19.54 INFO::Fitting model to feature number 1091, ASV1502
#> 2026-06-17 13:45:19.54 INFO::Fitting model to feature number 1092, ASV1503
#> 2026-06-17 13:45:19.54 INFO::Fitting model to feature number 1093, ASV1504
#> 2026-06-17 13:45:19.54 WARNING::Fitting problem for feature 1093 returning NA
#> 2026-06-17 13:45:19.55 INFO::Fitting model to feature number 1094, ASV1507
#> 2026-06-17 13:45:19.55 INFO::Fitting model to feature number 1095, ASV1509
#> 2026-06-17 13:45:19.55 INFO::Fitting model to feature number 1096, ASV1510
#> 2026-06-17 13:45:19.55 INFO::Fitting model to feature number 1097, ASV1513
#> 2026-06-17 13:45:19.56 INFO::Fitting model to feature number 1098, ASV1515
#> 2026-06-17 13:45:19.56 INFO::Fitting model to feature number 1099, ASV1516
#> 2026-06-17 13:45:19.56 INFO::Fitting model to feature number 1100, ASV1517
#> 2026-06-17 13:45:19.57 INFO::Fitting model to feature number 1101, ASV1518
#> 2026-06-17 13:45:19.57 INFO::Fitting model to feature number 1102, ASV1521
#> 2026-06-17 13:45:19.57 INFO::Fitting model to feature number 1103, ASV1523
#> 2026-06-17 13:45:19.57 INFO::Fitting model to feature number 1104, ASV1525
#> 2026-06-17 13:45:19.58 INFO::Fitting model to feature number 1105, ASV1526
#> 2026-06-17 13:45:19.58 INFO::Fitting model to feature number 1106, ASV1527
#> 2026-06-17 13:45:19.58 INFO::Fitting model to feature number 1107, ASV1528
#> 2026-06-17 13:45:19.58 INFO::Fitting model to feature number 1108, ASV1531
#> 2026-06-17 13:45:19.59 INFO::Fitting model to feature number 1109, ASV1532
#> 2026-06-17 13:45:19.59 INFO::Fitting model to feature number 1110, ASV1533
#> 2026-06-17 13:45:19.59 INFO::Fitting model to feature number 1111, ASV1534
#> 2026-06-17 13:45:19.60 WARNING::Fitting problem for feature 1111 returning NA
#> 2026-06-17 13:45:19.60 INFO::Fitting model to feature number 1112, ASV1536
#> 2026-06-17 13:45:19.60 INFO::Fitting model to feature number 1113, ASV1537
#> 2026-06-17 13:45:19.60 INFO::Fitting model to feature number 1114, ASV1542
#> 2026-06-17 13:45:19.61 INFO::Fitting model to feature number 1115, ASV1543
#> 2026-06-17 13:45:19.61 INFO::Fitting model to feature number 1116, ASV1544
#> 2026-06-17 13:45:19.61 INFO::Fitting model to feature number 1117, ASV1545
#> 2026-06-17 13:45:19.61 INFO::Fitting model to feature number 1118, ASV1548
#> 2026-06-17 13:45:19.62 INFO::Fitting model to feature number 1119, ASV1549
#> 2026-06-17 13:45:19.62 INFO::Fitting model to feature number 1120, ASV1550
#> 2026-06-17 13:45:19.62 INFO::Fitting model to feature number 1121, ASV1551
#> 2026-06-17 13:45:19.62 INFO::Fitting model to feature number 1122, ASV1552
#> 2026-06-17 13:45:19.63 INFO::Fitting model to feature number 1123, ASV1553
#> 2026-06-17 13:45:19.63 INFO::Fitting model to feature number 1124, ASV1554
#> 2026-06-17 13:45:19.63 INFO::Fitting model to feature number 1125, ASV1557
#> 2026-06-17 13:45:19.63 INFO::Fitting model to feature number 1126, ASV1558
#> 2026-06-17 13:45:19.64 INFO::Fitting model to feature number 1127, ASV1560
#> 2026-06-17 13:45:19.64 INFO::Fitting model to feature number 1128, ASV1562
#> 2026-06-17 13:45:19.64 INFO::Fitting model to feature number 1129, ASV1564
#> 2026-06-17 13:45:19.65 INFO::Fitting model to feature number 1130, ASV1565
#> 2026-06-17 13:45:19.65 INFO::Fitting model to feature number 1131, ASV1568
#> 2026-06-17 13:45:19.65 INFO::Fitting model to feature number 1132, ASV1569
#> 2026-06-17 13:45:19.65 INFO::Fitting model to feature number 1133, ASV1570
#> 2026-06-17 13:45:19.65 INFO::Fitting model to feature number 1134, ASV1571
#> 2026-06-17 13:45:19.66 INFO::Fitting model to feature number 1135, ASV1573
#> 2026-06-17 13:45:19.66 INFO::Fitting model to feature number 1136, ASV1574
#> 2026-06-17 13:45:19.66 INFO::Fitting model to feature number 1137, ASV1575
#> 2026-06-17 13:45:19.66 INFO::Fitting model to feature number 1138, ASV1576
#> 2026-06-17 13:45:19.66 INFO::Fitting model to feature number 1139, ASV1577
#> 2026-06-17 13:45:19.67 INFO::Fitting model to feature number 1140, ASV1578
#> 2026-06-17 13:45:19.67 INFO::Fitting model to feature number 1141, ASV1580
#> 2026-06-17 13:45:19.67 INFO::Fitting model to feature number 1142, ASV1581
#> 2026-06-17 13:45:19.68 INFO::Fitting model to feature number 1143, ASV1584
#> 2026-06-17 13:45:19.68 INFO::Fitting model to feature number 1144, ASV1585
#> 2026-06-17 13:45:19.68 INFO::Fitting model to feature number 1145, ASV1586
#> 2026-06-17 13:45:19.68 INFO::Fitting model to feature number 1146, ASV1588
#> 2026-06-17 13:45:19.68 INFO::Fitting model to feature number 1147, ASV1589
#> 2026-06-17 13:45:19.69 INFO::Fitting model to feature number 1148, ASV1592
#> 2026-06-17 13:45:19.69 INFO::Fitting model to feature number 1149, ASV1593
#> 2026-06-17 13:45:19.69 INFO::Fitting model to feature number 1150, ASV1594
#> 2026-06-17 13:45:19.69 INFO::Fitting model to feature number 1151, ASV1601
#> 2026-06-17 13:45:19.70 INFO::Fitting model to feature number 1152, ASV1602
#> 2026-06-17 13:45:19.70 INFO::Fitting model to feature number 1153, ASV1603
#> 2026-06-17 13:45:19.70 INFO::Fitting model to feature number 1154, ASV1604
#> 2026-06-17 13:45:19.71 WARNING::Fitting problem for feature 1154 returning NA
#> 2026-06-17 13:45:19.71 INFO::Fitting model to feature number 1155, ASV1605
#> 2026-06-17 13:45:19.71 INFO::Fitting model to feature number 1156, ASV1606
#> 2026-06-17 13:45:19.72 INFO::Fitting model to feature number 1157, ASV1607
#> 2026-06-17 13:45:19.72 INFO::Fitting model to feature number 1158, ASV1609
#> 2026-06-17 13:45:19.72 WARNING::Fitting problem for feature 1158 returning NA
#> 2026-06-17 13:45:19.72 INFO::Fitting model to feature number 1159, ASV1612
#> 2026-06-17 13:45:19.73 INFO::Fitting model to feature number 1160, ASV1613
#> 2026-06-17 13:45:19.73 INFO::Fitting model to feature number 1161, ASV1614
#> 2026-06-17 13:45:19.73 INFO::Fitting model to feature number 1162, ASV1615
#> 2026-06-17 13:45:19.73 INFO::Fitting model to feature number 1163, ASV1616
#> 2026-06-17 13:45:19.73 INFO::Fitting model to feature number 1164, ASV1620
#> 2026-06-17 13:45:19.74 INFO::Fitting model to feature number 1165, ASV1621
#> 2026-06-17 13:45:19.74 INFO::Fitting model to feature number 1166, ASV1622
#> 2026-06-17 13:45:19.74 WARNING::Fitting problem for feature 1166 returning NA
#> 2026-06-17 13:45:19.74 INFO::Fitting model to feature number 1167, ASV1624
#> 2026-06-17 13:45:19.75 INFO::Fitting model to feature number 1168, ASV1625
#> 2026-06-17 13:45:19.75 INFO::Fitting model to feature number 1169, ASV1628
#> 2026-06-17 13:45:19.75 INFO::Fitting model to feature number 1170, ASV1629
#> 2026-06-17 13:45:19.76 INFO::Fitting model to feature number 1171, ASV1630
#> 2026-06-17 13:45:19.76 INFO::Fitting model to feature number 1172, ASV1631
#> 2026-06-17 13:45:19.76 INFO::Fitting model to feature number 1173, ASV1632
#> 2026-06-17 13:45:19.76 INFO::Fitting model to feature number 1174, ASV1633
#> 2026-06-17 13:45:19.77 INFO::Fitting model to feature number 1175, ASV1634
#> 2026-06-17 13:45:19.77 INFO::Fitting model to feature number 1176, ASV1635
#> 2026-06-17 13:45:19.77 INFO::Fitting model to feature number 1177, ASV1636
#> 2026-06-17 13:45:19.77 INFO::Fitting model to feature number 1178, ASV1639
#> 2026-06-17 13:45:19.78 INFO::Fitting model to feature number 1179, ASV1641
#> 2026-06-17 13:45:19.78 WARNING::Fitting problem for feature 1179 returning NA
#> 2026-06-17 13:45:19.78 INFO::Fitting model to feature number 1180, ASV1644
#> 2026-06-17 13:45:19.78 INFO::Fitting model to feature number 1181, ASV1645
#> 2026-06-17 13:45:19.78 INFO::Fitting model to feature number 1182, ASV1649
#> 2026-06-17 13:45:19.79 INFO::Fitting model to feature number 1183, ASV1650
#> 2026-06-17 13:45:19.79 INFO::Fitting model to feature number 1184, ASV1651
#> 2026-06-17 13:45:19.79 WARNING::Fitting problem for feature 1184 returning NA
#> 2026-06-17 13:45:19.80 INFO::Fitting model to feature number 1185, ASV1653
#> 2026-06-17 13:45:19.80 INFO::Fitting model to feature number 1186, ASV1654
#> 2026-06-17 13:45:19.80 INFO::Fitting model to feature number 1187, ASV1656
#> 2026-06-17 13:45:19.80 INFO::Fitting model to feature number 1188, ASV1657
#> 2026-06-17 13:45:19.81 INFO::Fitting model to feature number 1189, ASV1658
#> 2026-06-17 13:45:19.81 INFO::Fitting model to feature number 1190, ASV1661
#> 2026-06-17 13:45:19.81 INFO::Fitting model to feature number 1191, ASV1663
#> 2026-06-17 13:45:19.82 WARNING::Fitting problem for feature 1191 returning NA
#> 2026-06-17 13:45:19.82 INFO::Fitting model to feature number 1192, ASV1664
#> 2026-06-17 13:45:19.82 INFO::Fitting model to feature number 1193, ASV1666
#> 2026-06-17 13:45:19.82 INFO::Fitting model to feature number 1194, ASV1667
#> 2026-06-17 13:45:19.82 INFO::Fitting model to feature number 1195, ASV1670
#> 2026-06-17 13:45:19.83 INFO::Fitting model to feature number 1196, ASV1671
#> 2026-06-17 13:45:19.83 INFO::Fitting model to feature number 1197, ASV1673
#> 2026-06-17 13:45:19.83 INFO::Fitting model to feature number 1198, ASV1674
#> 2026-06-17 13:45:19.83 INFO::Fitting model to feature number 1199, ASV1677
#> 2026-06-17 13:45:19.84 INFO::Fitting model to feature number 1200, ASV1681
#> 2026-06-17 13:45:19.84 INFO::Fitting model to feature number 1201, ASV1683
#> 2026-06-17 13:45:19.84 INFO::Fitting model to feature number 1202, ASV1684
#> 2026-06-17 13:45:19.85 INFO::Fitting model to feature number 1203, ASV1687
#> 2026-06-17 13:45:19.85 WARNING::Fitting problem for feature 1203 returning NA
#> 2026-06-17 13:45:19.85 INFO::Fitting model to feature number 1204, ASV1689
#> 2026-06-17 13:45:19.85 INFO::Fitting model to feature number 1205, ASV1690
#> 2026-06-17 13:45:19.86 INFO::Fitting model to feature number 1206, ASV1691
#> 2026-06-17 13:45:19.86 INFO::Fitting model to feature number 1207, ASV1694
#> 2026-06-17 13:45:19.86 INFO::Fitting model to feature number 1208, ASV1697
#> 2026-06-17 13:45:19.87 INFO::Fitting model to feature number 1209, ASV1704
#> 2026-06-17 13:45:19.87 INFO::Fitting model to feature number 1210, ASV1706
#> 2026-06-17 13:45:19.87 INFO::Fitting model to feature number 1211, ASV1707
#> 2026-06-17 13:45:19.87 WARNING::Fitting problem for feature 1211 returning NA
#> 2026-06-17 13:45:19.88 INFO::Fitting model to feature number 1212, ASV1712
#> 2026-06-17 13:45:19.88 WARNING::Fitting problem for feature 1212 returning NA
#> 2026-06-17 13:45:20.04 WARNING::Deleting existing residuals file: res_maaslin3/fits/residuals_linear.rds
#> 2026-06-17 13:45:20.04 INFO::Writing residuals to file res_maaslin3/fits/residuals_linear.rds
#> 2026-06-17 13:45:20.05 WARNING::Deleting existing fitted file: res_maaslin3/fits/fitted_linear.rds
#> 2026-06-17 13:45:20.05 INFO::Writing fitted values to file res_maaslin3/fits/fitted_linear.rds
#> 2026-06-17 13:45:20.07 WARNING::Deleting existing residuals file: res_maaslin3/fits/residuals_logistic.rds
#> 2026-06-17 13:45:20.07 INFO::Writing residuals to file res_maaslin3/fits/residuals_logistic.rds
#> 2026-06-17 13:45:20.19 WARNING::Deleting existing fitted file: res_maaslin3/fits/fitted_logistic.rds
#> 2026-06-17 13:45:20.19 INFO::Writing fitted values to file res_maaslin3/fits/fitted_logistic.rds
#> 2026-06-17 13:45:20.31 INFO::Writing all the results to file (ordered 
#>             by increasing individual q-values): res_maaslin3/all_results.tsv
#> 2026-06-17 13:45:20.35 INFO::Writing the significant results without errors (those which have joint q-values less than or equal to the threshold of 0.100000 ) to file (ordered by increasing individual q-values): res_maaslin3/significant_results.tsv
#> 2026-06-17 13:45:20.37 INFO::Writing summary plot of significant
#>                         results to file: res_maaslin3/figures/summary_plot.pdf
#> 2026-06-17 13:45:22.25 INFO::Writing association plots (one for each significant association) to output folder: res_maaslin3/figures
#> 2026-06-17 13:45:22.25 INFO::Plotting associations from most to least significant, grouped by metadata
#> 2026-06-17 13:45:22.26 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV13
#> 2026-06-17 13:45:22.58 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV187
#> 2026-06-17 13:45:22.89 INFO::Creating scatter plot for continuous 
#>                         data (linear), nb_seq vs ASV758
#> 2026-06-17 13:45:23.18 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV166
#> 2026-06-17 13:45:23.52 INFO::Creating scatter plot for continuous 
#>                         data (linear), nb_seq vs ASV1128
#> 2026-06-17 13:45:23.81 INFO::Creating scatter plot for continuous 
#>                         data (linear), nb_seq vs ASV481
#> 2026-06-17 13:45:24.10 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV509
#> 2026-06-17 13:45:24.41 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV378
#> 2026-06-17 13:45:24.73 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV340
#> 2026-06-17 13:45:25.05 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV56
#> 2026-06-17 13:45:25.38 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV71
#> 2026-06-17 13:45:25.70 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV217
#> 2026-06-17 13:45:26.02 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV2
#> 2026-06-17 13:45:26.33 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV445
#> 2026-06-17 13:45:26.65 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV438
#> 2026-06-17 13:45:26.99 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV263
#> 2026-06-17 13:45:27.31 INFO::Creating scatter plot for continuous 
#>                         data (linear), nb_seq vs ASV724
#> 2026-06-17 13:45:27.59 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV494
#> 2026-06-17 13:45:27.91 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV519
#> 2026-06-17 13:45:28.23 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV581
#> 2026-06-17 13:45:28.57 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV734
#> 2026-06-17 13:45:28.89 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV8
#> 2026-06-17 13:45:29.20 INFO::Creating scatter plot for continuous 
#>                         data (linear), nb_seq vs ASV1047
#> 2026-06-17 13:45:29.49 INFO::Creating scatter plot for continuous 
#>                         data (linear), nb_seq vs ASV1198
#> 2026-06-17 13:45:29.78 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV295
#> 2026-06-17 13:45:30.10 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV512
#> 2026-06-17 13:45:30.44 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV216
#> 2026-06-17 13:45:30.75 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV234
#> 2026-06-17 13:45:31.06 INFO::Creating boxplot for continuous data (logistic), nb_seq vs ASV727
#> 2026-06-17 13:45:31.38 INFO::Creating scatter plot for continuous 
#>                         data (linear), nb_seq vs ASV1427

gg_maaslin3_plot(res, type = "volcano")
#> Warning: Removed 1020 rows containing missing values or values outside the scale range
#> (`geom_point()`).



# Set to 0 the sequences numbers of 80% of taxa in "High" samples
 data_fungi_mini_with_less_otu_in_High =  multiply_counts_pq(data_fungi_mini,
      fact = "Height", prop=0.8,
      conditions = "High",
      multipliers = 0)
#> Modified 36 taxa in 28 matched samples

ggbetween_pq(data_fungi_mini_with_less_otu_in_High,
  "Height",
  one_plot=T)
#> Taxa are now in columns.
#> Warning: The mean number of sequences per samples vary across modalities of the variable 'Height' You should use rarefy_by_sample = TRUE or try hill_pq() with correction_for_sample_size = TRUE
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
```
