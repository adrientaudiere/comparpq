# S7 class for comparing phyloseq objects

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A class to store and compare multiple phyloseq objects. It contains:

- A list of phyloseq objects

- A summary table with computed characteristics for each phyloseq

- A list of comparison characteristics between phyloseq objects

An S7 class to store and compare multiple phyloseq objects.

## Usage

``` r
list_phyloseq(
  physeq_list,
  names = NULL,
  same_primer_seq_tech = TRUE,
  same_bioinfo_pipeline = TRUE
)
```

## Arguments

- physeq_list:

  (required) A named list of phyloseq objects.

- names:

  (character vector, default NULL) Optional names for the phyloseq
  objects. If NULL and the list is unnamed, names are generated
  automatically.

- same_primer_seq_tech:

  (logical, default TRUE) Whether the same primer and sequencing
  technology was used across all phyloseq objects. Set to FALSE when
  comparing different primers (e.g., ITS1 vs ITS2) or technologies
  (e.g., Illumina vs PacBio).

- same_bioinfo_pipeline:

  (logical, default TRUE) Whether the same bioinformatics pipeline was
  used across all phyloseq objects. Set to FALSE when comparing
  different clustering methods, taxonomic databases, or analysis
  parameters.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The combination of `same_primer_seq_tech` and `same_bioinfo_pipeline`
along with detected sample overlap determines the comparison type:

|  |  |  |  |
|----|----|----|----|
| same_samples | same_primer_seq_tech | same_bioinfo_pipeline | Type |
| TRUE | TRUE | TRUE | REPRODUCIBILITY |
| TRUE | TRUE | FALSE | ROBUSTNESS |
| TRUE | FALSE | \- | REPLICABILITY |
| nested | \- | \- | NESTED_ROBUSTNESS |
| FALSE | \- | \- | EXPLORATION or SEPARATE_ANALYSIS |

## Slots

- `phyloseq_list`:

  A named list of phyloseq objects

- `summary_table`:

  A tibble summarizing each phyloseq object

- `comparison`:

  A list of characteristics comparing the phyloseq objects

## Types of Comparison

The `list_phyloseq` class determines the type of comparison based on:

- Detected characteristics (sample overlap, nested samples, shared
  modalities)

- User-provided parameters (`same_primer_seq_tech` and
  `same_bioinfo_pipeline`)

There are six main types of comparisons:

- REPRODUCIBILITY:

  Same pipeline (`same_bioinfo_pipeline = TRUE`), same primer/technology
  (`same_primer_seq_tech = TRUE`), same samples. Used to test
  **reproducibility** of results when running the exact same analysis
  multiple times.

- ROBUSTNESS:

  Different pipeline (`same_bioinfo_pipeline = FALSE`, e.g., different
  clustering method, different assignment database) but same
  primer/technology (`same_primer_seq_tech = TRUE`), same samples. Used
  to test **robustness** of conclusions to methodological choices.

- NESTED_ROBUSTNESS:

  One phyloseq object is derived from another, with samples being a
  subset (e.g., rarefied version created with
  [`rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html)).
  Used to test **robustness** to data processing choices like
  rarefaction, filtering, or subsetting. Comparisons should focus on the
  common (nested) samples.

- REPLICABILITY:

  Different primer and/or technology (`same_primer_seq_tech = FALSE`,
  e.g., ITS1 vs ITS2, Illumina vs PacBio), same samples. Used to test
  **replicability** across taxonomic groups or sequencing technologies.

- EXPLORATION:

  Different samples but with shared modalities. Useful to **explore**
  differences among groups of samples. Note: consider merging samples
  into one phyloseq object for some analyses instead.

- SEPARATE_ANALYSIS:

  Different samples with no shared modalities. **Separate analysis** of
  each phyloseq object is recommended as direct comparison may not be
  meaningful.

## Examples

``` r
# REPRODUCIBILITY: Same samples, same pipeline (default)
lpq_repro <- list_phyloseq(list(run1 = data_fungi, run2 = data_fungi))
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: REPRODUCIBILITY
#> ℹ 185 common samples, 1420 common taxa
#> ✔ list_phyloseq created (REPRODUCIBILITY)

# ROBUSTNESS: Same samples, different pipeline
lpq_robust <- list_phyloseq(
  list(method_A = data_fungi, method_B = data_fungi),
  same_bioinfo_pipeline = FALSE
)
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: ROBUSTNESS
#> ℹ 185 common samples, 1420 common taxa
#> ✔ list_phyloseq created (ROBUSTNESS)

# REPLICABILITY: Same samples, different primer/technology
lpq_replic <- list_phyloseq(
  list(ITS1 = data_fungi, ITS2 = data_fungi),
  same_primer_seq_tech = FALSE
)
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: REPLICABILITY
#> ℹ 185 common samples, 1420 common taxa
#> ✔ list_phyloseq created (REPLICABILITY)
```
