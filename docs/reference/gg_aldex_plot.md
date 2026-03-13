# ggplot2 version of ALDEx2 diagnostic plots

Creates ggplot2 versions of the four ALDEx2 diagnostic plot types: MW
(Bland-Altman style), MA, volcano, and volcano.var. Supports faceting by
a `name` column, making it suitable for output from
[`aldex_lpq()`](https://adrientaudiere.github.io/comparpq/reference/aldex_lpq.md).

## Usage

``` r
gg_aldex_plot(
  x,
  type = c("MW", "MA", "volcano", "volcano.var"),
  test = c("welch", "wilcox", "effect", "both"),
  cutoff.pval = 0.05,
  cutoff.effect = 1,
  rare = 0,
  all.col = "grey30",
  called.col = "red",
  rare.col = "black",
  point.size = 1.5,
  point.alpha = 0.4
)
```

## Arguments

- x:

  (data.frame or tibble, required) ALDEx2 results, typically from
  [`ALDEx2::aldex()`](https://rdrr.io/pkg/ALDEx2/man/aldex.html) or
  [`aldex_lpq()`](https://adrientaudiere.github.io/comparpq/reference/aldex_lpq.md).
  Must contain columns `diff.btw`, `diff.win`, `rab.all`, and `effect`.
  When `test` is `"welch"`, column `we.eBH` is required; when
  `"wilcox"`, column `wi.eBH` is required.

- type:

  (character, default "MW") Plot type. One of `"MW"` (dispersion vs
  difference), `"MA"` (abundance vs difference), `"volcano"` (difference
  vs -log10 q-value), or `"volcano.var"` (dispersion vs -log10 q-value).

- test:

  (character, default "welch") Statistical test used for significance
  calling. One of `"welch"`, `"wilcox"`, `"effect"`, or `"both"`
  (effect + welch).

- cutoff.pval:

  (numeric, default 0.05) q-value threshold for significance.

- cutoff.effect:

  (numeric, default 1) Effect size threshold (must be \>= 0.5 when
  `test` is `"effect"` or `"both"`).

- rare:

  (numeric, default 0) Abundance threshold below which taxa are marked
  as rare (only used for MW and MA plots).

- all.col:

  (character, default "grey30") Color for non-significant points.

- called.col:

  (character, default "red") Color for significant points.

- rare.col:

  (character, default "black") Color for rare taxa points.

- point.size:

  (numeric, default 1.5) Size of points.

- point.alpha:

  (numeric, default 0.4) Alpha transparency of non-significant points.

## Value

A ggplot2 object. If `x` contains a `name` column (e.g., from
[`aldex_lpq()`](https://adrientaudiere.github.io/comparpq/reference/aldex_lpq.md)),
the plot is faceted by `name`.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function reimplements
[`ALDEx2::aldex.plot()`](https://rdrr.io/pkg/ALDEx2/man/aldex.plot.html)
using ggplot2, providing a more customizable and composable output. The
four plot types correspond to the original ALDEx2 types:

- MW:

  Bland-Altman style: within-condition dispersion (x) vs
  between-condition difference (y), with +/-1 effect size lines.

- MA:

  MA plot: median log2 relative abundance (x) vs median log2 difference
  (y).

- volcano:

  Volcano plot: median log2 difference (x) vs -log10 q-value (y).

- volcano.var:

  Variance volcano: median log2 dispersion (x) vs -log10 q-value (y).

## See also

[`aldex_lpq()`](https://adrientaudiere.github.io/comparpq/reference/aldex_lpq.md),
[`ALDEx2::aldex.plot()`](https://rdrr.io/pkg/ALDEx2/man/aldex.plot.html)

## Examples

``` r
data_fungi_high <- multiply_counts_pq(data_fungi, "Height", "High",
  4,
  prop_taxa = 0.1, seed = 42
)
#> Modified 142 taxa in 41 matched samples

aldex_pq(data_fungi_high,
  bifactor = "Height",
  modalities = c("Low", "High")
) |>
  gg_aldex_plot()
#> Taxa are now in rows.
#> aldex.clr: generating Monte-Carlo instances and clr values
#> conditions vector supplied
#> operating in serial mode
#> aldex.scaleSim: adjusting samples to reflect scale uncertainty.
#> aldex.ttest: doing t-test
#> aldex.effect: calculating effect sizes


lpq <- list_phyloseq(
  list(
    fungi = data_fungi,
    fungi_height = data_fungi_high
  ),
  same_bioinfo_pipeline = FALSE
)
#> ℹ Building summary table for 2 phyloseq objects...
#> ℹ Computing comparison characteristics...
#> ℹ Checking sample and taxa overlap...
#> ℹ Detected comparison type: ROBUSTNESS
#> ℹ 185 common samples, 1420 common taxa
#> ✔ list_phyloseq created (ROBUSTNESS)
# From aldex_lpq (faceted by name)
lpq_res <- aldex_lpq(lpq,
  bifactor = "Height",
  modalities = c("Low", "High")
)
#> Running ALDEx2 on 2 phyloseq objects
#> Bifactor: Height
#> → Processing: fungi
#> Taxa are now in rows.
#> aldex.clr: generating Monte-Carlo instances and clr values
#> conditions vector supplied
#> operating in serial mode
#> aldex.scaleSim: adjusting samples to reflect scale uncertainty.
#> aldex.ttest: doing t-test
#> aldex.effect: calculating effect sizes
#> → Processing: fungi_height
#> Taxa are now in rows.
#> aldex.clr: generating Monte-Carlo instances and clr values
#> conditions vector supplied
#> operating in serial mode
#> aldex.scaleSim: adjusting samples to reflect scale uncertainty.
#> aldex.ttest: doing t-test
#> aldex.effect: calculating effect sizes
gg_aldex_plot(lpq_res, type = "MA")

gg_aldex_plot(lpq_res, type = "MW", test = "wilcox")


gg_aldex_plot(lpq_res, type = "volcano")

gg_aldex_plot(lpq_res, type = "volcano.var")



gingival_pq <-
  MicrobiomeBenchmarkData::getBenchmarkData("HMP_2012_16S_gingival_V35_subset",
    dryrun = FALSE
  )[[1]] |>
  mia::convertToPhyloseq()
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "atomicVector" in cache; using the first, from namespace 'Matrix'
#> Also defined by ‘Rmpfr’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Finished HMP_2012_16S_gingival_V35_subset.
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’

aldex_res <- aldex_pq(gingival_pq,
  bifactor = "body_subsite",
  modalities = c("supragingival_plaque", "subgingival_plaque")
)
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘tidytree’
#> aldex.clr: generating Monte-Carlo instances and clr values
#> conditions vector supplied
#> operating in serial mode
#> aldex.scaleSim: adjusting samples to reflect scale uncertainty.
#> aldex.ttest: doing t-test
#> aldex.effect: calculating effect sizes

gg_aldex_plot(aldex_res, type = "volcano")

```
