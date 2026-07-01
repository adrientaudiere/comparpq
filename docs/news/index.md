# Changelog

## comparpq 0.3.0 (Development version)

## comparpq 0.2.1

- [`refseq_comp_lpq()`](https://adrientaudiere.github.io/comparpq/reference/refseq_comp_lpq.md)
  new function to compare `@refseq` sequences across all phyloseq
  objects in a `list_phyloseq` using k-mer Jaccard similarity and
  union-find connected components. Returns per-threshold Venn diagrams
  and shared-cluster counts. No igraph dependency.

- [`find_primers_pq()`](https://adrientaudiere.github.io/comparpq/reference/find_primers_pq.md)
  new function to detect taxa whose reference sequences match primer
  sequences (IUPAC-aware, forward and reverse complement). Returns a
  data frame suitable for use with
  [`tidypq::filter_taxa_pq()`](https://adrientaudiere.github.io/tidypq/reference/filter_taxa_pq.html)
  to prune contaminated taxa.

- [`community_sharing_barplot_pq()`](https://adrientaudiere.github.io/comparpq/reference/community_sharing_barplot_pq.md)
  new function to display pairwise community-sharing metrics as grouped
  bar charts, faceted by metric or by pair. Companion to
  [`community_sharing_pq()`](https://adrientaudiere.github.io/comparpq/reference/community_sharing_pq.md).

- [`community_sharing_pq()`](https://adrientaudiere.github.io/comparpq/reference/community_sharing_pq.md)
  new function to visualize community sharing between 2–4 modalities of
  a sample variable as a network figure: each node is a pie chart of
  taxonomic composition, and curved links encode multiple pairwise
  similarity metrics (Bray-Curtis, Jaccard, shared species, shared
  genera proportion). Supports label-permutation significance testing
  (`n_perm`). Requires packages `ggforce` and `RColorBrewer`.

- [`default_sharing_metrics()`](https://adrientaudiere.github.io/comparpq/reference/default_sharing_metrics.md)
  new function returning the 4 built-in metric definitions used by
  [`community_sharing_pq()`](https://adrientaudiere.github.io/comparpq/reference/community_sharing_pq.md)
  and
  [`community_sharing_barplot_pq()`](https://adrientaudiere.github.io/comparpq/reference/community_sharing_barplot_pq.md).

- [`make_sharing_metric()`](https://adrientaudiere.github.io/comparpq/reference/make_sharing_metric.md)
  new function to create custom metric definitions for
  [`community_sharing_pq()`](https://adrientaudiere.github.io/comparpq/reference/community_sharing_pq.md)
  and
  [`community_sharing_barplot_pq()`](https://adrientaudiere.github.io/comparpq/reference/community_sharing_barplot_pq.md).

- [`div_pq()`](https://adrientaudiere.github.io/comparpq/reference/div_pq.md)
  and `hill_samples_pq()` now use `divent` (via
  [`MiscMetabar::divent_hill_matrix_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/divent_hill_matrix_pq.html))
  instead of
  [`vegan::renyi()`](https://vegandevs.github.io/vegan/reference/renyi.html)
  for Hill number computation, and
  [`divent::ent_shannon()`](https://ericmarcon.github.io/divent/reference/ent_shannon.html)
  /
  [`divent::ent_simpson()`](https://ericmarcon.github.io/divent/reference/ent_simpson.html)
  instead of
  [`vegan::diversity()`](https://vegandevs.github.io/vegan/reference/diversity.html)
  for Shannon and Simpson indices. The default estimator is now
  `"UnveilJ"` (bias-corrected); pass `estimator = "naive"` to restore
  old numeric behavior.

- [`div_pq()`](https://adrientaudiere.github.io/comparpq/reference/div_pq.md):
  the `scales` parameter is deprecated in favour of `q`. The `hill`
  parameter is deprecated; only Hill numbers are now supported.

- [`add_external_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_external_seq_pq.md)
  now checks for a `refseq` slot upfront and emits a clear error when
  absent, instead of crashing with a cryptic message. It also strips the
  `phy_tree` slot before calling
  [`merge_phyloseq()`](https://rdrr.io/pkg/phyloseq/man/merge_phyloseq.html)
  to avoid tip-count mismatches on objects that carry a tree.

- [`add_shuffle_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_shuffle_seq_pq.md)
  now checks for a `refseq` slot upfront and emits a clear error when
  absent. It also strips the `phy_tree` slot before calling
  [`merge_phyloseq()`](https://rdrr.io/pkg/phyloseq/man/merge_phyloseq.html)
  to avoid tip-count mismatches.

- [`compare_refseq()`](https://adrientaudiere.github.io/comparpq/reference/compare_refseq.md)
  correctly handles `list_phyloseq` S7 objects by accessing
  `@phyloseq_list` directly instead of calling
  [`length()`](https://rdrr.io/r/base/length.html) on the S7 object
  itself.

- [`estim_cor_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_cor_pq.md)
  /
  [`estim_cor_lpq()`](https://adrientaudiere.github.io/comparpq/reference/estim_cor_lpq.md)
  bootstrap now passes `use = "complete.obs"` to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) and `na.rm = TRUE`
  to [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html),
  preventing NaN-induced crashes on degenerate resamples.

- [`estim_diff_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_diff_pq.md)
  now validates that each group has at least 3 samples before delegating
  to `dabestr`, providing an informative error message instead of a
  cryptic dabestr crash.

- [`rainplot_taxo_na()`](https://adrientaudiere.github.io/comparpq/reference/rainplot_taxo_na.md)
  now checks that requested rank columns exist in the
  [`psmelt()`](https://rdrr.io/pkg/phyloseq/man/psmelt.html) output
  before calling
  [`across()`](https://dplyr.tidyverse.org/reference/across.html),
  providing a clear error when all-NA rank columns are dropped.

- [`tc_heatmap()`](https://adrientaudiere.github.io/comparpq/reference/tc_heatmap.md)
  new function to visualize the correspondence between two taxonomic
  ranks as a heatmap, where each cell shows the number of taxa assigned
  to a given pair of rank values.

- [`taxtab_replace_pattern_by_NA()`](https://adrientaudiere.github.io/comparpq/reference/taxtab_replace_pattern_by_NA.md)
  fixes an inner-loop variable bug where patterns were applied to all
  `taxonomic_ranks` columns simultaneously instead of one at a time.

- [`tc_points_matrix()`](https://adrientaudiere.github.io/comparpq/reference/tc_points_matrix.md)
  now checks that requested rank columns exist in the
  [`psmelt()`](https://rdrr.io/pkg/phyloseq/man/psmelt.html) output
  before grouping, providing a clear error when all-NA rank columns are
  dropped.

- Add param `compute_dist` to
  [`list_phyloseq()`](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)

- [`length()`](https://rdrr.io/r/base/length.html),
  [`names()`](https://rdrr.io/r/base/names.html), `[()`, and `[[()` now
  work correctly on `list_phyloseq` objects. S7 stores the class
  attribute as `"comparpq::list_phyloseq"` (package-qualified), which
  prevented S3 dispatch from finding `length.list_phyloseq`. The
  constructor now prepends the bare `"list_phyloseq"` name to the class
  vector, enabling S3 dispatch. The four accessor methods are now
  documented and exported.

## comparpq 0.1.2

- Add params `significance`, `test` and `p_alpha` to
  [`div_pq()`](https://adrientaudiere.github.io/comparpq/reference/div_pq.md)
  to report tuckey hsd paired-test using letters.
- [`gg_hill_lpq()`](https://adrientaudiere.github.io/comparpq/reference/gg_hill_lpq.md)
  new function to visualize Hill diversity correlations across pairs of
  phyloseq objects in a `list_phyloseq`. Produces a faceted scatter plot
  (pairs × Hill orders) with optional 1:1 line, regression line, and
  per-panel correlation annotation, enabling visual assessment of
  REPRODUCIBILITY, ROBUSTNESS, and REPLICABILITY.
- [`gg_bubbles_pq()`](https://adrientaudiere.github.io/comparpq/reference/gg_bubbles_pq.md)
  now uses
  [`merge_lpq()`](https://adrientaudiere.github.io/comparpq/reference/merge_lpq.md)
  when a `list_phyloseq` is passed, merging into a single phyloseq and
  faceting by `source_name` instead of building separate patchwork
  panels. `diff_contour` now works with any `facet_by` variable (not
  just list_phyloseq), highlighting taxa unique to each facet level with
  a distinct contour color from `diff_contour_colors`. No longer limited
  to 2 or 3 objects. New `match_by` parameter controls how taxa are
  matched when merging list_phyloseq objects.
- [`gg_bubbles_pq()`](https://adrientaudiere.github.io/comparpq/reference/gg_bubbles_pq.md)
  new ggplot2-based circle-packed bubble plot of taxa abundances. Unlike
  [`bubbles_pq()`](https://adrientaudiere.github.io/comparpq/reference/bubbles_pq.md),
  it does not require d3js/Observable and supports faceting by a
  `@sam_data` variable to display one bubble chart per level. Uses
  `packcircles` for layout computation.
- [`compare_refseq()`](https://adrientaudiere.github.io/comparpq/reference/compare_refseq.md)
  new function to compare reference sequences (`refseq` slot) between
  two phyloseq objects, identifying shared and unique ASVs/OTUs by name
  and by DNA sequence content, including detection of
  same-name-different-sequence and same-sequence-different-name
  mismatches. Computes mean nearest-neighbor k-mer distance for unique
  sequences.
- Add function
  [`apply_to_lpq()`](https://adrientaudiere.github.io/comparpq/reference/apply_to_lpq.md)
  to apply a function to each phyloseq object in a list_phyloseq
- [`estim_cor_lpq()`](https://adrientaudiere.github.io/comparpq/reference/estim_cor_lpq.md)
  new function to compute bootstrap correlation/regression across a
  list_phyloseq
- [`estim_cor_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_cor_pq.md)
  new function to compute bootstrap correlation and regression CIs for
  diversity vs numeric variables
- [`estim_diff_lpq()`](https://adrientaudiere.github.io/comparpq/reference/estim_diff_lpq.md)
  new function to run estimation statistics (effect sizes + CIs) across
  a list_phyloseq
- [`estim_diff_pq()`](https://adrientaudiere.github.io/comparpq/reference/estim_diff_pq.md)
  new function for estimation statistics (Gardner-Altman/Cumming plots)
  comparing diversity across groups via dabestr
- [`merge_lpq()`](https://adrientaudiere.github.io/comparpq/reference/merge_lpq.md)
  new function to merge a list_phyloseq into a single phyloseq object
  where each original phyloseq becomes one sample. Taxa can be matched
  by reference sequences (`match_by = "refseq"`, default) or by taxa
  names (`match_by = "names"`).
- [`simple_venn_pq()`](https://adrientaudiere.github.io/comparpq/reference/simple_venn_pq.md)
  new function to draw Venn diagrams of shared taxa across 2-4 sample
  groups using pure ggplot2 (no external Venn package needed), with
  support for multiple taxonomic ranks and compact, clearly labeled
  circles/ellipses.
- [`simple_venn_pq()`](https://adrientaudiere.github.io/comparpq/reference/simple_venn_pq.md)
  now accepts a `list_phyloseq` object as input, automatically merging
  it via
  [`merge_lpq()`](https://adrientaudiere.github.io/comparpq/reference/merge_lpq.md)
  before drawing the Venn diagram.
