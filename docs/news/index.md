# Changelog

## comparpq 0.4.0 (Development version)

- [`tc_metrics_unit()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_unit.md)
  and
  [`tc_metrics_unit_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_unit_vec.md)
  new functions scoring a taxonomic assignation against a **per-unit
  truth** (one expected lineage per ASV or OTU) as Hleap et al. 2021 do
  in `optimize_n_score.py::score()`: every scored unit falls in exactly
  one cell of the confusion matrix, so TP + FP + FN + TN is the number
  of scored units. \[tc_metrics_mock()\] keeps its set-membership
  scoring for mocks that only come with a list of expected taxa. The
  external-control rate `ext_correct` is counted over the controls that
  have a truth at the rank. An optional logical `foreign` column of the
  truth marks the units foreign to the expected community: they leave
  the matrix at every rank, and `n_foreign` and `foreign_named` report
  how many there are and the share a method names.
- [`build_comparison_grid()`](https://adrientaudiere.github.io/comparpq/reference/build_comparison_grid.md)
  new function to enumerate, across a named list of phyloseq objects (or
  a `list_phyloseq` object), every pairwise combination of taxonomic
  databases per rank (`tax_table` columns named `"<Rank>_<Db>"`),
  producing a long-format grid (one row per phyloseq x rank x database
  pair) to drive systematic
  [`compare_taxo_db()`](https://adrientaudiere.github.io/comparpq/reference/compare_taxo_db.md)
  or
  [`tc_congruence_metrics()`](https://adrientaudiere.github.io/comparpq/reference/tc_congruence_metrics.md)
  comparisons.
- [`compare_taxo_db()`](https://adrientaudiere.github.io/comparpq/reference/compare_taxo_db.md)
  new function to assemble, in a single call, the congruence metrics, a
  filtered contingency table, and the
  [`tc_bar()`](https://adrientaudiere.github.io/comparpq/reference/tc_bar.md),
  [`tc_sankey()`](https://adrientaudiere.github.io/comparpq/reference/tc_sankey.md),
  [`tc_heatmap()`](https://adrientaudiere.github.io/comparpq/reference/tc_heatmap.md)
  and
  [`tc_circle()`](https://adrientaudiere.github.io/comparpq/reference/tc_circle.md)
  views comparing two taxonomic-assignment columns of a `phyloseq`
  object.
- [`count_taxo_congruence()`](https://adrientaudiere.github.io/comparpq/reference/count_taxo_congruence.md)
  new function to classify each taxon into a congruence category
  (`both_equal`, `both_na`, `only_<db1>`, `only_<db2>`, `different`)
  when comparing two `tax_table` columns, and to count both taxa and
  sequences (with percentages) per category.
- [`extract_species_epithet()`](https://adrientaudiere.github.io/comparpq/reference/extract_species_epithet.md)
  new function to extract the specific epithet from a `Species` value,
  handling binomials, underscores and infraspecific parts.
- [`harmonize_sp_names_pq()`](https://adrientaudiere.github.io/comparpq/reference/harmonize_sp_names_pq.md)
  new function to rewrite `Species_<db>` columns to their epithet and
  optionally verify names via
  [`taxinfo::gna_verifier_pq()`](https://adrientaudiere.github.io/taxinfo/reference/gna_verifier_pq.html)
  (offline epithet-only mode with `verify = FALSE`).
- [`plot_congruence_counts()`](https://adrientaudiere.github.io/comparpq/reference/plot_congruence_counts.md)
  new function to draw a stacked barplot of taxonomic-assignment
  congruence across several ranks, keeping each database’s color
  consistent across ranks.
- [`simple_venn_pq()`](https://adrientaudiere.github.io/comparpq/reference/simple_venn_pq.md)
  now orders the groups following the levels of `fact` when it is a
  factor, as its documentation states, instead of their order of first
  appearance in `sam_data`. The custom `labels` are therefore matched to
  the right groups: before, a label could be printed on the ellipse and
  sample count of another group. Character columns keep the order of
  first appearance.
- [`simple_venn_pq()`](https://adrientaudiere.github.io/comparpq/reference/simple_venn_pq.md)
  keeps long group names inside the figure: names on the sides of the
  diagram are pushed outward instead of running over the ellipses, and
  the plotting area is widened according to the length of the names, so
  they no longer cross the panel frame (or the next panel when
  `combine = TRUE`).
- [`tc_metrics_mock()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock.md)
  computes the MCC in double precision, so it is no longer `NA` (with an
  integer-overflow warning) when the product of the four
  confusion-matrix margins exceeds `.Machine$integer.max` (about 216
  taxa per margin), and returns an MCC of 0 instead of `NaN` when a
  margin is 0 (e.g. no taxon assigned at the rank, or every negative
  control assigned), following Chicco & Jurman (2020).
- Fix missing `Remotes` field in `DESCRIPTION` so that
  [`pak::pkg_install()`](https://pak.r-lib.org/reference/pkg_install.html)
  can resolve GitHub-only dependencies (`MiscMetabar`, `phylopq`,
  `taxinfo`) when installing comparpq as a transitive dependency of
  pqverse.

## comparpq 0.3.0

### Breaking changes

- `taxo2tree()` is removed from comparpq and relocated to the `phylopq`
  package, its natural home for phylogenetic tree construction from
  taxonomy tables. Calls to `comparpq::taxo2tree()` now fail with
  `could not find function`; use
  [`phylopq::taxo2tree()`](https://adrientaudiere.github.io/phylopq/reference/taxo2tree.html)
  instead (the interface is unchanged).

### New features

- [`merge_clust_lpq()`](https://adrientaudiere.github.io/comparpq/reference/merge_clust_lpq.md)
  new function to merge a `list_phyloseq` into a single phyloseq object
  while keeping every sample separate, unifying taxa across objects by
  clustering their reference sequences (vsearch by default, via
  [`MiscMetabar::postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.html)).
  Sample names are suffixed with the parent object name on collision,
  and a `source_name` column records the parent object of each sample.
  It complements
  [`merge_lpq()`](https://adrientaudiere.github.io/comparpq/reference/merge_lpq.md),
  which instead collapses each object into a single sample using exact
  sequence matching.

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
