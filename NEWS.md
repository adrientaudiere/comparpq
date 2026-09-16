# comparpq 0.4.0 (Development version)

* `build_comparison_grid()` new function to enumerate, across a named list of phyloseq objects (or a `list_phyloseq` object), every pairwise combination of taxonomic databases per rank (`tax_table` columns named `"<Rank>_<Db>"`), producing a long-format grid (one row per phyloseq x rank x database pair) to drive systematic `compare_taxo_db()` or `tc_congruence_metrics()` comparisons.
* `compare_taxo_db()` new function to assemble, in a single call, the congruence metrics, a filtered contingency table, and the `tc_bar()`, `tc_sankey()`, `tc_heatmap()` and `tc_circle()` views comparing two taxonomic-assignment columns of a `phyloseq` object.
* `count_taxo_congruence()` new function to classify each taxon into a congruence category (`both_equal`, `both_na`, `only_<db1>`, `only_<db2>`, `different`) when comparing two `tax_table` columns, and to count both taxa and sequences (with percentages) per category.
* `extract_species_epithet()` new function to extract the specific epithet from a `Species` value, handling binomials, underscores and infraspecific parts.
* `harmonize_sp_names_pq()` new function to rewrite `Species_<db>` columns to their epithet and optionally verify names via `taxinfo::gna_verifier_pq()` (offline epithet-only mode with `verify = FALSE`).
* `plot_congruence_counts()` new function to draw a stacked barplot of taxonomic-assignment congruence across several ranks, keeping each database's color consistent across ranks.
* `tc_metrics_mock()` computes the MCC in double precision, so it is no longer `NA` (with an integer-overflow warning) when the product of the four confusion-matrix margins exceeds `.Machine$integer.max` (about 216 taxa per margin), and returns an MCC of 0 instead of `NaN` when a margin is 0 (e.g. no taxon assigned at the rank, or every negative control assigned), following Chicco & Jurman (2020).
* Fix missing `Remotes` field in `DESCRIPTION` so that `pak::pkg_install()` can resolve GitHub-only dependencies (`MiscMetabar`, `phylopq`, `taxinfo`) when installing comparpq as a transitive dependency of pqverse.

# comparpq 0.3.0
## Breaking changes

* `taxo2tree()` is removed from comparpq and relocated to the `phylopq` package, its natural home for phylogenetic tree construction from taxonomy tables. Calls to `comparpq::taxo2tree()` now fail with `could not find function`; use `phylopq::taxo2tree()` instead (the interface is unchanged).

## New features

* `merge_clust_lpq()` new function to merge a `list_phyloseq` into a single phyloseq object while keeping every sample separate, unifying taxa across objects by clustering their reference sequences (vsearch by default, via `MiscMetabar::postcluster_pq()`). Sample names are suffixed with the parent object name on collision, and a `source_name` column records the parent object of each sample. It complements `merge_lpq()`, which instead collapses each object into a single sample using exact sequence matching.

# comparpq 0.2.1
* `refseq_comp_lpq()` new function to compare `@refseq` sequences across all phyloseq objects in a `list_phyloseq` using k-mer Jaccard similarity and union-find connected components. Returns per-threshold Venn diagrams and shared-cluster counts. No igraph dependency.
* `find_primers_pq()` new function to detect taxa whose reference sequences match primer sequences (IUPAC-aware, forward and reverse complement). Returns a data frame suitable for use with `tidypq::filter_taxa_pq()` to prune contaminated taxa.

* `community_sharing_barplot_pq()` new function to display pairwise community-sharing metrics as grouped bar charts, faceted by metric or by pair. Companion to `community_sharing_pq()`.
* `community_sharing_pq()` new function to visualize community sharing between 2–4 modalities of a sample variable as a network figure: each node is a pie chart of taxonomic composition, and curved links encode multiple pairwise similarity metrics (Bray-Curtis, Jaccard, shared species, shared genera proportion). Supports label-permutation significance testing (`n_perm`). Requires packages `ggforce` and `RColorBrewer`.
* `default_sharing_metrics()` new function returning the 4 built-in metric definitions used by `community_sharing_pq()` and `community_sharing_barplot_pq()`.
* `make_sharing_metric()` new function to create custom metric definitions for `community_sharing_pq()` and `community_sharing_barplot_pq()`.

* `div_pq()` and `hill_samples_pq()` now use `divent` (via `MiscMetabar::divent_hill_matrix_pq()`) instead of `vegan::renyi()` for Hill number computation, and `divent::ent_shannon()` / `divent::ent_simpson()` instead of `vegan::diversity()` for Shannon and Simpson indices. The default estimator is now `"UnveilJ"` (bias-corrected); pass `estimator = "naive"` to restore old numeric behavior.
* `div_pq()`: the `scales` parameter is deprecated in favour of `q`. The `hill` parameter is deprecated; only Hill numbers are now supported.

* `add_external_seq_pq()` now checks for a `refseq` slot upfront and emits a clear error when absent, instead of crashing with a cryptic message. It also strips the `phy_tree` slot before calling `merge_phyloseq()` to avoid tip-count mismatches on objects that carry a tree.
* `add_shuffle_seq_pq()` now checks for a `refseq` slot upfront and emits a clear error when absent. It also strips the `phy_tree` slot before calling `merge_phyloseq()` to avoid tip-count mismatches.
* `compare_refseq()` correctly handles `list_phyloseq` S7 objects by accessing `@phyloseq_list` directly instead of calling `length()` on the S7 object itself.
* `estim_cor_pq()` / `estim_cor_lpq()` bootstrap now passes `use = "complete.obs"` to `stats::cor()` and `na.rm = TRUE` to `stats::quantile()`, preventing NaN-induced crashes on degenerate resamples.
* `estim_diff_pq()` now validates that each group has at least 3 samples before delegating to `dabestr`, providing an informative error message instead of a cryptic dabestr crash.
* `rainplot_taxo_na()` now checks that requested rank columns exist in the `psmelt()` output before calling `across()`, providing a clear error when all-NA rank columns are dropped.
* `tc_heatmap()` new function to visualize the correspondence between two taxonomic ranks as a heatmap, where each cell shows the number of taxa assigned to a given pair of rank values.
* `taxtab_replace_pattern_by_NA()` fixes an inner-loop variable bug where patterns were applied to all `taxonomic_ranks` columns simultaneously instead of one at a time.
* `tc_points_matrix()` now checks that requested rank columns exist in the `psmelt()` output before grouping, providing a clear error when all-NA rank columns are dropped.
* Add param `compute_dist` to `list_phyloseq()`
* `length()`, `names()`, `[()`, and `[[()` now work correctly on `list_phyloseq` objects. S7 stores the class attribute as `"comparpq::list_phyloseq"` (package-qualified), which prevented S3 dispatch from finding `length.list_phyloseq`. The constructor now prepends the bare `"list_phyloseq"` name to the class vector, enabling S3 dispatch. The four accessor methods are now documented and exported.


# comparpq 0.1.2

* Add params `significance`, `test` and `p_alpha` to `div_pq()` to report tuckey hsd paired-test using letters.
* `gg_hill_lpq()` new function to visualize Hill diversity correlations across pairs of phyloseq objects in a `list_phyloseq`. Produces a faceted scatter plot (pairs × Hill orders) with optional 1:1 line, regression line, and per-panel correlation annotation, enabling visual assessment of REPRODUCIBILITY, ROBUSTNESS, and REPLICABILITY.
* `gg_bubbles_pq()` now uses `merge_lpq()` when a `list_phyloseq` is passed, merging into a single phyloseq and faceting by `source_name` instead of building separate patchwork panels. `diff_contour` now works with any `facet_by` variable (not just list_phyloseq), highlighting taxa unique to each facet level with a distinct contour color from `diff_contour_colors`. No longer limited to 2 or 3 objects. New `match_by` parameter controls how taxa are matched when merging list_phyloseq objects.
* `gg_bubbles_pq()` new ggplot2-based circle-packed bubble plot of taxa abundances. Unlike `bubbles_pq()`, it does not require d3js/Observable and supports faceting by a `@sam_data` variable to display one bubble chart per level. Uses `packcircles` for layout computation.
* `compare_refseq()` new function to compare reference sequences (`refseq` slot) between two phyloseq objects, identifying shared and unique ASVs/OTUs by name and by DNA sequence content, including detection of same-name-different-sequence and same-sequence-different-name mismatches. Computes mean nearest-neighbor k-mer distance for unique sequences.
* Add function `apply_to_lpq()` to apply a function to each phyloseq object in a list_phyloseq
* `estim_cor_lpq()` new function to compute bootstrap correlation/regression across a list_phyloseq
* `estim_cor_pq()` new function to compute bootstrap correlation and regression CIs for diversity vs numeric variables
* `estim_diff_lpq()` new function to run estimation statistics (effect sizes + CIs) across a list_phyloseq
* `estim_diff_pq()` new function for estimation statistics (Gardner-Altman/Cumming plots) comparing diversity across groups via dabestr
* `merge_lpq()` new function to merge a list_phyloseq into a single phyloseq object where each original phyloseq becomes one sample. Taxa can be matched by reference sequences (`match_by = "refseq"`, default) or by taxa names (`match_by = "names"`).
* `simple_venn_pq()` new function to draw Venn diagrams of shared taxa across 2-4 sample groups using pure ggplot2 (no external Venn package needed), with support for multiple taxonomic ranks and compact, clearly labeled circles/ellipses.
* `simple_venn_pq()` now accepts a `list_phyloseq` object as input, automatically merging it via `merge_lpq()` before drawing the Venn diagram.

# Initial comparpq (0.01)

* `taxo2tree()` gains `use_taxa_names` parameter to exclude taxa names (e.g., ASV_1) as terminal leaves and use lowest rank values instead
* `tc_linked_trees()` new function to plot two taxonomy trees facing each other with linked correspondences between matching taxa
* `tc_linked_trees()` gains `link_by_taxa` parameter to draw links based on taxa correspondence with line width proportional to ASV count
* `tc_linked_trees()` gains `physeq_2 = NULL` default allowing comparison of different rank columns from the same phyloseq object
* `tc_linked_trees()` labels are now scaled by depth (shallower nodes larger) and positioned above branches to avoid overlap


