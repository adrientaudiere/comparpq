# Changelog

## comparpq 0.1.3 (Development version)

- Add param `compute_dist` to
  [`list_phyloseq()`](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)

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
