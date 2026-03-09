# comparpq 0.1.0 (development version)

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


