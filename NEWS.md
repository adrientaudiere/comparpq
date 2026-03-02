# comparpq 0.1.0 (development version)

* Add function `apply_to_lpq()` to apply a function to each phyloseq object in a list_phyloseq
* `estim_cor_lpq()` new function to compute bootstrap correlation/regression across a list_phyloseq
* `estim_cor_pq()` new function to compute bootstrap correlation and regression CIs for diversity vs numeric variables
* `estim_diff_lpq()` new function to run estimation statistics (effect sizes + CIs) across a list_phyloseq
* `estim_diff_pq()` new function for estimation statistics (Gardner-Altman/Cumming plots) comparing diversity across groups via dabestr

# Initial comparpq (0.01)

* `taxo2tree()` gains `use_taxa_names` parameter to exclude taxa names (e.g., ASV_1) as terminal leaves and use lowest rank values instead
* `tc_linked_trees()` new function to plot two taxonomy trees facing each other with linked correspondences between matching taxa
* `tc_linked_trees()` gains `link_by_taxa` parameter to draw links based on taxa correspondence with line width proportional to ASV count
* `tc_linked_trees()` gains `physeq_2 = NULL` default allowing comparison of different rank columns from the same phyloseq object
* `tc_linked_trees()` labels are now scaled by depth (shallower nodes larger) and positioned above branches to avoid overlap


