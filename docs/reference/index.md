# Package index

## list_phyloseq Class

S7 class for storing and comparing multiple phyloseq objects

- [`list_phyloseq()`](https://adrientaudiere.github.io/comparpq/reference/list_phyloseq.md)
  : S7 class for comparing phyloseq objects

## list_phyloseq Utilities

Functions for manipulating list_phyloseq objects

- [`add_phyloseq()`](https://adrientaudiere.github.io/comparpq/reference/add_phyloseq.md)
  : Add a phyloseq object to a list_phyloseq
- [`filter_common_lpq()`](https://adrientaudiere.github.io/comparpq/reference/filter_common_lpq.md)
  : Filter phyloseq objects to keep only shared samples and/or taxa
- [`n_levels_lpq()`](https://adrientaudiere.github.io/comparpq/reference/n_levels_lpq.md)
  : Count unique taxonomic levels across phyloseq objects
- [`remove_phyloseq()`](https://adrientaudiere.github.io/comparpq/reference/remove_phyloseq.md)
  : Remove a phyloseq object from a list_phyloseq
- [`shared_mod_lpq()`](https://adrientaudiere.github.io/comparpq/reference/shared_mod_lpq.md)
  : Display shared sample_data modalities
- [`update_list_phyloseq()`](https://adrientaudiere.github.io/comparpq/reference/update_list_phyloseq.md)
  : Update the summary table and comparison characteristics

## list_phyloseq Visualization

Visualization functions for list_phyloseq objects

- [`formattable_lpq()`](https://adrientaudiere.github.io/comparpq/reference/formattable_lpq.md)
  : Formattable visualization for list_phyloseq summary
- [`formattable_lpq_full()`](https://adrientaudiere.github.io/comparpq/reference/formattable_lpq_full.md)
  : Extended formattable for list_phyloseq with comparison info
- [`upset_lpq()`](https://adrientaudiere.github.io/comparpq/reference/upset_lpq.md)
  : UpSet or Venn plot of shared taxonomic values across phyloseq
  objects

## list_phyloseq Analysis

Statistical analysis functions for list_phyloseq objects

- [`adonis_lpq()`](https://adrientaudiere.github.io/comparpq/reference/adonis_lpq.md)
  : PERMANOVA analysis on each phyloseq object in a list_phyloseq
- [`glmulti_lpq()`](https://adrientaudiere.github.io/comparpq/reference/glmulti_lpq.md)
  : Automated model selection for Hill diversity on each phyloseq in a
  list_phyloseq

## Taxonomic Utilities

Functions for utilities on taxonomic data

- [`taxo2tree()`](https://adrientaudiere.github.io/comparpq/reference/taxo2tree.md)
  : Convert taxonomy dataframe to phylogenetic tree

## Taxonomic Comparison Metrics

Functions for computing accuracy metrics comparing taxonomic assignments

- [`tc_congruence_metrics()`](https://adrientaudiere.github.io/comparpq/reference/tc_congruence_metrics.md)
  : Compute congruence metrics between two taxonomic assignments
- [`tc_metrics_mock()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock.md)
  : Compute accuracy metrics of multiple taxonomic assignations method
  using mock for multi-rank and multi assignation methods
- [`tc_metrics_mock_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock_vec.md)
  : Compute accuracy metrics of taxonomic assignation using a mock
  (known) community for one rank

## Taxonomic Comparison Plots

Visualization functions for comparing taxonomic assignments

- [`tc_df_pq()`](https://adrientaudiere.github.io/comparpq/reference/tc_df_pq.md)
  : Contingency table of two taxonomic ranks
- [`tc_bar()`](https://adrientaudiere.github.io/comparpq/reference/tc_bar.md)
  : Barchart of ratio to compare 2 taxonomic ranks
- [`tc_circle()`](https://adrientaudiere.github.io/comparpq/reference/tc_circle.md)
  : Circle of correspondence between two taxonomic levels
- [`tc_points_matrix()`](https://adrientaudiere.github.io/comparpq/reference/tc_points_matrix.md)
  : Matrix of point to compare two taxonomic ranks
- [`rainplot_taxo_na()`](https://adrientaudiere.github.io/comparpq/reference/rainplot_taxo_na.md)
  : Rainplot of the nb taxa assigned (not NA)
- [`tc_linked_trees()`](https://adrientaudiere.github.io/comparpq/reference/tc_linked_trees.md)
  : Plot two taxonomy trees with linked correspondences
- [`tc_sankey()`](https://adrientaudiere.github.io/comparpq/reference/tc_sankey.md)
  : Sankey diagram to compare two taxonomic ranks

## Phyloseq Visualization

Interactive plotting functions for phyloseq objects

- [`bubbles_pq()`](https://adrientaudiere.github.io/comparpq/reference/bubbles_pq.md)
  :

  Bubble plot of phyloseq object with observablehq
  [![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Fake Taxa Creation

Functions for adding fake taxa to test taxonomic assignment accuracy

- [`add_external_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_external_seq_pq.md)
  : Add external sequences to a phyloseq object
- [`add_shuffle_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_shuffle_seq_pq.md)
  : Add fake sequences by shuffling existing ones in a phyloseq object

## Tax Table Manipulation

Functions for modifying the tax_table slot of phyloseq objects

- [`resolve_taxo_conflict()`](https://adrientaudiere.github.io/comparpq/reference/resolve_taxo_conflict.md)
  : Resolve taxonomic conflict in the tax_table of a phyloseq object
- [`select_ranks_pq()`](https://adrientaudiere.github.io/comparpq/reference/select_ranks_pq.md)
  : Select taxonomic ranks in a phyloseq object
- [`rename_ranks_pq()`](https://adrientaudiere.github.io/comparpq/reference/rename_ranks_pq.md)
  : Rename names of ranks in the tax_table slot of a phyloseq object
- [`taxtab_replace_pattern_by_NA()`](https://adrientaudiere.github.io/comparpq/reference/taxtab_replace_pattern_by_NA.md)
  : Replace taxonomic value with a given pattern by NA

## Formatting Utilities

Helper functions for formatting tables

- [`factor_formatter()`](https://adrientaudiere.github.io/comparpq/reference/factor_formatter.md)
  : Format factor columns with funky colored backgrounds
