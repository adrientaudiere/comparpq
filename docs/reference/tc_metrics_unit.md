# Accuracy metrics of several assignations against a per-unit truth

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Run
[`tc_metrics_unit_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_unit_vec.md)
on every assignation of `ranks_df` and every rank of `truth`, and add
the metrics of the external controls, which the per-rank function cannot
compute on its own:

- `ext_fungi`: share of the external controls named Fungi at the kingdom
  rank, an error on every database;

- `ext_correct`: share of the external controls given their own lineage,
  among those with a truth at the rank (NA in `external_truth`: left
  out, as a real unit below its truth depth); needs `external_truth`.

## Usage

``` r
tc_metrics_unit(
  physeq,
  ranks_df,
  truth,
  truth_ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  external_truth = NULL,
  external_correct_ranks = NULL,
  accepted_suffix = "_accepted",
  external_scoring = c("matrix", "aside"),
  fake_taxa = TRUE,
  fake_pattern = "^fake_",
  external_pattern = "^external_",
  kingdom_rank = "Kingdom",
  fungi_name = "Fungi",
  verbose = FALSE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- ranks_df:

  (required) A data.frame of `tax_table` column names: one column per
  assignation (method x database x parameters), one row per rank of
  `truth_ranks`, in the same order.

- truth:

  (required) The per-unit truth, as in
  [`tc_metrics_unit_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_unit_vec.md).

- truth_ranks:

  (character vector) The ranks of `truth`, from the highest to the
  lowest, used to tell whether `rank` is deeper than `truth_depth`.

- external_truth:

  (data.frame, default NULL) The true lineage of the external controls:
  a `taxon` column and one column per rank. Without it `ext_correct` is
  not computed.

- external_correct_ranks:

  (character vector, default NULL) Ranks where `ext_correct` is
  computed. NULL means every rank of `truth_ranks`; a database holding
  only a few non-fungal representatives should restrict it to the
  kingdom.

- accepted_suffix:

  (character, default "\_accepted") Suffix of the optional second-name
  columns of `truth`.

- external_scoring:

  ("matrix" or "aside") Whether the external controls enter the
  confusion matrix. They should when the reference database cannot hold
  them (a Fungi-only database), and stay aside when it can, where naming
  them is a correct answer rather than an error.

- fake_taxa:

  (logical, default TRUE) If TRUE, the controls are identified by
  `fake_pattern` and `external_pattern` and scored as above. If FALSE,
  every unit is treated as a real one.

- fake_pattern:

  (character, default "^fake\_") Regular expression identifying the
  shuffled controls
  ([`add_shuffle_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_shuffle_seq_pq.md)).

- external_pattern:

  (character, default "^external\_") Regular expression identifying the
  external controls
  ([`add_external_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_external_seq_pq.md)).

- kingdom_rank:

  (character, default "Kingdom") Rank at which `ext_fungi` is counted.

- fungi_name:

  (character, default "Fungi") Value of `kingdom_rank` counted by
  `ext_fungi`.

- verbose:

  (logical, default TRUE) If TRUE, print informative messages.

## Value

A long-format data.frame with four columns: `method_db`, `tax_level`,
`metrics` and `values`, as
[`tc_metrics_mock()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock.md)
returns.

## See also

[`tc_metrics_unit_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_unit_vec.md),
[`tc_metrics_mock()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock.md)

## Author

Adrien Taudière

## Examples

``` r
truth <- data.frame(
  taxon = phyloseq::taxa_names(data_fungi_mini),
  Phylum = phyloseq::tax_table(data_fungi_mini)[, "Phylum"],
  Class = phyloseq::tax_table(data_fungi_mini)[, "Class"],
  truth_depth = "Class"
)
ranks_df <- data.frame(seed = c("Phylum", "Class"))
tc_metrics_unit(
  data_fungi_mini,
  ranks_df = ranks_df,
  truth = truth,
  truth_ranks = c("Phylum", "Class"),
  fake_taxa = FALSE,
  verbose = FALSE
)
#>    method_db tax_level            metrics      values
#> 1       seed    Phylum                 TP 45.00000000
#> 2       seed    Phylum                 FP  0.00000000
#> 3       seed    Phylum                 FN  0.00000000
#> 4       seed    Phylum                 TN  0.00000000
#> 5       seed    Phylum                FDR  0.00000000
#> 6       seed    Phylum                PPV  1.00000000
#> 7       seed    Phylum                TPR  1.00000000
#> 8       seed    Phylum                TNR         NaN
#> 9       seed    Phylum           F1_score  1.00000000
#> 10      seed    Phylum                ACC  1.00000000
#> 11      seed    Phylum                MCC  0.00000000
#> 12      seed    Phylum            NA_real  0.00000000
#> 13      seed    Phylum            NA_fake         NaN
#> 14      seed    Phylum ctrl_assigned_fake  0.00000000
#> 15      seed    Phylum  ctrl_assigned_ext  0.00000000
#> 16      seed    Phylum      misassign_seq  0.00000000
#> 17      seed    Phylum             n_real 45.00000000
#> 18      seed    Phylum         n_left_out  0.00000000
#> 19      seed    Phylum          n_foreign  0.00000000
#> 20      seed    Phylum      foreign_named         NaN
#> 21      seed     Class                 TP 44.00000000
#> 22      seed     Class                 FP  0.00000000
#> 23      seed     Class                 FN  1.00000000
#> 24      seed     Class                 TN  0.00000000
#> 25      seed     Class                FDR  0.00000000
#> 26      seed     Class                PPV  1.00000000
#> 27      seed     Class                TPR  0.97777778
#> 28      seed     Class                TNR         NaN
#> 29      seed     Class           F1_score  0.98876404
#> 30      seed     Class                ACC  0.97777778
#> 31      seed     Class                MCC  0.00000000
#> 32      seed     Class            NA_real  0.02222222
#> 33      seed     Class            NA_fake         NaN
#> 34      seed     Class ctrl_assigned_fake  0.00000000
#> 35      seed     Class  ctrl_assigned_ext  0.00000000
#> 36      seed     Class      misassign_seq  0.00000000
#> 37      seed     Class             n_real 45.00000000
#> 38      seed     Class         n_left_out  0.00000000
#> 39      seed     Class          n_foreign  0.00000000
#> 40      seed     Class      foreign_named         NaN
```
