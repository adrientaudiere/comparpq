# Accuracy metrics of a taxonomic assignation against a per-unit truth

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Score one assignation column at one rank as Hleap et al. (2021) do
(`optimize_n_score.py::score()`): every scored unit falls in exactly one
cell of the confusion matrix, so TP + FP + FN + TN is the number of
scored units.

|                                             |                   |                                                                |
|---------------------------------------------|-------------------|----------------------------------------------------------------|
| Unit                                        | Value at the rank | Cell                                                           |
| real, with a truth at this rank             | its own truth     | TP                                                             |
| real, with a truth at this rank             | another name      | FP                                                             |
| real, with a truth at this rank             | NA                | FN                                                             |
| real, without truth                         | any name / NA     | FP / FN                                                        |
| real, foreign to the mock (`truth$foreign`) | any name / NA     | left out, counted by `foreign_named`                           |
| real, truth shallower than the rank         | —                 | left out                                                       |
| shuffled control (`fake_pattern`)           | any name / NA     | FP / TN                                                        |
| external control (`external_pattern`)       | any name / NA     | FP / TN when `external_scoring = "matrix"`, left out otherwise |

This differs from
[`tc_metrics_mock_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock_vec.md),
which compares each value to the *set* of the expected taxa, counts FN
on the rows of the truth table and leaves the controls out of FP. Use
this one when every unit has its own truth (e.g. matched against the
Sanger sequences of the mock strains), and the other when the mock only
comes with a list of expected taxa.

## Usage

``` r
tc_metrics_unit_vec(
  physeq,
  taxonomic_rank,
  truth,
  rank = taxonomic_rank,
  truth_ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  accepted_suffix = "_accepted",
  external_scoring = c("matrix", "aside"),
  fake_taxa = TRUE,
  fake_pattern = "^fake_",
  external_pattern = "^external_",
  verbose = TRUE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- taxonomic_rank:

  (required) Name (or number) of the `tax_table` column holding the
  assignation to score.

- truth:

  (required) A data.frame with one row per unit: a `taxon` column naming
  the unit, one column per rank of `truth_ranks` holding its truth, and
  a `truth_depth` column naming the deepest rank with a truth (`NA` when
  the unit has none). A unit is scored down to `truth_depth` and left
  out of the matrix below it. Units absent from `truth` are scored as
  units without truth. Optional columns `<rank>_accepted` hold a second
  name that counts as correct too (e.g. the accepted name of a synonym).
  An optional logical column `foreign` marks the units foreign to the
  expected community (e.g. matching no Sanger sequence of the mock, even
  remotely): they leave the matrix at every rank, and `foreign_named`
  reports the share a method names.

- rank:

  (character, default `taxonomic_rank`) Rank of `truth` this column is
  compared to, when the column is named after the method rather than the
  rank (e.g. `"Genus_dada2__unite"` against `"Genus"`).

- truth_ranks:

  (character vector) The ranks of `truth`, from the highest to the
  lowest, used to tell whether `rank` is deeper than `truth_depth`.

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

- verbose:

  (logical, default TRUE) If TRUE, print informative messages.

## Value

A list of metrics:

- TP, FP, FN, TN: the four cells, in units;

- FDR = FP / (FP + TP), PPV = TP / (TP + FP), TPR = TP / (TP + FN), TNR
  = TN / (TN + FP) (every FP, controls and real misassignments alike),
  F1_score = 2 TP / (2 TP + FP + FN), ACC and MCC (0 when a margin is
  zero, Chicco & Jurman 2020);

- NA_real: share of the scored real units left NA;

- NA_fake: share of the shuffled controls left NA (1 is the correct
  answer);

- ctrl_assigned_fake, ctrl_assigned_ext: number of controls given a
  value;

- misassign_seq: share of the scored real units given a wrong name;

- n_real, n_left_out: number of real units scored, and of units left out
  of the matrix at this rank (truth shallower than the rank);

- n_foreign, foreign_named: number of real units marked `foreign` in
  `truth`, and the share of them given a value at this rank (NaN when
  there is none).

## See also

[`tc_metrics_unit()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_unit.md),
[`tc_metrics_mock_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock_vec.md)

## Author

Adrien Taudière

## Examples

``` r
truth <- data.frame(
  taxon = phyloseq::taxa_names(data_fungi_mini),
  Phylum = phyloseq::tax_table(data_fungi_mini)[, "Phylum"],
  truth_depth = "Phylum"
)
tc_metrics_unit_vec(
  data_fungi_mini,
  taxonomic_rank = "Phylum",
  truth = truth,
  truth_ranks = "Phylum",
  fake_taxa = FALSE,
  verbose = FALSE
)
#> $TP
#> [1] 45
#> 
#> $FP
#> [1] 0
#> 
#> $FN
#> [1] 0
#> 
#> $TN
#> [1] 0
#> 
#> $FDR
#> [1] 0
#> 
#> $PPV
#> [1] 1
#> 
#> $TPR
#> [1] 1
#> 
#> $TNR
#> [1] NaN
#> 
#> $F1_score
#> [1] 1
#> 
#> $ACC
#> [1] 1
#> 
#> $MCC
#> [1] 0
#> 
#> $NA_real
#> [1] 0
#> 
#> $NA_fake
#> [1] NaN
#> 
#> $ctrl_assigned_fake
#> [1] 0
#> 
#> $ctrl_assigned_ext
#> [1] 0
#> 
#> $misassign_seq
#> [1] 0
#> 
#> $n_real
#> [1] 45
#> 
#> $n_left_out
#> [1] 0
#> 
#> $n_foreign
#> [1] 0
#> 
#> $foreign_named
#> [1] NaN
#> 
```
