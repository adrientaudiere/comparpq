# Compute accuracy metrics of taxonomic assignation using a mock (known) community for one rank

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Compute numerous metrics comparing the computed taxonomic assignation to
a true assignation.

Note that to compute all metrics, one need to insert fake taxa (by
shuffling sequences and/or by adding external sequences). The user must
fake taxa using functions
[`add_external_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_external_seq_pq.md),
[`add_shuffle_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_shuffle_seq_pq.md))
before taxonomic assignation.

## Usage

``` r
tc_metrics_mock_vec(
  physeq,
  taxonomic_rank,
  true_values,
  fake_taxa = TRUE,
  fake_pattern = c("^fake_", "^external_"),
  verbose = TRUE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- taxonomic_rank:

  (required) Name (or number) of a taxonomic rank to count.

- true_values:

  (required) A vector with the true taxonomic assignation

- fake_taxa:

  (logical, default TRUE) If TRUE, the fake_pattern vector is used to
  identify fake taxa, i.e. taxa who are not in the reference database
  (see
  [`add_external_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_external_seq_pq.md))
  or taxa with fake sequences (see
  [`add_shuffle_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_shuffle_seq_pq.md)).

- fake_pattern:

  (character vector, default c("^fake\_", "^external\_")) A vector of
  patterns used to identify the fake taxa using a regex search in their
  name.

- verbose:

  (logical, default TRUE) If TRUE, print informative messages.

## Value

A list of metrics (see the confusion matrix
[article](https://en.wikipedia.org/wiki/Confusion_matrix) on wikipedia):

- TP (number of *true positive*)

- FP (number of *false positive*)

- FN (number of *false negative*)

- FDR (*false discovery rate*) = FP / (FP + TP)

- TPR (*true positive rate*, also named *recall* or *sensitivity*)

- PPV (*positive predictive value*, also named *precision*) = TP / (TP +
  FP)

- F1_score (*F1 score*) = 2 \* TP / (2 \* TP + FP + FN)

If fake taxa are present and `fake_taxa` is true, other metrics are
computed:

- TN (number of *true negative*)

- ACC (*Accuracy*) = (TP + TN) / (TP + TN + FP + FN)

- MCC (*Matthews correlation coefficient*) = (TP \* TN - FP \* FN) /
  sqrt((TP + FP) \* (TP + FN) \* (FP + TN) \* (TN + FN))

## See also

[`tc_metrics_mock()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock.md),
[`add_external_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_external_seq_pq.md),
[`add_shuffle_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_shuffle_seq_pq.md))

## Author

Adrien Taudière
