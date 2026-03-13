# Compute accuracy metrics of multiple taxonomic assignations method using mock for multi-rank and multi assignation methods

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
tc_metrics_mock(
  physeq,
  ranks_df,
  true_values_df,
  fake_taxa = TRUE,
  fake_pattern = c("^fake_", "^external_")
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- ranks_df:

  (required). A dataframe with at least one column (one database or one
  method) and a number of row equal to the column in true_values_df

- true_values_df:

  (required). A dataframe with the true taxonomic assignation. Note that
  the column names (names of taxonomic ranks) of the true_values_df
  defined the names present in the `tax_level` column of the resulting
  dataframe.

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

## Value

A long-format dataframe with 4 columns: (i) the name of the `method_db`
(ii) the name of the `tax_level` (taxonomic rank), (iii) the `metrics`
(see
[`tc_metrics_mock_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock_vec.md)
for more details) and (iv) the `values`.

## See also

[`tc_metrics_mock_vec()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock_vec.md)

## Author

Adrien Taudière
