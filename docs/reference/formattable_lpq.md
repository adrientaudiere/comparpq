# Formattable visualization for list_phyloseq summary

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Create a visualization table to display the summary_table from a
list_phyloseq object using the formattable package with colored bars.

Display the summary_table from a list_phyloseq object as a formattable
table with colored bars for numeric columns and colored indicators for
logical columns.

## Usage

``` r
formattable_lpq(
  x,
  columns = c("name", "n_samples", "n_taxa", "n_sequences", "mean_seq_per_sample",
    "mean_seq_per_taxon", "has_sam_data", "has_tax_table", "has_refseq", "has_phy_tree"),
  bar_colors = NULL,
  round_digits = 1,
  void_style = FALSE,
  log10_transform = TRUE,
  ...
)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- columns:

  (character vector, default selection of key columns) Character vector
  of column names to display. If NULL, displays a curated selection of
  columns.

- bar_colors:

  (named list, default NULL) Named list of colors for numeric columns
  with bars. Names should match column names. Default colors are
  provided for common columns.

- round_digits:

  (integer, default 1) Number of decimal places for rounding numeric
  columns.

- void_style:

  (logical, default FALSE) If TRUE, returns a formattable without any
  custom styling.

- log10_transform:

  (logical, default TRUE) If TRUE, applies log10 transformation to
  numeric columns with a range greater than 1000.

- ...:

  Additional arguments passed to
  [`formattable::formattable()`](https://renkun-ken.github.io/formattable/reference/formattable.html).

## Value

A formattable object

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function is inspired by
[`MiscMetabar::formattable_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/formattable_pq.html).
Numeric columns are displayed with proportional colored bars. Logical
columns (has_sam_data, has_tax_table, etc.) are displayed with
checkmarks or X marks with colored backgrounds.

## Examples

``` r
if (FALSE) { # \dontrun{
lpq <- list_phyloseq(list(data1 = data_fungi, data2 = data_fungi_mini))
formattable_lpq(lpq)

# Custom columns
formattable_lpq(lpq,
  columns = c("name", "n_samples", "n_taxa", "n_sequences")
)

# Custom colors
formattable_lpq(lpq, bar_colors = list(
  n_samples = "steelblue",
  n_taxa = "darkgreen",
  n_sequences = "purple"
))
} # }
```
