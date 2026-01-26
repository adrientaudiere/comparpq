# Extended formattable for list_phyloseq with comparison info

Create an extended formattable table that also displays comparison
characteristics from a list_phyloseq object.

## Usage

``` r
formattable_lpq_full(x, show_summary = TRUE, show_comparison = TRUE, ...)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- show_summary:

  (logical, default TRUE) If TRUE, show the summary table.

- show_comparison:

  (logical, default TRUE) If TRUE, show comparison info.

- ...:

  Additional arguments passed to
  [`formattable_lpq()`](https://adrientaudiere.github.io/comparpq/reference/formattable_lpq.md).

## Value

A list containing formattable objects for summary and comparison, or a
single formattable if only one is requested.

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)
