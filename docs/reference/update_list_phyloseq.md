# Update the summary table and comparison characteristics

Recomputes the summary_table and comparison slots. Useful after
modifying the phyloseq objects in place or to change the comparison
parameters.

## Usage

``` r
update_list_phyloseq(
  x,
  same_primer_seq_tech = NULL,
  same_bioinfo_pipeline = NULL
)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- same_primer_seq_tech:

  (logical, default NULL) Whether the same primer and sequencing
  technology was used. If NULL, preserves the original value.

- same_bioinfo_pipeline:

  (logical, default NULL) Whether the same bioinformatics pipeline was
  used. If NULL, preserves the original value.

## Value

An updated list_phyloseq object

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)
