# Glomeromycota OTU phyloseq dataset

An example
[phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object used to illustrate the taxonomic-comparison helpers (`tc_*()`) of
comparpq. It bundles an OTU table, taxonomy table and sample data for a
set of Glomeromycota (arbuscular mycorrhizal fungi) communities.

## Usage

``` r
data(Glom_otu)
```

## Format

A
[phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with 1147 taxa and 444 samples, containing an `otu_table`, a
`tax_table` and a `sample_data` slot.

## Source

Derived from an environmental metabarcoding dataset; bundled for
documentation and examples only.
