# Harmonize Species\_ columns and verify names via GNA Verifier

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

For each database suffix (e.g. `"_EUK"`, `"_Unite"`), (1) rewrites the
`Species<suffix>` column to keep only the specific epithet (see
[`extract_species_epithet()`](https://adrientaudiere.github.io/comparpq/reference/extract_species_epithet.md))
and, when `verify = TRUE`, (2) calls
[`taxinfo::gna_verifier_pq()`](https://adrientaudiere.github.io/taxinfo/reference/gna_verifier_pq.html)
on `c("Genus<suffix>", "Species<suffix>")` and appends the resulting
columns (`currentCanonicalSimple<suffix>`, `currentName<suffix>`, ...)
to the `tax_table`, renamed from `<db>_<col>` to `<col>_<db>` for
consistency with the `<Rank>_<Db>` naming.

The `verify = TRUE` path requires the suggested package taxinfo and
network access; the `verify = FALSE` path performs only the offline
epithet harmonization and needs no extra package.

## Usage

``` r
harmonize_sp_names_pq(
  physeq,
  suffixes,
  verify = TRUE,
  force_recompute = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq::phyloseq-class()`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object.

- suffixes:

  (character vector, required) Database suffixes (e.g.
  `c("_EUK", "_Unite")`). A leading underscore is added when missing.

- verify:

  (logical, default `TRUE`) When `FALSE`, only epithet harmonization is
  performed, without calling
  [`taxinfo::gna_verifier_pq()`](https://adrientaudiere.github.io/taxinfo/reference/gna_verifier_pq.html).

- force_recompute:

  (logical, default `FALSE`) When `TRUE`, drops existing result columns
  from a previous run before re-adding them; also forwarded to
  [`taxinfo::gna_verifier_pq()`](https://adrientaudiere.github.io/taxinfo/reference/gna_verifier_pq.html).

- ...:

  Additional arguments passed to
  [`taxinfo::gna_verifier_pq()`](https://adrientaudiere.github.io/taxinfo/reference/gna_verifier_pq.html)
  (e.g. `data_sources = 5` for Index Fungorum).

## Value

The `physeq` object with harmonized `Species<suffix>` columns and, when
`verify = TRUE`, the new `currentCanonicalSimple<suffix>`,
`currentName<suffix>`, ... columns.

## See also

[`extract_species_epithet()`](https://adrientaudiere.github.io/comparpq/reference/extract_species_epithet.md),
[`taxinfo::gna_verifier_pq()`](https://adrientaudiere.github.io/taxinfo/reference/gna_verifier_pq.html)

## Author

Adrien Taudière

## Examples

``` r
# Offline: keep only the epithet in a Species_<db> column.
pq <- Glom_otu
mat <- as(pq@tax_table, "matrix")
mat <- cbind(
  mat,
  Genus_EUK = rep("Quercus", nrow(mat)),
  Species_EUK = rep("Quercus_ilex", nrow(mat))
)
pq@tax_table <- phyloseq::tax_table(mat)
pq <- harmonize_sp_names_pq(pq, suffixes = "_EUK", verify = FALSE)

if (FALSE) { # \dontrun{
# verify = TRUE additionally validates names via taxinfo::gna_verifier_pq()
# (requires the 'taxinfo' package and network access).
pq <- harmonize_sp_names_pq(pq, suffixes = "_EUK", data_sources = 5)
} # }
```
