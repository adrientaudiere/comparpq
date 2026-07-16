# Extract the species epithet from a Species value

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Returns the specific epithet from a value of a `Species` column,
handling the many formats found in reference databases: bare epithet
(`"iris"`), binomial with a space or an underscore (`"Genus iris"`,
`"Genus_iris"`), and trailing infraspecific parts
(`"Genus_iris_var._occidentalis"`). Underscores are treated as spaces.
When the value holds only the genus (matching `genus_val`), `NA` is
returned.

## Usage

``` r
extract_species_epithet(species_val, genus_val)
```

## Arguments

- species_val:

  (character) A single value of the `Species` column.

- genus_val:

  (character) The corresponding `Genus` value (may be `NA`), used to
  decide whether the first word is the genus.

## Value

The species epithet as a length-one character, or `NA_character_`.

## See also

[`harmonize_sp_names_pq()`](https://adrientaudiere.github.io/comparpq/reference/harmonize_sp_names_pq.md)

## Author

Adrien Taudière

## Examples

``` r
extract_species_epithet("Quercus_ilex", "Quercus")
#> [1] "ilex"
extract_species_epithet("Quercus ilex var. rotundifolia", "Quercus")
#> [1] "ilex"
extract_species_epithet("ilex", "Quercus")
#> [1] "ilex"
```
