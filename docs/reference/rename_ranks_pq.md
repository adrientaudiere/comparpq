# Rename names of ranks in the tax_table slot of a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The users can use the function using a couple of vector
`old_names`/`new_names` or using a `pattern` to replace by `replacement`

## Usage

``` r
rename_ranks_pq(
  physeq,
  old_names = NULL,
  new_names = NULL,
  pattern = NULL,
  replacement = NULL,
  fixed = FALSE,
  perl = FALSE,
  useBytes = FALSE
)
```

## Arguments

- physeq:

  (required) A phyloseq object.

- old_names:

  (character vector, default NULL) Names of the columns to replace.

- new_names:

  (character vector, default NULL) New names for the columns.

- pattern:

  (character, default NULL) Pattern to replace by the `replacement`
  argument.

- replacement:

  (character, default NULL) Replacement string for the pattern.

- fixed:

  (logical, default FALSE) See ?grep.

- perl:

  (logical, default FALSE) See ?grep.

- useBytes:

  (logical, default FALSE) See ?grep.

## Value

An object of class phyloseq

## Author

Adrien Taudière

## Examples

``` r
rename_ranks_pq(data_fungi, c("Confidence.Ranking", "Phylum"), c("Conf.Rank.Guild", "Phyla"))@tax_table |>
  colnames()
#>  [1] "Domain"          "Phyla"           "Class"           "Order"          
#>  [5] "Family"          "Genus"           "Species"         "Trophic.Mode"   
#>  [9] "Guild"           "Trait"           "Conf.Rank.Guild" "Genus_species"  
rename_ranks_pq(data_fungi, pattern = ".", replacement = "_", fixed = TRUE)@tax_table |>
  colnames()
#>  [1] "Domain"             "Phylum"             "Class"             
#>  [4] "Order"              "Family"             "Genus"             
#>  [7] "Species"            "Trophic_Mode"       "Guild"             
#> [10] "Trait"              "Confidence_Ranking" "Genus_species"     
```
