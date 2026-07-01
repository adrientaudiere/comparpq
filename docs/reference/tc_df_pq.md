# Contingency table of two taxonomic ranks

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Creates a cross-tabulation (contingency table) comparing two taxonomic
ranks from a phyloseq object. Useful for comparing taxonomic assignments
from different databases, algorithms, or taxonomic levels.

## Usage

``` r
tc_df_pq(physeq, rank_1 = "Family", rank_2 = "Class", ...)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- rank_1:

  (character, default "Family") The name of the first taxonomic rank
  (column in tax_table slot) for the cross-tabulation rows.

- rank_2:

  (character, default "Class") The name of the second taxonomic rank
  (column in tax_table slot) for the cross-tabulation columns.

- ...:

  Additional arguments passed to
  [`gtsummary::tbl_cross()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_cross.html).

## Value

A gtsummary tbl_cross object displaying the cross-tabulation of the two
taxonomic ranks.

## Author

Adrien Taudière

## Examples

``` r
tc_df_pq(data_fungi_mini)


  

```

Class

Total

Agaricomycetes

Atractiellomycetes

Tremellomycetes

Unknown

Family

  

  

  

  

  

    Aporpiaceae

1

0

0

0

1

    Atractiellales_fam_Incertae_sedis

0

1

0

0

1

    Auriculariaceae

2

0

0

0

2

    Cantharellales_fam_Incertae_sedis

1

0

0

0

1

    Corticiaceae

1

0

0

0

1

    Entolomataceae

1

0

0

0

1

    Exidiaceae

3

0

0

0

3

    Hericiaceae

1

0

0

0

1

    Hymenochaetales_fam_Incertae_sedis

1

0

0

0

1

    Hyphodermataceae

2

0

0

0

2

    Lyophyllaceae

4

0

0

0

4

    Peniophoraceae

1

0

0

0

1

    Phanerochaetaceae

1

0

0

0

1

    Polyporaceae

5

0

0

0

5

    Pterulaceae

1

0

0

0

1

    Russulales_fam_Incertae_sedis

1

0

0

0

1

    Schizoporaceae

4

0

0

0

4

    Steccherinaceae

1

0

0

0

1

    Stereaceae

6

0

0

0

6

    Tricholomataceae

1

0

0

0

1

    Unknown

3

0

2

1

6

Total

41

1

2

1

45

tc_df_pq(data_fungi_mini, rank_1 = "Order", rank_2 = "Family")

[TABLE]

if (FALSE) { \# \dontrun{ \# Compare taxonomic assignments from
different methods ref_fasta \<-
[system.file](https://rdrr.io/r/base/system.file.html)("extdata",
"mini_UNITE_fungi.fasta.gz", package = "MiscMetabar", mustWork = TRUE )
data_fungi_mini2 \<- data_fungi_mini \|\>
[add_new_taxonomy_pq](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.html)(ref_fasta,
suffix = "\_sintax", method = "sintax") \|\>
[add_new_taxonomy_pq](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.html)(ref_fasta,
suffix = "\_lca", method = "lca") tc_df_pq(data_fungi_mini2, rank_1 =
"Class_lca", rank_2 = "Class_sintax") } \# }
