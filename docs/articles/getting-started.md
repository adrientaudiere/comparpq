# Getting Started with comparpq

## Introduction

**comparpq** is an extension package to
[{MiscMetabar}](https://adrientaudiere.github.io/MiscMetabar/) and
[{phyloseq}](https://joey711.github.io/phyloseq/) designed for comparing
phyloseq objects. This package is particularly useful for:

- Evaluating taxonomic assignment accuracy using mock communities
- Creating interactive visualizations of microbiome data
- Comparing multiple taxonomic databases or assignment methods
- Data preprocessing and quality control

## Installation

You can install comparpq from GitHub:

``` r
# Install from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("adrientaudiere/comparpq")
```

## Loading Required Libraries

``` r
library(comparpq)

# The package automatically loads dependencies including:
# - phyloseq for handling microbiome data
# - MiscMetabar for additional microbiome analysis tools
```

## Working with Example Data

The package includes the `Glom_otu` dataset, which contains glomalean
(arbuscular mycorrhizal fungi) OTU data:

``` r
Glom_otu
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1147 taxa and 444 samples ]
#> sample_data() Sample Data:       [ 444 samples by 118 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1147 taxa by 14 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1147 reference sequences ]
```

``` r
# Basic information about the dataset
cat("Number of taxa:", phyloseq::ntaxa(Glom_otu), "\n")
#> Number of taxa: 1147
cat("Number of samples:", phyloseq::nsamples(Glom_otu), "\n")
#> Number of samples: 444
cat("Taxonomic ranks:", paste(phyloseq::rank_names(Glom_otu), collapse = ", "), "\n")
#> Taxonomic ranks: Kingdom, Phyla, Class, Order, Family, Genus, Species, Kingdom__eukaryome_Glomero, Phylum__eukaryome_Glomero, Class__eukaryome_Glomero, Order__eukaryome_Glomero, Family__eukaryome_Glomero, Genus__eukaryome_Glomero, Species__eukaryome_Glomero
```

## Interactive Visualizations

### Bubble Plots with Observable HQ

One of the key features of comparpq is the ability to create interactive
bubble plots using Observable HQ:

``` r
# Basic bubble plot
bubbles_pq(Glom_otu)
```

``` r
# Bubble plot colored by Family
bubbles_pq(Glom_otu,
  rank_color = "Family",
  min_nb_seq = 1000
)
```

``` r
# Customized bubble plot with different color scheme
bubbles_pq(Glom_otu,
  rank_color = "Class",
  categorical_scheme = "d3.schemePastel1",
  min_nb_seq = 500,
  log1ptransform = TRUE
)
```

### Traditional Plots

For non-interactive visualizations, you can use:

``` r
# Create comparison matrices for visualization
matrix_data <- tc_points_matrix(Glom_otu)

# Bar chart comparison (requires appropriate data structure)
# tc_bar(matrix_data)

# Circular comparison plot
# tc_circle(matrix_data)
```

## Data Manipulation and Preprocessing

### Working with Taxonomic Tables

``` r
# Select specific taxonomic ranks
Glom_otu_subset <- select_ranks_pq(Glom_otu, Kingdom, Genus)
cat("Original ranks:", paste(phyloseq::rank_names(Glom_otu), collapse = ", "), "\n")
#> Original ranks: Kingdom, Phyla, Class, Order, Family, Genus, Species, Kingdom__eukaryome_Glomero, Phylum__eukaryome_Glomero, Class__eukaryome_Glomero, Order__eukaryome_Glomero, Family__eukaryome_Glomero, Genus__eukaryome_Glomero, Species__eukaryome_Glomero
cat("Selected ranks:", paste(phyloseq::rank_names(Glom_otu_subset), collapse = ", "), "\n")
#> Selected ranks: Kingdom, Genus
```

``` r
# Rename taxonomic ranks
renamed_physeq <- rename_ranks_pq(Glom_otu,
  new_names = c(
    "Domain", "Division", "Class_level",
    "Order_level", "Family_level", "Genus_level"
  )
)
cat("Renamed ranks:", paste(phyloseq::rank_names(renamed_physeq), collapse = ", "), "\n")
#> Renamed ranks: Kingdom, Phyla, Class, Order, Family, Genus, Species, Kingdom__eukaryome_Glomero, Phylum__eukaryome_Glomero, Class__eukaryome_Glomero, Order__eukaryome_Glomero, Family__eukaryome_Glomero, Genus__eukaryome_Glomero, Species__eukaryome_Glomero
```

### Pattern Replacement and Data Cleaning

``` r
# Replace unwanted patterns with NA
# This is useful for cleaning taxonomic assignments
cleaned_physeq <- taxtab_replace_pattern_by_NA(Glom_otu,
  pattern = "^uncultured",
  taxonomic_ranks = "Genus"
)

# Check how many taxa were affected
original_genus <- phyloseq::tax_table(Glom_otu)[, "Genus"]
cleaned_genus <- phyloseq::tax_table(cleaned_physeq)[, "Genus"]
cat(
  "Uncultured genera replaced with NA:",
  sum(is.na(cleaned_genus)) - sum(is.na(original_genus)), "\n"
)
#> Uncultured genera replaced with NA: 0
```

## Mock Community Analysis

### Creating Fake Taxa for Accuracy Assessment

To evaluate taxonomic assignment accuracy, you can add fake taxa to your
dataset:

``` r
# Add shuffled sequences (creates fake taxa from existing sequences)
physeq_with_shuffled <- add_shuffle_seq_pq(Glom_otu, n_fake = 10)
cat("Original taxa count:", phyloseq::ntaxa(Glom_otu), "\n")
#> Original taxa count: 1147
cat("With fake taxa:", phyloseq::ntaxa(physeq_with_shuffled), "\n")
#> With fake taxa: 1157

# Add external sequences (simulates contamination or database gaps)
external_seq <- Biostrings::readDNAStringSet(
  system.file("extdata/ex_little.fasta", package = "MiscMetabar")
)
physeq_with_external <- add_external_seq_pq(Glom_otu, ext_seqs = external_seq)
cat("With external sequences added:", phyloseq::ntaxa(physeq_with_external), "\n")
#> With external sequences added: 1149
```

### Accuracy Metrics Computation

``` r
# Example of computing accuracy metrics
# This requires properly formatted comparison data

# Step 1: Prepare your comparison data structure
# ranks_df should contain your taxonomic assignments to evaluate
# true_values_df should contain the known correct assignments

# Example structure (replace with your actual data):
# ranks_df <- data.frame(
#   method1 = c("Fungi", "Ascomycota", "Sordariomycetes"),
#   method2 = c("Fungi", "Basidiomycota", "Agaricomycetes")
# )
#
# true_values_df <- data.frame(
#   Kingdom = c("Fungi", "Fungi"),
#   Phylum = c("Glomeromycota", "Glomeromycota")
# )

# Step 2: Compute metrics
# accuracy_results <- tc_metrics_mock(physeq_with_shuffled,
#                                   ranks_df = ranks_df,
#                                   true_values_df = true_values_df,
#                                   fake_taxa = TRUE)
```

## Visualization of NA Patterns

The package includes specialized visualization for understanding missing
taxonomic data:

``` r
# Visualize NA patterns in taxonomic assignments
rainplot_taxo_na(Glom_otu)
```

## Resolving Taxonomic Conflicts

When working with multiple taxonomic databases, conflicts may arise:

``` r
# Resolve taxonomic conflicts between different assignment methods
# This function helps when you have conflicting taxonomic assignments
resolved_physeq <- resolve_taxo_conflict(Glom_otu,
  conflict_rank = "Genus",
  resolve_method = "majority"
)
```

## Advanced Usage Tips

### Workflow for Method Comparison

1.  **Prepare your data**: Start with a phyloseq object containing your
    OTU/ASV data
2.  **Add fake taxa**: Use
    [`add_shuffle_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_shuffle_seq_pq.md)
    and
    [`add_external_seq_pq()`](https://adrientaudiere.github.io/comparpq/reference/add_external_seq_pq.md)
    to create test cases
3.  **Perform taxonomic assignments**: Apply different methods/databases
    to your data
4.  **Compare results**: Use
    [`tc_metrics_mock()`](https://adrientaudiere.github.io/comparpq/reference/tc_metrics_mock.md)
    to quantify accuracy
5.  **Visualize**: Create plots using
    [`bubbles_pq()`](https://adrientaudiere.github.io/comparpq/reference/bubbles_pq.md),
    [`tc_bar()`](https://adrientaudiere.github.io/comparpq/reference/tc_bar.md),
    or
    [`tc_circle()`](https://adrientaudiere.github.io/comparpq/reference/tc_circle.md)

### Data Quality Control

``` r
# Example quality control workflow
qc_physeq <- Glom_otu

# 1. Clean up problematic annotations
qc_physeq <- taxtab_replace_pattern_by_NA(qc_physeq,
  pattern = "unknown|unidentified|uncultured",
  taxonomic_ranks = c("Family", "Genus")
)

# 2. Select relevant taxonomic ranks
qc_physeq <- select_ranks_pq(qc_physeq,
  taxonomic_ranks = c("Kingdom", "Phyla", "Class", "Order", "Family")
)

# 3. Standardize rank names
qc_physeq <- rename_ranks_pq(qc_physeq,
  old_names = c("Kingdom", "Phyla", "Class", "Order", "Family"),
  new_names = c("Kingdom__", "Phyla__", "Class__", "Order__", "Family__")
)

cat("Quality control completed. Final dataset has", phyloseq::ntaxa(qc_physeq), "taxa\n")
#> Quality control completed. Final dataset has 1147 taxa
```

## Getting Help

- **Function documentation**: Use `?function_name` for detailed help on
  any function
- **Package website**: Visit
  <https://adrientaudiere.github.io/comparpq/>
- **Report issues**: <https://github.com/adrientaudiere/comparpq/issues>
- **Related packages**:
  - [{MiscMetabar}](https://adrientaudiere.github.io/MiscMetabar/):
    Extended microbiome analysis
  - [{phyloseq}](https://joey711.github.io/phyloseq/): Core microbiome
    data handling

## Summary

The comparpq package provides a comprehensive toolkit for:

- **Interactive visualization** of phyloseq objects using modern web
  technologies
- **Accuracy assessment** of taxonomic assignment methods using mock
  communities  
- **Data preprocessing** and quality control for microbiome analyses
- **Comparative analysis** between different taxonomic databases or
  assignment approaches

This makes it particularly valuable for researchers working on taxonomic
assignment method development, database comparison, and microbiome data
quality assessment.
