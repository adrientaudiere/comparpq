# Automated model selection for Hill diversity on each phyloseq in a list_phyloseq

Performs automated model selection and multimodel inference using
[`MiscMetabar::glmutli_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/glmutli_pq.html)
on each phyloseq object in a list_phyloseq. Returns a summary table with
the results from all phyloseq objects.

## Usage

``` r
glmulti_lpq(
  x,
  formula,
  fitfunction = "lm",
  hill_scales = c(0, 1, 2),
  aic_step = 2,
  confsetsize = 100,
  plotty = FALSE,
  level = 1,
  method = "h",
  crit = "aicc",
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  (required) A list_phyloseq object.

- formula:

  (character, required) A model formula for glmulti. Variables must be
  present in the `sample_data` slot of **all** phyloseq objects. Hill
  numbers (Hill_0, Hill_1, Hill_2) and Abundance are automatically
  available.

- fitfunction:

  (character, default "lm") The model fitting function to use. Options
  include "lm" for linear models or "glm" for generalized linear models.

- hill_scales:

  (numeric vector, default c(0, 1, 2)) The q values for Hill number
  computation. Defaults to Hill numbers 0 (richness), 1 (Shannon
  exponential), and 2 (inverse Simpson).

- aic_step:

  (numeric, default 2) The AIC score threshold for model selection.
  Models within this threshold from the best model are included.

- confsetsize:

  (integer, default 100) The number of models to return in the
  confidence set.

- plotty:

  (logical, default FALSE) If TRUE, display IC profile during glmulti
  search.

- level:

  (integer, default 1) Model complexity level. 1 for main effects only,
  2 for pairwise interactions.

- method:

  (character, default "h") The search method for glmulti. Options: "h"
  (exhaustive), "g" (genetic algorithm), "l" (branch-and-bound), "d"
  (summary only).

- crit:

  (character, default "aicc") Information criterion for model selection.
  Options include "aic", "aicc" (small-sample corrected AIC), "bic".

- verbose:

  (logical, default TRUE) If TRUE, print progress messages.

- ...:

  Additional arguments passed to
  [`MiscMetabar::glmutli_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/glmutli_pq.html).

## Value

A tibble with the combined results from all phyloseq objects, containing
the following columns:

- name:

  Name of the phyloseq object

- variable:

  The variable name from the model

- estimates:

  The model coefficient estimate

- unconditional_interval:

  Confidence interval from model averaging

- nb_model:

  Number of models containing this variable

- importance:

  Relative importance of the variable (sum of Akaike weights)

- alpha:

  Significance level

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function requires that the list_phyloseq type is NOT
`SEPARATE_ANALYSIS`, as the formula must contain variables that are
common across all phyloseq objects.

The function wraps
[`MiscMetabar::glmutli_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/glmutli_pq.html),
which itself wraps the glmulti package for automated model selection.
For each phyloseq object, Hill diversity indices are computed and used
as response variables in the model selection process.

## See also

[`MiscMetabar::glmutli_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/glmutli_pq.html),
[`MiscMetabar::hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.html)

## Examples

``` r
lpq <- list_phyloseq(
  list(
    fungi = data_fungi,
    fungi_clust = postcluster_pq(data_fungi)
  ),
  same_bioinfo_pipeline = FALSE
)
#> Partitioning sequences by 6-mer similarity:
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
#> 
#> Time difference of 0.34 secs
#> 
#> Sorting by relatedness within 779 groups:
#> iteration 1 of up to 15 (100.0% stability) iteration 2 of up to 15 (100.0% stability) iteration 3 of up to 15 (100.0% stability) iteration 4 of up to 15 (100.0% stability) iteration 4 of up to 15 (100.0% stability) 
#> 
#> Time difference of 0.27 secs
#> 
#> Clustering sequences by 9-mer similarity:
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
#> 
#> Time difference of 3.42 secs
#> 
#> Clusters via relatedness sorting: 100% (0% exclusively)
#> Clusters via rare 6-mers: 100% (0% exclusively)
#> Estimated clustering effectiveness: 100%
#> 

results <- glmulti_lpq(lpq, formula = "Hill_0 ~ Height + Time")
#> Running glmulti on 2 phyloseq objects
#> Formula: Hill_0 ~ Height + Time
#> Hill scales: 0, 1, 2
#>   Processing: fungi
#> Taxa are now in rows.
#> Joining with `by = join_by(Sample)`
#> Warning: Error processing 'fungi': .onLoad failed in loadNamespace() for 'rJava', details:
#>   call: dyn.load(file, DLLpath = DLLpath, ...)
#>   error: unable to load shared object '/home/adrien/R/x86_64-pc-linux-gnu-library/4.5/rJava/libs/rJava.so':
#>   libjvm.so: cannot open shared object file: No such file or directory
#>   Processing: fungi_clust
#> Taxa are now in rows.
#> Joining with `by = join_by(Sample)`
#> Warning: Error processing 'fungi_clust': .onLoad failed in loadNamespace() for 'rJava', details:
#>   call: dyn.load(file, DLLpath = DLLpath, ...)
#>   error: unable to load shared object '/home/adrien/R/x86_64-pc-linux-gnu-library/4.5/rJava/libs/rJava.so':
#>   libjvm.so: cannot open shared object file: No such file or directory
#> Error in glmulti_lpq(lpq, formula = "Hill_0 ~ Height + Time"): All glmulti analyses failed. Check your data and formula.
results
#> Error: object 'results' not found

# With interactions
results_int <- glmulti_lpq(lpq, formula = "Hill_1 ~ Height * Time", level = 2)
#> Running glmulti on 2 phyloseq objects
#> Formula: Hill_1 ~ Height * Time
#> Hill scales: 0, 1, 2
#>   Processing: fungi
#> Taxa are now in rows.
#> Joining with `by = join_by(Sample)`
#> Warning: Error processing 'fungi': .onLoad failed in loadNamespace() for 'rJava', details:
#>   call: dyn.load(file, DLLpath = DLLpath, ...)
#>   error: unable to load shared object '/home/adrien/R/x86_64-pc-linux-gnu-library/4.5/rJava/libs/rJava.so':
#>   libjvm.so: cannot open shared object file: No such file or directory
#>   Processing: fungi_clust
#> Taxa are now in rows.
#> Joining with `by = join_by(Sample)`
#> Warning: Error processing 'fungi_clust': .onLoad failed in loadNamespace() for 'rJava', details:
#>   call: dyn.load(file, DLLpath = DLLpath, ...)
#>   error: unable to load shared object '/home/adrien/R/x86_64-pc-linux-gnu-library/4.5/rJava/libs/rJava.so':
#>   libjvm.so: cannot open shared object file: No such file or directory
#> Error in glmulti_lpq(lpq, formula = "Hill_1 ~ Height * Time", level = 2): All glmulti analyses failed. Check your data and formula.
```
