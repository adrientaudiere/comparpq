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
  q = c(0, 1, 2),
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

- q:

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
if (FALSE) { # \dontrun{
lpq <- list_phyloseq(
  list(
    fungi = data_fungi,
    fungi_clust = postcluster_pq(data_fungi)
  ),
  same_bioinfo_pipeline = FALSE
)

results <- glmulti_lpq(lpq, formula = "Hill_0 ~ Height + Time")
results

# With interactions
results_int <- glmulti_lpq(lpq, formula = "Hill_1 ~ Height * Time", level = 2)
} # }
```
