#' Extract variable names from a formula string
#' @param formula_string A character string representing the right part of a formula
#' @return A character vector of variable names
#' @noRd
extract_formula_vars <- function(formula_string) {
  # Remove spaces and split by common formula operators
  clean_formula <- gsub("\\s+", "", formula_string)

  # Split by operators: +, *, :, |
  vars <- unlist(strsplit(clean_formula, "[+*:|]"))

  # Remove empty strings and duplicates
  vars <- unique(vars[vars != ""])

  # Remove numbers (for things like Condition(1))
  vars <- vars[!grepl("^[0-9]+$", vars)]

  # Handle Condition() wrapper

  vars <- gsub("^Condition\\((.+)\\)$", "\\1", vars)

  return(vars)
}

#' Validate formula variables against shared modalities
#'
#' @param x A list_phyloseq object
#' @param formula A formula string
#' @return TRUE if valid, stops with error otherwise
#' @noRd
validate_formula_vars <- function(x, formula) {
  # Extract variables from formula
  vars <- extract_formula_vars(formula)

  # Get all shared sample_data column names
  shared_mod <- x@comparison$shared_sam_data_modalities
  shared_cols <- names(shared_mod)

  # Also get common columns that exist in all phyloseq objects
  # (some columns may exist but not have "shared modalities" if they are numeric)
  sam_data_list <- purrr::map(
    x@phyloseq_list,
    ~ {
      sam <- phyloseq::sample_data(.x, errorIfNULL = FALSE)
      if (is.null(sam)) {
        return(NULL)
      }
      colnames(sam)
    }
  )
  sam_data_list <- purrr::compact(sam_data_list)
  common_cols <- Reduce(intersect, sam_data_list)

  # Check each variable
  missing_vars <- vars[!(vars %in% common_cols)]

  if (length(missing_vars) > 0) {
    stop(
      "Formula variable(s) not found in common sample_data columns: ",
      paste(missing_vars, collapse = ", "),
      "\n",
      "Available common columns: ",
      paste(common_cols, collapse = ", ")
    )
  }

  return(TRUE)
}

#' PERMANOVA analysis on each phyloseq object in a list_phyloseq
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Performs a PERMANOVA analysis using [MiscMetabar::adonis_pq()] on each
#' phyloseq object in a list_phyloseq and returns a summary table with the
#' results.
#'
#' @param x (required) A list_phyloseq object.
#' @param formula (character, required) The right part of a formula for
#'   [vegan::adonis2()]. Variables must be present in the `sample_data` slot
#'   of **all** phyloseq objects. The formula should contain variables that
#'   are in the shared modalities.
#' @param dist_method (character, default "bray") The distance method to use.
#'   See [phyloseq::distance()] for available methods.
#' @param na_remove (logical, default FALSE) If TRUE, samples with NA values
#'   in the formula variables are removed before analysis.
#' @param correction_for_sample_size (logical, default FALSE) If TRUE, adds
#'   library size to the formula following Weiss et al. 2017 recommendations.
#' @param rarefy_nb_seqs (logical, default FALSE) If TRUE, rarefy each sample
#'   before computing distances.
#' @param by (character, default "margin") The `by` argument passed to
#'   [vegan::adonis2()]. Options are "terms", "margin", or NULL.
#' @param verbose (logical, default TRUE) If TRUE, print progress messages.
#' @param ... Additional arguments passed to [MiscMetabar::adonis_pq()].
#'
#' @return A tibble with one row per phyloseq object and the following columns:
#'   \describe{
#'     \item{name}{Name of the phyloseq object}
#'     \item{term}{The term(s) from the formula (one row per term if multiple)}
#'     \item{Df}{Degrees of freedom}
#'     \item{SumOfSqs}{Sum of squares}
#'     \item{R2}{R-squared value}
#'     \item{F}{F statistic}
#'     \item{Pr(>F)}{P-value}
#'     \item{partial_R2}{Partial R-squared (if available)}
#'   }
#'
#' @details
#' This function requires that the list_phyloseq type is NOT `SEPARATE_ANALYSIS`,
#' as the formula must contain variables that are common across all phyloseq
#' objects.
#'
#' The function is a wrapper around [MiscMetabar::adonis_pq()], which itself
#' wraps [vegan::adonis2()].
#'
#' @examples
#' \dontrun{
#' lpq <- list_phyloseq(
#'   list(
#'     fungi = data_fungi,
#'     fungi_clust = postcluster_pq(data_fungi)
#'   ),
#'   same_bioinfo_pipeline = FALSE
#' )
#'
#' # Run PERMANOVA on each phyloseq object
#' results <- adonis_lpq(lpq, formula = "Height+Time", na_remove = TRUE)
#' results
#'
#' # With Jaccard distance
#' results_jaccard <- adonis_lpq(lpq, formula = "Height", dist_method = "jaccard")
#' }
#'
#' @seealso [MiscMetabar::adonis_pq()], [vegan::adonis2()]
#' @export
adonis_lpq <- function(
  x,
  formula,
  dist_method = "bray",
  na_remove = FALSE,
  correction_for_sample_size = FALSE,
  rarefy_nb_seqs = FALSE,
  by = "margin",
  verbose = TRUE,
  ...
) {
  # Validate input
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  # Check comparison type
  if (x@comparison$type_of_comparison == "SEPARATE_ANALYSIS") {
    stop(
      "adonis_lpq cannot be used with SEPARATE_ANALYSIS type list_phyloseq.\n",
      "The formula must contain variables common to all phyloseq objects."
    )
  }

  # Validate formula variables
  validate_formula_vars(x, formula)

  if (verbose) {
    message(
      "Running PERMANOVA on ",
      length(x@phyloseq_list),
      " phyloseq objects"
    )
    message("Formula: ~ ", formula)
    message("Distance method: ", dist_method)
  }

  # Run adonis_pq on each phyloseq object
  results_list <- purrr::imap(x@phyloseq_list, function(pq, name) {
    if (verbose) {
      cli::cli_alert(cli::col_red("Processing: {name}"))
    }

    tryCatch(
      {
        res <- MiscMetabar::adonis_pq(
          physeq = pq,
          formula = formula,
          dist_method = dist_method,
          na_remove = na_remove,
          correction_for_sample_size = correction_for_sample_size,
          rarefy_nb_seqs = rarefy_nb_seqs,
          verbose = FALSE,
          by = by,
          ...
        )

        # Convert result to data frame
        res_df <- as.data.frame(res)

        # Add term names as a column
        res_df$term <- rownames(res_df)

        # Add phyloseq name
        res_df$name <- name

        # Reorder columns
        res_df <- res_df[, c(
          "name",
          "term",
          setdiff(colnames(res_df), c("name", "term"))
        )]

        rownames(res_df) <- NULL

        return(res_df)
      },
      error = function(e) {
        warning("Error processing '", name, "': ", e$message)
        return(NULL)
      }
    )
  })

  # Combine results
  results_list <- purrr::compact(results_list)

  if (length(results_list) == 0) {
    stop("All PERMANOVA analyses failed. Check your data and formula.")
  }

  results <- dplyr::bind_rows(results_list)

  # Convert to tibble
  results <- tibble::as_tibble(results)

  # Add class for potential method dispatch
  class(results) <- c("adonis_lpq_result", class(results))

  return(results)
}


#' Automated model selection for Hill diversity on each phyloseq in a list_phyloseq
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Performs automated model selection and multimodel inference using
#' [MiscMetabar::glmutli_pq()] on each phyloseq object in a list_phyloseq.
#' Returns a summary table with the results from all phyloseq objects.
#'
#' @param x (required) A list_phyloseq object.
#' @param formula (character, required) A model formula for glmulti. Variables
#'   must be present in the `sample_data` slot of **all** phyloseq objects.
#'   Hill numbers (Hill_0, Hill_1, Hill_2) and Abundance are automatically
#'   available.
#' @param fitfunction (character, default "lm") The model fitting function to
#'   use. Options include "lm" for linear models or "glm" for generalized
#'   linear models.
#' @param hill_scales (numeric vector, default c(0, 1, 2)) The q values for
#'   Hill number computation. Defaults to Hill numbers 0 (richness), 1
#'   (Shannon exponential), and 2 (inverse Simpson).
#' @param aic_step (numeric, default 2) The AIC score threshold for model
#'   selection. Models within this threshold from the best model are included.
#' @param confsetsize (integer, default 100) The number of models to return in
#'   the confidence set.
#' @param plotty (logical, default FALSE) If TRUE, display IC profile during
#'   glmulti search.
#' @param level (integer, default 1) Model complexity level. 1 for main effects
#'   only, 2 for pairwise interactions.
#' @param method (character, default "h") The search method for glmulti.
#'   Options: "h" (exhaustive), "g" (genetic algorithm), "l" (branch-and-bound),
#'   "d" (summary only).
#' @param crit (character, default "aicc") Information criterion for model
#'   selection. Options include "aic", "aicc" (small-sample corrected AIC),
#'   "bic".
#' @param verbose (logical, default TRUE) If TRUE, print progress messages.
#' @param ... Additional arguments passed to [MiscMetabar::glmutli_pq()].
#'
#' @return A tibble with the combined results from all phyloseq objects,
#'   containing the following columns:
#'   \describe{
#'     \item{name}{Name of the phyloseq object}
#'     \item{variable}{The variable name from the model}
#'     \item{estimates}{The model coefficient estimate}
#'     \item{unconditional_interval}{Confidence interval from model averaging}
#'     \item{nb_model}{Number of models containing this variable}
#'     \item{importance}{Relative importance of the variable (sum of Akaike
#'       weights)}
#'     \item{alpha}{Significance level}
#'   }
#'
#' @details
#' This function requires that the list_phyloseq type is NOT `SEPARATE_ANALYSIS`,
#' as the formula must contain variables that are common across all phyloseq
#' objects.
#'
#' The function wraps [MiscMetabar::glmutli_pq()], which itself wraps the
#' \pkg{glmulti} package for automated model selection. For each phyloseq
#' object, Hill diversity indices are computed and used as response variables
#' in the model selection process.
#'
#' @examples
#' lpq <- list_phyloseq(
#'   list(
#'     fungi = data_fungi,
#'     fungi_clust = postcluster_pq(data_fungi)
#'   ),
#'   same_bioinfo_pipeline = FALSE
#' )
#'
#' results <- glmulti_lpq(lpq, formula = "Hill_0 ~ Height + Time")
#' results
#'
#' # With interactions
#' results_int <- glmulti_lpq(lpq, formula = "Hill_1 ~ Height * Time", level = 2)
#'
#' @seealso [MiscMetabar::glmutli_pq()], [MiscMetabar::hill_pq()]
#' @export
glmulti_lpq <- function(
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
) {
  # Validate input
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  # Check comparison type
  if (x@comparison$type_of_comparison == "SEPARATE_ANALYSIS") {
    stop(
      "glmulti_lpq cannot be used with SEPARATE_ANALYSIS type list_phyloseq.\n",
      "The formula must contain variables common to all phyloseq objects."
    )
  }

  # Extract the right-hand side variables from the formula for validation
  formula_parts <- strsplit(formula, "~")[[1]]
  if (length(formula_parts) == 2) {
    rhs <- trimws(formula_parts[2])
    validate_formula_vars(x, rhs)
  }

  if (verbose) {
    message("Running glmulti on ", length(x@phyloseq_list), " phyloseq objects")
    message("Formula: ", formula)
    message("Hill scales: ", paste(hill_scales, collapse = ", "))
  }

  # Run glmutli_pq on each phyloseq object
  results_list <- purrr::imap(x@phyloseq_list, function(pq, name) {
    if (verbose) {
      cli::cli_alert(cli::col_red("Processing: {name}"))
    }

    tryCatch(
      {
        res <- MiscMetabar::glmutli_pq(
          physeq = pq,
          formula = formula,
          fitfunction = fitfunction,
          hill_scales = hill_scales,
          aic_step = aic_step,
          confsetsize = confsetsize,
          plotty = plotty,
          level = level,
          method = method,
          crit = crit,
          ...
        )

        # Add phyloseq name as first column
        res$name <- name
        res <- res[, c("name", setdiff(colnames(res), "name"))]

        return(res)
      },
      error = function(e) {
        warning("Error processing '", name, "': ", e$message)
        return(NULL)
      }
    )
  })

  # Combine results
  results_list <- purrr::compact(results_list)

  if (length(results_list) == 0) {
    stop("All glmulti analyses failed. Check your data and formula.")
  }

  results <- dplyr::bind_rows(results_list)

  # Convert to tibble
  results <- tibble::as_tibble(results)

  # Add class for potential method dispatch
  class(results) <- c("glmulti_lpq_result", class(results))

  return(results)
}


#' Suffix taxa names to make them distinguishable across phyloseq objects
#'
#' When phyloseq objects share no common taxa, appends `_<name>` to each
#' taxon name so that taxa from different objects are visually distinguishable
#' in combined result tables.
#'
#' @param phyloseq_list A named list of phyloseq objects
#' @param n_common_taxa Integer, number of common taxa
#' @return A named list of phyloseq objects (possibly with suffixed taxa names)
#' @noRd
maybe_suffix_taxa <- function(phyloseq_list, n_common_taxa) {
  if (n_common_taxa > 0) {
    return(phyloseq_list)
  }
  purrr::imap(phyloseq_list, function(pq, name) {
    new_names <- paste0(phyloseq::taxa_names(pq), "_", name)
    phyloseq::taxa_names(pq) <- new_names
    pq
  })
}


#' Validate that a factor column exists in common sample_data columns
#'
#' @param x A list_phyloseq object
#' @param fact Character, the column name to check
#' @param func_name Character, the calling function name (for error messages)
#' @noRd
validate_fact_column <- function(x, fact, func_name) {
  if (x@comparison$type_of_comparison == "SEPARATE_ANALYSIS") {
    stop(
      func_name, " cannot be used with SEPARATE_ANALYSIS type list_phyloseq.\n",
      "The factor must be a variable common to all phyloseq objects."
    )
  }

  sam_data_list <- purrr::map(
    x@phyloseq_list,
    ~ colnames(phyloseq::sample_data(.x, errorIfNULL = FALSE))
  )
  sam_data_list <- purrr::compact(sam_data_list)
  common_cols <- Reduce(intersect, sam_data_list)

  if (!fact %in% common_cols) {
    stop(
      "Factor '", fact, "' not found in common sample_data columns.\n",
      "Available common columns: ", paste(common_cols, collapse = ", ")
    )
  }
}


#' ANCOM-BC analysis on each phyloseq object in a list_phyloseq
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Performs ANCOM-BC differential abundance analysis using
#' [MiscMetabar::ancombc_pq()] on each phyloseq object in a list_phyloseq
#' and returns a combined result table.
#'
#' @param x (required) A list_phyloseq object.
#' @param fact (character, required) The name of a column in `sample_data`
#'   to use as the grouping factor. Must be present in **all** phyloseq
#'   objects.
#' @param levels_fact (character vector, default NULL) Levels of the factor
#'   to include. If NULL, all levels are used.
#' @param tax_level (character, default "Class") Taxonomic level for
#'   agglomeration before analysis.
#' @param verbose (logical, default TRUE) If TRUE, print progress messages.
#' @param ... Additional arguments passed to [MiscMetabar::ancombc_pq()].
#'
#' @return A tibble with the combined ANCOM-BC results from all phyloseq
#'   objects. Contains the `$res` data frame from each ANCOM-BC run, with
#'   an additional `name` column identifying the source phyloseq object.
#'
#' @details
#' This function requires that the list_phyloseq type is NOT
#' `SEPARATE_ANALYSIS`, as the factor must be common across all phyloseq
#' objects.
#'
#' When no common taxa exist across the phyloseq objects, taxa names are
#' suffixed with the phyloseq object name to make them distinguishable.
#'
#' @examples
#' data_fungi_high <- multiply_counts_pq(data_fungi, "Height", "High", 2)
#'
#' lpq <- list_phyloseq(
#'   list(
#'     fungi = data_fungi,
#'     fungi_height = data_fungi_high
#'   )
#' )
#'
#' results <- ancombc_lpq(lpq, fact = "Height")
#' results
#'
#' results |>
#'   filter(diff_HeightLow) |>
#'   ggplot(aes(y = taxon, x = lfc_HeightLow, color = name, shape = q_HeightLow < 0.1)) +
#'   geom_point(size = 4) +
#'   geom_vline(xintercept = 0) +
#'   theme_minimal()
#'
#' results_Genus <- ancombc_lpq(lpq, fact = "Height", tax_level = "Genus")
#'
#' results_Genus |>
#'   filter(diff_HeightLow) |>
#'   ggplot(aes(y = taxon, x = lfc_HeightLow, color = name, shape = q_HeightLow < 0.1)) +
#'   geom_point(size = 4) +
#'   geom_vline(xintercept = 0) +
#'   theme_minimal()
#'
#' @seealso [MiscMetabar::ancombc_pq()], [aldex_lpq()], [multipatt_lpq()]
#' @export
ancombc_lpq <- function(
  x,
  fact,
  levels_fact = NULL,
  tax_level = NULL,
  verbose = TRUE,
  ...
) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))
  validate_fact_column(x, fact, "ancombc_lpq")

  pq_list <- maybe_suffix_taxa(x@phyloseq_list, x@comparison$n_common_taxa)

  if (verbose) {
    message(
      "Running ANCOM-BC on ", length(pq_list), " phyloseq objects"
    )
    message("Factor: ", fact, " | Tax level: ", tax_level)
  }

  results_list <- purrr::imap(pq_list, function(pq, name) {
    if (verbose) cli::cli_alert(cli::col_red("Processing: {name}"))

    tryCatch(
      {
        res <- MiscMetabar::ancombc_pq(
          physeq = pq,
          fact = fact,
          levels_fact = levels_fact,
          tax_level = tax_level,
          ...
        )

        res_df <- res$res
        res_df$name <- name
        res_df <- res_df[, c("name", setdiff(colnames(res_df), "name"))]
        return(res_df)
      },
      error = function(e) {
        warning("Error processing '", name, "': ", e$message)
        return(NULL)
      }
    )
  })

  results_list <- purrr::compact(results_list)

  if (length(results_list) == 0) {
    stop("All ANCOM-BC analyses failed. Check your data and factor.")
  }

  results <- dplyr::bind_rows(results_list)
  results <- tibble::as_tibble(results)
  class(results) <- c("ancombc_lpq_result", class(results))
  return(results)
}


#' ALDEx2 analysis on each phyloseq object in a list_phyloseq
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Performs ALDEx2 differential abundance analysis using
#' [MiscMetabar::aldex_pq()] on each phyloseq object in a list_phyloseq
#' and returns a combined result table.
#'
#' @param x (required) A list_phyloseq object.
#' @param bifactor (character, required) The name of a dichotomous column
#'   in `sample_data` to use as the grouping factor. Must be present in
#'   **all** phyloseq objects.
#' @param modalities (character vector, default NULL) Two modalities of
#'   `bifactor` to compare. If NULL, uses the two levels present in the
#'   data.
#' @param gamma (numeric, default 0.5) The gamma parameter for ALDEx2.
#' @param verbose (logical, default TRUE) If TRUE, print progress messages.
#' @param ... Additional arguments passed to [MiscMetabar::aldex_pq()].
#'
#' @return A tibble with the combined ALDEx2 results from all phyloseq
#'   objects. Each row corresponds to one taxon in one phyloseq object,
#'   with columns from ALDEx2 output plus `taxon` (from rownames) and
#'   `name` (identifying the source phyloseq object).
#'
#' @details
#' This function requires that the list_phyloseq type is NOT
#' `SEPARATE_ANALYSIS`, as the bifactor must be common across all phyloseq
#' objects.
#'
#' When no common taxa exist across the phyloseq objects, taxa names are
#' suffixed with the phyloseq object name to make them distinguishable.
#'
#' @examples
#' \dontrun{
#' lpq <- list_phyloseq(
#'   list(
#'     fungi = data_fungi,
#'     fungi_clust = postcluster_pq(data_fungi)
#'   ),
#'   same_bioinfo_pipeline = FALSE
#' )
#'
#' results <- aldex_lpq(lpq,
#'   bifactor = "Height",
#'   modalities = c("Low", "High")
#' )
#' results
#'
#' ALDEx2::aldex.plot(filter(results, name == "fungi"), type = "volcano")
#' ALDEx2::aldex.plot(filter(results, name == "fungi_clust"), type = "volcano")
#'
#' ggplot(results, aes(y = taxon, x = effect, col = wi.eBH)) +
#'   geom_point()
#' }
#'
#' @seealso [MiscMetabar::aldex_pq()], [ancombc_lpq()], [multipatt_lpq()]
#' @export
aldex_lpq <- function(
  x,
  bifactor,
  modalities = NULL,
  gamma = 0.5,
  verbose = TRUE,
  ...
) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))
  validate_fact_column(x, bifactor, "aldex_lpq")

  pq_list <- maybe_suffix_taxa(x@phyloseq_list, x@comparison$n_common_taxa)

  if (verbose) {
    message("Running ALDEx2 on ", length(pq_list), " phyloseq objects")
    message("Bifactor: ", bifactor)
  }

  results_list <- purrr::imap(pq_list, function(pq, name) {
    if (verbose) cli::cli_alert(cli::col_red("Processing: {name}"))

    tryCatch(
      {
        res <- MiscMetabar::aldex_pq(
          physeq = pq,
          bifactor = bifactor,
          modalities = modalities,
          gamma = gamma,
          ...
        )

        res_df <- tibble::as_tibble(res, rownames = "taxon")
        res_df$name <- name
        res_df <- res_df[, c("name", setdiff(colnames(res_df), "name"))]
        return(res_df)
      },
      error = function(e) {
        warning("Error processing '", name, "': ", e$message)
        return(NULL)
      }
    )
  })

  results_list <- purrr::compact(results_list)

  if (length(results_list) == 0) {
    stop("All ALDEx2 analyses failed. Check your data and bifactor.")
  }

  results <- dplyr::bind_rows(results_list)
  results <- tibble::as_tibble(results)
  class(results) <- c("aldex_lpq_result", class(results))
  return(results)
}


#' Indicator species analysis on each phyloseq object in a list_phyloseq
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Performs indicator species analysis using [indicspecies::multipatt()] on
#' each phyloseq object in a list_phyloseq and returns a combined result
#' table of significant indicator taxa.
#'
#' @param x (required) A list_phyloseq object.
#' @param fact (character, required) The name of a column in `sample_data`
#'   to use as the grouping factor. Must be present in **all** phyloseq
#'   objects.
#' @param p_adjust_method (character, default "BH") The p-value adjustment
#'   method. See [stats::p.adjust()] for available methods.
#' @param pval (numeric, default 0.05) The significance threshold for
#'   adjusted p-values.
#' @param control (list, default `permute::how(nperm = 999)`) Permutation
#'   control settings for the permutation test.
#' @param verbose (logical, default TRUE) If TRUE, print progress messages.
#' @param ... Additional arguments passed to [indicspecies::multipatt()].
#'
#' @return A tibble with the combined significant indicator taxa from all
#'   phyloseq objects. Contains columns from `multipatt()$sign` output
#'   plus `taxon` (taxon name), `p.adj` (adjusted p-value), and `name`
#'   (identifying the source phyloseq object). Only taxa with
#'   `p.adj < pval` are included.
#'
#' @details
#' This function requires that the list_phyloseq type is NOT
#' `SEPARATE_ANALYSIS`, as the factor must be common across all phyloseq
#' objects.
#'
#' Unlike [MiscMetabar::multipatt_pq()] which returns a plot, this function
#' returns the underlying data as a tibble, making it easier to compare
#' results across phyloseq objects.
#'
#' When no common taxa exist across the phyloseq objects, taxa names are
#' suffixed with the phyloseq object name to make them distinguishable.
#'
#' @examples
#' \dontrun{
#' lpq <- list_phyloseq(
#'   list(
#'     fungi = data_fungi,
#'     fungi_clust = postcluster_pq(data_fungi)
#'   ),
#'   same_bioinfo_pipeline = FALSE
#' )
#'
#' results <- multipatt_lpq(lpq, fact = "Height")
#' results
#' }
#'
#' @seealso [indicspecies::multipatt()], [MiscMetabar::multipatt_pq()],
#'   [ancombc_lpq()], [aldex_lpq()]
#' @export
multipatt_lpq <- function(
  x,
  fact,
  p_adjust_method = "BH",
  pval = 0.05,
  control = permute::how(nperm = 999),
  verbose = TRUE,
  ...
) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))
  validate_fact_column(x, fact, "multipatt_lpq")

  pq_list <- maybe_suffix_taxa(x@phyloseq_list, x@comparison$n_common_taxa)

  if (verbose) {
    message("Running multipatt on ", length(pq_list), " phyloseq objects")
    message("Factor: ", fact)
  }

  results_list <- purrr::imap(pq_list, function(pq, name) {
    if (verbose) cli::cli_alert(cli::col_red("Processing: {name}"))

    tryCatch(
      {
        pq <- MiscMetabar::taxa_as_columns(pq)
        res <- indicspecies::multipatt(
          as.matrix(phyloseq::otu_table(pq)),
          phyloseq::sample_data(pq)[[fact]],
          control = control,
          ...
        )

        res_df <- res$sign
        res_df$p.adj <- stats::p.adjust(
          res_df$p.value,
          method = p_adjust_method
        )
        res_df$taxon <- rownames(res_df)
        res_df <- res_df[res_df$p.adj < pval, , drop = FALSE]

        if (nrow(res_df) == 0) {
          if (verbose) {
            message("    No significant indicators found for '", name, "'")
          }
          return(NULL)
        }

        res_df <- tibble::as_tibble(res_df)
        res_df$name <- name
        res_df <- res_df[, c("name", setdiff(colnames(res_df), "name"))]
        return(res_df)
      },
      error = function(e) {
        warning("Error processing '", name, "': ", e$message)
        return(NULL)
      }
    )
  })

  results_list <- purrr::compact(results_list)

  if (length(results_list) == 0) {
    stop(
      "No significant indicator taxa found in any phyloseq object. ",
      "Consider adjusting pval threshold or checking your data."
    )
  }

  results <- dplyr::bind_rows(results_list)
  results <- tibble::as_tibble(results)
  class(results) <- c("multipatt_lpq_result", class(results))
  return(results)
}
