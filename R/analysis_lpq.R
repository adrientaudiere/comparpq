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
      message("  Processing: ", name)
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
      message("  Processing: ", name)
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
