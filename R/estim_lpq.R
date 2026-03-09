#' Estimation statistics for categorical comparisons on a list_phyloseq
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#' alt="lifecycle-experimental"></a>
#'
#' @description
#' Applies [estim_diff_pq()] to each phyloseq object in a list_phyloseq
#' and combines the results. This allows comparing estimation statistics
#' across different bioinformatic pipelines or parameter settings.
#'
#' @param x (required) A list_phyloseq object.
#' @param fact (character, required) The name of a categorical column in
#'   `sample_data` to use as the grouping factor. Must be present in
#'   **all** phyloseq objects.
#' @param verbose (logical, default TRUE) If TRUE, print progress
#'   messages.
#' @param ... Additional arguments passed to [estim_diff_pq()].
#'
#' @return A list of class `"estim_diff_lpq_result"` with components:
#'   \describe{
#'     \item{results}{A named list of `estim_diff_pq_result` objects
#'       (one per phyloseq)}
#'     \item{summary}{A tibble combining all summaries with an
#'       additional `name` column identifying the source phyloseq}
#'   }
#'
#' @examples
#' lpq <- list_phyloseq(
#'   list(
#'     fungi = data_fungi,
#'     fungi_clust = postcluster_pq(data_fungi),
#'     fungi_rarefy = rarefy_even_depth(data_fungi),
#'     fungi_with_less_otu_in_High =  multiply_counts_pq(data_fungi,
#'       fact = "Height", prop=0.8,
#'       conditions = "High",
#'       multipliers = 0)
#'   ),
#'   same_bioinfo_pipeline = FALSE
#' )
#'
#' results <- estim_diff_lpq(lpq, fact = "Height")
#' results$summary
#'
#' # Plot results for two phyloseq objects
#'  ggplot(results$summary, aes(x = metric, y = effect_size, color=name)) +
#'    facet_wrap(~comparison) +
#'    geom_point(position=position_dodge(width=0.5)) +
#'    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position=position_dodge(width=0.5 ))
#' @seealso [estim_diff_pq()], [estim_cor_lpq()], [adonis_lpq()]
#' @export
estim_diff_lpq <- function(x, fact, ..., verbose = TRUE) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))
  validate_fact_column(x, fact, "estim_diff_lpq")

  if (verbose) {
    message(
      "Running estimation statistics (categorical) on ",
      length(x@phyloseq_list),
      " phyloseq objects"
    )
    message("Factor: ", fact)
  }

  results_list <- purrr::imap(x@phyloseq_list, function(pq, name) {
    if (verbose) {
      cli::cli_alert(cli::col_red("Processing: {name}"))
    }

    tryCatch(
      {
        estim_diff_pq(physeq = pq, fact = fact, ...)
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
      "All estimation statistics analyses failed. Check your data and factor."
    )
  }

  # Combine summaries
  summary_list <- purrr::imap(results_list, function(res, name) {
    tbl <- res$summary
    tbl$name <- name
    tbl <- tbl[, c("name", setdiff(colnames(tbl), "name"))]
    tbl
  })
  summary_tbl <- dplyr::bind_rows(summary_list)

  result <- list(
    results = results_list,
    summary = summary_tbl
  )
  class(result) <- "estim_diff_lpq_result"
  result
}


#' Estimation statistics for numeric correlation on a list_phyloseq
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#' alt="lifecycle-experimental"></a>
#'
#' @description
#' Applies [estim_cor_pq()] to each phyloseq object in a list_phyloseq
#' and combines the results.
#'
#' @param x (required) A list_phyloseq object.
#' @param variable (character, required) The name of a numeric column
#'   in `sample_data`. Must be present in **all** phyloseq objects.
#' @param verbose (logical, default TRUE) If TRUE, print progress
#'   messages.
#' @param ... Additional arguments passed to [estim_cor_pq()].
#'
#' @return A list of class `"estim_cor_lpq_result"` with components:
#'   \describe{
#'     \item{results}{A named list of `estim_cor_pq_result` objects
#'       (one per phyloseq)}
#'     \item{correlations}{A tibble combining all correlations with an
#'       additional `name` column}
#'     \item{regressions}{A tibble combining all regressions with an
#'       additional `name` column}
#'   }
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
#' # Assuming a numeric variable exists in sample_data
#' results <- estim_cor_lpq(lpq, variable = "lib_size")
#' results$correlations
#' }
#'
#' @seealso [estim_cor_pq()], [estim_diff_lpq()]
#' @export
estim_cor_lpq <- function(x, variable, ..., verbose = TRUE) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  # Validate variable exists in all phyloseq objects
  # Reuse validate_fact_column logic (works for any column name)
  validate_fact_column(x, variable, "estim_cor_lpq")

  if (verbose) {
    message(
      "Running estimation statistics (correlation) on ",
      length(x@phyloseq_list),
      " phyloseq objects"
    )
    message("Variable: ", variable)
  }

  results_list <- purrr::imap(x@phyloseq_list, function(pq, name) {
    if (verbose) {
      cli::cli_alert(cli::col_red("Processing: {name}"))
    }

    tryCatch(
      {
        estim_cor_pq(physeq = pq, variable = variable, ...)
      },
      error = function(e) {
        warning("Error processing '", name, "': ", e$message)
        return(NULL)
      }
    )
  })

  results_list <- purrr::compact(results_list)

  if (length(results_list) == 0) {
    stop("All correlation analyses failed. Check your data and variable.")
  }

  # Combine correlations
  cor_list <- purrr::imap(results_list, function(res, name) {
    tbl <- res$correlations
    tbl$name <- name
    tbl <- tbl[, c("name", setdiff(colnames(tbl), "name"))]
    tbl
  })
  cor_tbl <- dplyr::bind_rows(cor_list)

  # Combine regressions
  reg_list <- purrr::imap(results_list, function(res, name) {
    tbl <- res$regressions
    tbl$name <- name
    tbl <- tbl[, c("name", setdiff(colnames(tbl), "name"))]
    tbl
  })
  reg_tbl <- dplyr::bind_rows(reg_list)

  result <- list(
    results = results_list,
    correlations = cor_tbl,
    regressions = reg_tbl
  )
  class(result) <- "estim_cor_lpq_result"
  result
}


# ==============================================================================
# Print methods
# ==============================================================================

#' @export
print.estim_diff_lpq_result <- function(x, ...) {
  cat("Estimation statistics: categorical comparison (list_phyloseq)\n")
  cat("Phyloseq objects:", length(x$results), "\n\n")
  print(x$summary)
  cat(
    "\nNote: P-values (permutation test, Welch t-test, Mann-Whitney U)",
    "are provided\nfor legacy purposes only. Focus on the effect size",
    "and its confidence interval\nfor interpretation.",
    "See Ho et al. (2019) Nature Methods 16:565-566.\n"
  )
  invisible(x)
}

#' @export
print.estim_cor_lpq_result <- function(x, ...) {
  cat("Estimation statistics: numeric correlation (list_phyloseq)\n")
  cat("Phyloseq objects:", length(x$results), "\n\n")
  cat("Correlations:\n")
  print(x$correlations)
  cat("\nRegressions:\n")
  print(x$regressions)
  cat(
    "\nNote: P-values (from cor.test) are provided for legacy purposes",
    "only.\nFocus on the effect size (correlation estimate, regression",
    "slope) and its\nbootstrap confidence interval for interpretation.",
    "See Ho et al. (2019)\nNature Methods 16:565-566.\n"
  )
  invisible(x)
}
