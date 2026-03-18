# ==============================================================================
# Internal helpers for estimation statistics
# ==============================================================================

#' Compute Hill numbers per sample
#'
#' @param physeq A phyloseq object
#' @param q Numeric vector of Hill number orders (q values)
#' @return A data.frame with Hill columns and all sample_data columns
#' @noRd
hill_samples_pq <- function(physeq, q = c(0, 1, 2)) {
  otu <- as.data.frame(phyloseq::otu_table(taxa_as_columns(physeq)))

  hill_df <- MiscMetabar::divent_hill_matrix_pq(otu, q = q)
  colnames(hill_df) <- paste0("Hill_", q)
  sam <- as.data.frame(phyloseq::sample_data(physeq))
  cbind(hill_df, sam)
}


#' Compute diversity metrics per sample
#'
#' Delegates to hill_samples_pq or uses a custom function.
#'
#' @param physeq A phyloseq object
#' @param q Numeric vector of Hill number orders
#' @param custom_fn A function taking a phyloseq and returning a named numeric
#'   vector (names = sample names) or a data.frame with one row per sample
#' @return A data.frame with metric columns and sample_data columns
#' @noRd
diversity_samples_pq <- function(
  physeq,
  q = c(0, 1, 2),
  custom_fn = NULL
) {
  if (is.null(custom_fn)) {
    return(hill_samples_pq(physeq, q = q))
  }

  result <- custom_fn(physeq)

  sam <- as.data.frame(phyloseq::sample_data(physeq))
  samp_names <- phyloseq::sample_names(physeq)

  if (is.numeric(result) && !is.null(names(result))) {
    # Named numeric vector
    if (!all(samp_names %in% names(result))) {
      stop(
        "custom_fn must return values for all samples. Missing: ",
        paste(setdiff(samp_names, names(result)), collapse = ", ")
      )
    }
    df <- data.frame(
      custom_metric = result[samp_names],
      row.names = samp_names
    )
    return(cbind(df, sam))
  }

  if (is.data.frame(result)) {
    if (nrow(result) != length(samp_names)) {
      stop(
        "custom_fn data.frame must have one row per sample (",
        length(samp_names),
        " expected, got ",
        nrow(result),
        ")"
      )
    }
    if (
      !is.null(rownames(result)) &&
        all(samp_names %in% rownames(result))
    ) {
      result <- result[samp_names, , drop = FALSE]
    }
    return(cbind(result, sam))
  }

  stop("custom_fn must return a named numeric vector or a data.frame")
}


#' Bootstrap percentile CI for a correlation coefficient
#'
#' @param x,y Numeric vectors
#' @param method Correlation method ("pearson", "spearman", "kendall")
#' @param resamples Number of bootstrap resamples
#' @param ci Confidence level (0-100)
#' @return A list with estimate, ci_lower, ci_upper, boot_distribution
#' @noRd
bootstrap_cor <- function(
  x,
  y,
  method = "pearson",
  resamples = 5000,
  ci = 95,
  use = "complete.obs"
) {
  n <- length(x)
  estimate <- stats::cor(x, y, method = method, use = use)

  alpha <- (100 - ci) / 200
  boot_vals <- vapply(
    seq_len(resamples),
    \(i) {
      idx <- sample.int(n, replace = TRUE)
      stats::cor(x[idx], y[idx], method = method, use = use)
    },
    numeric(1)
  )

  list(
    estimate = estimate,
    ci_lower = stats::quantile(boot_vals, alpha, names = FALSE, na.rm = TRUE),
    ci_upper = stats::quantile(
      boot_vals,
      1 - alpha,
      names = FALSE,
      na.rm = TRUE
    ),
    boot_distribution = boot_vals
  )
}


#' Bootstrap percentile CI for linear regression coefficients
#'
#' @param x,y Numeric vectors
#' @param resamples Number of bootstrap resamples
#' @param ci Confidence level (0-100)
#' @return A list with intercept, slope, slope_ci_lower, slope_ci_upper,
#'   boot_slopes
#' @noRd
bootstrap_lm <- function(x, y, resamples = 5000, ci = 95) {
  fit <- stats::lm(y ~ x)
  intercept <- stats::coef(fit)[1]
  slope <- stats::coef(fit)[2]

  n <- length(x)
  alpha <- (100 - ci) / 200
  boot_slopes <- vapply(
    seq_len(resamples),
    \(i) {
      idx <- sample.int(n, replace = TRUE)
      stats::coef(stats::lm(y[idx] ~ x[idx]))[2]
    },
    numeric(1)
  )

  list(
    intercept = unname(intercept),
    slope = unname(slope),
    slope_ci_lower = stats::quantile(
      boot_slopes,
      alpha,
      names = FALSE,
      na.rm = TRUE
    ),
    slope_ci_upper = stats::quantile(
      boot_slopes,
      1 - alpha,
      names = FALSE,
      na.rm = TRUE
    ),
    boot_slopes = boot_slopes
  )
}


# ==============================================================================
# Exported functions
# ==============================================================================

#' Estimation statistics for categorical comparisons on a phyloseq object
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#' alt="lifecycle-experimental"></a>
#'
#' @description
#' Computes diversity metrics (Hill numbers by default) per sample and
#' compares them across groups defined by a categorical variable using
#' estimation statistics (effect sizes + bootstrap confidence intervals)
#' via the \pkg{dabestr} package.
#'
#' This approach replaces traditional p-value-based hypothesis testing
#' with Gardner-Altman or Cumming estimation plots, following the
#' estimation statistics framework (Ho et al. 2019).
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param fact (character, required) The name of a categorical column in
#'   `sample_data` to use as the grouping factor.
#' @param q (numeric vector, default `c(0, 1, 2)`) The q values
#'   for Hill number computation: 0 = richness, 1 = Shannon exponential,
#'   2 = inverse Simpson.
#' @param custom_fn (function, default NULL) An optional custom diversity
#'   function. Must take a phyloseq object and return a named numeric
#'   vector (names = sample names) or a data.frame with one row per
#'   sample. If provided, `q` is ignored.
#' @param effect_type (character, default `"mean_diff"`) The type of
#'   effect size to compute. One of: `"mean_diff"`, `"median_diff"`,
#'   `"cohens_d"`, `"hedges_g"`, `"cliffs_delta"`.
#' @param idx (list or character vector, default NULL) The group ordering
#'   for comparisons. If NULL, uses factor levels with first level as
#'   control. For 2 groups: `c("Control", "Treatment")`. For 3+ groups:
#'   `list(c("Ctrl", "T1", "T2"))`.
#' @param ci (numeric, default 95) Confidence interval level (0-100).
#' @param resamples (integer, default 5000) Number of bootstrap
#'   resamples.
#' @param na_remove (logical, default TRUE) If TRUE, samples with NA in
#'   `fact` are removed before analysis.
#' @param ... Additional arguments passed to dabestr plotting functions.
#'
#' @return A list of class `"estim_diff_pq_result"` with components:
#'   \describe{
#'     \item{data}{The diversity data.frame used for analysis}
#'     \item{dabest_objects}{A named list of dabestr objects (one per
#'       metric)}
#'     \item{plots}{A named list of dabestr plots (one per metric)}
#'     \item{summary}{A tibble summarizing all effect sizes and CIs with
#'       columns: `metric`, `comparison`, `effect_size`, `ci_lower`,
#'       `ci_upper`, `pvalue_permtest`, `pvalue_welch`,
#'       `pvalue_mann_whitney`}
#'     \item{effect_type}{The effect size type used}
#'   }
#'
#' @details
#' The function uses the \pkg{dabestr} package to produce estimation
#' plots. For two groups, Gardner-Altman plots are produced
#' (`float_contrast = TRUE`). For three or more groups, Cumming plots are
#' used (`float_contrast = FALSE`).
#'
#' @references
#' Ho, J., Tumkaya, T., Aryal, S., Choi, H., & Claridge-Chang, A.
#' (2019). Moving beyond P values: data analysis with estimation
#' graphics. *Nature Methods*, 16(7), 565-566.
#'
#' @examples
#' library(phyloseq)
#' data("data_fungi", package = "MiscMetabar")
#'
#' pq <- subset_samples(data_fungi, !is.na(Height))
#' pq <- clean_pq(pq)
#'
#' res <- estim_diff_pq(pq, fact = "Height")
#' res
#' res$plots$Hill_0
#' res$summary
#'
#' @seealso [estim_cor_pq()], [estim_diff_lpq()], [adonis_lpq()]
#' @export
estim_diff_pq <- function(
  physeq,
  fact,
  q = c(0, 1, 2),
  custom_fn = NULL,
  effect_type = "cohens_d",
  idx = NULL,
  ci = 95,
  resamples = 5000,
  na_remove = TRUE,
  ...
) {
  # Validate inputs
  stopifnot(inherits(physeq, "phyloseq"))
  rlang::check_installed("dabestr", reason = "to compute estimation statistics")

  valid_effects <- c(
    "mean_diff",
    "median_diff",
    "cohens_d",
    "hedges_g",
    "cliffs_delta"
  )
  if (!effect_type %in% valid_effects) {
    stop(
      "effect_type must be one of: ",
      paste(valid_effects, collapse = ", ")
    )
  }

  sam <- phyloseq::sample_data(physeq)
  if (!fact %in% colnames(sam)) {
    stop("Factor '", fact, "' not found in sample_data")
  }

  # Compute diversity
  div_df <- diversity_samples_pq(physeq, q, custom_fn)

  # Handle NAs
  if (na_remove) {
    keep <- !is.na(div_df[[fact]])
    div_df <- div_df[keep, , drop = FALSE]
  }

  if (nrow(div_df) == 0) {
    stop("No samples remaining after NA removal")
  }

  # Ensure factor
  div_df[[fact]] <- as.factor(div_df[[fact]])
  levs <- levels(div_df[[fact]])

  if (length(levs) < 2) {
    stop("Factor '", fact, "' must have at least 2 levels, got ", length(levs))
  }

  # Validate minimum group size (dabestr requires >= 3 samples per group)
  group_sizes <- table(div_df[[fact]])
  if (any(group_sizes < 3)) {
    small_groups <- names(group_sizes)[group_sizes < 3]
    stop(
      "Each group needs at least 3 samples for estimation statistics. ",
      "Group(s) with insufficient samples: ",
      paste0(
        "'",
        small_groups,
        "' (n=",
        group_sizes[small_groups],
        ")",
        collapse = ", "
      )
    )
  }

  # Build idx if NULL
  if (is.null(idx)) {
    if (length(levs) == 2) {
      idx <- levs
    } else {
      idx <- list(levs)
    }
  }

  # Identify metric columns
  metric_cols <- if (is.null(custom_fn)) {
    paste0("Hill_", q)
  } else {
    setdiff(colnames(div_df), colnames(sam))
  }

  # Add sample_id column for dabestr

  div_df$.sample_id <- rownames(div_df)

  dabest_objects <- list()
  plots <- list()
  summary_rows <- list()

  fact_sym <- rlang::sym(fact)

  for (metric in metric_cols) {
    metric_sym <- rlang::sym(metric)

    dab <- dabestr::load(
      data = div_df,
      x = !!fact_sym,
      y = !!metric_sym,
      idx = idx,
      resamples = resamples,
      ci = ci
    )

    # Apply effect size function
    dab_effect <- switch(
      effect_type,
      mean_diff = dabestr::mean_diff(dab),
      median_diff = dabestr::median_diff(dab),
      cohens_d = dabestr::cohens_d(dab),
      hedges_g = dabestr::hedges_g(dab),
      cliffs_delta = dabestr::cliffs_delta(dab)
    )

    dabest_objects[[metric]] <- dab_effect

    float <- length(levs) == 2
    p <- dabestr::dabest_plot(dab_effect, float_contrast = float, ...)
    plots[[metric]] <- p

    # Extract summary from the dabestr result
    res_tbl <- dab_effect$boot_result
    perm_tbl <- dab_effect$permtest_pvals
    for (i in seq_len(nrow(res_tbl))) {
      # Extract p-values from permtest_pvals
      pval_perm <- perm_tbl$pval_permtest[i]
      pvals_list <- perm_tbl$pvalues[[i]]
      pval_welch <- pvals_list$pvalue_welch
      pval_mann_whitney <- pvals_list$pvalue_mann_whitney

      summary_rows <- c(
        summary_rows,
        list(tibble::tibble(
          metric = metric,
          comparison = paste0(
            res_tbl$control_group[i],
            " vs ",
            res_tbl$test_group[i]
          ),
          effect_size = res_tbl$difference[i],
          ci_lower = res_tbl$bca_ci_low[i],
          ci_upper = res_tbl$bca_ci_high[i],
          pvalue_permtest = pval_perm,
          pvalue_welch = pval_welch %||% NA_real_,
          pvalue_mann_whitney = pval_mann_whitney %||% NA_real_
        ))
      )
    }
  }

  summary_tbl <- dplyr::bind_rows(summary_rows)

  result <- list(
    data = div_df,
    dabest_objects = dabest_objects,
    plots = plots,
    summary = summary_tbl,
    effect_type = effect_type
  )
  class(result) <- "estim_diff_pq_result"
  result
}


#' Estimation statistics for numeric variable correlation on a phyloseq
#' object
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#' alt="lifecycle-experimental"></a>
#'
#' @description
#' Computes diversity metrics (Hill numbers by default) per sample and
#' assesses their relationship with a numeric variable using bootstrap
#' confidence intervals for correlation coefficients and regression
#' slopes.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param variable (character, required) The name of a numeric column in
#'   `sample_data`.
#' @param q (numeric vector, default `c(0, 1, 2)`) The q
#'   values for Hill number computation.
#' @param custom_fn (function, default NULL) An optional custom diversity
#'   function (see [estim_diff_pq()] for details).
#' @param method (character, default `"pearson"`) Correlation method.
#'   One of `"pearson"`, `"spearman"`, `"kendall"`.
#' @param resamples (integer, default 5000) Number of bootstrap
#'   resamples.
#' @param ci (numeric, default 95) Confidence interval level (0-100).
#' @param na_remove (logical, default TRUE) If TRUE, samples with NA in
#'   `variable` are removed.
#'
#' @return A list of class `"estim_cor_pq_result"` with components:
#'   \describe{
#'     \item{data}{The diversity data.frame used for analysis}
#'     \item{correlations}{A tibble with columns: `metric`, `estimate`,
#'       `ci_lower`, `ci_upper`, `method`, `pvalue`}
#'     \item{regressions}{A tibble with columns: `metric`, `intercept`,
#'       `slope`, `slope_ci_lower`, `slope_ci_upper`}
#'     \item{plots}{A named list of ggplot2 scatter plots with
#'       regression line and bootstrap CI ribbon (one per metric)}
#'   }
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data("data_fungi", package = "MiscMetabar")
#'
#' # Add a numeric variable for demonstration
#' sam <- sample_data(data_fungi)
#' sam$lib_size <- sample_sums(data_fungi)
#' sample_data(data_fungi) <- sam
#'
#' res <- estim_cor_pq(data_fungi, variable = "lib_size")
#' res
#' res$plots$Hill_0
#' res$correlations
#' }
#'
#' @seealso [estim_diff_pq()], [estim_cor_lpq()]
#' @export
estim_cor_pq <- function(
  physeq,
  variable,
  q = c(0, 1, 2),
  custom_fn = NULL,
  method = "pearson",
  resamples = 5000,
  ci = 95,
  na_remove = TRUE
) {
  stopifnot(inherits(physeq, "phyloseq"))

  sam <- phyloseq::sample_data(physeq)
  if (!variable %in% colnames(sam)) {
    stop("Variable '", variable, "' not found in sample_data")
  }
  if (!is.numeric(sam[[variable]])) {
    stop("Variable '", variable, "' must be numeric")
  }

  valid_methods <- c("pearson", "spearman", "kendall")
  if (!method %in% valid_methods) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # Compute diversity
  div_df <- diversity_samples_pq(physeq, q, custom_fn)

  # Handle NAs
  if (na_remove) {
    keep <- !is.na(div_df[[variable]])
    div_df <- div_df[keep, , drop = FALSE]
  }

  if (nrow(div_df) < 3) {
    stop("Need at least 3 non-NA samples for correlation")
  }

  # Identify metric columns
  metric_cols <- if (is.null(custom_fn)) {
    paste0("Hill_", q)
  } else {
    setdiff(
      colnames(div_df),
      colnames(as.data.frame(phyloseq::sample_data(physeq)))
    )
  }

  x_vals <- div_df[[variable]]

  cor_rows <- list()
  reg_rows <- list()
  plots <- list()

  for (metric in metric_cols) {
    y_vals <- div_df[[metric]]

    # Bootstrap correlation
    cor_res <- bootstrap_cor(
      x_vals,
      y_vals,
      method = method,
      resamples = resamples,
      ci = ci
    )
    # Classical cor.test p-value for legacy purposes
    cor_test <- stats::cor.test(x_vals, y_vals, method = method)
    cor_rows <- c(
      cor_rows,
      list(tibble::tibble(
        metric = metric,
        estimate = cor_res$estimate,
        ci_lower = cor_res$ci_lower,
        ci_upper = cor_res$ci_upper,
        method = method,
        pvalue = cor_test$p.value
      ))
    )

    # Bootstrap regression
    reg_res <- bootstrap_lm(x_vals, y_vals, resamples = resamples, ci = ci)
    reg_rows <- c(
      reg_rows,
      list(tibble::tibble(
        metric = metric,
        intercept = reg_res$intercept,
        slope = reg_res$slope,
        slope_ci_lower = reg_res$slope_ci_lower,
        slope_ci_upper = reg_res$slope_ci_upper
      ))
    )

    # Scatter plot with regression line and CI ribbon
    plot_df <- data.frame(x = x_vals, y = y_vals)

    # Compute prediction ribbon from bootstrap slopes
    x_seq <- seq(min(x_vals), max(x_vals), length.out = 100)
    boot_preds <- vapply(
      reg_res$boot_slopes,
      \(s) {
        reg_res$intercept + s * x_seq
      },
      numeric(100)
    )
    alpha <- (100 - ci) / 200
    ribbon_lower <- apply(
      boot_preds,
      1,
      stats::quantile,
      probs = alpha,
      names = FALSE
    )
    ribbon_upper <- apply(
      boot_preds,
      1,
      stats::quantile,
      probs = 1 - alpha,
      names = FALSE
    )
    ribbon_df <- data.frame(
      x = x_seq,
      lower = ribbon_lower,
      upper = ribbon_upper,
      fitted = reg_res$intercept + reg_res$slope * x_seq
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_ribbon(
        data = ribbon_df,
        ggplot2::aes(
          x = .data$x,
          ymin = .data$lower,
          ymax = .data$upper,
          y = NULL
        ),
        alpha = 0.2,
        fill = "steelblue"
      ) +
      ggplot2::geom_line(
        data = ribbon_df,
        ggplot2::aes(x = .data$x, y = .data$fitted),
        color = "steelblue",
        linewidth = 1
      ) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::labs(
        x = variable,
        y = metric,
        title = paste0(metric, " vs ", variable),
        subtitle = paste0(
          method,
          " r = ",
          round(cor_res$estimate, 3),
          " [",
          round(cor_res$ci_lower, 3),
          ", ",
          round(cor_res$ci_upper, 3),
          "]"
        )
      ) +
      ggplot2::theme_minimal()

    plots[[metric]] <- p
  }

  result <- list(
    data = div_df,
    correlations = dplyr::bind_rows(cor_rows),
    regressions = dplyr::bind_rows(reg_rows),
    plots = plots
  )
  class(result) <- "estim_cor_pq_result"
  result
}


# ==============================================================================
# Print methods
# ==============================================================================

#' @export
print.estim_diff_pq_result <- function(x, ...) {
  cat("Estimation statistics: categorical comparison\n")
  cat("Effect size type:", x$effect_type, "\n")
  cat("Metrics:", paste(unique(x$summary$metric), collapse = ", "), "\n")
  cat("Comparisons:", length(unique(x$summary$comparison)), "\n\n")
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
print.estim_cor_pq_result <- function(x, ...) {
  cat("Estimation statistics: numeric correlation\n")
  cat("Metrics:", paste(unique(x$correlations$metric), collapse = ", "), "\n\n")
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
