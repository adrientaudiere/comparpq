#' Scatter plots of Hill diversity across pairs of phyloseq objects
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#' alt="lifecycle-experimental"></a>
#'
#' @description
#' For each pair of phyloseq objects in a [list_phyloseq], computes Hill
#' numbers (richness, Shannon exponential, inverse Simpson) per sample and
#' draws scatter plots comparing the diversity of the **same samples** across
#' objects. This is a visual diagnostic for REPRODUCIBILITY, ROBUSTNESS, and
#' REPLICABILITY comparisons.
#'
#' Each point represents one sample. A 1:1 reference line and optional
#' regression line help assess how similar (or different) the diversity
#' estimates are between two pipelines, runs, or primers.
#'
#' @param x (required) A [list_phyloseq] object.
#' @param hill_scales (numeric vector, default `c(0, 1, 2)`) Hill number
#'   orders to compute: 0 = richness, 1 = Shannon exponential, 2 = inverse
#'   Simpson.
#' @param pairs (list of integer pairs or NULL, default NULL) Which pairs of
#'   phyloseq objects to compare. Each element must be a length-2 integer
#'   vector of indices. If NULL, all pairwise combinations are used.
#' @param add_1to1 (logical, default TRUE) If TRUE, adds a dashed 1:1
#'   identity line to each panel.
#' @param add_smooth (logical, default TRUE) If TRUE, adds a linear
#'   regression line with confidence interval via [ggplot2::geom_smooth()].
#' @param cor_method (character, default `"pearson"`) Correlation method for
#'   the annotation label. One of `"pearson"`, `"spearman"`, `"kendall"`.
#' @param point_color (character, default `"black"`) Color for scatter
#'   points.
#' @param point_alpha (numeric, default 0.6) Transparency for scatter points.
#' @param point_size (numeric, default 2) Size for scatter points.
#' @param smooth_color (character, default `"steelblue"`) Color for the
#'   regression line and confidence ribbon.
#' @param verbose (logical, default TRUE) If TRUE, print a message when
#'   common samples are used instead of all samples.
#'
#' @return A [ggplot2::ggplot()] object with panels arranged in a grid:
#'   rows = pairs of phyloseq objects, columns = Hill number orders.
#'
#' @details
#' The function works on **common samples** across all phyloseq objects. When
#' the [list_phyloseq] has `same_samples = TRUE`, all samples are used. When
#' samples differ (e.g., NESTED_ROBUSTNESS), only the intersection is kept.
#'
#' Comparison type context:
#' \describe{
#'   \item{REPRODUCIBILITY}{High correlation expected for all Hill orders.}
#'   \item{ROBUSTNESS}{Moderate to high correlation is desirable; deviations
#'     indicate pipeline sensitivity.}
#'   \item{REPLICABILITY}{Correlation may vary across Hill orders depending on
#'     primer or technology differences.}
#' }
#'
#' @examples
#' lpq <- list_phyloseq(
#'   list(run1 = data_fungi, run2 = data_fungi),
#'   same_bioinfo_pipeline = FALSE
#' )
#'
#' gg_hill_lpq(lpq)
#' gg_hill_lpq(lpq, hill_scales = c(0, 1))
#' gg_hill_lpq(lpq, add_smooth = FALSE, add_1to1 = TRUE)
#'
#' @seealso [estim_cor_pq()], [estim_diff_lpq()]
#' @export
gg_hill_lpq <- function(
  x,
  hill_scales = c(0, 1, 2),
  pairs = NULL,
  add_1to1 = TRUE,
  add_smooth = TRUE,
  cor_method = "pearson",
  point_color = "black",
  point_alpha = 0.6,
  point_size = 2,
  smooth_color = "steelblue",
  verbose = TRUE
) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  valid_methods <- c("pearson", "spearman", "kendall")
  if (!cor_method %in% valid_methods) {
    stop("cor_method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  if (length(x@phyloseq_list) < 2) {
    stop(
      "list_phyloseq must contain at least 2 phyloseq objects to compare pairs"
    )
  }

  # ---- Determine common samples ------------------------------------------------
  common_samps <- x@comparison$common_samples

  if (length(common_samps) < 2) {
    stop(
      "At least 2 common samples are required across all phyloseq objects. ",
      "Found: ",
      length(common_samps)
    )
  }

  if (verbose && !x@comparison$same_samples) {
    message(
      "Samples differ across phyloseq objects. Using ",
      length(common_samps),
      " common samples."
    )
  }

  # ---- Compute Hill numbers per phyloseq (on common samples) ------------------
  hill_list <- purrr::imap(x@phyloseq_list, function(pq, name) {
    pq_sub <- phyloseq::prune_samples(
      phyloseq::sample_names(pq) %in% common_samps,
      pq
    )
    hill_df <- hill_samples_pq(pq_sub, hill_scales)
    hill_df <- hill_df[, paste0("Hill_", hill_scales), drop = FALSE]
    hill_df$sample_name <- rownames(hill_df)
    hill_df
  })

  # ---- Determine pairs to plot -------------------------------------------------
  nms <- names(x@phyloseq_list)
  n_pq <- length(nms)

  if (is.null(pairs)) {
    pair_indices <- utils::combn(seq_len(n_pq), 2, simplify = FALSE)
  } else {
    pair_indices <- pairs
    bad <- purrr::map_lgl(pair_indices, ~ length(.x) != 2)
    if (any(bad)) {
      stop("Each element of 'pairs' must be an integer vector of length 2")
    }
  }

  if (length(pair_indices) == 0) {
    stop("No pairs to compare")
  }

  # ---- Build long data frame for plotting -------------------------------------
  hill_cols <- paste0("Hill_", hill_scales)

  plot_data <- purrr::map_dfr(pair_indices, function(idx) {
    i <- idx[1]
    j <- idx[2]
    df_i <- hill_list[[i]]
    df_j <- hill_list[[j]]

    merged <- merge(
      df_i,
      df_j,
      by = "sample_name",
      suffixes = c("_x", "_y")
    )

    purrr::map_dfr(hill_cols, function(col) {
      data.frame(
        sample_name = merged$sample_name,
        pair = paste0(nms[i], "\nvs\n", nms[j]),
        pq_x = nms[i],
        pq_y = nms[j],
        hill_scale = col,
        x_value = merged[[paste0(col, "_x")]],
        y_value = merged[[paste0(col, "_y")]],
        stringsAsFactors = FALSE
      )
    })
  })

  # Order pair factor to preserve the original object ordering
  pair_levels <- purrr::map_chr(pair_indices, function(idx) {
    paste0(nms[idx[1]], "\nvs\n", nms[idx[2]])
  })
  plot_data$pair <- factor(plot_data$pair, levels = unique(pair_levels))
  plot_data$hill_scale <- factor(plot_data$hill_scale, levels = hill_cols)

  # ---- Compute per-panel correlation labels -----------------------------------
  cor_labels <- plot_data |>
    dplyr::group_by(.data$pair, .data$hill_scale) |>
    dplyr::summarise(
      n = dplyr::n(),
      r = stats::cor(.data$x_value, .data$y_value, method = cor_method),
      label = paste0(
        cor_method,
        " r = ",
        round(.data$r, 3),
        "\nn = ",
        .data$n
      ),
      .groups = "drop"
    )

  # ---- Build the ggplot -------------------------------------------------------
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$x_value, y = .data$y_value)
  )

  if (add_1to1) {
    p <- p +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "grey50",
        linewidth = 0.6
      )
  }

  if (add_smooth) {
    p <- p +
      ggplot2::geom_smooth(
        method = "lm",
        formula = y ~ x,
        color = smooth_color,
        fill = smooth_color,
        alpha = 0.2,
        linewidth = 0.8
      )
  }

  p <- p +
    ggplot2::geom_point(
      color = point_color,
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::geom_text(
      data = cor_labels,
      ggplot2::aes(label = .data$label),
      x = Inf,
      y = -Inf,
      hjust = 1.1,
      vjust = -0.3,
      size = 3,
      color = "grey30",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$pair),
      cols = ggplot2::vars(.data$hill_scale),
      scales = "free"
    ) +
    ggplot2::labs(
      x = "Hill diversity (phyloseq x)",
      y = "Hill diversity (phyloseq y)",
      title = paste0(
        "Hill diversity across pairs\n(",
        x@comparison$type_of_comparison,
        ")"
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 9),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey80")
    )

  p
}
