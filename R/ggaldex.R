#' ggplot2 version of ALDEx2 diagnostic plots
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Creates ggplot2 versions of the four ALDEx2 diagnostic plot types:
#' MW (Bland-Altman style), MA, volcano, and volcano.var. Supports
#' faceting by a `name` column, making it suitable for output from
#' [aldex_lpq()].
#'
#' @param x (data.frame or tibble, required) ALDEx2 results, typically
#'   from [ALDEx2::aldex()] or [aldex_lpq()]. Must contain columns
#'   `diff.btw`, `diff.win`, `rab.all`, and `effect`. When `test` is
#'   `"welch"`, column `we.eBH` is required; when `"wilcox"`, column
#'   `wi.eBH` is required.
#' @param type (character, default "MW") Plot type. One of `"MW"`
#'   (dispersion vs difference), `"MA"` (abundance vs difference),
#'   `"volcano"` (difference vs -log10 q-value), or `"volcano.var"`
#'   (dispersion vs -log10 q-value).
#' @param test (character, default "welch") Statistical test used for
#'   significance calling. One of `"welch"`, `"wilcox"`, `"effect"`,
#'   or `"both"` (effect + welch).
#' @param cutoff.pval (numeric, default 0.05) q-value threshold for
#'   significance.
#' @param cutoff.effect (numeric, default 1) Effect size threshold
#'   (must be >= 0.5 when `test` is `"effect"` or `"both"`).
#' @param rare (numeric, default 0) Abundance threshold below which
#'   taxa are marked as rare (only used for MW and MA plots).
#' @param all.col (character, default "grey30") Color for non-significant
#'   points.
#' @param called.col (character, default "red") Color for significant
#'   points.
#' @param rare.col (character, default "black") Color for rare taxa
#'   points.
#' @param point.size (numeric, default 1.5) Size of points.
#' @param point.alpha (numeric, default 0.4) Alpha transparency of
#'   non-significant points.
#'
#' @return A ggplot2 object. If `x` contains a `name` column (e.g.,
#'   from [aldex_lpq()]), the plot is faceted by `name`.
#'
#' @details
#' This function reimplements [ALDEx2::aldex.plot()] using ggplot2,
#' providing a more customizable and composable output. The four plot
#' types correspond to the original ALDEx2 types:
#'
#' \describe{
#'   \item{MW}{Bland-Altman style: within-condition dispersion (x) vs
#'     between-condition difference (y), with +/-1 effect size lines.}
#'   \item{MA}{MA plot: median log2 relative abundance (x) vs median
#'     log2 difference (y).}
#'   \item{volcano}{Volcano plot: median log2 difference (x) vs
#'     -log10 q-value (y).}
#'   \item{volcano.var}{Variance volcano: median log2 dispersion (x)
#'     vs -log10 q-value (y).}
#' }
#'
#' @examples
#'
#' data_fungi_high <- multiply_counts_pq(data_fungi, "Height", "High",
#'   4,
#'   prop_taxa = 0.1, seed = 42
#' )
#'
#' aldex_pq(data_fungi_high,
#'   bifactor = "Height",
#'   modalities = c("Low", "High")
#' ) |>
#'   gg_aldex_plot()
#'
#' lpq <- list_phyloseq(
#'   list(
#'     fungi = data_fungi,
#'     fungi_height = data_fungi_high
#'   ),
#'   same_bioinfo_pipeline = FALSE
#' )
#' # From aldex_lpq (faceted by name)
#' lpq_res <- aldex_lpq(lpq,
#'   bifactor = "Height",
#'   modalities = c("Low", "High")
#' )
#' gg_aldex_plot(lpq_res, type = "MA")
#' gg_aldex_plot(lpq_res, type = "MW", test = "wilcox")
#'
#' gg_aldex_plot(lpq_res, type = "volcano")
#' gg_aldex_plot(lpq_res, type = "volcano.var")
#'
#'
#' gingival_pq <-
#'   MicrobiomeBenchmarkData::getBenchmarkData("HMP_2012_16S_gingival_V35_subset",
#'     dryrun = FALSE
#'   )[[1]] |>
#'   mia::convertToPhyloseq()
#'
#' aldex_res <- aldex_pq(gingival_pq,
#'   bifactor = "body_subsite",
#'   modalities = c("supragingival_plaque", "subgingival_plaque")
#' )
#'
#' gg_aldex_plot(aldex_res, type = "volcano")
#'
#' @seealso [aldex_lpq()], [ALDEx2::aldex.plot()]
#' @export
gg_aldex_plot <- function(
  x,
  type = c("MW", "MA", "volcano", "volcano.var"),
  test = c("welch", "wilcox", "effect", "both"),
  cutoff.pval = 0.05,
  cutoff.effect = 1,
  rare = 0,
  all.col = "grey30",
  called.col = "red",
  rare.col = "black",
  point.size = 1.5,
  point.alpha = 0.4
) {
  type <- match.arg(type)
  test <- match.arg(test)

  if (length(x$effect) == 0) {
    stop("Column 'effect' not found. Please run aldex.effect before plotting.")
  }

  # Determine significance
  if (test == "welch") {
    if (is.null(x$we.eBH)) {
      stop("Welch test results (we.eBH) not in dataset.")
    }
    p.add <- min(x$we.eBH[x$we.eBH > 0], na.rm = TRUE) / 10
    x$.called <- x$we.eBH <= cutoff.pval
    x$.qval <- x$we.eBH + p.add
  } else if (test == "wilcox") {
    if (is.null(x$wi.eBH)) {
      stop("Wilcoxon test results (wi.eBH) not in dataset.")
    }
    p.add <- min(x$wi.eBH[x$wi.eBH > 0], na.rm = TRUE) / 10
    x$.called <- x$wi.eBH <= cutoff.pval
    x$.qval <- x$wi.eBH + p.add
  } else if (test == "effect") {
    if (cutoff.effect < 0.5) {
      stop("Please set cutoff.effect to at least 0.5.")
    }
    x$.called <- abs(x$effect) >= cutoff.effect
  } else if (test == "both") {
    if (cutoff.effect < 0.5) {
      stop("Please set cutoff.effect to at least 0.5.")
    }
    if (is.null(x$we.eBH)) {
      stop("Welch test results (we.eBH) not in dataset.")
    }
    p.add <- min(x$we.eBH[x$we.eBH > 0], na.rm = TRUE) / 10
    x$.called <- abs(x$effect) >= cutoff.effect & x$we.eBH <= cutoff.pval
    x$.qval <- x$we.eBH + p.add
  }

  # Determine rare status
  x$.rare <- x$rab.all < rare

  # Build status column for coloring
  x$.status <- "not significant"
  x$.status[x$.rare] <- "rare"
  x$.status[x$.called] <- "significant"
  x$.status <- factor(
    x$.status,
    levels = c("not significant", "rare", "significant")
  )

  color_values <- c(
    "not significant" = all.col,
    "rare" = rare.col,
    "significant" = called.col
  )

  has_facet <- "name" %in% colnames(x)

  if (type == "MW") {
    p <- ggplot2::ggplot(
      x,
      ggplot2::aes(x = diff.win, y = diff.btw, color = .status)
    ) +
      ggplot2::geom_point(size = point.size, alpha = point.alpha) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "darkgrey"
      ) +
      ggplot2::geom_abline(
        slope = -1,
        intercept = 0,
        linetype = "dashed",
        color = "darkgrey"
      ) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::labs(
        x = expression("Median" ~ Log[2] ~ "Dispersion"),
        y = expression("Median" ~ Log[2] ~ "Difference"),
        color = NULL
      ) +
      ggplot2::theme_minimal()
  } else if (type == "MA") {
    p <- ggplot2::ggplot(
      x,
      ggplot2::aes(x = rab.all, y = diff.btw, color = .status)
    ) +
      ggplot2::geom_point(size = point.size, alpha = point.alpha) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::labs(
        x = expression("Median" ~ Log[2] ~ "Relative Abundance"),
        y = expression("Median" ~ Log[2] ~ "Difference"),
        color = NULL
      ) +
      ggplot2::theme_minimal()
  } else if (type == "volcano") {
    x$.neg_log10_q <- -1 * log10(x$.qval)
    p <- ggplot2::ggplot(
      x,
      ggplot2::aes(x = diff.btw, y = .neg_log10_q, color = .status)
    ) +
      ggplot2::geom_point(size = point.size, alpha = point.alpha) +
      ggplot2::geom_hline(
        yintercept = -1 * log10(cutoff.pval),
        linetype = "dashed",
        color = "darkgrey"
      ) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::labs(
        x = expression("Median" ~ Log[2] ~ "Difference"),
        y = expression(-Log[10] ~ "q-value"),
        color = NULL
      ) +
      ggplot2::theme_minimal()
  } else if (type == "volcano.var") {
    x$.neg_log10_q <- -1 * log10(x$.qval)
    p <- ggplot2::ggplot(
      x,
      ggplot2::aes(x = diff.win, y = .neg_log10_q, color = .status)
    ) +
      ggplot2::geom_point(size = point.size, alpha = point.alpha) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::labs(
        x = expression("Median" ~ Log[2] ~ "Dispersion"),
        y = expression(-Log[10] ~ "q-value"),
        color = NULL
      ) +
      ggplot2::theme_minimal()
  }

  if (has_facet) {
    p <- p + ggplot2::facet_wrap(~name)
  }

  p
}
