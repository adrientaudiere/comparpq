################################################################################
#' Circle-packed bubble plot of a phyloseq object using ggplot2
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Creates a static circle-packed bubble plot of taxa abundances from a
#' phyloseq object using ggplot2. Circles can be packed in a circular layout
#' (tight, default) or a square layout. Optionally facets the plot by a sample
#' data variable, producing one bubble chart per level.
#'
#' @inheritParams tc_points_matrix
#' @param rank_label (character, default "Taxa") The name of the column in the
#'   `@tax_table` slot to label the circles. If set to "Taxa", the taxa names
#'   are used.
#' @param rank_color (character, default "Family") The name of the column in
#'   the `@tax_table` slot to color the circles.
#' @param layout (character, default "circle") The packing layout. `"circle"`
#'   produces a tight circular packing. `"square"` constrains circles inside a
#'   square boundary, with large circles placed centrally.
#' @param facet_by (character, default NULL) A column name from `@sam_data` to
#'   facet the plot. When set, one bubble chart is produced per level of the
#'   variable, with taxa abundances computed within each level.
#' @param log1ptransform (logical, default FALSE) If TRUE, the number of
#'   sequences is log1p transformed before computing circle sizes.
#' @param min_nb_seq (integer, default 0) Minimum number of sequences to filter
#'   out taxa with low abundance.
#' @param label_size (numeric, default 2) Font size for the labels inside
#'   circles.
#' @param label_color (character, default "grey10") Color for the label text.
#' @param show_labels (logical, default TRUE) If TRUE, labels are displayed
#'   inside circles. Only circles large enough to fit text are labeled.
#' @param border_color (character, default "white") Color for circle borders.
#' @param border_width (numeric, default 0.5) Width of circle borders.
#' @param alpha (numeric, default 0.8) Transparency of circle fill.
#' @param npoints (integer, default 50) Number of vertices used to approximate
#'   each circle polygon. Higher values produce smoother circles.
#' @param ncol_facet (integer, default NULL) Number of columns for facet layout.
#'   Passed to [ggplot2::facet_wrap()].
#' @param return_dataframe (logical, default FALSE) If TRUE, the plot is not
#'   returned, but the resulting dataframe to plot is returned.
#'
#' @return A ggplot2 object, or a data.frame if `return_dataframe = TRUE`.
#' @export
#' @author Adrien Taudiere
#'
#' @examples
#' gg_bubbles_pq(physeq = data_fungi_mini, rank_color = "Class")
#' gg_bubbles_pq(
#'   physeq = data_fungi_mini, rank_color = "Class",
#'   layout = "square"
#' )
#' \dontrun{
#' gg_bubbles_pq(
#'   physeq = data_fungi, rank_color = "Order",
#'   facet_by = "Height", min_nb_seq = 100
#' ) + no_legend()
#' }
gg_bubbles_pq <- function(
  physeq,
  rank_label = "Taxa",
  rank_color = "Family",
  layout = "circle",
  facet_by = NULL,
  log1ptransform = FALSE,
  min_nb_seq = 0,
  label_size = 2,
  label_color = "grey10",
  show_labels = TRUE,
  border_color = "white",
  border_width = 0.5,
  alpha = 0.8,
  npoints = 50,
  ncol_facet = NULL,
  return_dataframe = FALSE
) {
  verify_pq(physeq)
  layout <- match.arg(layout, c("circle", "square"))

  if (!is.null(facet_by)) {
    if (!facet_by %in% colnames(physeq@sam_data)) {
      stop(
        paste0(
          "'",
          facet_by,
          "' not found in sample_data of the phyloseq object."
        )
      )
    }
  }

  build_bubble_df <- function(pq, facet_label = NULL) {
    if (rank_label == "Taxa") {
      label <- taxa_names(pq)
    } else {
      label <- as.vector(pq@tax_table[, rank_label])
    }

    df <- data.frame(
      value = taxa_sums(pq),
      label = label,
      rank_value_color = as.vector(pq@tax_table[, rank_color]),
      stringsAsFactors = FALSE
    )

    if (min_nb_seq > 0) {
      df <- df[df$value > min_nb_seq, ]
    }

    if (nrow(df) == 0) {
      return(NULL)
    }

    if (log1ptransform) {
      df$value <- log1p(df$value)
    }

    if (!is.null(facet_label)) {
      df$facet <- facet_label
    }

    df
  }

  if (is.null(facet_by)) {
    df_all <- build_bubble_df(physeq)
  } else {
    levels_var <- unique(as.character(physeq@sam_data[[facet_by]]))
    levels_var <- levels_var[!is.na(levels_var)]
    df_list <- lapply(levels_var, function(lvl) {
      sam_values <- as.character(physeq@sam_data[[facet_by]])
      idx <- !is.na(sam_values) & sam_values == lvl
      if (sum(idx) == 0) {
        return(NULL)
      }
      pq_sub <- prune_samples(idx, physeq)
      pq_sub <- prune_taxa(taxa_sums(pq_sub) > 0, pq_sub)
      if (ntaxa(pq_sub) == 0) {
        return(NULL)
      }
      build_bubble_df(pq_sub, facet_label = lvl)
    })
    df_all <- do.call(rbind, df_list)
  }

  if (is.null(df_all) || nrow(df_all) == 0) {
    stop("No taxa remaining after filtering. Try lowering `min_nb_seq`.")
  }

  if (return_dataframe) {
    return(df_all)
  }

  compute_layout <- function(df) {
    lay <- packcircles::circleProgressiveLayout(
      df$value,
      sizetype = "area"
    )

    if (layout == "square") {
      # Map circle centers from disc to square using radial stretching.
      # Each center is moved outward along its angle so that the overall
      # circular silhouette fills a square. Radii stay unchanged.
      cx <- mean(range(lay$x))
      cy <- mean(range(lay$y))
      dx <- lay$x - cx
      dy <- lay$y - cy
      r <- sqrt(dx^2 + dy^2)
      theta <- atan2(dy, dx)

      # For each angle, the square boundary is farther out than the circle
      # boundary by factor 1 / max(|cos θ|, |sin θ|).
      stretch <- ifelse(
        r > 0,
        1 / pmax(abs(cos(theta)), abs(sin(theta))),
        1
      )
      lay$x <- cx + dx * stretch
      lay$y <- cy + dy * stretch

      # Resolve any overlaps introduced by the stretching
      init_dat <- data.frame(x = lay$x, y = lay$y, radius = lay$radius)
      xrng <- range(lay$x - lay$radius, lay$x + lay$radius)
      yrng <- range(lay$y - lay$radius, lay$y + lay$radius)
      pad <- max(lay$radius) * 0.1
      lay <- packcircles::circleRepelLayout(
        init_dat,
        xlim = c(xrng[1] - pad, xrng[2] + pad),
        ylim = c(yrng[1] - pad, yrng[2] + pad),
        xysizecols = c(1, 2, 3),
        sizetype = "radius",
        maxiter = 1000,
        wrap = FALSE
      )$layout
    }

    vertices <- packcircles::circleLayoutVertices(lay, npoints = npoints)
    vertices$label <- rep(df$label, each = npoints + 1)
    vertices$rank_value_color <- rep(df$rank_value_color, each = npoints + 1)

    center <- data.frame(
      x = lay$x,
      y = lay$y,
      radius = lay$radius,
      label = df$label
    )

    list(vertices = vertices, center = center)
  }

  if (is.null(facet_by)) {
    result <- compute_layout(df_all)
    verts <- result$vertices
    centers <- result$center
  } else {
    facet_levels <- unique(df_all$facet)
    verts_list <- list()
    centers_list <- list()
    for (fl in facet_levels) {
      df_sub <- df_all[df_all$facet == fl, ]
      result <- compute_layout(df_sub)
      result$vertices$facet <- fl
      result$center$facet <- fl
      verts_list[[fl]] <- result$vertices
      centers_list[[fl]] <- result$center
    }
    verts <- do.call(rbind, verts_list)
    centers <- do.call(rbind, centers_list)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = verts,
      ggplot2::aes(
        x = x,
        y = y,
        group = id,
        fill = rank_value_color
      ),
      color = border_color,
      linewidth = border_width,
      alpha = alpha
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  if (show_labels) {
    max_radius <- max(centers$radius, na.rm = TRUE)
    label_threshold <- max_radius * 0.05
    centers_labeled <- centers[centers$radius > label_threshold, ]

    if (nrow(centers_labeled) > 0) {
      p <- p +
        ggplot2::geom_text(
          data = centers_labeled,
          ggplot2::aes(x = x, y = y, label = label),
          size = label_size,
          color = label_color
        )
    }
  }

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(~facet, ncol = ncol_facet)
  }

  p
}
################################################################################
