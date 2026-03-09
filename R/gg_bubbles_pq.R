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
#' When a [list_phyloseq] is passed as `physeq`, it is first merged into a
#' single phyloseq object using [merge_lpq()] (each original phyloseq becomes
#' one sample) and the plot is automatically faceted by `source_name`.
#'
#' When `diff_contour = TRUE` together with `facet_by` (or a list_phyloseq
#' input), all pairwise comparisons between facet levels are shown side by
#' side using \pkg{patchwork}. For each pair (A vs B), taxa unique to A are
#' highlighted with A's color and taxa unique to B with B's color. Shared
#' taxa receive a transparent contour. This makes it easy to spot which taxa
#' are exclusive to each group in every pairwise comparison.
#'
#' @inheritParams tc_points_matrix
#' @param rank_label (character, default "Taxa") The name of the column in the
#'   `@tax_table` slot to label the circles. If set to "Taxa", the taxa names
#'   are used.
#' @param rank_color (character, default "Family") The name of the column in
#'   the `@tax_table` slot to color the circles.
#' @param rank_contour (character, default NULL) The name of a column in the
#'   `@tax_table` slot to color the circle borders (contours). When NULL, the
#'   fixed `border_color` is used for all borders. Ignored when
#'   `diff_contour = TRUE`.
#' @param layout (character, default "circle") The packing layout. `"circle"`
#'   produces a tight circular packing. `"square"` constrains circles inside a
#'   square boundary, with large circles placed centrally.
#' @param facet_by (character, default NULL) A column name from `@sam_data` to
#'   facet the plot. When set, one bubble chart is produced per level of the
#'   variable, with taxa abundances computed within each level. When `physeq`
#'   is a [list_phyloseq], this is automatically set to `"source_name"`.
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
#'   Passed to [ggplot2::facet_wrap()]. Not used if diff_contour is TRUE since
#'   patchwork is used for layout instead.
#' @param return_dataframe (logical, default FALSE) If TRUE, the plot is not
#'   returned, but the resulting dataframe to plot is returned. Ignored when
#'   `diff_contour = TRUE`.
#' @param diff_contour (logical, default FALSE) If TRUE and `facet_by` is set
#'   (or `physeq` is a [list_phyloseq]), produces pairwise comparison panels
#'   for all pairs of facet levels using \pkg{patchwork}. For each pair, taxa
#'   unique to each side are highlighted with a distinct contour color from
#'   `diff_contour_colors`. Shared taxa get a transparent contour. When TRUE,
#'   `rank_contour` is ignored.
#' @param diff_contour_colors (character vector, default
#'   `c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")`) Border colors for taxa
#'   unique to each facet level in `diff_contour` mode. Recycled if shorter
#'   than the number of facet levels. Each level gets a distinct color so
#'   unique taxa from different groups are visually distinguishable.
#' @param diff_border_width (numeric, default 1.5) Border width used in
#'   `diff_contour` mode.
#' @param show_title (logical, default TRUE) If TRUE, adds an informative
#'   title describing what the fill color, contour color, circle size, and
#'   labels represent.
#' @param match_by (character, default `"refseq"`) How to match taxa when
#'   `physeq` is a [list_phyloseq]. Passed to [merge_lpq()]. One of
#'   `"refseq"` (match by reference sequences) or `"names"` (match by taxa
#'   names).
#'
#' @return A ggplot2 object, a patchwork object (when `diff_contour = TRUE`),
#'   or a data.frame if `return_dataframe = TRUE`.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' gg_bubbles_pq(physeq = data_fungi_mini, rank_color = "Class")
#' gg_bubbles_pq(
#'   physeq = data_fungi_mini, rank_color = "Class",
#'   rank_contour = "Order"
#' )
#' gg_bubbles_pq(
#'   physeq = data_fungi_mini, rank_color = "Class",
#'   layout = "square"
#' )
#'
#' # Faceted by sample variable
#' gg_bubbles_pq(
#'   physeq = data_fungi, rank_color = "Order",
#'   facet_by = "Height", min_nb_seq = 100
#' ) + no_legend()
#'
#' # Pairwise diff_contour on a faceted phyloseq
#' gg_bubbles_pq(
#'   physeq = data_fungi, rank_color = "Order",
#'   facet_by = "Height", min_nb_seq = 100,
#'   diff_contour = TRUE, show_labels = FALSE
#' ) & no_legend()
#'
#' # list_phyloseq: automatically merged and faceted
#' mini2 <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) < 10000)
#' lpq <- list_phyloseq(
#'   list(full = data_fungi, mini = data_fungi_mini, mini2 = mini2),
#' )
#' gg_bubbles_pq(lpq, rank_color = "Class")
#'
#' # list_phyloseq with diff_contour: pairwise panels
#' gg_bubbles_pq(lpq, rank_color = "Class", diff_contour = TRUE,
#'   show_labels = FALSE, diff_border_width = 1) & no_legend()
#'
gg_bubbles_pq <- function(
  physeq,
  rank_label = "Taxa",
  rank_color = "Family",
  rank_contour = NULL,
  layout = c("circle", "square"),
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
  return_dataframe = FALSE,
  diff_contour = FALSE,
  diff_contour_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
  diff_border_width = 1.5,
  show_title = TRUE,
  match_by = c("refseq", "names")
) {
  # ---- Branch: list_phyloseq -> merge and facet ------------------------------
  if (inherits(physeq, "comparpq::list_phyloseq")) {
    match_by <- match.arg(match_by)
    physeq <- merge_lpq(physeq, match_by = match_by, verbose = FALSE)
    facet_by <- "source_name"
  }

  # ---- Single phyloseq input ------------------------------------------------
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

  if (diff_contour && is.null(facet_by)) {
    warning(
      "'diff_contour' requires 'facet_by' to be set. ",
      "Ignoring diff_contour."
    )
    diff_contour <- FALSE
  }

  if (diff_contour && !is.null(rank_contour)) {
    message("'rank_contour' is ignored when 'diff_contour = TRUE'.")
    rank_contour <- NULL
  }

  # ---- diff_contour: pairwise comparison panels ------------------------------
  if (diff_contour) {
    rlang::check_installed(
      "patchwork",
      reason = "to combine pairwise diff_contour panels"
    )

    levels_var <- unique(as.character(physeq@sam_data[[facet_by]]))
    levels_var <- levels_var[!is.na(levels_var)]

    diff_colors <- stats::setNames(
      rep_len(diff_contour_colors, length(levels_var)),
      levels_var
    )

    # Shared fill scale across all panels
    all_color_vals <- sort(unique(na.omit(
      as.vector(physeq@tax_table[, rank_color])
    )))

    # Helper: subset physeq to a facet level, prune empty taxa
    subset_to_level <- function(lvl) {
      sam_values <- as.character(physeq@sam_data[[facet_by]])
      idx <- !is.na(sam_values) & sam_values == lvl
      pq_sub <- prune_samples(idx, physeq)
      prune_taxa(taxa_sums(pq_sub) > 0, pq_sub)
    }

    # Helper: add .cmpq_diff tag and build one bubble plot
    make_diff_panel <- function(pq, unique_taxa, contour_color, title) {
      phyloseq::tax_table(pq) <- cbind(
        phyloseq::tax_table(pq),
        .cmpq_diff = ifelse(
          phyloseq::taxa_names(pq) %in% unique_taxa,
          "unique",
          "shared"
        )
      )
      p <- gg_bubbles_pq(
        physeq = pq,
        rank_label = rank_label,
        rank_color = rank_color,
        rank_contour = ".cmpq_diff",
        layout = layout,
        log1ptransform = log1ptransform,
        min_nb_seq = min_nb_seq,
        label_size = label_size,
        label_color = label_color,
        show_labels = show_labels,
        border_width = diff_border_width,
        alpha = alpha,
        npoints = npoints,
        show_title = FALSE
      ) +
        ggplot2::scale_color_manual(
          values = c("unique" = contour_color, "shared" = "transparent"),
          guide = "none"
        )
      if (length(all_color_vals) > 0) {
        p <- p + ggplot2::scale_fill_discrete(limits = all_color_vals)
      }
      p + ggplot2::ggtitle(title)
    }

    # Build all pairwise panels
    pair_indices <- utils::combn(seq_along(levels_var), 2, simplify = FALSE)

    pair_plots <- lapply(pair_indices, function(idx) {
      i <- idx[1]
      j <- idx[2]
      lvl_i <- levels_var[i]
      lvl_j <- levels_var[j]

      pq_i <- subset_to_level(lvl_i)
      pq_j <- subset_to_level(lvl_j)

      unique_i <- setdiff(
        phyloseq::taxa_names(pq_i),
        phyloseq::taxa_names(pq_j)
      )
      unique_j <- setdiff(
        phyloseq::taxa_names(pq_j),
        phyloseq::taxa_names(pq_i)
      )

      p_i <- make_diff_panel(pq_i, unique_i, diff_colors[lvl_i], lvl_i)
      p_j <- make_diff_panel(pq_j, unique_j, diff_colors[lvl_j], lvl_j)

      (p_i | p_j) +
        patchwork::plot_annotation(
          title = paste0(lvl_i, " vs ", lvl_j),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(size = 11, face = "bold")
          )
        )
    })

    pw <- patchwork::wrap_plots(pair_plots, ncol = 1) +
      patchwork::plot_layout(guides = "collect")

    if (show_title) {
      size_desc <- if (log1ptransform) {
        "log1p(nb sequences)"
      } else {
        "Nb sequences"
      }
      label_desc <- if (rank_label == "Taxa") "taxa names" else rank_label
      title_parts <- paste0(
        "Fill: ",
        rank_color,
        " | Contour: unique taxa in each comparison",
        " | Size: ",
        size_desc,
        if (show_labels) paste0(" | Label: ", label_desc) else ""
      )
      pw <- pw +
        patchwork::plot_annotation(
          title = title_parts,
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(size = 10)
          )
        )
    }

    return(pw)
  }

  # ---- Build data frame (standard mode) --------------------------------------
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

    if (!is.null(rank_contour)) {
      df$rank_value_contour <- as.vector(pq@tax_table[, rank_contour])
    }

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

  # ---- Compute circle layout ------------------------------------------------
  compute_layout <- function(df) {
    lay <- packcircles::circleProgressiveLayout(
      df$value,
      sizetype = "area"
    )

    if (layout == "square") {
      cx <- mean(range(lay$x))
      cy <- mean(range(lay$y))
      dx <- lay$x - cx
      dy <- lay$y - cy
      r <- sqrt(dx^2 + dy^2)
      theta <- atan2(dy, dx)

      stretch <- ifelse(
        r > 0,
        1 / pmax(abs(cos(theta)), abs(sin(theta))),
        1
      )
      lay$x <- cx + dx * stretch
      lay$y <- cy + dy * stretch

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

    if (!is.null(rank_contour)) {
      vertices$rank_value_contour <- rep(
        df$rank_value_contour,
        each = npoints + 1
      )
    }

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

  # ---- Build plot ------------------------------------------------------------
  if (!is.null(rank_contour)) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data = verts,
        ggplot2::aes(
          x = x,
          y = y,
          group = id,
          fill = rank_value_color,
          color = rank_value_contour
        ),
        linewidth = border_width,
        alpha = alpha
      )
  } else {
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
      )
  }

  p <- p +
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

  if (show_title) {
    size_desc <- if (log1ptransform) {
      "log1p(nb sequences)"
    } else {
      "Nb sequences"
    }
    contour_part <- if (!is.null(rank_contour)) {
      paste0(" | Contour: ", rank_contour)
    } else {
      ""
    }
    label_desc <- if (rank_label == "Taxa") "taxa names" else rank_label
    title_parts <- paste0(
      "Fill: ",
      rank_color,
      contour_part,
      " | Size: ",
      size_desc,
      if (show_labels) paste0(" | Label: ", label_desc) else ""
    )
    p <- p + ggplot2::ggtitle(title_parts)
  }

  p
}
################################################################################
