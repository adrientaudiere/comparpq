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
#' When a [list_phyloseq] is passed as `physeq`, one bubble chart is produced
#' per phyloseq object and combined via \pkg{patchwork} with a shared legend.
#' If `diff_contour = TRUE` and there are at least 3 phyloseq objects, one
#' additional comparison panel is appended per pair, highlighting taxa unique
#' to each object in red.
#'
#' @inheritParams tc_points_matrix
#' @param rank_label (character, default "Taxa") The name of the column in the
#'   `@tax_table` slot to label the circles. If set to "Taxa", the taxa names
#'   are used.
#' @param rank_color (character, default "Family") The name of the column in
#'   the `@tax_table` slot to color the circles.
#' @param rank_contour (character, default NULL) The name of a column in the
#'   `@tax_table` slot to color the circle borders (contours). When NULL, the
#'   fixed `border_color` is used for all borders.
#' @param layout (character, default "circle") The packing layout. `"circle"`
#'   produces a tight circular packing. `"square"` constrains circles inside a
#'   square boundary, with large circles placed centrally.
#' @param facet_by (character, default NULL) A column name from `@sam_data` to
#'   facet the plot. When set, one bubble chart is produced per level of the
#'   variable, with taxa abundances computed within each level. Ignored when
#'   `physeq` is a [list_phyloseq].
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
#'   returned, but the resulting dataframe to plot is returned. Ignored when
#'   `physeq` is a [list_phyloseq].
#' @param diff_contour (logical, default FALSE) Only used when `physeq` is a
#'   [list_phyloseq] with at least 3 phyloseq objects. If TRUE, appends one
#'   pairwise comparison panel per pair of phyloseq objects, where taxa unique
#'   to each object (by taxa names) are highlighted with a red border contour.
#'   Matching taxa (present in both) receive a transparent contour. Taxa
#'   uniqueness is determined per pair by comparing `taxa_names()`.
#' @param diff_unique_color (character, default `"red"`) Border color for taxa
#'   that are unique to one object in a pairwise `diff_contour` comparison.
#' @param diff_border_width (numeric, default 1.5) Border width used in
#'   `diff_contour` comparison panels.
#'
#' @return A ggplot2 object (or patchwork when `physeq` is a [list_phyloseq]),
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
#' gg_bubbles_pq(
#'   physeq = data_fungi, rank_color = "Order",
#'   facet_by = "Height", min_nb_seq = 100
#' ) + no_legend()
#'
#' # list_phyloseq: one bubble chart per phyloseq, shared legend
#' lpq <- list_phyloseq(
#'   list(full = data_fungi, mini = data_fungi_mini)
#' )
#' gg_bubbles_pq(lpq, rank_color = "Class") & no_legend()
#'
#' # Highlight unique sequences when comparing two phyloseq objects.
#' # Here we contour taxa found only in data_fungi_mini,
#' # using transparent borders
#' # for shared taxa so that only unique ones stand out.
#'
#' lpq3 <- list_phyloseq(
#'   list(full = data_fungi, mini = data_fungi_mini, mini2 = mini2)
#' )
#' gg_bubbles_pq(lpq3, rank_color = "Class", diff_contour=T, show_labels=F,
#'  diff_border_width=1, diff_unique_color="#284028") & no_legend()
#'
#'
#' pq_list <- list_phyloseq(list("full" = data_fungi_mini, "mini" = mini2),
#'   same_bioinfo_pipeline = FALSE)
#' unique_seqs <-
#'   pq_list@comparison$refseq_comparison$full_vs_mini$unique_seqs_1
#'
#' tax_table(data_fungi_mini) <- cbind(
#'   tax_table(data_fungi_mini),
#'   unique_to_full = ifelse(
#'     as.character(data_fungi_mini@refseq) %in% unique_seqs, "Only_1", "both"
#'   )
#' )
#'
#' (gg_bubbles_pq(data_fungi_mini, rank_color = "Class", rank_contour = "unique_to_full", border_width = 1) +
#'   ggplot2::scale_color_manual(
#'     values = c("Only_1" = "red", "both" = "transparent")
#'   )) /
#'   gg_bubbles_pq(mini2, rank_color = "Class")
#'
gg_bubbles_pq <- function(
  physeq,
  rank_label = "Taxa",
  rank_color = "Family",
  rank_contour = NULL,
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
  return_dataframe = FALSE,
  diff_contour = FALSE,
  diff_unique_color = "red",
  diff_border_width = 1.5
) {
  # ---- Branch: list_phyloseq input -------------------------------------------
  if (inherits(physeq, "comparpq::list_phyloseq")) {
    rlang::check_installed(
      "patchwork",
      reason = "to combine bubble plots for list_phyloseq"
    )

    pq_list <- physeq
    n_pq <- length(pq_list@phyloseq_list)
    pq_names <- names(pq_list@phyloseq_list)

    # Collect all rank_color values across all objects for a shared legend
    all_color_vals <- sort(unique(na.omit(unlist(purrr::map(
      pq_list@phyloseq_list,
      ~ as.vector(.x@tax_table[, rank_color])
    )))))

    # Build individual plots with a shared fill scale
    individual_plots <- purrr::imap(pq_list@phyloseq_list, function(pq, name) {
      p <- gg_bubbles_pq(
        physeq = pq,
        rank_label = rank_label,
        rank_color = rank_color,
        rank_contour = rank_contour,
        layout = layout,
        log1ptransform = log1ptransform,
        min_nb_seq = min_nb_seq,
        label_size = label_size,
        label_color = label_color,
        show_labels = show_labels,
        border_color = border_color,
        border_width = border_width,
        alpha = alpha,
        npoints = npoints
      )
      if (length(all_color_vals) > 0) {
        p <- p + ggplot2::scale_fill_discrete(limits = all_color_vals)
      }
      p + ggplot2::ggtitle(name)
    })

    if (!diff_contour) {
      return(
        patchwork::wrap_plots(individual_plots) +
          patchwork::plot_layout(guides = "collect")
      )
    }

    if (!n_pq %in% c(2L, 3L)) {
      message(
        "`diff_contour` is only applied when the list_phyloseq contains ",
        "exactly 2 or 3 phyloseq objects. Returning individual plots only."
      )
      return(
        patchwork::wrap_plots(individual_plots) +
          patchwork::plot_layout(guides = "collect")
      )
    }

    # Helper: stamp ".cmpq_diff" onto a phyloseq tax_table
    add_diff_col <- function(pq, unique_taxa) {
      phyloseq::tax_table(pq) <- cbind(
        phyloseq::tax_table(pq),
        .cmpq_diff = ifelse(
          phyloseq::taxa_names(pq) %in% unique_taxa,
          "unique",
          "shared"
        )
      )
      pq
    }

    # Helper: one bubble plot with unique taxa outlined in diff_unique_color
    make_diff_half <- function(pq, title) {
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
        npoints = npoints
      ) +
        ggplot2::scale_color_manual(
          values = c("unique" = diff_unique_color, "shared" = "transparent"),
          guide = "none"
        )
      if (length(all_color_vals) > 0) {
        p <- p + ggplot2::scale_fill_discrete(limits = all_color_vals)
      }
      p + ggplot2::ggtitle(title)
    }

    # n = 2: replace the two individual plots with diff-highlighted versions
    if (n_pq == 2L) {
      pq_1 <- pq_list@phyloseq_list[[1]]
      pq_2 <- pq_list@phyloseq_list[[2]]
      unique_1 <- setdiff(
        phyloseq::taxa_names(pq_1),
        phyloseq::taxa_names(pq_2)
      )
      unique_2 <- setdiff(
        phyloseq::taxa_names(pq_2),
        phyloseq::taxa_names(pq_1)
      )

      p1 <- make_diff_half(add_diff_col(pq_1, unique_1), pq_names[1])
      p2 <- make_diff_half(add_diff_col(pq_2, unique_2), pq_names[2])

      return(
        patchwork::wrap_plots(list(p1, p2)) +
          patchwork::plot_layout(guides = "collect")
      )
    }

    # n = 3: three pairwise comparison panels stacked in a single column
    pair_indices <- utils::combn(seq_len(n_pq), 2, simplify = FALSE)

    diff_panels <- purrr::map(pair_indices, function(idx) {
      i <- idx[1]
      j <- idx[2]
      pq_i <- pq_list@phyloseq_list[[i]]
      pq_j <- pq_list@phyloseq_list[[j]]
      name_i <- pq_names[i]
      name_j <- pq_names[j]

      unique_i <- setdiff(
        phyloseq::taxa_names(pq_i),
        phyloseq::taxa_names(pq_j)
      )
      unique_j <- setdiff(
        phyloseq::taxa_names(pq_j),
        phyloseq::taxa_names(pq_i)
      )

      p_i <- make_diff_half(
        add_diff_col(pq_i, unique_i),
        paste0(name_i, "\n(unique vs ", name_j, ")")
      )
      p_j <- make_diff_half(
        add_diff_col(pq_j, unique_j),
        paste0(name_j, "\n(unique vs ", name_i, ")")
      )

      (p_i | p_j) +
        patchwork::plot_annotation(
          title = paste0(name_i, " vs ", name_j),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(size = 11, face = "bold")
          )
        )
    })

    diff_combined <- patchwork::wrap_plots(diff_panels, ncol = 1) +
      patchwork::plot_layout(guides = "collect")

    return(diff_combined)
  }

  # ---- Single phyloseq input (original logic) --------------------------------
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

  p
}
################################################################################
