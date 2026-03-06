utils::globalVariables(c("x", "y"))

#' Venn diagram of shared taxa across sample groups
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#'   alt="lifecycle-experimental"></a>
#'
#' Draws a Venn diagram showing shared and unique taxa (or higher-rank
#' groups) across 2 to 4 levels of a sample variable, using only
#' ggplot2 (no external Venn diagram package needed).
#'
#' When `taxonomic_rank` is a character vector of length > 1 (the
#' default), all ranks are displayed in a single combined figure using
#' the patchwork package (must be installed). Set `combine = FALSE` to
#' get a named list of individual plots instead.
#'
#' For 2 and 3 groups, circles are used. For 4 groups, ellipses are
#' used to ensure all intersection regions are representable.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param fact (character, required) Name of a variable in
#'   `sample_data(physeq)` defining the groups (2-4 levels).
#' @param min_nb_seq (integer, default 0) Minimum total read count for
#'   a taxon to be considered present in a group. A taxon must have
#'   strictly more than `min_nb_seq` reads in a group to be included.
#' @param taxonomic_rank (character or NULL) Taxonomic rank(s) at
#'   which to aggregate (via [phyloseq::tax_glom()]) before computing
#'   the Venn diagram. Defaults to all standard ranks (Kingdom through
#'   Species). Use `NULL` to skip aggregation and work at ASV/OTU
#'   level.
#' @param na_remove (logical, default TRUE) Remove samples with NA in
#'   `fact` and, when aggregating, taxa with NA at `taxonomic_rank`.
#' @param add_nb_samples (logical, default TRUE) Append sample count
#'   to group labels.
#' @param fill_alpha (numeric, default 0.3) Fill transparency for
#'   shapes.
#' @param border_size (numeric, default 0.8) Border line width.
#' @param count_type (character, default `"rank"`) What to count in
#'   each Venn region. One of:
#'   - `"rank"`: number of unique taxonomic levels (e.g. number of
#'     shared Classes). This is the default.
#'   - `"taxa"`: number of ASVs/OTUs assigned to the shared
#'     taxonomic levels.
#'   - `"sequences"`: total number of reads for ASVs/OTUs assigned
#'     to the shared taxonomic levels.
#'   - `"rank_taxa"`: shows both rank and taxa counts as
#'     `"nb_rank (nb_taxa)"`.
#'   Ignored when `taxonomic_rank` is `NULL` (ASV-level), where
#'   `"rank"` and `"taxa"` are equivalent.
#' @param text_size (numeric, default 4) Base size of count labels
#'   inside regions.
#' @param scale_text (logical, default FALSE) If `TRUE`, scale the
#'   size of count labels proportionally to the count value. The
#'   `text_size` parameter then acts as the base (minimum) size.
#' @param hide_zero (logical, default TRUE) If `TRUE`, hide count
#'   labels that are zero (or `"0 (0)"` when `count_type = "rank_taxa"`).
#' @param label_size (numeric, default 4.5) Size of group name labels.
#' @param colors (character or NULL) Vector of colors, one per group.
#'   Defaults to a 4-color qualitative palette.
#' @param combine (logical, default TRUE) When `taxonomic_rank` has
#'   length > 1, combine plots into a single patchwork figure. Set to
#'   `FALSE` to return a named list of individual ggplot objects.
#'   Requires the patchwork package.
#' @param verbose (logical, default TRUE) Print a message when no taxa
#'   meet the criteria.
#'
#' @return A ggplot2 object (single rank), a patchwork object (multiple
#'   ranks with `combine = TRUE`), or a named list of ggplot2 objects
#'   (multiple ranks with `combine = FALSE`).
#'
#' @export
#' @author Adrien Taudiere
#'
#' @examples
#' # Default: all ranks combined in one figure
#' simple_venn_pq(data_fungi_mini, "Height")
#'
#' # At genus level only
#' simple_venn_pq(data_fungi_mini, "Height", taxonomic_rank = "Genus")
#'
#' # Multiple ranks as a list
#' plots <- simple_venn_pq(
#'   data_fungi_mini, "Height",
#'   taxonomic_rank = c("Family", "Genus"),
#'   combine = FALSE
#' )
#' plots[["Family"]]
#'
#' # Count ASVs instead of rank levels
#' simple_venn_pq(data_fungi_mini, "Height",
#'   taxonomic_rank = "Genus", count_type = "taxa"
#' )
#'
#' # Scale text by count value
#' simple_venn_pq(data_fungi_mini, "Height",
#'   taxonomic_rank = "Genus", scale_text = TRUE
#' )
simple_venn_pq <- function(
  physeq,
  fact,
  min_nb_seq = 0,
  taxonomic_rank = c(
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  ),
  na_remove = TRUE,
  count_type = c("rank", "taxa", "sequences", "rank_taxa"),
  add_nb_samples = TRUE,
  fill_alpha = 0.3,
  border_size = 0.8,
  text_size = 4,
  scale_text = FALSE,
  hide_zero = TRUE,
  label_size = 4.5,
  colors = NULL,
  combine = TRUE,
  verbose = TRUE
) {
  # Multiple ranks -> combine or return list of plots
  if (length(taxonomic_rank) > 1) {
    plots <- lapply(taxonomic_rank, \(rank) {
      simple_venn_pq(
        physeq,
        fact,
        min_nb_seq = min_nb_seq,
        taxonomic_rank = rank,
        na_remove = na_remove,
        count_type = count_type,
        add_nb_samples = add_nb_samples,
        fill_alpha = fill_alpha,
        border_size = border_size,
        text_size = text_size,
        scale_text = scale_text,
        hide_zero = hide_zero,
        label_size = label_size,
        colors = colors,
        combine = combine,
        verbose = verbose
      )
    })
    names(plots) <- taxonomic_rank
    # Remove NULL entries (ranks with no taxa)
    plots <- Filter(Negate(is.null), plots)
    if (combine) {
      if (!requireNamespace("patchwork", quietly = TRUE)) {
        stop(
          "The 'patchwork' package is required when combine = TRUE. ",
          "Install it with install.packages('patchwork') or use ",
          "combine = FALSE.",
          call. = FALSE
        )
      }
      # Add light grey border around each panel
      plots <- lapply(plots, \(pl) {
        pl +
          ggplot2::theme(
            plot.background = ggplot2::element_rect(
              color = "grey80",
              fill = NA,
              linewidth = 0.5
            ),
            plot.margin = ggplot2::margin(12, 12, 12, 12)
          )
      })
      combined <- patchwork::wrap_plots(plots) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10))
      return(combined)
    }
    return(plots)
  }

  count_type <- match.arg(count_type)

  if (!inherits(physeq, "phyloseq")) {
    stop("'physeq' must be a phyloseq object.")
  }

  sam <- phyloseq::sample_data(physeq)
  if (!fact %in% colnames(sam)) {
    stop("'", fact, "' not found in sample_data.")
  }

  groups <- as.character(sam[[fact]])

  if (na_remove) {
    keep <- !is.na(groups)
    if (sum(keep) < phyloseq::nsamples(physeq)) {
      physeq <- phyloseq::prune_samples(keep, physeq)
      groups <- groups[keep]
    }
  }

  levels_fact <- unique(groups)
  n_groups <- length(levels_fact)

  if (n_groups < 2 || n_groups > 4) {
    stop(
      "simple_venn_pq() requires 2 to 4 groups. '",
      fact,
      "' has ",
      n_groups,
      " unique level(s)."
    )
  }

  # Keep original data for taxa/sequences counting
  physeq_original <- physeq

  # OTU table before glom (rows = taxa)
  otu_orig <- as.data.frame(phyloseq::otu_table(physeq_original))
  if (!phyloseq::taxa_are_rows(physeq_original)) {
    otu_orig <- as.data.frame(t(otu_orig))
  }

  # Aggregate at taxonomic rank if requested
  if (!is.null(taxonomic_rank)) {
    physeq <- phyloseq::tax_glom(
      physeq,
      taxrank = taxonomic_rank,
      NArm = na_remove
    )
  }

  # Build sets: taxa with > min_nb_seq reads per group
  otu <- as.data.frame(phyloseq::otu_table(physeq))
  if (!phyloseq::taxa_are_rows(physeq)) {
    otu <- as.data.frame(t(otu))
  }

  sets <- lapply(
    stats::setNames(levels_fact, levels_fact),
    \(lev) {
      samp_idx <- which(groups == lev)
      if (length(samp_idx) == 1) {
        taxa_sums <- otu[, samp_idx]
      } else {
        taxa_sums <- rowSums(otu[, samp_idx])
      }
      names(which(taxa_sums > min_nb_seq))
    }
  )

  # Map to rank values for display clarity
  if (!is.null(taxonomic_rank)) {
    tt <- as.data.frame(phyloseq::tax_table(physeq))
    rank_map <- stats::setNames(tt[[taxonomic_rank]], rownames(tt))
    sets <- lapply(sets, \(s) unique(rank_map[s]))
  }

  # Sample counts per group
  n_samples <- vapply(
    levels_fact,
    \(lev) sum(groups == lev),
    integer(1)
  )

  # Group labels
  if (add_nb_samples) {
    group_labels <- paste0(levels_fact, "\n(n=", n_samples, ")")
  } else {
    group_labels <- levels_fact
  }

  # Compute region memberships
  all_items <- unique(unlist(sets))
  if (length(all_items) == 0) {
    if (verbose) {
      message("No taxa found meeting criteria.")
    }
    return(invisible(NULL))
  }

  membership <- vapply(
    sets,
    \(s) all_items %in% s,
    logical(length(all_items))
  )
  codes <- apply(
    membership,
    1,
    \(row) paste(as.integer(row), collapse = "")
  )

  # Compute region counts based on count_type
  region_rank_counts <- as.list(table(codes))

  if (
    count_type %in%
      c("taxa", "sequences", "rank_taxa") &&
      !is.null(taxonomic_rank)
  ) {
    # Map original ASVs to their rank value
    tt_orig <- as.data.frame(phyloseq::tax_table(physeq_original))
    asv_rank <- tt_orig[[taxonomic_rank]]
    names(asv_rank) <- rownames(tt_orig)
    # For each Venn code, find rank values in that region
    item_codes <- stats::setNames(codes, all_items)
    unique_codes <- unique(codes)
    region_taxa_counts <- lapply(
      stats::setNames(unique_codes, unique_codes),
      \(cd) {
        rank_vals <- names(item_codes[item_codes == cd])
        asv_idx <- which(asv_rank %in% rank_vals)
        if (count_type == "sequences") {
          sum(otu_orig[asv_idx, , drop = FALSE])
        } else {
          # "taxa" or "rank_taxa"
          length(asv_idx)
        }
      }
    )
  }

  # Build region_counts depending on count_type
  if (count_type == "rank" || is.null(taxonomic_rank)) {
    region_counts <- region_rank_counts
  } else if (count_type == "rank_taxa") {
    region_counts <- region_rank_counts
    region_counts_secondary <- region_taxa_counts
  } else {
    region_counts <- region_taxa_counts
  }

  # Define shapes and compute centroids
  shapes <- venn_shapes(n_groups)
  centroids <- venn_centroids(shapes, n_groups)

  # Default colors
  if (is.null(colors)) {
    default_colors <- c("#44AA99", "#CC6677", "#DDCC77", "#88CCEE")
    colors <- default_colors[seq_len(n_groups)]
  }

  # Build ggplot
  p <- ggplot2::ggplot() +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 10, 20, 10))

  # Draw shapes
  for (i in seq_len(n_groups)) {
    shape_df <- ellipse_polygon(shapes[[i]])
    p <- p +
      ggplot2::geom_polygon(
        data = shape_df,
        ggplot2::aes(x = x, y = y),
        fill = colors[i],
        alpha = fill_alpha,
        color = colors[i],
        linewidth = border_size
      )
  }

  # Add region count labels
  all_counts <- vapply(
    names(centroids),
    \(code) {
      if (code %in% names(region_counts)) {
        as.integer(region_counts[[code]])
      } else {
        0L
      }
    },
    integer(1)
  )
  max_count <- max(all_counts, 1L)

  # Build label strings
  if (count_type == "rank_taxa" && !is.null(taxonomic_rank)) {
    all_labels <- vapply(
      names(centroids),
      \(code) {
        primary <- all_counts[[code]]
        secondary <- if (code %in% names(region_counts_secondary)) {
          as.integer(region_counts_secondary[[code]])
        } else {
          0L
        }
        paste0(primary, " (", secondary, ")")
      },
      character(1)
    )
  } else {
    all_labels <- as.character(all_counts)
    names(all_labels) <- names(centroids)
  }

  for (code in names(centroids)) {
    count <- all_counts[[code]]
    if (hide_zero && count == 0L) {
      next
    }
    pos <- centroids[[code]]
    if (scale_text) {
      scaled_size <- text_size + text_size * (count / max_count)
    } else {
      scaled_size <- text_size
    }
    p <- p +
      ggplot2::annotate(
        "text",
        x = pos[1],
        y = pos[2],
        label = all_labels[[code]],
        size = scaled_size,
        fontface = "bold"
      )
  }

  # Add group name labels (darkened colors for readability)
  label_colors <- vapply(
    colors,
    \(col) {
      rgb_val <- grDevices::col2rgb(col)[, 1] / 255
      grDevices::rgb(rgb_val[1] * 0.55, rgb_val[2] * 0.55, rgb_val[3] * 0.55)
    },
    character(1)
  )
  label_pos <- venn_label_positions(n_groups)
  for (i in seq_len(n_groups)) {
    p <- p +
      ggplot2::annotate(
        "text",
        x = label_pos[[i]][1],
        y = label_pos[[i]][2],
        label = group_labels[i],
        size = label_size,
        fontface = "bold",
        color = label_colors[i]
      )
  }

  # Title
  if (!is.null(taxonomic_rank)) {
    count_label <- switch(
      count_type,
      rank = paste("nb.", taxonomic_rank),
      taxa = "nb. taxa",
      sequences = "nb. sequences",
      rank_taxa = paste0("nb. ", taxonomic_rank, " (nb. taxa)")
    )
    p <- p +
      ggplot2::labs(title = taxonomic_rank, caption = count_label) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 14,
          face = "bold",
          hjust = 0.05
        ),
        plot.caption = ggplot2::element_text(
          size = 10,
          hjust = 0.5,
          color = "grey40"
        )
      )
  }
  p
}


# --- Internal helpers ---------------------------------------------------

#' Define shape geometry for Venn layouts
#'
#' Returns a list of shape parameter lists (cx, cy, a, b, angle).
#' Circles have a == b and angle == 0.
#' @param n Integer, number of sets (2, 3, or 4).
#' @return List of shape parameter lists.
#' @noRd
venn_shapes <- function(n) {
  if (n == 2) {
    list(
      list(cx = -0.5, cy = 0, a = 1, b = 1, angle = 0),
      list(cx = 0.5, cy = 0, a = 1, b = 1, angle = 0)
    )
  } else if (n == 3) {
    list(
      list(cx = 0, cy = 0.58, a = 1, b = 1, angle = 0),
      list(cx = -0.5, cy = -0.29, a = 1, b = 1, angle = 0),
      list(cx = 0.5, cy = -0.29, a = 1, b = 1, angle = 0)
    )
  } else {
    # 4-set Venn with ellipses (based on VennDiagram layout)
    list(
      list(cx = -0.30, cy = -0.06, a = 1.44, b = 0.90, angle = pi / 6),
      list(cx = 0.00, cy = 0.14, a = 1.44, b = 0.90, angle = 5 * pi / 6),
      list(cx = 0.00, cy = -0.14, a = 1.44, b = 0.90, angle = 5 * pi / 6),
      list(cx = 0.30, cy = 0.06, a = 1.44, b = 0.90, angle = pi / 6)
    )
  }
}


#' Generate polygon vertices for an ellipse
#' @param shape List with cx, cy, a, b, angle.
#' @param n Number of vertices.
#' @return data.frame with x and y columns.
#' @noRd
ellipse_polygon <- function(shape, n = 200) {
  theta <- seq(0, 2 * pi, length.out = n)
  lx <- shape$a * cos(theta)
  ly <- shape$b * sin(theta)
  cos_a <- cos(shape$angle)
  sin_a <- sin(shape$angle)
  data.frame(
    x = shape$cx + lx * cos_a - ly * sin_a,
    y = shape$cy + lx * sin_a + ly * cos_a
  )
}


#' Test whether points lie inside an ellipse
#' @param px,py Numeric vectors of point coordinates.
#' @param shape List with cx, cy, a, b, angle.
#' @return Logical vector.
#' @noRd
point_in_ellipse <- function(px, py, shape) {
  dx <- px - shape$cx
  dy <- py - shape$cy
  cos_a <- cos(shape$angle)
  sin_a <- sin(shape$angle)
  lx <- dx * cos_a + dy * sin_a
  ly <- -dx * sin_a + dy * cos_a
  (lx / shape$a)^2 + (ly / shape$b)^2 <= 1
}


#' Compute centroids of each Venn region via grid sampling
#' @param shapes List of shape parameter lists.
#' @param n_sets Integer, number of sets.
#' @param n_grid Integer, grid resolution per axis.
#' @return Named list of numeric(2) centroids, keyed by binary code.
#' @noRd
venn_centroids <- function(shapes, n_sets, n_grid = 500) {
  all_pts <- do.call(rbind, lapply(shapes, ellipse_polygon))
  xlim <- range(all_pts$x) * 1.1
  ylim <- range(all_pts$y) * 1.1

  xseq <- seq(xlim[1], xlim[2], length.out = n_grid)
  yseq <- seq(ylim[1], ylim[2], length.out = n_grid)
  grid <- expand.grid(x = xseq, y = yseq)

  membership <- vapply(
    shapes,
    \(s) point_in_ellipse(grid$x, grid$y, s),
    logical(nrow(grid))
  )

  codes <- apply(
    membership,
    1,
    \(row) paste(as.integer(row), collapse = "")
  )

  all_zero <- paste(rep("0", n_sets), collapse = "")

  # All possible non-empty codes
  combos <- expand.grid(rep(list(0:1), n_sets))
  possible <- apply(combos, 1, \(row) paste(row, collapse = ""))
  possible <- setdiff(possible, all_zero)

  # Use spatial median (medoid) instead of mean for centroid,
  # which is more robust for thin/curved regions
  centroids <- list()
  for (code in possible) {
    idx <- which(codes == code)
    if (length(idx) > 0) {
      centroids[[code]] <- c(
        x = stats::median(grid$x[idx]),
        y = stats::median(grid$y[idx])
      )
    }
  }

  centroids
}


#' Fixed label positions for group names (outside the shapes)
#' @param n Integer, number of sets.
#' @return List of numeric(2) positions.
#' @noRd
venn_label_positions <- function(n) {
  if (n == 2) {
    list(c(-1.1, 1.1), c(1.1, 1.1))
  } else if (n == 3) {
    list(c(0, 2.0), c(-1.3, -1.0), c(1.3, -1.0))
  } else {
    list(c(-1.7, -0.7), c(-0.8, 1.5), c(0.8, -1.5), c(1.7, 0.7))
  }
}
