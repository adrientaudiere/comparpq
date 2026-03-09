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
#' @param physeq (phyloseq or list_phyloseq, required) A phyloseq
#'   object, or a [list_phyloseq] object. When a list_phyloseq is
#'   provided, it is first merged into a single phyloseq using
#'   [merge_lpq()] (each original phyloseq becomes one sample) and
#'   the `fact` parameter is automatically set to `"source_name"`.
#' @param fact (character, required when `physeq` is a phyloseq)
#'   Name of a variable in `sample_data(physeq)` defining the groups
#'   (2-4 levels). Ignored when `physeq` is a list_phyloseq.
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
#' @param show_na_count (logical, default FALSE) If `TRUE`, display
#'   the number of taxa with `NA` at the chosen `taxonomic_rank` in
#'   the bottom-left corner of the plot. When `count_type = "taxa"`,
#'   the sum of all Venn region counts plus the NA count equals
#'   `ntaxa(physeq)`. Ignored when `taxonomic_rank` is `NULL`.
#' @param count_taxa (logical, default TRUE) If `TRUE`, append a
#'   `"Taxa"` panel to the Venn diagram showing shared and unique
#'   individual taxa (ASVs/OTUs) alongside the aggregated taxonomic
#'   ranks. A temporary `Taxa` column is added to the tax_table
#'   with each taxon's name as its value. Ignored when
#'   `taxonomic_rank` is `NULL`.
#' @param match_by (character, default `"refseq"`) Passed to
#'   [merge_lpq()] when `physeq` is a list_phyloseq. One of
#'   `"refseq"` or `"names"`.
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
#'
#' # From a list_phyloseq object
#' lpq <- list_phyloseq(list(
#'   fungi = data_fungi_mini,
#'   fungi2 = data_fungi
#' ))
#' simple_venn_pq(lpq, taxonomic_rank = "Genus", count_taxa)
simple_venn_pq <- function(
  physeq,
  fact = NULL,
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
  show_na_count = FALSE,
  count_taxa = TRUE,
  match_by = c("refseq", "names"),
  combine = TRUE,
  verbose = TRUE,
  .lpq_n_samples = NULL
) {
  # Handle list_phyloseq input
  lpq_n_samples <- .lpq_n_samples
  if (inherits(physeq, "comparpq::list_phyloseq")) {
    match_by <- match.arg(match_by)
    # Capture original sample counts before merging
    lpq_n_samples <- vapply(
      physeq@phyloseq_list,
      phyloseq::nsamples,
      integer(1)
    )
    physeq <- merge_lpq(physeq, match_by = match_by, verbose = verbose)
    fact <- "source_name"
  } else if (is.null(fact)) {
    stop("'fact' is required when 'physeq' is a phyloseq object.")
  }

  # Add taxa-level count as a pseudo-rank
  if (count_taxa && !is.null(taxonomic_rank) && inherits(physeq, "phyloseq")) {
    tt <- as.data.frame(phyloseq::tax_table(physeq))
    tt[["Taxa (OTU/ASV)"]] <- phyloseq::taxa_names(physeq)
    phyloseq::tax_table(physeq) <- phyloseq::tax_table(as.matrix(tt))
    taxonomic_rank <- c(taxonomic_rank, "Taxa (OTU/ASV)")
  }

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
        show_na_count = show_na_count,
        count_taxa = FALSE,
        match_by = match_by,
        combine = combine,
        verbose = verbose,
        .lpq_n_samples = lpq_n_samples
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
    stop("'physeq' must be a phyloseq or list_phyloseq object.")
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

  if (n_groups < 2) {
    stop(
      "simple_venn_pq() requires at least 2 groups. '",
      fact,
      "' has ",
      n_groups,
      " unique level(s)."
    )
  }
  if (n_groups > 4) {
    stop(
      "simple_venn_pq() supports at most 4 groups. '",
      fact,
      "' has ",
      n_groups,
      " unique level(s). Consider using MiscMetabar::upset_pq() for more groups."
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
        names(taxa_sums) <- rownames(otu)
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
  if (!is.null(lpq_n_samples)) {
    n_samples <- lpq_n_samples[levels_fact]
  } else {
    n_samples <- vapply(
      levels_fact,
      \(lev) sum(groups == lev),
      integer(1)
    )
  }

  # Group labels
  group_labels <- levels_fact
  sample_labels <- if (add_nb_samples) {
    paste0("(n=", n_samples, ")")
  } else {
    NULL
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

  # Define shapes and look up precomputed centroids
  # Shrink shapes when sample counts are shown to leave room for labels
  shapes <- venn_shapes(n_groups)
  if (add_nb_samples) {
    shapes <- lapply(shapes, \(s) {
      s$cx <- s$cx * 0.85
      s$cy <- s$cy * 0.85
      s$a <- s$a * 0.85
      s$b <- s$b * 0.85
      s
    })
  }
  centroids <- venn_centroids_precomputed(n_groups, add_nb_samples)

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
    if (add_nb_samples) {
      p <- p +
        ggplot2::annotate(
          "text",
          x = label_pos[[i]][1],
          y = label_pos[[i]][2] - 0.22,
          label = sample_labels[i],
          size = label_size * 0.7,
          color = "grey50"
        )
    }
  }

  # NA count annotation
  if (show_na_count && !is.null(taxonomic_rank)) {
    if (count_type == "taxa" || count_type == "rank_taxa") {
      # Count excluded taxa so Venn + NA = ntaxa(physeq_original)
      na_count <- phyloseq::ntaxa(physeq_original) - sum(all_counts)
    } else {
      tt_na <- as.data.frame(phyloseq::tax_table(physeq_original))
      na_count <- sum(
        is.na(tt_na[[taxonomic_rank]]) |
          tt_na[[taxonomic_rank]] == ""
      )
    }
    # Position in bottom-left, outside the shapes
    na_x <- if (n_groups == 4) -1.8 else -1.5
    na_y <- if (n_groups == 3) -1.5 else -1.3
    if (add_nb_samples) {
      na_x <- na_x * 0.85
      na_y <- na_y * 0.85
    }
    p <- p +
      ggplot2::annotate(
        "text",
        x = na_x,
        y = na_y,
        label = paste0("NA: ", na_count),
        size = text_size * 0.9,
        color = "grey40",
        fontface = "italic",
        hjust = 0
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
    title_hjust <- if (n_groups == 3) 0.05 else 0.5
    p <- p +
      ggplot2::labs(title = taxonomic_rank, caption = count_label) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 14,
          face = "bold",
          hjust = title_hjust
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


#' Look up precomputed centroids for standard Venn layouts
#'
#' Returns centroids for the standard 2-, 3-, or 4-set Venn shapes,
#' with or without the 0.85 scaling applied for sample-count labels.
#' These were computed once via `venn_centroids()` with n_grid = 200
#' and hardcoded to avoid repeated grid sampling.
#' @param n_sets Integer, number of sets (2, 3, or 4).
#' @param scaled Logical, whether the 0.85 scaling is applied.
#' @return Named list of numeric(2) centroids, keyed by binary code.
#' @noRd
venn_centroids_precomputed <- function(n_sets, scaled = FALSE) {
  key <- paste0(n_sets, "_", if (scaled) "s" else "u")
  .venn_centroid_cache[[key]]
}

# Precomputed centroid positions for all standard Venn configurations.
# Keys: "{n_sets}_{u|s}" where u = unscaled, s = scaled (×0.85).
.venn_centroid_cache <- list(
  "2_u" = list(
    "01" = c(1.003293, -0.005527),
    "10" = c(-1.003156, -0.005527),
    "11" = c(-0.008223, -0.005527)
  ),
  "2_s" = list(
    "01" = c(0.852799, -0.004698),
    "10" = c(-0.852683, -0.004698),
    "11" = c(-0.006989, -0.004698)
  ),
  "3_u" = list(
    "001" = c(0.953547, -0.546447),
    "010" = c(-0.953409, -0.546447),
    "011" = c(0.008360, -0.736815),
    "100" = c(-0.008223, 1.103407),
    "101" = c(0.638484, 0.373664),
    "110" = c(-0.638347, 0.373664),
    "111" = c(-0.008223, -0.007072)
  ),
  "3_s" = list(
    "001" = c(0.810515, -0.464480),
    "010" = c(-0.810398, -0.464480),
    "011" = c(0.007106, -0.626293),
    "100" = c(-0.006989, 0.937896),
    "101" = c(0.542712, 0.317614),
    "110" = c(-0.542595, 0.317614),
    "111" = c(-0.006989, -0.006011)
  ),
  "4_u" = list(
    "0001" = c(1.321071, 0.471399),
    "0010" = c(0.404464, -1.055556),
    "0011" = c(0.098929, -0.922777),
    "0100" = c(-0.404306, 1.055626),
    "0101" = c(1.105399, 0.179286),
    "0110" = c(0.925672, -0.710331),
    "0111" = c(1.015535, -0.298717),
    "1000" = c(-1.320912, -0.471329),
    "1001" = c(-0.763759, -0.763442),
    "1010" = c(-1.105240, -0.179216),
    "1011" = c(-0.404306, -0.803276),
    "1100" = c(-0.098770, 0.922847),
    "1101" = c(0.404464, 0.803346),
    "1110" = c(-1.033350, 0.285509),
    "1111" = c(-0.008907, -0.033160)
  ),
  "4_s" = list(
    "0001" = c(1.122910, 0.400689),
    "0010" = c(0.343795, -0.897222),
    "0011" = c(0.084090, -0.784360),
    "0100" = c(-0.343660, 0.897282),
    "0101" = c(0.939589, 0.152393),
    "0110" = c(0.786821, -0.603781),
    "0111" = c(0.863205, -0.253909),
    "1000" = c(-1.122775, -0.400630),
    "1001" = c(-0.649195, -0.648926),
    "1010" = c(-0.939454, -0.152334),
    "1011" = c(-0.343660, -0.682785),
    "1100" = c(-0.083955, 0.784420),
    "1101" = c(0.343795, 0.682844),
    "1110" = c(-0.878347, 0.242683),
    "1111" = c(-0.007571, -0.028186)
  )
)


#' Compute centroids of each Venn region via grid sampling
#'
#' For each region, finds the grid point with the highest "depth
#' score": how deeply inside the included sets and how far outside
#' the excluded sets. This avoids placing labels at the center for
#' symmetric regions that span both sides of the diagram.
#' @param shapes List of shape parameter lists.
#' @param n_sets Integer, number of sets.
#' @param n_grid Integer, grid resolution per axis.
#' @return Named list of numeric(2) centroids, keyed by binary code.
#' @noRd
venn_centroids <- function(shapes, n_sets, n_grid = 200) {
  all_pts <- do.call(rbind, lapply(shapes, ellipse_polygon))
  xlim <- range(all_pts$x) * 1.1
  ylim <- range(all_pts$y) * 1.1

  gx <- seq(xlim[1], xlim[2], length.out = n_grid)
  gy <- seq(ylim[1], ylim[2], length.out = n_grid)
  px <- rep(gx, n_grid)
  py <- rep(gy, each = n_grid)
  n_pts <- length(px)

  # Compute normalized ellipse distance for each shape
  # (< 1 = inside, > 1 = outside)
  ellipse_dist <- matrix(0, nrow = n_pts, ncol = n_sets)
  for (s in seq_len(n_sets)) {
    sh <- shapes[[s]]
    dx <- px - sh$cx
    dy <- py - sh$cy
    cos_a <- cos(sh$angle)
    sin_a <- sin(sh$angle)
    lx <- dx * cos_a + dy * sin_a
    ly <- -dx * sin_a + dy * cos_a
    ellipse_dist[, s] <- sqrt((lx / sh$a)^2 + (ly / sh$b)^2)
  }

  # Encode region as integer (avoid slow apply+paste)
  membership <- ellipse_dist <= 1
  powers <- 10L^((n_sets - 1L):0L)
  codes_int <- as.integer(membership %*% powers)

  # All possible non-empty codes
  combos <- expand.grid(rep(list(0:1), n_sets))
  possible_int <- as.integer(as.matrix(combos) %*% powers)
  possible_int <- possible_int[possible_int > 0L]

  centroids <- list()
  for (code_int in possible_int) {
    idx <- which(codes_int == code_int)
    if (length(idx) == 0L) {
      next
    }
    # Decode bits
    bits <- as.integer(strsplit(
      formatC(code_int, width = n_sets, flag = "0"),
      ""
    )[[1]])
    # Depth score: min across sets of distance to boundary
    depth <- rep(Inf, length(idx))
    for (s in seq_len(n_sets)) {
      d <- ellipse_dist[idx, s]
      if (bits[s] == 1L) {
        depth <- pmin(depth, 1 - d)
      } else {
        depth <- pmin(depth, d - 1)
      }
    }
    best <- which.max(depth)
    code_str <- formatC(code_int, width = n_sets, flag = "0")
    centroids[[code_str]] <- c(x = px[idx[best]], y = py[idx[best]])
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
