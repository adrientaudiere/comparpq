################################################################################
#' Compute congruence metrics between two taxonomic assignments
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Computes metrics quantifying the congruence between two taxonomic assignments
#' for the same set of taxa (ASVs/OTUs). This is useful for comparing different
#' taxonomic databases, assignment methods, or reference versions.
#'
#' @param physeq_1 (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   with the first taxonomic assignment.
#' @param physeq_2 (phyloseq, default NULL) A \code{\link[phyloseq]{phyloseq-class}}
#'   object with the second taxonomic assignment. If NULL, uses physeq_1
#'   (useful for comparing different rank columns from the same object).
#' @param ranks_1 (character vector, required) Taxonomic rank names to use for
#'   the first assignment. Must match column names in physeq_1's tax_table.
#' @param ranks_2 (character vector, default NULL) Taxonomic rank names to use for
#'   the second assignment. Must match column names in physeq_2's tax_table. If
#'   NULL, uses the same ranks as ranks_1 vector.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{summary}{A data frame with counts and percentages for each category,
#'     including `leaf_match_congruence` which counts taxa where the deepest
#'     classification matches regardless of the path}
#'   \item{only_in_1}{Character vector of taxa names present only in physeq_1}
#'   \item{only_in_2}{Character vector of taxa names present only in physeq_2}
#'   \item{classified_only_1}{Character vector of common taxa classified only
#'     in physeq_1 (all NA in physeq_2 for the specified ranks)}
#'   \item{classified_only_2}{Character vector of common taxa classified only
#'     in physeq_2 (all NA in physeq_1 for the specified ranks)}
#'   \item{unclassified_both}{Character vector of common taxa with all NA in
#'     both physeq objects}
#'   \item{total_congruent}{Character vector of taxa with identical taxonomic paths
#'     (same values at all ranks, same NA positions)}
#'   \item{partial_1_deeper}{Character vector of taxa where physeq_1 classifies
#'     to a deeper rank but agrees on shared ranks}
#'   \item{partial_2_deeper}{Character vector of taxa where physeq_2 classifies
#'     to a deeper rank but agrees on shared ranks}
#'   \item{incongruent_leaves}{Character vector of taxa with disagreement at
#'     the deepest (leaf) level where both have assignments}
#'   \item{incongruent_nodes}{Character vector of taxa with disagreement at
#'     higher (internal node) levels, even if leaves match}
#'   \item{details}{A data frame with per-taxon details including: taxon name,
#'     depth_1, depth_2, leaf_1, leaf_2, leaf_match (TRUE if leaves are equal),
#'     and category}
#' }
#'
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' # Compare two taxonomic assignments from the same phyloseq object
#' metrics <- tc_congruence_metrics(
#'   subset_taxa(Glom_otu, Phyla == "Fungi" | Phylum__eukaryome_Glomero == "Fungi"),
#'   ranks_1 = c("Phyla", "Class", "Order", "Family"),
#'   ranks_2 = c(
#'     "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
#'     "Order__eukaryome_Glomero", "Family__eukaryome_Glomero"
#'   )
#' )
#'
#' # View summary
#' metrics$summary
#'
#' # Get taxa with leaf-level incongruence
#' metrics$incongruent_leaves
tc_congruence_metrics <- function(
  physeq_1,
  physeq_2 = NULL,
  ranks_1,
  ranks_2 = NULL
) {
  verify_pq(physeq_1)

  if (is.null(physeq_2)) {
    physeq_2 <- physeq_1
  } else {
    verify_pq(physeq_2)
  }

  if (is.null(ranks_2)) {
    ranks_2 <- ranks_1
  }

  if (is.null(ranks_2) && is.null(physeq_2)) {
    stop("Your are comparing two equivalent object.
    Either ranks_2 or physeq_2 must be provided.")
  }

  common_taxa <- intersect(taxa_names(physeq_1), taxa_names(physeq_2))
  only_in_1_taxa <- setdiff(taxa_names(physeq_1), taxa_names(physeq_2))
  only_in_2_taxa <- setdiff(taxa_names(physeq_2), taxa_names(physeq_1))

  if (length(common_taxa) == 0) {
    warning("No common taxa names between physeq_1 and physeq_2")
    return(list(
      summary = data.frame(
        category = c("only_in_1", "only_in_2", "common"),
        count = c(length(only_in_1_taxa), length(only_in_2_taxa), 0),
        stringsAsFactors = FALSE
      ),
      only_in_1 = only_in_1_taxa,
      only_in_2 = only_in_2_taxa,
      total_congruent = character(0),
      partial_1_deeper = character(0),
      partial_2_deeper = character(0),
      incongruent_leaves = character(0),
      incongruent_nodes = character(0),
      details = data.frame()
    ))
  }

  tax_1 <- as.data.frame(physeq_1@tax_table[common_taxa, ranks_1, drop = FALSE])
  tax_2 <- as.data.frame(physeq_2@tax_table[common_taxa, ranks_2, drop = FALSE])

  clean_na <- function(x) {
    x[x == "" | x == "NA_NA" | x == "NA"] <- NA
    x
  }
  tax_1 <- as.data.frame(lapply(tax_1, clean_na), stringsAsFactors = FALSE)
  tax_2 <- as.data.frame(lapply(tax_2, clean_na), stringsAsFactors = FALSE)
  rownames(tax_1) <- common_taxa
  rownames(tax_2) <- common_taxa

  total_congruent <- character(0)
  partial_1_deeper <- character(0)
  partial_2_deeper <- character(0)
  incongruent_leaves <- character(0)
  incongruent_nodes <- character(0)
  unclassified_both <- character(0)
  classified_only_1 <- character(0)
  classified_only_2 <- character(0)

  details_list <- vector("list", length(common_taxa))

  for (i in seq_along(common_taxa)) {
    taxon <- common_taxa[i]
    path_1 <- as.character(tax_1[taxon, ])
    path_2 <- as.character(tax_2[taxon, ])

    depth_1 <- sum(!is.na(path_1))
    depth_2 <- sum(!is.na(path_2))

    leaf_pos_1 <- if (depth_1 > 0) {
      max(which(!is.na(path_1)))
    } else {
      NA
    }
    leaf_pos_2 <- if (depth_2 > 0) {
      max(which(!is.na(path_2)))
    } else {
      NA
    }
    last_valid_1 <- if (depth_1 > 0) {
      path_1[leaf_pos_1]
    } else {
      NA
    }
    last_valid_2 <- if (depth_2 > 0) {
      path_2[leaf_pos_2]
    } else {
      NA
    }

    leaf_match <- !is.na(last_valid_1) && !is.na(last_valid_2) &&
      last_valid_1 == last_valid_2

    details_list[[i]] <- data.frame(
      taxon = taxon,
      depth_1 = depth_1,
      depth_2 = depth_2,
      leaf_1 = last_valid_1,
      leaf_2 = last_valid_2,
      leaf_match = leaf_match,
      stringsAsFactors = FALSE
    )

    if (depth_1 == 0 && depth_2 == 0) {
      unclassified_both <- c(unclassified_both, taxon)
      details_list[[i]]$category <- "unclassified_both"
      next
    }

    if (depth_1 == 0) {
      classified_only_2 <- c(classified_only_2, taxon)
      details_list[[i]]$category <- "classified_only_2"
      next
    }

    if (depth_2 == 0) {
      classified_only_1 <- c(classified_only_1, taxon)
      details_list[[i]]$category <- "classified_only_1"
      next
    }

    max_leaf_pos <- max(leaf_pos_1, leaf_pos_2)
    min_leaf_pos <- min(leaf_pos_1, leaf_pos_2)

    paths_identical <- all(path_1 == path_2, na.rm = TRUE) &&
      all(is.na(path_1) == is.na(path_2))

    if (paths_identical) {
      total_congruent <- c(total_congruent, taxon)
      details_list[[i]]$category <- "congruent"
    } else if (leaf_pos_1 != leaf_pos_2) {
      shared_range <- seq_len(min_leaf_pos)
      shared_values_match <- all(
        path_1[shared_range] == path_2[shared_range],
        na.rm = TRUE
      )
      shared_na_match <- all(is.na(path_1[shared_range]) == is.na(path_2[shared_range]))

      if (shared_values_match && shared_na_match) {
        if (leaf_pos_1 > leaf_pos_2) {
          partial_1_deeper <- c(partial_1_deeper, taxon)
          details_list[[i]]$category <- "partial_1_deeper"
        } else {
          partial_2_deeper <- c(partial_2_deeper, taxon)
          details_list[[i]]$category <- "partial_2_deeper"
        }
      } else {
        first_diff <- which(
          path_1[shared_range] != path_2[shared_range] |
            xor(is.na(path_1[shared_range]), is.na(path_2[shared_range]))
        )[1]
        if (!is.na(first_diff) && first_diff == min_leaf_pos) {
          incongruent_leaves <- c(incongruent_leaves, taxon)
          details_list[[i]]$category <- "incongruent_leaves"
        } else {
          incongruent_nodes <- c(incongruent_nodes, taxon)
          details_list[[i]]$category <- "incongruent_nodes"
        }
      }
    } else {
      first_diff <- which(
        path_1 != path_2 | xor(is.na(path_1), is.na(path_2))
      )[1]
      if (!is.na(first_diff) && first_diff == leaf_pos_1) {
        incongruent_leaves <- c(incongruent_leaves, taxon)
        details_list[[i]]$category <- "incongruent_leaves"
      } else {
        incongruent_nodes <- c(incongruent_nodes, taxon)
        details_list[[i]]$category <- "incongruent_nodes"
      }
    }
  }

  details <- do.call(rbind, details_list)

  n_leaf_match <- sum(details$leaf_match, na.rm = TRUE)
  n_total <- length(taxa_names(physeq_1)) + length(only_in_2_taxa)
  summary_df <- data.frame(
    category = c(
      "only_in_1", "only_in_2", "classified_only_1", "classified_only_2",
      "unclassified_both", "total_congruent", "leaf_match_congruence",
      "partial_1_deeper", "partial_2_deeper",
      "incongruent_leaves", "incongruent_nodes"
    ),
    count = c(
      length(only_in_1_taxa), length(only_in_2_taxa), length(classified_only_1),
      length(classified_only_2), length(unclassified_both), length(total_congruent),
      n_leaf_match, length(partial_1_deeper), length(partial_2_deeper),
      length(incongruent_leaves), length(incongruent_nodes)
    ),
    stringsAsFactors = FALSE
  )
  summary_df$percentage <- round(100 * summary_df$count / n_total, 2)

  list(
    summary = summary_df,
    only_in_1 = only_in_1_taxa,
    only_in_2 = only_in_2_taxa,
    classified_only_1 = classified_only_1,
    classified_only_2 = classified_only_2,
    unclassified_both = unclassified_both,
    total_congruent = total_congruent,
    partial_1_deeper = partial_1_deeper,
    partial_2_deeper = partial_2_deeper,
    incongruent_leaves = incongruent_leaves,
    incongruent_nodes = incongruent_nodes,
    details = details
  )
}


################################################################################
#' Plot two taxonomy trees with linked correspondences
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Creates a visualization of two taxonomy trees facing each other, with
#' matching taxa aligned horizontally. Optionally draws lines connecting
#' matching leaf taxa. Useful for comparing taxonomic assignments from
#' different databases or classification methods.
#'
#' @param physeq_1 (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   for the left tree.
#' @param physeq_2 (phyloseq, default NULL) A \code{\link[phyloseq]{phyloseq-class}}
#'   object for the right tree. If NULL, uses physeq_1 (useful for comparing
#'   different rank columns from the same object).
#' @param ranks_1 (character vector, required) Taxonomic rank names to use for
#'   the first tree. Must match column names in physeq_1's tax_table.
#' @param ranks_2 (character vector, required) Taxonomic rank names to use for
#'   the second tree. Must match column names in physeq_2's tax_table.
#' @param tree_distance (numeric, default 1) Distance between the two trees.
#' @param show_tip_links (logical, default TRUE) Whether to draw lines between
#'   matching leaf taxa (tips). Internal nodes are never linked.
#' @param link_by_taxa (logical, default FALSE) If TRUE, links are drawn based
#'   on taxa (ASV/OTU) correspondence: for each taxon in physeq_1, a link is
#'   drawn to its corresponding tip in tree_2 based on its taxonomy. Line width
#'   is proportional to the number of taxa sharing each link. Requires that
#'   physeq_1 and physeq_2 have the same taxa names (rownames).
#' @param link_alpha (numeric, default 0.5) Transparency of the correspondence
#'   lines.
#' @param link_color (character, default "grey50") Color of the correspondence
#'   lines. Use NULL to color by label.
#' @param show_labels (logical, default TRUE) Whether to show node/tip labels.
#' @param label_size (numeric, default 2) Base size for labels. Labels are
#'   scaled by depth: shallower nodes (closer to root) get larger text (up to
#'   1.5x), deeper nodes get smaller text (down to 0.6x). All labels are
#'   positioned above their branch to avoid overlap.
#' @param internal_node_singletons (logical, default FALSE) Whether to create
#'   internal nodes for singletons. Passed to \code{\link{taxo2tree}}.
#' @param use_taxa_names (logical, default FALSE) Whether to use taxa names
#'   (e.g., ASV_1, ASV_2) as terminal leaves. If FALSE (default), collapses
#'   identical taxonomy paths and uses the lowest rank value as tip labels.
#'   Passed to \code{\link{taxo2tree}}.
#'
#' @return A ggplot2 object that can be further customized.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' # Compare two taxonomic assignments from the same phyloseq object
#' # using different rank columns (physeq_2 defaults to physeq_1)
#' tc_linked_trees(
#'   Glom_otu,
#'   ranks_1 = c("Kingdom", "Class", "Order", "Family"),
#'   ranks_2 = c(
#'     "Kingdom__eukaryome_Glomero", "Class__eukaryome_Glomero",
#'     "Order__eukaryome_Glomero", "Family__eukaryome_Glomero"
#'   )
#' )
#'
#' \dontrun{
#' # Link by taxa with line width proportional to ASV count
#' tc_linked_trees(
#'   Glom_otu,
#'   ranks_1 = c("Kingdom", "Class", "Order", "Family"),
#'   ranks_2 = c(
#'     "Kingdom__eukaryome_Glomero", "Class__eukaryome_Glomero",
#'     "Order__eukaryome_Glomero", "Family__eukaryome_Glomero"
#'   ),
#'   link_by_taxa = TRUE
#' )
#'
#' # Without tip links
#' tc_linked_trees(
#'   Glom_otu,
#'   ranks_1 = c("Kingdom", "Class", "Order", "Family"),
#'   ranks_2 = c(
#'     "Kingdom__eukaryome_Glomero", "Class__eukaryome_Glomero",
#'     "Order__eukaryome_Glomero", "Family__eukaryome_Glomero"
#'   ),
#'   show_tip_links = FALSE
#' )
#'
#' # Color links by label
#' tc_linked_trees(
#'   Glom_otu,
#'   ranks_1 = c("Kingdom", "Class", "Order", "Family"),
#'   ranks_2 = c(
#'     "Kingdom__eukaryome_Glomero", "Class__eukaryome_Glomero",
#'     "Order__eukaryome_Glomero", "Family__eukaryome_Glomero"
#'   ),
#'   link_color = NULL
#' ) +
#'   theme(legend.position = "none")
#'
#' # Using one phyloseq object for both trees
#' # and different rank levels with link_by_taxa
#' tc_linked_trees(
#'   Glom_otu,
#'   ranks_1 = c("Kingdom", "Class", "Order", "Family"),
#'   ranks_2 = c(
#'     "Kingdom__eukaryome_Glomero", "Class__eukaryome_Glomero",
#'     "Order__eukaryome_Glomero", "Family__eukaryome_Glomero"
#'   ), link_color = NULL,
#'   link_by_taxa = TRUE,
#' ) +
#'   theme(legend.position = "none")
#' }
tc_linked_trees <- function(
  physeq_1,
  physeq_2 = NULL,
  ranks_1,
  ranks_2,
  tree_distance = 1,
  show_tip_links = TRUE,
  link_by_taxa = FALSE,
  link_alpha = 0.5,
  link_color = "grey50",
  show_labels = TRUE,
  label_size = 2,
  internal_node_singletons = FALSE,
  use_taxa_names = FALSE
) {
  rlang::check_installed("ggtree", reason = "to plot taxonomy trees")

  verify_pq(physeq_1)

  if (is.null(physeq_2)) {
    physeq_2 <- physeq_1
  } else {
    verify_pq(physeq_2)
  }

  tree_1 <- taxo2tree(
    physeq_1,
    ranks = ranks_1,
    internal_node_singletons = internal_node_singletons,
    use_taxa_names = use_taxa_names
  )
  tree_2 <- taxo2tree(
    physeq_2,
    ranks = ranks_2,
    internal_node_singletons = internal_node_singletons,
    use_taxa_names = use_taxa_names
  )

  clean_label <- function(x) {
    x <- gsub("^NA$", "", x, fixed = TRUE)
    x <- gsub("NA_NA", "", x, fixed = TRUE)
    x[is.na(x)] <- ""
    x
  }

  tips_1 <- tree_1$tip.label
  tips_2 <- tree_2$tip.label
  common_tips <- intersect(tips_1, tips_2)

  if (length(common_tips) > 1) {
    constraint_order <- tips_1[tips_1 %in% common_tips]
    tree_2 <- ape::rotateConstr(tree_2, constraint_order)
  }

  p1 <- ggtree::ggtree(tree_1, layout = "rectangular")
  d1 <- p1$data
  d1$label <- clean_label(d1$label)

  p2 <- ggtree::ggtree(tree_2, layout = "rectangular")
  d2 <- p2$data
  d2$label <- clean_label(d2$label)

  n_tips_1 <- sum(d1$isTip)
  n_tips_2 <- sum(d2$isTip)
  max_tips <- max(n_tips_1, n_tips_2)

  scale_y <- function(data, max_tips) {
    tips <- data[data$isTip, ]
    if (nrow(tips) == 0) {
      return(data)
    }

    tip_y_order <- order(tips$y)
    new_tip_y <- seq_len(nrow(tips))
    tips$y[tip_y_order] <- new_tip_y

    data[data$isTip, "y"] <- tips$y

    nodes <- data[!data$isTip, ]
    for (i in seq_len(nrow(nodes))) {
      node_id <- nodes$node[i]
      children_idx <- which(data$parent == node_id)
      if (length(children_idx) > 0) {
        nodes$y[i] <- mean(data$y[children_idx])
      }
    }
    data[!data$isTip, "y"] <- nodes$y

    data
  }

  d1 <- scale_y(d1, max_tips)
  d2 <- scale_y(d2, max_tips)

  d1$tree <- "tree_1"
  d2$tree <- "tree_2"

  max_x1 <- max(d1$x, na.rm = TRUE)
  max_x2 <- max(d2$x, na.rm = TRUE)
  d2$x <- max_x2 - d2$x + max_x1 + tree_distance

  labels_1 <- d1$label[d1$label != ""]
  labels_2 <- d2$label[d2$label != ""]
  common_labels <- intersect(labels_1, labels_2)

  p <- ggplot() +
    ggtree::geom_tree(data = d1, aes(x = x, y = y), layout = "rectangular") +
    ggtree::geom_tree(data = d2, aes(x = x, y = y), layout = "rectangular")

  if (show_tip_links) {
    if (link_by_taxa) {
      common_taxa <- intersect(taxa_names(physeq_1), taxa_names(physeq_2))
      if (length(common_taxa) == 0) {
        warning("No common taxa names between physeq_1 and physeq_2")
      } else {
        get_tip_label <- function(tax_row, ranks) {
          vals <- tax_row[ranks]
          vals <- vals[!is.na(vals) & vals != "" & vals != "NA_NA"]
          if (length(vals) > 0) vals[length(vals)] else "Unknown"
        }

        tax_1 <- as.data.frame(physeq_1@tax_table[common_taxa, , drop = FALSE])
        tax_2 <- as.data.frame(physeq_2@tax_table[common_taxa, , drop = FALSE])

        taxa_links <- data.frame(
          taxon = common_taxa,
          tip_1 = apply(tax_1, 1, get_tip_label, ranks = ranks_1),
          tip_2 = apply(tax_2, 1, get_tip_label, ranks = ranks_2),
          stringsAsFactors = FALSE
        )

        link_counts <- taxa_links |>
          dplyr::count(tip_1, tip_2, name = "n_taxa")

        tips_d1 <- d1[d1$isTip, c("label", "x", "y")]
        tips_d2 <- d2[d2$isTip, c("label", "x", "y")]

        link_data <- link_counts |>
          dplyr::inner_join(tips_d1, by = c("tip_1" = "label")) |>
          dplyr::rename(x_1 = x, y_1 = y) |>
          dplyr::inner_join(tips_d2, by = c("tip_2" = "label")) |>
          dplyr::rename(x_2 = x, y_2 = y)

        if (nrow(link_data) > 0) {
          max_n <- max(link_data$n_taxa)
          link_data$line_width <- 0.5 + (link_data$n_taxa / max_n) * 3

          if (is.null(link_color)) {
            p <- p +
              geom_segment(
                data = link_data,
                aes(
                  x = x_1,
                  y = y_1,
                  xend = x_2,
                  yend = y_2,
                  linewidth = line_width,
                  colour = tip_1
                ),
                alpha = link_alpha
              ) +
              scale_linewidth_identity()
          } else {
            p <- p +
              geom_segment(
                data = link_data,
                aes(
                  x = x_1,
                  y = y_1,
                  xend = x_2,
                  yend = y_2,
                  linewidth = line_width
                ),
                alpha = link_alpha,
                colour = link_color
              ) +
              scale_linewidth_identity()
          }
        }
      }
    } else if (length(common_labels) > 0) {
      tips_d1 <- d1[
        d1$isTip & d1$label %in% common_labels,
        c("label", "x", "y")
      ]
      tips_d2 <- d2[
        d2$isTip & d2$label %in% common_labels,
        c("label", "x", "y")
      ]

      if (nrow(tips_d1) > 0 && nrow(tips_d2) > 0) {
        link_data <- merge(
          tips_d1,
          tips_d2,
          by = "label",
          suffixes = c("_1", "_2")
        )

        if (nrow(link_data) > 0) {
          if (is.null(link_color)) {
            p <- p +
              geom_segment(
                data = link_data,
                aes(
                  x = x_1,
                  y = y_1,
                  xend = x_2,
                  yend = y_2,
                  colour = label
                ),
                alpha = link_alpha
              )
          } else {
            p <- p +
              geom_segment(
                data = link_data,
                aes(x = x_1, y = y_1, xend = x_2, yend = y_2),
                alpha = link_alpha,
                colour = link_color
              )
          }
        }
      }
    }
  }

  if (show_labels) {
    calc_depth <- function(data) {
      max_x <- max(data$x, na.rm = TRUE)
      min_x <- min(data$x, na.rm = TRUE)
      x_range <- max_x - min_x
      if (x_range == 0) {
        x_range <- 1
      }
      data$depth <- (data$x - min_x) / x_range
      data
    }

    d1 <- calc_depth(d1)
    d2 <- calc_depth(d2)
    d2$depth <- 1 - d2$depth

    d1_labels <- d1[d1$label != "", ]
    d2_labels <- d2[d2$label != "", ]

    min_size <- label_size * 0.6
    max_size <- label_size * 1.5
    d1_labels$text_size <- min_size +
      (1 - d1_labels$depth) * (max_size - min_size)
    d2_labels$text_size <- min_size +
      (1 - d2_labels$depth) * (max_size - min_size)

    d1_tips <- d1_labels[d1_labels$isTip, ]
    d1_nodes <- d1_labels[!d1_labels$isTip, ]
    d2_tips <- d2_labels[d2_labels$isTip, ]
    d2_nodes <- d2_labels[!d2_labels$isTip, ]

    if (nrow(d1_tips) > 0) {
      p <- p +
        geom_text(
          data = d1_tips,
          aes(x = x, y = y, label = label, size = text_size),
          hjust = 1,
          vjust = -0.3,
          nudge_x = -0.05
        )
    }

    if (nrow(d2_tips) > 0) {
      p <- p +
        geom_text(
          data = d2_tips,
          aes(x = x, y = y, label = label, size = text_size),
          hjust = 0,
          vjust = -0.3,
          nudge_x = 0.05
        )
    }

    if (nrow(d1_nodes) > 0) {
      p <- p +
        geom_text(
          data = d1_nodes,
          aes(x = x, y = y, label = label, size = text_size),
          hjust = 1,
          vjust = -0.5,
          nudge_x = -0.05
        )
    }

    if (nrow(d2_nodes) > 0) {
      p <- p +
        geom_text(
          data = d2_nodes,
          aes(x = x, y = y, label = label, size = text_size),
          hjust = 0,
          vjust = -0.5,
          nudge_x = 0.05
        )
    }

    p <- p + scale_size_identity()

    # add a title with the name of the ranks and phyloseq objects
    if (!is.null(physeq_2) &&
      !identical(physeq_1, physeq_2)) {
      name_1 <- deparse(substitute(physeq_1))
      name_2 <- deparse(substitute(physeq_2))
    } else {
      name_1 <- deparse(substitute(physeq_1))
      name_2 <- deparse(substitute(physeq_1))
    }
    p <- p + labs(
      title =
        paste0(name_1, " (left) vs ", name_2, " (right)"),
      subtitle =
        paste0(
          paste(ranks_1, collapse = " > "), "\n",
          paste(ranks_2, collapse = " > ")
        )
    )
  }

  p + theme_void()
}
