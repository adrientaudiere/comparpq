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
#' @param physeq_2 (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   for the right tree.
#' @param ranks_1 (character vector, required) Taxonomic rank names to use for
#'   the first tree. Must match column names in physeq_1's tax_table.
#' @param ranks_2 (character vector, required) Taxonomic rank names to use for
#'   the second tree. Must match column names in physeq_2's tax_table.
#' @param tree_distance (numeric, default 1) Distance between the two trees.
#' @param show_tip_links (logical, default TRUE) Whether to draw lines between
#'   matching leaf taxa (tips). Internal nodes are never linked.
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
#' # using different rank columns
#' Glom_otu_ab <- subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 10000)
#'
#' tc_linked_trees(
#'   Glom_otu, Glom_otu_ab,
#'   ranks_1 = c("Kingdom", "Class", "Order", "Family"),
#'   ranks_2 = c(
#'     "Kingdom__eukaryome_Glomero", "Class__eukaryome_Glomero",
#'     "Order__eukaryome_Glomero", "Family__eukaryome_Glomero"
#'   )
#' )
#'
#' \dontrun{
#' # Without tip links
#' tc_linked_trees(
#'   Glom_otu, Glom_otu_ab,
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
#'   Glom_otu, Glom_otu_ab,
#'   ranks_1 = c("Kingdom", "Class", "Order", "Family"),
#'   ranks_2 = c(
#'     "Kingdom__eukaryome_Glomero", "Class__eukaryome_Glomero",
#'     "Order__eukaryome_Glomero", "Family__eukaryome_Glomero"
#'   ),
#'   link_color = NULL
#' ) +
#'   theme(legend.position = "none")
#' }
tc_linked_trees <- function(physeq_1,
                            physeq_2,
                            ranks_1,
                            ranks_2,
                            tree_distance = 1,
                            show_tip_links = TRUE,
                            link_alpha = 0.5,
                            link_color = "grey50",
                            show_labels = TRUE,
                            label_size = 2,
                            internal_node_singletons = FALSE,
                            use_taxa_names = FALSE) {
  rlang::check_installed("ggtree", reason = "to plot taxonomy trees")

  verify_pq(physeq_1)
  verify_pq(physeq_2)

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
    x <- gsub("^NA$", "", x)
    x <- gsub("NA_NA", "", x)
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
    if (nrow(tips) == 0) return(data)

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

  if (show_tip_links && length(common_labels) > 0) {
    tips_d1 <- d1[d1$isTip & d1$label %in% common_labels, c("label", "x", "y")]
    tips_d2 <- d2[d2$isTip & d2$label %in% common_labels, c("label", "x", "y")]

    if (nrow(tips_d1) > 0 && nrow(tips_d2) > 0) {
      link_data <- merge(tips_d1, tips_d2, by = "label", suffixes = c("_1", "_2"))

      if (nrow(link_data) > 0) {
        if (is.null(link_color)) {
          p <- p +
            geom_segment(
              data = link_data,
              aes(
                x = x_1, y = y_1,
                xend = x_2, yend = y_2,
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

  if (show_labels) {
    calc_depth <- function(data) {
      max_x <- max(data$x, na.rm = TRUE)
      min_x <- min(data$x, na.rm = TRUE)
      x_range <- max_x - min_x
      if (x_range == 0) x_range <- 1
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
    d1_labels$text_size <- min_size + (1 - d1_labels$depth) * (max_size - min_size)
    d2_labels$text_size <- min_size + (1 - d2_labels$depth) * (max_size - min_size)

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
  }

  p + theme_void()
}
