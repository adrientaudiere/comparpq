#' Convert taxonomy dataframe to phylogenetic tree
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"> <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description Creates a phylo object from a taxonomy
#' table with hierarchical taxonomic ranks as columns.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#' @param ranks Character vector specifying the column names to use as
#'  taxonomic ranks, ordered from highest to lowest.
#'  By default: c("Domain", "Phylum", "Class", "Order", "Family",
#'   "Genus", "Species").
#' @param internal_node_singletons Logical, if TRUE, create internal nodes
#'   for singleton. If FALSE, internal nodes with only one descendant
#'   are discarded.
#' @param use_taxa_names Logical, if TRUE (default), use the taxa names
#'   (rownames, e.g., ASV_1, ASV_2) as terminal leaves. If FALSE, collapse
#'   identical taxonomy paths and use the lowest rank value as tip labels.
#'   This is useful for creating cleaner trees that show only taxonomic
#'   structure without individual ASV/OTU names.
#' @return A phylo object (ape package) representing the taxonomic tree.
#'
#' @author Adrien Taudière
#'
#' @examples
#' data_fungi_mini@phy_tree <-
#'   phyloseq::phy_tree(taxo2tree(data_fungi_mini,
#'     ranks = c(
#'       "Domain", "Phylum", "Class", "Order",
#'       "Family", "Genus", "Genus_species"
#'     )
#'   ))
#'
#' library(ggtree)
#'
#' ggtree(data_fungi_mini@phy_tree) +
#'   geom_nodelab(size = 2, nudge_x = -0.2, nudge_y = 0.6) +
#'   geom_tiplab()
#'
#' psm <- psmelt(data_fungi_mini) |>
#'   group_by(OTU) |>
#'   summarize(
#'     Mol_Abund = sum(Abundance),
#'     Frequency = sum(Abundance > 0),
#'     Class = tidyr::replace_na(unique(Class), "Unknown")
#'   )
#'
#'
#' ggtree(data_fungi_mini@phy_tree) %<+% psm +
#'   geom_tiplab(offset = .6, hjust = .5) +
#'   geom_tippoint(aes(
#'     shape = Class,
#'     color = log10(1 + Mol_Abund),
#'     size = Frequency
#'   )) +
#'   geom_nodelab(nudge_x = -0.2, nudge_y = 0.6) +
#'   scale_color_viridis_c()
#'
#'
#' # Without internal node singletons and only until the Genus level
#' tree_wo_singletons <-
#'   phyloseq::phy_tree(taxo2tree(data_fungi_mini,
#'     ranks = c(
#'       "Domain", "Phylum", "Class", "Order",
#'       "Family", "Genus"
#'     ),
#'     internal_node_singletons = FALSE
#'   ))
#'
#' ggtree::ggtree(tree_wo_singletons) +
#'   ggtree::geom_nodelab(size = 2, nudge_x = -0.2, nudge_y = 0.6) +
#'   ggtree::geom_tiplab()
#'
#' # Without taxa names (collapse identical paths)
#' tree_no_taxa <- taxo2tree(data_fungi_mini,
#'   ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
#'   use_taxa_names = FALSE
#' )
#'
#' ggtree::ggtree(tree_no_taxa) +
#'   ggtree::geom_nodelab(size = 2, nudge_x = -0.2, nudge_y = 0.6) +
#'   ggtree::geom_tiplab()
#' @export
taxo2tree <- function(
  physeq,
  ranks = c(
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  ),
  internal_node_singletons = TRUE,
  use_taxa_names = TRUE
) {
  verify_pq(physeq)

  missing_ranks <- setdiff(ranks, colnames(physeq@tax_table))
  if (length(missing_ranks) > 0) {
    stop(
      "The following ranks are not present in the tax_table: ",
      paste(missing_ranks, collapse = ", ")
    )
  }

  tax_df <- physeq@tax_table |>
    as.data.frame() |>
    select(all_of(ranks)) |>
    mutate(across(everything(), ~ na_if(.x, "NA_NA")))

  if (!use_taxa_names) {
    tax_df <- tax_df |>
      dplyr::distinct() |>
      dplyr::mutate(
        .tip_label = apply(
          dplyr::pick(dplyr::everything()),
          1,
          \(x) {
            non_na <- x[!is.na(x) & x != ""]
            if (length(non_na) > 0) non_na[length(non_na)] else "Unknown"
          }
        )
      )

    dup_labels <- tax_df$.tip_label[duplicated(tax_df$.tip_label)]
    if (length(dup_labels) > 0) {
      for (lbl in unique(dup_labels)) {
        idx <- which(tax_df$.tip_label == lbl)
        tax_df$.tip_label[idx] <- paste0(lbl, "_", seq_along(idx))
      }
    }

    rownames(tax_df) <- tax_df$.tip_label
    tax_df$.tip_label <- NULL
  }

  tax_mat <- as.matrix(tax_df[, ranks, drop = FALSE])

  newick <- build_newick(
    tax_mat = tax_mat,
    ranks = ranks,
    internal_node_singletons = internal_node_singletons
  )

  tree <- ape::read.tree(text = newick)

  return(tree)
}


#' Build Newick format string from taxonomy matrix
#'
#' @description
#' Internally used by \code{taxo2tree()} to build the Newick format string.
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">  <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @inheritParams taxo2tree
#' @param internal_node_singletons Logical, if TRUE, create internal nodes
#'   for singleton. If FALSE, internal nodes with only one descendant
#'   are discarded.
#' @return Newick format string
#' @noRd
#' @keywords internal
build_newick <- function(tax_mat, ranks, internal_node_singletons = TRUE) {
  n_ranks <- length(ranks)
  tip_labels <- rownames(tax_mat)

  clades <- lapply(1:nrow(tax_mat), function(i) {
    list(rows = i, newick = tip_labels[i], label = "")
  })

  for (rank_idx in n_ranks:1) {
    current_rank <- ranks[rank_idx]
    new_clades <- list()

    remaining_clades <- clades

    while (length(remaining_clades) > 0) {
      # Take first clade
      first_clade <- remaining_clades[[1]]
      first_rows <- first_clade$rows
      first_val <- tax_mat[first_rows[1], rank_idx]

      if (is.na(first_val)) {
        same_val_idx <- 1
      } else {
        same_val_idx <- which(sapply(remaining_clades, function(clade) {
          val <- tax_mat[clade$rows[1], rank_idx]
          !is.na(val) && val == first_val
        }))
      }

      # Extract clades with same value
      grouped_clades <- remaining_clades[same_val_idx]
      remaining_clades <- remaining_clades[-same_val_idx]

      if (length(grouped_clades) == 1) {
        merged_clade <- grouped_clades[[1]]
        if (internal_node_singletons && !is.na(first_val)) {
          merged_clade$newick <- paste0(
            "(",
            merged_clade$newick,
            ")",
            first_val
          )
          merged_clade$label <- first_val
        } else if (!is.na(first_val)) {
          merged_clade$label <- first_val
        }
      } else {
        all_rows <- unlist(lapply(grouped_clades, function(x) x$rows))
        all_newick <- sapply(grouped_clades, function(x) x$newick)

        combined_newick <- paste0("(", paste(all_newick, collapse = ","), ")")

        node_label <- ""
        if (rank_idx > 1 && !is.na(first_val)) {
          node_label <- first_val
          combined_newick <- paste0(combined_newick, node_label)
        }

        merged_clade <- list(
          rows = all_rows,
          newick = combined_newick,
          label = node_label
        )
      }

      new_clades[[length(new_clades) + 1]] <- merged_clade
    }

    clades <- new_clades
  }

  if (length(clades) == 1) {
    newick <- paste0(clades[[1]]$newick, ";")
  } else {
    all_newick <- sapply(clades, function(x) x$newick)
    newick <- paste0("(", paste(all_newick, collapse = ","), ");")
  }

  return(newick)
}
