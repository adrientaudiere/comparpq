# cross_tax_pq(data_fungi_mini)
cross_tax_pq <- function(physeq, rank_1 = "Family", rank_2 = "Class", ...) {
  taxatab <- dplyr::select(as.data.frame(unclass(physeq@tax_table)), one_of(c(rank_1, rank_2)))
  tbl_sum <- gtsummary::tbl_cross(taxatab, ...)
  return(tbl_sum)
}





## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## IN PROGRESS #todo
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Sankey for comparing two ranks
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

tc_sankey <- function(physeq,
                      rank_1 = NULL,
                      rank_2 = NULL) {
  require(ggalluvial)
  df_sank <- as.data.frame(unclass(physeq@tax_table)) |>
    dplyr::select(rank_1, rank_2) |>
    dplyr::count(.data[[rank_1]], .data[[rank_2]])

  ggplot(df_sank, aes(
    axis1 = .data[[rank_1]],
    axis2 = .data[[rank_2]],
    y = n
  )) +
    geom_alluvium() +
    geom_stratum(width = 1 / 12, color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(
      limits = c("Gender", "Dept"),
      expand = c(.05, .05)
    ) +
    theme_void()
}

data("Glom_otu")
tc_sankey(Glom_otu, "Family__eukaryome_Glomero", "Family")
















## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Taxonomy tree with linked correspondence
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# todo : add a width of branch in function of the number of ASV/sequences

ranks_1 <- c(1:6)
ranks_2 <- c(9:14)
collapse_taxa <- TRUE

data("Glom_otu")

physeq_interm <- subset_samples(Glom_otu, sample_sums(Glom_otu) > 50000)
physeq_interm <- subset_taxa_pq(physeq_interm, taxa_sums(physeq_interm) >
  10000)
physeq_interm <- clean_pq(subset_taxa_pq(physeq_interm, taxa_sums(as_binary_otu_table(physeq_interm)) >
  30))
physeq <- physeq_interm




require(ggtree)

formula_taxo_1 <- formula(paste0("~", paste0(colnames(physeq@tax_table)[ranks_1],
  collapse =
    "/"
)))

data_taxo_1 <- as.data.frame(physeq@tax_table[, ranks_1]) |>
  dplyr::mutate_if(is.character, as.factor)

if (collapse_taxa) {
  data_taxo_1 <- data_taxo_1 |>
    group_by(across(everything())) |>
    summarise(n = n())
}

phy_tree_1 <- ape::as.phylo.formula(formula_taxo_1,
  data = data_taxo_1, collapse =
    T
)
p1 <- ggtree(phy_tree_1, layout = "roundrect")


formula_taxo_2 <- formula(paste0("~", paste0(colnames(physeq@tax_table)[ranks_2],
  collapse =
    "/"
)))

data_taxo_2 <- as.data.frame(physeq@tax_table[, ranks_2]) |>
  dplyr::mutate_if(is.character, as.factor)

if (collapse_taxa) {
  data_taxo_2 <- data_taxo_2 |>
    group_by(across(everything())) |>
    summarise(n = n())
}

phy_tree_2 <- ape::as.phylo.formula(formula_taxo_2,
  data = data_taxo_2, collapse =
    T
)
p2 <- ggtree(phy_tree_2, layout = "roundrect")

# https://yulab-smu.top/treedata-book/chapter2.html#ggtree-fortify
d1 <- p1$data
d1$label <- gsub("NA", "", d1$label)

d2 <- p2$data
d2$label <- gsub("NA", "", d2$label)

## reverse x-axis and
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1 + geom_tree(data = d2, layout = "ellipse") +
  ggnewscale::new_scale_fill()

dd <- bind_rows(d1, d2) %>%
  dplyr::filter(label %in% names(table(label))[table(label) > 1]) %>%
  dplyr::filter(!is.na(label)) %>%
  dplyr::filter(label != "")

pp +
  ggrepel::geom_label_repel(aes(label = label), size = 1, data = d2) +
  ggrepel::geom_label_repel(aes(label = label), size = 1, data = d1) +
  geom_line(aes(x, y, group = label, colour = label),
    alpha = 0.5,
    data = dd
  )
