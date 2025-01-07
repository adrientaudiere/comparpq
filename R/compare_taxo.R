################################################################################
################################################################################
#' Barchart of ratio to compare 2 taxonomic ranks
#'
#'
#' Useful to compare taxonomy from two different source/db/algo side-by-side
#'
#' @inheritParams tc_points_matrix
#' @param color_rank Define the taxonomic rank for color
#'   as the number or the name of the column in tax_table slot
#' @param point_size
#' @param point_alpha
#'
#' @return A ggplot2 object
#' @export
#' @author Adrien Taudière
#' @examples
#'
#' tc_bar(Glom_otu,
#'             rank_1 = 5,
#'             rank_2 = 13,
#'             color_rank = 3)

tc_bar <- function(physeq,
                         rank_1,
                         rank_2,
                         color_rank,
                         point_size = 0.3,
                         point_alpha = 0.3,
                         merge_sample_by = NULL) {
  verify_pq(physeq)

  if (!is.null(merge_sample_by)) {
    physeq <- clean_pq(merge_samples2(physeq, merge_sample_by))
  }

  psm <- psmelt(physeq)

  ranks1 <- colnames(physeq@tax_table[, c(rank_1, color_rank)])
  ranks2 <- colnames(physeq@tax_table[, c(rank_2, color_rank)])

  psm2 <- psm |>
    group_by(across(all_of(c(
      "Sample", ranks1, ranks2
    )))) |>
    summarise(Abundance = sum(Abundance), n = n()) |>
    dplyr::filter(Abundance > 0)

  p <- ggplot(psm2, aes(
    x = .data[[ranks2[1]]],
    y = log10(Abundance),
    fill = .data[[ranks2[[2]]]],
    label = .data[[ranks2[1]]]
  )) +
    stat_summary(fun.y = mean,
                 position = "dodge",
                 geom = "bar") +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = .2,
      position = position_dodge(.9)
    ) +
    geom_jitter(size = point_size, alpha = point_alpha) +
    stat_summary(
      aes(
        x = .data[[ranks1[1]]],
        y = -log10(Abundance),
        fill = .data[[ranks1[[2]]]]
      ),
      fun.y = mean,
      position = "dodge",
      geom = "bar"
    ) +
    stat_summary(
      aes(x = .data[[ranks1[1]]], y = -log10(Abundance)),
      fun.data = mean_se,
      geom = "errorbar",
      width = .2,
      position = position_dodge(.9)
    ) +
    theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )) +
    geom_hline(aes(yintercept = 0)) +
    geom_jitter(aes(x = .data[[ranks1[1]]], y = -log10(Abundance)), size =
                  point_size, alpha = point_alpha)

  return(p)
}
################################################################################
################################################################################


################################################################################
################################################################################
#' Matrix of point to compare two taxonomic ranks
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @param physeq (required): a \code{\link[phyloseq]{phyloseq-class}} object obtained
#'   using the `phyloseq` package.
#' @param rank_1 Define the first taxonomic rank as the number or the name
#'    of the column in tax_table slot
#' @param rank_2 Define the second taxonomic rank as the number or the name
#'    of the column in tax_table slot
#' @param color_1 Color for rank_1 values
#' @param color_2 Color for rank_1 values
#' @param stat_across_sample Either "mean" or "sum". How the abundance is
#'   computed across samples.
#' @param merge_sample_by a vector to determine
#'   which samples to merge using [merge_samples2()] function.
#'   Need to be in \code{physeq@sam_data}
#'
#' @return A ggplot2 object
#' @export
#' @author Adrien Taudière
#' @examples
#' tc_points_matrix(Glom_otu, 6, 14)
#' tc_points_matrix(Glom_otu, 4, 12)
#' tc_points_matrix(Glom_otu, 4, 12, stat_across_sample = "mean")
#'
#' Glom_otu@sam_data$unique_value <- rep("samp", nsamples(Glom_otu))
#' tc_points_matrix(as_binary_otu_table(Glom_otu), 5, 13,
#'                        stat_across_sample = "sum", merge_sample_by="unique_value")
#' tc_points_matrix(as_binary_otu_table(Glom_otu), 5, 13,
#'                        stat_across_sample = "mean")
#' tc_points_matrix(Glom_otu, 5, 13,
#'                        stat_across_sample = "mean")

tc_points_matrix <- function(physeq,
                                   rank_1,
                                   rank_2,
                                   color_1 = "#dc863b",
                                   color_2 = "#2e7891",
                                   stat_across_sample = "sum",
                                   merge_sample_by = NULL) {
  verify_pq(physeq)

  if (!is.null(merge_sample_by)) {
    physeq <- clean_pq(merge_samples2(physeq, merge_sample_by))
  }

  psm <- psmelt(physeq)

  ranks1 <- colnames(physeq@tax_table[, c(rank_1)])
  ranks2 <- colnames(physeq@tax_table[, c(rank_2)])

  psm2 <- psm |>
    group_by(across(all_of(c(
      "Sample", ranks1, ranks2
    )))) |>
    summarise(Abundance = sum(Abundance), n = n()) |>
    dplyr::filter(Abundance > 0)

  psm2_summary_taxo1 <- psm2 %>%
    group_by(.data[[ranks1[[1]]]]) %>%
    summarise(mean_ab = mean(Abundance),
              sum_ab = sum(Abundance)) |>
    mutate(Rank = .data[[ranks1[[1]]]])

  psm2_summary_taxo2 <- psm2 %>%
    group_by(.data[[ranks2[[1]]]]) %>%
    summarise(mean_ab = mean(Abundance),
              sum_ab = sum(Abundance)) |>
    mutate(Rank = .data[[ranks2[[1]]]])

  psm2_summary <- full_join(psm2_summary_taxo1,
                            psm2_summary_taxo2,
                            by = join_by(Rank == Rank))

  if (stat_across_sample == "mean") {
    p <- ggplot(psm2_summary, aes(x = .data[[ranks1[[1]]]],
                                  y = .data[[ranks2[[1]]]],
                                  size =  mean_ab.x)) +
      geom_point(shape = 21,
                 color = color_1,
                 fill = alpha(color_1, 0.2)) +
      geom_point(
        aes(size = mean_ab.y),
        color = color_2,
        shape = 21,
        ,
        fill = alpha(color_2, 0.2)
      )

  } else if (stat_across_sample == "sum") {
    p <- ggplot(psm2_summary, aes(x = .data[[ranks1[[1]]]],
                                  y = .data[[ranks2[[1]]]],
                                  size =  sum_ab.x)) +
      geom_point(shape = 21,
                 color = color_1,
                 fill = alpha(color_1, 0.2)) +
      geom_point(
        aes(size = sum_ab.y),
        color = color_2,
        shape = 21,
        ,
        fill = alpha(color_2, 0.2)
      )

  } else {
    stop("Param stat_across_sample must be set to mean or sum !")
  }
  p <- p+theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))

  return(p)
}
################################################################################
################################################################################

################################################################################
################################################################################
#' Rainplot of the nb taxa assigned (not NA)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @inheritParams tc_points_matrix
#' @param ranks The ranks to include in the rainplot.
#' @param min_nb_seq Minimum number of sequences to filter out taxa with low
#'   abundance
#' @param merge_sample_by a vector to determine
#'   which samples to merge using [merge_samples2()] function.
#'   Need to be in \code{physeq@sam_data}
#' @param sample_colored (logical, default FALSE) Do points are colored by
#'   samples?
#' @param sample_linked (logical, default FALSE) Do points are linked by
#'   samples?
#' @param ...
#'
#' @return A ggplot2 object
#' @export
#' @author Adrien Taudière
#' @examples
#' Glom_otu@sam_data$tmt_type <- paste0(Glom_otu@sam_data$Tmt, "_", Glom_otu@sam_data$Type)
#' rainplot_taxo_na(Glom_otu)
#' rainplot_taxo_na(
#'   Glom_otu,
#'   merge_sample_by = "tmt_type",
#'   sample_colored = TRUE,
#'   sample_linked = TRUE
#' )
#' rainplot_taxo_na(Glom_otu, ranks = c(4, 12), rain.side = 'f1x1')
#' rainplot_taxo_na(
#'   Glom_otu,
#'   ranks = c(6, 14),
#'   rain.side = 'f1x1',
#'   sample_linked = TRUE
#' ) +
#'   theme(legend.position = "none")

rainplot_taxo_na <- function(physeq,
                             ranks = c(1:7),
                             min_nb_seq = 0,
                             merge_sample_by = NULL,
                             sample_colored = FALSE,
                             sample_linked = FALSE,
                             ...) {
  verify_pq(physeq)

  if (!is.null(merge_sample_by)) {
    physeq <- clean_pq(merge_samples2(physeq, merge_sample_by))
  }


  rank_names <- colnames(physeq@tax_table)[ranks]

  psm <- psmelt(physeq)

  sum_not_na <- function(x) {
    sum(!is.na(x))
  }

  prop_not_na <- function(x) {
    sum(!is.na(x)) / length(x)
  }

  psm2 <- psm |>
    dplyr::filter(Abundance > min_nb_seq) |>
    group_by(Sample) |>
    summarise(across(all_of(rank_names), prop_not_na)) |>
    tidyr::pivot_longer(!Sample)

  psm2$name <- factor(psm2$name, levels = rank_names)

  p <- ggplot(psm2, aes(x = name, y = value, fill = name)) +
    coord_flip()

  if (sample_colored && sample_linked) {
    p <- p + ggrain::geom_rain(cov = "Sample", id.long.var = "Sample", ...)
  } else if (sample_colored && !sample_linked) {
    p <- p + ggrain::geom_rain(cov = "Sample", ...)
  } else if (!sample_colored && sample_linked) {
    p <- p + ggrain::geom_rain(id.long.var = "Sample", ...)
  } else if (!sample_colored && !sample_linked) {
    p <- p + ggrain::geom_rain(...)
  }

  return(p)
}

################################################################################
################################################################################


################################################################################
#' Circle of correspondance between two taxonomic levels
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#' @inheritParams tc_points_matrix
#' @param suffix_1 A suffix to add to rank_1 values (default "_1")
#' @param suffix_2 A suffix to add to rank_2 values (default "_2")
#'
#' @return A ggplot2 object
#' @export
#' @author Adrien Taudière
#' @examples
#' tc_circle(
#' Glom_otu,
#' "Genus__eukaryome_Glomero",
#' "Genus",
#' suffix_1 = "_Euk",
#' suffix_2 = "_Marjaam"
#' )
#'
#' tc_circle(
#'   Glom_otu,
#'   "Family__eukaryome_Glomero",
#'   "Family",
#'   suffix_1 = "_Euk",
#'   suffix_2 = "_Marjaam"
#' )

tc_circle <- function(physeq,
                            rank_1 = NULL,
                            rank_2 = NULL,
                            suffix_1 = "_1",
                            suffix_2 = "_2") {
  tab <- table(physeq@tax_table[, rank_1], physeq@tax_table[, rank_2])
  df_circle <- data.frame(
    from = rep(paste0(rownames(tab), suffix_1), times = ncol(tab)),
    to = rep(paste0(colnames(tab), suffix_2), each = nrow(tab)),
    value = as.vector(tab),
    stringsAsFactors = FALSE
  )

  suffix_name <- c(paste0(rownames(tab), suffix_1), paste0(colnames(tab), suffix_2))

  uniq_names <- unique(c(rownames(tab), colnames(tab)))

  grid.col <- circlize::rand_color(length(uniq_names))

  circlize::circos.par(gap.degree = 1)

  circlize::chordDiagram(df_circle)
  circlize::circos.clear()
}
################################################################################









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
tc_sankey(Glom_otu, "Family__eukaryome_Glomero", "Family")
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
    geom_stratum(width = 1 / 12, color = "grey")  +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Gender", "Dept"),
                     expand = c(.05, .05)) +
    theme_void()
}


















## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Taxonomy tree with linked correspondance
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# todo : add a width of branch in fonction of the number of ASV/sequences

physeq

ranks_1 = c(1:6)
ranks_2 = c(9:14)
collapse_taxa = TRUE

data("Glom_otu")

physeq_interm <- subset_samples(Glom_otu, sample_sums(Glom_otu) > 50000)
physeq_interm <- subset_taxa_pq(physeq_interm, taxa_sums(physeq_interm) >
                                  10000)
physeq_interm <- clean_pq(subset_taxa_pq(physeq_interm, taxa_sums(as_binary_otu_table(physeq_interm)) >
                                           30))
physeq = physeq_interm




require(ggtree)

formula_taxo_1 <- formula(paste0("~", paste0(colnames(physeq@tax_table)[ranks_1], collapse =
                                               "/")))

data_taxo_1 <- as.data.frame(physeq@tax_table[, ranks_1]) |>
  dplyr::mutate_if(is.character, as.factor)

if (collapse_taxa) {
  data_taxo_1 <- data_taxo_1 |>
    group_by(across(everything())) |>
    summarise(n = n())
}

phy_tree_1 <- ape::as.phylo.formula(formula_taxo_1, data = data_taxo_1, collapse =
                                      T)
p1 <- ggtree(phy_tree_1, layout = 'roundrect')


formula_taxo_2 <- formula(paste0("~", paste0(colnames(physeq@tax_table)[ranks_2], collapse =
                                               "/")))

data_taxo_2 <- as.data.frame(physeq@tax_table[, ranks_2]) |>
  dplyr::mutate_if(is.character, as.factor)

if (collapse_taxa) {
  data_taxo_2 <- data_taxo_2 |>
    group_by(across(everything())) |>
    summarise(n = n())
}

phy_tree_2 <- ape::as.phylo.formula(formula_taxo_2, data = data_taxo_2, collapse =
                                      T)
p2 <- ggtree(phy_tree_2, layout = 'roundrect')

# https://yulab-smu.top/treedata-book/chapter2.html#ggtree-fortify
d1 <- p1$data
d1$label <- gsub("NA", "", d1$label)

d2 <- p2$data
d2$label <- gsub("NA", "", d2$label)

## reverse x-axis and
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1 + geom_tree(data = d2, layout = 'ellipse') +
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
            data = dd)






