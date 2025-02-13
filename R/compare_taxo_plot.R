################################################################################
#' Barchart of ratio to compare 2 taxonomic ranks
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful to compare taxonomy from two different source/db/algo side-by-side
#'
#' @inheritParams tc_points_matrix
#' @param color_rank Define the taxonomic rank for color
#'   as the number or the name of the column in tax_table slot
#' @param point_size Size of points
#' @param point_alpha Transparency of points
#' @param log10trans (Logical, default TRUE) Do the abundancy is log10
#'   transformed?
#'
#' @return A ggplot2 object
#' @export
#' @author Adrien Taudière
#' @examples
#'
#' tc_bar(subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000), rank_1 = 5, rank_2 = 13, color_rank = 3)
#' \dontrun{
#' tc_bar(Glom_otu, rank_1 = 5, rank_2 = 13, color_rank = 3)
#' tc_bar(as_binary_otu_table(Glom_otu), rank_1 = 5, rank_2 = 13, color_rank = 3, log10trans = FALSE)
#' tc_bar(Glom_otu,
#'   rank_1 = "Genus",
#'   rank_2 = "Genus__eukaryome_Glomero",
#'   color_rank = "Family"
#' )
#' }
tc_bar <- function(physeq,
                   rank_1,
                   rank_2,
                   color_rank,
                   point_size = 0.3,
                   point_alpha = 0.3,
                   merge_sample_by = NULL,
                   log10trans = TRUE) {
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

  if (log10trans) {
    psm2 <- psm2 |>
      mutate(Abundance = log10(Abundance))
  }
  p <- ggplot(psm2, aes(
    x = .data[[ranks2[1]]],
    y = Abundance,
    fill = .data[[ranks2[[2]]]],
    label = .data[[ranks2[1]]]
  )) +
    stat_summary(
      fun.y = mean,
      position = "dodge",
      geom = "bar"
    ) +
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
        y = -Abundance,
        fill = .data[[ranks1[[2]]]]
      ),
      fun.y = mean,
      position = "dodge",
      geom = "bar"
    ) +
    stat_summary(
      aes(x = .data[[ranks1[1]]], y = -Abundance),
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
    geom_jitter(aes(x = .data[[ranks1[1]]], y = -Abundance),
      size =
        point_size, alpha = point_alpha
    )

  if (log10trans) {
    p <- p + xlab("Abundance (log10)")
  }
  return(p)
}
################################################################################

################################################################################
#' Matrix of point to compare two taxonomic ranks
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'  Useful to compare taxonomy from two different source/db/algo side-by-side
#'
#' @param physeq (required): a \code{\link[phyloseq]{phyloseq-class}} object obtained
#'   using the `phyloseq` package.
#' @param rank_1 Define the first taxonomic rank as the number or the name
#'    of the column in tax_table slot
#' @param rank_2 Define the second taxonomic rank as the number or the name
#'    of the column in tax_table slot
#' @param color_1 Color for rank_1 values
#' @param color_2 Color for rank_1 values
#' @param stat_across_sample Either "mean" or "sum". Set how the abundance is
#'   computed across samples.
#' @param merge_sample_by a vector to determine
#'   which samples to merge using [merge_samples2()] function.
#'   Need to be in \code{physeq@sam_data}
#'
#' @return A ggplot2 object
#' @export
#' @author Adrien Taudière
#' @examples
#'
#' tc_points_matrix(
#'   subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000),
#'   "Order", "Order__eukaryome_Glomero"
#' )
#' \dontrun{
#' tc_points_matrix(Glom_otu, 6, 14)
#' tc_points_matrix(Glom_otu, 4, 12)
#' tc_points_matrix(Glom_otu, 4, 12, stat_across_sample = "mean")
#'
#' Glom_otu@sam_data$unique_value <- rep("samp", nsamples(Glom_otu))
#' tc_points_matrix(as_binary_otu_table(Glom_otu), 5, 13,
#'   stat_across_sample = "sum", merge_sample_by = "unique_value"
#' )
#' tc_points_matrix(as_binary_otu_table(Glom_otu), 5, 13,
#'   stat_across_sample = "mean"
#' )
#' tc_points_matrix(Glom_otu, 5, 13,
#'   stat_across_sample = "mean"
#' )
#' }
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
    summarise(
      mean_ab = mean(Abundance),
      sum_ab = sum(Abundance)
    ) |>
    mutate(Rank = .data[[ranks1[[1]]]])

  psm2_summary_taxo2 <- psm2 %>%
    group_by(.data[[ranks2[[1]]]]) %>%
    summarise(
      mean_ab = mean(Abundance),
      sum_ab = sum(Abundance)
    ) |>
    mutate(Rank = .data[[ranks2[[1]]]])

  psm2_summary <- full_join(psm2_summary_taxo1,
    psm2_summary_taxo2,
    by = join_by(Rank == Rank)
  )

  if (stat_across_sample == "mean") {
    p <- ggplot(psm2_summary, aes(
      x = .data[[ranks1[[1]]]],
      y = .data[[ranks2[[1]]]],
      size = mean_ab.x
    )) +
      geom_point(
        shape = 21,
        color = color_1,
        fill = alpha(color_1, 0.2)
      ) +
      geom_point(
        aes(size = mean_ab.y),
        color = color_2,
        shape = 21, ,
        fill = alpha(color_2, 0.2)
      )
  } else if (stat_across_sample == "sum") {
    p <- ggplot(psm2_summary, aes(
      x = .data[[ranks1[[1]]]],
      y = .data[[ranks2[[1]]]],
      size = sum_ab.x
    )) +
      geom_point(
        shape = 21,
        color = color_1,
        fill = alpha(color_1, 0.2)
      ) +
      geom_point(
        aes(size = sum_ab.y),
        color = color_2,
        shape = 21, ,
        fill = alpha(color_2, 0.2)
      )
  } else {
    stop("Param stat_across_sample must be set to mean or sum !")
  }
  p <- p + theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))

  return(p)
}
################################################################################

################################################################################
#' Rainplot of the nb taxa assigned (not NA)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful to compare taxonomy from two different source/db/algo
#'
#' @inheritParams tc_points_matrix
#' @param ranks The ranks to include in the rainplot. If left to NULL, all
#'   ranks are used. Each rank can be defined either by integer for the index
#'   or by its full names (exacly matching the colnames of the `tax_table` slot).
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
#' rainplot_taxo_na(subset_taxa_pq(Glom_otu, taxa_sums(Glom_otu) > 5000),
#'   ranks = c("Family", "Family__eukaryome_Glomero")
#' )
#' \dontrun{
#' rainplot_taxo_na(Glom_otu)
#'
#' Glom_otu@sam_data$tmt_type <- paste0(Glom_otu@sam_data$Tmt, "_", Glom_otu@sam_data$Type)
#' rainplot_taxo_na(
#'   Glom_otu,
#'   merge_sample_by = "tmt_type",
#'   sample_colored = TRUE,
#'   sample_linked = TRUE
#' )
#' rainplot_taxo_na(Glom_otu, ranks = c(4, 12), rain.side = "f1x1")
#' rainplot_taxo_na(
#'   Glom_otu,
#'   ranks = c(6, 14),
#'   rain.side = "f1x1",
#'   sample_linked = TRUE
#' ) +
#'   theme(legend.position = "none")
#' }
rainplot_taxo_na <- function(physeq,
                             ranks = NULL,
                             min_nb_seq = 0,
                             merge_sample_by = NULL,
                             sample_colored = FALSE,
                             sample_linked = FALSE,
                             ...) {
  verify_pq(physeq)

  if (is.null(ranks)) {
    ranks <- seq(1, ncol(physeq@tax_table))
  }
  if (!is.null(merge_sample_by)) {
    physeq <- clean_pq(merge_samples2(physeq, merge_sample_by))
  }

  rank_names <- colnames(physeq@tax_table[, ranks])

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
    coord_flip() +
    ylab("Proportion of known (non NA) value") +
    xlab("Rank")

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
#' Circle of correspondence between two taxonomic levels
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful to compare taxonomy from two different source/db/algo side-by-side
#'
#' @inheritParams tc_points_matrix
#' @param suffix_1 A suffix to add to rank_1 values (default "_1")
#' @param suffix_2 A suffix to add to rank_2 values (default "_2")
#'
#' @return A ggplot2 object
#' @export
#' @author Adrien Taudière
#' @examples
#' tc_circle(
#'   Glom_otu,
#'   "Genus__eukaryome_Glomero",
#'   "Genus",
#'   suffix_1 = "_Euk",
#'   suffix_2 = "_Marjaam"
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
