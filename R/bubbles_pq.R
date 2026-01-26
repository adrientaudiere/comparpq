#' Bubble plot of phyloseq object with observablehq
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#'
#' Creates an interactive bubble plot visualization of taxa abundances from a
#' phyloseq object using an Observable notebook.
#'
#' @inheritParams tc_points_matrix
#' @param rank_label (character, default "Taxa") The name of the column in the
#'   @tax_table slot to label the points. If set to "Taxa", the taxa name is used.
#' @param rank_color (character, default "Family") The name of the column in the
#'   @tax_table slot to color the points.
#' @param categorical_scheme (character, default "d3.schemeCategory10") A color
#'   scheme from [d3js](https://d3js.org/d3-scale-chromatic/categorical).
#' @param label_color (character, default "grey10") Color for the label text.
#' @param value_color (character, default "grey20") Color for the value text.
#' @param label_size (integer, default 8) Font size for the label.
#' @param log1ptransform (logical, default FALSE) If TRUE, the number of
#'   sequences is log1p transformed.
#' @param min_nb_seq (integer, default 0) Minimum number of sequences to filter
#'   out taxa with low abundance.
#' @param randomize (logical, default FALSE) If TRUE, shuffles the order of the
#'   taxa.
#' @param seed (integer, default 32) Seed for the randomization.
#' @param width (integer, default 600) The notebook width for visualization in
#'   pixels.
#' @param include (character vector, default c("TitleCell", "key", "chart"))
#'   You can modify include to remove the title and/or the legend, for example
#'   using `include = c("chart")`.
#' @param notebook (character, default "https://observablehq.com/d/d755af3197af2320")
#'   You can change the notebook url if you know what you are doing.
#' @param return_dataframe (logical, default FALSE) If TRUE, the plot is not
#'   returned, but the resulting dataframe to plot is returned.
#' @param title (character, default "") Title of the plot.
#' @param title_size (character, default "22px") Font size for the title.
#'
#' @return A htmlwidget
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' bubbles_pq(physeq = data_fungi, rank_color = "Class", randomize = TRUE)
#' bubbles_pq(physeq = data_fungi, min_nb_seq = 10000, categorical_scheme = "d3.schemePastel1")
#' bubbles_pq(physeq = data_fungi, min_nb_seq = 1000, categorical_scheme = "d3.schemeTableau10", label_color = "purple", value_color = "white")
#' bubbles_pq(physeq = data_fungi, rank_color = "Order", include = c("chart"), log1ptransform = TRUE)
#' bubbles_pq(physeq = data_fungi, rank_label = "Genus", rank_color = "Class", include = c("chart"), randomize = TRUE, seed = 3)
#' bubbles_pq(physeq = as_binary_otu_table(data_fungi), min_nb_seq = 20, title = "Nb of samples per OTU. <br> Only OTU present in more than 20 samples are shown", title_size = "12px")
bubbles_pq <- function(physeq,
                       rank_label = "Taxa",
                       rank_color = "Family",
                       categorical_scheme = "d3.schemeCategory10",
                       label_color = "grey10",
                       value_color = "grey20",
                       label_size = 8,
                       log1ptransform = FALSE,
                       min_nb_seq = 0,
                       randomize = FALSE,
                       seed = 32,
                       width = 600,
                       include = c("TitleCell", "key", "chart"),
                       notebook = "https://observablehq.com/d/d755af3197af2320",
                       return_dataframe = FALSE,
                       title = "",
                       title_size = "22px") {
  if (rank_label == "Taxa") {
    label <- taxa_names(physeq)
  } else {
    label <- as.vector(physeq@tax_table[, rank_label])
  }

  df <- data.frame(
    "value" = taxa_sums(physeq),
    "label" = label,
    "color" = fac2col(physeq@tax_table[, rank_color]),
    "rank_value_color" = as.vector(physeq@tax_table[, rank_color])
  )

  df$id <- paste0(df$label, ".", df$rank_value_color, ".", df$label)

  if (min_nb_seq > 0) {
    df <- df |>
      dplyr::filter(value > min_nb_seq)
  }

  if (log1ptransform) {
    df$value <- log1p(df$value)
  }

  if (randomize) {
    set.seed(seed)
    df <- df[sample(nrow(df)), ]
  }

  if (return_dataframe) {
    return(df)
  }
  robservable::robservable(
    notebook = notebook,
    include = include,
    input = list(
      data = df,
      categoricalScheme = categorical_scheme,
      colorLabel = label_color,
      colorValue = value_color,
      title = title,
      titleSize = title_size
    ),
    width = width
  )
}
