#' @title Formattable visualization for list_phyloseq summary
#'
#' @description
#' 
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' 
#' Create a visualization table to display the summary_table from a list_phyloseq
#' object using the formattable package with colored bars.
#'
#' @name formattable_lpq
#' @importFrom formattable formattable formatter style csscolor gradient percent proportion
#' @export
NULL


#' Create a color bar formatter for numeric columns
#'
#' @param color The color for the bar gradient
#' @param alpha Transparency for the lighter end of the gradient
#' @param na_color Color for NA values
#' @return A formattable formatter
#' @noRd
color_bar_formatter <- function(
  color = "#1a9641",
  alpha = 0.2,
  na_color = "white"
) {
  color_transp <- transp(color, alpha)
  formattable::formatter(
    "span",
    style = x ~ formattable::style(
      "font-size" = "85%",
      "display" = "inline-block",
      direction = "rtl",
      `border-radius` = "3px",
      `padding-right` = "4px",
      `background-color` = ifelse(
        is.na(x) | x == 0,
        na_color,
        formattable::csscolor(
          formattable::gradient(as.numeric(x), color_transp, color)
        )
      ),
      width = formattable::percent(
        formattable::proportion(as.numeric(x), na.rm = TRUE)
      )
    )
  )
}
#' Format factor columns with funky colored backgrounds
#'
#' @param x A factor or character vector
#' @return A formattable object with colored backgrounds
#' @export
factor_formatter <- function(x) {
  if (!is.factor(x) && !is.character(x)) {
    return(x)
  }

  # Get unique levels
  levels_x <- if (is.factor(x)) levels(x) else unique(x)

  # Generate funky colors for each level
  colors <- MiscMetabar::funky_color(length(levels_x))

  # Create a named vector of colors
  color_map <- setNames(colors, levels_x)

  # Apply formatting
  formattable::formatter("span", style = function(y) {
    formattable::style(
      display = "block",
      padding = "0 4px",
      `border-radius` = "4px",
      `background-color` = color_map[as.character(y)]
    )
  })
}
#' Create a logical formatter with colored background
#'
#' @param true_color Color for TRUE values
#' @param false_color Color for FALSE values
#' @param na_color Color for NA values
#' @return A formattable formatter
#' @noRd
logical_formatter <- function(
  true_color = "#53a134ff",
  false_color = "#c23133c2",
  na_color = "lightgray"
) {
  formattable::formatter(
    "span",
    style = x ~ formattable::style(
      display = "block",
      `border-radius` = "5px",
      padding = "2px 4px",
      `text-align` = "center",
      `background-color` = dplyr::case_when(
        is.na(x) ~ na_color,
        x == TRUE ~ transp(true_color, 0.4),
        x == FALSE ~ transp(false_color, 0.4),
        TRUE ~ na_color
      ),
      color = dplyr::case_when(
        is.na(x) ~ "gray",
        x == TRUE ~ true_color,
        x == FALSE ~ false_color,
        TRUE ~ "gray"
      ),
      `font-weight` = "bold"
    ),
    x ~ ifelse(is.na(x), "NA", ifelse(x, "✓", "✗"))
  )
}

#' Create a formattable table from list_phyloseq summary
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Display the summary_table from a list_phyloseq object as a formattable
#' table with colored bars for numeric columns and colored indicators
#' for logical columns.
#'
#' @param x (required) A list_phyloseq object.
#' @param columns (character vector, default selection of key columns) Character
#'   vector of column names to display. If NULL, displays a curated selection
#'   of columns.
#' @param bar_colors (named list, default NULL) Named list of colors for numeric
#'   columns with bars. Names should match column names. Default colors are
#'   provided for common columns.
#' @param round_digits (integer, default 1) Number of decimal places for rounding
#'   numeric columns.
#' @param void_style (logical, default FALSE) If TRUE, returns a formattable
#'   without any custom styling.
#' @param log10_transform (logical, default TRUE) If TRUE, applies log10
#'   transformation to numeric columns with a range greater than 1000.
#' @param ... Additional arguments passed to `formattable::formattable()`.
#'
#' @return A formattable object
#'
#' @details
#' This function is inspired by `MiscMetabar::formattable_pq()`.
#' Numeric columns are displayed with proportional colored bars.
#' Logical columns (has_sam_data, has_tax_table, etc.) are displayed
#' with checkmarks or X marks with colored backgrounds.
#'
#' @examples
#' \dontrun{
#' lpq <- list_phyloseq(list(data1 = data_fungi, data2 = data_fungi_mini))
#' formattable_lpq(lpq)
#'
#' # Custom columns
#' formattable_lpq(lpq,
#'   columns = c("name", "n_samples", "n_taxa", "n_sequences"))
#'
#' # Custom colors
#' formattable_lpq(lpq, bar_colors = list(
#'   n_samples = "steelblue",
#'   n_taxa = "darkgreen",
#'   n_sequences = "purple"
#' ))
#' }
#'
#' @export
formattable_lpq <- function(
  x,
  columns = c(
    "name",
    "n_samples",
    "n_taxa",
    "n_sequences",
    "mean_seq_per_sample",
    "mean_seq_per_taxon",
    "has_sam_data",
    "has_tax_table",
    "has_refseq",
    "has_phy_tree"
  ),
  bar_colors = NULL,
  round_digits = 1,
  void_style = FALSE,
  log10_transform = TRUE,
  ...
) {
  if (!requireNamespace("formattable", quietly = TRUE)) {
    stop(
      "Package 'formattable' is required. Install it with: install.packages('formattable')"
    )
  }

  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  df <- x@summary_table

  # Subset to selected columns
  df <- df[, columns, drop = FALSE]

  # Round numeric columns
  numeric_cols <- purrr::map_lgl(df, is.numeric)
  df <- df |>
    dplyr::mutate(dplyr::across(
      dplyr::where(is.numeric),
      ~ round(.x, round_digits)
    ))

  if (log10_transform) {
    # Apply log10 transformation with dplyr to numeric columns with a range > 1000
    df <- df |>
      dplyr::mutate(dplyr::across(
        dplyr::where(is.numeric),
        ~ if ((max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)) > 1000) {
          log10(.x + 1)
        } else {
          .x
        }
      ))
  }
  if (void_style) {
    return(formattable::formattable(df, ...))
  }

  # Default bar colors
  default_bar_colors <- list(
    n_samples = "#3288bd",
    n_taxa = "#b34e41ff",
    n_sequences = "#7065a1ff",
    n_occurence = "#a55c7eff",
    mean_seq_length = "#e6f598",
    sd_seq_length = "#e6f598",
    min_seq_length = "#e6f598",
    max_seq_length = "#e6f598",
    mean_seq_per_sample = "#f46d43",
    sd_seq_per_sample = "#f46d43",
    min_seq_per_sample = "#f46d43",
    max_seq_per_sample = "#f46d43",
    mean_seq_per_taxon = "#3288bd",
    sd_seq_per_taxon = "#3288bd"
  )

  # Merge with user-provided colors
  if (!is.null(bar_colors)) {
    default_bar_colors <- modifyList(default_bar_colors, bar_colors)
  }

  # Add grey for all columns not in default
  for (col in columns) {
    if (is.numeric(df[[col]]) && !(col %in% names(default_bar_colors))) {
      default_bar_colors[[col]] <- "gray"
    }
  }

  # Build formatter list
  formatter_list <- list()

  # Name column formatter (bold)
  if ("name" %in% columns) {
    formatter_list[["name"]] <- formattable::formatter(
      "span",
      style = ~ formattable::style(
        `font-weight` = "bold",
        `padding-right` = "8px"
      )
    )
  }

  # Numeric columns with color bars
  numeric_bar_cols <- intersect(
    names(default_bar_colors),
    columns[purrr::map_lgl(df[columns], is.numeric)]
  )

  for (col in numeric_bar_cols) {
    formatter_list[[col]] <- color_bar_formatter(
      color = default_bar_colors[[col]],
      alpha = 0.15
    )
  }

  # Logical columns with checkmarks
  logical_cols <- columns[purrr::map_lgl(df[columns], is.logical)]
  for (col in logical_cols) {
    formatter_list[[col]] <- logical_formatter()
  }

  # Factor/character columns with funky colors
  factor_cols <- columns[purrr::map_lgl(df[columns], function(x) {
    is.factor(x) || is.character(x)
  })]
  for (col in factor_cols) {
    formatter_list[[col]] <- factor_formatter(df[[col]])
  }
  formattable::formattable(df, formatter_list, ...)
}

#' Extended formattable for list_phyloseq with comparison info
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Create an extended formattable table that also displays comparison
#' characteristics from a list_phyloseq object.
#'
#' @param x (required) A list_phyloseq object.
#' @param show_summary (logical, default TRUE) If TRUE, show the summary table.
#' @param show_comparison (logical, default TRUE) If TRUE, show comparison info.
#' @param ... Additional arguments passed to `formattable_lpq()`.
#'
#' @return A list containing formattable objects for summary and comparison,
#'   or a single formattable if only one is requested.
#'
#' @export
formattable_lpq_full <- function(
  x,
  show_summary = TRUE,
  show_comparison = TRUE,
  ...
) {
  if (!requireNamespace("formattable", quietly = TRUE)) {
    stop(
      "Package 'formattable' is required. Install it with: install.packages('formattable')"
    )
  }

  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  result <- list()

  # Summary table
  if (show_summary) {
    result$summary <- formattable_lpq(x, ...)
  }

  # Comparison characteristics as a table
  if (show_comparison) {
    comp <- x@comparison

    # Build a comparison summary tibble
    comp_df <- tibble::tibble(
      Characteristic = c(
        "Same sample_data structure",
        "Same samples",
        "Same taxa",
        "Number of common samples",
        "Number of common taxa",
        "All have sample_data",
        "All have tax_table",
        "All have refseq",
        "All have phy_tree"
      ),
      Value = c(
        as.character(comp$same_sam_data_structure),
        as.character(comp$same_samples),
        as.character(comp$same_taxa),
        as.character(comp$n_common_samples),
        as.character(comp$n_common_taxa),
        as.character(comp$all_have_sam_data),
        as.character(comp$all_have_tax_table),
        as.character(comp$all_have_refseq),
        as.character(comp$all_have_phy_tree)
      )
    )

    # Apply formatting
    result$comparison <- formattable::formattable(
      comp_df,
      list(
        Characteristic = formattable::formatter(
          "span",
          style = ~ formattable::style(`font-weight` = "bold")
        ),
        Value = formattable::formatter(
          "span",
          style = x ~ formattable::style(
            `background-color` = dplyr::case_when(
              x == "TRUE" ~ transp("#1a9641", 0.4),
              x == "FALSE" ~ transp("#d7191c", 0.4),
              x == "NA" ~ transp("gray", 0.3),
              TRUE ~ "transparent"
            ),
            `border-radius` = "4px",
            `padding` = "2px 6px",
            `display` = "inline-block"
          )
        )
      )
    )
  }

  # Return based on what was requested
  if (show_summary && !show_comparison) {
    return(result$summary)
  } else if (!show_summary && show_comparison) {
    return(result$comparison)
  } else {
    return(result)
  }
}

#' Display shared sample_data modalities
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Create a table showing which sample_data variable modalities
#' are shared across phyloseq objects in a list_phyloseq.
#'
#' @param x (required) A list_phyloseq object.
#' @param max_modalities (integer, default 10) Maximum number of modalities to
#'   display per variable.
#'
#' @return A tibble or NULL if no shared modalities exist
#'
#' @export
#' @examples
#' lpq <- list_phyloseq(list(data1 = data_fungi, data2 = data_fungi_mini))
#' 
#' shared_mod_lpq(lpq)
#' shared_mod_lpq(lpq, 10)
shared_mod_lpq <- function(x, max_modalities = NULL) {
 
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  shared_mod <- x@comparison$shared_sam_data_modalities

  if (length(shared_mod) == 0) {
    message("No shared sample_data modalities found.")
    return(NULL)
  }

  # Build a table
  df_list <- purrr::imap(shared_mod, function(mods, var_name) {
    n_mods <- length(mods)
    if (n_mods > max_modalities && !is.null(max_modalities)) {
      mods_display <- c(
        mods[seq_len(max_modalities)],
        paste0("... (+", n_mods - max_modalities, " more)")
      )
    } else {
      mods_display <- mods
    }
    tibble::tibble(
      Variable = var_name,
      `N_shared` = n_mods,
      `Shared_modalities` = paste(mods_display, collapse = ", ")
    )
  })

  df <- dplyr::bind_rows(df_list)
  return(df)
}

#' UpSet or Venn plot of shared taxonomic values across phyloseq objects
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Create an UpSet plot (or Venn diagram) showing the shared 
#' taxonomic values at a specified rank across all phyloseq objects 
#' in a list_phyloseq.
#'
#' @param x (required) A list_phyloseq object.
#' @param tax_rank (character, required) The name of the taxonomic rank column
#'   present in the `@tax_table` slot of each phyloseq object. For example,
#'   "Genus", "Family", or "Species".
#' @param plot_type (character, default "auto") Type of plot to generate.
#'   One of "auto", "upset", or "venn". If "auto", uses Venn diagram for
#'   4 or fewer phyloseq objects, UpSet plot otherwise.
#' @param remove_na (logical, default TRUE) If TRUE, remove NA values from
#'   the taxonomic rank before computing intersections.
#' @param ... Additional arguments passed to [ComplexUpset::upset()] or
#'   [ggVennDiagram::ggVennDiagram()].
#'
#' @return A ggplot2 object (both UpSet and Venn diagrams)
#'
#' @details
#' This function extracts the unique values for the specified taxonomic rank
#' from each phyloseq object and creates a visualization showing the
#' intersections between them. UpSet plots are generally better for visualizing
#' complex intersections with more than 4 sets, while Venn diagrams work well
#' for 2-4 sets.
#'
#' @examples
#' data("enterotype", package = "phyloseq")
#' lpq <- list_phyloseq(list(fung = data_fungi, 
#'   fung_mini= data_fungi_mini, 
#'   fung_rarefy = rarefy_even_depth(data_fungi),
#'   enterotype = enterotype))
#' upset_lpq(lpq, plot_type = "upset")
#' lpq2 <- list_phyloseq(list(fung = data_fungi,
#'  fung_mini= data_fungi_mini)) 
#' upset_lpq(lpq2, tax_rank = "Family")
#'
#' @export
upset_lpq <- function(
    x,
    tax_rank= "Genus",
    plot_type = "auto",
    remove_na = TRUE,
    ...) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  plot_type <- match.arg(plot_type, c("auto", "upset", "venn"))

  # Extract unique taxonomic values for each phyloseq object
  tax_values_list <- purrr::map(x@phyloseq_list, function(pq) {
    if (is.null(phyloseq::tax_table(pq, errorIfNULL = FALSE))) {
      warning("One phyloseq object has no tax_table, skipping.")
      return(character(0))
    }

    tax_tab <- phyloseq::tax_table(pq)

    if (!(tax_rank %in% colnames(tax_tab))) {
      stop(
        "Taxonomic rank '", tax_rank, "' not found in tax_table. ",
        "Available ranks: ", paste(colnames(tax_tab), collapse = ", ")
      )
    }

    values <- as.vector(tax_tab[, tax_rank])

    if (remove_na) {
      values <- values[!is.na(values)]
    }

    unique(values)
  })

  names(tax_values_list) <- names(x@phyloseq_list)
  n_sets <- length(tax_values_list)

  if (plot_type == "auto") {
    plot_type <- if (n_sets <= 4) "venn" else "upset"
  }

  if (plot_type == "venn") {
    if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
      stop(
        "Package 'ggVennDiagram' is required for Venn diagrams. ",
        "Install it with: install.packages('ggVennDiagram')"
      )
    }
    if (n_sets > 7) {
      warning(
        "Venn diagrams with more than 7 sets may not display well. ",
        "Consider using plot_type = 'upset'."
      )
    }

    p <- ggVennDiagram::ggVennDiagram(tax_values_list, ...) +
      ggplot2::labs(title = paste("Shared", tax_rank, "across phyloseq objects"))

    return(p)
  } else {
    # UpSet plot using ComplexUpset
    if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
      stop(
        "Package 'ComplexUpset' is required for UpSet plots. ",
        "Install it with: install.packages('ComplexUpset')"
      )
    }

    # ComplexUpset requires a data frame with binary columns for each set
    # Get all unique values across all sets
    all_values <- unique(unlist(tax_values_list))

    # Create a data frame with one row per unique value
    # and one column per phyloseq object (TRUE/FALSE membership)
    upset_data <- data.frame(
      value = all_values,
      stringsAsFactors = FALSE
    )

    # Add binary columns for each phyloseq object
    for (set_name in names(tax_values_list)) {
      upset_data[[set_name]] <- all_values %in% tax_values_list[[set_name]]
    }

    # Get the set names for the intersect parameter
    set_names <- names(tax_values_list)

    # Create the plot
    p <- ComplexUpset::upset(
      upset_data,
      intersect = set_names,
      name = tax_rank,
      sort_sets = "descending",
      sort_intersections = "descending",
      ...
    )

    return(p)
  }
}
