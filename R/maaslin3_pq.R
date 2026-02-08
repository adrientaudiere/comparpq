#' Run MaAsLin3 differential abundance analysis on a phyloseq object
#'
#' #' #TODO VERY experimental
#' 
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' A wrapper around `maaslin3::maaslin3()` for phyloseq objects. Optionally
#' includes the number of reads (library size) as a covariate in the model
#' to account for differences in sequencing depth.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param formula (character, required) A formula string for the model
#'   (e.g., `"~ Height"` or `"~ Height + Time"`).
#' @param reference (named list, default NULL) Reference levels for categorical
#'   variables with more than 2 levels. Format: `list(varname = "ref_level")`.
#'   For example: `list(Height = "Low")` sets "Low" as the reference for Height.
#'   Variables not in this list will use their first level as reference.
#' @param correction_for_sample_size (logical, default TRUE) If TRUE, adds
#'   `nb_seq` (library size) to the formula to control for sequencing depth.
#' @param output (character, default "res_maaslin3") Output directory for
#'   MaAsLin3 results.
#' @param ... Additional arguments passed to `maaslin3::maaslin3()`.
#'
#' @return The result object from `maaslin3::maaslin3()`, containing:
#'   - `results`: Data frame with differential abundance results
#'   - `results_ordered`: Results ordered by significance
#'   - Other MaAsLin3 output components
#'
#' @details
#' MaAsLin3 requires categorical variables with more than 2 levels to have
#' a defined reference level. This function handles this by:
#' 1. Converting character variables to factors
#' 2. Setting reference levels via the `reference` parameter
#' 3. Using `relevel()` to set the specified reference as the first level
#'
#' The `correction_for_sample_size` option adds library size as a covariate,
#' which can help account for compositional effects and sequencing depth
#' differences between samples.
#'
#' @export
#' @author Adrien Taudière
#'
#' @seealso [MiscMetabar::aldex_pq()], [ancombc_lpq()], [gg_maaslin3_plot()]
#'
#' @references
#' Nickols WA, et al. (2024). MaAsLin 3: Refining and extending generalized
#' multivariable linear models for meta-omic association discovery. bioRxiv.
#'
#' @examples
#' # Basic usage with a binary factor
#' data_fungi_high_low <- subset_samples(
#' data_fungi, Height %in% c("High", "Low")
#' )
#' res <- maaslin3_pq(data_fungi_high_low, formula = "~ Height")
#' res$results
#'
#' # Specify reference level for multi-level factor
#' res <- maaslin3_pq(
#'   data_fungi,
#'   formula = "~ Height",
#'   reference = list(Height = "Low")
#' )
#'
#' # Multiple variables with references
#' res <- maaslin3_pq(
#'   data_fungi,
#'   formula = "~ Height + Site",
#'   reference = list(Height = "Low", Site = "A")
#' )
#'
#' # Without library size correction
#' res <- maaslin3_pq(
#'   data_fungi,
#'   formula = "~ Height",
#'   correction_for_sample_size = FALSE
#' )
#'
#' # Plot results
#' gg_maaslin3_plot(res, type = "volcano")
#'             
#' # Full example with GlobalPatterns dataset                                                                                                
#'  data("GlobalPatterns")                                                                          
#'                                                                                                  
#'  # Subset to two very different environments: Feces vs Soil                                      
#'  gp_subset <- subset_samples(GlobalPatterns, SampleType %in% c("Feces", "Soil"))                 
#'                                                                                                  
#'  # Agglomerate at Phylum level for clearer results                                               
#'  gp_phylum <- tax_glom(gp_subset, taxrank = "Phylum", NArm = FALSE)                              
#'                                                                                                  
#'  # Run MaAsLin3 with Soil as reference                                                           
#'  res <- maaslin3_pq(                                                                             
#'    gp_phylum,                                                                                    
#'    formula = "~ SampleType",                                                                     
#'    reference = list(SampleType = "Soil"),                                                        
#'    output = "output/maaslin3_example"  ,
#'   correction_for_sample_size = FALSE                                                          
#'  )                                                                                               
#'                                                                                                  
#'  # Plot results
#'  gg_maaslin3_plot(res, type = "volcano", signif_threshold = 0.1)
#'  gg_maaslin3_plot(res, type = "forest", top_n = 15) 
maaslin3_pq <- function(
  physeq,
  formula,
  reference = NULL,
  correction_for_sample_size = TRUE,
  output = "res_maaslin3",
  ...
) {
  if (!requireNamespace("maaslin3", quietly = TRUE)) {
    stop(
      "Package 'maaslin3' is required. Install it with:\n",
      "  BiocManager::install('biobakery/maaslin3')"
    )
  }

  verify_pq(physeq)
  physeq <- taxa_as_columns(physeq)

  # Get metadata with library size info
  metadata <- phyloseq::sample_data(
    add_info_to_sam_data(physeq, add_nb_otu = FALSE)
  ) |>
    unclass() |>
    as.data.frame()

  rownames(metadata) <- phyloseq::sample_names(physeq)

  # Extract variable names from formula
  formula_vars <- all.vars(as.formula(formula))

  # Convert character columns to factors and set reference levels
  for (var in formula_vars) {
    if (var %in% colnames(metadata)) {
      col <- metadata[[var]]

      # Convert character to factor
      if (is.character(col)) {
        metadata[[var]] <- as.factor(col)
        col <- metadata[[var]]
      }

      # Set reference level if specified
      if (is.factor(col) && !is.null(reference) && var %in% names(reference)) {
        ref_level <- reference[[var]]
        if (!ref_level %in% levels(col)) {
          stop(
            "Reference level '", ref_level, "' not found in variable '", var,
            "'. Available levels: ", paste(levels(col), collapse = ", ")
          )
        }
        metadata[[var]] <- relevel(col, ref = ref_level)
      }
    }
  }

  # Get OTU table
  otutable <- phyloseq::otu_table(physeq) |>
    as.data.frame()

  # Add library size to formula if requested
  if (correction_for_sample_size) {
    formula <- paste0(formula, " + nb_seq")
  }

  # Run MaAsLin3
  res_maaslin3 <- maaslin3::maaslin3(
    input_data = otutable,
    input_metadata = metadata,
    formula = formula,
    output = output,
    ...
  )

  return(res_maaslin3)
}


#' Plot MaAsLin3 results
#' 
#' #TODO VERY experimental
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Creates ggplot2 visualizations of MaAsLin3 differential abundance results.
#' Supports summary plots (default), volcano plots, and forest plots.
#'
#' @param res (list, required) The full result object from `maaslin3_pq()`.
#'   For types other than "summary", can also be directly the `$results` data frame.
#' @param type (character, default "summary") Type of plot. One of:
#'   - `"summary"`: Default MaAsLin3 summary plot with coefficients and heatmap
#'     (uses `maaslin3::maaslin_plot_results()` internally)
#'   - `"volcano"`: Effect size (coef) vs -log10(p-value)
#'   - `"forest"`: Effect sizes with confidence intervals
#' @param model (character, default "abundance") Which model results to plot.
#'   One of `"abundance"`, `"prevalence"`, or `"both"`. Used for volcano and
#'   forest plots.
#' @param metadata_filter (character, default NULL) Filter results to a
#'   specific metadata variable. If NULL, uses the first non-intercept variable.
#'   Used for volcano and forest plots.
#' @param signif_threshold (numeric, default 0.1) Significance threshold
#'   for q-value (FDR-corrected). For summary plot, this is passed to
#'   `max_significance`. For other plots, used to color significant features.
#' @param use_qval (logical, default TRUE) If TRUE, use q-values (FDR-corrected)
#'   for significance. If FALSE, use raw p-values. Used for volcano and forest.
#' @param top_n (integer, default 25) For summary plot, number of top features
#'   to display (`summary_plot_first_n`). For forest plot, show only the top N
#'   features by absolute effect size.
#' @param show_labels (logical, default FALSE) For volcano plots, show
#'   feature labels for significant points.
#' @param point_size (numeric, default 2) Size of points in volcano plot.
#' @param colors (character vector, default NULL) Colors for significant/
#'   non-significant points. If NULL, uses c("grey60", "firebrick3").
#' @param show_annotations (logical, default TRUE) If TRUE, adds annotations
#'   indicating which group has higher values on each side of the plot.
#'   For simple formulas like `~ Height`, shows "More in High" vs "More in Low".
#' @param reference_label (character, default NULL) Label for the reference
#'   group. If NULL, attempts to extract from the metadata column name.
#' @param coef_plot_vars (character vector, default NULL) For summary plot only.
#'   Variables to include in coefficient plot. If NULL, uses all non-intercept
#'   variables.
#' @param heatmap_vars (character vector, default NULL) For summary plot only.
#'   Variables to include in heatmap. If NULL, determined automatically.
#' @param normalization (character, default "TSS") Normalization method used.
#'   Only needed for summary plot. One of "TSS", "CLR", or "NONE".
#' @param transform (character, default "LOG") Transform method used.
#'   Only needed for summary plot. One of "LOG", "PLOG", or "NONE".
#'
#' @return A ggplot2 object.
#'
#' @details
#' **Summary plot** (default): Uses `maaslin3::maaslin_plot_results()` to create
#' the standard MaAsLin3 visualization with sorted per-feature coefficients
#' plotted with standard errors, and a heatmap for additional variables. This
#' requires the full maaslin3 result object (not just the results data frame).
#'
#' **Volcano plot**: Shows the relationship between effect size (coefficient)
#' and statistical significance. Points above the horizontal dashed line
#' are significant at the specified threshold. Points are colored by
#' significance.
#'
#' **Forest plot**: Shows effect sizes with 95
#' percent confidence intervals for the top features. Features are ordered by effect size. Significant
#' associations are highlighted.
#'
#' The `coef` in MaAsLin3 abundance models represents log2 fold change:
#' a one-unit change in the variable corresponds to a 2^coef fold change
#' in relative abundance.
#'
#' @export
#' @author Adrien Taudière
#'
#' @seealso [maaslin3_pq()], [gg_aldex_plot()]
#'
#' @examples
#' \dontrun{
#' # Run MaAsLin3
#' res <- maaslin3_pq(
#'   data_fungi,
#'   formula = "~ Height",
#'   reference = list(Height = "Low")
#' )
#'
#' # Summary plot (default) - uses maaslin3's native visualization
#' gg_maaslin3_plot(res)
#'
#' # Volcano plot
#' gg_maaslin3_plot(res, type = "volcano")
#'
#' # Forest plot of top 15 features
#' gg_maaslin3_plot(res, type = "forest", top_n = 15)
#'
#' # Plot prevalence model results
#' gg_maaslin3_plot(res, type = "volcano", model = "prevalence")
#'
#' # Customize significance threshold
#' gg_maaslin3_plot(res, type = "volcano", signif_threshold = 0.05, show_labels = TRUE)
#' 
#'  # Complete Example with HMP2 Dataset                                                              
#'                                                                                                  
#'  library(maaslin3)                                                                                                              
#'  # ============================================================                                  
#'  # Convert maaslin3 default HMP2 dataset to phyloseq                                             
#'  # ============================================================                                  
#'                                                                                                  
#'  # Read the HMP2 default dataset from maaslin3 package                                           
#'  taxa_table_name <- system.file("extdata", "HMP2_taxonomy.tsv", package = "maaslin3")            
#'  taxa_table <- read.csv(taxa_table_name, sep = "\t", row.names = 1)                              
#'                                                                                                  
#'  metadata_name <- system.file("extdata", "HMP2_metadata.tsv", package = "maaslin3")              
#'  metadata <- read.csv(metadata_name, sep = "\t", row.names = 1)                                  
#'                                                                                                  
#'  # Set factor levels                                                                             
#'  metadata$antibiotics <- factor(metadata$antibiotics, levels = c("No", "Yes"))                   
#'                                                                                                  
#'  # Create phyloseq object (HMP2 data has samples as rows, taxa as columns)                       
#'  otu <- otu_table(as.matrix(taxa_table), taxa_are_rows = FALSE)                                  
#'  sam <- sample_data(metadata)                                                                    
#'  species_names <- colnames(taxa_table)                                                           
#'  tax_df <- data.frame(                                                                           
#'    Species = species_names,                                                                      
#'    Genus = sapply(strsplit(species_names, "_"), \(x) x[1]),                                      
#'    row.names = species_names                                                                     
#'  )                                                                                               
#'  tax <- tax_table(as.matrix(tax_df))                                                             
#'  physeq_hmp2 <- phyloseq(otu, sam, tax)                                                          
#'                                                                                                  
#'  # ============================================================                                  
#'  # Run MaAsLin3 analysis                                                                         
#'  # ============================================================                                  
#'                                                                                                  
#'  res <- maaslin3_pq(                                                                             
#'    physeq_hmp2,                                                                                  
#'    formula = "~ antibiotics",                                                                    
#'    reference = list(antibiotics = "No"),                                                         
#'    output = "output/maaslin3_hmp2",                                                              
#'    correction_for_sample_size = FALSE                                                            
#'  )                                                                                               
#'                                                                                                  
#'  # ============================================================                                  
#'  # Compare plots                                                                                 
#'  # ============================================================                                  
#'                                                                                                  
#'  # Summary plot (NEW DEFAULT) - uses maaslin3's native visualization                             
#'  gg_maaslin3_plot(res)                                                                           
#'                                                                                                  
#'  # Alternative plot types                                                                        
#'  gg_maaslin3_plot(res, type = "volcano")                                                         
#'  gg_maaslin3_plot(res, type = "forest", top_n = 15)
#' }
#' 
gg_maaslin3_plot <- function(
  res,
  type = c("summary", "volcano", "forest"),
  model = c("abundance", "prevalence", "both"),
  metadata_filter = NULL,
  signif_threshold = 0.1,
  use_qval = TRUE,
  top_n = 25,
  show_labels = FALSE,
  point_size = 2,
  colors = NULL,
  show_annotations = TRUE,
  reference_label = NULL,
  coef_plot_vars = NULL,
  heatmap_vars = NULL,
  normalization = "TSS",
  transform = "LOG"
) {
  type <- match.arg(type)
  model <- match.arg(model)

  # Handle summary plot type using maaslin3's native plotting

  if (type == "summary") {
    # Check required components - maaslin3 uses "metadata" for unstandardized metadata
    required_components <- c("transformed_data", "metadata", "fit_data_abundance")
    if (!is.list(res) || !all(required_components %in% names(res))) {
      stop(
        "For type='summary', 'res' must be the full maaslin3 result object from maaslin3_pq().\n",
        "Required components: ", paste(required_components, collapse = ", ")
      )
    }

    # Create temporary output directory
    temp_output <- tempfile(pattern = "maaslin3_plot_")
    dir.create(temp_output, showWarnings = FALSE)
    dir.create(file.path(temp_output, "figures"), showWarnings = FALSE)

    # Call maaslin_plot_results
    # Note: maaslin3 result stores unstandardized metadata in $metadata
    plot_res <- maaslin3::maaslin_plot_results(
      output = temp_output,
      transformed_data = res$transformed_data,
      unstandardized_metadata = res$metadata,
      fit_data_abundance = res$fit_data_abundance,
      fit_data_prevalence = res$fit_data_prevalence,
      normalization = normalization,
      transform = transform,
      max_significance = signif_threshold,
      plot_summary_plot = TRUE,
      summary_plot_first_n = top_n,
      coef_plot_vars = coef_plot_vars,
      heatmap_vars = heatmap_vars,
      plot_associations = FALSE
    )

    # Clean up temp directory
    unlink(temp_output, recursive = TRUE)

    return(plot_res$summary_plot$final)
  }


  # Extract results data frame from maaslin3 output
  # MaAsLin3 stores results in nested structure:
  #   $fit_data_abundance$results and $fit_data_prevalence$results
  df <- NULL

  if (is.data.frame(res)) {
    # Direct data frame passed
    df <- res
  } else if (is.list(res)) {
    # Try different maaslin3 output structures
    if (model == "abundance" || model == "both") {
      if ("fit_data_abundance" %in% names(res) &&
          "results" %in% names(res$fit_data_abundance)) {
        df <- res$fit_data_abundance$results
        if (!is.null(df)) df$model <- "abundance"
      }
    }

    if (model == "prevalence" || model == "both") {
      prev_df <- NULL
      if ("fit_data_prevalence" %in% names(res) &&
          "results" %in% names(res$fit_data_prevalence)) {
        prev_df <- res$fit_data_prevalence$results
        if (!is.null(prev_df)) prev_df$model <- "prevalence"
      }

      if (model == "prevalence") {
        df <- prev_df
      } else if (model == "both" && !is.null(prev_df)) {
        df <- rbind(df, prev_df)
      }
    }

    # Fallback: try $results directly (older structure or custom)
    if (is.null(df) && "results" %in% names(res)) {
      df <- res$results
    }
  }

  if (is.null(df) || nrow(df) == 0) {
    stop(
      "'res' must be a maaslin3 result object or a data frame.\n",
      "Expected structure: res$fit_data_abundance$results or res$fit_data_prevalence$results"
    )
  }

  # Determine p-value column
  pval_col <- if (use_qval) "qval_individual" else "pval_individual"
  if (!pval_col %in% colnames(df)) {
    # Fallback for different column names
    pval_col <- if (use_qval) {
      grep("qval", colnames(df), value = TRUE)[1]
    } else {
      grep("pval", colnames(df), value = TRUE)[1]
    }
  }

  if (is.na(pval_col) || !pval_col %in% colnames(df)) {
    stop("Cannot find p-value/q-value column in results.")
  }

  # Filter by metadata variable
  if (!is.null(metadata_filter)) {
    df <- df[df$metadata == metadata_filter, ]
  } else if ("metadata" %in% colnames(df)) {
    # Exclude intercept and nb_seq (library size covariate)
    available_meta <- unique(df$metadata)
    available_meta <- available_meta[!available_meta %in% c("(Intercept)", "nb_seq")]
    if (length(available_meta) > 0) {
      metadata_filter <- available_meta[1]
      df <- df[df$metadata == metadata_filter, ]
    }
  }

  if (nrow(df) == 0) {
    stop("No results after filtering.")
  }

  # Add significance column
  df$significant <- df[[pval_col]] < signif_threshold
  df$neg_log10_p <- -log10(df[[pval_col]])

  # Handle infinite values
  max_finite <- max(df$neg_log10_p[is.finite(df$neg_log10_p)], na.rm = TRUE)
  df$neg_log10_p[!is.finite(df$neg_log10_p)] <- max_finite * 1.1

  # Set colors
  if (is.null(colors)) {
    colors <- c("FALSE" = "grey60", "TRUE" = "firebrick3")
  } else {
    names(colors) <- c("FALSE", "TRUE")
  }

  # Parse metadata to get comparison and reference labels for annotations
  comparison_label <- NULL
  ref_label <- reference_label

  if (show_annotations && !is.null(metadata_filter)) {
    # metadata_filter is like "HeightHigh" or "SampleTypeFeces"
    # Try to extract variable name and comparison level
    # Common patterns: VariableLevel (e.g., HeightHigh, SampleTypeFeces)

    # Check if there's a 'value' column (full maaslin3 output)
    if ("value" %in% colnames(df)) {
      comparison_label <- unique(df$value)[1]
      # Extract variable name by removing the value from metadata
      var_name <- gsub(paste0(comparison_label, "$"), "", metadata_filter)
      if (is.null(ref_label)) {
        ref_label <- paste0("ref. ", var_name)
      }
    } else {
      # Simplified output - try to parse metadata column
      # Look for common patterns like "HeightHigh", "SampleTypeFeces"
      meta_str <- metadata_filter

      # Try to find where the variable name ends and level begins

      # Common variable names: Height, SampleType, Treatment, Group, etc.
      common_vars <- c(
        "SampleType", "Treatment", "Group", "Condition", "Status",
        "Height", "Time", "Site", "Location", "Diet", "Age", "Sex"
      )

      for (var in common_vars) {
        if (startsWith(meta_str, var) && nchar(meta_str) > nchar(var)) {
          comparison_label <- substring(meta_str, nchar(var) + 1)
          if (is.null(ref_label)) {
            ref_label <- paste0("ref. ", var)
          }
          break
        }
      }

      # If no match found, use generic labels
      if (is.null(comparison_label)) {
        comparison_label <- meta_str
        if (is.null(ref_label)) {
          ref_label <- "reference"
        }
      }
    }
  }

  # Build plot based on type
  if (type == "volcano") {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$coef, y = .data$neg_log10_p, color = .data$significant)
    ) +
      ggplot2::geom_point(size = point_size, alpha = 0.7) +
      ggplot2::geom_hline(
        yintercept = -log10(signif_threshold),
        linetype = "dashed",
        color = "grey40"
      ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey40") +
      ggplot2::scale_color_manual(
        values = colors,
        labels = c("FALSE" = "Not significant", "TRUE" = "Significant"),
        name = paste0(if (use_qval) "q" else "p", "-value < ", signif_threshold)
      ) +
      ggplot2::labs(
        x = "Coefficient (log2 fold change)",
        y = paste0("-log10(", if (use_qval) "q" else "p", "-value)"),
        title = if (show_annotations && !is.null(comparison_label)) {
          paste0("MaAsLin3: ", comparison_label, " vs ", ref_label)
        } else {
          paste0("MaAsLin3 Volcano Plot", if (!is.null(metadata_filter)) paste0(" - ", metadata_filter) else "")
        }
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")

    if (show_labels) {
      df_signif <- df[df$significant, ]
      if (nrow(df_signif) > 0) {
        p <- p + ggplot2::geom_text(
          data = df_signif,
          ggplot2::aes(label = .data$feature),
          hjust = -0.1,
          vjust = 0.5,
          size = 3,
          check_overlap = TRUE
        )
      }
    }

    # Add annotations for volcano plot
    if (show_annotations && !is.null(comparison_label)) {
      x_range <- range(df$coef, na.rm = TRUE)
      y_max <- max(df$neg_log10_p, na.rm = TRUE)

      p <- p + ggplot2::annotate(
        "text",
        x = x_range[2] * 0.7,
        y = y_max * 0.95,
        label = paste0("\u2191 in ", comparison_label),
        hjust = 0.5,
        size = 3.5,
        color = "darkred",
        fontface = "bold"
      ) + ggplot2::annotate(
        "text",
        x = x_range[1] * 0.7,
        y = y_max * 0.95,
        label = paste0("\u2191 in ", ref_label),
        hjust = 0.5,
        size = 3.5,
        color = "darkblue",
        fontface = "bold"
      )
    }

  } else if (type == "forest") {
    # Select top N by absolute coefficient
    df <- df[order(abs(df$coef), decreasing = TRUE), ]
    df <- utils::head(df, top_n)
    df$feature <- factor(df$feature, levels = df$feature[order(df$coef)])

    # Calculate CI
    df$ci_low <- df$coef - 1.96 * df$stderr
    df$ci_high <- df$coef + 1.96 * df$stderr

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$coef,
        y = .data$feature,
        color = .data$significant
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = .data$ci_low, xmax = .data$ci_high),
        height = 0.2
      ) +
      ggplot2::geom_point(size = point_size + 1) +
      ggplot2::scale_color_manual(
        values = colors,
        labels = c("FALSE" = "Not significant", "TRUE" = "Significant"),
        name = paste0(if (use_qval) "q" else "p", "-value < ", signif_threshold)
      ) +
      ggplot2::labs(
        x = "Coefficient (log2 fold change)",
        y = NULL,
        title = paste0(
          "MaAsLin3 Forest Plot - Top ", nrow(df), " features",
          if (!is.null(metadata_filter)) paste0(" (", metadata_filter, ")") else ""
        ),
        subtitle = if (show_annotations && !is.null(comparison_label)) {
          paste0("Positive = more in ", comparison_label, " | Negative = more in ", ref_label)
        } else {
          NULL
        }
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")
  }

  return(p)
}
