#' @title S7 class for comparing phyloseq objects
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' A class to store and compare multiple phyloseq objects. It contains:
#' - A list of phyloseq objects
#' - A summary table with computed characteristics for each phyloseq
#' - A list of comparison characteristics between phyloseq objects
#'
#' @section Types of Comparison:
#' The `list_phyloseq` class determines the type of comparison based on:
#' - Detected characteristics (sample overlap, nested samples, shared modalities)
#' - User-provided parameters (`same_primer_seq_tech` and `same_bioinfo_pipeline`)
#'
#' There are six main types of comparisons:
#'
#' \describe{
#'   \item{REPRODUCIBILITY}{Same pipeline (`same_bioinfo_pipeline = TRUE`),
#'     same primer/technology (`same_primer_seq_tech = TRUE`), same samples.
#'     Used to test **reproducibility** of results when running the exact same
#'     analysis multiple times.}
#'
#'   \item{ROBUSTNESS}{Different pipeline (`same_bioinfo_pipeline = FALSE`,
#'     e.g., different clustering method, different assignment database) but
#'     same primer/technology (`same_primer_seq_tech = TRUE`), same samples.
#'     Used to test **robustness** of conclusions to methodological choices.}
#'
#'   \item{NESTED_ROBUSTNESS}{One phyloseq object is derived from another, with
#'     samples being a subset (e.g., rarefied version created with
#'     `rarefy_even_depth()`). Used to test **robustness** to data processing
#'     choices like rarefaction, filtering, or subsetting. Comparisons should
#'     focus on the common (nested) samples.}
#'
#'   \item{REPLICABILITY}{Different primer and/or technology
#'     (`same_primer_seq_tech = FALSE`, e.g., ITS1 vs ITS2, Illumina vs PacBio),
#'     same samples. Used to test **replicability** across taxonomic groups or
#'     sequencing technologies.}
#'
#'   \item{EXPLORATION}{Different samples but with shared modalities.
#'     Useful to **explore** differences among groups of samples.
#'     Note: consider merging samples into one phyloseq object for
#'     some analyses instead.}
#'
#'   \item{SEPARATE_ANALYSIS}{Different samples with no shared modalities.
#'     **Separate analysis** of each phyloseq object is recommended as direct
#'     comparison may not be meaningful.}
#' }
#'
#' @name list_phyloseq
#' @importFrom S7 new_class new_property class_list class_any
#' @importFrom tibble tibble
#' @importFrom phyloseq nsamples ntaxa sample_sums taxa_sums sample_data
#' @export
NULL

# Helper functions for phyloseq summary ----------------------------------------

#' Compute summary statistics for a single phyloseq object
#' @param physeq (required) A phyloseq object.
#' @param name (character, default NULL) Optional name for the phyloseq object.
#' @param compute_kmer_dist (logical, default FALSE) If TRUE, compute mean k-mer
#'   distance.
#' @param kmer_dist (integer, default 5) k value for k-mer distance calculation.
#' @param verbose (logical, default TRUE) If TRUE, print progress messages.
#' @return A tibble row with computed characteristics
#' @noRd
compute_phyloseq_summary <- function(
  physeq,
  name = NULL,
  compute_kmer_dist = FALSE,
  kmer_dist = 5,
  verbose = TRUE
) {
  tib <- tibble::tibble(
    name = name %||% "unnamed",
    n_samples = phyloseq::nsamples(physeq),
    n_taxa = phyloseq::ntaxa(physeq),
    n_sequences = sum(phyloseq::sample_sums(physeq)),
    n_occurence = sum(phyloseq::otu_table(physeq) > 0),
    mean_seq_length = if (
      !is.null(phyloseq::refseq(physeq, errorIfNULL = FALSE))
    ) {
      mean(Biostrings::width(phyloseq::refseq(physeq)))
    } else {
      NA_real_
    },
    min_seq_length = if (
      !is.null(phyloseq::refseq(physeq, errorIfNULL = FALSE))
    ) {
      min(Biostrings::width(phyloseq::refseq(physeq)))
    } else {
      NA_real_
    },
    mean_seq_per_sample = mean(phyloseq::sample_sums(physeq)),
    sd_seq_per_sample = stats::sd(phyloseq::sample_sums(physeq)),
    min_seq_per_sample = min(phyloseq::sample_sums(physeq)),
    max_seq_per_sample = max(phyloseq::sample_sums(physeq)),
    mean_seq_per_taxon = mean(phyloseq::taxa_sums(physeq)),
    sd_seq_per_taxon = stats::sd(phyloseq::taxa_sums(physeq)),
    has_sam_data = !is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE)),
    has_tax_table = !is.null(phyloseq::tax_table(physeq, errorIfNULL = FALSE)),
    has_refseq = !is.null(phyloseq::refseq(physeq, errorIfNULL = FALSE)),
    has_phy_tree = !is.null(phyloseq::phy_tree(physeq, errorIfNULL = FALSE))
  )

  if (compute_kmer_dist) {
    if (verbose) {
      cli::cli_alert_info(
        "Computing mean {kmer_dist}-mer distance for phyloseq object '{name}'"
      )
    }
    tib$mean_kmer_dist <- if (
      !is.null(phyloseq::refseq(physeq, errorIfNULL = FALSE))
    ) {
      seqs <- ape::as.DNAbin(phyloseq::refseq(physeq))
      mean(kmer::kdistance(seqs, k = kmer_dist))
    } else {
      NA_real_
    }
  }
  return(tib)
}

#' Build summary table for a list of phyloseq objects
#' @param physeq_list A named list of phyloseq objects
#' @return A tibble with one row per phyloseq object
#' @noRd
build_summary_table <- function(physeq_list) {
  names_list <- names(physeq_list)
  if (is.null(names_list)) {
    names_list <- paste0("physeq_", seq_along(physeq_list))
  }
  purrr::map2_dfr(physeq_list, names_list, compute_phyloseq_summary)
}

# Helper functions for comparison characteristics ------------------------------

#' Check if all phyloseq objects share the same sample_data structure
#' @param physeq_list A list of phyloseq objects
#' @return Logical
#' @noRd
check_same_sam_data_structure <- function(physeq_list) {
  sam_data_list <- purrr::map(
    physeq_list,
    ~ {
      sam <- phyloseq::sample_data(.x, errorIfNULL = FALSE)
      if (is.null(sam)) {
        return(NULL)
      }
      colnames(sam)
    }
  )
  sam_data_list <- purrr::compact(sam_data_list)
  if (length(sam_data_list) < 2) {
    return(NA)
  }
  reference <- sam_data_list[[1]]
  all(purrr::map_lgl(sam_data_list[-1], ~ identical(sort(.x), sort(reference))))
}

#' Check if all phyloseq objects share the same samples
#' @param physeq_list A list of phyloseq objects
#' @return Logical
#' @noRd
check_same_samples <- function(physeq_list) {
  sample_names_list <- purrr::map(physeq_list, phyloseq::sample_names)
  reference <- sample_names_list[[1]]
  all(purrr::map_lgl(
    sample_names_list[-1],
    ~ identical(sort(.x), sort(reference))
  ))
}

#' Check if samples are nested (one phyloseq's samples are subset of another's)
#'
#' @description
#' Detects if phyloseq objects have nested sample relationships, which is common
#' when one object is derived from another (e.g., rarefied version of original).
#'
#' @param physeq_list A list of phyloseq objects
#' @return A list with:
#'   - `is_nested`: Logical, TRUE if samples are nested
#'   - `nesting_structure`: Character describing the nesting relationship
#' @noRd
check_nested_samples <- function(physeq_list) {
  sample_names_list <- purrr::map(physeq_list, phyloseq::sample_names)
  names_list <- names(physeq_list)
  n <- length(sample_names_list)

  if (n < 2) {
    return(list(is_nested = FALSE, nesting_structure = NULL))
  }

  # Check all pairwise combinations for subset relationships
  nested_pairs <- list()

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) {
        samples_i <- sample_names_list[[i]]
        samples_j <- sample_names_list[[j]]

        # Check if i is a proper subset of j (i contained in j, but not equal)
        if (
          all(samples_i %in% samples_j) && length(samples_i) < length(samples_j)
        ) {
          nested_pairs <- c(
            nested_pairs,
            list(list(
              subset = names_list[i],
              superset = names_list[j],
              n_subset = length(samples_i),
              n_superset = length(samples_j)
            ))
          )
        }
      }
    }
  }

  if (length(nested_pairs) > 0) {
    # Build description of nesting structure
    descriptions <- purrr::map_chr(nested_pairs, function(pair) {
      paste0(
        pair$subset,
        " (",
        pair$n_subset,
        " samples) is nested in ",
        pair$superset,
        " (",
        pair$n_superset,
        " samples)"
      )
    })
    nesting_structure <- paste(descriptions, collapse = "; ")

    return(list(is_nested = TRUE, nesting_structure = nesting_structure))
  }

  return(list(is_nested = FALSE, nesting_structure = NULL))
}

#' Check if all phyloseq objects share the same taxa
#' @param physeq_list A list of phyloseq objects
#' @return Logical
#' @noRd
check_same_taxa <- function(physeq_list) {
  taxa_names_list <- purrr::map(physeq_list, phyloseq::taxa_names)
  reference <- taxa_names_list[[1]]
  all(purrr::map_lgl(
    taxa_names_list[-1],
    ~ identical(sort(.x), sort(reference))
  ))
}

#' Get shared modalities across sample_data variables
#' @param physeq_list A list of phyloseq objects
#' @return A named list with shared modalities per variable
#' @noRd
get_shared_sam_data_modalities <- function(physeq_list) {
  sam_data_list <- purrr::map(
    physeq_list,
    ~ {
      sam <- phyloseq::sample_data(.x, errorIfNULL = FALSE)
      if (is.null(sam)) {
        return(NULL)
      }
      as.data.frame(sam)
    }
  )
  sam_data_list <- purrr::compact(sam_data_list)
  if (length(sam_data_list) < 2) {
    return(list())
  }

  common_cols <- Reduce(intersect, purrr::map(sam_data_list, colnames))
  if (length(common_cols) == 0) {
    return(list())
  }

  result <- purrr::map(common_cols, function(col) {
    values_list <- purrr::map(
      sam_data_list,
      ~ {
        val <- .x[[col]]
        if (is.factor(val) || is.character(val)) {
          unique(as.character(val))
        } else {
          NULL
        }
      }
    )
    values_list <- purrr::compact(values_list)
    if (length(values_list) < 2) {
      return(NULL)
    }
    Reduce(intersect, values_list)
  })
  names(result) <- common_cols
  purrr::compact(result)
}

#' Get common samples across all phyloseq objects
#' @param physeq_list A list of phyloseq objects
#' @return Character vector of shared sample names
#' @noRd
get_common_samples <- function(physeq_list) {
  sample_names_list <- purrr::map(physeq_list, phyloseq::sample_names)
  Reduce(intersect, sample_names_list)
}

#' Get common taxa across all phyloseq objects
#' @param physeq_list A list of phyloseq objects
#' @return Character vector of shared taxa names
#' @noRd
get_common_taxa <- function(physeq_list) {
  taxa_names_list <- purrr::map(physeq_list, phyloseq::taxa_names)
  Reduce(intersect, taxa_names_list)
}

#' Determine the type of comparison based on detected characteristics
#'
#' @description
#' Determines the comparison type based on sample overlap, shared modalities,
#' and user-provided information about primer/sequencing technology and
#' bioinformatics pipeline.
#'
#' The six comparison types are:
#'
#' 1. **REPRODUCIBILITY**: Same pipeline, same primer/technology, same samples.
#'    Used to test reproducibility of results.
#'
#' 2. **ROBUSTNESS**: Different pipeline (e.g., different clustering method,
#'    different assignment database) but same primer/technology, same samples.
#'    Used to test robustness to methodological choices.
#'
#' 3. **NESTED_ROBUSTNESS**: One phyloseq object is derived from another (e.g.,
#'    rarefied version). Samples from one are a subset of samples from another.
#'    Used to test robustness to data processing choices like rarefaction.
#'
#' 4. **REPLICABILITY**: Different primer (ITS1, ITS2) and/or technology
#'    (Illumina, PacBio), same samples.
#'    Used to test replicability across taxonomic groups or sequencing technologies.
#'
#' 5. **EXPLORATION**: Different samples but with shared modalities.
#'    Useful to explore differences among groups of samples.
#'    Consider merging samples into one phyloseq object instead.
#'
#' 6. **SEPARATE_ANALYSIS**: Different samples with no shared modalities.
#'    Separate analysis of each phyloseq object is recommended.
#'
#' @param same_samples Logical, whether all phyloseq objects have the same samples
#' @param nested_info List from check_nested_samples() with is_nested and nesting_structure
#' @param has_shared_modalities Logical, whether there are shared sample_data modalities
#' @param same_primer_seq_tech Logical, whether the same primer and sequencing
#'   technology was used across all phyloseq objects
#' @param same_bioinfo_pipeline Logical, whether the same bioinformatics pipeline
#'   was used across all phyloseq objects
#' @return A list with `type` (character) and `description` (character)
#' @noRd
determine_comparison_type <- function(
  same_samples,
  nested_info,
  has_shared_modalities,
  same_primer_seq_tech = TRUE,
  same_bioinfo_pipeline = TRUE
) {
  if (same_samples) {
    # Same samples: differentiate based on pipeline and primer/technology
    if (!same_primer_seq_tech) {
      # Different primer/technology: REPLICABILITY
      list(
        type = "REPLICABILITY",
        description = paste0(
          "Same samples with different primer/sequencing technology.\n",
          "  Used to test replicability across taxonomic groups or\n",
          "  sequencing technologies (e.g., ITS1 vs ITS2, Illumina vs PacBio)."
        )
      )
    } else if (!same_bioinfo_pipeline) {
      # Same primer/technology but different pipeline: ROBUSTNESS
      list(
        type = "ROBUSTNESS",
        description = paste0(
          "Same samples, same primer/technology, different bioinformatics pipeline.\n",
          "  Used to test robustness to methodological choices\n",
          "  (e.g., different clustering, different taxonomic database)."
        )
      )
    } else {
      # Same pipeline, same primer/technology: REPRODUCIBILITY
      list(
        type = "REPRODUCIBILITY",
        description = paste0(
          "Same samples, same primer/technology, same bioinformatics pipeline.\n",
          "  Used to test reproducibility of results when running\n",
          "  the exact same analysis multiple times."
        )
      )
    }
  } else if (nested_info$is_nested) {
    # Nested samples: one phyloseq is derived from another (e.g., rarefied)
    list(
      type = "NESTED_ROBUSTNESS",
      description = paste0(
        "Nested samples detected (one phyloseq derived from another).\n",
        "  Nesting: ",
        nested_info$nesting_structure,
        "\n",
        "  Useful to test robustness to data processing (e.g., rarefaction).\n",
        "  Comparisons should focus on the common (nested) samples."
      )
    )
  } else if (has_shared_modalities) {
    # Different samples but shared modalities: EXPLORATION
    list(
      type = "EXPLORATION",
      description = paste0(
        "Different samples but shared modalities detected.\n",
        "  Useful to explore differences among groups of samples.\n",
        "  Consider merging samples into one phyloseq object for some analyses."
      )
    )
  } else {
    # Different samples, no shared modalities: SEPARATE_ANALYSIS
    list(
      type = "SEPARATE_ANALYSIS",
      description = paste0(
        "Different samples with no shared modalities detected.\n",
        "  Separate analysis of each phyloseq object is recommended.\n",
        "  Direct comparison may not be meaningful."
      )
    )
  }
}

#' Build comparison characteristics for a list of phyloseq objects
#' @param physeq_list A named list of phyloseq objects
#' @param same_primer_seq_tech Logical, whether the same primer and sequencing
#'   technology was used across all phyloseq objects
#' @param same_bioinfo_pipeline Logical, whether the same bioinformatics pipeline
#'   was used across all phyloseq objects
#' @param compute_dist Logical, whether to compute k-mer distances between
#'   refseq slots
#' @return A list of comparison characteristics
#' @noRd
build_comparison_characteristics <- function(
  physeq_list,
  same_primer_seq_tech = TRUE,
  same_bioinfo_pipeline = TRUE,
  compute_dist = TRUE,
  verbose = TRUE
) {
  if (verbose) {
    cli::cli_alert_info("Checking sample and taxa overlap...")
  }
  common_samples <- get_common_samples(physeq_list)
  common_taxa <- get_common_taxa(physeq_list)
  same_samples <- check_same_samples(physeq_list)
  nested_info <- check_nested_samples(physeq_list)
  shared_modalities <- get_shared_sam_data_modalities(physeq_list)
  has_shared_modalities <- length(shared_modalities) > 0

  # Determine comparison type
  comparison_type <- determine_comparison_type(
    same_samples,
    nested_info,
    has_shared_modalities,
    same_primer_seq_tech,
    same_bioinfo_pipeline
  )

  if (verbose) {
    cli::cli_alert_info(
      "Detected comparison type: {comparison_type$type}"
    )
    cli::cli_alert_info(
      "{length(common_samples)} common sample{?s}, {length(common_taxa)} common tax{?on/a}"
    )
  }

  list(
    type_of_comparison = comparison_type$type,
    type_description = comparison_type$description,
    same_primer_seq_tech = same_primer_seq_tech,
    same_bioinfo_pipeline = same_bioinfo_pipeline,
    same_sam_data_structure = check_same_sam_data_structure(physeq_list),
    same_samples = same_samples,
    nested_samples = nested_info$is_nested,
    nesting_structure = nested_info$nesting_structure,
    same_taxa = check_same_taxa(physeq_list),
    n_common_samples = length(common_samples),
    common_samples = common_samples,
    n_common_taxa = length(common_taxa),
    common_taxa = common_taxa,
    shared_sam_data_modalities = shared_modalities,
    all_have_sam_data = all(purrr::map_lgl(
      physeq_list,
      ~ !is.null(phyloseq::sample_data(.x, errorIfNULL = FALSE))
    )),
    all_have_tax_table = all(purrr::map_lgl(
      physeq_list,
      ~ !is.null(phyloseq::tax_table(.x, errorIfNULL = FALSE))
    )),
    all_have_refseq = all(purrr::map_lgl(
      physeq_list,
      ~ !is.null(phyloseq::refseq(.x, errorIfNULL = FALSE))
    )),
    all_have_phy_tree = all(purrr::map_lgl(
      physeq_list,
      ~ !is.null(phyloseq::phy_tree(.x, errorIfNULL = FALSE))
    )),
    refseq_comparison = build_refseq_comparison(
      physeq_list,
      compute_dist,
      verbose
    )
  )
}

#' Build refseq comparison for all pairs of phyloseq objects
#'
#' Returns NULL if not all objects have a refseq slot, or a named list of
#' `compare_refseq` results (one per pair).
#' @param physeq_list A named list of phyloseq objects
#' @param compute_dist Logical, whether to compute k-mer distances
#' @param verbose Logical, whether to print progress messages
#' @return A named list of compare_refseq results, or NULL
#' @noRd
build_refseq_comparison <- function(
  physeq_list,
  compute_dist = TRUE,
  verbose = TRUE
) {
  all_have_refseq <- all(purrr::map_lgl(
    physeq_list,
    ~ !is.null(phyloseq::refseq(.x, errorIfNULL = FALSE))
  ))
  if (!all_have_refseq || length(physeq_list) < 2) {
    if (verbose && !all_have_refseq) {
      cli::cli_alert_info(
        "Skipping refseq comparison (not all objects have refseq)"
      )
    }
    return(NULL)
  }

  if (!compute_dist) {
    if (verbose) {
      cli::cli_alert_info(
        "Skipping refseq distance computation (compute_dist = FALSE)"
      )
    }
    return(NULL)
  }

  nms <- names(physeq_list)
  pairs <- utils::combn(seq_along(physeq_list), 2, simplify = FALSE)

  results <- purrr::map(pairs, function(idx) {
    compare_refseq(
      physeq_list[[idx[1]]],
      physeq_list[[idx[2]]],
      name1 = nms[idx[1]],
      name2 = nms[idx[2]],
      verbose = FALSE
    )
  })
  names(results) <- purrr::map_chr(
    pairs,
    ~ paste0(nms[.x[1]], "_vs_", nms[.x[2]])
  )
  results
}

# S7 Class Definition ----------------------------------------------------------

#' list_phyloseq S7 class
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' An S7 class to store and compare multiple phyloseq objects.
#'
#' @param physeq_list (required) A named list of phyloseq objects.
#' @param names (character vector, default NULL) Optional names for the phyloseq
#'   objects. If NULL and the list is unnamed, names are generated automatically.
#' @param same_primer_seq_tech (logical, default TRUE) Whether the same primer
#'   and sequencing technology was used across all phyloseq objects. Set to FALSE
#'   when comparing different primers (e.g., ITS1 vs ITS2) or technologies
#'   (e.g., Illumina vs PacBio).
#' @param same_bioinfo_pipeline (logical, default TRUE) Whether the same
#'   bioinformatics pipeline was used across all phyloseq objects. Set to FALSE
#'   when comparing different clustering methods, taxonomic databases, or
#'   analysis parameters.
#' @param compute_dist (logical, default TRUE) Whether to compute pairwise
#'   k-mer distances between refseq slots. Set to FALSE to skip this
#'   potentially slow computation, especially when phyloseq objects have many
#'   taxa or very different sequence sets.
#' @param verbose (logical, default TRUE) Whether to print progress messages
#'   via [cli::cli_alert_info()].
#'
#' @slot phyloseq_list A named list of phyloseq objects
#' @slot summary_table A tibble summarizing each phyloseq object
#' @slot comparison A list of characteristics comparing the phyloseq objects
#'
#' @details
#' The combination of `same_primer_seq_tech` and `same_bioinfo_pipeline` along
#' with detected sample overlap determines the comparison type:
#'
#' | same_samples | same_primer_seq_tech | same_bioinfo_pipeline | Type |
#' |--------------|----------------------|-----------------------|------|
#' | TRUE | TRUE | TRUE | REPRODUCIBILITY |
#' | TRUE | TRUE | FALSE | ROBUSTNESS |
#' | TRUE | FALSE | - | REPLICABILITY |
#' | nested | - | - | NESTED_ROBUSTNESS |
#' | FALSE | - | - | EXPLORATION or SEPARATE_ANALYSIS |
#'
#' @examples
#' # REPRODUCIBILITY: Same samples, same pipeline (default)
#' lpq_repro <- list_phyloseq(list(run1 = data_fungi, run2 = data_fungi))
#'
#' # ROBUSTNESS: Same samples, different pipeline
#' lpq_robust <- list_phyloseq(
#'   list(method_A = data_fungi, method_B = data_fungi),
#'   same_bioinfo_pipeline = FALSE
#' )
#'
#' # REPLICABILITY: Same samples, different primer/technology
#' lpq_replic <- list_phyloseq(
#'   list(ITS1 = data_fungi, ITS2 = data_fungi),
#'   same_primer_seq_tech = FALSE
#' )
#'
#' @export
list_phyloseq <- S7::new_class(
  name = "list_phyloseq",
  properties = list(
    phyloseq_list = S7::new_property(
      class = S7::class_list,
      validator = function(value) {
        if (length(value) == 0) {
          return("phyloseq_list must contain at least one phyloseq object")
        }
        are_phyloseq <- purrr::map_lgl(value, ~ inherits(.x, "phyloseq"))
        if (!all(are_phyloseq)) {
          return("All elements in phyloseq_list must be phyloseq objects")
        }
        NULL
      }
    ),
    summary_table = S7::new_property(
      class = S7::class_any,
      validator = function(value) {
        if (!inherits(value, "tbl_df")) {
          return("summary_table must be a tibble")
        }
        NULL
      }
    ),
    comparison = S7::new_property(
      class = S7::class_list
    )
  ),
  constructor = function(
    physeq_list,
    names = NULL,
    same_primer_seq_tech = TRUE,
    same_bioinfo_pipeline = TRUE,
    compute_dist = TRUE,
    verbose = TRUE
  ) {
    if (!is.null(names)) {
      if (length(names) != length(physeq_list)) {
        stop("Length of 'names' must match length of 'physeq_list'")
      }
      names(physeq_list) <- names
    } else if (is.null(names(physeq_list))) {
      names(physeq_list) <- paste0("physeq_", seq_along(physeq_list))
    }

    if (verbose) {
      cli::cli_alert_info(
        "Building summary table for {length(physeq_list)} phyloseq object{?s}..."
      )
    }
    summary_table <- build_summary_table(physeq_list)

    if (verbose) {
      cli::cli_alert_info("Computing comparison characteristics...")
    }
    comparison <- build_comparison_characteristics(
      physeq_list,
      same_primer_seq_tech,
      same_bioinfo_pipeline,
      compute_dist,
      verbose
    )

    if (verbose) {
      cli::cli_alert_success(
        "list_phyloseq created ({comparison$type_of_comparison})"
      )
    }

    S7::new_object(
      S7::S7_object(),
      phyloseq_list = physeq_list,
      summary_table = summary_table,
      comparison = comparison
    )
  }
)

# Methods ----------------------------------------------------------------------

# Print method for list_phyloseq
S7::method(print, list_phyloseq) <- function(x, ...) {
  cat(
    "list_phyloseq object with",
    length(x@phyloseq_list),
    "phyloseq objects\n"
  )
  cat("\n--- Summary ---\n")
  print(x@summary_table)
  cat("\n--- Comparison characteristics ---\n")
  cat("Type of comparison:", x@comparison$type_of_comparison, "\n")
  cat(x@comparison$type_description, "\n\n")

  cat("Same primer/seq tech:", x@comparison$same_primer_seq_tech, "\n")
  cat("Same bioinfo pipeline:", x@comparison$same_bioinfo_pipeline, "\n")
  cat("Same sample_data structure:", x@comparison$same_sam_data_structure, "\n")
  cat("Same samples:", x@comparison$same_samples, "\n")
  cat("Nested samples:", x@comparison$nested_samples, "\n")
  cat("Same taxa:", x@comparison$same_taxa, "\n")
  cat("Common samples:", x@comparison$n_common_samples, "\n")
  cat("Common taxa:", x@comparison$n_common_taxa, "\n")

  if (!is.null(x@comparison$refseq_comparison)) {
    cat("\n--- Reference sequence comparison ---\n")
    for (pair_name in names(x@comparison$refseq_comparison)) {
      rc <- x@comparison$refseq_comparison[[pair_name]]
      cat(
        pair_name,
        ": ",
        length(rc$shared_seqs),
        " shared seqs, ",
        length(rc$unique_seqs_1),
        " unique in ",
        rc$name1,
        ", ",
        length(rc$unique_seqs_2),
        " unique in ",
        rc$name2,
        "\n",
        sep = ""
      )
    }
  }
  invisible(x)
}

# Length method for list_phyloseq
S7::method(length, list_phyloseq) <- function(x) {
  length(x@phyloseq_list)
}

# Names method for list_phyloseq
S7::method(names, list_phyloseq) <- function(x) {
  names(x@phyloseq_list)
}

# Subset method for list_phyloseq
S7::method(`[`, list_phyloseq) <- function(x, i) {
  list_phyloseq(
    x@phyloseq_list[i],
    same_primer_seq_tech = x@comparison$same_primer_seq_tech,
    same_bioinfo_pipeline = x@comparison$same_bioinfo_pipeline,
    verbose = FALSE
  )
}

# Extract method for list_phyloseq (single element)
S7::method(`[[`, list_phyloseq) <- function(x, i) {
  x@phyloseq_list[[i]]
}

# Utility functions ------------------------------------------------------------

#' Update the summary table and comparison characteristics
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Recomputes the summary_table and comparison slots. Useful after
#' modifying the phyloseq objects in place or to change the comparison
#' parameters.
#'
#' @param x (required) A list_phyloseq object.
#' @param same_primer_seq_tech (logical, default NULL) Whether the same primer
#'   and sequencing technology was used. If NULL, preserves the original value.
#' @param same_bioinfo_pipeline (logical, default NULL) Whether the same
#'   bioinformatics pipeline was used. If NULL, preserves the original value.
#' @param compute_dist (logical, default TRUE) Whether to compute pairwise
#'   k-mer distances between refseq slots.
#' @return An updated list_phyloseq object
#' @export
update_list_phyloseq <- function(
  x,
  same_primer_seq_tech = NULL,
  same_bioinfo_pipeline = NULL,
  compute_dist = TRUE,
  verbose = TRUE
) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  # Preserve original values if not specified
  if (is.null(same_primer_seq_tech)) {
    same_primer_seq_tech <- x@comparison$same_primer_seq_tech
  }
  if (is.null(same_bioinfo_pipeline)) {
    same_bioinfo_pipeline <- x@comparison$same_bioinfo_pipeline
  }

  list_phyloseq(
    x@phyloseq_list,
    same_primer_seq_tech = same_primer_seq_tech,
    same_bioinfo_pipeline = same_bioinfo_pipeline,
    compute_dist = compute_dist,
    verbose = verbose
  )
}

#' Add a phyloseq object to a list_phyloseq
#'
#' @param x (required) A list_phyloseq object.
#' @param physeq (required) A phyloseq object to add.
#' @param name (character, default NULL) Optional name for the new phyloseq object.
#'   If NULL, a name is generated automatically.
#' @return A new list_phyloseq object with the added phyloseq
#' @export
add_phyloseq <- function(x, physeq, name = NULL, verbose = TRUE) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))
  stopifnot(inherits(physeq, "phyloseq"))

  new_list <- x@phyloseq_list
  idx <- length(new_list) + 1

  if (is.null(name)) {
    name <- paste0("physeq_", idx)
  }

  new_list[[name]] <- physeq
  list_phyloseq(
    new_list,
    same_primer_seq_tech = x@comparison$same_primer_seq_tech,
    same_bioinfo_pipeline = x@comparison$same_bioinfo_pipeline,
    verbose = verbose
  )
}

#' Remove a phyloseq object from a list_phyloseq
#'
#' @param x (required) A list_phyloseq object.
#' @param name (character or integer) Name or index of the phyloseq object to remove.
#' @return A new list_phyloseq object without the removed phyloseq
#' @export
remove_phyloseq <- function(x, name, verbose = TRUE) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  new_list <- x@phyloseq_list
  if (is.character(name)) {
    new_list[[name]] <- NULL
  } else if (is.numeric(name)) {
    new_list <- new_list[-name]
  }

  if (length(new_list) == 0) {
    stop("Cannot remove the last phyloseq object from list_phyloseq")
  }

  list_phyloseq(
    new_list,
    same_primer_seq_tech = x@comparison$same_primer_seq_tech,
    same_bioinfo_pipeline = x@comparison$same_bioinfo_pipeline,
    verbose = verbose
  )
}

#' Filter phyloseq objects to keep only shared samples and/or taxa
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Filters each phyloseq object in a list_phyloseq to retain only the samples
#' and/or taxa that are common across all objects. This is useful for making
#' direct comparisons on a common basis.
#'
#' @param x (required) A list_phyloseq object.
#' @param filter_samples (logical, default TRUE) If TRUE, filter to keep only
#'   samples present in all phyloseq objects.
#' @param filter_taxa (logical, default FALSE) If TRUE, filter to keep only
#'   taxa present in all phyloseq objects.
#' @param clean (logical, default TRUE) If TRUE, apply [MiscMetabar::clean_pq()]
#'   after filtering to remove empty samples/taxa.
#' @param verbose (logical, default TRUE) If TRUE, print information about the
#'   filtering process.
#'
#' @return A new list_phyloseq object with filtered phyloseq objects
#'
#' @details
#' This function is particularly useful for:
#' - **NESTED_ROBUSTNESS** comparisons: filter to common samples when comparing
#'   original vs rarefied data
#' - **EXPLORATION** comparisons: filter to common samples/taxa when comparing
#'   different sample groups
#' - Any comparison where you need to ensure all phyloseq objects contain the
#'   same samples or taxa
#'
#' @examples
#' lpq <- list_phyloseq(list(
#'   original = data_fungi,
#'   mini = data_fungi_mini,
#'   rarefied = rarefy_even_depth(data_fungi)
#' ))
#'
#' # Filter to keep only common samples (useful for nested comparisons)
#' lpq_filtered <- filter_common_lpq(lpq, filter_samples = TRUE, verbose = FALSE)
#'
#' # Filter to keep only common taxa
#' lpq_filtered <- filter_common_lpq(lpq, filter_samples = FALSE, filter_taxa = TRUE)
#'
#' # Filter both samples and taxa
#' lpq_filtered <- filter_common_lpq(lpq, filter_samples = TRUE, filter_taxa = TRUE)
#'
#' @export
filter_common_lpq <- function(
  x,
  filter_samples = TRUE,
  filter_taxa = FALSE,
  clean = TRUE,
  verbose = TRUE
) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  if (!filter_samples && !filter_taxa) {
    if (verbose) {
      message("No filtering requested. Returning original list_phyloseq.")
    }
    return(x)
  }
  common_samples <- x@comparison$common_samples
  common_taxa <- x@comparison$common_taxa

  if (filter_samples && length(common_samples) == 0) {
    stop("No common samples found across phyloseq objects. Cannot filter.")
  }
  if (filter_taxa && length(common_taxa) == 0) {
    stop("No common taxa found across phyloseq objects. Cannot filter.")
  }

  if (verbose) {
    if (filter_samples) {
      message("Filtering to ", length(common_samples), " common samples.")
    }
    if (filter_taxa) {
      message("Filtering to ", length(common_taxa), " common taxa.")
    }
  }

  filtered_list <- purrr::imap(x@phyloseq_list, function(pq, name) {
    original_nsamples <- phyloseq::nsamples(pq)
    original_ntaxa <- phyloseq::ntaxa(pq)

    if (filter_samples) {
      samples_to_keep <- phyloseq::sample_names(pq) %in% common_samples
      pq <- phyloseq::prune_samples(samples_to_keep, pq)
    }

    if (filter_taxa) {
      taxa_to_keep <- phyloseq::taxa_names(pq) %in% common_taxa
      pq <- phyloseq::prune_taxa(taxa_to_keep, pq)
    }

    if (clean) {
      pq <- MiscMetabar::clean_pq(pq)
    }

    if (verbose) {
      message(
        "  ",
        name,
        ": ",
        original_nsamples,
        " -> ",
        phyloseq::nsamples(pq),
        " samples, ",
        original_ntaxa,
        " -> ",
        phyloseq::ntaxa(pq),
        " taxa"
      )
    }
    pq
  })

  list_phyloseq(
    filtered_list,
    same_primer_seq_tech = x@comparison$same_primer_seq_tech,
    same_bioinfo_pipeline = x@comparison$same_bioinfo_pipeline,
    verbose = verbose
  )
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


#' Count unique taxonomic levels across phyloseq objects
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Creates a summary table showing the number of unique taxonomic values
#' (levels) for each taxonomic rank across all phyloseq objects in a
#' list_phyloseq. This is useful for comparing taxonomic resolution and
#' diversity across different datasets or classification methods.
#'
#' @param x (required) A list_phyloseq object.
#' @param taxonomic_ranks (character vector, required) Names of taxonomic ranks
#'   to count. Must be present in the tax_table of ALL phyloseq objects in the
#'   list.
#' @param na.rm (logical, default TRUE) If TRUE, NA values are excluded when
#'   counting unique levels.
#'
#' @return A data frame with:
#' \describe{
#'   \item{Rows}{One row per phyloseq object (named by the phyloseq name)}
#'   \item{Columns}{One column per taxonomic rank, containing the count of
#'     unique values for that rank in that phyloseq object}
#' }
#'
#' @export
#' @author Adrien Taudière
#'
#' @seealso [upset_lpq()]
#' @examples
#' lpq <- list_phyloseq(list(fungi = data_fungi, fungi_mini = data_fungi_mini))
#'
#' n_levels_lpq(lpq, c("Phylum", "Class", "Order", "Family", "Genus"))
n_levels_lpq <- function(x, taxonomic_ranks, na.rm = TRUE) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  pq_names <- names(x@phyloseq_list)

  missing_ranks_per_pq <- purrr::imap(x@phyloseq_list, function(pq, name) {
    available_ranks <- colnames(pq@tax_table)
    missing <- setdiff(taxonomic_ranks, available_ranks)
    if (length(missing) > 0) {
      return(list(name = name, missing = missing))
    }
    NULL
  })

  missing_ranks_per_pq <- purrr::compact(missing_ranks_per_pq)

  if (length(missing_ranks_per_pq) > 0) {
    msg <- purrr::map_chr(missing_ranks_per_pq, function(m) {
      paste0("  - ", m$name, ": missing ", paste(m$missing, collapse = ", "))
    })
    stop(
      "Some taxonomic_ranks are missing in the following phyloseq objects:\n",
      paste(msg, collapse = "\n")
    )
  }

  count_levels <- function(pq, ranks, na.rm) {
    tax_df <- as.data.frame(pq@tax_table[, ranks, drop = FALSE])
    purrr::map_int(tax_df, function(col) {
      if (na.rm) {
        length(unique(col[!is.na(col) & col != "" & col != "NA_NA"]))
      } else {
        length(unique(col))
      }
    })
  }

  counts_list <- purrr::map(
    x@phyloseq_list,
    count_levels,
    ranks = taxonomic_ranks,
    na.rm = na.rm
  )

  result <- do.call(rbind, counts_list)
  result <- as.data.frame(result)
  rownames(result) <- pq_names

  result
}


#' Apply a function to all phyloseq objects in a list_phyloseq
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Applies a function to each phyloseq object in a list_phyloseq and returns
#' a new list_phyloseq with the transformed objects. This is useful for
#' applying cleaning, filtering, or transformation functions uniformly across
#' all phyloseq objects.
#'
#' @param x (required) A list_phyloseq object.
#' @param .f (function, required) A function to apply to each phyloseq object.
#'   The function must take a phyloseq object as its first argument and return
#'   a phyloseq object.
#' @param ... Additional arguments passed to `.f`.
#' @param verbose (logical, default TRUE) If TRUE, print information about the
#'   transformation applied to each phyloseq object.
#'
#' @return A new list_phyloseq object with the transformed phyloseq objects.
#'   The comparison parameters (`same_primer_seq_tech`, `same_bioinfo_pipeline`)
#'   are preserved from the original object.
#'
#' @details
#' Common functions to apply include:
#' - `MiscMetabar::clean_pq()`: Remove empty samples and taxa
#' - `phyloseq::taxa_as_rows()`: Ensure taxa are rows in otu_table
#' - `phyloseq::rarefy_even_depth()`: Rarefy to even depth
#' - `phyloseq::transform_sample_counts()`: Transform counts (e.g., relative abundance)
#' - `phyloseq::subset_taxa()`: Filter taxa based on criteria
#' - `phyloseq::subset_samples()`: Filter samples based on criteria
#'
#' @export
#' @author Adrien Taudière
#'
#' @seealso [filter_common_lpq()], [update_list_phyloseq()]
#'
#' @examples
#' lpq <- list_phyloseq(list(fungi = data_fungi, fungi_mini = data_fungi_mini))
#'
#' # Apply clean_pq to all phyloseq objects
#' lpq_clean <- apply_to_lpq(lpq, MiscMetabar::clean_pq)
#'
#' # Apply taxa_as_rows
#' lpq_rows <- apply_to_lpq(lpq, MiscMetabar::taxa_as_rows)
#'
#' # Apply rarefy_even_depth with a specific rngseed
#' lpq_rar <- apply_to_lpq(lpq, rarefy_even_depth, rngseed = 21)
#'
#' lpq_rar
#'
#' # Transform to relative abundance
#' lpq_rel <- apply_to_lpq(
#'   lpq,
#'   phyloseq::transform_sample_counts,
#'   function(x) x / sum(x)
#' )
#'
#' # Chain multiple transformations
#' lpq_processed <- lpq |>
#'   apply_to_lpq(MiscMetabar::clean_pq) |>
#'   apply_to_lpq(MiscMetabar::taxa_as_rows)
apply_to_lpq <- function(x, .f, ..., verbose = TRUE) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  if (!is.function(.f)) {
    stop("`.f` must be a function")
  }

  func_name <- deparse(substitute(.f))

  transformed_list <- purrr::imap(x@phyloseq_list, function(pq, name) {
    if (verbose) {
      original_nsamples <- phyloseq::nsamples(pq)
      original_ntaxa <- phyloseq::ntaxa(pq)
    }

    result <- .f(pq, ...)

    if (!inherits(result, "phyloseq")) {
      stop(
        "Function `",
        func_name,
        "` did not return a phyloseq object for '",
        name,
        "'. Got class: ",
        paste(class(result), collapse = ", ")
      )
    }

    if (verbose) {
      new_nsamples <- phyloseq::nsamples(result)
      new_ntaxa <- phyloseq::ntaxa(result)
      message(
        "  ",
        name,
        ": ",
        original_nsamples,
        " -> ",
        new_nsamples,
        " samples, ",
        original_ntaxa,
        " -> ",
        new_ntaxa,
        " taxa"
      )
    }

    result
  })

  if (verbose) {
    message(
      "Applied `",
      func_name,
      "` to ",
      length(transformed_list),
      " phyloseq objects."
    )
  }

  list_phyloseq(
    transformed_list,
    same_primer_seq_tech = x@comparison$same_primer_seq_tech,
    same_bioinfo_pipeline = x@comparison$same_bioinfo_pipeline,
    verbose = verbose
  )
}
