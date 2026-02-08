#' Add fake sequences by shuffling existing ones in a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful to compute true negative values with functions
#'   [tc_metrics_mock()] and [tc_metrics_mock_vec()].
#'
#' Note the that the tax_table for additional sequences is full of NA and
#'   that the corresponding otu_table is full of 0.
#'
#' @inheritParams tc_points_matrix
#' @param n_fake (integer, default NULL) A number of fake taxa to add.
#'   Must specify either `n_fake` or `prop_fake`, not both.
#' @param prop_fake (numeric, default NULL) A proportion of fake taxa to add.
#'   Must specify either `n_fake` or `prop_fake`, not both.
#' @param prefix (character, default "fake_") A prefix to add to the taxa name.
#' @returns A phyloseq object
#'
#' @export
#' @author Adrien Taudière
#' @examples
#' d_fake_F <- data_fungi_mini |>
#'   add_shuffle_seq_pq(prop_fake = 0.1)
#' ntaxa(d_fake_F) - ntaxa(data_fungi_mini)
add_shuffle_seq_pq <- function(
  physeq,
  n_fake = NULL,
  prop_fake = NULL,
  prefix = "fake_"
) {
  taxasrow <- taxa_are_rows(physeq)
  if (taxasrow) {
    physeq <- MiscMetabar::taxa_as_columns(physeq)
  }
  if (is.null(n_fake) && is.null(prop_fake)) {
    stop("You must specify either n_fake or prop_fake param.")
  }
  if (!is.null(n_fake) && !is.null(prop_fake)) {
    stop("You must specify either n_fake or prop_fake param, not both !")
  }
  if (is.null(n_fake) && !is.null(prop_fake)) {
    n_fake <- prop_fake * ntaxa(physeq)
  }
  taxa_to_subset <- taxa_names(physeq) %in% sample(taxa_names(physeq), n_fake)
  names(taxa_to_subset) <- taxa_names(physeq)

  subset_physeq <- subset_taxa_pq(physeq, taxa_to_subset)
  refseq <- subset_physeq@refseq
  fake_refseq <- Biostrings::DNAStringSet(lapply(refseq, sample))
  subset_physeq@refseq <- refseq(fake_refseq)

  taxa_names(subset_physeq) <- paste0(prefix, 1:n_fake)
  new_physeq <- merge_phyloseq(physeq, subset_physeq)

  new_tax_tab <- as.matrix(unclass(new_physeq@tax_table))
  new_tax_tab[taxa_names(subset_physeq), ] <- NA
  new_physeq@tax_table <- tax_table(new_tax_tab)

  new_otu_tab <- as.matrix(unclass(new_physeq@otu_table))
  new_otu_tab[, taxa_names(subset_physeq)] <- 0
  new_physeq@otu_table <- otu_table(new_otu_tab, taxa_are_rows = FALSE)
  if (taxasrow) {
    newphyseq <- taxa_as_rows(new_physeq)
  }
  return(new_physeq)
}

#' Add external sequences to a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful to compute true negative values with functions
#'   [tc_metrics_mock()] and [tc_metrics_mock_vec()].
#'
#'   Note the that the tax_table for additional sequences is full of NA and
#'   that the corresponding otu_table is full of 0.
#'
#' @inheritParams tc_points_matrix
#' @param ext_seqs (DNAStringSet, required) A DNAStringSet object containing
#'   external sequences to add.
#' @param prefix (character, default "external_") A prefix to add to the taxa name.
#' @returns A phyloseq object
#'
#' @export
#' @author Adrien Taudière
#' @examples
#' add_external_seq_pq(
#'   data_fungi_mini,
#'   Biostrings::readDNAStringSet(
#'     system.file("extdata/ex_little.fasta", package = "MiscMetabar")
#'   )
#' )
#'
add_external_seq_pq <- function(physeq, ext_seqs, prefix = "external_") {
  refseq <- physeq@refseq
  names(ext_seqs) <- paste0(prefix, names(ext_seqs))
  external_samtab <- physeq@sam_data

  external_otu_tab <- otu_table(
    matrix(0, nrow = nsamples(physeq), ncol = length(ext_seqs)),
    taxa_are_rows = FALSE
  )
  colnames(external_otu_tab) <- names(ext_seqs)
  rownames(external_otu_tab) <- sample_names(physeq)

  fake_tax_tab <- matrix(
    "NA",
    nrow = length(ext_seqs),
    ncol = ncol(physeq@tax_table)
  )
  colnames(fake_tax_tab) <- colnames(physeq@tax_table)
  rownames(fake_tax_tab) <- names(ext_seqs)

  fake_pq <- phyloseq(
    refseq(ext_seqs),
    otu_table(external_otu_tab),
    sample_data(external_samtab),
    tax_table(fake_tax_tab)
  )

  new_physeq <- merge_phyloseq(physeq, fake_pq)
  new_tax_tab <- as.matrix(unclass(new_physeq@tax_table))
  new_tax_tab[taxa_names(fake_pq), ] <- NA
  new_physeq@tax_table <- tax_table(new_tax_tab)

  return(new_physeq)
}


#' Multiply OTU counts conditionally based on sample metadata
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Multiplies OTU table counts by different values depending on the
#' level of a factor column in `sample_data`. Samples whose factor
#' value does not match any condition (including `NA` samples) are
#' left unchanged unless explicitly targeted.
#'
#' @inheritParams tc_points_matrix
#' @param fact (character, required) The name of a column in
#'   `sample_data(physeq)`.
#' @param conditions (character vector, required) Values of `fact` to
#'   match. Use `NA` to target samples with missing values.
#' @param multipliers (numeric vector, required) Multiplying factors
#'   corresponding to each element of `conditions`. Must be the same
#'   length as `conditions`.
#' @param prop_taxa (numeric, default 0.5) Proportion of taxa (between
#'   0 and 1) to which the multiplication is applied. Taxa are
#'   randomly sampled. Set to 1 to apply to all taxa.
#' @param seed (integer, default NULL) Random seed for reproducible
#'   taxa sampling when `prop_taxa < 1`.
#' @param compensate (logical, default FALSE) If TRUE, scale down the
#'   non-selected taxa so that total library size per sample is
#'   preserved. This creates a pure compositional shift: selected
#'   taxa gain relative abundance while the rest lose it. This mode
#'   is essential for testing differential abundance methods (ALDEx2,
#'   ANCOM-BC) which normalize for library size.
#' @param min_prevalence (numeric, default 0) Minimum prevalence
#'   (proportion of matched samples with non-zero counts) for a taxon
#'   to be eligible for selection. Setting this > 0 (e.g. 0.5) ensures
#'   that only taxa actually present in the target samples are
#'   inflated, producing a stronger and more realistic DA signal.
#' @param round (logical, default TRUE) If TRUE, round the resulting
#'   counts to integers (since OTU tables typically contain integer
#'   counts).
#' @param verbose (logical, default TRUE) If TRUE, print a message
#'   with the number of modified taxa.
#'
#' @return A phyloseq object with modified OTU counts. The selected taxa
#'   are stored as an attribute accessible via `attr(result, "taxa_modified")`.
#'
#' @details
#' Each sample is checked against `conditions` in order. When the
#' sample's value in `fact` matches a condition, its counts are
#' multiplied by the corresponding multiplier. Samples that do not
#' match any condition are left unchanged.
#'
#' When `prop_taxa < 1`, only a random subset of taxa is affected by
#' the multiplication. This is useful to simulate differential
#' abundance where only some taxa respond to a condition.
#'
#' **Simulating differential abundance for DA methods:**
#' DA methods like ALDEx2 and ANCOM-BC work on compositions (relative
#' abundances). Simply multiplying counts changes library size, which
#' these methods normalize away. Use `compensate = TRUE` to create a
#' real compositional shift: selected taxa are inflated, then non-
#' selected taxa are scaled down so total counts per sample stay
#' constant. Combine with `min_prevalence > 0` and small `prop_taxa`
#' (e.g. 0.05) for a realistic DA signal on few, prevalent taxa.
#'
#' **Compensation limits:** When the multiplier is very high and
#' selected taxa already have high counts, the non-selected taxa may
#' not have enough counts to compensate. In such cases, non-selected
#' taxa are zeroed and selected taxa are scaled down to preserve
#' library size. A warning is issued when this occurs. To avoid this,
#' use moderate multipliers (2-3) or select taxa with lower initial
#' abundance.
#'
#' To target `NA` values in the factor column, include `NA` in
#' `conditions` (e.g., `conditions = c("High", NA)`).
#'
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' # Multiply counts by 1.2 for "High" samples in half of taxa
#' d <- multiply_counts_pq(data_fungi,
#'   fact = "Height",
#'   conditions = "High", multipliers = 2
#' )
#'
#' # Multiple conditions with different multipliers
#' d2 <- multiply_counts_pq(data_fungi,
#'   fact = "Height",
#'   conditions = c("High", "Low"), multipliers = c(1.5, 0.8)
#' )
#'
#' # Target NA samples
#' d3 <- multiply_counts_pq(data_fungi,
#'   fact = "Height",
#'   conditions = c("High", NA), multipliers = c(1.2, 0.5)
#' )
#'
#' sum(data_fungi@otu_table)
#' sum(d@otu_table)
#' sum(d2@otu_table)
#' sum(d3@otu_table)
#'
#' # Apply only to 30% of taxa with a reproducible seed set
#' d4 <- multiply_counts_pq(data_fungi,
#'   fact = "Height",
#'   conditions = "High", multipliers = 2, prop_taxa = 0.3, seed = 42
#' )
#'
#' # Simulate DA for compositional methods (ALDEx2, ANCOM-BC):
#' # Use compensate=TRUE to create a real compositional shift
#' # Use min_prevalence to select only prevalent taxa
#' # Use small prop_taxa to affect few taxa strongly
#' d_da <- multiply_counts_pq(
#'   data_fungi,
#'   fact = "Height",
#'   conditions = "High",
#'   multipliers = 5,
#'   prop_taxa = 0.05,
#'   min_prevalence = 0.5,
#'   compensate = TRUE,
#'   seed = 123
#' )
#'
#' # Library sizes are preserved
#' sum(sample_sums(data_fungi))
#' sum(sample_sums(d_da))
#'
#' # Test with maaslin3
#' res <- maaslin3_pq(d_da,
#'   formula = "~Height"
#' )
#' 
#' gg_maaslin3_plot(res, type = "volcano")
#'
multiply_counts_pq <- function(
  physeq,
  fact,
  conditions,
  multipliers,
  prop_taxa = 0.5,
  seed = NULL,
  compensate = FALSE,
  min_prevalence = 0,
  round = TRUE,
  verbose = TRUE
) {
  if (length(conditions) != length(multipliers)) {
    stop(
      "'conditions' and 'multipliers' must have the same length. Got ",
      length(conditions),
      " and ",
      length(multipliers),
      "."
    )
  }

  verify_pq(physeq)

  if (prop_taxa <= 0 || prop_taxa > 1) {
    stop("'prop_taxa' must be between 0 (exclusive) and 1 (inclusive).")
  }

  if (min_prevalence < 0 || min_prevalence > 1) {
    stop("'min_prevalence' must be between 0 and 1.")
  }

  sd <- as.data.frame(phyloseq::sample_data(physeq))
  if (!fact %in% colnames(sd)) {
    stop("Column '", fact, "' not found in sample_data.")
  }

  fact_values <- sd[[fact]]
  taxasrow <- phyloseq::taxa_are_rows(physeq)
  otu <- as.matrix(phyloseq::otu_table(physeq))

  # Identify all matched samples across all conditions (for prevalence filter)
  all_matched <- rep(FALSE, length(fact_values))
  for (i in seq_along(conditions)) {
    cond <- conditions[i]
    if (is.na(cond)) {
      all_matched <- all_matched | is.na(fact_values)
    } else {
      all_matched <- all_matched | (!is.na(fact_values) & fact_values == cond)
    }
  }

  # Filter taxa by prevalence in matched samples
  all_taxa <- phyloseq::taxa_names(physeq)
  if (min_prevalence > 0 && any(all_matched)) {
    if (taxasrow) {
      matched_otu <- otu[, all_matched, drop = FALSE]
      prevalence <- rowSums(matched_otu > 0) / sum(all_matched)
    } else {
      matched_otu <- otu[all_matched, , drop = FALSE]
      prevalence <- colSums(matched_otu > 0) / sum(all_matched)
    }
    eligible_taxa <- all_taxa[prevalence >= min_prevalence]
    if (length(eligible_taxa) == 0) {
      warning(
        "No taxa meet min_prevalence threshold of ",
        min_prevalence,
        ". Using all taxa."
      )
      eligible_taxa <- all_taxa
    }
  } else {
    eligible_taxa <- all_taxa
  }

  # Select taxa subset from eligible taxa
  n_taxa <- ceiling(length(eligible_taxa) * prop_taxa)
  n_taxa <- min(n_taxa, length(eligible_taxa))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  selected_taxa <- sample(eligible_taxa, n_taxa)

  # Store original library sizes for compensation
  if (compensate) {
    if (taxasrow) {
      original_libsize <- colSums(otu)
    } else {
      original_libsize <- rowSums(otu)
    }
  }

  for (i in seq_along(conditions)) {
    cond <- conditions[i]
    mult <- multipliers[i]

    if (is.na(cond)) {
      matched <- is.na(fact_values)
    } else {
      matched <- !is.na(fact_values) & fact_values == cond
    }

    if (taxasrow) {
      taxa_idx <- rownames(otu) %in% selected_taxa
      otu[taxa_idx, matched] <- otu[taxa_idx, matched] * mult
    } else {
      taxa_idx <- colnames(otu) %in% selected_taxa
      otu[matched, taxa_idx] <- otu[matched, taxa_idx] * mult
    }
  }

  # Compensate: scale down non-selected taxa to preserve library size
  if (compensate) {
    n_uncompensated <- 0

    for (i in seq_along(conditions)) {
      cond <- conditions[i]

      if (is.na(cond)) {
        matched <- is.na(fact_values)
      } else {
        matched <- !is.na(fact_values) & fact_values == cond
      }

      matched_idx <- which(matched)

      for (j in matched_idx) {
        if (taxasrow) {
          taxa_idx <- rownames(otu) %in% selected_taxa
          new_libsize <- sum(otu[, j])
          target_libsize <- original_libsize[j]

          if (new_libsize > target_libsize) {
            nonselected_sum <- sum(otu[!taxa_idx, j])
            excess <- new_libsize - target_libsize

            if (nonselected_sum > 0 && excess <= nonselected_sum) {
              scale_factor <- (nonselected_sum - excess) / nonselected_sum
              otu[!taxa_idx, j] <- otu[!taxa_idx, j] * scale_factor
            } else if (excess > nonselected_sum) {
              # Cannot fully compensate - zero out non-selected and cap selected
              otu[!taxa_idx, j] <- 0
              # Scale down selected taxa to reach target
              selected_sum <- sum(otu[taxa_idx, j])
              if (selected_sum > 0) {
                otu[taxa_idx, j] <- otu[taxa_idx, j] *
                  (target_libsize / selected_sum)
              }
              n_uncompensated <- n_uncompensated + 1
            }
          }
        } else {
          taxa_idx <- colnames(otu) %in% selected_taxa
          new_libsize <- sum(otu[j, ])
          target_libsize <- original_libsize[j]

          if (new_libsize > target_libsize) {
            nonselected_sum <- sum(otu[j, !taxa_idx])
            excess <- new_libsize - target_libsize

            if (nonselected_sum > 0 && excess <= nonselected_sum) {
              scale_factor <- (nonselected_sum - excess) / nonselected_sum
              otu[j, !taxa_idx] <- otu[j, !taxa_idx] * scale_factor
            } else if (excess > nonselected_sum) {
              # Cannot fully compensate - zero out non-selected and cap selected
              otu[j, !taxa_idx] <- 0
              # Scale down selected taxa to reach target
              selected_sum <- sum(otu[j, taxa_idx])
              if (selected_sum > 0) {
                otu[j, taxa_idx] <- otu[j, taxa_idx] *
                  (target_libsize / selected_sum)
              }
              n_uncompensated <- n_uncompensated + 1
            }
          }
        }
      }
    }

    if (n_uncompensated > 0) {
      warning(
        n_uncompensated,
        " sample(s) could not be fully compensated ",
        "(multiplier too high relative to non-selected taxa abundance). ",
        "Non-selected taxa were zeroed and selected taxa were scaled down ",
        "to preserve library size."
      )
    }
  }

  if (round) {
    otu <- round(otu)
  }

  phyloseq::otu_table(physeq) <- phyloseq::otu_table(
    otu,
    taxa_are_rows = taxasrow
  )

  # Store selected taxa as attribute for reference
  attr(physeq, "taxa_modified") <- selected_taxa

  if (verbose) {
    message(
      "Modified ",
      length(selected_taxa),
      " taxa in ",
      sum(all_matched),
      " matched samples"
    )
  }

  return(physeq)
}


#' Simulate differential abundance by redistributing OTU counts
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Creates a fake differential abundance signal by redistributing counts
#' within matched samples: selected "DA" taxa receive counts taken from
#' non-selected taxa. Library sizes (row sums) are preserved, but taxa
#' totals (column sums) are allowed to change to create a detectable DA
#' signal.
#'
#' This approach is suitable for testing DA methods (ALDEx2, ANCOM-BC)
#' because it creates a true compositional shift while maintaining the
#' library size normalization these methods rely on.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param fact (character, required) Column name in `sample_data` for
#'   grouping.
#' @param conditions (character vector, required) Levels of `fact` where
#'   DA should be simulated.
#' @param effect_size (numeric, default 3) The fold-change in relative
#'   abundance for selected taxa in matched samples. Values > 1 increase
#'   abundance.
#' @param prop_taxa (numeric, default 0.05) Proportion of taxa to be
#'   selected as DA taxa.
#' @param min_prevalence (numeric, default 0.1) Minimum prevalence in
#'   matched samples for a taxon to be eligible.
#' @param seed (integer, default NULL) Random seed for reproducibility.
#' @param verbose (logical, default FALSE) Print progress info.
#'
#' @return A phyloseq object with modified OTU counts creating a DA signal.
#'   The selected DA taxa are stored in `attr(result, "da_taxa")`.
#'
#' @details
#' The algorithm works as follows:
#' 1. Select a subset of prevalent taxa as "DA" taxa
#' 2. For matched samples only:
#'    - Multiply selected taxa counts by `effect_size`
#'    - Scale down non-selected taxa to preserve library size
#' 3. Round to integers
#'
#' Unlike [multiply_counts_pq()], this function guarantees that library
#' sizes are exactly preserved (except for rounding), creating a pure
#' compositional shift that DA methods can detect.
#'
#' @export
#' @author Adrien Taudière
#'
#' @seealso [multiply_counts_pq()]
#'
#' @examples
#' \dontrun{
#' d_da <- permute_da_pq(
#'   data_fungi,
#'   fact = "Height",
#'   conditions = "High",
#'   effect_size = 3,
#'   prop_taxa = 0.05,
#'   seed = 123
#' )
#'
#' # Library sizes are preserved
#' all.equal(sample_sums(data_fungi), sample_sums(d_da))
#'
#' # Test with ALDEx2
#' res <- MiscMetabar::aldex_pq(d_da,
#'   bifactor = "Height",
#'   modalities = c("Low", "High")
#' )
#' gg_aldex_plot(res, type = "volcano")
#'
#' # Check which taxa were made DA
#' attr(d_da, "da_taxa")
#' }
permute_da_pq <- function(
  physeq,
  fact,
  conditions,
  effect_size = 3,
  prop_taxa = 0.05,
  min_prevalence = 0.1,
  seed = NULL,
  verbose = FALSE
) {
  verify_pq(physeq)

  if (effect_size <= 0) {
    stop("'effect_size' must be positive.")
  }

  sd <- as.data.frame(phyloseq::sample_data(physeq))
  if (!fact %in% colnames(sd)) {
    stop("Column '", fact, "' not found in sample_data.")
  }

  fact_values <- sd[[fact]]
  taxasrow <- phyloseq::taxa_are_rows(physeq)
  otu <- unclass(as.matrix(phyloseq::otu_table(physeq)))

  # Ensure samples are rows for easier manipulation
  if (taxasrow) {
    otu <- t(otu)
  }

  # Original library sizes
  orig_libsize <- rowSums(otu)

  # Identify matched samples
  all_matched <- rep(FALSE, nrow(otu))
  for (cond in conditions) {
    if (is.na(cond)) {
      all_matched <- all_matched | is.na(fact_values)
    } else {
      all_matched <- all_matched | (!is.na(fact_values) & fact_values == cond)
    }
  }

  if (sum(all_matched) == 0) {
    warning("No samples match the specified conditions.")
    return(physeq)
  }

  # Filter taxa by prevalence in matched samples
  all_taxa <- phyloseq::taxa_names(physeq)
  matched_otu <- otu[all_matched, , drop = FALSE]
  prevalence <- colSums(matched_otu > 0) / sum(all_matched)
  eligible_taxa <- all_taxa[prevalence >= min_prevalence]

  if (length(eligible_taxa) == 0) {
    warning(
      "No taxa meet min_prevalence threshold of ",
      min_prevalence,
      ". Using top 10% most prevalent taxa."
    )
    eligible_taxa <- all_taxa[prevalence >= quantile(prevalence, 0.9)]
  }

  # Select DA taxa
  n_taxa <- max(1, ceiling(length(eligible_taxa) * prop_taxa))
  n_taxa <- min(n_taxa, length(eligible_taxa))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  selected_taxa <- sample(eligible_taxa, n_taxa)
  taxa_idx <- colnames(otu) %in% selected_taxa

  if (verbose) {
    message(
      "Selected ",
      length(selected_taxa),
      " DA taxa from ",
      length(eligible_taxa),
      " eligible taxa"
    )
    message("Matched samples: ", sum(all_matched))
  }

  # Apply effect to matched samples and compensate to preserve library size
  otu_result <- otu
  matched_idx <- which(all_matched)

  for (j in matched_idx) {
    # Original counts
    selected_counts <- otu[j, taxa_idx]
    nonselected_counts <- otu[j, !taxa_idx]
    target_libsize <- orig_libsize[j]

    # Apply effect to selected taxa
    new_selected <- selected_counts * effect_size

    # Calculate how much we need from non-selected to keep library size
    new_selected_sum <- sum(new_selected)
    remaining <- target_libsize - new_selected_sum

    if (remaining <= 0) {
      # Selected taxa alone exceed library size - scale them down
      otu_result[j, taxa_idx] <- new_selected *
        (target_libsize / new_selected_sum)
      otu_result[j, !taxa_idx] <- 0
    } else {
      # Scale down non-selected to fit remaining budget
      nonselected_sum <- sum(nonselected_counts)
      if (nonselected_sum > 0) {
        scale_factor <- remaining / nonselected_sum
        otu_result[j, taxa_idx] <- new_selected
        otu_result[j, !taxa_idx] <- nonselected_counts * scale_factor
      } else {
        otu_result[j, taxa_idx] <- new_selected
      }
    }
  }

  # Round
  otu_result <- round(otu_result)

  # Restore orientation
  if (taxasrow) {
    otu_result <- t(otu_result)
  }

  phyloseq::otu_table(physeq) <- phyloseq::otu_table(
    otu_result,
    taxa_are_rows = taxasrow
  )

  # Store selected taxa as attribute for reference
  attr(physeq, "da_taxa") <- selected_taxa

  return(physeq)
}


#' Simulate microbiome data with differential abundance using MIDASim
#'
#' #TODO : NOT WORKING with aldex2
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @description
#' Uses the MIDASim package to generate realistic simulated microbiome count
#' data with known differential abundance effects. The function learns the
#' structure (taxon correlations, sparsity patterns, library size distribution)
#' from a template phyloseq object and generates new data with specified DA
#' taxa.
#'
#' @param physeq (phyloseq, required) A phyloseq object to use as template.
#'   MIDASim will learn taxon-taxon correlations and abundance distributions
#'   from this data.
#' @param fact (character, required) Column name in `sample_data` defining
#'   the binary factor for DA simulation.
#' @param condition (character, required) The level of `fact` where taxa
#'   should be differentially abundant (case group).
#' @param n_samples (integer, default NULL) Number of samples to simulate.
#'   If NULL, uses the same number as the template.
#' @param prop_case (numeric, default 0.5) Proportion of simulated samples
#'   that belong to the case group (where DA taxa are elevated).
#' @param n_da_taxa (integer, default 10) Number of taxa to make
#'   differentially abundant.
#' @param effect_size (numeric, default 1) Log-fold change for DA taxa in

#'   case samples. Positive values increase abundance; negative decrease.
#'   Typical values: 0.5-2 for moderate effects.
#' @param da_taxa_idx (integer vector, default NULL) Specific taxon indices
#'   to make DA. If NULL, randomly selects `n_da_taxa` from prevalent taxa.
#' @param min_prevalence (numeric, default 0.1) Minimum prevalence for a
#'   taxon to be eligible for DA selection.
#' @param mode (character, default "nonparametric") MIDASim mode. One of
#'   "nonparametric" (more realistic) or "parametric" (more controllable).
#' @param seed (integer, default NULL) Random seed for reproducibility.
#' @param verbose (logical, default FALSE) Print progress messages.
#'
#' @return A phyloseq object with simulated counts. The object includes:
#'   - Simulated OTU table with DA effects
#'   - Tax table from template (taxa are preserved)
#'   - New sample_data with binary factor column
#'   - Attributes: `da_taxa_idx` (indices), `da_taxa_names` (names),
#'     `effect_size`, `condition`
#'
#' @details
#' The simulation workflow follows MIDASim's three-step process:
#'
#' 1. **Setup**: `MIDASim.setup()` learns the template's taxon abundance
#'    distributions, sparsity patterns, and taxon-taxon correlations.
#'
#' 2. **Modify**: `MIDASim.modify()` introduces DA effects by creating
#'    sample-specific relative abundances. For case samples (where
#'    `fact == condition`), selected taxa are multiplied by
#'    `exp(effect_size)`, then renormalized to sum to 1.
#'
#' 3. **Simulate**: `MIDASim()` generates realistic count data preserving
#'    the template's correlation structure.
#'
#' This approach is superior to simple count multiplication because:
#' - Taxon-taxon correlations are preserved
#' - Sparsity patterns are realistic
#' - Library sizes follow realistic distributions
#' - The DA signal is embedded in the data generation process
#'
#' @export
#' @author Adrien Taudière
#'
#' @seealso [multiply_counts_pq()], [permute_da_pq()]
#'
#' @references
#' He M, et al. (2024). MIDASim: a fast and simple simulator for realistic
#' microbiome data. Microbiome. doi:10.1186/s40168-024-01822-z
#'
#' @examples
#' # Requires MIDASim package
#' # install.packages("MIDASim")
#'
#' # Simulate data with 10 DA taxa having log-fold change of 2
#' sim_pq <- midasim_pq(
#'   data_fungi_mini,
#'   fact = "Height",
#'   condition = "High",
#'   n_da_taxa = 10,
#'   effect_size = 2,
#'   seed = 42
#' )
#'
#' # Test with ALDEx2
#' sim_pq@sam_data$Height <-
#'   as.character(sim_pq@sam_data$Height)
#' res <- MiscMetabar::aldex_pq(
#'   sim_pq,
#'   bifactor = "Height",
#'   modalities = c("Low", "High")
#' )
#' gg_aldex_plot(res, type = "volcano")
#'
#' # Check which taxa were set as DA
#' attr(sim_pq, "da_taxa_names")
#'
#' res_height <- ancombc_pq(
#'   sim_pq,
#'   fact = "Height",
#'   levels_fact = c("Low", "High"),
#'   verbose = TRUE, tax_level = NULL
#' )
#'
#' ggplot(
#'   res_height$res,
#'   aes(
#'     y = reorder(taxon, lfc_HeightHigh),
#'     x = lfc_HeightHigh,
#'     color = diff_HeightHigh
#'   )
#' ) +
#'   geom_vline(xintercept = 0) +
#'   geom_segment(aes(
#'     xend = 0, y = reorder(taxon, lfc_HeightHigh),
#'     yend = reorder(taxon, lfc_HeightHigh)
#'   ), color = "darkgrey") +
#'   geom_point()
midasim_pq <- function(
  physeq,
  fact,
  condition,
  n_samples = NULL,
  prop_case = 0.5,
  n_da_taxa = 10,
  effect_size = 1,
  da_taxa_idx = NULL,
  min_prevalence = 0.1,
  mode = c("nonparametric", "parametric"),
  seed = NULL,
  verbose = TRUE
) {
  mode <- match.arg(mode)

  verify_pq(physeq)
  physeq <- taxa_as_columns(physeq)

  # Check MIDASim availability
  if (!requireNamespace("MIDASim", quietly = TRUE)) {
    stop(
      "Package 'MIDASim' is required. Install it with:\n",
      "  install.packages('MIDASim')"
    )
  }

  sd <- as.data.frame(phyloseq::sample_data(physeq))
  if (!fact %in% colnames(sd)) {
    stop("Column '", fact, "' not found in sample_data.")
  }

  if (!condition %in% sd[[fact]]) {
    stop("Condition '", condition, "' not found in column '", fact, "'.")
  }

  # Extract OTU table (samples as rows, taxa as columns)
  taxasrow <- phyloseq::taxa_are_rows(physeq)
  otu <- as.matrix(phyloseq::otu_table(physeq))
  if (taxasrow) {
    otu <- t(otu)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Number of samples and taxa
  if (is.null(n_samples)) {
    n_samples <- nrow(otu)
  } else {
    set.seed(seed)
    physeq <- prune_samples(
      sample_names(physeq)[sample(1:nrow(otu), n_samples, replace = FALSE)],
      physeq
    )
  }
  n_taxa <- ncol(otu)
  taxa_names_vec <- colnames(otu)

  if (verbose) {
    message("Template: ", nrow(otu), " samples, ", n_taxa, " taxa")
    message("Simulating: ", n_samples, " samples")
  }

  # Step 1: MIDASim.setup
  if (verbose) {
    message("Running MIDASim.setup (mode = ", mode, ")...")
  }

  setup_args <- list(otu.tab = otu, mode = mode)
  if (mode == "nonparametric") {
    setup_args$n.break.ties <- 10
  }
  midas_setup <- do.call(MIDASim::MIDASim.setup, setup_args)

  # Determine DA taxa
  if (is.null(da_taxa_idx)) {
    # Select based on prevalence
    prevalence <- colSums(otu > 0) / nrow(otu)
    eligible_idx <- which(prevalence >= min_prevalence)

    if (length(eligible_idx) < n_da_taxa) {
      warning(
        "Only ",
        length(eligible_idx),
        " taxa meet prevalence threshold. ",
        "Using all eligible taxa."
      )
      da_taxa_idx <- eligible_idx
    } else {
      da_taxa_idx <- sample(eligible_idx, n_da_taxa)
    }
  }

  if (verbose) {
    message("Selected ", length(da_taxa_idx), " DA taxa")
  }

  # Step 2: Create sample-specific relative abundances with DA effect
  # Y = 0 for control, Y = 1 for case
  n_case <- round(n_samples * prop_case)
  n_control <- n_samples - n_case
  Y <- c(rep(0, n_control), rep(1, n_case))

  # Beta vector: effect_size for DA taxa, 0 otherwise
  beta <- rep(0, n_taxa)
  beta[da_taxa_idx] <- effect_size

  # Create individual relative abundances
  mean_rel_abund <- midas_setup$mean.rel.abund
  individual_rel_abund <- matrix(
    mean_rel_abund,
    nrow = n_samples,
    ncol = n_taxa,
    byrow = TRUE
  )

  # Apply DA effect: multiply by exp(Y * beta)
  # For controls (Y=0): no change

  # For cases (Y=1): multiply DA taxa by exp(effect_size)
  effect_matrix <- exp(outer(Y, beta))
  individual_rel_abund <- individual_rel_abund * effect_matrix
  individual_rel_abund <- individual_rel_abund / rowSums(individual_rel_abund)

  # Step 3: MIDASim.modify
  if (verbose) {
    message("Running MIDASim.modify...")
  }

  midas_modified <- MIDASim::MIDASim.modify(
    midas_setup,
    individual.rel.abund = individual_rel_abund
  )

  # Step 4: Simulate
  if (verbose) {
    message("Running MIDASim...")
  }

  sim_result <- MIDASim::MIDASim(midas_modified)

  # Extract simulated counts
  sim_counts <- sim_result$sim_count
  colnames(sim_counts) <- taxa_names_vec
  rownames(sim_counts) <- paste0("sim_sample_", seq_len(n_samples))

  # Create sample_data
  other_level <- setdiff(unique(sd[[fact]]), condition)[1]
  if (is.na(other_level)) {
    other_level <- "control"
  }

  sim_fact <- ifelse(Y == 1, condition, other_level)
  sim_sam_data <- data.frame(
    row.names = rownames(sim_counts),
    stringsAsFactors = FALSE
  )
  sim_sam_data[[fact]] <- factor(sim_fact, levels = c(other_level, condition))

  # Build phyloseq object
  new_otu <- phyloseq::otu_table(sim_counts, taxa_are_rows = FALSE)

  # Keep tax_table if available
  if (!is.null(physeq@tax_table)) {
    new_tax <- physeq@tax_table
  } else {
    new_tax <- NULL
  }

  new_pq <- phyloseq::phyloseq(
    new_otu,
    phyloseq::sample_data(sim_sam_data)
  )

  if (!is.null(new_tax)) {
    new_pq <- phyloseq::merge_phyloseq(new_pq, new_tax)
  }

  # Store DA info as attributes
  attr(new_pq, "da_taxa_idx") <- da_taxa_idx
  attr(new_pq, "da_taxa_names") <- taxa_names_vec[da_taxa_idx]
  attr(new_pq, "effect_size") <- effect_size
  attr(new_pq, "condition") <- condition
  attr(new_pq, "Y") <- Y

  if (verbose) {
    message(
      "Done. Simulated ",
      n_samples,
      " samples with ",
      length(da_taxa_idx),
      " DA taxa (effect = ",
      effect_size,
      ")"
    )
  }

  return(new_pq)
}
