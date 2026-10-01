################################################################################
#' Accuracy metrics of a taxonomic assignation against a per-unit truth
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Score one assignation column at one rank as Hleap et al. (2021) do
#' (`optimize_n_score.py::score()`): every scored unit falls in exactly one
#' cell of the confusion matrix, so TP + FP + FN + TN is the number of scored
#' units.
#'
#' | Unit | Value at the rank | Cell |
#' | --- | --- | --- |
#' | real, with a truth at this rank | its own truth | TP |
#' | real, with a truth at this rank | another name | FP |
#' | real, with a truth at this rank | NA | FN |
#' | real, without truth | any name / NA | FP / FN |
#' | real, truth shallower than the rank | — | left out |
#' | shuffled control (`fake_pattern`) | any name / NA | FP / TN |
#' | external control (`external_pattern`) | any name / NA | FP / TN when `external_scoring = "matrix"`, left out otherwise |
#'
#' This differs from [tc_metrics_mock_vec()], which compares each value to the
#' *set* of the expected taxa, counts FN on the rows of the truth table and
#' leaves the controls out of FP. Use this one when every unit has its own
#' truth (e.g. matched against the Sanger sequences of the mock strains), and
#' the other when the mock only comes with a list of expected taxa.
#'
#' @inheritParams tc_points_matrix
#' @param taxonomic_rank (required) Name (or number) of the `tax_table` column
#'   holding the assignation to score.
#' @param truth (required) A data.frame with one row per unit: a `taxon` column
#'   naming the unit, one column per rank of `truth_ranks` holding its truth,
#'   and a `truth_depth` column naming the deepest rank with a truth (`NA` when
#'   the unit has none). A unit is scored down to `truth_depth` and left out of
#'   the matrix below it. Units absent from `truth` are scored as units without
#'   truth. Optional columns `<rank>_accepted` hold a second name that counts
#'   as correct too (e.g. the accepted name of a synonym).
#' @param rank (character, default `taxonomic_rank`) Rank of `truth` this
#'   column is compared to, when the column is named after the method rather
#'   than the rank (e.g. `"Genus_dada2__unite"` against `"Genus"`).
#' @param truth_ranks (character vector) The ranks of `truth`, from the highest
#'   to the lowest, used to tell whether `rank` is deeper than `truth_depth`.
#' @param accepted_suffix (character, default "_accepted") Suffix of the
#'   optional second-name columns of `truth`.
#' @param external_scoring ("matrix" or "aside") Whether the external controls
#'   enter the confusion matrix. They should when the reference database cannot
#'   hold them (a Fungi-only database), and stay aside when it can, where
#'   naming them is a correct answer rather than an error.
#' @param fake_taxa (logical, default TRUE) If TRUE, the controls are
#'   identified by `fake_pattern` and `external_pattern` and scored as above.
#'   If FALSE, every unit is treated as a real one.
#' @param fake_pattern (character, default "^fake_") Regular expression
#'   identifying the shuffled controls ([add_shuffle_seq_pq()]).
#' @param external_pattern (character, default "^external_") Regular expression
#'   identifying the external controls ([add_external_seq_pq()]).
#' @param verbose (logical, default TRUE) If TRUE, print informative messages.
#'
#' @returns A list of metrics:
#'
#'  - TP, FP, FN, TN: the four cells, in units;
#'
#'  - FDR = FP / (FP + TP), PPV = TP / (TP + FP), TPR = TP / (TP + FN),
#'    TNR = TN / (TN + FP) (every FP, controls and real misassignments alike),
#'    F1_score = 2 TP / (2 TP + FP + FN), ACC and MCC (0 when a margin is zero,
#'    Chicco & Jurman 2020);
#'
#'  - NA_real: share of the scored real units left NA;
#'
#'  - NA_fake: share of the shuffled controls left NA (1 is the correct answer);
#'
#'  - ctrl_assigned_fake, ctrl_assigned_ext: number of controls given a value;
#'
#'  - misassign_seq: share of the scored real units given a wrong name;
#'
#'  - n_real, n_left_out: number of real units scored, and of units left out of
#'    the matrix at this rank (truth shallower than the rank).
#'
#' @export
#' @seealso [tc_metrics_unit()], [tc_metrics_mock_vec()]
#' @author Adrien Taudière
#' @examples
#' truth <- data.frame(
#'   taxon = phyloseq::taxa_names(data_fungi_mini),
#'   Phylum = phyloseq::tax_table(data_fungi_mini)[, "Phylum"],
#'   truth_depth = "Phylum"
#' )
#' tc_metrics_unit_vec(
#'   data_fungi_mini,
#'   taxonomic_rank = "Phylum",
#'   truth = truth,
#'   truth_ranks = "Phylum",
#'   fake_taxa = FALSE,
#'   verbose = FALSE
#' )
tc_metrics_unit_vec <- function(
  physeq,
  taxonomic_rank,
  truth,
  rank = taxonomic_rank,
  truth_ranks = c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  ),
  accepted_suffix = "_accepted",
  external_scoring = c("matrix", "aside"),
  fake_taxa = TRUE,
  fake_pattern = "^fake_",
  external_pattern = "^external_",
  verbose = TRUE
) {
  external_scoring <- match.arg(external_scoring)
  if (!rank %in% truth_ranks) {
    stop("`rank` (", rank, ") is not one of `truth_ranks`.")
  }
  if (!all(c("taxon", "truth_depth", rank) %in% colnames(truth))) {
    stop(
      "`truth` must have the columns taxon, truth_depth and ",
      rank,
      "."
    )
  }

  units <- taxa_names(physeq)
  values <- as.character(physeq@tax_table[, taxonomic_rank])
  names(values) <- units

  is_fake <- fake_taxa & grepl(fake_pattern, units)
  is_ext <- fake_taxa & grepl(external_pattern, units)
  is_real <- !is_fake & !is_ext
  if (verbose) {
    message(
      sum(is_fake),
      " shuffled and ",
      sum(is_ext),
      " external controls found; the external ones ",
      ifelse(
        external_scoring == "matrix",
        "enter the matrix",
        "stay aside"
      ),
      "."
    )
  }

  rows <- match(units, truth$taxon)
  depth <- truth$truth_depth[rows]
  expected <- truth[[rank]][rows]
  accepted_col <- paste0(rank, accepted_suffix)
  accepted <- if (accepted_col %in% colnames(truth)) {
    truth[[accepted_col]][rows]
  } else {
    rep(NA_character_, length(units))
  }

  # A unit is left out of the matrix at the ranks below its truth (a tie
  # between strains of the same genus says nothing about the species).
  deeper_than_truth <- is_real &
    !is.na(depth) &
    match(rank, truth_ranks) > match(depth, truth_ranks)
  scored_real <- is_real & !deeper_than_truth
  has_truth <- scored_real & !is.na(expected)

  correct <- !is.na(values) &
    ((!is.na(expected) & values == expected) |
      (!is.na(accepted) & values == accepted))

  TP <- sum(has_truth & correct)
  FP_real <- sum(scored_real & !is.na(values) & !correct)
  FN <- sum(scored_real & is.na(values))

  scored_ctrl <- is_fake | (is_ext & external_scoring == "matrix")
  FP_ctrl <- sum(scored_ctrl & !is.na(values))
  TN <- sum(scored_ctrl & is.na(values))
  FP <- FP_real + FP_ctrl

  FDR <- FP / (FP + TP)
  PPV <- TP / (TP + FP)
  TPR <- TP / (TP + FN)
  TNR <- TN / (TN + FP)
  F1_score <- 2 * TP / (2 * TP + FP + FN)
  ACC <- (TP + TN) / (TP + TN + FP + FN)
  # as.numeric(): the product of the four margins overflows beyond
  # .Machine$integer.max. A zero margin leaves the MCC undefined: 0 by
  # convention (Chicco & Jurman 2020).
  mcc_denominator <- as.numeric(TP + FP) * (TP + FN) * (FP + TN) * (TN + FN)
  MCC <- if (mcc_denominator == 0) {
    0
  } else {
    (as.numeric(TP) * TN - as.numeric(FP) * FN) / sqrt(mcc_denominator)
  }

  list(
    "TP" = TP,
    "FP" = FP,
    "FN" = FN,
    "TN" = TN,
    "FDR" = FDR,
    "PPV" = PPV,
    "TPR" = TPR,
    "TNR" = TNR,
    "F1_score" = F1_score,
    "ACC" = ACC,
    "MCC" = MCC,
    "NA_real" = sum(scored_real & is.na(values)) / sum(scored_real),
    # NaN rather than a share when the object holds no shuffled control
    # (fake_taxa = FALSE, or an input without controls).
    "NA_fake" = sum(is_fake & is.na(values)) / sum(is_fake),
    "ctrl_assigned_fake" = sum(is_fake & !is.na(values)),
    "ctrl_assigned_ext" = sum(is_ext & !is.na(values)),
    "misassign_seq" = FP_real / sum(scored_real),
    "n_real" = sum(scored_real),
    "n_left_out" = sum(deeper_than_truth)
  )
}
################################################################################

################################################################################
#' Accuracy metrics of several assignations against a per-unit truth
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Run [tc_metrics_unit_vec()] on every assignation of `ranks_df` and every
#' rank of `truth`, and add the metrics of the external controls, which the
#' per-rank function cannot compute on its own:
#'
#'  - `ext_fungi`: share of the external controls named Fungi at the kingdom
#'    rank, an error on every database;
#'
#'  - `ext_correct`: share of the external controls given their own lineage,
#'    among those with a truth at the rank (NA in `external_truth`: left
#'    out, as a real unit below its truth depth); needs `external_truth`.
#'
#' @inheritParams tc_metrics_unit_vec
#' @param ranks_df (required) A data.frame of `tax_table` column names: one
#'   column per assignation (method x database x parameters), one row per rank
#'   of `truth_ranks`, in the same order.
#' @param truth (required) The per-unit truth, as in [tc_metrics_unit_vec()].
#' @param external_truth (data.frame, default NULL) The true lineage of the
#'   external controls: a `taxon` column and one column per rank. Without it
#'   `ext_correct` is not computed.
#' @param external_correct_ranks (character vector, default NULL) Ranks where
#'   `ext_correct` is computed. NULL means every rank of `truth_ranks`; a
#'   database holding only a few non-fungal representatives should restrict it
#'   to the kingdom.
#' @param kingdom_rank (character, default "Kingdom") Rank at which
#'   `ext_fungi` is counted.
#' @param fungi_name (character, default "Fungi") Value of `kingdom_rank`
#'   counted by `ext_fungi`.
#'
#' @returns A long-format data.frame with four columns: `method_db`,
#'   `tax_level`, `metrics` and `values`, as [tc_metrics_mock()] returns.
#' @export
#' @seealso [tc_metrics_unit_vec()], [tc_metrics_mock()]
#' @author Adrien Taudière
#' @examples
#' truth <- data.frame(
#'   taxon = phyloseq::taxa_names(data_fungi_mini),
#'   Phylum = phyloseq::tax_table(data_fungi_mini)[, "Phylum"],
#'   Class = phyloseq::tax_table(data_fungi_mini)[, "Class"],
#'   truth_depth = "Class"
#' )
#' ranks_df <- data.frame(seed = c("Phylum", "Class"))
#' tc_metrics_unit(
#'   data_fungi_mini,
#'   ranks_df = ranks_df,
#'   truth = truth,
#'   truth_ranks = c("Phylum", "Class"),
#'   fake_taxa = FALSE,
#'   verbose = FALSE
#' )
tc_metrics_unit <- function(
  physeq,
  ranks_df,
  truth,
  truth_ranks = c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  ),
  external_truth = NULL,
  external_correct_ranks = NULL,
  accepted_suffix = "_accepted",
  external_scoring = c("matrix", "aside"),
  fake_taxa = TRUE,
  fake_pattern = "^fake_",
  external_pattern = "^external_",
  kingdom_rank = "Kingdom",
  fungi_name = "Fungi",
  verbose = FALSE
) {
  external_scoring <- match.arg(external_scoring)
  if (nrow(ranks_df) != length(truth_ranks)) {
    stop(
      "`ranks_df` must have one row per rank of `truth_ranks` (",
      nrow(ranks_df),
      " rows against ",
      length(truth_ranks),
      " ranks)."
    )
  }
  if (is.null(external_correct_ranks)) {
    external_correct_ranks <- truth_ranks
  }

  units <- taxa_names(physeq)
  is_ext <- fake_taxa & grepl(external_pattern, units)

  res <- lapply(seq_len(ncol(ranks_df)), function(nc) {
    columns <- as.character(ranks_df[, nc])
    per_rank <- lapply(seq_along(truth_ranks), function(i) {
      metrics <- tc_metrics_unit_vec(
        physeq,
        taxonomic_rank = columns[i],
        truth = truth,
        rank = truth_ranks[i],
        truth_ranks = truth_ranks,
        accepted_suffix = accepted_suffix,
        external_scoring = external_scoring,
        fake_taxa = fake_taxa,
        fake_pattern = fake_pattern,
        external_pattern = external_pattern,
        verbose = verbose
      )
      values <- as.character(physeq@tax_table[, columns[i]])
      names(values) <- units

      if (sum(is_ext) > 0 && truth_ranks[i] == kingdom_rank) {
        metrics[["ext_fungi"]] <- sum(
          !is.na(values[is_ext]) & values[is_ext] == fungi_name
        ) /
          sum(is_ext)
      }
      if (
        sum(is_ext) > 0 &&
          !is.null(external_truth) &&
          truth_ranks[i] %in% external_correct_ranks
      ) {
        expected <- external_truth[[truth_ranks[i]]][
          match(units, external_truth$taxon)
        ]
        metrics[["ext_correct"]] <- sum(
          is_ext & !is.na(values) & !is.na(expected) & values == expected
        ) /
          sum(is_ext & !is.na(expected))
      }

      data.frame(
        "method_db" = rep(colnames(ranks_df)[nc], length(metrics)),
        "tax_level" = rep(truth_ranks[i], length(metrics)),
        "metrics" = names(metrics),
        "values" = as.vector(unlist(metrics))
      )
    })
    do.call(rbind, per_rank)
  })
  do.call(rbind, res)
}
################################################################################
