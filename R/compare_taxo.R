################################################################################
#' Compute accuracy metrics of multiple taxonomic assignations method using
#'  mock for multi-rank and multi assignation methods
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Compute numerous metrics comparing the computed taxonomic assignation to a
#' true assignation.
#'
#' Note that to compute all metrics, one need to insert fake
#' taxa (by shuffling sequences and/or by adding external sequences). The user
#' must fake taxa using functions [add_external_seq_pq()],
#' [add_shuffle_seq_pq()]) before taxonomic assignation.
#'
#' @inheritParams tc_points_matrix
#' @param ranks_df (required). A dataframe with at least one column (one database or one
#'   method) and a number of row equal to the column in true_values_df
#' @param true_values_df  (required). A dataframe with the true taxonomic assignation.
#'   Note that the column names (names of taxonomic ranks) of the true_values_df defined
#'   the names present in the `tax_level` column of the resulting dataframe.
#' @param fake_taxa (logical, default TRUE). If TRUE, the fake_pattern vector
#'   is used to identify fake taxa, i.e. taxa who are not in the reference
#'   database (see [add_external_seq_pq()]) or taxa with fake sequences
#'   (see [add_shuffle_seq_pq()]).
#' @param fake_pattern (vector of pattern): A vector used to identify the fake
#'   taxa using a regex search in their name.
#'
#' @returns A long-format dataframe with 4 columns: (i) the name of the `method_db`
#' (ii) the name of the `tax_level` (taxonomic rank), (iii) the `metrics` (see
#' [tc_metrics_mock_vec()] for more details) and (iv) the `values`.
#' @export
#' @seealso [tc_metrics_mock_vec()]
#' @author Adrien Taudière
tc_metrics_mock <- function(physeq,
                            ranks_df,
                            true_values_df,
                            fake_taxa = TRUE,
                            fake_pattern = c("^fake_", "^external_")) {
  if (nrow(ranks_df) != ncol(true_values_df)) {
    stop("The number of rows of ranks_df must be equal to the number of column in true_values_df")
  }

  res_df <- data.frame(matrix(NA, nrow = 0, ncol = 4))
  colnames(res_df) <- c("method_db", "tax_level", "metrics", "values")


  for (nc in 1:ncol(ranks_df)) {
    for (i in 1:nrow(ranks_df)) {
      res_list <- tc_metrics_mock_vec(
        physeq,
        taxonomic_rank = ranks_df[i, nc],
        true_values = true_values_df[, i],
        fake_taxa = fake_taxa,
        fake_pattern = fake_pattern
      )

      res_df <- rbind(
        res_df,
        data.frame(
          "method_db" = rep(colnames(ranks_df)[nc], length(res_list)),
          "tax_level" = rep(colnames(true_values_df)[i], length(res_list)),
          "metrics" = names(res_list),
          "values" = as.vector(unlist(res_list))
        )
      )
      # print(paste0(colnames(ranks_df)[nc], " - ", colnames(true_values_df)[i]))
    }
  }
  return(res_df)
}
################################################################################

################################################################################
#' Compute accuracy metrics of taxonomic assignation using a mock (known)
#'  community for one rank
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Compute numerous metrics comparing the computed taxonomic assignation to a
#' true assignation.
#'
#'
#' Note that to compute all metrics, one need to insert fake
#' taxa (by shuffling sequences and/or by adding external sequences). The user
#' must fake taxa using functions [add_external_seq_pq()],
#' [add_shuffle_seq_pq()]) before taxonomic assignation.
#'
#'
#' @inheritParams tc_points_matrix
#' @param taxonomic_rank (required) Name (or number) of a taxonomic rank
#'   to count.
#' @param true_values (required) A vector with the true taxonomic assignation
#' @param fake_taxa (logical, default TRUE). If TRUE, the fake_pattern vector
#'   is used to identify fake taxa, i.e. taxa who are not in the reference
#'   database (see [add_external_seq_pq()]) or taxa with fake sequences
#'   (see [add_shuffle_seq_pq()]).
#' @param fake_pattern (vector of pattern): A vector used to identify the fake
#'   taxa using a regex search in their name.
#'
#' @returns A list of metrics (see the confusion matrix
#'   [article](https://en.wikipedia.org/wiki/Confusion_matrix) on wikipedia):
#'
#'  - TP (number of *true positive*)
#'
#'  - FP (number of *false positive*)
#'
#'  - FN (number of *false negative*)
#'
#'  - FDR (*false discovery rate*) = FP / (FP + TP)
#'
#'  - TPR (*true positive rate*, also named *recall* or *sensitivity*)
#'
#'  - PPV (*positive predictive value*, also named *precision*) = TP / (TP + FP)
#'
#'  - F1_score (*F1 score*) = 2 * TP / (2 * TP + FP + FN)
#'
#'  If fake taxa are present and `fake_taxa` is true, other metrics are computed:
#'
#'  - TN (number of *true negative*)
#'
#'  - ACC (*Accuracy*) = (TP + TN) / (TP + TN + FP + FN)
#'
#'  - MCC (*Matthews correlation coefficient*) = (TP * TN - FP * FN) /
#'      sqrt((TP + FP) * (TP + FN) * (FP + TN) * (TN + FN))
#'
#' @export
#' @seealso [tc_metrics_mock()], [add_external_seq_pq()], [add_shuffle_seq_pq()])
#'
#' @author Adrien Taudière
tc_metrics_mock_vec <- function(physeq,
                                taxonomic_rank,
                                true_values,
                                fake_taxa = TRUE,
                                fake_pattern = c("^fake_", "^external_"),
                                verbose = TRUE) {
  if (fake_taxa) {
    fake_pattern <- paste(fake_pattern, collapse = "|")

    fake_taxa_names <- taxa_names(physeq)[grepl(fake_pattern, taxa_names(physeq))]

    fake_taxa_cond <- taxa_names(physeq) %in% fake_taxa_names
    names(fake_taxa_cond) <- taxa_names(physeq)

    if (verbose) {
      message(
        length(fake_taxa_names),
        "fake taxa were found using pattern",
        fake_pattern,
        "."
      )
    }

    # Hleap et al.  2021 [https://doi.org/10.1111/1755-0998.13407]
    rank_pq <- physeq@tax_table[!fake_taxa_cond, taxonomic_rank]

    TP <- sum(rank_pq %in% true_values)
    if (sum(is.na(rank_pq)) == length(rank_pq)) {
      FP <- 0
    } else {
      FP <- sum(!(rank_pq[!is.na(rank_pq)] %in% true_values))
    }

    FN <- sum(!(true_values %in% rank_pq))

    TN <- sum(is.na(
      subset_taxa_pq(physeq, fake_taxa_cond, clean_pq = FALSE)@tax_table[, taxonomic_rank]
    ))

    FDR <- FP / (FP + TP)
    TPR <- TP / (TP + FN)
    TNR <- TN / (TN + FN)
    PPV <- TP / (TP + FP)
    # PPV+FDR == 1

    # F1_score <- 2 * (TPR * PPV) / TPR + PPV
    F1_score <- 2 * TP / (2 * TP + FP + FN)

    MCC <- (TP * TN - FP * FN) /
      sqrt((TP + FP) * (TP + FN) * (FP + TN) * (TN + FN))

    ACC <- (TP + TN) / (TP + TN + FP + FN)

    res <- list(
      "TP" = TP,
      "FP" = FP,
      "FN" = FN,
      "FDR" = FDR,
      "TPR" = TPR,
      "TNR" = TNR,
      "PPV" = PPV,
      "F1_score" = F1_score,
      "TN" = TN,
      "MCC" = MCC,
      "ACC" = ACC
    )
  } else {
    # Hleap et al.  2021 [https://doi.org/10.1111/1755-0998.13407]
    rank_pq <- physeq@tax_table[, taxonomic_rank]

    TP <- sum(rank_pq %in% true_values)
    FP <- sum(!(rank_pq[!is.na(rank_pq)] %in% true_values))
    FN <- sum(!(true_values %in% rank_pq))

    FDR <- FP / (FP + TP)
    TPR <- TP / (TP + FN)
    PPV <- TP / (TP + FP)
    # PPV+FDR == 1

    # F1_score <- 2 * (TPR * PPV) / TPR + PPV
    F1_score <- 2 * TP / (2 * TP + FP + FN)

    res <- list(
      "TP" = TP,
      "FP" = FP,
      "FN" = FN,
      "FDR" = FDR,
      "TPR" = TPR,
      "PPV" = PPV,
      "F1_score" = F1_score
    )
  }
  return(res)
}
################################################################################
