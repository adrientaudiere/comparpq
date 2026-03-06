#' Compare refseq slots between two phyloseq objects
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#'   alt="lifecycle-experimental"></a>
#'
#' Performs a detailed comparison of the reference sequences (`refseq` slot)
#' between two phyloseq objects. Identifies shared and unique ASVs/OTUs both
#' by taxa name and by actual DNA sequence content.
#'
#' This is useful to detect subtle differences between two phyloseq objects
#' that may share the same samples but differ in their ASV/OTU composition,
#' e.g. after different bioinformatics pipelines or filtering steps.
#'
#' When unique (non-shared) sequences exist, the mean nearest-neighbor k-mer
#' distance from each unique sequence to the other set is computed. This
#' gives a sense of how different the unmatched sequences are from their
#' closest counterpart in the other object.
#'
#' @param physeq1 (phyloseq or list_phyloseq, required) First phyloseq
#'   object, which must have a `refseq` slot. Alternatively, a
#'   [list_phyloseq] object; in that case the first two phyloseq objects
#'   are used (with their names) and `physeq2` is ignored.
#' @param physeq2 (phyloseq, default NULL) Second phyloseq object. Must
#'   have a `refseq` slot. Ignored when `physeq1` is a [list_phyloseq].
#' @param name1 (character, default NULL) Label for the first phyloseq
#'   object in the output. If NULL, inferred from the [list_phyloseq]
#'   names or defaults to `"physeq1"`.
#' @param name2 (character, default NULL) Label for the second phyloseq
#'   object in the output. If NULL, inferred from the [list_phyloseq]
#'   names or defaults to `"physeq2"`.
#' @param k (integer, default 5) k-mer size for nearest-neighbor distance
#'   computation on unique sequences.
#' @param max_seqs (integer, default 500) Maximum number of unique sequences
#'   per object to use for distance computation. If there are more, a random
#'   sample is drawn. Set to `Inf` to use all (may be slow).
#' @param seed (integer, default NULL) Random seed for reproducible sampling
#'   when `max_seqs` is exceeded.
#' @param verbose (logical, default TRUE) If TRUE, print a summary of the
#'   comparison.
#'
#' @return A list of class `"compare_refseq"` with components:
#' \describe{
#'   \item{`n_taxa`}{Named integer vector with the number of taxa in each
#'     object.}
#'   \item{`shared_names`}{Character vector of taxa names present in both
#'     objects.}
#'   \item{`unique_names_1`}{Taxa names only in `physeq1`.}
#'   \item{`unique_names_2`}{Taxa names only in `physeq2`.}
#'   \item{`shared_seqs`}{Character vector of unique DNA sequences found in
#'     both objects (regardless of taxa name).}
#'   \item{`unique_seqs_1`}{DNA sequences only in `physeq1`.}
#'   \item{`unique_seqs_2`}{DNA sequences only in `physeq2`.}
#'   \item{`same_name_diff_seq`}{A data frame of taxa with the same name but
#'     different sequences across objects. Columns: `taxa_name`,
#'     `seq_physeq1`, `seq_physeq2`.}
#'   \item{`same_seq_diff_name`}{A data frame of sequences shared between
#'     objects under different names. Columns: `sequence`, `name_physeq1`,
#'     `name_physeq2`.}
#'   \item{`mean_nn_kmer_dist_1`}{Mean nearest-neighbor k-mer distance from
#'     sequences unique to `physeq1` to the closest sequence in `physeq2`.
#'     `NA` if no unique sequences.}
#'   \item{`mean_nn_kmer_dist_2`}{Mean nearest-neighbor k-mer distance from
#'     sequences unique to `physeq2` to the closest sequence in `physeq1`.
#'     `NA` if no unique sequences.}
#' }
#'
#' @export
#' @author Adrien Taudiere
#'
#' @examples
#' # Compare a phyloseq object with itself (all shared)
#' res <- compare_refseq(data_fungi_mini, data_fungi_mini)
#' res
#'
#' # Compare with a subset
#' sub <- prune_taxa(taxa_names(data_fungi_mini)[1:20], data_fungi_mini)
#' res2 <- compare_refseq(data_fungi_mini, sub,
#'   name1 = "full", name2 = "subset"
#' )
#' res2
#'
#' # From a list_phyloseq (uses the first two objects and their names)
#' lpq <- list_phyloseq(list(full = data_fungi_mini, subset = sub))
#' res3 <- compare_refseq(lpq)
#' res3
compare_refseq <- function(
  physeq1,
  physeq2 = NULL,
  name1 = NULL,
  name2 = NULL,
  k = 5,
  max_seqs = 500,
  seed = NULL,
  verbose = TRUE
) {
  # --- list_phyloseq dispatch ---
  if (inherits(physeq1, "comparpq::list_phyloseq")) {
    lpq <- physeq1
    if (length(lpq) < 2) {
      stop("list_phyloseq must contain at least 2 phyloseq objects.")
    }
    if (length(lpq) > 2) {
      message(
        "list_phyloseq contains ",
        length(lpq),
        " objects. Using the first two: '",
        names(lpq)[1],
        "' and '",
        names(lpq)[2],
        "'."
      )
    }
    physeq1 <- lpq[[1]]
    physeq2 <- lpq[[2]]
    if (is.null(name1)) {
      name1 <- names(lpq)[1]
    }
    if (is.null(name2)) {
      name2 <- names(lpq)[2]
    }
  }

  if (is.null(physeq2)) {
    stop("'physeq2' is required when 'physeq1' is a phyloseq object.")
  }
  if (is.null(name1)) {
    name1 <- "physeq1"
  }
  if (is.null(name2)) {
    name2 <- "physeq2"
  }

  rs1 <- phyloseq::refseq(physeq1, errorIfNULL = FALSE)
  rs2 <- phyloseq::refseq(physeq2, errorIfNULL = FALSE)

  if (is.null(rs1)) {
    stop("'", name1, "' has no refseq slot.")
  }
  if (is.null(rs2)) {
    stop("'", name2, "' has no refseq slot.")
  }

  names1 <- names(rs1)
  names2 <- names(rs2)
  seqs1 <- as.character(rs1)
  seqs2 <- as.character(rs2)

  # --- By name ---
  shared_names <- intersect(names1, names2)
  unique_names_1 <- setdiff(names1, names2)
  unique_names_2 <- setdiff(names2, names1)

  # --- By sequence ---
  seqs_set1 <- unique(seqs1)
  seqs_set2 <- unique(seqs2)
  shared_seqs <- intersect(seqs_set1, seqs_set2)
  unique_seqs_1 <- setdiff(seqs_set1, seqs_set2)
  unique_seqs_2 <- setdiff(seqs_set2, seqs_set1)

  # --- Same name, different sequence ---
  if (length(shared_names) > 0) {
    seq_in_1 <- seqs1[shared_names]
    seq_in_2 <- seqs2[shared_names]
    differ <- seq_in_1 != seq_in_2
    if (any(differ)) {
      diff_names <- shared_names[differ]
      same_name_diff_seq <- data.frame(
        taxa_name = diff_names,
        seq_physeq1 = unname(seq_in_1[diff_names]),
        seq_physeq2 = unname(seq_in_2[diff_names]),
        stringsAsFactors = FALSE
      )
    } else {
      same_name_diff_seq <- data.frame(
        taxa_name = character(0),
        seq_physeq1 = character(0),
        seq_physeq2 = character(0),
        stringsAsFactors = FALSE
      )
    }
  } else {
    same_name_diff_seq <- data.frame(
      taxa_name = character(0),
      seq_physeq1 = character(0),
      seq_physeq2 = character(0),
      stringsAsFactors = FALSE
    )
  }

  # --- Same sequence, different name ---
  seq_to_names1 <- split(names1, seqs1)
  seq_to_names2 <- split(names2, seqs2)

  same_seq_diff_name_rows <- list()
  for (seq in shared_seqs) {
    n1 <- seq_to_names1[[seq]]
    n2 <- seq_to_names2[[seq]]
    only_n1 <- setdiff(n1, n2)
    only_n2 <- setdiff(n2, n1)
    if (length(only_n1) > 0 || length(only_n2) > 0) {
      same_seq_diff_name_rows <- c(
        same_seq_diff_name_rows,
        list(data.frame(
          sequence = seq,
          name_physeq1 = paste(n1, collapse = ", "),
          name_physeq2 = paste(n2, collapse = ", "),
          stringsAsFactors = FALSE
        ))
      )
    }
  }

  if (length(same_seq_diff_name_rows) > 0) {
    same_seq_diff_name <- do.call(rbind, same_seq_diff_name_rows)
    rownames(same_seq_diff_name) <- NULL
  } else {
    same_seq_diff_name <- data.frame(
      sequence = character(0),
      name_physeq1 = character(0),
      name_physeq2 = character(0),
      stringsAsFactors = FALSE
    )
  }

  # --- Mean nearest-neighbor k-mer distance for unique sequences ---
  mean_nn_kmer_dist_1 <- NA_real_
  mean_nn_kmer_dist_2 <- NA_real_

  if (length(unique_seqs_1) > 0) {
    mean_nn_kmer_dist_1 <- compute_mean_nn_kmer(
      unique_seqs_1,
      rs2,
      k,
      max_seqs,
      seed
    )
  }
  if (length(unique_seqs_2) > 0) {
    mean_nn_kmer_dist_2 <- compute_mean_nn_kmer(
      unique_seqs_2,
      rs1,
      k,
      max_seqs,
      seed
    )
  }

  result <- structure(
    list(
      name1 = name1,
      name2 = name2,
      k = k,
      n_taxa = stats::setNames(
        c(length(names1), length(names2)),
        c(name1, name2)
      ),
      shared_names = shared_names,
      unique_names_1 = unique_names_1,
      unique_names_2 = unique_names_2,
      shared_seqs = shared_seqs,
      unique_seqs_1 = unique_seqs_1,
      unique_seqs_2 = unique_seqs_2,
      same_name_diff_seq = same_name_diff_seq,
      same_seq_diff_name = same_seq_diff_name,
      mean_nn_kmer_dist_1 = mean_nn_kmer_dist_1,
      mean_nn_kmer_dist_2 = mean_nn_kmer_dist_2
    ),
    class = "compare_refseq"
  )

  if (verbose) {
    print(result)
  }

  invisible(result)
}


#' @export
print.compare_refseq <- function(x, ...) {
  cat("== Reference Sequence Comparison ==\n")
  cat(x$name1, ":", x$n_taxa[[1]], "taxa\n")
  cat(x$name2, ":", x$n_taxa[[2]], "taxa\n\n")

  cat("-- By taxa name --\n")
  cat("  Shared :", length(x$shared_names), "\n")
  cat("  Only in", x$name1, ":", length(x$unique_names_1), "\n")
  cat("  Only in", x$name2, ":", length(x$unique_names_2), "\n\n")

  cat("-- By DNA sequence --\n")
  cat("  Shared :", length(x$shared_seqs), "\n")
  cat("  Only in", x$name1, ":", length(x$unique_seqs_1), "\n")
  cat("  Only in", x$name2, ":", length(x$unique_seqs_2), "\n\n")

  if (nrow(x$same_name_diff_seq) > 0) {
    cat(
      "  ! ",
      nrow(x$same_name_diff_seq),
      " taxa share a name but have DIFFERENT sequences\n"
    )
  }
  if (nrow(x$same_seq_diff_name) > 0) {
    cat(
      "  ! ",
      nrow(x$same_seq_diff_name),
      " sequences are shared but under DIFFERENT names\n"
    )
  }

  has_dist <- !is.na(x$mean_nn_kmer_dist_1) || !is.na(x$mean_nn_kmer_dist_2)
  if (has_dist) {
    cat(
      "\n-- Mean nearest-neighbor k-mer distance (k=",
      x$k,
      ") --\n",
      sep = ""
    )
    if (!is.na(x$mean_nn_kmer_dist_1)) {
      cat(
        "  Unique in",
        x$name1,
        "->",
        x$name2,
        ":",
        round(x$mean_nn_kmer_dist_1, 4),
        "\n"
      )
    }
    if (!is.na(x$mean_nn_kmer_dist_2)) {
      cat(
        "  Unique in",
        x$name2,
        "->",
        x$name1,
        ":",
        round(x$mean_nn_kmer_dist_2, 4),
        "\n"
      )
    }
  }

  invisible(x)
}


#' Compute mean nearest-neighbor k-mer distance from query sequences to a
#' target DNAStringSet
#'
#' @param query_seqs Character vector of DNA sequences.
#' @param target_refseq A DNAStringSet (the target set to search against).
#' @param k Integer, k-mer size.
#' @param max_seqs Integer, max query sequences to use.
#' @param seed Integer or NULL, random seed for sampling.
#' @return Numeric scalar: mean of minimum k-mer distances.
#' @noRd
compute_mean_nn_kmer <- function(query_seqs, target_refseq, k, max_seqs, seed) {
  if (length(query_seqs) > max_seqs) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    query_seqs <- sample(query_seqs, max_seqs)
  }

  query_dss <- Biostrings::DNAStringSet(query_seqs)
  names(query_dss) <- paste0("q_", seq_along(query_dss))
  combined <- c(query_dss, target_refseq)

  dnabin <- ape::as.DNAbin(combined)
  dist_mat <- as.matrix(kmer::kdistance(dnabin, k = k))

  n_query <- length(query_dss)
  target_idx <- n_query + seq_len(length(target_refseq))

  nn_dists <- apply(
    dist_mat[seq_len(n_query), target_idx, drop = FALSE],
    1,
    min
  )
  mean(nn_dists)
}
