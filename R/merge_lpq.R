#' Merge a list_phyloseq into a single phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#'   alt="lifecycle-experimental"></a>
#'
#' Merges all phyloseq objects from a [list_phyloseq] into a single
#' phyloseq object where each original phyloseq becomes one sample.
#' Abundances are summed across samples within each phyloseq object.
#'
#' Taxa are matched across phyloseq objects either by reference
#' sequences (`refseq` slot, default) or by taxa names. Matching by
#' refseq is preferred because taxa names are often inconsistent
#' across independently built phyloseq objects (e.g., ASV_1 in one
#' object is not the same taxon as ASV_1 in another).
#'
#' @param x (list_phyloseq, required) A [list_phyloseq] object.
#' @param match_by (character, default `"refseq"`) How to match taxa
#'   across phyloseq objects. One of:
#'   - `"refseq"`: match by DNA sequence in the `refseq` slot
#'     (recommended). All phyloseq objects must have a `refseq` slot.
#'   - `"names"`: match by taxa names. Use only when taxa names are
#'     consistent across objects (e.g., same pipeline, same database).
#' @param tax_priority (character or integer, default 1L) Which
#'   phyloseq object's taxonomy to use when taxa are matched. Either
#'   a name from the list or an integer index. When a taxon appears
#'   in multiple objects, the taxonomy from the priority object is
#'   used; if absent there, the first available taxonomy is used.
#' @param verbose (logical, default TRUE) Print information about the
#'   merge.
#'
#' @return A phyloseq object with:
#' \describe{
#'   \item{`otu_table`}{One column per original phyloseq object
#'     (summed across its samples), one row per unique taxon.}
#'   \item{`sample_data`}{One row per original phyloseq, with a
#'     column `source_name` containing the list_phyloseq names.}
#'   \item{`tax_table`}{Taxonomy from the priority object (or first
#'     available).}
#'   \item{`refseq`}{Present when `match_by = "refseq"`.}
#' }
#'
#' @export
#' @author Adrien Taudière
#'
#' @seealso [list_phyloseq], [simple_venn_pq()]
#'
#' @examples
#' lpq <- list_phyloseq(list(
#'   fungi = data_fungi_mini,
#'   fungi2 = data_fungi_mini
#' ))
#'
#' # Merge by refseq (default)
#' merged <- merge_lpq(lpq)
#' merged
#'
#' # Merge by taxa names
#' merged_names <- merge_lpq(lpq, match_by = "names")
merge_lpq <- function(
  x,
  match_by = c("refseq", "names"),
  tax_priority = 1L,
  verbose = TRUE
) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))
  match_by <- match.arg(match_by)

  pq_list <- x@phyloseq_list
  pq_names <- names(pq_list)
  n_pq <- length(pq_list)

  if (match_by == "refseq") {
    all_have_refseq <- all(vapply(
      pq_list,
      \(pq) !is.null(phyloseq::refseq(pq, errorIfNULL = FALSE)),
      logical(1)
    ))
    if (!all_have_refseq) {
      stop(
        "All phyloseq objects must have a refseq slot when ",
        "match_by = \"refseq\". Use match_by = \"names\" instead ",
        "or add reference sequences to all objects.",
        call. = FALSE
      )
    }
  }

  # Resolve tax_priority to an index
  if (is.character(tax_priority)) {
    tax_priority <- match(tax_priority, pq_names)
    if (is.na(tax_priority)) {
      stop(
        "'tax_priority' name not found in list_phyloseq names.",
        call. = FALSE
      )
    }
  }

  if (match_by == "refseq") {
    result <- merge_by_refseq(pq_list, pq_names, tax_priority, verbose)
  } else {
    result <- merge_by_names(pq_list, pq_names, tax_priority, verbose)
  }

  result
}


#' Merge phyloseq objects by matching refseq sequences
#' @param pq_list Named list of phyloseq objects
#' @param pq_names Character vector of names
#' @param tax_priority Integer index for taxonomy priority
#' @param verbose Logical
#' @return A phyloseq object
#' @noRd
merge_by_refseq <- function(pq_list, pq_names, tax_priority, verbose) {
  n_pq <- length(pq_list)

  # Collect all unique sequences across all objects
  all_seqs_list <- lapply(pq_list, \(pq) {
    seqs <- phyloseq::refseq(pq)
    stats::setNames(as.character(seqs), phyloseq::taxa_names(pq))
  })

  all_unique_seqs <- unique(unlist(all_seqs_list, use.names = FALSE))
  n_taxa_merged <- length(all_unique_seqs)

  if (verbose) {
    n_per_pq <- vapply(all_seqs_list, length, integer(1))
    message(
      "Merging ",
      n_pq,
      " phyloseq objects by refseq: ",
      paste(n_per_pq, collapse = " + "),
      " taxa -> ",
      n_taxa_merged,
      " unique sequences."
    )
  }

  # Build merged OTU table: rows = unique sequences, columns = phyloseq names
  otu_mat <- matrix(
    0L,
    nrow = n_taxa_merged,
    ncol = n_pq,
    dimnames = list(
      paste0("seq_", seq_len(n_taxa_merged)),
      pq_names
    )
  )

  for (i in seq_len(n_pq)) {
    pq <- pq_list[[i]]
    otu <- as.data.frame(phyloseq::otu_table(pq))
    if (!phyloseq::taxa_are_rows(pq)) {
      otu <- as.data.frame(t(otu))
    }
    # Sum across samples for this phyloseq
    taxa_sums <- rowSums(otu)
    # Map each taxon to its sequence position in all_unique_seqs
    seqs <- all_seqs_list[[i]]
    seq_idx <- match(seqs, all_unique_seqs)
    otu_mat[seq_idx, i] <- as.integer(taxa_sums)
  }

  # Build taxonomy table using priority object
  # For each unique sequence, find taxonomy from priority object first,
  # then fall back to first available
  tax_tables <- lapply(pq_list, \(pq) {
    tt <- phyloseq::tax_table(pq, errorIfNULL = FALSE)
    if (is.null(tt)) {
      return(NULL)
    }
    as.data.frame(tt)
  })

  # Determine common tax ranks
  all_ranks <- unique(unlist(lapply(
    Filter(Negate(is.null), tax_tables),
    colnames
  )))

  tax_mat <- matrix(
    NA_character_,
    nrow = n_taxa_merged,
    ncol = length(all_ranks),
    dimnames = list(rownames(otu_mat), all_ranks)
  )

  for (seq_i in seq_len(n_taxa_merged)) {
    seq_val <- all_unique_seqs[seq_i]
    filled <- FALSE

    # Try priority object first
    if (!is.null(tax_tables[[tax_priority]])) {
      seqs_priority <- all_seqs_list[[tax_priority]]
      taxa_match <- names(seqs_priority)[seqs_priority == seq_val]
      if (length(taxa_match) > 0) {
        tt <- tax_tables[[tax_priority]]
        if (taxa_match[1] %in% rownames(tt)) {
          available_ranks <- intersect(all_ranks, colnames(tt))
          tax_mat[seq_i, available_ranks] <- as.character(
            tt[taxa_match[1], available_ranks]
          )
          filled <- TRUE
        }
      }
    }

    # Fall back to first available
    if (!filled) {
      for (j in seq_len(n_pq)) {
        if (is.null(tax_tables[[j]])) {
          next
        }
        seqs_j <- all_seqs_list[[j]]
        taxa_match <- names(seqs_j)[seqs_j == seq_val]
        if (length(taxa_match) > 0) {
          tt <- tax_tables[[j]]
          if (taxa_match[1] %in% rownames(tt)) {
            available_ranks <- intersect(all_ranks, colnames(tt))
            tax_mat[seq_i, available_ranks] <- as.character(
              tt[taxa_match[1], available_ranks]
            )
            break
          }
        }
      }
    }
  }

  # Build refseq
  refseq <- Biostrings::DNAStringSet(all_unique_seqs)
  names(refseq) <- rownames(otu_mat)

  # Build sample_data
  sam_df <- data.frame(
    source_name = pq_names,
    row.names = pq_names,
    stringsAsFactors = FALSE
  )

  # Assemble phyloseq
  otu_obj <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  tax_obj <- phyloseq::tax_table(as.matrix(tax_mat))
  sam_obj <- phyloseq::sample_data(sam_df)

  phyloseq::phyloseq(otu_obj, tax_obj, sam_obj, refseq)
}


#' Merge phyloseq objects by matching taxa names
#' @param pq_list Named list of phyloseq objects
#' @param pq_names Character vector of names
#' @param tax_priority Integer index for taxonomy priority
#' @param verbose Logical
#' @return A phyloseq object
#' @noRd
merge_by_names <- function(pq_list, pq_names, tax_priority, verbose) {
  n_pq <- length(pq_list)

  # Collect all unique taxa names
  all_taxa_names <- unique(unlist(
    lapply(pq_list, phyloseq::taxa_names)
  ))
  n_taxa_merged <- length(all_taxa_names)

  if (verbose) {
    n_per_pq <- vapply(pq_list, phyloseq::ntaxa, integer(1))
    message(
      "Merging ",
      n_pq,
      " phyloseq objects by names: ",
      paste(n_per_pq, collapse = " + "),
      " taxa -> ",
      n_taxa_merged,
      " unique names."
    )
  }

  # Build merged OTU table
  otu_mat <- matrix(
    0L,
    nrow = n_taxa_merged,
    ncol = n_pq,
    dimnames = list(all_taxa_names, pq_names)
  )

  for (i in seq_len(n_pq)) {
    pq <- pq_list[[i]]
    otu <- as.data.frame(phyloseq::otu_table(pq))
    if (!phyloseq::taxa_are_rows(pq)) {
      otu <- as.data.frame(t(otu))
    }
    taxa_sums <- rowSums(otu)
    taxa_in_pq <- phyloseq::taxa_names(pq)
    otu_mat[taxa_in_pq, i] <- as.integer(taxa_sums[taxa_in_pq])
  }

  # Build taxonomy table
  tax_tables <- lapply(pq_list, \(pq) {
    tt <- phyloseq::tax_table(pq, errorIfNULL = FALSE)
    if (is.null(tt)) {
      return(NULL)
    }
    as.data.frame(tt)
  })

  all_ranks <- unique(unlist(lapply(
    Filter(Negate(is.null), tax_tables),
    colnames
  )))

  tax_mat <- matrix(
    NA_character_,
    nrow = n_taxa_merged,
    ncol = length(all_ranks),
    dimnames = list(all_taxa_names, all_ranks)
  )

  for (taxon in all_taxa_names) {
    filled <- FALSE

    # Try priority object first
    if (!is.null(tax_tables[[tax_priority]])) {
      tt <- tax_tables[[tax_priority]]
      if (taxon %in% rownames(tt)) {
        available_ranks <- intersect(all_ranks, colnames(tt))
        tax_mat[taxon, available_ranks] <- as.character(
          tt[taxon, available_ranks]
        )
        filled <- TRUE
      }
    }

    if (!filled) {
      for (j in seq_len(n_pq)) {
        if (is.null(tax_tables[[j]])) {
          next
        }
        tt <- tax_tables[[j]]
        if (taxon %in% rownames(tt)) {
          available_ranks <- intersect(all_ranks, colnames(tt))
          tax_mat[taxon, available_ranks] <- as.character(
            tt[taxon, available_ranks]
          )
          break
        }
      }
    }
  }

  # Build refseq if all objects have one
  refseq <- NULL
  all_have_refseq <- all(vapply(
    pq_list,
    \(pq) !is.null(phyloseq::refseq(pq, errorIfNULL = FALSE)),
    logical(1)
  ))
  if (all_have_refseq) {
    all_seqs <- character(n_taxa_merged)
    names(all_seqs) <- all_taxa_names
    # Priority object first, then fill gaps
    for (j in c(tax_priority, setdiff(seq_len(n_pq), tax_priority))) {
      seqs_j <- as.character(phyloseq::refseq(pq_list[[j]]))
      names_j <- phyloseq::taxa_names(pq_list[[j]])
      missing <- all_seqs[all_taxa_names] == ""
      to_fill <- intersect(
        names_j[names_j %in% all_taxa_names],
        names(missing[missing])
      )
      if (length(to_fill) > 0) {
        idx <- match(to_fill, names_j)
        all_seqs[to_fill] <- seqs_j[idx]
      }
    }
    if (all(all_seqs != "")) {
      refseq <- Biostrings::DNAStringSet(all_seqs)
    }
  }

  # Build sample_data
  sam_df <- data.frame(
    source_name = pq_names,
    row.names = pq_names,
    stringsAsFactors = FALSE
  )

  # Assemble phyloseq
  otu_obj <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  tax_obj <- phyloseq::tax_table(as.matrix(tax_mat))
  sam_obj <- phyloseq::sample_data(sam_df)

  if (!is.null(refseq)) {
    phyloseq::phyloseq(otu_obj, tax_obj, sam_obj, refseq)
  } else {
    phyloseq::phyloseq(otu_obj, tax_obj, sam_obj)
  }
}
