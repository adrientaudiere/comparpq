#' Merge a list_phyloseq into one phyloseq by clustering reference sequences
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#'   alt="lifecycle-experimental"></a>
#'
#' Merges all phyloseq objects from a [list_phyloseq] into a **single**
#' phyloseq object while **keeping every sample separate**. Taxa are
#' unified *across* objects by clustering their reference sequences
#' (`refseq` slot) with a clustering algorithm (vsearch by default), so
#' that sequences that are similar enough (controlled by `id`) become a
#' single merged taxon whose counts are summed within each original
#' sample.
#'
#' This differs from [merge_lpq()], which collapses each phyloseq object
#' into a single sample and matches taxa by **exact** sequence identity.
#' Use `merge_clust_lpq()` when you want to pool the raw samples of
#' several independently built datasets into one object and reconcile
#' their ASVs/OTUs through fuzzy (identity-threshold) clustering.
#'
#' Samples keep their original names; when the same sample name occurs in
#' more than one object, the parent object name is appended as a suffix
#' (`<sample>_<object>`) to keep sample names unique. A column
#' (`source_name` by default) recording the parent object is added to the
#' `sample_data`.
#'
#' @param x (list_phyloseq or list, required) A [list_phyloseq] object, or
#'   a (preferably named) list of phyloseq objects, all of which must have
#'   a `refseq` slot.
#' @param method (character, default `"vsearch"`) Clustering method passed
#'   to [MiscMetabar::postcluster_pq()]. One of `"clusterize"`,
#'   `"vsearch"`, `"swarm"` or `"mmseqs2"`.
#' @param id (numeric, default 0.97) Sequence identity threshold for the
#'   clustering (used by the `vsearch`/`mmseqs2` methods). Ignored by
#'   `swarm`.
#' @param source_col (character, default `"source_name"`) Name of the
#'   `sample_data` column that will store the parent object name for each
#'   sample.
#' @param verbose (logical, default TRUE) Print information about the
#'   merge.
#' @param ... Further arguments passed to [MiscMetabar::postcluster_pq()]
#'   (e.g. `tax_adjust`, `rank_propagation`, `nproc`).
#'
#' @return A phyloseq object with:
#' \describe{
#'   \item{`otu_table`}{One column per original sample (across all
#'     objects, suffixed on name collision), one row per clustered taxon.}
#'   \item{`sample_data`}{One row per original sample, with a `source_col`
#'     column giving the parent object name.}
#'   \item{`tax_table`}{Taxonomy carried from the clustering
#'     representative of each merged taxon.}
#'   \item{`refseq`}{Representative sequence of each clustered taxon.}
#' }
#'
#' @export
#' @author Adrien Taudière
#'
#' @seealso [merge_lpq()], [list_phyloseq], [MiscMetabar::postcluster_pq()]
#'
#' @examples
#' \dontrun{
#' # Requires vsearch to be installed by default (MiscMetabar::install_vsearch())
#' library(MiscMetabar)
#' pq1 <- postcluster_pq(data_fungi_mini, method = "vsearch", id = 0.97)
#' pq2 <- clean_pq(prune_samples(sample_names(pq1)[1:4], pq1))
#'
#' lpq <- list_phyloseq(list(run_a = pq2, run_b = pq1))
#'
#' merged <- merge_clust_lpq(lpq, id = 0.97)
#' merged
#' table(sample_data(merged)$source_name)
#' }
merge_clust_lpq <- function(
  x,
  method = "vsearch",
  id = 0.97,
  source_col = "source_name",
  verbose = TRUE,
  ...
) {
  if (inherits(x, "comparpq::list_phyloseq")) {
    pq_list <- x@phyloseq_list
  } else if (
    is.list(x) && all(vapply(x, \(p) inherits(p, "phyloseq"), logical(1)))
  ) {
    pq_list <- x
  } else {
    cli::cli_abort(
      "{.arg x} must be a {.cls list_phyloseq} or a list of phyloseq objects."
    )
  }

  if (length(pq_list) < 2) {
    cli::cli_abort("{.arg x} must contain at least two phyloseq objects.")
  }

  pq_names <- names(pq_list)
  if (is.null(pq_names) || any(pq_names == "")) {
    pq_names <- paste0("physeq_", seq_along(pq_list))
    names(pq_list) <- pq_names
  }

  all_have_refseq <- all(vapply(
    pq_list,
    \(pq) !is.null(phyloseq::refseq(pq, errorIfNULL = FALSE)),
    logical(1)
  ))
  if (!all_have_refseq) {
    cli::cli_abort(
      "All phyloseq objects must have a {.field refseq} slot to be clustered."
    )
  }

  union_pq <- build_union_pq(pq_list, pq_names, source_col, verbose)

  n_before <- phyloseq::ntaxa(union_pq)
  merged <- MiscMetabar::postcluster_pq(
    union_pq,
    method = method,
    id = id,
    ...
  )

  if (verbose) {
    message(
      "Clustered ",
      n_before,
      " pooled taxa (from ",
      length(pq_list),
      " objects) into ",
      phyloseq::ntaxa(merged),
      " taxa across ",
      phyloseq::nsamples(merged),
      " samples."
    )
  }

  merged
}


#' Build a naive union phyloseq keeping all samples, deduplicating by
#' exact sequence
#'
#' All samples of all objects are kept (suffixed on name collision).
#' Taxa sharing an **identical** reference sequence across objects are
#' merged into a single row (counts summed within their own samples), so
#' the resulting `refseq` slot has no exact duplicates and can be fed to
#' [MiscMetabar::postcluster_pq()] for fuzzy (identity-threshold)
#' clustering of the remaining near-duplicates.
#'
#' @param pq_list Named list of phyloseq objects
#' @param pq_names Character vector of names
#' @param source_col Name of the sample_data column storing the parent name
#' @param verbose Logical
#' @return A phyloseq object
#' @noRd
build_union_pq <- function(pq_list, pq_names, source_col, verbose) {
  n_pq <- length(pq_list)

  all_sample_names <- unlist(lapply(pq_list, phyloseq::sample_names))
  dup_names <- unique(all_sample_names[duplicated(all_sample_names)])

  otu_list <- vector("list", n_pq)
  seq_of_taxon <- list()
  tax_list <- vector("list", n_pq)
  sam_list <- vector("list", n_pq)

  for (i in seq_len(n_pq)) {
    pq <- pq_list[[i]]
    nm <- pq_names[i]

    otu <- as(phyloseq::otu_table(pq), "matrix")
    if (!phyloseq::taxa_are_rows(pq)) {
      otu <- t(otu)
    }
    storage.mode(otu) <- "double"

    smp <- colnames(otu)
    new_smp <- ifelse(smp %in% dup_names, paste0(smp, "_", nm), smp)
    colnames(otu) <- new_smp
    otu_list[[i]] <- otu

    rs <- as.character(phyloseq::refseq(pq))
    seq_of_taxon[[i]] <- rs[rownames(otu)]

    tt <- phyloseq::tax_table(pq, errorIfNULL = FALSE)
    if (!is.null(tt)) {
      tax_list[[i]] <- as.data.frame(tt, stringsAsFactors = FALSE)
    }

    sd <- phyloseq::sample_data(pq, errorIfNULL = FALSE)
    if (is.null(sd)) {
      sd <- data.frame(row.names = new_smp)
    } else {
      sd <- data.frame(sd, check.names = FALSE, stringsAsFactors = FALSE)
      rownames(sd) <- new_smp
    }
    sd[[source_col]] <- nm
    sam_list[[i]] <- sd
  }

  # Unique sequences become the merged taxa (exact-match deduplication).
  # Zero-pad ids so alphabetical order matches creation order (keeps the
  # otu_table / tax_table / refseq slots consistently ordered).
  unique_seqs <- unique(unlist(seq_of_taxon, use.names = FALSE))
  seq_ids <- sprintf(
    paste0("union_%0", nchar(length(unique_seqs)), "d"),
    seq_along(unique_seqs)
  )
  names(seq_ids) <- unique_seqs

  global_samples <- unlist(lapply(otu_list, colnames), use.names = FALSE)

  otu_mat <- matrix(
    0,
    nrow = length(unique_seqs),
    ncol = length(global_samples),
    dimnames = list(unname(seq_ids), global_samples)
  )
  for (i in seq_len(n_pq)) {
    rows <- seq_ids[seq_of_taxon[[i]]]
    block <- otu_list[[i]]
    # Sum in case several taxa within one object share a sequence
    agg <- rowsum(block, group = rows, reorder = FALSE)
    otu_mat[rownames(agg), colnames(block)] <-
      otu_mat[rownames(agg), colnames(block)] + agg
  }

  all_ranks <- unique(unlist(lapply(
    Filter(Negate(is.null), tax_list),
    colnames
  )))
  tax_mat <- matrix(
    NA_character_,
    nrow = length(unique_seqs),
    ncol = length(all_ranks),
    dimnames = list(unname(seq_ids), all_ranks)
  )
  # First object providing a taxonomy for a given sequence wins
  for (i in seq_len(n_pq)) {
    tt <- tax_list[[i]]
    if (is.null(tt)) {
      next
    }
    taxa_i <- names(seq_of_taxon[[i]])
    ids_i <- unname(seq_ids[seq_of_taxon[[i]]])
    ranks_i <- intersect(all_ranks, colnames(tt))
    for (k in seq_along(taxa_i)) {
      id <- ids_i[k]
      if (all(is.na(tax_mat[id, ranks_i]))) {
        tax_mat[id, ranks_i] <- as.character(tt[taxa_i[k], ranks_i])
      }
    }
  }

  sam_df <- dplyr::bind_rows(lapply(sam_list, function(d) {
    d$.rowname <- rownames(d)
    d
  }))
  sam_df <- as.data.frame(sam_df, stringsAsFactors = FALSE)
  rownames(sam_df) <- sam_df$.rowname
  sam_df$.rowname <- NULL

  refseq <- Biostrings::DNAStringSet(stats::setNames(unique_seqs, seq_ids))

  if (verbose) {
    n_per_pq <- vapply(otu_list, nrow, integer(1))
    message(
      "Pooling ",
      n_pq,
      " phyloseq objects: ",
      paste(n_per_pq, collapse = " + "),
      " taxa -> ",
      length(unique_seqs),
      " unique sequences across ",
      length(global_samples),
      " samples before clustering."
    )
  }

  phyloseq::phyloseq(
    phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(tax_mat)),
    phyloseq::sample_data(sam_df),
    refseq
  )
}
