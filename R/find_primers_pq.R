#' Find taxa whose reference sequences match primer sequences
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#'   alt="lifecycle-experimental"></a>
#'
#' @description
#' Searches every sequence in `@refseq` for occurrences of the supplied
#' primers (forward **and** reverse complement) using IUPAC-aware matching via
#' [Biostrings::vcountPattern()]. Returns a data frame of taxa that match at
#' least one primer, which can be passed directly to
#' `tidypq::filter_taxa_pq()` to prune them from the phyloseq object.
#'
#' @param physeq (required) A [phyloseq::phyloseq] object with a populated
#'   `@refseq` slot.
#' @param primers (required) A named character vector of primer sequences.
#'   IUPAC ambiguity codes (M, R, Y, S, W, K, B, D, H, V, N) are supported.
#' @param verbose (logical, default `TRUE`) If `TRUE`, print a summary message.
#'
#' @return A `data.frame` (or `NULL` if no matches) with columns:
#'   \describe{
#'     \item{`taxon`}{Character. Taxa name as in `taxa_names(physeq)`.}
#'     \item{`matched_primers`}{Character. Comma-separated names of matching
#'       primers.}
#'     \item{`n_reads`}{Numeric. Total read count across all samples
#'       (`phyloseq::taxa_sums(physeq)`).}
#'   }
#'
#' @importFrom Biostrings DNAString reverseComplement vcountPattern
#'
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' primers <- c(
#'   mcrA_fwd = "GGTGGTGTMGGDTTCACMCARTA",
#'   mcrA_rev = "CGTTCATBGCGTAGTTVGGRTAGT"
#' )
#' bad <- find_primers_pq(data_fungi_mini, primers)
#' bad
#'
#' # Prune contaminated taxa (requires tidypq):
#' # if (!is.null(bad)) {
#' #   tidypq::filter_taxa_pq(
#' #     data_fungi_mini,
#' #     !taxa_names(data_fungi_mini) %in% bad$taxon,
#' #     clean_phyloseq_object = TRUE
#' #   )
#' # }
find_primers_pq <- function(physeq, primers, verbose = TRUE) {
  refseq <- physeq@refseq
  if (is.null(refseq)) {
    stop("physeq has no @refseq slot.")
  }
  if (length(primers) == 0L) {
    stop("primers must be a non-empty named character vector.")
  }

  match_mat <- vapply(
    names(primers),
    function(p_name) {
      p <- DNAString(primers[[p_name]])
      rc_p <- reverseComplement(p)
      fwd <- vcountPattern(p, refseq, fixed = FALSE) > 0L
      rev <- vcountPattern(rc_p, refseq, fixed = FALSE) > 0L
      fwd | rev
    },
    logical(length(refseq))
  )

  any_match <- rowSums(match_mat) > 0L
  n_hits <- sum(any_match)

  if (verbose) {
    message(
      n_hits,
      " taxa matched at least one primer out of ",
      length(refseq),
      "."
    )
  }

  if (!any(any_match)) {
    return(NULL)
  }

  bad_idx <- which(any_match)
  data.frame(
    taxon = names(refseq)[bad_idx],
    matched_primers = apply(
      match_mat[bad_idx, , drop = FALSE],
      1,
      \(row) paste(names(primers)[row], collapse = ", ")
    ),
    n_reads = phyloseq::taxa_sums(physeq)[bad_idx],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
