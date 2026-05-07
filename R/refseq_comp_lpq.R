uf_components <- function(n_nodes, edges) {
  parent <- seq_len(n_nodes)

  find_root <- function(x) {
    while (parent[x] != x) {
      x <- parent[x]
    }
    x
  }

  if (!is.null(edges) && nrow(edges) > 0L) {
    for (i in seq_len(nrow(edges))) {
      ra <- find_root(edges[i, 1L])
      rb <- find_root(edges[i, 2L])
      if (ra != rb) {
        parent[max(ra, rb)] <- min(ra, rb)
      }
    }
  }

  vapply(seq_len(n_nodes), find_root, integer(1L))
}


#' Compare refseq sequences across all objects in a list_phyloseq
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#'   alt="lifecycle-experimental"></a>
#'
#' @description
#' For each pair of phyloseq objects in `x`, computes pairwise sequence
#' similarity between `@refseq` slots using k-mer Jaccard similarity, then
#' groups near-identical sequences across all objects via connected components
#' (union-find). Returns a list of results per threshold, and one Venn
#' diagram per threshold (requires \pkg{ggVennDiagram}).
#'
#' Similarity is approximated as `p ≈ (2j / (1 + j))^(1/k)` where `j` is the
#' Jaccard index of shared k-mer presence and `k` is the k-mer size. The
#' approximation is accurate within ~2 percentage points for identities in the
#' 90–100% range.
#'
#' @param x (required) A [list_phyloseq] object. All phyloseq objects must have
#'   a populated `@refseq` slot. Supports 2–7 objects (Venn diagrams are
#'   suppressed for > 7).
#' @param thresholds (numeric vector, default `c(0, 1, 2, 3, 5, 10)`) Maximum
#'   dissimilarity thresholds in percent. `0` means exact string matches only;
#'   `10` means sequences within 10% divergence are grouped together.
#' @param k (integer, default `7`) K-mer size for building presence/absence
#'   profiles. Larger `k` increases specificity; smaller `k` increases
#'   sensitivity for divergent sequences.
#' @param verbose (logical, default `TRUE`) If `TRUE`, print progress messages.
#'
#' @return A list of class `"refseq_comp_lpq_result"` with:
#'   \describe{
#'     \item{`venn_plots`}{Named list of ggVennDiagram objects (one per
#'       threshold, names follow `"<threshold>_pct"`), or `NULL` if
#'       \pkg{ggVennDiagram} is not installed or N > 7.}
#'     \item{`results`}{List of per-threshold results, each containing
#'       `threshold`, `venn_data` (named list of component IDs per set),
#'       and `n_shared` (named list of shared cluster counts per pair).}
#'     \item{`labels`}{Character vector of phyloseq object names.}
#'     \item{`thresholds`}{The thresholds used.}
#'     \item{`k`}{The k-mer size used.}
#'   }
#'
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency
#'
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' lpq <- list_phyloseq(list(
#'   a = data_fungi_mini,
#'   b = prune_taxa(taxa_names(data_fungi_mini)[1:20], data_fungi_mini)
#' ))
#' res <- refseq_comp_lpq(lpq)
refseq_comp_lpq <- function(
  x,
  thresholds = c(0, 1, 2, 3, 5, 10),
  k = 7,
  verbose = TRUE
) {
  stopifnot(inherits(x, "comparpq::list_phyloseq"))

  pq_list <- x@phyloseq_list
  labels <- names(pq_list)
  n_obj <- length(pq_list)

  if (n_obj < 2) {
    stop("list_phyloseq must contain at least 2 phyloseq objects.")
  }

  missing_refseq <- vapply(
    pq_list,
    \(ps) is.null(phyloseq::refseq(ps, errorIfNULL = FALSE)),
    logical(1)
  )
  if (any(missing_refseq)) {
    stop(
      "The following phyloseq objects have no @refseq slot: ",
      paste(labels[missing_refseq], collapse = ", ")
    )
  }

  seqs_list <- lapply(pq_list, \(ps) as.character(ps@refseq))
  ns <- lengths(seqs_list)
  offsets <- c(0L, cumsum(ns[-n_obj]))

  if (verbose) {
    message("Building k-mer profiles (k = ", k, ") for ", n_obj, " objects...")
  }
  km_list <- lapply(seqs_list, function(seqs) {
    dna <- DNAStringSet(seqs)
    (oligonucleotideFrequency(dna, width = k) > 0) * 1L
  })
  nk_list <- lapply(km_list, rowSums)

  jaccard_to_pid <- function(j) 100 * (2 * j / (1 + j))^(1 / k)

  pairs <- utils::combn(n_obj, 2, simplify = FALSE)
  pair_names <- vapply(
    pairs,
    \(p) paste(labels[p[1]], labels[p[2]], sep = "_vs_"),
    character(1)
  )

  if (verbose) {
    message("Computing ", length(pairs), " pairwise similarity matrices...")
  }
  pid_mats <- lapply(pairs, function(p) {
    i <- p[1]
    j <- p[2]
    if (verbose) {
      message("  ", labels[i], " vs ", labels[j])
    }
    shared <- tcrossprod(km_list[[i]], km_list[[j]])
    union_size <- outer(nk_list[[i]], nk_list[[j]], "+") - shared
    jaccard_to_pid(shared / pmax(union_size, 1))
  })
  names(pid_mats) <- pair_names

  process_threshold <- function(thresh) {
    if (thresh == 0) {
      all_seqs <- unlist(seqs_list, use.names = FALSE)
      seq_groups <- split(seq_along(all_seqs), all_seqs)
      edges <- do.call(
        rbind,
        lapply(seq_groups, function(idx) {
          if (length(idx) < 2) {
            return(NULL)
          }
          t(utils::combn(idx, 2))
        })
      )
    } else {
      min_pid_val <- 100 - thresh
      edges <- do.call(
        rbind,
        lapply(seq_along(pairs), function(p_idx) {
          p <- pairs[[p_idx]]
          e <- which(pid_mats[[p_idx]] >= min_pid_val, arr.ind = TRUE)
          if (nrow(e) == 0) {
            return(NULL)
          }
          cbind(e[, 1] + offsets[p[1]], e[, 2] + offsets[p[2]])
        })
      )
    }

    comp_id <- uf_components(sum(ns), edges)

    node_set <- rep(labels, ns)
    ct <- table(comp_id, node_set)
    in_set <- lapply(labels, \(lbl) ct[, lbl] > 0)
    names(in_set) <- labels
    comp_names <- paste0("c", rownames(ct))

    venn_data <- lapply(labels, \(lbl) comp_names[in_set[[lbl]]])
    names(venn_data) <- labels

    n_shared <- stats::setNames(
      lapply(pairs, \(p) sum(in_set[[labels[p[1]]]] & in_set[[labels[p[2]]]])),
      pair_names
    )

    list(threshold = thresh, venn_data = venn_data, n_shared = n_shared)
  }

  results <- lapply(thresholds, process_threshold)

  venn_plots <- NULL
  if (n_obj <= 7 && requireNamespace("ggVennDiagram", quietly = TRUE)) {
    venn_plots <- lapply(results, function(r) {
      thresh <- r$threshold
      ggVennDiagram::ggVennDiagram(
        r$venn_data,
        label_alpha = 0,
        label = "count"
      ) +
        ggplot2::scale_fill_gradient(low = "white", high = "#4CAF50") +
        ggplot2::labs(
          title = if (thresh == 0) {
            "Exact matches"
          } else {
            sprintf("Within %d%% dissim.", thresh)
          }
        ) +
        ggplot2::theme(legend.position = "none")
    })
    names(venn_plots) <- paste0(thresholds, "_pct")
  } else if (n_obj > 7 && verbose) {
    message(
      "Venn diagrams suppressed: ggVennDiagram supports at most 7 sets (",
      n_obj,
      " provided)."
    )
  }

  structure(
    list(
      venn_plots = venn_plots,
      results = results,
      labels = labels,
      thresholds = thresholds,
      k = k
    ),
    class = "refseq_comp_lpq_result"
  )
}


#' @export
print.refseq_comp_lpq_result <- function(x, ...) {
  cat(
    "== Refseq similarity comparison across",
    length(x$labels),
    "phyloseq objects ==\n"
  )
  cat("Objects    :", paste(x$labels, collapse = ", "), "\n")
  cat("K-mer size :", x$k, "\n\n")
  cat(sprintf(
    "  %-6s  %s\n",
    "Dissim.",
    paste(
      vapply(
        names(x$results[[1]]$n_shared),
        \(pn) {
          parts <- strsplit(pn, "_vs_")[[1]]
          sprintf("%s vs %s", parts[1], parts[2])
        },
        character(1)
      ),
      collapse = "  "
    )
  ))
  for (r in x$results) {
    counts <- vapply(r$n_shared, as.character, character(1))
    cat(sprintf(
      "  %-6s  %s\n",
      paste0(r$threshold, "%"),
      paste(counts, collapse = "        ")
    ))
  }
  invisible(x)
}
