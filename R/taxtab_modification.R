################################################################################
#' Resolve taxonomic conflict in the tax_table of a phyloseq object
#'
#' @param physeq A phyloseq object
#' @param pattern_tax_ranks A vector of pattern to aggregate taxonomic ranks.
#'   For example "^Genus" stand for all taxonomic ranks (columns in tax_table
#'   slot) starting by Genus
#' @param method One of "consensus", "rel_majority", "abs_majority",
#'   "preference" or "unanimity". See details in the documentation of the
#'   function [MiscMetabar::resolve_vector_ranks()].
#' @param keep_tax_ranks (logical, default TRUE) Do we keep the old taxonomic
#'   ranks?
#' @param new_names A vector of new names for the taxonomic columns
#' @param strict (logical, default FALSE). If TRUE, NA are considered as
#'   informative in resolving conflict (i.e. NA are taking into account in vote).
#'   See details for more informations.
#' @param second_method One of "consensus", "rel_majority", "abs_majority",
#'   or "unanimity". Only used if method = "preference". See details.
#' @param nb_agree_threshold (Int, default 1) The minimum number of times a
#'   value must arise to be selected using vote. If 2, we only kept
#'   taxonomic value present at least 2 times in the vector.
#' @param preference_pattern A pattern to match the only column used as prefered
#'   one if method = "preference".
#' @param collapse_string (default '/'). The character to collapse taxonomic names
#'   when multiple assignment is done.
#' @param replace_collapsed_rank_by_NA (logical, default FALSE). If set to TRUE,
#'   all multiple assignments (all taxonomic rank including the 'collapse_string'
#'   parameter) are replaced by NA.
#'
#' @returns A phyloseq object
#' @export
#'
#' @author Adrien Taudière
#'
#' @examples
#'
#' data_fungi_mini_new <- assign_sintax(data_fungi_mini,
#'   ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz",
#'     package = "MiscMetabar"
#'   ),
#'   behavior = "add_to_phyloseq"
#' )
#'
#' resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "consensus", new_names = "Genus_consensus")@tax_table
#' resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "consensus", new_names = "Genus_consensus", replace_collapsed_rank_by_NA = TRUE)@tax_table
#'
#' resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "preference", preference_pattern = ".x$")@tax_table
#' resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "abs_majority")@tax_table
#' resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "rel_majority")@tax_table
#' resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "unanimity")@tax_table
#'
#' resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\.", "^Family\\.", "^Species\\."), method = "consensus")@tax_table
resolve_taxo_conflict <- function(physeq,
                                  pattern_tax_ranks = NULL,
                                  method = c(
                                    "consensus",
                                    "rel_majority",
                                    "abs_majority",
                                    "preference",
                                    "unanimity"
                                  ),
                                  strict = FALSE,
                                  second_method = c(
                                    "consensus",
                                    "rel_majority",
                                    "abs_majority",
                                    "preference",
                                    "unanimity"
                                  ),
                                  nb_agree_threshold = 1,
                                  keep_tax_ranks = TRUE,
                                  new_names = NULL,
                                  preference_pattern = NULL,
                                  collapse_string = "/",
                                  replace_collapsed_rank_by_NA = FALSE) {
  method <- match.arg(method)

  verify_pq(physeq)
  new_physeq <- physeq
  taxtab <- as_tibble(data.frame(physeq@tax_table))

  if (is.null(new_names)) {
    new_names <- paste0(pattern_tax_ranks, "_", method[1])
  }


  for (i in seq_along(pattern_tax_ranks)) {
    new_tax_ranks <- taxtab |>
      dplyr::select(matches(pattern_tax_ranks[[i]]))

    if (method == "preference") {
      preference_index <- grep(preference_pattern, colnames(new_tax_ranks))
      if (length(preference_index) > 1) {
        stop(
          "Your preference_pattern match multiple rank name
           (column in the tax_table slot) and must only match one."
        )
      } else if (length(preference_index) == 0) {
        stop(
          "Your preference_pattern do not match any rank name
           (column in the tax_table slot)"
        )
      }
    } else {
      preference_index <- NULL
    }

    new_tax_r <- apply(new_tax_ranks, 1, function(x) {
      resolve_vector_ranks(
        x,
        method = method,
        preference_index = preference_index,
        collapse_string = collapse_string,
        replace_collapsed_rank_by_NA = replace_collapsed_rank_by_NA
      )
    })

    taxtab <- tibble::add_column(taxtab, !!new_names[[i]] := new_tax_r)
  }

  if (!keep_tax_ranks) {
    taxtab <- taxtab |>
      dplyr::select(!matches(pattern_tax_ranks))
  }

  new_physeq@tax_table <- tax_table(as.matrix(taxtab))
  taxa_names(new_physeq@tax_table) <- taxa_names(physeq)

  return(new_physeq)
}
################################################################################





################################################################################
#' Select taxonomic ranks in a phyloseq object
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful in tidyverse-like pipeline.
#'
#' @param physeq A phyloseq object
#' @param ... <tidy-select> One or more unquoted expressions separated by commas. Variable names can be used as if they were positions in the data frame, so expressions like x:y can be used to select a range of variables. See ?dplyr::select
#' @returns A phyloseq object
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' select_ranks_pq(data_fungi, Order, Family)@tax_table |>
#'   dim()
#' select_ranks_pq(data_fungi, !Order)@tax_table |>
#'   colnames()
#' #
select_ranks_pq <- function(physeq, ...) {
  verify_pq(physeq)
  new_physeq <- physeq

  taxtab <- as_tibble(data.frame(physeq@tax_table)) |>
    select(...)

  new_physeq@tax_table <- tax_table(as.matrix(taxtab))
  taxa_names(new_physeq@tax_table) <- taxa_names(physeq)

  verify_pq(new_physeq)
  return(new_physeq)
}

################################################################################

################################################################################
#' Rename names of ranks in the tax_table slot of a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' The users can use the function using a couple of vector `old_names`/`new_names`
#' or using a `pattern` to replace by `replacement`
#'
#' @param physeq A phyloseq object
#' @param old_names Names of the names to replace
#' @param new_names Names of the new names
#' @param pattern Pattern to replace by the args replacement
#' @param replacement Pattern of replacement
#' @param fixed see ?grep
#' @param perl see ?grep
#' @param useBytes see ?grep
#'
#' @returns An object of class phyloseq
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' rename_ranks_pq(data_fungi, c("Confidence.Ranking", "Phylum"), c("Conf.Rank.Guild", "Phyla"))@tax_table |>
#'   colnames()
#' rename_ranks_pq(data_fungi, pattern = ".", replacement = "_", fixed = TRUE)@tax_table |>
#'   colnames()
rename_ranks_pq <- function(physeq,
                            old_names = NULL,
                            new_names = NULL,
                            pattern = NULL,
                            replacement = NULL,
                            fixed = FALSE,
                            perl = FALSE,
                            useBytes = FALSE) {
  verify_pq(physeq)
  new_physeq <- physeq

  if (!is.null(pattern)) {
    if (!is.null(old_names) || !is.null(new_names)) {
      stop(
        "You must specify either (old_names and new_names) OR (pattern and replacement) parameters."
      )
    }
    old_names <- grep(
      pattern,
      colnames(new_physeq@tax_table),
      value = TRUE,
      fixed = fixed,
      perl = perl,
      useBytes = useBytes
    )
    new_names <- gsub(
      pattern,
      replacement,
      old_names,
      fixed = fixed,
      perl = perl,
      useBytes = useBytes
    )
  }

  rename_rank_pq <- function(physeq,
                             old_name = NULL,
                             new_name = NULL) {
    verify_pq(physeq)
    new_physeq <- physeq

    taxtab <- as_tibble(data.frame(physeq@tax_table))
    colnames(taxtab) <- gsub(old_name, new_name, colnames(taxtab), fixed = TRUE)

    new_physeq@tax_table <- tax_table(as.matrix(taxtab))
    taxa_names(new_physeq@tax_table) <- taxa_names(physeq)

    verify_pq(new_physeq)
    return(new_physeq)
  }

  for (i in seq_along(old_names)) {
    new_physeq <- rename_rank_pq(new_physeq, old_names[i], new_names[i])
  }
  return(new_physeq)
}
################################################################################




################################################################################
#' Replace taxonomic value with a given pattern by NA
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Use [base::gsub()] for substitution.
#'
#'
#' @inheritParams tc_points_matrix
#' @param patterns A vector of patterns to select taxonomic value to convert to NA
#' @param taxonomic_ranks A list of taxonomic ranks where we want the
#'   substitution occur. If left to NULL, all taxonomic ranks are modified.
#' @param progress_bar (logical, default FALSE) Do we print progress during
#'   the calculation?
#' @param ... Other arguments to be passed on to [base::gsub()]
#' @return A phyloseq object
#' @export
#'
#' @author Adrien Taudière
#' @examples
#' data_fungi@tax_table["ASV85", "Family"]
#' data_fungi2 <- taxtab_replace_pattern_by_NA(data_fungi, "fam_Incertae_sedis", taxonomic_ranks = "Family")
#' data_fungi2@tax_table["ASV85", "Family"]
#'
#' # By default patterns ".*_incertae_sedis" and "unclassified.*" are replaced by NA
#' data_fungi3 <-
#'   taxtab_replace_pattern_by_NA(data_fungi, ignore.case = TRUE, progress_bar = TRUE)
#' data_fungi3@tax_table["ASV85", "Family"]
taxtab_replace_pattern_by_NA <- function(physeq,
                                         patterns = c(".*_incertae_sedis", "unclassified.*"),
                                         taxonomic_ranks = NULL,
                                         progress_bar = FALSE,
                                         ...) {
  if (is.null(taxonomic_ranks)) {
    taxonomic_ranks <- colnames(physeq@tax_table)
  }
  taxtab <- as.matrix(unclass(physeq@tax_table))

  if (progress_bar) {
    pb <- txtProgressBar(
      min = 0,
      max = length(taxonomic_ranks) * length(patterns),
      style = 3,
      width = 50,
      char = "="
    )
    pb_val <- 0
  }

  for (pat in patterns) {
    for (taxrank in taxonomic_ranks) {
      if (progress_bar) {
        pb_val <- pb_val + 1
        setTxtProgressBar(pb, pb_val)
      }
      taxtab[, taxonomic_ranks] <-
        gsub(pat, NA, taxtab[, taxonomic_ranks], ...)
    }
  }

  physeq@tax_table <- tax_table(taxtab)
  return(physeq)
}
################################################################################
