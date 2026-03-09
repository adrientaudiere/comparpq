#' Diversity indices per sample, optionally grouped by a modality
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange"
#' alt="lifecycle-experimental"></a>
#'
#' @description
#' Computes alpha-diversity indices (via [vegan::diversity()]) and/or
#' Hill numbers / Rényi entropy (via [vegan::renyi()]) for each sample in
#' a phyloseq object. When `modality` is supplied the computation is done
#' separately for each level of that grouping variable and the level label
#' is appended as an extra column. `NA` values in the modality column are
#' kept as a distinct group so that samples with missing metadata are never
#' silently dropped. When `modality = NULL` all samples are analysed
#' together and no grouping column is added to the result.
#'
#' @param physeq (phyloseq, required) A phyloseq object.
#' @param modality (character or NULL, default NULL) Name of a column in
#'   `sample_data(physeq)` used to split samples into groups. `NA` values
#'   are treated as a separate group. When `NULL`, indices are computed for
#'   the whole dataset without grouping.
#' @param indices (character, default `"shannon"`) One or more index names
#'   accepted by [vegan::diversity()] (e.g. `"shannon"`, `"simpson"`,
#'   `"invsimpson"`). Set to `NULL` to skip classical diversity indices.
#' @param scales (numeric or NULL, default NULL) Scale values passed to
#'   [vegan::renyi()]. When `NULL`, Hill / Rényi computation is skipped.
#' @param hill (logical, default TRUE) If `TRUE`, return Hill numbers;
#'   if `FALSE`, return Rényi entropy. Passed to [vegan::renyi()]. Only
#'   relevant when `scales` is not `NULL`.
#' @param aggregate (logical, default FALSE) If `TRUE` and `modality` is
#'   not `NULL`, summarise per-group results using `funs`.
#' @param funs (named list, default `list(mean = mean, sd = sd)`) Summary
#'   functions applied when `aggregate = TRUE`. Passed to [dplyr::across()].
#'
#' @return A data frame with one row per sample. Columns correspond to the
#'   requested `indices` and/or Hill / Rényi scales. When `modality` is
#'   supplied, an additional column named after `modality` identifies the
#'   group (possibly `NA`). When `aggregate = TRUE` rows are collapsed to
#'   one per group with summary statistics.
#' @author Adrien Taudière
#' 
#' @examples
#' div_pq(data_fungi_mini, indices = c("shannon", "simpson"))
#'
#' div_pq(
#'   data_fungi_mini,
#'   modality = "Height",
#'   indices = c("shannon", "simpson"),
#'   scales = c(0, 1, 2),
#'   hill = TRUE
#' )
#'
#' div_pq(
#'   data_fungi_mini,
#'   modality = "Height",
#'   indices = "shannon",
#'   aggregate = TRUE
#' )
#' @seealso [vegan::diversity()], [vegan::renyi()]
#' @export
div_pq <- function(
  physeq,
  modality = NULL,
  indices = "shannon",
  scales = NULL,
  hill = TRUE,
  aggregate = FALSE,
  funs = list(mean = mean, sd = sd)
) {
  .compute_div <- function(sub) {
    comm <- as.data.frame(otu_table(sub))
    if (taxa_are_rows(sub)) {
      comm <- t(comm)
    }
    df <- data.frame(row.names = seq_len(nrow(comm)))
    if (!is.null(indices)) {
      df_div <- data.frame(lapply(indices, \(idx) {
        vegan::diversity(comm, index = idx)
      }))
      names(df_div) <- indices
      df <- cbind(df, df_div)
    }
    if (!is.null(scales)) {
      df_renyi <- as.data.frame(rbind(vegan::renyi(
        comm,
        scales = scales,
        hill = hill
      )))
      names(df_renyi) <- paste0(if (hill) "hill_" else "renyi_", scales)
      df <- cbind(df, df_renyi)
    }
    df
  }

  if (is.null(modality)) {
    return(.compute_div(physeq))
  }

  mod_col <- sample_data(physeq)[[modality]]
  mod_values <- unique(mod_col)
  res <- lapply(mod_values, \(val) {
    if (is.na(val)) {
      mask <- is.na(mod_col)
    } else {
      mask <- !is.na(mod_col) & mod_col == val
    }
    df <- .compute_div(prune_samples(mask, physeq))
    df[[modality]] <- val
    df
  }) |>
    do.call(what = rbind)
  if (aggregate) {
    div_cols <- setdiff(names(res), modality)
    res <- res |>
      group_by(.data[[modality]]) |>
      summarise(
        n_samples = n(),
        across(all_of(div_cols), funs, .names = "{.col}_({.fn})")
      )
  }
  res
}
