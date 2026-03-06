# to document

div_pq <- function(physeq, modality, indices = "shannon", scales = NULL, hill = FALSE, aggregate = FALSE, funs = list(mean = mean, sd = sd)) {
  mod_values <- unique(as.character(sample_data(physeq)[[modality]]))
  res <- lapply(mod_values, \(val) {
    sub <- prune_samples(sample_data(physeq)[[modality]] == val, physeq)
    comm <- as.data.frame(otu_table(sub))
    if (taxa_are_rows(sub)) {
      comm <- t(comm)
    }
    df <- data.frame(row.names = seq_len(nrow(comm)))
    if (!is.null(indices)) {
      df_div <- data.frame(lapply(indices, \(idx) vegan::diversity(comm, index = idx)))
      names(df_div) <- indices
      df <- cbind(df, df_div)
    }
    if (!is.null(scales)) {
      df_renyi <- as.data.frame(vegan::renyi(comm, scales = scales, hill = hill))
      names(df_renyi) <- paste0(if (hill) { "hill_" } else { "renyi_" }, scales)
      df <- cbind(df, df_renyi)
    }
    df[[modality]] <- val
    df
  }) |>
    do.call(what = rbind)
  if (aggregate) {
    div_cols <- setdiff(names(res), modality)
    res <- res |>
      group_by(.data[[modality]]) |>
      summarise(n_samples = n(),
      across(all_of(div_cols), funs, .names = "{.col}_({.fn})"))
  }
  res
}