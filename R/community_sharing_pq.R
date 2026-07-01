################################################################################
# ── Community sharing plot ────────────────────────────────────────────────────

# Internal helper: aggregate an OTU table by modality (sum of reads per group).
# Returns a matrix with taxa as rows and modalities as columns.
.agg_by_mod <- function(ps, fact, mods) {
  otu <- as.data.frame(as.matrix(
    if (phyloseq::taxa_are_rows(ps)) {
      phyloseq::otu_table(ps)
    } else {
      t(phyloseq::otu_table(ps))
    }
  ))
  sd <- as.data.frame(phyloseq::sample_data(ps))
  sapply(mods, \(mod) {
    samps <- intersect(rownames(sd)[sd[[fact]] == mod], colnames(otu))
    rowSums(otu[, samps, drop = FALSE])
  })
}

# Default node angles (degrees) for 2-4 modalities.
.sharing_angles <- function(n) {
  switch(
    as.character(n),
    "2" = c(180, 0),
    "3" = c(90, 210, 330),
    "4" = c(135, 225, 315, 45)
  )
}

#' Build a single metric definition for `community_sharing_pq()`
#'
#' Build a single metric definition for `community_sharing_pq()`.
#'
#' @description
#' Build a single metric definition for `community_sharing_pq()`.
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @param label Character. Human-readable label used in the legend.
#' @param color Character. Color of the curve connecting modalities.
#' @param fn Function with signature `function(a, b, otu_sp, cache)` returning
#'   a single numeric value. `a`, `b` are modality names; `otu_sp` is the
#'   species-level OTU table aggregated by modality; `cache` is the object
#'   returned by `prep` (or `NULL`).
#' @param fmt `sprintf` format string for legend stats (min / mean / max).
#'   Default `"%.2f"`; use `"%.0f"` for integer counts.
#' @param prep Optional. Function `function(physeq, fact, modalities)` called
#'   once to precompute helper data shared across all pairs (e.g. genus-level
#'   aggregation). Returns an object passed as `cache` to `fn`.
#' @param bounds Either `NULL` (default) or a numeric vector `c(lower, upper)`
#'   giving the theoretical range of the metric. When non-`NULL`, linewidth is
#'   scaled to `linewidth_range` using these fixed bounds, so the same metric
#'   value always maps to the same linewidth regardless of the observed data.
#'   When `NULL`, linewidth is rescaled within the observed range of values.
#'   Use `bounds = c(0, 1)` for metrics naturally expressed as proportions or
#'   similarities (e.g. Jaccard, Bray-Curtis). All metrics should be expressed
#'   as similarities: higher value = more similar communities = thicker link.
#'
#' @return A list with the metric definition.
#' @examples
#' \dontrun{
#' # Custom metric: robust Aitchison similarity (vegan >= 2.6).
#' # The method uses matrix completion and needs the full matrix, so we
#' # compute the full distance once in `prep` and look up each pair in `fn`.
#' rob_aitch <- make_sharing_metric(
#'   label = "Robust Aitchison similarity",
#'   color = "#FF7F00",
#'   fmt   = "%.2f",
#'   prep  = function(physeq, fact, modalities) {
#'     otu_sp <- .agg_by_mod(physeq, fact, modalities)
#'     as.matrix(vegan::vegdist(t(otu_sp), method = "robust.aitchison"))
#'   },
#'   fn    = \(a, b, otu_sp, cache) 1 - cache[a, b]
#' )
#'
#' # Custom metric: number of unique (non-shared) species
#' unique_sp <- make_sharing_metric(
#'   label = "Unique species count",
#'   color = "#A65628",
#'   fmt   = "%.0f",
#'   fn    = \(a, b, otu_sp, cache) {
#'     sum(xor(otu_sp[, a] > 0, otu_sp[, b] > 0))
#'   }
#' )
#'
#' # Custom metric with prep step: Sorensen similarity at Family rank
#' sor_family <- make_sharing_metric(
#'   label = "Sorensen at Family rank",
#'   color = "#F781BF",
#'   fmt   = "%.2f",
#'   prep  = function(physeq, fact, modalities) {
#'     d_fam <- phyloseq::tax_glom(physeq, "Family", NArm = FALSE)
#'     .agg_by_mod(d_fam, fact, modalities)
#'   },
#'   fn    = \(a, b, otu_sp, cache) {
#'     1 - as.numeric(vegan::vegdist(
#'       t(cache[, c(a, b)]),
#'       method = "bray", binary = TRUE
#'     ))
#'   }
#' )
#' }
#' @export
#' @author Adrien Taudière
make_sharing_metric <- function(
  label,
  color,
  fn,
  fmt = "%.2f",
  prep = NULL,
  bounds = NULL
) {
  list(
    label = label,
    color = color,
    fn = fn,
    fmt = fmt,
    prep = prep,
    bounds = bounds
  )
}

#' Default metrics for `community_sharing_pq()`
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Returns the 4 default metrics: shared species count, Bray-Curtis similarity,
#' Jaccard binary similarity, and proportion of shared genera.
#'
#' `bray_sim`, `jac_sim`, and `genus_prop` use `bounds = c(0, 1)` for globally
#' consistent linewidth scaling. `shared_sp` uses `bounds = NULL` (count metric
#' without a fixed upper bound; rescaled within the observed range).
#'
#' @return A named list of metric definitions.
#' @examples
#' \dontrun{
#' # Default set of 4 metrics
#' mets <- default_sharing_metrics()
#' community_sharing_pq(data_fungi, fact = "Height", metrics = mets)
#'
#' # Extend defaults with a robust Aitchison similarity metric
#' mets_plus <- c(
#'   default_sharing_metrics(),
#'   list(
#'     rob_aitch = make_sharing_metric(
#'       label = "Robust Aitchison similarity",
#'       color = "#FF7F00",
#'       prep  = function(physeq, fact, modalities) {
#'         otu_sp <- .agg_by_mod(physeq, fact, modalities)
#'         as.matrix(vegan::vegdist(t(otu_sp), method = "robust.aitchison"))
#'       },
#'       fn    = \(a, b, otu_sp, cache) 1 - cache[a, b]
#'     )
#'   )
#' )
#' community_sharing_pq(data_fungi, fact = "Height", metrics = mets_plus)
#'
#' # Keep only a subset of the defaults
#' community_sharing_pq(
#'   data_fungi,
#'   fact    = "Height",
#'   metrics = default_sharing_metrics()[c("shared_sp", "bray_sim")]
#' )
#' }
#' @export
#' @author Adrien Taudière
default_sharing_metrics <- function() {
  list(
    shared_sp = make_sharing_metric(
      label = "Shared species",
      color = "#E41A1C",
      fmt = "%.0f",
      bounds = NULL,
      fn = \(a, b, otu_sp, cache) sum(otu_sp[, a] > 0 & otu_sp[, b] > 0)
    ),
    bray_sim = make_sharing_metric(
      label = "Bray-Curtis similarity",
      color = "#377EB8",
      fmt = "%.2f",
      bounds = c(0, 1),
      fn = \(a, b, otu_sp, cache) {
        1 - as.numeric(vegan::vegdist(t(otu_sp[, c(a, b)]), method = "bray"))
      }
    ),
    jac_sim = make_sharing_metric(
      label = "Jaccard similarity",
      color = "#4DAF4A",
      fmt = "%.2f",
      bounds = c(0, 1),
      fn = \(a, b, otu_sp, cache) {
        1 -
          as.numeric(vegan::vegdist(
            t(otu_sp[, c(a, b)]),
            method = "jaccard",
            binary = TRUE
          ))
      }
    ),
    genus_prop = make_sharing_metric(
      label = "Shared genera (prop.)",
      color = "#984EA3",
      fmt = "%.2f",
      bounds = c(0, 1),
      prep = function(physeq, fact, modalities) {
        d_gen <- phyloseq::tax_glom(physeq, "Genus", NArm = FALSE)
        .agg_by_mod(d_gen, fact, modalities)
      },
      fn = \(a, b, otu_sp, cache) {
        ga <- rownames(cache)[cache[, a] > 0]
        gb <- rownames(cache)[cache[, b] > 0]
        length(intersect(ga, gb)) / length(union(ga, gb))
      }
    )
  )
}

# Internal helper: compute pairwise metrics for a given physeq + modalities.
# Returns a list with $pairs_df, $otu_sp, $caches.
.compute_pairs_df <- function(physeq, fact, modalities, metrics) {
  otu_sp <- .agg_by_mod(physeq, fact, modalities)
  caches <- lapply(metrics, \(m) {
    if (is.null(m$prep)) NULL else m$prep(physeq, fact, modalities)
  })
  pairs_df <- purrr::map_dfr(
    combn(modalities, 2, simplify = FALSE),
    function(pair) {
      row <- dplyr::tibble(from = pair[1], to = pair[2])
      for (nm in names(metrics)) {
        row[[nm]] <- metrics[[nm]]$fn(pair[1], pair[2], otu_sp, caches[[nm]])
      }
      row
    }
  )
  list(pairs_df = pairs_df, otu_sp = otu_sp, caches = caches)
}

#' Community sharing plot: modalities as pie nodes with multi-metric links
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Draws a figure with one node per modality of `fact` (2 to 4 supported),
#' positioned on a regular polygon. Each node is a pie chart showing the
#' taxonomic composition at rank `pie_taxrank`. Between each pair of nodes,
#' one curved link per metric is drawn; linewidth is rescaled within each
#' metric to `linewidth_range`. A legend below the figure gives each metric's
#' min / mean / max over all pairs.
#'
#' @param physeq A `phyloseq` object.
#' @param fact Character. Name of a `sample_data(physeq)` column used to group
#'   samples into modalities. Must have 2 to 4 unique values.
#' @param metrics Named list of metric definitions. See
#'   [default_sharing_metrics()] and [make_sharing_metric()] for the format.
#'   Default uses 4 metrics: shared species, Bray-Curtis similarity, Jaccard
#'   similarity, and proportion of common genera.
#' @param pie_taxrank Character. Taxonomic rank for the pie charts. Default
#'   `"Class"`.
#' @param pie_r Numeric. Pie radius in data units. Default `0.28`.
#' @param label_offset Numeric. Distance between pie edge and node label.
#'   Default `0.18`.
#' @param curvature_range Numeric of length 2. Range of `geom_curve` curvatures
#'   used to fan out metrics between a pair of nodes. Default `c(-0.35, 0.35)`.
#'   The sign is flipped per pair so that the first metric always fans toward
#'   the plot centre (and the last toward the border), regardless of which
#'   side of the plot the pair is on.
#' @param linewidth_range Numeric of length 2. `linewidth_range[1]` (minimum
#'   linewidth) corresponds to the weakest similarity or sharing; `linewidth_range[2]`
#'   (maximum linewidth) corresponds to the strongest. All built-in metrics are
#'   expressed as similarities (higher value = more similar communities = thicker
#'   link). For bounded metrics (`bounds = c(0, 1)`), the mapping is global
#'   across runs; for unbounded metrics it is scaled within the observed range.
#'   Default `c(0.8, 4.5)`.
#' @param n_perm Integer. Number of label-permutation iterations for significance
#'   testing. Default `0` (no test). See Details.
#' @param sig_threshold Numeric in `(0, 1)`. p-value threshold below which a
#'   link is considered significant. Default `0.05`. Only used when `n_perm > 0`.
#' @param nonsig_alpha Numeric in `[0, 1]`. Alpha (transparency) applied to
#'   non-significant links. Default `0.15`. Significant links keep alpha `0.75`.
#'   Only used when `n_perm > 0`.
#' @param seed Integer or `NULL`. Passed to [set.seed()] before the permutation
#'   loop for reproducibility. Default `NULL` (no seed).
#' @param palette Either the name of an RColorBrewer qualitative palette (e.g.
#'   `"Set3"`) or a character vector of fill colors used for the top taxa.
#' @param max_taxa Integer. Keep colors for the `max_taxa` most abundant taxa
#'   (summed across all modalities) and collapse the rest into an `"Other"`
#'   category filled with `other_color`. Default `12`.
#' @param other_color Fill color for the `"Other"` category. Default
#'   `"grey70"`.
#' @param show_na_modality Logical. If `TRUE`, samples whose `fact` value is
#'   `NA` are grouped into an extra `"NA"` modality and displayed as an
#'   additional node. The total number of modalities (including `"NA"`) must
#'   still be between 2 and 4. Default `FALSE`.
#' @param show_na Logical. If `TRUE` (default), taxa with an `NA`/empty value
#'   at `pie_taxrank` are shown as a dedicated `"NA"` category filled with
#'   `na_color`. If `FALSE`, they are dropped.
#' @param na_color Fill color for the `"NA"` category. Default `"grey40"`.
#' @param pie_border_color Color of the circle drawn around each pie. Default
#'   `"black"`.
#' @param pie_border_width Linewidth of the pie border. Default `0.6`.
#' @param title Plot title (optional).
#' @param base_size Base font size. Default `12`.
#'
#' @details
#' **Permutation null model.** When `n_perm > 0`, significance is assessed by
#' label permutation: the `fact` column in `sample_data(physeq)` is shuffled
#' uniformly at random among samples (preserving group sizes), then all `prep()`
#' functions and metric computations are re-run on the permuted data. This is
#' repeated `n_perm` times. The empirical p-value for each (pair, metric)
#' combination is the proportion of permuted values greater than or equal to
#' the observed value (one-sided upper-tail test). This null model is
#' appropriate for all built-in metrics because they compare aggregated
#' community profiles between groups, and the aggregation depends on which
#' samples are assigned to each group. Non-significant links are drawn with
#' alpha `nonsig_alpha` (faded); significant links use alpha `0.75`.
#'
#' **Performance.** Each permutation re-runs all `prep()` functions. For
#' metrics that call `phyloseq::tax_glom()` (e.g. `genus_prop`), this can be
#' slow on large datasets. Recommended range: `n_perm = 99` to `n_perm = 199`.
#'
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' # Default: 4 metrics, pie charts at Class rank
#' community_sharing_pq(data_fungi, fact = "Height")
#'
#' # Custom metrics: add robust Aitchison similarity (vegan >= 2.6)
#' mets <- c(
#'   default_sharing_metrics()[1],
#'   list(
#'     rob_aitch = make_sharing_metric(
#'       label = "Robust Aitchison similarity",
#'       color = "#FF7F00",
#'       prep  = function(physeq, fact, modalities) {
#'         otu_sp <- .agg_by_mod(physeq, fact, modalities)
#'         as.matrix(1 - vegan::vegdist(t(otu_sp), method = "robust.aitchison"))
#'       },
#'       fn    = \(a, b, otu_sp, cache) 1 - cache[a, b]
#'     )
#'   )
#' )
#' community_sharing_pq(data_fungi, fact = "Height", metrics = mets)
#'
#' # Only show a single metric (Jaccard)
#' community_sharing_pq(
#'   data_fungi,
#'   fact    = "Height",
#'   metrics = default_sharing_metrics()["jac_sim"]
#' )
#'
#' # Use Phylum rank for the pie charts and a different palette
#' community_sharing_pq(
#'   data_fungi,
#'   fact        = "Height",
#'   pie_taxrank = "Phylum",
#'   palette     = "Set2",
#'   max_taxa    = 8
#' )
#'
#' # Hide NA taxa in pies and add a title
#' community_sharing_pq(
#'   data_fungi,
#'   fact    = "Height",
#'   show_na = FALSE,
#'   title   = "Community sharing across heights"
#' )
#'
#' # Include samples with NA height as a fourth modality
#' community_sharing_pq(
#'   data_fungi,
#'   fact             = "Height",
#'   show_na_modality = TRUE
#' )
#'
#' # Permutation significance test (99 permutations, faded non-sig links)
#' community_sharing_pq(
#'   data_fungi,
#'   fact   = "Height",
#'   n_perm = 99,
#'   seed   = 42
#' )
#' }
#' @importFrom vegan vegdist
#' @importFrom scales rescale
#' @importFrom tidyr pivot_longer
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom purrr map_dfr pmap_chr map_lgl
#' @export
#' @author Adrien Taudière
community_sharing_pq <- function(
  physeq,
  fact,
  metrics = default_sharing_metrics(),
  pie_taxrank = "Class",
  pie_r = 0.28,
  label_offset = 0.18,
  curvature_range = c(-0.35, 0.35),
  linewidth_range = c(0.8, 4.5),
  n_perm = 0,
  sig_threshold = 0.05,
  nonsig_alpha = 0.15,
  seed = NULL,
  palette = "Set3",
  max_taxa = 12,
  other_color = "grey70",
  show_na_modality = FALSE,
  show_na = TRUE,
  na_color = "grey40",
  pie_border_color = "black",
  pie_border_width = 0.6,
  title = NULL,
  base_size = 12
) {
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop(
      "Package 'ggforce' is required for community_sharing_pq(). ",
      "Install it with: install.packages('ggforce')"
    )
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop(
      "Package 'RColorBrewer' is required for community_sharing_pq(). ",
      "Install it with: install.packages('RColorBrewer')"
    )
  }

  sd <- as.data.frame(phyloseq::sample_data(physeq))
  if (!fact %in% colnames(sd)) {
    stop("Column '", fact, "' not found in sample_data.")
  }

  fact_values <- sd[[fact]]

  # Promote NA samples to an explicit "NA" modality when requested
  if (show_na_modality && anyNA(fact_values)) {
    fact_values <- ifelse(is.na(fact_values), "NA", as.character(fact_values))
    sd[[fact]] <- fact_values
    phyloseq::sample_data(physeq) <- phyloseq::sample_data(sd)
  }

  modalities <- if (is.factor(fact_values)) {
    levels(droplevels(fact_values))
  } else {
    sort(unique(stats::na.omit(fact_values)))
  }
  n_mod <- length(modalities)
  if (n_mod < 2 || n_mod > 4) {
    stop("Number of modalities must be between 2 and 4, got ", n_mod, ".")
  }

  # ── Node positions ─────────────────────────────────────────────────────────
  angles <- .sharing_angles(n_mod) * pi / 180
  n_per_mod <- as.data.frame(table(fact_values))
  colnames(n_per_mod) <- c("name", "n")
  n_per_mod$name <- as.character(n_per_mod$name)

  nodes <- dplyr::tibble(
    name = modalities,
    x = cos(angles),
    y = sin(angles)
  ) |>
    dplyr::left_join(n_per_mod, by = "name") |>
    dplyr::mutate(
      label_r = 1 + pie_r + label_offset,
      x_label = x * label_r,
      y_label = y * label_r,
      hjust = dplyr::case_when(abs(x) < 0.1 ~ 0.5, x > 0 ~ 0, TRUE ~ 1),
      vjust = dplyr::case_when(y > 0.3 ~ 0, y < -0.3 ~ 1, TRUE ~ 0.5),
      node_label = paste0(name, "\n(n=", n, ")")
    )

  # ── Pairwise metric computation ────────────────────────────────────────────
  res <- .compute_pairs_df(physeq, fact, modalities, metrics)
  pairs_df <- res$pairs_df

  # ── Permutation significance test ──────────────────────────────────────────
  if (n_perm > 0) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    obs_mat <- as.matrix(pairs_df[, names(metrics)])
    perm_vals <- array(NA_real_, c(n_perm, nrow(pairs_df), length(metrics)))

    cli::cli_progress_bar("Permutations", total = n_perm, clear = FALSE)
    for (i in seq_len(n_perm)) {
      ps_perm <- physeq
      sd_perm <- as.data.frame(phyloseq::sample_data(physeq))
      sd_perm[[fact]] <- sample(sd_perm[[fact]])
      phyloseq::sample_data(ps_perm) <- phyloseq::sample_data(sd_perm)
      res_perm <- .compute_pairs_df(ps_perm, fact, modalities, metrics)
      perm_vals[i, , ] <- as.matrix(res_perm$pairs_df[, names(metrics)])
      cli::cli_progress_update()
    }
    cli::cli_progress_done()

    pval_mat <- matrix(
      NA_real_,
      nrow = nrow(pairs_df),
      ncol = length(metrics),
      dimnames = list(NULL, paste0("pval_", names(metrics)))
    )
    for (j in seq_along(metrics)) {
      for (k in seq_len(nrow(pairs_df))) {
        pval_mat[k, j] <- mean(perm_vals[, k, j] >= obs_mat[k, j])
      }
    }
    pairs_df <- dplyr::bind_cols(pairs_df, as.data.frame(pval_mat))
  }

  # ── Metric metadata (color, curvature, format, bounds) ─────────────────────
  n_met <- length(metrics)
  curvatures <- if (n_met == 1) {
    0
  } else {
    seq(curvature_range[1], curvature_range[2], length.out = n_met)
  }
  metric_meta <- dplyr::tibble(
    metric = names(metrics),
    label = vapply(metrics, \(m) m$label, character(1)),
    color = vapply(metrics, \(m) m$color, character(1)),
    fmt = vapply(metrics, \(m) m$fmt, character(1)),
    curvature = curvatures,
    bounds = lapply(metrics, \(m) m$bounds)
  )

  # ── Link data ──────────────────────────────────────────────────────────────
  # bounds-aware rescale: when from is supplied use global range, otherwise
  # fall back to the observed range (or midpoint for constant values).
  safe_rescale <- function(x, to, from = NULL) {
    if (!is.null(from)) {
      scales::rescale(x, to = to, from = from)
    } else if (length(unique(x)) < 2) {
      rep(mean(to), length(x))
    } else {
      scales::rescale(x, to = to)
    }
  }
  link_df <- pairs_df |>
    tidyr::pivot_longer(
      dplyr::all_of(names(metrics)),
      names_to = "metric",
      values_to = "raw_value"
    ) |>
    dplyr::left_join(metric_meta, by = "metric") |>
    dplyr::left_join(
      dplyr::select(nodes, name, x, y),
      by = c("from" = "name")
    ) |>
    dplyr::rename(x_from = x, y_from = y) |>
    dplyr::left_join(dplyr::select(nodes, name, x, y), by = c("to" = "name")) |>
    dplyr::rename(x_to = x, y_to = y) |>
    dplyr::mutate(
      met_idx = match(metric, names(metrics)),
      # Flip curvature per pair so that metric 1 always sits closest to the
      # plot centre (and metric n closest to the border), regardless of which
      # side of the plot the pair is on.
      pair_flip = {
        dx <- x_to - x_from
        dy <- y_to - y_from
        mx <- (x_from + x_to) / 2
        my <- (y_from + y_to) / 2
        ifelse(dx * my - dy * mx > 0, -1, 1)
      },
      curvature_pair = curvature * pair_flip
    ) |>
    dplyr::group_by(metric) |>
    dplyr::mutate(
      lw = safe_rescale(raw_value, linewidth_range, from = bounds[[1]])
    ) |>
    dplyr::ungroup()

  # Add significance flag
  if (n_perm > 0) {
    pval_long <- pairs_df |>
      dplyr::select(from, to, dplyr::starts_with("pval_")) |>
      tidyr::pivot_longer(
        dplyr::starts_with("pval_"),
        names_to = "metric",
        values_to = "pval",
        names_prefix = "pval_"
      )
    link_df <- dplyr::left_join(
      link_df,
      pval_long,
      by = c("from", "to", "metric")
    ) |>
      dplyr::mutate(significant = pval < sig_threshold)
  } else {
    link_df <- dplyr::mutate(link_df, significant = TRUE)
  }

  # ── Legend (min / mean / max) ──────────────────────────────────────────────
  legend_y0 <- -1.4
  legend_step <- -0.18
  # Row 0 is reserved for the "Metrics" header; metric rows start at row 1.
  legend_df <- pairs_df |>
    tidyr::pivot_longer(
      dplyr::all_of(names(metrics)),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::group_by(metric) |>
    dplyr::summarise(
      vmin = min(value),
      vmean = mean(value),
      vmax = max(value),
      .groups = "drop"
    ) |>
    dplyr::left_join(metric_meta, by = "metric") |>
    dplyr::arrange(match(metric, metric_meta$metric)) |>
    dplyr::mutate(
      stats_label = purrr::pmap_chr(
        list(label, fmt, vmin, vmean, vmax),
        function(lbl, f, mn, me, mx) {
          sprintf(
            paste0("%s  [min: ", f, "  mean: ", f, "  max: ", f, "]"),
            lbl,
            mn,
            me,
            mx
          )
        }
      ),
      x = -1.65,
      xend = -1.35,
      y = legend_y0 + legend_step * dplyr::row_number(),
      yend = y
    )

  # ── Pie charts (composition at rank `pie_taxrank`) ─────────────────────────
  d_rank <- phyloseq::tax_glom(physeq, pie_taxrank, NArm = FALSE)
  otu_rank <- .agg_by_mod(d_rank, fact, modalities)
  taxa_lab <- as.character(phyloseq::tax_table(d_rank)[
    rownames(otu_rank),
    pie_taxrank
  ])

  na_mask <- is.na(taxa_lab) | taxa_lab == ""
  if (show_na) {
    taxa_lab[na_mask] <- "NA"
  } else {
    otu_rank <- otu_rank[!na_mask, , drop = FALSE]
    taxa_lab <- taxa_lab[!na_mask]
  }
  rownames(otu_rank) <- taxa_lab

  # Collapse rare taxa into "Other" (NA kept aside as its own category)
  real_taxa <- setdiff(rownames(otu_rank), "NA")
  total_abund <- rowSums(otu_rank[real_taxa, , drop = FALSE])
  n_keep <- min(max_taxa, length(real_taxa))
  top_n <- names(sort(total_abund, decreasing = TRUE))[seq_len(n_keep)]

  is_other <- !(rownames(otu_rank) %in% c(top_n, "NA"))
  if (any(is_other)) {
    other_row <- colSums(otu_rank[is_other, , drop = FALSE])
    otu_rank <- rbind(otu_rank[!is_other, , drop = FALSE], Other = other_row)
  }

  tax_order <- c(
    top_n,
    if ("Other" %in% rownames(otu_rank)) "Other",
    if ("NA" %in% rownames(otu_rank)) "NA"
  )
  otu_rank <- otu_rank[tax_order, , drop = FALSE]

  pie_df <- as.data.frame(otu_rank) |>
    tibble::rownames_to_column(pie_taxrank) |>
    tidyr::pivot_longer(
      -dplyr::all_of(pie_taxrank),
      names_to = "modality",
      values_to = "count"
    ) |>
    dplyr::filter(count > 0) |>
    dplyr::group_by(modality) |>
    dplyr::mutate(prop = count / sum(count)) |>
    dplyr::ungroup() |>
    dplyr::left_join(
      dplyr::select(nodes, name, x, y),
      by = c("modality" = "name")
    ) |>
    dplyr::mutate(
      !!pie_taxrank := factor(.data[[pie_taxrank]], levels = tax_order)
    )

  # geom_curve layers (one per metric x pair_flip x significant)
  make_curve_layers <- function(m) {
    sub <- dplyr::filter(link_df, metric == m)
    meta <- dplyr::filter(metric_meta, metric == m)
    combos <- unique(sub[, c("pair_flip", "significant")])
    lapply(seq_len(nrow(combos)), function(i) {
      flip <- combos$pair_flip[i]
      sig <- combos$significant[i]
      d <- dplyr::filter(sub, pair_flip == flip, significant == sig)
      ggplot2::geom_curve(
        data = d,
        ggplot2::aes(
          x = x_from,
          y = y_from,
          xend = x_to,
          yend = y_to,
          linewidth = lw
        ),
        color = meta$color,
        curvature = meta$curvature * flip,
        alpha = if (sig) 0.75 else nonsig_alpha,
        lineend = "round"
      )
    })
  }

  top_colors <- if (
    length(palette) == 1 &&
      palette %in% rownames(RColorBrewer::brewer.pal.info)
  ) {
    max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    n_use <- max(3, min(length(top_n), max_n))
    cols <- RColorBrewer::brewer.pal(n_use, palette)[seq_len(min(
      length(top_n),
      max_n
    ))]
    if (length(top_n) > max_n) {
      cols <- grDevices::colorRampPalette(cols)(length(top_n))
    }
    cols
  } else {
    rep(palette, length.out = length(top_n))
  }
  names(top_colors) <- top_n
  fill_values <- top_colors
  if ("Other" %in% tax_order) {
    fill_values <- c(fill_values, Other = other_color)
  }
  if ("NA" %in% tax_order) {
    fill_values <- c(fill_values, "NA" = na_color)
  }
  fill_scale <- ggplot2::scale_fill_manual(
    values = fill_values,
    breaks = tax_order,
    name = pie_taxrank
  )

  # ── Significance legend (only when permutations were run) ──────────────────
  if (n_perm > 0) {
    perm_header_y <- legend_y0 + legend_step * (n_met + 1.5)
    perm_desc_y <- legend_y0 + legend_step * (n_met + 2.5)
    sig_rows_y <- legend_y0 + legend_step * (n_met + c(3.5, 4.5))

    perm_desc <- sprintf(
      "%d label permutations (group sizes preserved); p = prop(perm \u2265 obs)",
      n_perm
    )

    sig_legend_df <- dplyr::tibble(
      label = c(
        sprintf("p < %g \u2014 significant", sig_threshold),
        sprintf("p \u2265 %g \u2014 not significant", sig_threshold)
      ),
      y = sig_rows_y,
      sig = c(TRUE, FALSE)
    )
    y_bottom <- sig_rows_y[2] - 0.15
  } else {
    perm_header_y <- NULL
    perm_desc_y <- NULL
    perm_desc <- NULL
    sig_legend_df <- NULL
    y_bottom <- legend_y0 + legend_step * n_met - 0.15
  }
  y_top <- 1 + pie_r + label_offset + 0.25

  p <- ggplot2::ggplot()
  for (nm in names(metrics)) {
    for (layer in make_curve_layers(nm)) {
      p <- p + layer
    }
  }

  p +
    ggforce::geom_arc_bar(
      data = pie_df,
      ggplot2::aes(
        x0 = x,
        y0 = y,
        r0 = 0,
        r = pie_r,
        amount = prop,
        fill = .data[[pie_taxrank]]
      ),
      stat = "pie",
      color = "white",
      linewidth = 0.3
    ) +
    ggforce::geom_circle(
      data = nodes,
      ggplot2::aes(x0 = x, y0 = y, r = pie_r),
      color = pie_border_color,
      fill = NA,
      linewidth = pie_border_width,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(
        x = x_label,
        y = y_label,
        label = node_label,
        hjust = hjust,
        vjust = vjust
      ),
      fontface = "bold",
      size = 4,
      lineheight = 0.9
    ) +
    ggplot2::annotate(
      "text",
      x = -1.65,
      y = legend_y0,
      label = "Metrics of similarity (higher = more similar)",
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = 3.2,
      color = "gray20"
    ) +
    ggplot2::geom_segment(
      data = legend_df,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      color = legend_df$color,
      linewidth = 1.8
    ) +
    ggplot2::geom_text(
      data = legend_df,
      ggplot2::aes(x = xend + 0.06, y = y, label = stats_label),
      hjust = 0,
      size = 3,
      color = "gray30"
    ) +
    (if (!is.null(sig_legend_df)) {
      list(
        ggplot2::annotate(
          "text",
          x = -1.65,
          y = perm_header_y,
          label = "Permutations",
          hjust = 0,
          vjust = 0.5,
          fontface = "bold",
          size = 3.2,
          color = "gray20"
        ),
        ggplot2::annotate(
          "text",
          x = -1.65,
          y = perm_desc_y,
          label = perm_desc,
          hjust = 0,
          vjust = 0.5,
          fontface = "italic",
          size = 2.6,
          color = "gray40"
        ),
        ggplot2::geom_segment(
          data = dplyr::filter(sig_legend_df, sig),
          ggplot2::aes(x = -1.65, xend = -1.35, y = y, yend = y),
          color = "grey30",
          linewidth = 1.8,
          alpha = 0.75
        ),
        ggplot2::geom_segment(
          data = dplyr::filter(sig_legend_df, !sig),
          ggplot2::aes(x = -1.65, xend = -1.35, y = y, yend = y),
          color = "grey30",
          linewidth = 1.8,
          alpha = nonsig_alpha
        ),
        ggplot2::geom_text(
          data = sig_legend_df,
          ggplot2::aes(x = -1.29, y = y, label = label),
          hjust = 0,
          size = 3,
          color = "gray30"
        )
      )
    } else {
      list()
    }) +
    fill_scale +
    ggplot2::coord_equal(xlim = c(-1.9, 1.9), ylim = c(y_bottom, y_top)) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold", size = 9),
      legend.text = ggplot2::element_text(size = 8),
      legend.key.size = ggplot2::unit(0.4, "cm")
    )
}

#' Companion bar chart for `community_sharing_pq()`
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Computes the same pairwise metrics as [community_sharing_pq()] and displays
#' them as grouped bars. Useful for precise numerical comparison alongside the
#' network figure.
#'
#' @param physeq A `phyloseq` object.
#' @param fact Character. Name of a `sample_data(physeq)` column. Must have
#'   2 to 4 unique values.
#' @param metrics Named list of metric definitions from [make_sharing_metric()].
#'   Default: [default_sharing_metrics()].
#' @param facet_by Character. One of `"metric"` (default: one panel per metric,
#'   x-axis = pair label) or `"pair"` (one panel per pair, x-axis = metric
#'   label). Bar fill colors come from each metric's `color` field.
#' @param show_na_modality Logical. Same meaning as in [community_sharing_pq()].
#'   Default `FALSE`.
#' @param base_size Base font size. Default `12`.
#' @param title Plot title. Default `NULL`.
#'
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' # One panel per metric
#' community_sharing_barplot_pq(data_fungi, fact = "Height")
#'
#' # One panel per pair
#' community_sharing_barplot_pq(data_fungi, fact = "Height", facet_by = "pair")
#' }
#' @export
#' @author Adrien Taudière
community_sharing_barplot_pq <- function(
  physeq,
  fact,
  metrics = default_sharing_metrics(),
  facet_by = c("metric", "pair"),
  show_na_modality = FALSE,
  base_size = 12,
  title = NULL
) {
  facet_by <- match.arg(facet_by)

  sd <- as.data.frame(phyloseq::sample_data(physeq))
  if (!fact %in% colnames(sd)) {
    stop("Column '", fact, "' not found in sample_data.")
  }
  fact_values <- sd[[fact]]
  if (show_na_modality && anyNA(fact_values)) {
    fact_values <- ifelse(is.na(fact_values), "NA", as.character(fact_values))
    sd[[fact]] <- fact_values
    phyloseq::sample_data(physeq) <- phyloseq::sample_data(sd)
  }
  modalities <- if (is.factor(fact_values)) {
    levels(droplevels(fact_values))
  } else {
    sort(unique(stats::na.omit(fact_values)))
  }
  n_mod <- length(modalities)
  if (n_mod < 2 || n_mod > 4) {
    stop("Number of modalities must be between 2 and 4, got ", n_mod, ".")
  }

  pairs_df <- .compute_pairs_df(physeq, fact, modalities, metrics)$pairs_df

  metric_colors <- stats::setNames(
    vapply(metrics, \(m) m$color, character(1)),
    names(metrics)
  )
  metric_labels <- stats::setNames(
    vapply(metrics, \(m) m$label, character(1)),
    names(metrics)
  )

  plot_df <- pairs_df |>
    dplyr::mutate(pair = paste0(from, " \u2014 ", to)) |>
    tidyr::pivot_longer(
      dplyr::all_of(names(metrics)),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric_label = factor(metric_labels[metric], levels = metric_labels)
    )

  color_scale <- ggplot2::scale_fill_manual(
    values = stats::setNames(metric_colors, metric_labels[names(metric_colors)])
  )

  base_theme <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )

  if (facet_by == "metric") {
    ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = pair, y = value, fill = metric_label)
    ) +
      ggplot2::geom_col(alpha = 0.85) +
      ggplot2::facet_wrap(~metric_label, scales = "free_y") +
      color_scale +
      ggplot2::labs(x = NULL, y = "Metric value", title = title) +
      base_theme
  } else {
    ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = metric_label, y = value, fill = metric_label)
    ) +
      ggplot2::geom_col(alpha = 0.85) +
      ggplot2::facet_wrap(~pair, scales = "free_y") +
      color_scale +
      ggplot2::labs(x = NULL, y = "Metric value", title = title) +
      base_theme
  }
}
