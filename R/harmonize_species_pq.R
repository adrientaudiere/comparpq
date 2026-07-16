################################################################################
#' Extract the species epithet from a Species value
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Returns the specific epithet from a value of a `Species` column, handling the
#' many formats found in reference databases: bare epithet (`"iris"`), binomial
#' with a space or an underscore (`"Genus iris"`, `"Genus_iris"`), and trailing
#' infraspecific parts (`"Genus_iris_var._occidentalis"`). Underscores are
#' treated as spaces. When the value holds only the genus (matching `genus_val`),
#' `NA` is returned.
#'
#' @param species_val (character) A single value of the `Species` column.
#' @param genus_val (character) The corresponding `Genus` value (may be `NA`),
#'   used to decide whether the first word is the genus.
#'
#' @returns The species epithet as a length-one character, or `NA_character_`.
#' @export
#' @author Adrien Taudière
#'
#' @seealso [harmonize_sp_names_pq()]
#'
#' @examples
#' extract_species_epithet("Quercus_ilex", "Quercus")
#' extract_species_epithet("Quercus ilex var. rotundifolia", "Quercus")
#' extract_species_epithet("ilex", "Quercus")
extract_species_epithet <- function(species_val, genus_val) {
  if (is.na(species_val) || species_val == "") {
    return(NA_character_)
  }

  s <- gsub("_", " ", species_val)
  s <- trimws(s)
  if (s == "") {
    return(NA_character_)
  }

  words <- strsplit(s, "\\s+")[[1]]

  if (length(words) == 1) {
    if (!is.na(genus_val) && tolower(words[1]) == tolower(genus_val)) {
      return(NA_character_)
    }
    return(words[1])
  }

  if (!is.na(genus_val) && tolower(words[1]) == tolower(genus_val)) {
    return(words[2])
  }

  if (grepl("^[A-Z]", words[1])) {
    return(words[2])
  }

  words[1]
}

################################################################################
#' Harmonize Species_<db> columns and verify names via GNA Verifier
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' For each database suffix (e.g. `"_EUK"`, `"_Unite"`), (1) rewrites the
#' `Species<suffix>` column to keep only the specific epithet (see
#' [extract_species_epithet()]) and, when `verify = TRUE`, (2) calls
#' [taxinfo::gna_verifier_pq()] on `c("Genus<suffix>", "Species<suffix>")` and
#' appends the resulting columns (`currentCanonicalSimple<suffix>`,
#' `currentName<suffix>`, ...) to the `tax_table`, renamed from
#' `<db>_<col>` to `<col>_<db>` for consistency with the `<Rank>_<Db>` naming.
#'
#' The `verify = TRUE` path requires the suggested package \pkg{taxinfo} and
#' network access; the `verify = FALSE` path performs only the offline epithet
#' harmonization and needs no extra package.
#'
#' @param physeq (required) A [phyloseq::phyloseq-class()] object.
#' @param suffixes (character vector, required) Database suffixes
#'   (e.g. `c("_EUK", "_Unite")`). A leading underscore is added when missing.
#' @param verify (logical, default `TRUE`) When `FALSE`, only epithet
#'   harmonization is performed, without calling [taxinfo::gna_verifier_pq()].
#' @param force_recompute (logical, default `FALSE`) When `TRUE`, drops existing
#'   result columns from a previous run before re-adding them; also forwarded to
#'   [taxinfo::gna_verifier_pq()].
#' @param ... Additional arguments passed to [taxinfo::gna_verifier_pq()]
#'   (e.g. `data_sources = 5` for Index Fungorum).
#'
#' @returns The `physeq` object with harmonized `Species<suffix>` columns and,
#'   when `verify = TRUE`, the new `currentCanonicalSimple<suffix>`,
#'   `currentName<suffix>`, ... columns.
#' @export
#' @author Adrien Taudière
#'
#' @seealso [extract_species_epithet()], [taxinfo::gna_verifier_pq()]
#'
#' @examples
#' # Offline: keep only the epithet in a Species_<db> column.
#' pq <- Glom_otu
#' mat <- as(pq@tax_table, "matrix")
#' mat <- cbind(
#'   mat,
#'   Genus_EUK = rep("Quercus", nrow(mat)),
#'   Species_EUK = rep("Quercus_ilex", nrow(mat))
#' )
#' pq@tax_table <- phyloseq::tax_table(mat)
#' pq <- harmonize_sp_names_pq(pq, suffixes = "_EUK", verify = FALSE)
#'
#' \dontrun{
#' # verify = TRUE additionally validates names via taxinfo::gna_verifier_pq()
#' # (requires the 'taxinfo' package and network access).
#' pq <- harmonize_sp_names_pq(pq, suffixes = "_EUK", data_sources = 5)
#' }
harmonize_sp_names_pq <- function(
  physeq,
  suffixes,
  verify = TRUE,
  force_recompute = FALSE,
  ...
) {
  verify_pq(physeq)

  suffixes <- vapply(
    suffixes,
    function(s) {
      if (!startsWith(s, "_")) {
        paste0("_", s)
      } else {
        s
      }
    },
    character(1)
  )

  # --- Step 1: harmonize the Species<suffix> columns ---
  for (suf in suffixes) {
    species_col <- paste0("Species", suf)
    genus_col <- paste0("Genus", suf)

    if (!(species_col %in% colnames(physeq@tax_table))) {
      next
    }

    species_vals <- as.character(physeq@tax_table[, species_col])
    genus_vals <- if (genus_col %in% colnames(physeq@tax_table)) {
      as.character(physeq@tax_table[, genus_col])
    } else {
      rep(NA_character_, length(species_vals))
    }

    physeq@tax_table[, species_col] <- mapply(
      extract_species_epithet,
      species_vals,
      genus_vals,
      USE.NAMES = FALSE
    )
  }

  # --- Step 2: GNA Verifier verification for each database ---
  if (!verify) {
    return(physeq)
  }

  if (!requireNamespace("taxinfo", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg taxinfo} is required when {.code verify = TRUE}.",
      "i" = "Install it, or call with {.code verify = FALSE} for offline epithet harmonization only."
    ))
  }

  for (suf in suffixes) {
    genus_col <- paste0("Genus", suf)
    species_col <- paste0("Species", suf)

    if (
      !(genus_col %in%
        colnames(physeq@tax_table) &&
        species_col %in% colnames(physeq@tax_table))
    ) {
      next
    }

    db_name <- gsub("^_", "", suf)
    prefix <- paste0(db_name, "_")

    if (force_recompute) {
      result_cols <- c(
        "currentName",
        "currentCanonicalSimple",
        "genusEpithet",
        "specificEpithet",
        "genusSpeciesEpithet",
        "namePublishedInYear",
        "authorship",
        "bracketauthorship",
        "scientificNameAuthorship"
      )
      stale_cols <- paste0(result_cols, "_", db_name)
      stale_cols <- intersect(stale_cols, colnames(physeq@tax_table))
      if (length(stale_cols) > 0) {
        tax_mat <- as(physeq@tax_table, "matrix")
        tax_mat <- tax_mat[,
          !(colnames(tax_mat) %in% stale_cols),
          drop = FALSE
        ]
        physeq@tax_table <- phyloseq::tax_table(tax_mat)
      }
    }

    physeq <- taxinfo::gna_verifier_pq(
      physeq,
      taxonomic_rank = c(genus_col, species_col),
      col_prefix = prefix,
      force_recompute = force_recompute,
      ...
    )

    pattern <- paste0("^", db_name, "_")
    cols_to_rename <- grep(
      pattern,
      colnames(physeq@tax_table),
      value = TRUE
    )
    if (length(cols_to_rename) > 0) {
      new_names <- sub(pattern, "", cols_to_rename)
      new_names <- paste0(new_names, "_", db_name)
      idx <- match(cols_to_rename, colnames(physeq@tax_table))
      colnames(physeq@tax_table)[idx] <- new_names
    }
  }

  physeq
}
################################################################################
