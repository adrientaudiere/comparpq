#' Add fake sequences by shuffling existing ones in a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful to compute true negative values with functions
#'   [tc_metrics_mock()] and [tc_metrics_mock_vec()].
#'
#' Note the that the tax_table for additional sequences is full of NA and
#'   that the corresponding otu_table is full of 0.
#'
#' @inheritParams tc_points_matrix
#' @param n_fake A number of fake taxa added.
#' @param prop_fake A proportion of fake taxa added.
#' @param prefix A prefix add to the taxa name
#' @returns A phyloseq object
#'
#' @export
#' @author Adrien Taudière
#' @examples
#' d_fake_F <- data_fungi_mini |>
#'   add_shuffle_seq_pq(prop_fake = 0.1)
#' ntaxa(d_fake_F) - ntaxa(data_fungi_mini)
add_shuffle_seq_pq <- function(physeq,
                               n_fake = NULL,
                               prop_fake = NULL,
                               prefix = "fake_") {
  taxasrow <- taxa_are_rows(physeq)
  if (taxasrow) {
    physeq <- MiscMetabar::taxa_as_columns(physeq)
  }
  if (is.null(n_fake) && is.null(prop_fake)) {
    stop("You must specify either n_fake or prop_fake param.")
  }
  if (!is.null(n_fake) && !is.null(prop_fake)) {
    stop("You must specify either n_fake or prop_fake param, not both !")
  }
  if (is.null(n_fake) && !is.null(prop_fake)) {
    n_fake <- prop_fake * ntaxa(physeq)
  }
  taxa_to_subset <- taxa_names(physeq) %in% sample(taxa_names(physeq), n_fake)
  names(taxa_to_subset) <- taxa_names(physeq)

  subset_physeq <- subset_taxa_pq(physeq, taxa_to_subset)
  refseq <- subset_physeq@refseq
  fake_refseq <- Biostrings::DNAStringSet(lapply(refseq, sample))
  subset_physeq@refseq <- refseq(fake_refseq)

  taxa_names(subset_physeq) <- paste0(prefix, 1:n_fake)
  new_physeq <- merge_phyloseq(physeq, subset_physeq)

  new_tax_tab <- as.matrix(unclass(new_physeq@tax_table))
  new_tax_tab[taxa_names(subset_physeq), ] <- NA
  new_physeq@tax_table <- tax_table(new_tax_tab)

  new_otu_tab <- as.matrix(unclass(new_physeq@otu_table))
  new_otu_tab[, taxa_names(subset_physeq)] <- 0
  new_physeq@otu_table <- otu_table(new_otu_tab, taxa_are_rows = FALSE)
  if (taxasrow) {
    newphyseq <- taxa_as_rows(new_physeq)
  }
  return(new_physeq)
}

#' Add external sequences to a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful to compute true negative values with functions
#'   [tc_metrics_mock()] and [tc_metrics_mock_vec()].
#'
#'   Note the that the tax_table for additional sequences is full of NA and
#'   that the corresponding otu_table is full of 0.
#'
#' @inheritParams tc_points_matrix
#' @param ext_seqs A DNAStringSet object
#' @param prefix A prefix add to the taxa name
#' @returns A phyloseq object
#'
#' @export
#' @author Adrien Taudière
#' @examples
#' add_external_seq_pq(
#'   data_fungi_mini,
#'   Biostrings::readDNAStringSet(
#'     system.file("extdata/ex_little.fasta", package = "MiscMetabar")
#'   )
#' )
#'
add_external_seq_pq <- function(physeq,
                                ext_seqs,
                                prefix = "external_") {
  refseq <- physeq@refseq
  names(ext_seqs) <- paste0(prefix, names(ext_seqs))
  external_samtab <- physeq@sam_data

  external_otu_tab <- otu_table(matrix(0,
    nrow = nsamples(physeq), ncol =
      length(ext_seqs)
  ), taxa_are_rows = FALSE)
  colnames(external_otu_tab) <- names(ext_seqs)
  rownames(external_otu_tab) <- sample_names(physeq)

  fake_tax_tab <- matrix("NA",
    nrow = length(ext_seqs),
    ncol = ncol(physeq@tax_table)
  )
  colnames(fake_tax_tab) <- colnames(physeq@tax_table)
  rownames(fake_tax_tab) <- names(ext_seqs)

  fake_pq <- phyloseq(
    refseq(ext_seqs),
    otu_table(external_otu_tab),
    sample_data(external_samtab),
    tax_table(fake_tax_tab)
  )

  new_physeq <- merge_phyloseq(physeq, fake_pq)
  new_tax_tab <- as.matrix(unclass(new_physeq@tax_table))
  new_tax_tab[taxa_names(fake_pq), ] <- NA
  new_physeq@tax_table <- tax_table(new_tax_tab)

  return(new_physeq)
}
