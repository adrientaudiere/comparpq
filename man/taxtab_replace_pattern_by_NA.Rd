% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxtab_modification.R
\name{taxtab_replace_pattern_by_NA}
\alias{taxtab_replace_pattern_by_NA}
\title{Replace taxonomic value with a given pattern by NA}
\usage{
taxtab_replace_pattern_by_NA(
  physeq,
  patterns = c(".*_incertae_sedis", "unclassified.*"),
  taxonomic_ranks = NULL,
  progress_bar = FALSE,
  ...
)
}
\arguments{
\item{physeq}{(required): a \code{\link[phyloseq]{phyloseq-class}} object obtained
using the \code{phyloseq} package.}

\item{patterns}{A vector of patterns to select taxonomic value to convert to NA}

\item{taxonomic_ranks}{A list of taxonomic ranks where we want the
substitution occur. If left to NULL, all taxonomic ranks are modified.}

\item{progress_bar}{(logical, default FALSE) Do we print progress during
the calculation?}

\item{...}{Other arguments to be passed on to \code{\link[base:grep]{base::gsub()}}}
}
\value{
A phyloseq object
}
\description{
\if{html}{\out{
<a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
<img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
}}


Use \code{\link[base:grep]{base::gsub()}} for substitution.
}
\examples{
data_fungi@tax_table["ASV85", "Family"]
data_fungi2 <- taxtab_replace_pattern_by_NA(data_fungi, "fam_Incertae_sedis", taxonomic_ranks = "Family")
data_fungi2@tax_table["ASV85", "Family"]

# By default patterns ".*_incertae_sedis" and "unclassified.*" are replaced by NA
data_fungi3 <-
  taxtab_replace_pattern_by_NA(data_fungi, ignore.case = TRUE, progress_bar = TRUE)
data_fungi3@tax_table["ASV85", "Family"]
}
\author{
Adrien Taudière
}
