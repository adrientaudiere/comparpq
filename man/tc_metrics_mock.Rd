% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_taxo.R
\name{tc_metrics_mock}
\alias{tc_metrics_mock}
\title{Compute accuracy metrics of multiple taxonomic assignations method using
mock for multi-rank and multi assignation methods}
\usage{
tc_metrics_mock(
  physeq,
  ranks_df,
  true_values_df,
  fake_taxa = TRUE,
  fake_pattern = c("^fake_", "^external_")
)
}
\arguments{
\item{physeq}{(required): a \code{\link[phyloseq]{phyloseq-class}} object obtained
using the \code{phyloseq} package.}

\item{ranks_df}{(required). A dataframe with at least one column (one database or one
method) and a number of row equal to the column in true_values_df}

\item{true_values_df}{(required). A dataframe with the true taxonomic assignation.
Note that the column names (names of taxonomic ranks) of the true_values_df defined
the names present in the \code{tax_level} column of the resulting dataframe.}

\item{fake_taxa}{(logical, default TRUE). If TRUE, the fake_pattern vector
is used to identify fake taxa, i.e. taxa who are not in the reference
database (see \code{\link[=add_external_seq_pq]{add_external_seq_pq()}}) or taxa with fake sequences
(see \code{\link[=add_shuffle_seq_pq]{add_shuffle_seq_pq()}}).}

\item{fake_pattern}{(vector of pattern): A vector used to identify the fake
taxa using a regex search in their name.}
}
\value{
A long-format dataframe with 4 columns: (i) the name of the \code{method_db}
(ii) the name of the \code{tax_level} (taxonomic rank), (iii) the \code{metrics} (see
\code{\link[=tc_metrics_mock_vec]{tc_metrics_mock_vec()}} for more details) and (iv) the \code{values}.
}
\description{
\if{html}{\out{
<a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
<img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
}}


Compute numerous metrics comparing the computed taxonomic assignation to a
true assignation.

Note that to compute all metrics, one need to insert fake
taxa (by shuffling sequences and/or by adding external sequences). The user
must fake taxa using functions \code{\link[=add_external_seq_pq]{add_external_seq_pq()}},
\code{\link[=add_shuffle_seq_pq]{add_shuffle_seq_pq()}}) before taxonomic assignation.
}
\seealso{
\code{\link[=tc_metrics_mock_vec]{tc_metrics_mock_vec()}}
}
\author{
Adrien Taudière
}
