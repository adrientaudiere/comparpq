% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_taxo.R
\name{tc_metrics_mock_vec}
\alias{tc_metrics_mock_vec}
\title{Compute accuracy metrics of taxonomic assignation using a mock (known)
community for one rank}
\usage{
tc_metrics_mock_vec(
  physeq,
  taxonomic_rank,
  true_values,
  fake_taxa = TRUE,
  fake_pattern = c("^fake_", "^external_"),
  verbose = TRUE
)
}
\arguments{
\item{physeq}{(required): a \code{\link[phyloseq]{phyloseq-class}} object obtained
using the \code{phyloseq} package.}

\item{taxonomic_rank}{(required) Name (or number) of a taxonomic rank
to count.}

\item{true_values}{(required) A vector with the true taxonomic assignation}

\item{fake_taxa}{(logical, default TRUE). If TRUE, the fake_pattern vector
is used to identify fake taxa, i.e. taxa who are not in the reference
database (see \code{\link[=add_external_seq_pq]{add_external_seq_pq()}}) or taxa with fake sequences
(see \code{\link[=add_shuffle_seq_pq]{add_shuffle_seq_pq()}}).}

\item{fake_pattern}{(vector of pattern): A vector used to identify the fake
taxa using a regex search in their name.}
}
\value{
A list of metrics (see the confusion matrix
\href{https://en.wikipedia.org/wiki/Confusion_matrix}{article} on wikipedia):
\itemize{
\item TP (number of \emph{true positive})
\item FP (number of \emph{false positive})
\item FN (number of \emph{false negative})
\item FDR (\emph{false discovery rate}) = FP / (FP + TP)
\item TPR (\emph{true positive rate}, also named \emph{recall} or \emph{sensitivity})
\item PPV (\emph{positive predictive value}, also named \emph{precision}) = TP / (TP + FP)
\item F1_score (\emph{F1 score}) = 2 * TP / (2 * TP + FP + FN)
}

If fake taxa are present and \code{fake_taxa} is true, other metrics are computed:
\itemize{
\item TN (number of \emph{true negative})
\item ACC (\emph{Accuracy}) = (TP + TN) / (TP + TN + FP + FN)
\item MCC (\emph{Matthews correlation coefficient}) = (TP * TN - FP * FN) /
sqrt((TP + FP) * (TP + FN) * (FP + TN) * (TN + FN))
}
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
\code{\link[=tc_metrics_mock]{tc_metrics_mock()}}, \code{\link[=add_external_seq_pq]{add_external_seq_pq()}}, \code{\link[=add_shuffle_seq_pq]{add_shuffle_seq_pq()}})
}
\author{
Adrien Taudière
}
