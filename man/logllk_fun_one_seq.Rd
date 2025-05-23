% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LOGLLK_FUN_2.R
\name{logllk_fun_one_seq}
\alias{logllk_fun_one_seq}
\title{logllk_fun_one_seq}
\usage{
logllk_fun_one_seq(
  Ri,
  wi,
  gi,
  ai,
  bi,
  dict,
  theta_0_samp,
  theta_samp,
  theta_til_samp
)
}
\arguments{
\item{Ri}{A string representing the sequence.}

\item{wi}{Indicator for the first motif (1 if present, 0 otherwise).}

\item{gi}{Indicator for the second motif (1 if present, 0 otherwise).}

\item{ai}{The positions of the first motif within the sequence.}

\item{bi}{The positions of the second motif within the sequence.}

\item{dict}{A dictionary of unique characters in the sequence.}

\item{theta_0_samp}{A vector of sampled probabilities for the background model.}

\item{theta_samp}{A matrix of sampled probabilities for the first motif.}

\item{theta_til_samp}{A matrix of sampled probabilities for the second motif.}
}
\value{
A numeric value representing the log-likelihood of the sequence.
}
\description{
This function calculates the log-likelihood of a single sequence based on the presence or absence of motifs.
}
