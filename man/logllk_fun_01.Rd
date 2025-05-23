% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LOGLLK_FUN_2.R
\name{logllk_fun_01}
\alias{logllk_fun_01}
\title{logllk_fun_01}
\usage{
logllk_fun_01(Ri, bi, theta_0_samp, theta_til_samp, dict)
}
\arguments{
\item{Ri}{A string representing the sequence.}

\item{bi}{The positions of the second motif within the sequence.}

\item{theta_0_samp}{A vector of sampled probabilities for the background model.}

\item{theta_til_samp}{A matrix of sampled probabilities for the second motif.}

\item{dict}{A dictionary of unique characters in the sequence.}
}
\value{
A numeric value representing the log-likelihood of the sequence.
}
\description{
This function calculates the log-likelihood of a sequence when the second motif is present (gi = 1) and the first motif is absent (wi = 0).
}
