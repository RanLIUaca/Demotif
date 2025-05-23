% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/THETA_MH_FUN_2.R
\name{theta_samp_fun}
\alias{theta_samp_fun}
\title{theta_samp_fun}
\usage{
theta_samp_fun(
  R,
  W_samp,
  G_samp,
  A_samp,
  B_samp,
  motif_len_w,
  motif_len_g,
  dict,
  pri_alw
)
}
\arguments{
\item{R}{A vector of string sequences.}

\item{W_samp}{A vector indicating the presence of the first motif (sampled).}

\item{G_samp}{A vector indicating the presence of the second motif (sampled).}

\item{A_samp}{A vector of positions for the first motif (sampled).}

\item{B_samp}{A vector of positions for the second motif (sampled).}

\item{motif_len_w}{Length of the first motif.}

\item{motif_len_g}{Length of the second motif.}

\item{dict}{A dictionary of possible characters in the sequences.}

\item{pri_alw}{A prior Dirichlet parameter for the probabilities of the first motif.}
}
\value{
A list containing:
\itemize{
\item \code{theta_samp}: A matrix of sampled probabilities for each position in the first motif.
\item \code{Post_alpha_w}: A matrix of posterior Dirichlet parameters for each position in the first motif.
}
}
\description{
theta_samp_fun
}
