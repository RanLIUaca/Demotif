% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/res_ana_RL_MH_FUN.R
\name{res_ana}
\alias{res_ana}
\title{res_ana}
\usage{
res_ana(motif_len_w, motif_len_g, burnin, N)
}
\arguments{
\item{motif_len_w}{The length of the first binding motif.}

\item{motif_len_g}{The length of the second binding motif.}

\item{burnin}{The number of burn-in.}

\item{N}{The number of sampling iterations.}
}
\value{
A list containing:
\itemize{
\item Samples for \code{theta_0} (background probabilities),
\item Samples for \code{theta} (the first binding motif probabilities),
\item Samples for \verb{tilde\{theta\}} (the second binding motif probabilities),
\item Samples for \code{W} and \code{G} (indicators for motif presence),
\item Samples for \code{A} and \code{B} (positions of motifs),
\item Log-likelihood values,
\item Estimated parameters (\code{theta_0}, \code{theta}, \verb{tilde\{theta\}}),
\item The number of jumps for both motifs,
\item a plot for the curve of the loglikelihood, the estimated theta logo and the estimated tilde{theta} logo,
\item a plot for the logos before and after a jump.
@export
}
}
\description{
This function performs the analysis of motif estimation and sampling results. It includes reading the sampled data, estimating parameters, plotting results, and analyzing motif jumps.
}
