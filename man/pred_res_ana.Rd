\name{pred_res_ana}
\alias{pred_res_ana}
\title{Predictive Residue Analysis}
\description{
This function performs predictive residue analysis using an MCMC-based approach.
}
\usage{
pred_res_ana(Data, theta_0_samp, theta_samp, theta_til_samp, 
             motif_len_w, motif_len_g, burnin, N)
}
\arguments{
  \item{Data}{A list containing observed sequence data and related information.}
  \item{theta_0_samp}{A parameter matrix for theta_0.}
  \item{theta_samp}{A parameter matrix for theta.}
  \item{theta_til_samp}{A parameter matrix for theta_til.}
  \item{motif_len_w}{The length of the W motif.}
  \item{motif_len_g}{The length of the G motif.}
  \item{burnin}{The number of burn-in iterations for MCMC sampling.}
  \item{N}{The total number of MCMC iterations.}
}
\value{
A data frame containing:
\item{est_W}{Estimated W values from the best MCMC iteration.}
\item{est_G}{Estimated G values from the best MCMC iteration.}
\item{prob_00}{Posterior probability of (W = 0, G = 0).}
\item{prob_01}{Posterior probability of (W = 0, G = 1).}
\item{prob_10}{Posterior probability of (W = 1, G = 0).}
\item{prob_11}{Posterior probability of (W = 1, G = 1).}
}
\examples{
# Example usage:
# Data <- list(seq = ..., W_obs = ..., G_obs = ...)
# result <- pred_res_ana(Data, theta_0_samp, theta_samp, theta_til_samp, 
#                        motif_len_w, motif_len_g, burnin, N)
}
