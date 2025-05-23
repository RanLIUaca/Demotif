% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_demotif_FUN.R
\name{predict_demotif}
\alias{predict_demotif}
\title{Predict Motif}
\usage{
predict_demotif(
  Data_name,
  Data_one_row,
  dict,
  theta_0,
  theta,
  theta_til,
  motif_len_w,
  motif_len_g,
  burnin,
  num_max
)
}
\arguments{
\item{Data_name}{A string representing the dataset name.}

\item{Data_one_row}{A dataframe containing one row of sequence data with columns:
\enumerate{
\item \code{seq} (sequence string)
\item \code{W_obs} (observed W values, may contain NAs)
\item \code{G_obs} (observed G values, may contain NAs)
}}

\item{dict}{A dictionary containing motif-related information.}

\item{theta_0}{Estimated parameter for theta_0.}

\item{theta}{Estimated parameter for theta.}

\item{theta_til}{Estimated parameter for theta_til.}

\item{motif_len_w}{Length of motif W.}

\item{motif_len_g}{Length of motif G.}

\item{burnin}{Number of initial iterations to discard as burn-in.}

\item{num_max}{Total number of MCMC iterations.}
}
\value{
A list containing:
\itemize{
\item \code{res_frame}: A dataframe with predicted motif states and probabilities.
\item \code{est_A}: Estimated start and end positions for motif W.
\item \code{est_B}: Estimated start and end positions for motif G.
}
}
\description{
This function predicts motif locations in a given sequence using MCMC sampling.
It processes the input data, generates samples using \code{generate_samp()},
reads the sampled data, and determines the most probable motif locations.
}
