% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_samp_FUN.R
\name{generate_samp}
\alias{generate_samp}
\title{Generate Sample Data for Motif Analysis}
\usage{
generate_samp(
  Data_name,
  Data_pre,
  dict,
  theta_0,
  theta,
  theta_til,
  motif_len_w,
  motif_len_g,
  num_max
)
}
\arguments{
\item{Data_name}{A string indicating the dataset name.}

\item{Data_pre}{A data frame containing the input sequence (\code{seq}), observed \code{W} (\code{W_obs}), and observed \code{G} (\code{G_obs}).}

\item{dict}{A dictionary containing motif-related information.}

\item{theta_0}{Estimated parameter for theta_0.}

\item{theta}{Estimated parameter for theta.}

\item{theta_til}{Estimated parameter for theta_til.}

\item{motif_len_w}{The length of motif W.}

\item{motif_len_g}{The length of motif G.}

\item{num_max}{The number of sampling iterations.}
}
\value{
This function does not return a value but saves sampled data as CSV files in './result/temp_data'.
}
\description{
This function generates sample motif data by iteratively updating
the binding locations (A, B) and binding labels (W, G) using probabilistic sampling.
}
