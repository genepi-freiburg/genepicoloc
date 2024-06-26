% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_coloc.R
\name{process_wrapper}
\alias{process_wrapper}
\title{Process results of coloc_wrapper function.}
\usage{
process_wrapper(coloc_output, N_top_SNPs = 5, remove_full_results = T)
}
\arguments{
\item{coloc_output}{output of coloc_wrapper function.}

\item{N_top_SNPs}{Number of SNPs with highest PP.H4 to output.}

\item{remove_full_results}{Should data.frame with full coloc output be removed?
Usually TRUE, in this case only first SNPs are used in output.}
}
\value{
results of coloc.abf function
}
\description{
Process results of coloc_wrapper function.
}
