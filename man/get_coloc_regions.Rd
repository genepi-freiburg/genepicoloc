% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_coloc.R
\name{get_coloc_regions}
\alias{get_coloc_regions}
\title{Get coloc regions}
\usage{
get_coloc_regions(
  sumstats,
  CHR_name = "CHR",
  POS_name = "POS",
  p_value_name = "P",
  p_threshold = 5e-08,
  halfwindow = 5e+05
)
}
\arguments{
\item{sumstats}{data frame read with read_sumstats_1(). Mandatory columns: CHR, BP, P}

\item{p_threshold}{search for regions until no more variants below this threshold remains}

\item{log_name}{iteration log will be written to this file}
}
\value{
data frame with extracted regions
}
\description{
Find significant regions in sumstats and merge close regions if needed
from a list of external sumstats.
}
\examples{
under development
}
