% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query_sumstats_custom.R
\name{tabix_GCKD_pGWAS}
\alias{tabix_GCKD_pGWAS}
\title{Query GCKD pGWAS summary statistics using tabix}
\usage{
tabix_GCKD_pGWAS(sumstats_file, coloc_regions_PASS)
}
\arguments{
\item{sumstats_file}{Path to the GCKD pGWAS summary statistics file}

\item{coloc_regions_PASS}{Data frame with columns CHR_var, BP_START_var, BP_STOP_var}
}
\value{
Data frame with formatted summary statistics
}
\description{
Retrieves GCKD pGWAS summary statistics for specified genomic regions
}
\keyword{internal}
