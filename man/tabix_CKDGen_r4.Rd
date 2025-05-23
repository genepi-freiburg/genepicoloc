% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query_sumstats_custom.R
\name{tabix_CKDGen_r4}
\alias{tabix_CKDGen_r4}
\title{Query CKDGen r4 summary statistics using tabix}
\usage{
tabix_CKDGen_r4(sumstats_file, coloc_regions_PASS)
}
\arguments{
\item{sumstats_file}{Path to the CKDGen r4 summary statistics file}

\item{coloc_regions_PASS}{Data frame with columns CHR_var, BP_START_var, BP_STOP_var}
}
\value{
Data frame with formatted summary statistics
}
\description{
Retrieves CKDGen r4 summary statistics for specified genomic regions
}
\keyword{internal}
