% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_query_sumstats.R
\name{query_UKB_GWAS}
\alias{query_UKB_GWAS}
\title{Query UKB GWAS}
\usage{
query_UKB_GWAS(sumstats_file, CHR_var, BP_START_var, BP_STOP_var, ...)
}
\arguments{
\item{sumstats_file}{path tabix-indexed sumstats.}

\item{CHR_var}{chromosome (as.character "1", "2", ..., "X").}

\item{BP_START_var}{start of region, integer}

\item{BP_STOP_var}{end of region, integer}
}
\value{
data frame with extracted sumstats.
}
\description{
Query UK Biobank GWAS data to extract a region of interest
}
