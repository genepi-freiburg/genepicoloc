% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/development.R
\name{query_GTEXv8_GWAS}
\alias{query_GTEXv8_GWAS}
\title{Query GTEx v8 GWAS}
\usage{
query_GTEXv8_GWAS(
  sumstats_file,
  CHR_var,
  BP_START_var,
  BP_STOP_var,
  phenotype_id_var
)
}
\arguments{
\item{sumstats_file}{path to GTEx v7 sumstats.}

\item{CHR_var}{CHR (as.character "1", "2", ..., "X").}

\item{BP_START_var}{start of region, integer}

\item{BP_STOP_var}{end of region, integer}

\item{annotation_file}{annotation file with number of samples for each GWAS.}
}
\value{
data frame with extracted sumstats.
}
\description{
Query GTEx v8 GWAS data to extract a region of interest
}