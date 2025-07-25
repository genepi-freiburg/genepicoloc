% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocalization.R
\name{create_result_template}
\alias{create_result_template}
\title{Create a template for colocalization results with gene annotation}
\usage{
create_result_template(
  CHR_var,
  BP_START_var,
  BP_STOP_var,
  sumstats_1_file,
  sumstats_2_file,
  sumstats_1_tabix,
  sumstats_2_tabix,
  sumstats_1_max_nlog10P = NA,
  sumstats_2_max_nlog10P = NA,
  nsnps = NA,
  sumstats_1_QC = NA,
  sumstats_2_QC = NA,
  coloc = NA
)
}
\arguments{
\item{CHR_var}{Character or numeric. Chromosome identifier.}

\item{BP_START_var}{Numeric. Start position of the region.}

\item{BP_STOP_var}{Numeric. End position of the region.}

\item{sumstats_1_file}{Character. Path to primary summary statistics file.}

\item{sumstats_2_file}{Character. Path to secondary summary statistics file.}

\item{sumstats_1_tabix}{Character. Tabix status for primary data.}

\item{sumstats_2_tabix}{Character. Tabix status for secondary data.}

\item{sumstats_1_max_nlog10P}{Numeric. Maximum -log10(p) in primary data.}

\item{sumstats_2_max_nlog10P}{Numeric. Maximum -log10(p) in secondary data.}

\item{nsnps}{Numeric. Number of SNPs in the analysis.}

\item{sumstats_1_QC}{Character. QC status of primary data.}

\item{sumstats_2_QC}{Character. QC status of secondary data.}

\item{coloc}{Character. Colocalization analysis status.}
}
\value{
A data.frame with one row containing all metadata fields and
placeholder NA values for results to be filled later.
}
\description{
Creates a standardized data.frame template for storing colocalization
results. This ensures consistent output format across all analyses.
}
\details{
The template includes fields for:
\itemize{
\item Region coordinates and metadata
\item Gene annotation (nearest genes, region center annotation)
\item File paths and data quality indicators
\item Colocalization results (PP.H0-H4)
\item Directionality analysis results
\item Top SNP information
}

This standardized format facilitates downstream analysis and reporting.
}
\seealso{
\code{\link{run_region}} which creates and populates this template
\code{\link{format_coloc_output}} which adds results to the template
}
