% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocalization.R
\name{get_directionality}
\alias{get_directionality}
\title{Analyze effect directionality between datasets at lead SNP}
\usage{
get_directionality(sumstats_1, sumstats_2)
}
\arguments{
\item{sumstats_1}{A sumstats object containing the primary dataset.}

\item{sumstats_2}{A sumstats object containing the secondary dataset.}
}
\value{
A data.frame with one row containing:
\itemize{
\item sumstats_1_ind_Name: SNP identifier of lead variant
\item sumstats_1_ind_A1/A2: Alleles in primary dataset
\item sumstats_1_ind_nlog10P: Significance of lead variant
\item sumstats_1_ind_BETA: Effect size in primary dataset
\item sumstats_1_ind_BETA_sign: Sign of effect (+1 or -1)
\item sumstats_2_ind_A1/A2: Alleles in secondary dataset
\item sumstats_2_ind_nlog10P: Significance in secondary dataset
\item sumstats_2_ind_BETA: Effect size in secondary dataset
\item sumstats_2_ind_BETA_sign: Sign of effect (+1 or -1)
\item directionality: 1 if same direction, -1 if opposite, NA if undetermined
}
}
\description{
Identifies the most significant SNP in the primary dataset and analyzes
whether the effects in both datasets have the same direction. This helps
determine if the colocalized signal represents the same or opposite effects.
}
\details{
The function:
\enumerate{
\item Finds the intersection of SNPs between datasets
\item Identifies the most significant SNP in the primary dataset
\item Extracts information for this SNP from both datasets
\item Determines effect directionality accounting for allele coding
}

Directionality is determined by comparing effect signs while accounting
for potential allele flips between datasets.
}
\seealso{
\code{\link{detect_directionality}} for the directionality logic
\code{\link{run_region}} which calls this function
}
