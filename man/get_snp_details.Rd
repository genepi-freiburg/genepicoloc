% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebi_interface.R
\name{get_snp_details}
\alias{get_snp_details}
\title{Get genomic coordinates for an rsID}
\usage{
get_snp_details(rsid)
}
\arguments{
\item{rsid}{Character string. The rsID to look up (e.g., "rs12345").}
}
\value{
Character string in "chr:position" format, or NA if lookup fails.
}
\description{
Queries the GWAS catalog API to retrieve chromosomal coordinates
for a given rsID.
}
\details{
This helper function is used internally by get_gwas_associations to
convert rsIDs to genomic coordinates. It handles compound rsIDs
(e.g., "rs12345-A") by extracting just the rsID portion.
}
