% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_coloc.R
\name{create_coloc_params_df}
\alias{create_coloc_params_df}
\title{Create data.frame with input parameters for coloc}
\usage{
create_coloc_params_df(
  sumstats_1_args,
  EXPERIMENT,
  sumstats_path,
  files = NULL,
  sumstats_pattern = "gz$",
  grep_invert = NULL,
  sumstats_function,
  hyprcoloc = F,
  sumstats_type,
  sumstats_sdY = NA,
  extra_args = NULL,
  do_annotate = F,
  annotation_function,
  annotation_function_args,
  do_annotate_sumstats_1 = F,
  annotation_function_sumstats_1,
  annotation_function_args_sumstats_1
)
}
\arguments{
\item{sumstats_path}{path to folder with indexed sumstats.}

\item{sumstats_pattern}{usually "gz$".}

\item{grep_invert}{expression to exclude some found files
(e.g., annotation files in the same folder).}

\item{sumstats_function}{function name to query sumstats.}

\item{sumstats_type}{quant or cc.}

\item{sumstats_sdY}{sdY for quant traits,
by default NA and estimated from BETA and MAF.}
}
\value{
data.frame with input parameters for coloc
}
\description{
Create data.frame with input parameters for coloc
}
\examples{
create_coloc_params_df()
}
