#' Package initialization
#'
#' @description Function called when the package is loaded. Registers S3 methods
#' to ensure proper method dispatch.
#'
#' @param libname Library path where package is installed
#' @param pkgname Package name
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Register S3 methods that might have dispatch issues
  registerS3method("print", "sumstats_2_form", print.sumstats_2_form)
}

#' Package startup message
#'
#' @description Function called when the package is attached. Displays package
#' version and citation information.
#'
#' @param libname Library path where package is installed
#' @param pkgname Package name
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  # Get package version
  version <- utils::packageDescription(pkgname, fields = "Version")
  
  # Create startup message
  packageStartupMessage(
    "genepicoloc v", version, " - R package to facilitate genetic colocalization analyses\n"
    # "For citation information, use: citation('", pkgname, "')"
  )
}

# Suppress R CMD check NOTEs about undefined global variables
# These are typically column names used with data.table syntax
utils::globalVariables(c(
  "..required_cols",
  "..std_cols",
  "CHR",
  "Name",
  "POS",
  "Phenotype",
  "copy",
  "effectAlleleFreq",
  "region",
  "study_id",
  "sumstats_2_study"
))