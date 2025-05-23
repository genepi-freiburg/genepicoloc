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