.onLoad <- function(libname, pkgname) {
  # Register S3 methods that might have dispatch issues
  registerS3method("print", "sumstats_2_form", print.sumstats_2_form)
}