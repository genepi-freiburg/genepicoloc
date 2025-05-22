#' @export
print.sumstats_2_form <- function(x, ...) {
  cat("sumstats_2_form object with", length(x), "datasets\n")
  cat("sumstats_pheno type:", attr(x, "sumstats_pheno"), "\n")
  invisible(x)
}
