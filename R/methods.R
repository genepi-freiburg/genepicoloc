#' Print method for sumstats_2_form objects
#'
#' @description Displays a summary of the sumstats_2_form object including
#' the number of datasets and phenotype type.
#'
#' @param x A sumstats_2_form object
#' @param ... Additional arguments (ignored)
#' @return The input object (invisibly)
#' @export
print.sumstats_2_form <- function(x, ...) {
  cat("sumstats_2_form object with", length(x), "datasets\n")
  cat("sumstats_pheno type:", attr(x, "sumstats_pheno"), "\n")
  invisible(x)
}
