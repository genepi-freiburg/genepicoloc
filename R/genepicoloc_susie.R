#' Perform SuSiE-RSS for fine-mapping using summary statistics and LD matrix
#'
#' @description
#' This function runs the SuSiE-RSS algorithm for fine-mapping genetic loci
#' using summary statistics and a linkage disequilibrium (LD) matrix.
#'
#' @param sumstats_data A data frame containing summary statistics for each SNP.
#'   Must include columns for SNP identifiers, effect sizes, standard errors,
#'   and sample sizes.
#' @param R A matrix of linkage disequilibrium (LD) values between SNPs.
#'   Rows and columns must correspond to the SNPs in `sumstats_data`.
#' @param snp_col_name The name of the column in `sumstats_data` containing SNP identifiers.
#'   Default: "Name".
#' @param bhat_col_name The name of the column in `sumstats_data` containing effect sizes.
#'   Default: "BETA".
#' @param shat_col_name The name of the column in `sumstats_data` containing standard errors.
#'   Default: "SE".
#' @param n_col_name The name of the column in `sumstats_data` containing sample sizes.
#'   Default: "N".
#' @param n The sample size to use for fine-mapping. If not provided, the median
#'   sample size from `sumstats_data` is used.
#' @param L Maximum number of components (causal SNPs) to include. Default: 10.
#' @param coverage The credible set coverage probability. Default: 0.95.
#' @param refine Logical, whether to refine the estimated PIPs. Default: FALSE.
#' @param max_iter Maximum number of iterations for the algorithm. Default: 100.
#' @param min_abs_corr Minimum absolute correlation allowed in a credible set. Default: 0.5.
#' @param tol Convergence tolerance for the IBSS fitting procedure. Default: 0.001.
#' @param ... Additional arguments passed to `susie_rss`.
#'
#' @details
#' The function checks that SNP names in `sumstats_data` match the rownames and colnames
#' of `R`. It then runs `susie_rss` with the provided parameters.
#'
#' @return
#' A fitted SuSiE-RSS object.
#'
#' @examples
#' # Example usage:
#' # sumstats_data <- data.frame(Name = c("rs1", "rs2"), BETA = c(0.1, 0.2),
#' #                            SE = c(0.01, 0.02), N = c(1000, 1000))
#' # R <- matrix(c(1, 0.8, 0.8, 1), nrow = 2,
#' #             dimnames = list(c("rs1", "rs2"), c("rs1", "rs2")))
#' # result <- genepicoloc_susie_rss(sumstats_data, R)
#'
#' @seealso
#' \code{\link{susie_rss}} for details on the underlying algorithm.
#'
#' @importFrom susieR susie_rss
#'
#' @export
genepicoloc_susie_rss <- function(
    sumstats_data,
    R,
    snp_col_name = "Name",
    bhat_col_name = "BETA",
    shat_col_name = "SE",
    n_col_name = "N",
    n,
    L = 10,
    coverage = 0.95,
    refine = FALSE,
    max_iter = 100,
    min_abs_corr = 0.5,
    tol = 0.001,
    ...
) {
    # Check that SNP names match between sumstats_data and R
    stopifnot(all(sumstats_data[[snp_col_name]] == colnames(R)))
    stopifnot(all(sumstats_data[[snp_col_name]] == rownames(R)))

    # Use median sample size if n is not provided
    if (missing(n)) {
        n <- median(sumstats_data[[n_col_name]], na.rm = TRUE)
    }
    if (is.na(n) || is.null(n)) {
        stop("Sample size is missing")
    }

    # Run SuSiE-RSS
    fitted_rss <- susie_rss(
        bhat = sumstats_data[[bhat_col_name]],
        shat = sumstats_data[[shat_col_name]],
        n = n,
        R = R,
        L = L,
        coverage = coverage,
        refine = refine,
        max_iter = max_iter,
        min_abs_corr = min_abs_corr,
        tol = tol,
        ...
    )

    return(fitted_rss)
}

