# Approximate Bayes Factor Colocalization
#
# Clean-room implementation of the ABF colocalization method from:
#   Giambartolomei et al. (2014) PLoS Genetics 10(5): e1004383
#   Wallace (2020) PLoS Genetics 16(4): e1008720
#
# Computes posterior probabilities for 5 hypotheses:
#   H0: Neither trait associated
#   H1: Only trait 1 associated
#   H2: Only trait 2 associated
#   H3: Both traits associated, different causal variants
#   H4: Both traits associated, shared causal variant


#' Log-sum-exp: log(sum(exp(x))) with numerical stability
#' @param x Numeric vector of log values
#' @return Scalar: log(sum(exp(x)))
#' @keywords internal
logsum <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}


#' Log-difference-exp: log(exp(a) - exp(b)) with numerical stability
#' Requires a > b (i.e., exp(a) > exp(b))
#' @param a,b Scalars (log-scale)
#' @return Scalar: log(exp(a) - exp(b))
#' @keywords internal
logdiff <- function(a, b) {
  m <- max(a, b)
  m + log(exp(a - m) - exp(b - m))
}


#' Compute per-SNP log approximate Bayes factors from beta/se
#'
#' @param beta Numeric vector of effect sizes
#' @param se Numeric vector of standard errors
#' @param prior_var Numeric scalar: prior variance on true effect size (W)
#' @return Numeric vector of log ABFs
#' @keywords internal
compute_labf <- function(beta, se, prior_var) {
  V <- se^2
  z <- beta / se
  r <- prior_var / (prior_var + V)
  0.5 * (log(1 - r) + r * z^2)
}


#' Estimate trait standard deviation from GWAS summary statistics
#'
#' Uses the relationship between beta, MAF, N, and trait variance
#' to estimate sdY when not provided directly.
#' From Wakefield (2009) approximation.
#'
#' @param beta Numeric vector of effect sizes
#' @param se Numeric vector of standard errors
#' @param maf Numeric vector of minor allele frequencies
#' @param n Numeric vector (or scalar) of sample sizes
#' @return Estimated sdY (scalar)
#' @keywords internal
estimate_sdY <- function(beta, se, maf, n) {
  # V(beta-hat) = sigma^2 / (2 * N * MAF * (1-MAF))
  # => sigma^2 = V(beta-hat) * 2 * N * MAF * (1-MAF)
  # Use regression: varbeta ~ 1/(2*N*MAF*(1-MAF))
  oneover <- 1 / (2 * n * maf * (1 - maf))
  m <- stats::lm(se^2 ~ oneover - 1)
  sqrt(stats::coefficients(m))
}


#' Combine per-SNP ABFs under the 5 colocalization hypotheses
#'
#' @param labf1 Numeric vector of log ABFs for trait 1
#' @param labf2 Numeric vector of log ABFs for trait 2
#' @param p1 Prior probability a SNP is associated with trait 1 only
#' @param p2 Prior probability a SNP is associated with trait 2 only
#' @param p12 Prior probability a SNP is associated with both traits
#' @return Named numeric vector of posterior probabilities (PP.H0-H4)
#' @keywords internal
combine_abf <- function(labf1, labf2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5) {
  # Log-evidence for each hypothesis
  lH0 <- 0
  lH1 <- log(p1) + logsum(labf1)
  lH2 <- log(p2) + logsum(labf2)
  # H3: both associated, different variants = product of marginals minus diagonal
  lH3 <- log(p1) + log(p2) + logdiff(logsum(labf1) + logsum(labf2),
                                       logsum(labf1 + labf2))
  lH4 <- log(p12) + logsum(labf1 + labf2)

  # Posterior probabilities via softmax
  all_H <- c(lH0, lH1, lH2, lH3, lH4)
  pp <- exp(all_H - logsum(all_H))
  names(pp) <- paste0("PP.H", 0:4, ".abf")
  pp
}


#' Approximate Bayes Factor colocalization analysis
#'
#' Tests whether two traits share a causal variant at a genomic locus using
#' approximate Bayes factors computed from GWAS summary statistics.
#'
#' @param dataset1 List with elements: beta, varbeta, snp, type ("quant"/"cc"),
#'   and optionally sdY (quantitative) or s (case-control proportion), or MAF + N.
#' @param dataset2 Same format as dataset1.
#' @param p1 Prior probability a SNP is associated with trait 1 only (default: 1e-4)
#' @param p2 Prior probability a SNP is associated with trait 2 only (default: 1e-4)
#' @param p12 Prior probability a SNP is associated with both traits (default: 1e-5)
#'
#' @return List with elements:
#'   \describe{
#'     \item{summary}{Named vector: nsnps, PP.H0.abf through PP.H4.abf}
#'     \item{results}{Data frame with per-SNP posterior probabilities}
#'     \item{priors}{List with p1, p2, p12}
#'   }
#'
#' @references
#' Giambartolomei et al. (2014) PLoS Genetics 10(5): e1004383.
#' Wallace (2020) PLoS Genetics 16(4): e1008720.
#'
#' @export
coloc_abf <- function(dataset1, dataset2,
                      p1 = 1e-4, p2 = 1e-4, p12 = 1e-5) {

  # --- Resolve prior variance (W) for each dataset ---
  resolve_prior_var <- function(d) {
    if (d$type == "cc") {
      return(0.2^2)
    }
    # Quantitative: W = (0.15 * sdY)^2
    if (!is.null(d$sdY)) {
      return((0.15 * d$sdY)^2)
    }
    # Estimate sdY from MAF + N
    if (!is.null(d$MAF) && !is.null(d$N)) {
      sdY <- estimate_sdY(d$beta, sqrt(d$varbeta), d$MAF, d$N)
      return((0.15 * sdY)^2)
    }
    # Fallback: assume sdY = 1
    return(0.15^2)
  }

  W1 <- resolve_prior_var(dataset1)
  W2 <- resolve_prior_var(dataset2)

  # --- Find shared SNPs ---
  shared <- intersect(dataset1$snp, dataset2$snp)
  if (length(shared) == 0) {
    stop("No shared SNPs between datasets")
  }

  idx1 <- match(shared, dataset1$snp)
  idx2 <- match(shared, dataset2$snp)

  # --- Compute per-SNP log ABFs ---
  labf1 <- compute_labf(dataset1$beta[idx1], sqrt(dataset1$varbeta[idx1]), W1)
  labf2 <- compute_labf(dataset2$beta[idx2], sqrt(dataset2$varbeta[idx2]), W2)

  # --- Combine under 5 hypotheses ---
  pp <- combine_abf(labf1, labf2, p1, p2, p12)

  # --- Per-SNP posterior probabilities ---
  # Under H4, the posterior for each SNP being the shared causal variant
  labf_combined <- labf1 + labf2
  snp_pp_h4 <- exp(labf_combined - logsum(labf_combined))

  # Under H1/H2/H3, compute per-SNP contributions
  snp_pp_h0 <- rep(pp["PP.H0.abf"] / length(shared), length(shared))
  snp_pp_h1 <- pp["PP.H1.abf"] * exp(labf1 - logsum(labf1))
  snp_pp_h2 <- pp["PP.H2.abf"] * exp(labf2 - logsum(labf2))

  # H3 per-SNP: marginal contribution from each trait independently
  # (simplified: distribute proportionally)
  snp_pp_h3 <- rep(pp["PP.H3.abf"] / length(shared), length(shared))

  snp_pp_h4_scaled <- pp["PP.H4.abf"] * snp_pp_h4

  results <- data.frame(
    snp = shared,
    V.df1 = dataset1$varbeta[idx1],
    z.df1 = dataset1$beta[idx1] / sqrt(dataset1$varbeta[idx1]),
    V.df2 = dataset2$varbeta[idx2],
    z.df2 = dataset2$beta[idx2] / sqrt(dataset2$varbeta[idx2]),
    lABF.df1 = labf1,
    lABF.df2 = labf2,
    internal.sum.lABF = labf1 + labf2,
    SNP.PP.H0 = as.numeric(snp_pp_h0),
    SNP.PP.H1 = as.numeric(snp_pp_h1),
    SNP.PP.H2 = as.numeric(snp_pp_h2),
    SNP.PP.H3 = as.numeric(snp_pp_h3),
    SNP.PP.H4 = as.numeric(snp_pp_h4_scaled),
    stringsAsFactors = FALSE
  )

  summary_vec <- c(nsnps = length(shared), pp)

  list(
    summary = summary_vec,
    results = results,
    priors = list(p1 = p1, p2 = p2, p12 = p12)
  )
}
