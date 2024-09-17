#' LD_matrix
#' @export
create_LD_matrix <- function(sumstats, bfile, fixref=F, keep=NULL,
                             plink2="plink2", Name_var="Name", A1_var="A1") {
  # TODO add functionality to accept precalculated LD matrix
  tmp <- tempfile()
  plink2_cmd <- paste0(plink2, " --bfile ", bfile)
  if (!is.null(keep)) plink2_cmd <- paste0(plink2_cmd, " --keep ", keep, " ")
  if (fixref) {
    write.table(sumstats[,c(Name_var, A1_var)], tmp, row.names = F, col.names = F, quote=F)
    plink2_cmd2 <- paste0(plink2_cmd, " --extract ", tmp, " --ref-allele ", tmp,
                         " --make-bed --out ", paste0(tmp, "_ref"))
    system(plink2_cmd2)
    bfile <- paste0(tmp, "_ref")
  } else {
    write.table(sumstats[,Name_var], tmp, row.names = F, col.names = F, quote=F)
  }
  plink2_cmd <- paste0(plink2_cmd, " --extract ", tmp,
                       " --r-unphased square ref-based --out ", tmp, "_r_mat")
  system(plink2_cmd)
  # read resulting tables
  r_mat <- read.table(paste0(tmp, "_r_mat.unphased.vcor1"), comment.char = "", header=F)
  r_mat.vars <- read.table(paste0(tmp, "_r_mat.unphased.vcor1.vars"), comment.char = "", header=F)
  r_mat <- as.matrix(r_mat)
  dimnames(r_mat) <- list(r_mat.vars[["V1"]], r_mat.vars[["V1"]])
  # remove NAs
  to_remove_row <- apply(r_mat, 1, function(x) any(is.na(x)))
  message(paste0("Removed ", sum(to_remove_row), " rows"))
  r_mat <- r_mat[!to_remove_row,!to_remove_row]
  # Cleaning
  message("Cleaning tmp files")
  unlink(paste0(tmp, "*"))
  # output
  return(r_mat)
}

#' genepicoloc_susie_rss
#' @export
genepicoloc_susie <- function(raw_gt, pheno, pheno_name, scale_genotypes=F) {
  # input check
  no_IID <- "No IID column found in "
  if (! "IID" %in% colnames(raw_gt)) stop(paste0(no_IID, "the genotype input"))
  if (! "IID" %in% colnames(pheno)) stop(paste0(no_IID, "the phenotype input"))
  if (! pheno_name %in% colnames(pheno)) stop(paste0(pheno_name, " not found among columns in pheno"))
  # data check
  pheno <- pheno[, c("IID", pheno_name)]
  pheno <- subset(pheno, !is.na(pheno[[pheno_name]]))
  raw_gt <- as.data.frame(raw_gt)
  raw_gt <- raw_gt[,!sapply(raw_gt, function(x) {any(is.na(x))})]
  # subset and match
  raw_gt <- raw_gt[raw_gt$IID %in% pheno$IID,]
  pheno <- subset(pheno, IID %in% raw_gt$IID)
  pheno <- pheno[match(raw_gt$IID, pheno$IID),]
  stopifnot(all(pheno$IID == raw_gt$IID))
  # format
  pheno$IID <- NULL
  pheno_mat <- as.matrix(pheno)
  raw_gt_mat <- as.matrix(raw_gt[,7:ncol(raw_gt)])
  if (scale_genotypes) raw_gt_mat <- scale_genotypes(raw_gt_mat)
  # run susie()
  susie_out <- susie(raw_gt_mat, pheno_mat, verbose = TRUE)
  return(susie_out)
}

#' scale_genotypes
scale_genotypes <- function(raw_gt_mat) {
  if (!is.data.frame(raw_gt_mat)) raw_gt_mat <- as.data.frame(raw_gt_mat)
  raw_gt_mat <- sapply(raw_gt_mat, function(x) as.numeric(scale(x)))
  as.matrix(raw_gt_mat)
}


#' genepicoloc_susie_rss
#' @export
genepicoloc_susie_rss <- function(sumstats, LD_matrix,
                                  Name_name="Name", BETA_name="BETA",
                                  SE_name="SE", N_name="N", ...) {
  # checks
  stopifnot(all(sumstats[[Name_name]] == colnames(LD_matrix_mat)))
  stopifnot(all(sumstats[[Name_name]] == rownames(LD_matrix_mat)))
  n <- median(sumstats[[N_name]])
  # run susieR
  susie_rss_params <- list(bhat = sumstats[[BETA_name]],
                           shat = sumstats[[SE_name]],
                           n = n,
                           R = LD_matrix,
                           verbose = TRUE)
  if (length(list(...)) > 0) susie_rss_params <- c(susie_rss_params, list(...))
  fitted_rss <- do.call(susie_rss, susie_rss_params)
  # #   some additional options:
  # refine = T
  # L = 10
  # max_iter = 100000
  # estimate_residual_variance = T # for in-sample LD matrix
  # # QC - under development
  # attr(LD_matrix_mat, "eigen") = eigen(LD_matrix_mat, symmetric = TRUE)
  # lambda <- estimate_s_rss(z, LD_matrix_mat, n=n)
  # condz_in <- kriging_rss(z, LD_matrix_mat, n=n); condz_in$plot
  return(fitted_rss)
}


