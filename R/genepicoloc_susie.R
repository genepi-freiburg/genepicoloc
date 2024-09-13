# genepicoloc_susie ----
#' LD_matrix
#' @export
LD_matrix <- function(sumstats, fixref=F, tokeep=NULL,
                      plink2="plink2", bfile=NULL, rm_tmp=T,
                      Name_var="Name", A1_var="A1") {
  # TODO add functionality to accept precalculated LD matrix
  tmp <- tempfile()
  # bfile for plink2 LD estimation
  if (is.null(bfile)) bfile <- get_bfile(CHR=CHR_var)
  if (fixref) {
    write.table(sumstats[,c(Name_var, A1_var)], tmp, row.names = F, col.names = F, quote=F)
    plink2_cmd <- paste0(plink2, " --bfile ", bfile,
                         " --extract ", tmp,
                         " --ref-allele ", tmp,
                         " --make-bed --out ", paste0(tmp, "_ref"))
    if (!is.null(tokeep)) {
      plink2_cmd <- paste0(plink2_cmd, " --keep ", tokeep, " ")
    }
    system(plink2_cmd)
    bfile <- paste0(tmp, "_ref")
  } else {
    write.table(sumstats[,"Name"], tmp, row.names = F, col.names = F, quote=F)
  }
  plink2_cmd <- paste0(plink2, " --bfile ", bfile, " --extract ", tmp,
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
  # output
  if (rm_tmp) { message("Cleaning tmp files"); unlink(paste0(tmp, "*"))}
  return(r_mat)
}

#' genepicoloc_susie_rss
#' @export
genepicoloc_susie_rss <- function(sumstats, LD_matrix_mat,
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
                           R = LD_matrix_mat,
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
