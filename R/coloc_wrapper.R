#' Query sumstats and run coloc, supports multithreading / slurm.
#' @param args_list list of input arguments for coloc
#' usually created using 'create_coloc_params_df' function 
#' @param ellipsis with additional arguments (e.g., annotation files).
#' @return results of coloc.abf function
#' @export
coloc_wrapper <- function(CHR_var, BP_START_var, BP_STOP_var,
                          sumstats_1_file, sumstats_1_function,
                          sumstats_1_type, sumstats_1_sdY,
                          sumstats_2_file, sumstats_2_function,
                          sumstats_2_type, sumstats_2_sdY,
                          min_nlog10P = -log10(1e-5),
                          silent=T, ...) {
  # Declare nested function to run coloc needed for sumstats_2 which is a list
  process_list_coloc <- function(sumstats_1_df, sumstats_2_df) {
    # Initialize output template
    coloc_template <- out_template(CHR_var, BP_START_var, BP_STOP_var,
                                 sumstats_1_file, sumstats_1_max_nlog10P=sumstats_1_max_nlog10P,
                                 sumstats_2_file, sumstats_2_max_nlog10P=NA,
                                 nsnps=NA)
    ### process sumstats_2
    # checks
    if (!is.data.frame(sumstats_2_df)) {
      coloc_template$comment <- "sumstats_2_not_df"
      return(coloc_template)
    }
    if (nrow(sumstats_2_df) == 0) {
      coloc_template$comment <- "sumstats_2_zero_rows"
      return(coloc_template)
    }
    no_intersect <- all(!(sumstats_1_df$Name %in% sumstats_2_df$Name))
    if (no_intersect) {
      coloc_template$comment <- "no_SNP_intersect"
      return(coloc_template)
    }
    # temporarily disabled - needs more checks
    # sumstats_2_df <- check_sumstats(sumstats_2_df)
    # Add Phenotype column if present
    if ("Phenotype" %in% colnames(sumstats_2_df)) {
      coloc_template$sumstats_2_file <- paste0(coloc_template$sumstats_2_file, "_",
                                             unique(sumstats_2_df$Phenotype))
    }
    # calculate max nlog10P
    if ("nlog10P" %in% colnames(sumstats_2_df)) {
      sumstats_2_max_nlog10P <- max(sumstats_2_df[["nlog10P"]], na.rm=T) 
    } else if ("P" %in% colnames(sumstats_2_df)) {
      sumstats_2_max_nlog10P <- max(-log10(sumstats_2_df[["P"]]), na.rm=T)
    } else { stop("sumstats_2 has neither nlog10P nor P columns") }
    coloc_template$sumstats_2_max_nlog10P <- sumstats_2_max_nlog10P
    # Do not run coloc if there are no significant SNP in one of the sumstats
    if (sumstats_1_max_nlog10P < min_nlog10P) {
      coloc_template$comment <- "no_signif_sumstats1"
      return(coloc_template)
    }
    if (sumstats_2_max_nlog10P < min_nlog10P) {
      coloc_template$comment <- "no_signif_sumstats2"
      return(coloc_template)
    }
    # run coloc if both sumstats have significant SNPs
    coloc_output <- run_coloc(sumstats_1_df=sumstats_1_df,
                              sumstats_1_type=sumstats_1_type,
                              sumstats_1_sdY=sumstats_1_sdY,
                              sumstats_2_df=sumstats_2_df,
                              sumstats_2_type=sumstats_2_type,
                              sumstats_2_sdY=sumstats_2_sdY)
    coloc_template$comment <- "coloc_done"
    coloc_template <- format_coloc_output(coloc_output, coloc_template)
    return(coloc_template)
  }
  
  ### run code
  # process sumstats_1
  sumstats_1_df <- retrieve_sumstats(sumstats_1_function, sumstats_1_file,
                                     CHR_var = CHR_var, BP_START_var = BP_START_var,
                                     BP_STOP_var = BP_STOP_var)
  # checks
  if (!is.data.frame(sumstats_1_df)) stop("sumstats_1 is not a data.frame, please check the input.")
  if (nrow(sumstats_1_df) == 0) stop("sumstats_1 has 0 rows, please check the input.")
  # calculate max nlog10 or min P
  sumstats_1_max_nlog10P <- max(sumstats_1_df[["nlog10P"]], na.rm=T)
  if (sumstats_1_max_nlog10P < min_nlog10P) warning("sumstats_1_max_nlog10P is lower than min_nlog10P")
  # process sumstats_2
  sumstats_2_obj <- retrieve_sumstats(sumstats_2_function, sumstats_2_file,
                                      CHR_var = CHR_var, BP_START_var = BP_START_var,
                                      BP_STOP_var = BP_STOP_var)
  # sumstats_2_obj should be a list of data.frame(s), convert if needed
  if (is.data.frame(sumstats_2_obj)) {sumstats_2_obj <- list(sumstats_2_obj)}
  ### run coloc
  coloc_output <- lapply(sumstats_2_obj, function(sumstats_2_df) {
    process_list_coloc(sumstats_1_df = sumstats_1_df, sumstats_2_df = sumstats_2_df)
  })
  coloc_output <- do.call(rbind, coloc_output)
  return(coloc_output)
}

# helpers ----
#' Create a template data.frame for coloc_wrapper output
out_template <- function(CHR_var, BP_START_var, BP_STOP_var,
                         sumstats_1_file, sumstats_1_max_nlog10P,
                         sumstats_2_file, sumstats_2_max_nlog10P, nsnps=NA) {
  data.frame(CHR_var = CHR_var,
             BP_START_var = BP_START_var,
             BP_STOP_var = BP_STOP_var,
             sumstats_1_file = sumstats_1_file,
             sumstats_1_max_nlog10P = sumstats_1_max_nlog10P,
             sumstats_2_file = sumstats_2_file,
             sumstats_2_max_nlog10P	= sumstats_2_max_nlog10P,
             nsnps = nsnps,
             PP.H0.abf = NA, PP.H1.abf = NA, PP.H2.abf = NA,
             PP.H3.abf = NA, PP.H4.abf = NA,
             Top_coloc_SNP = NA, Top_coloc_SNP.PP.H4 = NA, priors = NA,
             comment=NA)
}

#' Format sumstats data.frame for the coloc.abf input
#' @param sumstats_df data frame.
#' @param sumstats_type quant or cc.
#' @param sumstats_sdY standard deviation for quantitative traits
format_for_coloc <- function(sumstats_df, sumstats_type, sumstats_sdY) {
  sumstats_df_coloc <- list(beta=sumstats_df$BETA,
                            varbeta=(sumstats_df$SE)^2,
                            snp=sumstats_df$Name,
                            type=sumstats_type)
  if (sumstats_type == "quant") {
    if (!is.na(sumstats_sdY)) {
      sumstats_df_coloc$sdY <- sumstats_sdY
    } else {
      sumstats_df_coloc$MAF <- sumstats_df$AF
      sumstats_df_coloc$N <- sumstats_df$N
    }
  }
  return(sumstats_df_coloc)
}

#' Run coloc.abf given either two data.frames:sumstats_1_df, sumstats_2_df,
#' or parameters to query these data.frames: CHR, BP, file, and function
#' sumstats_type and sumstats_sdY are mandatory
#' @importFrom coloc coloc.abf
run_coloc <- function(sumstats_1_df, sumstats_2_df,
                      CHR_var, BP_START_var, BP_STOP_var,
                      sumstats_1_file, sumstats_1_function,
                      sumstats_2_file, sumstats_2_function,
                      sumstats_1_type, sumstats_1_sdY,
                      sumstats_2_type, sumstats_2_sdY,
                      silent=T) {
  tmp_file <- tempfile()
  # Create sumstats data.frames (if sumstats_1_df or sumstats_2_df were not provided)
  if (missing(sumstats_1_df)) {
    sumstats_1_df <- retrieve_sumstats(sumstats_1_function, sumstats_1_file,
                                       CHR_var = CHR_var, BP_START_var = BP_START_var,
                                       BP_STOP_var = BP_STOP_var)
  }
  if (missing(sumstats_2_df)) {
    sumstats_2_df <- retrieve_sumstats(sumstats_2_function, sumstats_2_file,
                                       CHR_var = CHR_var, BP_START_var = BP_START_var,
                                       BP_STOP_var = BP_STOP_var)
  }
  # format both data.frames for coloc.abf input
  sumstats_df_1_coloc <- format_for_coloc(sumstats_1_df, sumstats_1_type, sumstats_1_sdY)
  sumstats_df_2_coloc <- format_for_coloc(sumstats_2_df, sumstats_2_type, sumstats_2_sdY)
  # run coloc.abf
  if (silent) {
    sink(tmp_file, type = "out")
    on.exit(sink())
    coloc_output <- invisible(force(suppressWarnings(coloc.abf(sumstats_df_1_coloc, sumstats_df_2_coloc))))
  } else {
    coloc_output <- coloc.abf(sumstats_df_1_coloc, sumstats_df_2_coloc)
  }
  unlink(tmp_file)
  return(coloc_output)
}


#' Format the output of coloc_wrapper function.
#' @param coloc_output output of coloc::coloc.abf function.
#' @param coloc_template template created by genepicoloc::out_template function.
#' @param N_top_SNPs Number of SNPs with highest PP.H4 to output.
#' @return formatted results of coloc.abf function
format_coloc_output <- function(coloc_output, coloc_template, N_top_SNPs = 5) {
  ### get information from coloc_output
  # SNPs with top PP.H4
  results_top <- coloc_output$results[,c("snp", "SNP.PP.H4")]
  results_top <- head(results_top[order(results_top$SNP.PP.H4, decreasing = T),], N_top_SNPs)
  Top_coloc_SNP <- paste(results_top$snp, collapse = ", ")
  SNP.PP.H4 <- results_top$SNP.PP.H4
  SNP.PP.H4 <- ifelse(SNP.PP.H4 < 0.01, format(SNP.PP.H4, scientific = TRUE, digits = 2),
                      sprintf("%.2f", round(SNP.PP.H4, 2)))
  SNP.PP.H4 <- paste(SNP.PP.H4, collapse = ", ")
  # priors
  priors <- paste(coloc_output$priors, collapse=", ")
  ### annotate coloc_template
  coloc_template$Top_coloc_SNP <- Top_coloc_SNP
  coloc_template$Top_coloc_SNP.PP.H4 <- SNP.PP.H4
  coloc_template$priors <- priors
  coloc_template$nsnps <- coloc_output$summary["nsnps"]
  coloc_template$PP.H0.abf <- coloc_output$summary["PP.H0.abf"]
  coloc_template$PP.H1.abf <- coloc_output$summary["PP.H1.abf"]
  coloc_template$PP.H2.abf <- coloc_output$summary["PP.H2.abf"]
  coloc_template$PP.H3.abf <- coloc_output$summary["PP.H3.abf"]
  coloc_template$PP.H4.abf <- coloc_output$summary["PP.H4.abf"]
  return(coloc_template)
}

