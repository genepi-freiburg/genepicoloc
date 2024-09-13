#' @title Query Kidney eQTL
#' @description Query GTEx v8 GWAS data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @examples
#' under development
#' @export
query_kidney_eQTL <- function(sumstats_file,
                              CHR_var, BP_START_var, BP_STOP_var,
                              ...) {
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T), header = T)
  if (nrow(sumstats) == 0) { sumstats <- data.frame(); return(sumstats) }
  sumstats <- unique(sumstats)
  # format by phenotype ID
  sumstats_list <- lapply(unique(sumstats$GeneID), function(x) {
    sumstats <- subset(sumstats, GeneID == x)
    sumstats$AF <- NA
    sumstats$N <- NA
    sumstats <- sumstats[,c("Name", "rsID", "CHR", "POS_hg38", "Alt", "Ref", "Beta", "Std", "Pvalue", "AF", "N", "GeneID")]
    colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N", "Phenotype")
    sumstats <- subset(sumstats, (!is.na(BETA)) & (!is.na(SE)))
    sumstats <- subset(sumstats, (! BETA %in% c(Inf, -Inf)) & (! SE %in% c(Inf, -Inf)))
    if (length(unique(sumstats$Phenotype)) > 1) {stop("Phenotype not unique in output query")}
    return(sumstats)
  })
  return(sumstats_list)
}

#' @title query UKB PPP pGWAS
#' @description Query UKB PPP pGWAS data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_UKB_PPP_EUR <- function(sumstats_file,
                              CHR_var, BP_START_var, BP_STOP_var, ...,
                              handle_underflow=F,
                              colClasses_int=c(4L,5L), ncol_sumstats=14) {
  if (CHR_var == "X") {CHR_var <- "23"}
  colClasses <- readtable_colCl(ncol_sumstats, colClasses_int)
  sumstats <- tabix_fun(sumstats_file, CHR_var, BP_START_var, BP_STOP_var, colClasses)
  sumstats$CHROM[sumstats$CHROM == "23"] <- "X"
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$ID <- paste0("chr", sumstats$CHROM, ":", sumstats$GENPOS, ":", sumstats$ALLELE0, ":", sumstats$ALLELE1)
  sumstats$rsID <- NA
  if (handle_underflow) {
    sumstats$P <- 10^(-handle_underflow(sumstats[["LOG10P"]]))
  } else {
    sumstats$P <- 10^(-sumstats$LOG10P)
  }
  sumstats <- sumstats[,c("ID", "rsID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "BETA", "SE", "P", "A1FREQ", "N")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  return(sumstats)
}

#' @title Query finngen GWAS
#' @description Query finngen GWAS data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @examples
#' query_finngen_GWAS(sumstats_file = "finngen_R9_E4_DM2REN.gz", CHR_var = "1", BP_START_var = 100000, BP_STOP_var = 110000)
#' @export
query_finngen_GWAS <- function(sumstats_file,
                               CHR_var, BP_START_var, BP_STOP_var, ...) {
  if (CHR_var == "X") {CHR_var <- "23"}
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T),  header = T, stringsAsFactors = FALSE, 
                         colClasses = c(alt = "character", ref= "character"))
  
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$X.chrom[sumstats$X.chrom == "23"] <- "X"
  sumstats$Name <- paste0("chr", sumstats$X.chrom, ":", sumstats$pos, ":", sumstats$ref, ":", sumstats$alt)
  sumstats <- sumstats[,c("Name", "rsids", "X.chrom", "pos", "alt", "ref", "beta", "sebeta", "pval", "af_alt")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF")
  return(sumstats)
}


### Helper functions ----

#' tabix_fun
#' Helper function to create tabix_cmd
tabix_fun <- function(sumstats_file, CHR_var, BP_START_var, BP_STOP_var,
                      colClasses=NULL, sep = "\t", header = T,
                      show_cmd=F) {
  tabix_cmd <- paste0("tabix -h ", sumstats_file, " ", CHR_var, ":", BP_START_var, "-", BP_STOP_var)
  if (show_cmd) {message(tabix_cmd)}
  tabix_txt <- system(tabix_cmd, intern = T)
  if (!identical(tabix_txt, character(0))) {
    if (!is.null(colClasses)) {
      sumstats <- read.table(text=tabix_txt, sep = sep, header = header, colClasses = colClasses)
    } else {
      sumstats <- read.table(text=tabix_txt, sep = sep, header = header)
    }
  } else {
    sumstats <- data.frame()
  }
  return(sumstats)
}

#' readtable_colCl
#' Helper function to create tabix_cmd
readtable_colCl <- function(ncol_sumstats, colClasses_int) {
  colClasses <- as.character(rep(NA, ncol_sumstats))
  colClasses[colClasses_int] <- "character"
  colClasses
}


#' Save coloc regions
#' @export
save_coloc_regions <- function(coloc_regions_list, sumstats_name, max_row=100000,
                               SKIP_name="Name", CHR_place=3, POS_place=4,
                               bgzip_bin="bgzip", tabix_bin="tabix",
                               regions_log="regions_log", coloc_regions="coloc_regions",
                               coloc_regions_PASS="coloc_regions_PASS",
                               sumstats_filt="sumstats_filt") {
  writeLines(coloc_regions_list[[regions_log]], con = paste0(sumstats_name, "_log.txt"))
  write.table(coloc_regions_list[[coloc_regions]], paste0(sumstats_name, "_", coloc_regions, ".tsv"),
              sep="\t", row.names = F, col.names = T, quote = F)
  write.table(coloc_regions_list[[coloc_regions_PASS]], paste0(sumstats_name, "_", coloc_regions_PASS, ".tsv"),
              sep="\t", row.names = F, col.names = T, quote = F)
  # if (nrow(coloc_regions_list[[sumstats_filt]]) > max_row) {
  write.table(coloc_regions_list[[sumstats_filt]], paste0(sumstats_name, "_subset.tsv"),
              sep="\t", row.names = F, col.names = T, quote = F)
  system(paste0(bgzip_bin, " -f ", sumstats_name, "_subset.tsv"))
  system(paste0("tabix -f -s", CHR_place, " -b", POS_place, " -e", POS_place, " ", sumstats_name, "_subset.tsv.gz -c ", SKIP_name))
  # } else {
  #   saveRDS(coloc_regions_list[[sumstats_filt]], paste0(sumstats_name, "_subset.RDS"))
  # }
}

# new line
#' Create data.frame with input parameters for coloc
#' @param sumstats_path path to folder with indexed sumstats.
#' @param sumstats_pattern usually "gz$".
#' @param grep_invert expression to exclude some found files
#' (e.g., annotation files in the same folder).
#' @param sumstats_function function name to query sumstats.
#' @param sumstats_type quant or cc.
#' @param sumstats_sdY sdY for quant traits,
#' by default NA and estimated from BETA and MAF.
#' @return data.frame with input parameters for coloc
#' @examples
#' create_coloc_params_df()
#' @export
create_coloc_params_df <- function(sumstats_1_args,
                                   EXPERIMENT,
                                   sumstats_path, files = NULL,
                                   sumstats_pattern = "gz$",
                                   grep_invert = NULL,
                                   sumstats_function,
                                   hyprcoloc = F,
                                   sumstats_type,
                                   sumstats_sdY = NA,
                                   extra_args = NULL,
                                   do_annotate = F, annotation_function, annotation_function_args,
                                   do_annotate_sumstats_1 = F, annotation_function_sumstats_1, annotation_function_args_sumstats_1
) {
  # get all files
  if (!is.null(sumstats_path)) {
    files <- list.files(sumstats_path, pattern = sumstats_pattern, full.names = T)
    if (length(files) == 0) {stop("No sumstats found under given path")}
    if (!is.null(grep_invert)) {
      files <- grep(grep_invert, files, value = T, invert = T)
    }
  } else {
    if (is.null(files)) { stop("files cannot be NULL when sumstats_path is NULL") }
    files <- files
  }
  if (hyprcoloc) {
    df_out <- data.frame(sumstats_file = files,
                         sumstats_function = sumstats_function)
    if (!is.null(extra_args)) {
      df_out$extra_args <- extra_args
    } else {
      df_out$extra_args <- NA
    }
    return(df_out) }
  sumstats_2_args <- data.frame(sumstats_2_file = files,
                                sumstats_2_function = sumstats_function,
                                sumstats_2_type = sumstats_type,
                                sumstats_2_sdY = sumstats_sdY)
  params_df <- merge(sumstats_1_args, sumstats_2_args)
  # This block is used for GTEXv8 study where sumstats are split by chr
  # Potentially relevant for other studies
  if (EXPERIMENT == "GTEXv8") {
    GTEXv8_CHR_match <- params_df$CHR_var ==
      gsub(".*chr([0-9X]+).*", "\\1", params_df$sumstats_2_file)
    params_df <- params_df[GTEXv8_CHR_match,]
  }
  list_out <- list(
    EXPERIMENT = EXPERIMENT,
    params_df = params_df)
  if (!is.null(extra_args)) {
    list_out$extra_args <- extra_args
  }
  if (do_annotate) {
    list_out$annotate <- list(
      do_annotate = do_annotate,
      annotation_function = annotation_function,
      annotation_function_args = c(list(study = EXPERIMENT),
                                   annotation_function_args))
  }
  if (do_annotate_sumstats_1) {
    list_out$annotate_sumstats_1 <- list(
      do_annotate_sumstats_1 = do_annotate_sumstats_1,
      annotation_function_sumstats_1 = annotation_function_sumstats_1,
      annotation_function_args_sumstats_1 = annotation_function_args_sumstats_1)
  }
  return(list_out)
}

#' Format input df for the coloc input
#' @param sumstats_df data frame.
#' @param sumstats_type quant or cc.
#' @param sumstats_sdY standard deviation for quantitative traits
#' @export
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


#' Template for coloc output
#' under development
#' @export
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
             Top_coloc_SNP = NA, Top_coloc_SNP.PP.H4 = NA, priors = NA)
}


#' Query sumstats and run coloc, supports multithreading / slurm.
#' @param args_list list of input arguments for coloc
#' usually created using 'create_coloc_params_df' function 
#' @param ellipsis with additional arguments (e.g., annotation files).
#' @return results of coloc.abf function
#' @export
coloc_wrapper <- function(CHR_var, BP_START_var, BP_STOP_var,
                          sumstats_1_file, sumstats_1_function,
                          ..., sumstats_1_type, sumstats_1_sdY,
                          sumstats_2_file, sumstats_2_function,
                          sumstats_2_type, sumstats_2_sdY,
                          hyprcoloc_mode = F, check_sumstats_2 = F,
                          do_process_wrapper = T, min_nlog10P = -log10(1e-5)) {
  # Declare nested function
  process_sumstats_2_df <- function(sumstats_1_df, sumstats_2_df) {
    # if sumstats_1/2 queries return 0 rows or there is no SNP intersect
    no_intersect <- all(!(sumstats_1_df$Name %in% sumstats_2_df$Name))
    if (nrow(sumstats_1_df) == 0 | nrow(sumstats_2_df) == 0 | no_intersect) {
      if (do_process_wrapper == F) {stop("Output is not consistent when nrow=0 and do_process_wrapper=F")}
      coloc_output <- out_template(CHR_var, BP_START_var, BP_STOP_var,
                                   sumstats_1_file, sumstats_1_max_nlog10P=NA,
                                   sumstats_2_file, sumstats_2_max_nlog10P=NA,
                                   nsnps=0)
    } else {
      # calculate max nlog10 or min P
      if (is.data.frame(sumstats_1_df)) {
        sumstats_1_max_nlog10P <- max(sumstats_1_df[["nlog10P"]], na.rm=T)
      } else { sumstats_1_max_nlog10P <- NULL }
      if (is.data.frame(sumstats_2_df)) {
        sumstats_2_max_nlog10P <- max(-log10(sumstats_2_df[["P"]]), na.rm=T)
      } else { sumstats_2_max_nlog10P <- NULL }
      # Do not run coloc if there are no significant SNP
      if (sumstats_1_max_nlog10P < min_nlog10P | sumstats_2_max_nlog10P < min_nlog10P) {
        coloc_output <- out_template(CHR_var, BP_START_var, BP_STOP_var,
                                     sumstats_1_file, sumstats_1_max_nlog10P,
                                     sumstats_2_file, sumstats_2_max_nlog10P)
      } else {
        # run coloc if both sumstats have significant SNPs
        sumstats_df_1_coloc <- format_for_coloc(sumstats_1_df, sumstats_1_type, sumstats_1_sdY)
        sumstats_df_2_coloc <- format_for_coloc(sumstats_2_df, sumstats_2_type, sumstats_2_sdY)
        coloc_output <- coloc.abf(dataset1=sumstats_df_1_coloc, dataset2=sumstats_df_2_coloc)
        coloc_output$region <- data.frame(CHR_var = CHR_var,
                                          BP_START_var = BP_START_var,
                                          BP_STOP_var = BP_STOP_var,
                                          sumstats_1_file = sumstats_1_file,
                                          sumstats_1_max_nlog10P = sumstats_1_max_nlog10P,
                                          sumstats_2_file = sumstats_2_file,
                                          sumstats_2_max_nlog10P = sumstats_2_max_nlog10P)
        
        if (do_process_wrapper) {
          coloc_output <- process_wrapper(coloc_output)
        }
      }
    }
    return(coloc_output)
  }
  # run code
  extra_args <- list(...)
  args_list <- list(CHR_var = CHR_var, BP_START_var = BP_START_var,
                    BP_STOP_var = BP_STOP_var)
  if (length(extra_args) > 0) {
    args_list <- c(args_list, extra_args)
  }
  sumstats_1_df <- do.call(sumstats_1_function,
                           c(list(sumstats_file = sumstats_1_file), args_list))
  if (hyprcoloc_mode) {
    return(sumstats_1_df)
  }
  # handle character P in case of underflow - already done at preprocessing
  # sumstats_1_df[["P"]] <- as.numeric(sumstats_1_df[["P"]])
  sumstats_2_obj <- do.call(sumstats_2_function,
                            c(list(sumstats_file = sumstats_2_file), args_list))
  # sumstats_2_obj can be either a list of data.frames or a data.frame
  # next block will process sumstats_2_obj as a list, so convert first if needed
  if (is.data.frame(sumstats_2_obj)) {sumstats_2_obj <- list(sumstats_2_obj)}
  # process_sumstats_2_df
  coloc_output <- lapply(sumstats_2_obj, function(sumstats_2_df) {
    # handle character P in case of underflow - already done at preprocessing
    # sumstats_2_df[["P"]] <- as.numeric(sumstats_2_df[["P"]])
    # sumstats_1_df <- check_sumstats(sumstats_1_df)
    if (check_sumstats_2) {
      sumstats_2_df <- check_sumstats(sumstats_2_df)
    }
    df_out <- process_sumstats_2_df(sumstats_1_df = sumstats_1_df,
                                    sumstats_2_df = sumstats_2_df)
    if ("Phenotype" %in% colnames(sumstats_2_df)) {
      df_out[["sumstats_2_file"]] <- paste0(df_out[["sumstats_2_file"]], "_", unique(sumstats_2_df[["Phenotype"]]))
    }
    return(df_out)
  })
  return(coloc_output)
}

#' Process results of coloc_wrapper function.
#' @param coloc_output output of coloc_wrapper function.
#' @param N_top_SNPs Number of SNPs with highest PP.H4 to output.
#' @param remove_full_results Should data.frame with full coloc output be removed?
#' Usually TRUE, in this case only first SNPs are used in output.
#' @return results of coloc.abf function
#' @export
process_wrapper <- function(coloc_output,
                            N_top_SNPs = 5,
                            remove_full_results = T) {
  results_top <- coloc_output$results[,c("snp", "SNP.PP.H4")]
  results_top <- head(results_top[order(results_top$SNP.PP.H4, decreasing = T),], N_top_SNPs)
  Top_coloc_SNP <- paste(results_top$snp, collapse = ", ")
  SNP.PP.H4 <- results_top$SNP.PP.H4
  SNP.PP.H4 <- ifelse(SNP.PP.H4 < 0.01, format(SNP.PP.H4, scientific = TRUE, digits = 2),
                      sprintf("%.2f", round(SNP.PP.H4, 2)))
  SNP.PP.H4 <- paste(SNP.PP.H4, collapse = ", ")
  priors <- paste(coloc_output$priors, collapse=", ")
  annotation_df <- data.frame(Top_coloc_SNP = Top_coloc_SNP,
                              Top_coloc_SNP.PP.H4 = SNP.PP.H4,
                              priors = priors)
  coloc_output$summary_df <- data.frame(coloc_output$region,
                                        data.frame(t(coloc_output$summary)),
                                        annotation_df)
  if (remove_full_results) {
    coloc_output$results <- NULL
  }
  return(coloc_output$summary_df)
}

# annotation_function = NULL,
# annotation_function_args = NULL,
# global_objects = NULL,
# N_nodes = NULL,
#' parallel wrapper
#' @description under development
#' @export
parallel_wrapper <- function(args_df, N_cpus_per_node = 10, output_folder="output",
                             do_rbind = T, do_annotate = NULL, do_annotate_sumstats_1 = NULL,
                             save_RDS = T, save_RDS_no_annotation = F,
                             dry_run = T, debug_mode = F,
                             min_nlog10P = -log10(1e-5),
                             run_slurm = NULL, check_sumstats_2 = F) {
  # Setup
  if (!is.null(run_slurm)) {
    stop("Slurm functionality is (temporarily) disabled.",
         " Please remove 'run_slurm' parameter and the function will execute the standard run.",
         " Wrap in slurm separately, if necessary (e.g., for each sumstats_2).")
  }
  if (!"parallel" %in% rownames(installed.packages())) {
    stop("'parallel' is currently required to run this function'")
  }
  system(paste0("mkdir -p ", output_folder))
  # Arguments
  EXPERIMENT <- args_df$EXPERIMENT
  print(EXPERIMENT)
  extra_args <- args_df$extra_args
  if (!is.null(extra_args)) {
    extra_args <- list(extra_args, min_nlog10P = min_nlog10P)
  } else {
    extra_args <- list(min_nlog10P = min_nlog10P)
  }
  params_df <- args_df$params_df
  if (is.null(do_annotate)) {
    if (!is.null(args_df$annotate)) {
      do_annotate <- args_df$annotate$do_annotate
      annotation_function <- args_df$annotate$annotation_function
      annotation_function_args <- args_df$annotate$annotation_function_args
    } else {
      do_annotate <- F
    }
  }
  if (is.null(do_annotate_sumstats_1)) {
    if (!is.null(args_df$annotate_sumstats_1)) {
      do_annotate_sumstats_1 <- args_df$annotate_sumstats_1$do_annotate_sumstats_1
      annotation_function_sumstats_1 <- args_df$annotate_sumstats_1$annotation_function_sumstats_1
      annotation_function_args_sumstats_1 <- args_df$annotate_sumstats_1$annotation_function_args_sumstats_1
    } else {
      do_annotate_sumstats_1 <- F
    }
  }
  if (do_annotate | do_annotate_sumstats_1) {  if (!do_rbind) { stop("do_annotate=T or do_annotate_sumstats_1=T cannot be used with do_rbind=F")}  }
  if (dry_run & !(debug_mode)) {params_df <- params_df[1:2,]; EXPERIMENT <- paste0(EXPERIMENT, "_dryrun")}
  if (debug_mode) {
    coloc_out <- lapply(1:nrow(params_df),
                        function(i) {
                          print(i); do.call(coloc_wrapper, c(params_df[i,], extra_args, check_sumstats_2=check_sumstats_2))
                        })
  } else {
    coloc_out <- parallel::mclapply(1:nrow(params_df),
                                    function(i) {
                                      do.call(coloc_wrapper, c(params_df[i,], extra_args, check_sumstats_2=check_sumstats_2))
                                    }, mc.cores = N_cpus_per_node)
  }
  if (do_rbind) {
    coloc_out <- do.call(rbind, lapply(coloc_out, function(x) {do.call(rbind, x)}))
  }
  if (do_annotate) {
    coloc_out <- do.call(annotation_function, c(annotation_function_args, list(coloc_out = coloc_out)))
  }
  if (do_annotate_sumstats_1) {
    coloc_out <- do.call(annotation_function_sumstats_1, c(annotation_function_args_sumstats_1, list(coloc_out = coloc_out)))
  }
  if (save_RDS_no_annotation | save_RDS) {
    EXPERIMENT <- paste0(output_folder, "/", EXPERIMENT)
  }
  if (save_RDS_no_annotation) {
    saveRDS(coloc_out, paste0(EXPERIMENT, "_no_annotation.RDS"))
  }
  if (save_RDS) {
    saveRDS(coloc_out, paste0(EXPERIMENT, ".RDS"))
  }
  return(coloc_out)
}

#' Summarize coloc results
#' @description under development
#' @export
summarize_coloc <- function(selected_studies,
                            output_folder = "output",
                            remove_dirname = T,
                            PP.H4.abf_filt=0.5,
                            PP.H3.abf_filt=NULL,
                            do_summary=T, do_xlsx=T) {
  if (! "data.table" %in% rownames(installed.packages())) { do_summary <- F }
  if (! "writexl" %in% rownames(installed.packages())) { do_xlsx <- F }
  files <- list.files(path = output_folder, pattern = ".RDS", full.names = T)
  files <- grep("dryrun|no_annotation|summary", files, invert = T, value = T)
  coloc_out <- sapply(files, readRDS, simplify = F)
  if (length(coloc_out) > 0) {
    names(coloc_out) <- gsub(".RDS", "", basename(names(coloc_out)))
    # filter before save
    coloc_out_filt <- sapply(coloc_out, function(x){
      x <- subset(x, !is.na(PP.H4.abf))
      if (remove_dirname) {
        x[["sumstats_1_file"]] <- basename(x[["sumstats_1_file"]])
        x[["sumstats_2_file"]] <- basename(x[["sumstats_2_file"]])
      }
      x
    }, simplify = F)
    # create summary table
    coloc_out_combined <- sapply(coloc_out_filt, function(x){
      if (!is.null(PP.H3.abf_filt)) {
        x <- subset(x, PP.H3.abf >= PP.H3.abf_filt | PP.H4.abf >= PP.H4.abf_filt)
      } else {
        x <- subset(x, PP.H4.abf >= PP.H4.abf_filt)
      }
      if (remove_dirname) {
        x[["sumstats_1_file"]] <- basename(x[["sumstats_1_file"]])
        x[["sumstats_2_file"]] <- basename(x[["sumstats_2_file"]])
      }
      x
    }, simplify = F)
    # remove empty results
    for (i in names(coloc_out_combined)) {
      if (nrow(coloc_out_combined[[i]]) == 0)
        coloc_out_combined[[i]] <- NULL
    }
    if (do_summary) {
      coloc_out_filt[["summary"]] <- data.table::rbindlist(coloc_out_combined, fill=TRUE, idcol = "Dataset")
      saveRDS(coloc_out_filt[["summary"]], paste0(output_folder, "/summary.RDS"))
    }
    sapply(names(coloc_out_filt), function(x) {
      if (do_xlsx) {
        writexl::write_xlsx(coloc_out_filt[[x]], paste0(output_folder, "/", x, ".xlsx"))
      } else {
        write.table(coloc_out_filt[[x]], paste0(output_folder, "/", x, ".csv"), sep = "\t", row.names = F, quote = F)
      }
    })
  }
}


list_to_create_args_list <- list(
  Kidney_eQTL =
    list(EXPERIMENT = "Kidney_eQTL",
         sumstats_path = NULL,
         files = c("genepicoloc/data/Kidney_eQTL.txt.gz"),
         sumstats_pattern = "gz$",
         sumstats_function = "query_kidney_eQTL",
         sumstats_type = "quant",
         sumstats_sdY = 1,
         do_annotate = T,
         annotation_function = "transcriptomics_annotation",
         annotation_function_args = data.frame(annotation_file = "genepicoloc/data/gencode.v26.GRCh38.genes.gtf_genes_format.txt")
    ),
  UKB_PPP_EUR =
    list(EXPERIMENT = "UKB_PPP_EUR",
         sumstats_path = NULL,
         files = c("genepicoloc/data/UMOD_P07911_OID20237_v1_Cardiometabolic.txt.gz"),
         sumstats_pattern = "gz$",
         sumstats_function = "query_UKB_PPP_EUR",
         sumstats_type = "quant",
         sumstats_sdY = NA,
         do_annotate = T,
         annotation_function = "proteomics_annotation",
         annotation_function_args = data.frame(annotation_file = "genepicoloc/data/olink_protein_map_3k_v1.tsv")
         ),
  FinnGen_r9 =
    list(EXPERIMENT = "FinnGen_r9",
         sumstats_path = NULL,
         files = c("genepicoloc/data/finngen_R9_N14_CHRONKIDNEYDIS.gz"),
         sumstats_pattern = "gz$",
         sumstats_function = "query_finngen_GWAS",
         sumstats_type = "cc",
         sumstats_sdY = NA,
         do_annotate = T,
         annotation_function = "standard_annotation",
         annotation_function_args = list(annotation_file = "genepicoloc/data/endpoints.tsv")
    )
)

#' transcriptomics annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
transcriptomics_annotation <- function(study, annotation_file,
                                       sumstats_file = "sumstats_2_file",
                                       coloc_out,
                                       CHR_var = "CHR_var", BP_START_var = "BP_START_var",
                                       BP_STOP_var = "BP_STOP_var") {
  annotation_df <- read.delim(annotation_file)
  colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
  merge_column <- paste0(study, "_gene_id_no_dot")
  annotation_df[[merge_column]] <- gsub("(ENSG[0-9]+).?.*", "\\1", annotation_df[[paste0(study, "_gene_id")]])
  coloc_out[[merge_column]] <- gsub(".*(ENSG[0-9]+).?.*", "\\1", basename(coloc_out[[sumstats_file]]))
  coloc_out[[paste0(study, "_Tissue")]] <- gsub("_ENSG.*", "", gsub("Formated", "" , gsub("_hg38.txt.gz", "", basename(coloc_out[[sumstats_file]]))))
  coloc_out <- merge(coloc_out, annotation_df, by = merge_column, all.x=T,
                     sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  coloc_out[[paste0(study, "_cis_trans")]] <- cis_trans_annotation(region_CHR_vec = coloc_out[[CHR_var]],
                                                                   region_BP_START_vec = coloc_out[[BP_START_var]],
                                                                   region_BP_STOP_vec = coloc_out[[BP_STOP_var]],
                                                                   gene_chr_vec = coloc_out[[paste0(study, "_chr")]],
                                                                   gene_start_vec = coloc_out[[paste0(study, "_gene_start")]],
                                                                   suggestive_window = 1e6)
  return(coloc_out)
}

#' proteomics annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
proteomics_annotation <- function(study, annotation_file,
                                  sumstats_file = "sumstats_2_file",
                                  coloc_out, CHR_var = "CHR_var", BP_START_var = "BP_START_var",
                                  BP_STOP_var = "BP_STOP_var") {
  merge_column <- paste0(study, "_OlinkID")
  annotation_df <- read.delim(annotation_file, sep="\t")
  colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
  annotation_df[[paste0(study, "_multiple_genes_per_OID")]] <- 0
  selected_cols <- paste0(study, "_", c("OlinkID", "olink_target_fullname", "UniProt", "Assay",
                                        "HGNC.symbol", "ensembl_id", "chr", "gene_start", "gene_end",
                                        "multiple_genes_per_OID"))
  annotation_df <- annotation_df[,selected_cols]
  if (sumstats_file == "sumstats_2_file") {
    coloc_out[[merge_column]] <- gsub(".*(OID[0-9]+).*", "\\1", coloc_out[[sumstats_file]])
  }
  cis_trans_column <- paste0(study, "_cis_trans")
  gene_chr_vec_name <- paste0(study, "_chr")
  gene_start_vec_name <- paste0(study, "_gene_start")
  coloc_out <- merge(coloc_out, annotation_df, by=merge_column, all.x=T,
                     sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  coloc_out[[cis_trans_column]] <- cis_trans_annotation(region_CHR_vec = coloc_out[[CHR_var]],
                                                        region_BP_START_vec = coloc_out[[BP_START_var]],
                                                        region_BP_STOP_vec = coloc_out[[BP_STOP_var]],
                                                        gene_chr_vec = coloc_out[[gene_chr_vec_name]],
                                                        gene_start_vec = coloc_out[[gene_start_vec_name]],
                                                        suggestive_window = 1e6)
  return(coloc_out)
}


#' standard annotation
#' mGWAS, UKB TOPMed, FinnGen r9
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
standard_annotation <- function(study, annotation_file,
                                sumstats_file = "sumstats_2_file",
                                coloc_out) {
  merge_column <- paste0(study, "_phenocode")
  annotation_df <- read.delim(annotation_file, colClasses = "character")
  coloc_out[[merge_column]] <- gsub("finngen_R9_(.*).gz", "\\1", basename(coloc_out[[sumstats_file]]))
  colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
  nrow_before <- nrow(coloc_out)
  coloc_out <- merge(coloc_out, by = merge_column, all.x=T,
                     annotation_df, sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
 if(nrow_before != nrow(coloc_out)) { stop("Merge produced different number of rows, check duplicates or missing annotations") }
  return(coloc_out)
}

cis_trans_annotation <- function(region_CHR_vec, region_BP_START_vec, region_BP_STOP_vec,
                                 gene_chr_vec, gene_start_vec,
                                 suggestive_window = 1e6) {
  cis_condition <- (gene_chr_vec == region_CHR_vec) &
    (gene_start_vec >= region_BP_START_vec & gene_start_vec <= region_BP_STOP_vec)
  suggestive_cis_condition <- (gene_chr_vec == region_CHR_vec) &
    (gene_start_vec >= region_BP_START_vec-suggestive_window &
       gene_start_vec <= region_BP_STOP_vec+suggestive_window)
  trans_condition <- ((gene_chr_vec != region_CHR_vec) | 
                        (gene_start_vec < region_BP_START_vec-suggestive_window) |
                        (gene_start_vec > region_BP_STOP_vec+suggestive_window))
  cis_trans <- ifelse(cis_condition, "cis",
                      ifelse(suggestive_cis_condition, "suggestive_cis",
                             ifelse(trans_condition, "trans", NA)))
  return(cis_trans)
}

