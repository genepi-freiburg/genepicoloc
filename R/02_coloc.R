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
create_coloc_params_df <- function(EXPERIMENT,
                                   sumstats_1_args,
                                   sumstats_path,
                                   sumstats_pattern = "gz$",
                                   grep_invert = NULL,
                                   sumstats_function,
                                   sumstats_type,
                                   sumstats_sdY = NA,
                                   extra_args = NULL,
                                   do_annotate = F, annotation_function, annotation_function_args
                                   ) {
  # get all files
  files <- list.files(sumstats_path, pattern = sumstats_pattern, full.names = T)
  if (length(files) == 0) {stop("No sumstats found under given path")}
  if (!is.null(grep_invert)) {
    files <- grep(grep_invert, files, value = T, invert = T)
  }
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
      annotation_function_args = annotation_function_args)
  }
  return(list_out)
}

#' Run coloc function using 2 extracted regions
#' @param sumstats_1_df data frame.
#' @param sumstats_1_type quant or cc.
#' @param sumstats_1_df data frame.
#' @param sumstats_2_type quant or cc.
#' @return results of coloc.abf function
#' @examples
#' run_coloc()
#' @export
run_coloc <- function(sumstats_1_df, sumstats_1_type, sumstats_1_sdY = NA,
                      sumstats_2_df, sumstats_2_type, sumstats_2_sdY = NA) {
  sumstats_df_1_coloc <- list(beta=sumstats_1_df$BETA,
                              varbeta=(sumstats_1_df$SE)^2,
                              snp=sumstats_1_df$Name,
                              type=sumstats_1_type)
  if (sumstats_1_type == "quant") {
    if (!is.na(sumstats_1_sdY)) {
      sumstats_df_1_coloc$sdY <- sumstats_1_sdY
    } else
      sumstats_df_1_coloc$MAF <- sumstats_1_df$AF
    sumstats_df_1_coloc$N <- sumstats_1_df$N
  }
  sumstats_df_2_coloc <- list(beta=sumstats_2_df$BETA,
                              varbeta=(sumstats_2_df$SE)^2,
                              snp=sumstats_2_df$Name,
                              type=sumstats_2_type)
  if (sumstats_2_type == "quant") {
    if (!is.na(sumstats_2_sdY)) {
      sumstats_df_2_coloc$sdY <- sumstats_2_sdY
    } else
      sumstats_df_2_coloc$MAF <- sumstats_2_df$AF
    sumstats_df_2_coloc$N <- sumstats_2_df$N
  }
  coloc_res <- coloc.abf(dataset1=sumstats_df_1_coloc, dataset2=sumstats_df_2_coloc)
  return(coloc_res)
}

#' Query sumstats and run coloc, supports multithreading / slurm.
#' @param args_list list of input arguments for coloc
#' usually created using 'create_coloc_params_df' function 
#' @param ellipsis with additional arguments (e.g., annotation files).
#' @return results of coloc.abf function
#' @examples
#' run_coloc()
#' @export
coloc_wrapper <- function(CHR_var, BP_START_var, BP_STOP_var,
                          sumstats_1_file, sumstats_1_function,
                          sumstats_1_type, sumstats_1_sdY,
                          sumstats_2_file, sumstats_2_function,
                          sumstats_2_type, sumstats_2_sdY,
                          ...,
                          do_process_wrapper = T) {
  # Declare nested function
  process_sumstats_2_df <- function(sumstats_1_df, sumstats_2_df) {
    if (nrow(sumstats_1_df) == 0 | nrow(sumstats_2_df) == 0) {
      if (do_process_wrapper == F) {stop("Output is not consistent when nrow=0 and do_process_wrapper=F")}
      if (is.data.frame(sumstats_1_df)) {
        suppressWarnings({sumstats_1_min_P <- min(sumstats_1_df$P, na.rm=T)})
      } else { sumstats_1_min_P <- NULL }
      if (is.data.frame(sumstats_2_df)) {
        suppressWarnings({sumstats_2_min_P <- min(sumstats_2_df$P, na.rm=T)})
      } else { sumstats_2_min_P <- NULL }
      coloc_output <- data.frame(CHR_var = CHR_var, BP_START_var = BP_START_var,
                                 BP_STOP_var = BP_STOP_var,
                                 sumstats_1_file = sumstats_1_file,
                                 sumstats_1_min_P = sumstats_1_min_P,
                                 sumstats_2_file = sumstats_2_file,
                                 sumstats_2_min_P	= sumstats_2_min_P,
                                 nsnps = 0,
                                 PP.H0.abf = NA, PP.H1.abf = NA, PP.H2.abf = NA,
                                 PP.H3.abf = NA, PP.H4.abf = NA, Top_coloc_SNP = NA,
                                 Top_coloc_SNP.PP.H4 = NA, priors = NA)
      return(coloc_output)
    }
    sumstats_1_min_P <- min(sumstats_1_df$P, na.rm=T)
    sumstats_2_min_P <- min(sumstats_2_df$P, na.rm=T)
    coloc_output <- run_coloc(sumstats_1_df = sumstats_1_df,
                              sumstats_1_type = sumstats_1_type,
                              sumstats_1_sdY = sumstats_1_sdY,
                              sumstats_2_df = sumstats_2_df,
                              sumstats_2_type = sumstats_2_type,
                              sumstats_2_sdY = sumstats_2_sdY)
    coloc_output$region <- data.frame(CHR_var = CHR_var,
                                      BP_START_var = BP_START_var,
                                      BP_STOP_var = BP_STOP_var,
                                      sumstats_1_file = sumstats_1_file,
                                      sumstats_1_min_P = sumstats_1_min_P,
                                      sumstats_2_file = sumstats_2_file,
                                      sumstats_2_min_P = sumstats_2_min_P)
    if (do_process_wrapper) {
      coloc_output <- process_wrapper(coloc_output)
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
  sumstats_2_obj <- do.call(sumstats_2_function,
                            c(list(sumstats_file = sumstats_2_file), args_list))
  # sumstats_2_obj can be either a list of data.frames or a data.frame
  # next block will process sumstats_2_obj as a list, so convert first if needed
  if (is.data.frame(sumstats_2_obj)) {sumstats_2_obj <- list(sumstats_2_obj)}
  # process_sumstats_2_df
  sumstats_2_list <- lapply(sumstats_2_obj, function(x) {
    df_out <- process_sumstats_2_df(sumstats_1_df = sumstats_1_df,
                                    sumstats_2_df = x)
    if ("Phenotype" %in% colnames(x)) {
      df_out$sumstats_2_file <- paste0(df_out$sumstats_2_file, "_", unique(x$Phenotype))
    }
    return(df_out)
  })
  # coloc_output <- do.call(rbind, sumstats_2_list)
  return(sumstats_2_list)
}


#' Process results of coloc_wrapper function.
#' @param coloc_output output of coloc_wrapper function.
#' @param N_top_SNPs Number of SNPs with highest PP.H4 to output.
#' @param remove_full_results Should data.frame with full coloc output be removed?
#' Usually TRUE, in this case only first SNPs are used in output.
#' @return results of coloc.abf function
#' @examples
#' run_coloc()
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

#' parallel wrapper
#' @description under development
#' @export
parallel_wrapper <- function(sumstats_2_args,
                             annotation_function = NULL,
                             annotation_function_args = NULL,
                             N_nodes = 10, N_cpus_per_node = 10,
                             do_rbind = T,
                             run_slurm = FALSE, global_objects = NULL,
                             dry_run = T, debug_mode = F,
                             save_RDS = T) {
  EXPERIMENT <- sumstats_2_args$EXPERIMENT
  extra_args <- sumstats_2_args$extra_args
  params_df <- sumstats_2_args$params_df
  if (!is.null(sumstats_2_args$annotate)) {
    do_annotate <- sumstats_2_args$annotate$do_annotate
    annotation_function <- sumstats_2_args$annotate$annotation_function
    annotation_function_args <- sumstats_2_args$annotate$annotation_function_args
  } else {
    do_annotate <- F
  }
  if (do_annotate) {  if (!do_rbind) { stop("do_annotate=T cannot be used with do_rbind=F")}  }
  if (dry_run) {params_df <- params_df[1:2,]; EXPERIMENT <- paste0(EXPERIMENT, "_dryrun")}
  if (debug_mode) {
    coloc_out <- lapply(1:nrow(params_df),
                        function(i) {
                          print(i); do.call(coloc_wrapper, c(params_df[i,], extra_args))
                        })
  }
  if (run_slurm) {
    coloc_slr_job <- slurm_apply(f = coloc_wrapper,
                                 params = params_df,
                                 extra_args,
                                 nodes = N_nodes, cpus_per_node = N_cpus_per_node,
                                 global_objects = global_objects,
                                 jobname = EXPERIMENT, submit = TRUE,
                                 slurm_options = list(time = "12:00:00", share = TRUE))
    coloc_out <- get_slurm_out(coloc_slr_job, outtype = "raw", wait = TRUE, ncores = NULL)
  } else if ("parallel" %in% rownames(installed.packages())) {
    if (debug_mode == F) {
      coloc_out <- parallel::mclapply(1:nrow(params_df),
                                      function(i) do.call(coloc_wrapper, c(params_df[i,], extra_args)),
                                      mc.cores = N_cpus_per_node)
    }
  } else {
    coloc_out <- lapply(1:nrow(params_df),
                        function(i) {
                          print(i); do.call(coloc_wrapper, c(params_df[i,], extra_args))
                        })
  }
  if (do_rbind) {
    coloc_out <- do.call(rbind, lapply(coloc_out, function(x) {do.call(rbind, x)}))
  }
  if (do_annotate) {
    coloc_out <- do.call(annotation_function, c(annotation_function_args, list(coloc_out = coloc_out)))
  }
  if (save_RDS) {
    saveRDS(coloc_out, paste0(EXPERIMENT, ".RDS"))
  }
  return(coloc_out)
  # cleanup_files(coloc_slr_job, wait = TRUE)
}


#' Read sumstats 1 and format columns
#' @param sumstats_file path to sumstats file
#' @param CHR_name name of the column with chromosomes
#' @param BP_name name of the column with positions
#' @param p_value_name name of the column with p-values
#' @param read_method either data.table (preferable, if available) or data.frame
#' @return data frame with formatted columns
#' @examples
#' under development
#' @export
read_sumstats_1 <- function(sumstats_file,
                            CHR_name = "CHR",
                            BP_name = "BP",
                            p_value_name = "P",
                            other_columns,
                            read_method = "data.table") {
  if (read_method == "data.table") {
    sumstats <- data.table::fread(sumstats_file)
    sumstats <- as.data.frame(sumstats)
  }
  if (read_method == "data.frame") {
    sumstats <- read.delim(sumstats_file, header = T)
  }
  sumstats <- sumstats[,c(CHR_name, BP_name, p_value_name, other_columns)]
  colnames(sumstats)[colnames(sumstats) == CHR_name] <- "CHR"
  colnames(sumstats)[colnames(sumstats) == BP_name] <- "BP"
  colnames(sumstats)[colnames(sumstats) == p_value_name] <- "P"
  return(sumstats)
}

#' Get coloc regions
#' @param sumstats data frame read with read_sumstats_1(). Mandatory columns: CHR, BP, P
#' @param p_threshold search for regions until no more variants below this threshold remains
#' @param log_name iteration log will be written to this file
#' @return data frame with extracted regions
#' @description
#' Find significant regions in sumstats and merge close regions if needed
#' from a list of external sumstats.
#' 
#' @examples
#' under development
#' @export
get_coloc_regions <- function(sumstats,
                              CHR_name = "CHR",
                              BP_name = "BP",
                              p_value_name = "P",
                              p_threshold = 5e-8,
                              halfwindow = 500000,
                              log_name = NA) {
  if(!all("CHR" %in% colnames(sumstats) &
          "BP" %in% colnames(sumstats) &
          "P" %in% colnames(sumstats))) {stop("Check column names")}
  # set up variables
  coloc_regions <- data.frame()
  # function-specific constants
  region_var <- 1
  comment_var <- "PASS"
  if (!is.na(log_name)) {
    fileConn <- file(log_name, "w")
    sink(fileConn)
  }
  # start iterations
  while(min(sumstats[[p_value_name]], na.rm = T) < p_threshold) {
    min_p_row <- sumstats[which.min(sumstats[[p_value_name]]),]
    CHR_var <- min_p_row[[CHR_name]]
    BP_var <- min_p_row[[BP_name]]
    BP_START_var <- BP_var - halfwindow
    BP_STOP_var <- BP_var + halfwindow
    print(paste0("Solving region ", region_var))
    print(paste0("Next most significant variant ", CHR_var, ":", BP_var))
    if (nrow(coloc_regions) > 0) {
      coloc_regions_filtered <- subset(coloc_regions, grepl("PASS", comment))
    } else {
      coloc_regions_filtered <- coloc_regions
    }
    if (CHR_var %in% coloc_regions_filtered[[CHR_name]]) {
      closest_region <- subset(coloc_regions_filtered, coloc_regions_filtered[[CHR_name]] == CHR_var)
      closest_region <- closest_region[which.min(sapply(closest_region$BP, function(x) abs(x - BP_var))),]
      print(paste0("Closest region so far: region=", closest_region$region, ", ", closest_region$CHR, ":", closest_region$BP_START, "-", closest_region$BP_STOP))
      if (abs(closest_region$BP_START - BP_var) < halfwindow |
          abs(closest_region$BP_STOP - BP_var) < halfwindow) {
        print(paste0("Next most significant variant is closer than ", halfwindow, " BP to the closest region, merging"))
        if (abs(closest_region$BP_START - BP_var) < halfwindow) {
          coloc_regions[coloc_regions$region == closest_region$region,]$BP_START <- BP_var
        } else if (abs(closest_region$BP_STOP - BP_var) < halfwindow) {
          coloc_regions[coloc_regions$region == closest_region$region,]$BP_STOP <- BP_var
        }
        updated_region <- subset(coloc_regions, region == closest_region$region)
        print(paste0("Updated region: region=", updated_region$region, ", ", updated_region$CHR, ":", updated_region$BP_START, "-", updated_region$BP_STOP))
        comment_var <- paste0("SKIP_merged_to_", updated_region$region)
      } else {
        print(paste0("Next most significant variant is further than ", halfwindow, " BP"))
      }
    } else {
      print(paste0("No hits on this chromosome so far"))
    }
    coloc_regions <- rbind(coloc_regions, data.frame(region = region_var, min_p_row, BP_START = BP_START_var, BP_STOP = BP_STOP_var, comment = comment_var))
    old_indeces <- 1:nrow(sumstats)
    indeces <- which(sumstats[[CHR_name]] == CHR_var & sumstats[[BP_name]] >= BP_START_var & sumstats[[BP_name]] <= BP_STOP_var)
    new_indeces <- old_indeces[! old_indeces %in% indeces]
    stopifnot(nrow(sumstats) - length(indeces) == length(new_indeces))
    sumstats <- sumstats[new_indeces,]
    region_var <- region_var + 1
    comment_var <- "PASS"
    print("----------------")
  }
  # close log
  if (!is.na(log_name)) {
    sink(); close(fileConn)
  }
  # return results
  colnames(coloc_regions)[colnames(coloc_regions) == "CHR"] <- "CHR_var"
  colnames(coloc_regions)[colnames(coloc_regions) == "BP_START"] <- "BP_START_var"
  colnames(coloc_regions)[colnames(coloc_regions) == "BP_STOP"] <- "BP_STOP_var"
  start_cols <- which(colnames(coloc_regions) %in% c("CHR_var", "BP_START_var", "BP_STOP_var"))
  end_cols <- which(!colnames(coloc_regions) %in% c("CHR_var", "BP_START_var", "BP_STOP_var"))
  coloc_regions <- coloc_regions[,c(start_cols, end_cols)]
  rownames(coloc_regions) <- NULL
  return(coloc_regions)
}

#' subset_sumstats_1
#' @description under development
#' @export
subset_sumstats_1 <- function(CHR_var, BP_START_var, BP_STOP_var,
                              sumstats, dbSNP_file = NULL,
                              remove_duplicates = T,
                              do_match_rs = F) {
  sumstats <- get(sumstats)
  if(!all("CHR" %in% colnames(sumstats) &
          "BP" %in% colnames(sumstats) &
          "P" %in% colnames(sumstats))) {stop("Check column names")}
  sumstats_filt <- subset(sumstats, CHR == CHR_var & BP >= BP_START_var & BP <= BP_STOP_var)
  if (do_match_rs) {
    sumstats_filt_rs_matched <- match_rs(dbSNP_file = dbSNP_file,
                                         CHR_var = CHR_var,
                                         BP_START_var = BP_START_var,
                                         BP_STOP_var = BP_STOP_var,
                                         sumstats = sumstats)
    # sumstats_filt_rs_unmatched <- subset(sumstats_filt_rs[!matched_alleles,], !SNP %in% sumstats_filt_rs_matched$SNP)
    sumstats_filt <- sumstats_filt_rs_matched
  }
  if (remove_duplicates) {
    sumstats_filt <- unique(sumstats_filt)
  }
  return(sumstats_filt)
}


#' match_rs
#' @description under development
#' @export
match_rs <- function(dbSNP_file, CHR_var, BP_START_var, BP_STOP_var,
                     sumstats) {
  rs_df <- query_dbSNP(dbSNP_file, CHR_var, BP_START_var, BP_STOP_var)
  stopifnot(all(names(table(sapply(strsplit(rs_df$V3, "rs"), length))) == "2"))
  rs_df <- subset(rs_df, V3 %in% sumstats$SNP)
  rs_df <- rs_df[,c(3,6)]
  rs_df$REF <- gsub("chr[0-9]+:[0-9]+:(.*):.*", "\\1", rs_df$V6)
  rs_df$ALT <- gsub("chr[0-9]+:[0-9]+:.*:(.*)", "\\1", rs_df$V6)
  rs_df <- unique(rs_df)
  sumstats <- merge(sumstats, by.x="SNP",
                    rs_df, by.y="V3")
  matched_alleles <- ((sumstats$REF == sumstats$A1) & (sumstats$ALT == sumstats$A2)) |
    ((sumstats$REF == sumstats$A2) & (sumstats$ALT == sumstats$A1))
  sumstats <- sumstats[matched_alleles,]
  sumstats <- sumstats[,c("V6", "SNP", "CHR", "BP", "A1", "A2", "b", "se", "P", "freq", "N")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  sumstats <- sumstats[order(sumstats$CHR, sumstats$POS),]
  return(sumstats)
}

