#' Create data.frame with input parameters for coloc
#' @param sumstats_number 1 or 2 (pair of sumstats for coloc).
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
create_coloc_params_df <- function(sumstats_number = 2,
                                   sumstats_path,
                                   sumstats_pattern = "gz$",
                                   grep_invert = NULL,
                                   sumstats_function,
                                   sumstats_type,
                                   sumstats_sdY = NA) {
  files <- list.files(sumstats_path, pattern = sumstats_pattern, full.names = T)
  if (length(files) == 0) {stop("No sumstats found under given path")}
  if (!is.null(grep_invert)) {
    files <- grep(grep_invert, files, value = T, invert = T)
  }
  params_df <- data.frame(matrix(nrow = length(files), ncol = 0))
  params_df[[paste0("sumstats_", sumstats_number, "_file")]] <- files
  params_df[[paste0("sumstats_", sumstats_number, "_function")]] <- sumstats_function
  params_df[[paste0("sumstats_", sumstats_number, "_type")]] <- sumstats_type
  params_df[[paste0("sumstats_", sumstats_number, "_sdY")]] <- sumstats_sdY
  return(params_df)
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
coloc_wrapper <- function(args_list,
                          ...,
                          do_process_wrapper = T) {
  extra_args <- list(...)
  for (i in names(args_list)) {
    assign(i, args_list[[i]])
  }
  region_list <- list(CHR_var, BP_START_var, BP_STOP_var)
  sumstats_1_df <- do.call(sumstats_1_function, c(list(sumstats_1_file), region_list, extra_args)) # extra_args
  sumstats_1_min_P <- min(sumstats_1_df$P, na.rm=T)
  sumstats_2_df <- do.call(sumstats_2_function, c(list(sumstats_2_file), region_list, extra_args))
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

coloc_wrapper_02 <- function(params_df_1,
                             params_df_2,
                             EXPERIMENT,
                             chunk_size = 10000,
                             do_process_wrapper = T,
                             extra_arguments) {
  params_df <- merge(params_df_1, params_df_2)
  system(paste0("mkdir -p resultsRDS"))
  params_chunks <- params_df_to_chunks(params_df, chunk_size = chunk_size)
  coloc_out_list <- lapply(1:length(params_chunks), function(i) {
    params_list <- params_df_to_list(params_chunks[[i]])
    coloc_slr_job <- slurm_map(x = params_list,
                               coloc_wrapper,
                               do_process_wrapper = do_process_wrapper,
                               extra_arguments,
                               nodes = 40, cpus_per_node = 5,
                               jobname = paste0("coloc_", i), submit = TRUE,
                               slurm_options = list(time = "12:00:00", share = TRUE))
    coloc_out <- get_slurm_out(coloc_slr_job, outtype = "raw", wait = TRUE, ncores = NULL)
    # cleanup_files(coloc_slr_job, wait = TRUE)
    # coloc_out <- mclapply(params_list, coloc_wrapper, mc.cores = 10)
    # coloc_out <- mclapply(coloc_out, process_wrapper, mc.cores = 10)
    saveRDS(coloc_out, paste0("resultsRDS/batch_", i, "_coloc.RDS"))
    # return(coloc_out)
  })
  # return(coloc_out_list)
  # coloc_out_df <- do.call(rbind, )
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

#' Find significant regions in sumstats and merge close regions if needed
#' from a list of external sumstats.
#'
#' @param sumstats data frame read with read_sumstats_1(). Mandatory columns: CHR, BP, P
#' @param p_threshold search for regions until no more variants below this threshold remains
#' @param log_name iteration log will be written to this file
#' @return data frame with extracted regions
#' @examples
#' under development
#' sumstats=read.table('CAD.colocalization.tsv.gz', sep='\t',header=T)
#' get_coloc_regions(sumstats, p_value_name = "P.value", p_threshold = 5e-100)
#' @export
get_coloc_regions <- function(sumstats,
                              CHR_name = "CHR",
                              BP_name = "BP",
                              p_value_name = "P",
                              p_threshold = 5e-8,
                              halfwindow = 500000,
                              log_name = "log.txt",
                              regions_name = "regions.txt") {
  if(!all("CHR" %in% colnames(sumstats) &
              "BP" %in% colnames(sumstats) &
              "P" %in% colnames(sumstats))) {stop("Check column names")}
  # set up variables
  coloc_regions <- data.frame()
  # function-specific constants
  region_var <- 1
  comment_var <- "PASS"
  fileConn <- file(log_name, "w")
  sink(fileConn)
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
  # close log and return results
  sink(); close(fileConn)
  coloc_regions_short <- coloc_regions[,c("CHR", "BP_START", "BP_STOP")]
  colnames(coloc_regions_short) <- paste0(colnames(coloc_regions_short), "_var")
  write.table(coloc_regions_short, regions_name, sep="\t", quote=F, row.names = F)
  write.table(coloc_regions, gsub(".txt", "_full.txt", regions_name), sep="\t", quote=F, row.names = F)
  return(coloc_regions_short)
}

#' Read sumstats 1 and format columns
#' @param sumstats sumstats file created by read_sumstats_1()
#' @param coloc_regions file created by get_coloc_regions()
#' @return data frame with formatted columns
#' @examples
#' under development
#' @export
subset_sumstats_1 <- function(sumstats,
                              coloc_regions) {
  if(!all("CHR" %in% colnames(sumstats) &
          "BP" %in% colnames(sumstats) &
          "P" %in% colnames(sumstats))) {stop("Check column names")}
  sumstats_filt_list <- lapply(1:nrow(coloc_regions), function(i) {
    coloc_regions_format_i <- coloc_regions[i,]
    subset(sumstats, CHR == coloc_regions_format_i$CHR_var &
             BP >= coloc_regions_format_i$BP_START_var &
             BP <= coloc_regions_format_i$BP_STOP_var)
  })
  sumstats_filt <- do.call(rbind, sumstats_filt_list)
  # write.table(sumstats_filt, "sumstats_filt.tsv", quote=F, row.names = F, sep="\t")
  return(sumstats_filt)
}
