#' Apply "coloc" to extracted regions from 2 summary statistics
#' @param sumstats_1_df data frame.
#' @param sumstats_1_type quant or cc.
#' @param sumstats_1_df data frame.
#' @param sumstats_2_type quant or cc.
#' @return results of coloc
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

wrapper_run_coloc <- function(...,
                              mclapply_use = F,
                              remove_full_results = T,
                              N_top_SNPs = 5,
                              add_annotation = T) {
  extra_args <- list(...)
  if (mclapply_use) {
    extra_args <- unlist(extra_args, recursive = FALSE)
  }
  if (length(extra_args) > 0) {
    for (i in 1:length(extra_args)) {
      if (names(extra_args)[i] == "") {
        assign(names(extra_args[[i]]), extra_args[[i]])
      } else {
        assign(names(extra_args)[i], extra_args[[i]])
      }
    }
  } else {
    stop("ellipsis arguments are empty")
  }
  region_list <- list(CHR_var, BP_START_var, BP_STOP_var)
  sumstats_1_df <- do.call(sumstats_1_function, c(list(sumstats_1_file), region_list, extra_args))
  sumstats_2_df <- do.call(sumstats_2_function, c(list(sumstats_2_file), region_list, extra_args))
  coloc_output <- run_coloc(sumstats_1_df = sumstats_1_df,
                            sumstats_1_type = sumstats_1_type,
                            sumstats_1_sdY = sumstats_1_sdY,
                            sumstats_2_df = sumstats_2_df,
                            sumstats_2_type = sumstats_2_type,
                            sumstats_2_sdY = sumstats_2_sdY)
  if (add_annotation) {
    results_top <- coloc_output$results[,c("snp", "SNP.PP.H4")]
    results_top <- head(results_top[order(results_top$SNP.PP.H4, decreasing = T),], N_top_SNPs)
    Top_coloc_SNP <- paste(results_top$snp, collapse = ", ")
    SNP.PP.H4 <- results_top$SNP.PP.H4
    SNP.PP.H4 <- ifelse(SNP.PP.H4 < 0.01, format(SNP.PP.H4, scientific = TRUE, digits = 2),
                        sprintf("%.2f", round(SNP.PP.H4, 2)))
    SNP.PP.H4 <- paste(SNP.PP.H4, collapse = ", ")
    priors <- paste(coloc_output$priors, collapse=", ")
    region_df <- data.frame(sumstats_1_file = sumstats_1_file,
                            sumstats_2_file = sumstats_2_file,
                            CHR_var = CHR_var,
                            BP_START_var = BP_START_var,
                            BP_STOP_var = BP_STOP_var)
    annotation_df <- data.frame(sumstats_1_min_P = min(sumstats_1_df$P, na.rm=T),
                                sumstats_2_min_P = min(sumstats_2_df$P, na.rm=T),
                                Top_coloc_SNP = Top_coloc_SNP,
                                Top_coloc_SNP.PP.H4 = SNP.PP.H4,
                                priors = priors)
    coloc_output$summary_df <- data.frame(region_df,
                                          data.frame(t(coloc_output$summary)),
                                          annotation_df)
  }
  if (remove_full_results) {
    coloc_output$results <- NULL
  }
  return(coloc_output$summary_df)
}



#' Find significant regions in sumstats and merge close regions if needed
#' from a list of external sumstats.
#'
#' @param sumstats data frame with sumstats. Mandatory columns: CHR, BP, P
#' @param p_threshold search for regions until no more variants below this threshold remains
#' @param log_name iteration log will be written to this file
#' @return data frame with extracted regions
#' @examples
#' under development
#' sumstats=read.table('CAD.colocalization.tsv.gz', sep='\t',header=T)
#' get_coloc_regions(sumstats, p_value_name = "P.value", p_threshold = 5e-100)
#' @export
get_coloc_regions <- function(sumstats,
                              p_value_name = "P",
                              CHR_name = "CHR",
                              BP_name = "BP",
                              p_threshold = 5e-8,
                              halfwindow = 500000,
                              log_name = "log.txt") {
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
  return(coloc_regions)
}

