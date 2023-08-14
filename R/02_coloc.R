#' Runs coloc between GCKD pQTL and every sumstats
#' from a list of external sumstats.
#'
#' @param sumstats_df_1 data frame.
#' @param sumstats_1_type quant by default.
#' @param sumstats_df_2 data frame.
#' @param sumstats_2_type quant by default.
#' @return data.table with results of coloc
#' @examples
#' run_coloc()
#' @export
run_coloc <- function(sumstats_df_1, sumstats_1_type,
                      sumstats_df_2, sumstats_2_type) {
  sumstats_df_1_coloc <- list(beta=sumstats_df_1$BETA,
                              varbeta=(sumstats_df_1$SE)^2,
                              snp=sumstats_df_1$Name,
                              type=sumstats_1_type)
  if (sumstats_1_type == "quant") {
    sumstats_df_1_coloc$MAF <- sumstats_df_1$AF
    sumstats_df_1_coloc$N <- sumstats_df_1$N
  }
  sumstats_df_2_coloc <- list(beta=sumstats_df_2$BETA,
                              varbeta=(sumstats_df_2$SE)^2,
                              snp=sumstats_df_2$Name,
                              type=sumstats_2_type)
  if (sumstats_2_type == "quant") {
    sumstats_df_2_coloc$MAF <- sumstats_df_2$AF
    sumstats_df_2_coloc$N <- sumstats_df_2$N
  }
  coloc_res <- coloc.abf(dataset1=sumstats_df_1_coloc, dataset2=sumstats_df_2_coloc)
  return(coloc_res)
}


#' Find significant regions in sumstats and merge close regions if needed
#' from a list of external sumstats.
#'
#' @param sumstats data frame with sumstats. Mandatory columns: CHR, BP, P
#' @param p_threshold search for regions until no more variants below this threshold remains
#' @param log_name iteration log will be written to this file
#' @return data frame with extracted regions
#' @examples
#' run_coloc()
#' @export
get_coloc_regions <- function(sumstats,
                              p_threshold = 5e-8,
                              log_name) {
  # set up variables
  coloc_regions <- data.frame()
  region_var <- 1
  halfwindow <- 500000
  comment_var <- "PASS"
  fileConn <- file(log_name, "w")
  sink(fileConn)
  # start iterations
  while(min(sumstats$P, na.rm = T) < 5e-8) {
    min_p_row <- sumstats[which.min(sumstats$P),]
    CHR_var <- min_p_row$CHR
    BP_var <- min_p_row$BP
    BP_START_var <- BP_var - halfwindow
    BP_STOP_var <- BP_var + halfwindow
    print(paste0("Solving region ", region_var))
    print(paste0("Next most significant variant ", CHR_var, ":", BP_var))
    if (nrow(coloc_regions) > 0) {
      coloc_regions_filtered <- subset(coloc_regions, grepl("PASS", comment))
    } else {
      coloc_regions_filtered <- coloc_regions
    }
    if (CHR_var %in% coloc_regions_filtered$CHR) {
      closest_region <- subset(coloc_regions_filtered, coloc_regions_filtered$CHR == CHR_var)
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
    indeces <- which(sumstats$CHR == CHR_var & sumstats$BP >= BP_START_var & sumstats$BP <= BP_STOP_var)
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

