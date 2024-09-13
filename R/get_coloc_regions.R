#' Get coloc regions
#' @param sumstats data frame read with read_sumstats(). Mandatory columns: CHR, BP, P
#' @param p_threshold search for regions until no more variants below this threshold remains
#' @param log_name iteration log will be written to this file
#' @return data frame with extracted regions
#' @description
#' Find significant regions in sumstats and merge close regions if needed
#' from a list of external sumstats.
#' @examples
#' under development
#' @export
get_coloc_regions <- function(sumstats,
                              CHR_name = "CHR",
                              POS_name = "POS",
                              p_value_name = NULL,
                              nlog10p_value_name = "nlog10P",
                              CHR_out = "CHR_var",
                              BP_START_var_out = "BP_START_var",
                              BP_STOP_var_out = "BP_STOP_var",
                              nlogP_threshold = 7.30103,
                              halfwindow = 500000) {
  if (!is.null(p_value_name)) {
    stop(paste0("This function now works with -log10(P), please add this column to sumstats first."))
  }
  # check if there are any significant regions
  if (max(sumstats[[nlog10p_value_name]], na.rm = T) < nlogP_threshold) {
    stop("No regions below the given threshold detected")
  } else {
    message(paste0("Some significant regions below the given threshold detected, ", 
            "starting iterations to identify them."))
  }
  # create a copy of sumstats object (for subsetting in the end)
  sumstats_backup <- sumstats
  # identify chromosomes without significant hits and remove them to speed up
  sumstats <- subset_chromosomes(sumstats=sumstats, CHR_name=CHR_name,
                                 nlog10p_value_name=nlog10p_value_name,
                                 nlogP_threshold=nlogP_threshold) 
  # set up variables
  coloc_regions <- data.frame()
  regions_log <- c()
  region_var <- 1
  comment_var <- "PASS"
  # start iterations
  while(max(sumstats[[nlog10p_value_name]], na.rm = T) > nlogP_threshold) {
    # Rmpfr has potential bug with which.min, therefore a fix
    which_max <- which(sumstats[[nlog10p_value_name]] == max(sumstats[[nlog10p_value_name]]))
    if (length(which_max) > 1) {which_max <- which_max[1]}
    min_p_row <- sumstats[which_max,]
    min_p_row[[nlog10p_value_name]] <- as.numeric(min_p_row[[nlog10p_value_name]])
    print(min_p_row)
    CHR_var <- min_p_row[[CHR_name]]
    BP_var <- min_p_row[[POS_name]]
    BP_START_var <- BP_var - halfwindow
    BP_STOP_var <- BP_var + halfwindow
    regions_log <- c(regions_log, paste0("Solving region ", region_var, ": Most significant variant ", CHR_var, ":", BP_var))
    if (nrow(coloc_regions) > 0) {
      coloc_regions_filtered <- subset(coloc_regions, grepl("PASS", comment))
    } else {
      coloc_regions_filtered <- coloc_regions
    }
    if (CHR_var %in% coloc_regions_filtered[[CHR_name]]) {
      closest_regions <- subset(coloc_regions_filtered, coloc_regions_filtered[[CHR_name]] == CHR_var)
      cr1 <- sapply(closest_regions[["BP_START"]], function(x) abs(x - BP_var))
      cr2 <- sapply(closest_regions[["BP_STOP"]], function(x) abs(x - BP_var))
      cr <- c(cr1,cr2)
      closest_region <- rbind(closest_regions, closest_regions)[which.min(cr),]
      regions_log <- c(regions_log, paste0("Closest region so far: region=", closest_region$region, ", ", closest_region[[CHR_name]], ":", closest_region$BP_START, "-", closest_region$BP_STOP))
      if (abs(closest_region$BP_START - BP_var) < halfwindow |
          abs(closest_region$BP_STOP - BP_var) < halfwindow) {
        regions_log <- c(regions_log, paste0("Next most significant variant is closer than ", halfwindow, " BP to the closest region, merging"))
        if (sum(cr < halfwindow) == 1) {
          if (abs(closest_region$BP_START - BP_var) < halfwindow) {
            coloc_regions[coloc_regions$region == closest_region$region,]$BP_START <- BP_var-halfwindow
          } else if (abs(closest_region$BP_STOP - BP_var) < halfwindow) {
            coloc_regions[coloc_regions$region == closest_region$region,]$BP_STOP <- BP_var+halfwindow
          }
        } else if (sum(cr < halfwindow) == 2) {
          if (abs(closest_region$BP_START - BP_var) < halfwindow) {
            coloc_regions[coloc_regions$region == closest_region$region,]$BP_START <- BP_var
          } else if (abs(closest_region$BP_STOP - BP_var) < halfwindow) {
            coloc_regions[coloc_regions$region == closest_region$region,]$BP_STOP <- BP_var
          }
        } else { stop ("More than 2 regions closer than halfwindow identified, please check the data.")}
        updated_region <- subset(coloc_regions, region == closest_region$region)
        regions_log <- c(regions_log, paste0("Updated region: region=", updated_region$region, ", ", updated_region$CHR, ":", updated_region$BP_START, "-", updated_region$BP_STOP))
        comment_var <- paste0("SKIP_merged_to_", updated_region$region)
      } else {
        regions_log <- c(regions_log, paste0("Next most significant variant is further than ", halfwindow, " BP"))
      }
    } else {
      regions_log <- c(regions_log, paste0("No hits on this chromosome so far"))
    }
    coloc_regions <- rbind(coloc_regions, data.frame(region = region_var, min_p_row, BP_START = BP_START_var, BP_STOP = BP_STOP_var, comment = comment_var))
    old_indeces <- 1:nrow(sumstats)
    indeces <- which(sumstats[[CHR_name]] == CHR_var & sumstats[[POS_name]] >= BP_START_var & sumstats[[POS_name]] <= BP_STOP_var)
    new_indeces <- old_indeces[! old_indeces %in% indeces]
    stopifnot(nrow(sumstats) - length(indeces) == length(new_indeces))
    sumstats <- sumstats[new_indeces,]
    region_var <- region_var + 1
    comment_var <- "PASS"
    regions_log <- c(regions_log, "----------------")
    if (nrow(sumstats) == 0) {
      break
    }
  }
  if (nrow(coloc_regions) > 0) {
    # fix negative BP
    coloc_regions$BP_START[coloc_regions$BP_START < 1] <- 1
    # naming
    colnames(coloc_regions)[colnames(coloc_regions) == CHR_name] <- CHR_out
    colnames(coloc_regions)[colnames(coloc_regions) == "BP_START"] <- BP_START_var_out
    colnames(coloc_regions)[colnames(coloc_regions) == "BP_STOP"] <- BP_STOP_var_out
    start_cols <- which(colnames(coloc_regions) %in% c(CHR_out, BP_START_var_out, BP_STOP_var_out))
    end_cols <- which(!colnames(coloc_regions) %in% c(CHR_out, BP_START_var_out, BP_STOP_var_out))
    coloc_regions <- coloc_regions[,c(start_cols, end_cols)]
    rownames(coloc_regions) <- NULL
  }
  # subset sumstats
  coloc_regions_PASS <- subset(coloc_regions, comment == "PASS", c("CHR_var", "BP_START_var", "BP_STOP_var"))
  sumstats_filt_list <- lapply(1:nrow(coloc_regions_PASS), function(i) {
    subset(sumstats_backup, sumstats_backup[[CHR_name]] == coloc_regions_PASS[i,][["CHR_var"]] & 
             sumstats_backup[[POS_name]] >= coloc_regions_PASS[i,][["BP_START_var"]] & 
             sumstats_backup[[POS_name]] <= coloc_regions_PASS[i,][["BP_STOP_var"]])
  })
  sumstats_filt <- do.call(rbind, sumstats_filt_list)
  # sort
  sumstats_filt <- sumstats_filt[with(sumstats_filt, order(sumstats_filt[[CHR_name]], sumstats_filt[[POS_name]])), ]
  coloc_regions_list <- list(coloc_regions = coloc_regions,
                             coloc_regions_PASS = coloc_regions_PASS,
                             regions_log = regions_log,
                             sumstats_filt = sumstats_filt)
  # return
  return(coloc_regions_list)
}


# helpers ----
#' subset_chromosomes
#' Subset chromosomes with significant signals to speed up computation
#' Used by get_coloc_regions()
subset_chromosomes <- function(sumstats, CHR_name, nlog10p_value_name, nlogP_threshold) {
  CHR_to_keep <- unique(sumstats[[CHR_name]][sumstats[[nlog10p_value_name]] > nlogP_threshold])
  sumstats <- sumstats[sumstats[[CHR_name]] %in% CHR_to_keep,]
  message(paste0("Found significant regions on chromosomes: ", paste(CHR_to_keep, collapse = ", ")))
  return(sumstats)
}
