#' Find significant genomic regions for colocalization analysis
#'
#' @description Identifies significant regions in GWAS summary statistics and
#' merges nearby regions if needed. The function iteratively finds variants
#' below the p-value threshold and creates regions around them.
#'
#' @param sumstats Data frame containing GWAS summary statistics
#' @param CHR_name Column name for chromosome in input data
#' @param POS_name Column name for position in input data
#' @param nlog10p_value_name Column name for -log10(p-value) in input data
#' @param CHR_out Column name for chromosome in output data
#' @param BP_START_var_out Column name for start position in output data
#' @param BP_STOP_var_out Column name for end position in output data
#' @param nlogP_threshold Significance threshold on -log10(p-value) scale (default: 7.30103, equals p=5e-8)
#' @param halfwindow Half-size of the window around significant variants in base pairs
#' @return A list containing:
#'   \item{coloc_regions}{Data frame with all identified regions}
#'   \item{coloc_regions_PASS}{Data frame with regions that passed filters}
#'   \item{regions_log}{Character vector with log messages}
#'   \item{sumstats_filt}{Filtered summary statistics in the identified regions}
#' @export
get_coloc_regions <- function(sumstats,
                              CHR_name = "CHR",
                              POS_name = "POS",
                              nlog10p_value_name = "nlog10P",
                              CHR_out = "CHR_var",
                              BP_START_var_out = "BP_START_var",
                              BP_STOP_var_out = "BP_STOP_var",
                              nlogP_threshold = 7.30103,
                              halfwindow = 500000) {
  # Validate input columns
  required_cols <- c(CHR_name, POS_name, nlog10p_value_name)
  missing_cols <- required_cols[!required_cols %in% colnames(sumstats)]
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for significant signals
  max_signal <- max(sumstats[[nlog10p_value_name]], na.rm = TRUE)
  if (max_signal < nlogP_threshold) {
    warning("No regions below the significance threshold (", nlogP_threshold, ") detected")
    return(list(
      coloc_regions = NA,
      coloc_regions_PASS = NA,
      regions_log = NA,
      sumstats_filt = NA
    ))
  }
  
  message("Found significant regions with -log10(p) > ", nlogP_threshold, 
          ", starting iterations to identify them.")
  
  # Keep original data for final filtering
  sumstats_backup <- sumstats
  
  # Subset to chromosomes with significant signals
  sig_chrs <- unique(sumstats[[CHR_name]][sumstats[[nlog10p_value_name]] > nlogP_threshold])
  sumstats <- sumstats[sumstats[[CHR_name]] %in% sig_chrs, ]
  message("Found significant regions on chromosomes: ", paste(sig_chrs, collapse = ", "))
  
  # Initialize variables
  coloc_regions <- data.frame()
  regions_log <- character()
  region_var <- 1
  comment_var <- "PASS"
  
  # Iteratively find regions
  while (nrow(sumstats) > 0 && max(sumstats[[nlog10p_value_name]], na.rm = TRUE) > nlogP_threshold) {
    # Find the most significant variant
    which_max <- which(sumstats[[nlog10p_value_name]] == max(sumstats[[nlog10p_value_name]]))
    if (length(which_max) > 1) which_max <- which_max[1]
    min_p_row <- sumstats[which_max, ]
    min_p_row[[nlog10p_value_name]] <- as.numeric(min_p_row[[nlog10p_value_name]])
    print(min_p_row)
    
    # Extract key information
    CHR_var <- min_p_row[[CHR_name]]
    BP_var <- min_p_row[[POS_name]]
    BP_START_var <- BP_var - halfwindow
    BP_STOP_var <- BP_var + halfwindow
    
    # Log the current step - KEEPING ORIGINAL LOG FORMAT
    regions_log <- c(regions_log, 
                     paste0("Solving region ", region_var, ": Most significant variant ", CHR_var, ":", BP_var))
    
    # Check for overlaps with existing regions
    if (nrow(coloc_regions) > 0) {
      coloc_regions_filtered <- subset(coloc_regions, grepl("PASS", comment))
      
      if (CHR_var %in% coloc_regions_filtered[[CHR_name]]) {
        closest_regions <- subset(coloc_regions_filtered, coloc_regions_filtered[[CHR_name]] == CHR_var)
        cr1 <- sapply(closest_regions[["BP_START"]], function(x) abs(x - BP_var))
        cr2 <- sapply(closest_regions[["BP_STOP"]], function(x) abs(x - BP_var))
        cr <- c(cr1, cr2)
        closest_region <- rbind(closest_regions, closest_regions)[which.min(cr), ]
        
        regions_log <- c(regions_log, 
                         paste0("Closest region so far: region=", closest_region$region, 
                                ", ", closest_region[[CHR_name]], ":", 
                                closest_region$BP_START, "-", closest_region$BP_STOP))
        
        if (abs(closest_region$BP_START - BP_var) < halfwindow || 
            abs(closest_region$BP_STOP - BP_var) < halfwindow) {
          regions_log <- c(regions_log, 
                           paste0("Next most significant variant is closer than ", 
                                  halfwindow, " BP to the closest region, merging"))
          
          if (sum(cr < halfwindow) == 1) {
            if (abs(closest_region$BP_START - BP_var) < halfwindow) {
              coloc_regions[coloc_regions$region == closest_region$region, ]$BP_START <- BP_var - halfwindow
            } else if (abs(closest_region$BP_STOP - BP_var) < halfwindow) {
              coloc_regions[coloc_regions$region == closest_region$region, ]$BP_STOP <- BP_var + halfwindow
            }
          } else if (sum(cr < halfwindow) == 2) {
            if (abs(closest_region$BP_START - BP_var) < halfwindow) {
              coloc_regions[coloc_regions$region == closest_region$region, ]$BP_START <- BP_var
            } else if (abs(closest_region$BP_STOP - BP_var) < halfwindow) {
              coloc_regions[coloc_regions$region == closest_region$region, ]$BP_STOP <- BP_var
            }
          } else {
            stop("More than 2 regions closer than halfwindow identified, please check the data.")
          }
          
          updated_region <- subset(coloc_regions, region == closest_region$region)
          regions_log <- c(regions_log, 
                           paste0("Updated region: region=", updated_region$region, 
                                  ", ", updated_region$CHR, ":", 
                                  updated_region$BP_START, "-", updated_region$BP_STOP))
          
          comment_var <- paste0("SKIP_merged_to_", updated_region$region)
        } else {
          regions_log <- c(regions_log, 
                           paste0("Next most significant variant is further than ", halfwindow, " BP"))
        }
      } else {
        regions_log <- c(regions_log, "No hits on this chromosome so far")
      }
    } else {
      regions_log <- c(regions_log, "No hits on this chromosome so far")
    }
    
    # Add the new region
    coloc_regions <- rbind(coloc_regions, 
                           data.frame(region = region_var, 
                                      min_p_row, 
                                      BP_START = BP_START_var, 
                                      BP_STOP = BP_STOP_var, 
                                      comment = comment_var))
    
    # Remove variants in this region from further consideration
    old_indeces <- 1:nrow(sumstats)
    indeces <- which(sumstats[[CHR_name]] == CHR_var & 
                       sumstats[[POS_name]] >= BP_START_var & 
                       sumstats[[POS_name]] <= BP_STOP_var)
    new_indeces <- old_indeces[!old_indeces %in% indeces]
    stopifnot(nrow(sumstats) - length(indeces) == length(new_indeces))
    sumstats <- sumstats[new_indeces, ]
    
    # Increment region counter and add separator to log
    region_var <- region_var + 1
    comment_var <- "PASS"
    regions_log <- c(regions_log, "----------------")
    
    # Check if we've processed all variants
    if (nrow(sumstats) == 0) {
      break
    }
  }
  
  # Process results
  if (nrow(coloc_regions) > 0) {
    # Ensure start positions are valid
    coloc_regions$BP_START[coloc_regions$BP_START < 1] <- 1
    
    # Rename columns to output names
    colnames(coloc_regions)[colnames(coloc_regions) == CHR_name] <- CHR_out
    colnames(coloc_regions)[colnames(coloc_regions) == "BP_START"] <- BP_START_var_out
    colnames(coloc_regions)[colnames(coloc_regions) == "BP_STOP"] <- BP_STOP_var_out
    
    # Reorder columns for better readability
    start_cols <- which(colnames(coloc_regions) %in% c(CHR_out, BP_START_var_out, BP_STOP_var_out))
    end_cols <- which(!colnames(coloc_regions) %in% c(CHR_out, BP_START_var_out, BP_STOP_var_out))
    coloc_regions <- coloc_regions[, c(start_cols, end_cols)]
    
    # Reset row names
    rownames(coloc_regions) <- NULL
  }
  
  # Filter to PASS regions only
  coloc_regions_PASS <- subset(coloc_regions, 
                               comment == "PASS", 
                               c(CHR_out, BP_START_var_out, BP_STOP_var_out))
  
  # Extract variants in regions that passed
  sumstats_filt_list <- lapply(1:nrow(coloc_regions_PASS), function(i) {
    subset(sumstats_backup, 
           sumstats_backup[[CHR_name]] == coloc_regions_PASS[i, ][[CHR_out]] & 
             sumstats_backup[[POS_name]] >= coloc_regions_PASS[i, ][[BP_START_var_out]] & 
             sumstats_backup[[POS_name]] <= coloc_regions_PASS[i, ][[BP_STOP_var_out]])
  })
  
  sumstats_filt <- do.call(rbind, sumstats_filt_list)
  
  # Sort the filtered data
  if (!is.null(sumstats_filt) && nrow(sumstats_filt) > 0) {
    sumstats_filt <- sumstats_filt[with(sumstats_filt, 
                                        order(sumstats_filt[[CHR_name]], 
                                              sumstats_filt[[POS_name]])), ]
  }
  
  # Return results
  coloc_regions_list <- list(
    coloc_regions = coloc_regions,
    coloc_regions_PASS = coloc_regions_PASS,
    regions_log = regions_log,
    sumstats_filt = sumstats_filt
  )
  
  return(coloc_regions_list)
}