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
  
  # Create results list
  coloc_regions_list <- list(
    coloc_regions = coloc_regions,
    coloc_regions_PASS = coloc_regions_PASS,
    regions_log = regions_log,
    sumstats_filt = sumstats_filt
  )
  
  # Add class to the results list
  class(coloc_regions_list) <- unique(c("coloc_regions_list", class(coloc_regions_list)))
  
  return(coloc_regions_list)
}

#' Save colocalization regions and results to files
#'
#' @description Writes the components of a colocalization regions list to
#' separate files: log messages, all regions, passing regions, and filtered summary
#' statistics. The summary statistics are compressed and indexed with bgzip/tabix.
#'
#' @param coloc_regions_list List containing colocalization regions data
#' @param sumstats_name Base name for output files
#' @return Path to the compressed summary statistics file if created, otherwise NULL
#' @export
write_regions <- function(coloc_regions_list, sumstats_name) {
  # Check that coloc_regions_list has the expected class
  if (!inherits(coloc_regions_list, "coloc_regions_list") && 
      !all(c("regions_log", "coloc_regions", "coloc_regions_PASS", "sumstats_filt") %in% names(coloc_regions_list))) {
    stop("Input must be a coloc_regions_list object containing required components")
  }
  
  # Initialize return value
  sumstats_file <- NULL
  
  # Check for NA values
  components_na <- sapply(coloc_regions_list, function(x) is.null(x) || all(is.na(x)))
  
  # Write log file if it exists
  if (!components_na["regions_log"]) {
    log_file <- paste0(sumstats_name, "_log.txt")
    writeLines(coloc_regions_list[["regions_log"]], con = log_file)
    message("Written: ", log_file)
  }
  
  # Write all regions if they exist
  if (!components_na["coloc_regions"]) {
    regions_file <- paste0(sumstats_name, "_coloc_regions.tsv")
    write.table(
      coloc_regions_list[["coloc_regions"]], 
      file = regions_file,
      sep = "\t", 
      row.names = FALSE, 
      col.names = TRUE, 
      quote = FALSE
    )
    message("Written: ", regions_file)
  }
  
  # Write passing regions if they exist
  if (!components_na["coloc_regions_PASS"]) {
    pass_file <- paste0(sumstats_name, "_coloc_regions_PASS.tsv")
    write.table(
      coloc_regions_list[["coloc_regions_PASS"]], 
      file = pass_file,
      sep = "\t", 
      row.names = FALSE, 
      col.names = TRUE, 
      quote = FALSE
    )
    message("Written: ", pass_file)
  }
  
  # Write and compress filtered summary statistics if they exist
  if (!components_na["sumstats_filt"]) {
    sumstats_file <- gc_bgzip_tabix(
      sumstats = coloc_regions_list[["sumstats_filt"]],
      sumstats_name = sumstats_name
    )
    message("Written and indexed: ", sumstats_file)
  } else {
    message("No filtered summary statistics to write")
  }
  
  return(sumstats_file)
}


#' Compress and index summary statistics with bgzip and tabix
#'
#' @description Writes summary statistics to a file, compresses it with bgzip,
#' and creates a tabix index for efficient region-based queries.
#'
#' @param sumstats Data frame containing summary statistics
#' @param sumstats_name Base name for output file
#' @param out_name Suffix to add to the base name
#' @param CHR Column name for chromosome
#' @param POS Column name for position
#' @param SKIP_name Value for tabix's -c parameter (comment character or column to skip)
#' @param order_sumstats Whether to sort the data by chromosome and position
#' @param bgzip_bin Path to bgzip executable
#' @param tabix_bin Path to tabix executable
#' @return Path to the compressed file
#' @importFrom data.table fwrite
#' @export
gc_bgzip_tabix <- function(sumstats, 
                           sumstats_name, 
                           out_name = "_subset",
                           CHR = "CHR", 
                           POS = "POS",
                           SKIP_name = "Name", 
                           order_sumstats = FALSE,
                           bgzip_bin = "bgzip", 
                           tabix_bin = "tabix") {
  # Check if required executables are available
  check_bin(bgzip_bin)
  check_bin(tabix_bin)
  
  # Check required columns
  if (!CHR %in% colnames(sumstats)) {
    stop("CHR column '", CHR, "' not found in sumstats")
  }
  if (!POS %in% colnames(sumstats)) {
    stop("POS column '", POS, "' not found in sumstats")
  }
  
  # Get column positions (needed for tabix)
  CHR_place <- which(colnames(sumstats) == CHR)
  POS_place <- which(colnames(sumstats) == POS)
  
  message("CHR column found at position ", CHR_place, 
          ". POS column found at position ", POS_place)
  
  # Order data if requested
  if (order_sumstats) {
    message("Ordering sumstats by CHR and POS")
    sumstats <- sumstats[order(sumstats[[CHR]], sumstats[[POS]]), ]
  }
  
  # Create file paths
  tsv_file <- paste0(sumstats_name, out_name, ".tsv")
  gz_file <- paste0(tsv_file, ".gz")
  
  # Write data to file
  message("Writing to ", tsv_file)
  data.table::fwrite(sumstats, tsv_file, sep = "\t", scipen = 999)
  
  # Compress with bgzip
  bgzip_cmd <- paste0(bgzip_bin, " -f ", tsv_file)
  message("Compressing with bgzip: ", bgzip_cmd)
  result <- system(bgzip_cmd)
  if (result != 0) {
    stop("bgzip compression failed with exit code ", result)
  }
  
  # Index with tabix
  tabix_cmd <- paste0(tabix_bin, " -f -s", CHR_place, " -b", POS_place, 
                      " -e", POS_place, " ", gz_file, " -c ", SKIP_name)
  message("Creating tabix index: ", tabix_cmd)
  result <- system(tabix_cmd)
  if (result != 0) {
    stop("tabix indexing failed with exit code ", result)
  }
  
  message("Successfully created ", gz_file, " and index")
  
  return(gz_file)
}

#' Check if a binary is available in the system path
#'
#' @param bin_name Name or path of the binary to check
#' @return TRUE if the binary is available, otherwise stops with an error
#' @keywords internal
check_bin <- function(bin_name) {
  # Use 'which' on Unix/Linux/Mac or 'where' on Windows
  cmd <- if (.Platform$OS.type == "windows") {
    paste0("where ", bin_name, " 2>NUL")
  } else {
    paste0("which ", bin_name, " 2>/dev/null")
  }
  
  # Run command and check result
  result <- system(cmd, intern = TRUE)
  
  if (length(result) == 0) {
    stop("Required binary '", bin_name, "' not found in system path. ",
         "Please install it or provide the full path.")
  }
  
  return(TRUE)
}
