#' Find significant genomic regions for colocalization analysis
#'
#' @description Identifies significant regions in GWAS summary statistics and
#' merges nearby regions if needed. The function iteratively finds variants
#' below the p-value threshold and creates regions around them.
#'
#' @param sumstats Data frame containing GWAS summary statistics
#' @param CHR_name Column name for chromosome in input data (default: "CHR")
#' @param POS_name Column name for position in input data (default: "POS")
#' @param nlog10p_value_name Column name for -log10(p-value) in input data (default: "nlog10P")
#' @param CHR_out Column name for chromosome in output data (default: "CHR_var")
#' @param BP_START_var_out Column name for start position in output data (default: "BP_START_var")
#' @param BP_STOP_var_out Column name for end position in output data (default: "BP_STOP_var")
#' @param nlogP_threshold Significance threshold on -log10(p-value) scale (default: 7.30103, equals p=5e-8)
#' @param halfwindow Half-size of the window around significant variants in base pairs (default: 500000)
#' @return A list of class "coloc_regions_list" containing:
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
  message("Found significant regions on chromosomes: ", paste(sort(sig_chrs), collapse = ", "))

  # --- Sorted-peak approach: O(S log S + S*K) instead of O(N*K) ---
  # S = number of significant variants (~84K for MVP), K = number of regions
  # (~200). The old approach scanned all N (~20M) variants K times.
  # 1. Extract only significant variants, sort by nlog10P descending.
  # 2. Walk through sorted peaks; for each, check if it falls inside an
  #    already-claimed region (fast, only K regions to check).
  # 3. If not, create a new region; if within halfwindow of an existing
  #    region, merge (extend).

  sig_mask <- sumstats[[nlog10p_value_name]] > nlogP_threshold
  sig_idx  <- which(sig_mask)
  sig_ord  <- sig_idx[order(sumstats[[nlog10p_value_name]][sig_idx], decreasing = TRUE)]

  chr_regions_list <- list()
  regions_log <- character()
  global_region <- 1L

  # Per-chromosome PASS region tracking (vectors, not data.frames)
  pass_chr   <- integer(0)    # chromosome of each PASS region
  pass_start <- numeric(0)    # BP_START
  pass_stop  <- numeric(0)    # BP_STOP
  pass_id    <- integer(0)    # region id

  # Claimed intervals: track all windows (PASS + merged) for skipping
  # variants already inside an existing window. Stored per-chromosome
  # in a list of 2-column matrices for fast lookup.
  claimed <- list()  # chr -> matrix(ncol=2) of [start, stop]

  for (si in seq_along(sig_ord)) {
    row_i   <- sig_ord[si]
    CHR_var <- sumstats[[CHR_name]][row_i]
    BP_var  <- sumstats[[POS_name]][row_i]

    # Skip if this variant is already inside a claimed window
    chr_key <- as.character(CHR_var)
    if (!is.null(claimed[[chr_key]])) {
      cmat <- claimed[[chr_key]]
      if (any(BP_var >= cmat[, 1] & BP_var <= cmat[, 2])) next
    }

    BP_START_var <- BP_var - halfwindow
    BP_STOP_var  <- BP_var + halfwindow

    regions_log <- c(regions_log,
                     paste0("Solving region ", global_region,
                            ": Most significant variant ", CHR_var, ":", BP_var))

    comment_var <- "PASS"

    # Check for overlaps with existing PASS regions on this chromosome
    same_chr <- which(pass_chr == CHR_var)
    if (length(same_chr) > 0) {
      d_start <- abs(pass_start[same_chr] - BP_var)
      d_stop  <- abs(pass_stop[same_chr]  - BP_var)
      cr <- c(d_start, d_stop)
      # Find closest region using the same rbind(regions, regions)[which.min(cr)]
      # logic as the original code
      closest_local <- ((which.min(cr) - 1L) %% length(same_chr)) + 1L
      closest_idx   <- same_chr[closest_local]

      regions_log <- c(regions_log,
                       paste0("Closest region so far: region=", pass_id[closest_idx],
                              ", ", CHR_var, ":", pass_start[closest_idx],
                              "-", pass_stop[closest_idx]))

      if (abs(pass_start[closest_idx] - BP_var) < halfwindow ||
          abs(pass_stop[closest_idx]  - BP_var) < halfwindow) {
        regions_log <- c(regions_log,
                         paste0("Next most significant variant is closer than ",
                                halfwindow, " BP to the closest region, merging"))

        n_close <- sum(cr < halfwindow)
        if (n_close == 1) {
          if (abs(pass_start[closest_idx] - BP_var) < halfwindow) {
            pass_start[closest_idx] <- BP_var - halfwindow
          } else {
            pass_stop[closest_idx]  <- BP_var + halfwindow
          }
        } else if (n_close == 2) {
          if (abs(pass_start[closest_idx] - BP_var) < halfwindow) {
            pass_start[closest_idx] <- BP_var
          } else {
            pass_stop[closest_idx]  <- BP_var
          }
        } else {
          stop("More than 2 regions closer than halfwindow identified, please check the data.")
        }

        regions_log <- c(regions_log,
                         paste0("Updated region: region=", pass_id[closest_idx],
                                ", ", CHR_var, ":", pass_start[closest_idx],
                                "-", pass_stop[closest_idx]))
        comment_var <- paste0("SKIP_merged_to_", pass_id[closest_idx])
      } else {
        regions_log <- c(regions_log,
                         paste0("Next most significant variant is further than ",
                                halfwindow, " BP"))
      }
    } else {
      regions_log <- c(regions_log, "No hits on this chromosome so far")
    }

    # Record region
    peak_row <- sumstats[row_i, , drop = FALSE]
    peak_row[[nlog10p_value_name]] <- as.numeric(peak_row[[nlog10p_value_name]])

    chr_regions_list[[length(chr_regions_list) + 1L]] <- data.frame(
      region = global_region, peak_row,
      BP_START = BP_START_var, BP_STOP = BP_STOP_var,
      comment = comment_var, stringsAsFactors = FALSE)

    if (comment_var == "PASS") {
      pass_chr   <- c(pass_chr,   CHR_var)
      pass_start <- c(pass_start, BP_START_var)
      pass_stop  <- c(pass_stop,  BP_STOP_var)
      pass_id    <- c(pass_id,    global_region)
    }

    # Add this window to the claimed set so subsequent variants are skipped
    if (is.null(claimed[[chr_key]])) {
      claimed[[chr_key]] <- matrix(c(BP_START_var, BP_STOP_var), ncol = 2)
    } else {
      claimed[[chr_key]] <- rbind(claimed[[chr_key]],
                                  c(BP_START_var, BP_STOP_var))
    }

    global_region <- global_region + 1L
    regions_log <- c(regions_log, "----------------")
  }

  # Combine per-chromosome results
  if (length(chr_regions_list) == 0) {
    warning("No regions found despite significant variants")
    return(list(coloc_regions = NA, coloc_regions_PASS = NA,
                regions_log = regions_log, sumstats_filt = NA))
  }

  coloc_regions <- do.call(rbind, chr_regions_list)
  rownames(coloc_regions) <- NULL

  # Apply updated PASS boundaries back (merges may have extended them)
  for (k in seq_along(pass_id)) {
    coloc_regions$BP_START[coloc_regions$region == pass_id[k]] <- pass_start[k]
    coloc_regions$BP_STOP[coloc_regions$region == pass_id[k]]  <- pass_stop[k]
  }

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

  # Filter to PASS regions only
  coloc_regions_PASS <- subset(coloc_regions,
                               comment == "PASS",
                               c(CHR_out, BP_START_var_out, BP_STOP_var_out))

  # Extract variants in regions that passed - one pass per chromosome
  # instead of one pass per region (543 regions x 16M rows -> 22 chr x ~750K)
  keep <- rep(FALSE, nrow(sumstats_backup))
  for (chr in unique(coloc_regions_PASS[[CHR_out]])) {
    chr_rows   <- which(sumstats_backup[[CHR_name]] == chr)
    chr_pos    <- sumstats_backup[[POS_name]][chr_rows]
    chr_regions <- coloc_regions_PASS[coloc_regions_PASS[[CHR_out]] == chr, , drop = FALSE]
    for (j in seq_len(nrow(chr_regions))) {
      keep[chr_rows[chr_pos >= chr_regions[[BP_START_var_out]][j] &
                     chr_pos <= chr_regions[[BP_STOP_var_out]][j]]] <- TRUE
    }
  }
  sumstats_filt <- sumstats_backup[keep, ]

  # Sort the filtered data
  if (nrow(sumstats_filt) > 0) {
    sumstats_filt <- sumstats_filt[order(sumstats_filt[[CHR_name]],
                                         sumstats_filt[[POS_name]]), ]
  }

  # Create results list
  coloc_regions_list <- list(
    coloc_regions = coloc_regions,
    coloc_regions_PASS = coloc_regions_PASS,
    regions_log = regions_log,
    sumstats_filt = sumstats_filt
  )

  class(coloc_regions_list) <- unique(c("coloc_regions_list", class(coloc_regions_list)))

  return(coloc_regions_list)
}

#' Save colocalization regions and results to files
#'
#' @description Writes the components of a colocalization regions list to
#' separate files: log messages, all regions, passing regions, and filtered summary
#' statistics. The summary statistics are compressed and indexed with bgzip/tabix.
#'
#' @param coloc_regions_list List of class "coloc_regions_list" containing colocalization regions data
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
#' @param out_name Suffix to add to the base name (default: "_subset")
#' @param CHR Column name for chromosome (default: "CHR")
#' @param POS Column name for position (default: "POS")
#' @param SKIP_name Value for tabix's -c parameter to skip the header line. If NULL (default), uses the first column name.
#' @param order_sumstats Whether to sort the data by chromosome and position (default: FALSE)
#' @param bgzip_bin Path to bgzip executable (default: "bgzip")
#' @param tabix_bin Path to tabix executable (default: "tabix")
#' @return Path to the compressed file (.tsv.gz)
#' @importFrom data.table fwrite
#' @export
gc_bgzip_tabix <- function(sumstats, 
                           sumstats_name, 
                           out_name = "_subset",
                           CHR = "CHR", 
                           POS = "POS",
                           SKIP_name = NULL,
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
  
  # Auto-detect SKIP_name from the first column if not provided
  if (is.null(SKIP_name)) {
    SKIP_name <- colnames(sumstats)[1]
  }

  # Get column positions (needed for tabix)
  CHR_place <- which(colnames(sumstats) == CHR)
  POS_place <- which(colnames(sumstats) == POS)

  message("CHR column found at position ", CHR_place,
          ". POS column found at position ", POS_place)
  
  # Order data if requested
  if (order_sumstats) {
    message("Ordering sumstats by CHR and POS")
    new_order <- order(sumstats[[CHR]], sumstats[[POS]])
    sumstats <- sumstats[new_order, ]
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

