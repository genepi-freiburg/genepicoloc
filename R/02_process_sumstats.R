#' Read and format input summary statistics (under development)
#' @param sumstats_file path to sumstats file.
#' @param sumstats sumstats data.frame loaded into R.
#' @param read_method "read.delim" or "data.table". If data.table package is
#' available, "data.table" can speed up reading (usually 2-3 times faster)
#' @param sep_var Separator type (',', '\t', ' ', etc). Should be provided only if 
#' read_method is read.delim (data.table usually determined separator automatically)
#' @param Name_name update development
#' @param CHR_name update development
#' @param POS_name update development
#' @param A1_name update development
#' @param A2_name update development
#' @param BETA_name update development
#' @param SE_name update development
#' @param p_value_name update development
#' @param AF_name update development
#' @param N_name update development
#' @param other_columns vector with names of other columns
#' @return data frame with formatted columns
#' @examples
#' under development
#' @export
read_sumstats <- function(sumstats_file, sumstats, sep_var,
                          Name_name = NULL, Name_name_new = "Name",
                          rsID_name = NULL, rsID_name_new = "rsID",
                          CHR_name = NULL, CHR_name_new = "CHR",
                          POS_name = NULL, POS_name_new = "POS",
                          A1_name = NULL, A1_name_new = "A1",
                          A2_name = NULL, A2_name_new = "A2",
                          BETA_name = NULL, BETA_name_new = "BETA",
                          SE_name = NULL, SE_name_new = "SE",
                          p_value_name = NULL, p_value_name_new = "P",
                          AF_name = NULL, AF_name_new = "AF",
                          N_name = NULL, N_name_new = "N",
                          other_columns = NULL) {
  if ("data.table" %in% rownames(installed.packages())) {
    read_method <- "data.table"
  } else { 
    read_method <- "read.delim"
    if (missing(sep_var)) { stop("Please provide separator type (',', '\t', ' ', etc)") }
  }
  ### Checks
  if (missing(CHR_name)) { warning("No chromosome column provided") }
  if (missing(POS_name)) { warning("No base pair position column provided") }
  if (missing(p_value_name)) { warning("No p-value column provided") }
  ### Read
  if (!missing(sumstats_file)) {
    print("'sumstats_file' provided, sumstats will be read from the file")
    if (read_method == "data.table") {
      sumstats <- data.table::fread(input = sumstats_file)
      sumstats <- as.data.frame(sumstats)
    } else if (read_method == "read.delim") {
      sumstats <- read.delim(file = sumstats_file, header = T, sep = sep_var)
    }
  } else if (!missing(sumstats)) {
    print(paste0("'sumstats' object provided, it will be used for formatting. ",
                 "If you like to read a file from the disk, please use 'sumstats_file' argument instead."))
    sumstats <- as.data.frame(sumstats)
  }
  ### Select
  all_cols <- c(Name_name, rsID_name, CHR_name, POS_name,
                A1_name, A2_name, BETA_name, SE_name,
                p_value_name, AF_name, N_name, other_columns)
  if (!all(all_cols %in% colnames(sumstats))) {
    stop("Not all provided colnames match colnames in the sumstats")
  }
  sumstats <- sumstats[,all_cols]
  ### Format
  if (!is.null(Name_name)) {
    colnames(sumstats)[colnames(sumstats) == Name_name] <- Name_name_new
  } else {
    warning("Name is empty, it will be created using chromosome, BP position, and alleles")
    sumstats[[Name_name_new]] <- paste(sumstats[[CHR_name]],
                                       sumstats[[POS_name]],
                                       sumstats[[A2_name]],
                                       sumstats[[A1_name]], sep = ":")
  }
  if (!is.null(rsID_name)) {
    colnames(sumstats)[colnames(sumstats) == rsID_name] <- rsID_name_new
  }
  if (!is.null(CHR_name)) {
    colnames(sumstats)[colnames(sumstats) == CHR_name] <- CHR_name_new
  }
  if (!is.null(POS_name)) {
    if (class(sumstats[[POS_name]]) != "integer") {
      warning("POS column is not of the 'integer' class, converting to integer")
      sumstats[[POS_name]] <- as.integer(sumstats[[POS_name]])
    }
    colnames(sumstats)[colnames(sumstats) == POS_name] <- POS_name_new
  }
  if (!is.null(A1_name)) {
    colnames(sumstats)[colnames(sumstats) == A1_name] <- A1_name_new
  }
  if (!is.null(A2_name)) {
    colnames(sumstats)[colnames(sumstats) == A2_name] <- A2_name_new
  }
  if (!is.null(BETA_name)) {
    if (class(sumstats[[BETA_name]]) != "numeric") {
      warning("Beta column is not of class 'numeric', converting to numeric")
      sumstats[[BETA_name]] <- as.numeric(sumstats[[BETA_name]])
    }
    colnames(sumstats)[colnames(sumstats) == BETA_name] <- BETA_name_new
  }
  if (!is.null(SE_name)) {
    if (class(sumstats[[SE_name]]) != "numeric") {
      warning("SE column is not of class 'numeric', converting to numeric")
      sumstats[[SE_name]] <- as.numeric(sumstats[[SE_name]])
    }
    colnames(sumstats)[colnames(sumstats) == SE_name] <- SE_name_new
  }
  if (!is.null(p_value_name)) {
    if (class(sumstats[[p_value_name]]) != "numeric") {
      warning("The p-value column is not of the 'numeric' class, converting to numeric")
      sumstats[[p_value_name]] <- as.numeric(sumstats[[p_value_name]])
    }
    colnames(sumstats)[colnames(sumstats) == p_value_name] <- p_value_name_new
  }
  if (!is.null(AF_name)) {
    colnames(sumstats)[colnames(sumstats) == AF_name] <- AF_name_new
  }
  if (!is.null(N_name)) {
    colnames(sumstats)[colnames(sumstats) == N_name] <- N_name_new
  }
  return(sumstats)
}

#' Get coloc regions
#' @param sumstats data frame read with read_sumstats(). Mandatory columns: CHR, BP, P
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
                              POS_name = "POS",
                              p_value_name = "P",
                              CHR_out = "CHR_var",
                              BP_START_var_out = "BP_START_var",
                              BP_STOP_var_out = "BP_STOP_var",
                              p_threshold = 5e-8,
                              halfwindow = 500000) {
  # creat a copy of sumstats object (for subsetting in the end)
  sumstats_backup <- sumstats
  if (is.character(sumstats[[p_value_name]])) {
    warning(paste0(p_value_name, " column is character, converting to numeric"))
    if (!"Rmpfr" %in% rownames(installed.packages())) {
      stop("Rmpfr is required to run this function")
    }
    sumstats[[p_value_name]] <- Rmpfr::mpfr(sumstats[[p_value_name]])
  }
  # set up variables
  coloc_regions <- data.frame()
  regions_log <- c()
  # function-specific constants
  region_var <- 1
  comment_var <- "PASS"
  if (min(sumstats[[p_value_name]], na.rm = T) > p_threshold) {
    stop("No regions below the given threshold detected")
  }
  # start iterations
  while(min(sumstats[[p_value_name]], na.rm = T) < p_threshold) {
    regions_log <- c(regions_log, paste0("Solving region ", region_var))
    print(paste0("Solving region ", region_var))
    min_p_row <- sumstats[which.min(sumstats[[p_value_name]]),]
    print(min_p_row)
    CHR_var <- min_p_row[[CHR_name]]
    BP_var <- min_p_row[[POS_name]]
    BP_START_var <- BP_var - halfwindow
    BP_STOP_var <- BP_var + halfwindow
    regions_log <- c(regions_log, paste0("Next most significant variant ", CHR_var, ":", BP_var))
    if (nrow(coloc_regions) > 0) {
      coloc_regions_filtered <- subset(coloc_regions, grepl("PASS", comment))
    } else {
      coloc_regions_filtered <- coloc_regions
    }
    if (CHR_var %in% coloc_regions_filtered[[CHR_name]]) {
      closest_region <- subset(coloc_regions_filtered, coloc_regions_filtered[[CHR_name]] == CHR_var)
      closest_region <- closest_region[which.min(sapply(closest_region[[POS_name]], function(x) abs(x - BP_var))),]
      regions_log <- c(regions_log, paste0("Closest region so far: region=", closest_region$region, ", ", closest_region[[CHR_name]], ":", closest_region$BP_START, "-", closest_region$BP_STOP))
      if (abs(closest_region$BP_START - BP_var) < halfwindow |
          abs(closest_region$BP_STOP - BP_var) < halfwindow) {
        regions_log <- c(regions_log, paste0("Next most significant variant is closer than ", halfwindow, " BP to the closest region, merging"))
        if (abs(closest_region$BP_START - BP_var) < halfwindow) {
          coloc_regions[coloc_regions$region == closest_region$region,]$BP_START <- BP_var-halfwindow
        } else if (abs(closest_region$BP_STOP - BP_var) < halfwindow) {
          coloc_regions[coloc_regions$region == closest_region$region,]$BP_STOP <- BP_var+halfwindow
        }
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
  }
  if (nrow(coloc_regions) > 0) {
    # fix negative BP
    coloc_regions$BP_START[coloc_regions$BP_START < 1] <- 1
    # naming
    colnames(coloc_regions)[colnames(coloc_regions) == "CHR"] <- CHR_out
    colnames(coloc_regions)[colnames(coloc_regions) == "BP_START"] <- BP_START_var_out
    colnames(coloc_regions)[colnames(coloc_regions) == "BP_STOP"] <- BP_STOP_var_out
    start_cols <- which(colnames(coloc_regions) %in% c(CHR_out, BP_START_var_out, BP_STOP_var_out))
    end_cols <- which(!colnames(coloc_regions) %in% c(CHR_out, BP_START_var_out, BP_STOP_var_out))
    coloc_regions <- coloc_regions[,c(start_cols, end_cols)]
    rownames(coloc_regions) <- NULL
  }
  # subset sumstats
  coloc_regions_PASS <- subset(coloc_regions, comment == "PASS")
  sumstats_filt_list <- lapply(1:nrow(coloc_regions_PASS), function(i) {
    subset(sumstats_backup, sumstats_backup[[CHR_name]] == coloc_regions_PASS[i,][["CHR_var"]] & 
             sumstats_backup[[POS_name]] >= coloc_regions_PASS[i,][["BP_START_var"]] & 
             sumstats_backup[[POS_name]] <= coloc_regions_PASS[i,][["BP_STOP_var"]])
  })
  sumstats_filt <- do.call(rbind, sumstats_filt_list)
  # return
  return(list(coloc_regions = coloc_regions,
              coloc_regions_PASS = coloc_regions_PASS,
              regions_log = regions_log,
              sumstats_filt = sumstats_filt))
}


#' subset_sumstats
#' @description under development
#' @export
subset_sumstats <- function(sumstats,
                            coloc_regions,
                            CHR_var, BP_START_var, BP_STOP_var,
                            CHR_name = "CHR",
                            POS_name = "POS",
                            dbSNP_file = NULL,
                            CHR_var_col_name = NULL,
                            BP_START_var_col_name = NULL,
                            BP_STOP_var_col_name = NULL,
                            remove_duplicates = T,
                            do_match_rs = F, N_mc.cores = 10, ...) {
  if (!missing(coloc_regions)) {
    sumstats_filt <- do.call(rbind, apply(coloc_regions, 1, function(row) {
      sumstats_filt <- subset(sumstats, sumstats[[CHR_name]] == row["CHR_var"] & 
                                sumstats[[POS_name]] >= row["BP_START_var"] & 
                                sumstats[[POS_name]] <= row["BP_STOP_var"])
    }))
  } else {
    sumstats_filt <- subset(sumstats, sumstats[[CHR_name]] == CHR_var & 
                              sumstats[[POS_name]] >= BP_START_var & 
                              sumstats[[POS_name]] <= BP_STOP_var)
  }
  if (do_match_rs) {
    if ("parallel" %in% rownames(installed.packages())) {
      sumstats_filt_list <- parallel::mclapply(1:nrow(sumstats_filt), function(i) {
        do.call(match_rs, list(dbSNP_file = "/data/public_resources/Ensembl_human_variation_b38_v109/dbSNP_v156_b38p14_rsid.vcf.gz",
                               sumstats = sumstats_filt[i,],
                               CHR_var_col_name = CHR_var_col_name,
                               BP_START_var_col_name = BP_START_var_col_name,
                               BP_STOP_var_col_name = BP_STOP_var_col_name,
                               A1_name = "A1",
                               A2_name = "A2",
                               SNP_name = "rsID"))
      }, mc.cores = N_mc.cores)
      sumstats_filt <- do.call(rbind, sumstats_filt_list)
    } else {
      stop("match_rs without parallel::mclappy not yet implemented ")
    }
  }
  if (remove_duplicates) {
    sumstats_filt <- unique(sumstats_filt)
  }
  return(sumstats_filt)
}







#' process sumstats
#' @description under development
#' @export
process_sumstats_1 <- function(sumstats_file, sumstats_name,
                             Name_name,
                             rsID_name = NULL,
                             CHR_name, POS_name,
                             A1_name, A2_name,
                             BETA_name, SE_name, p_value_name,
                             AF_name, N_name, p_threshold, format_columns=F) {
  sumstats_name_out <- paste0("subset/", sumstats_name)
  system(paste0("mkdir -p ", sumstats_name_out))
  if (format_columns) {
    if(grepl("gz$", sumstats_file)) {
      sumstats_file <- paste0("zcat ", sumstats_file, " | cut -f 1-16")
    } else {
      sumstats_file <- paste0("cut -f 1-16 ", sumstats_file)
    }
  }
  sumstats_1 <- read_sumstats(sumstats_file,
                              read_method = "data.table",
                              format_columns = format_columns,
                              Name_name = Name_name,
                              rsID_name = rsID_name,
                              CHR_name = CHR_name,
                              POS_name = POS_name,
                              A1_name = A1_name,
                              A2_name = A2_name,
                              BETA_name = BETA_name,
                              SE_name = SE_name,
                              p_value_name = p_value_name,
                              AF_name = AF_name,
                              N_name = N_name)
  coloc_regions_list <- get_coloc_regions(sumstats = sumstats_1,
                                          p_threshold = p_threshold,
                                          halfwindow = 500000)
  coloc_regions <- coloc_regions_list$coloc_regions
  if (nrow(coloc_regions) == 0) {
    writeLines("No significant regions found", con = paste0(sumstats_name_out, "/log_no_significant.txt"))
    return("No significant regions")
  }
  regions_log <- coloc_regions_list$regions_log
  extra_args <- list(do_match_rs = F,
                     sumstats = sumstats_1,
                     CHR_name = "CHR",
                     POS_name = "POS")
  sumstats_1_subset <- mclapply(1:nrow(coloc_regions), 
                                function(i) do.call(subset_sumstats, c(coloc_regions[i,], extra_args)),
                                mc.cores = 8)
  sumstats_1_subset <- unique(do.call(rbind, sumstats_1_subset))
  # save
  writeLines(regions_log, con = paste0(sumstats_name_out, "/", basename(sumstats_name_out), "_log.txt"))
  saveRDS(sumstats_1_subset, paste0(sumstats_name_out, "/", basename(sumstats_name_out), "_subset.RDS"))
  saveRDS(coloc_regions, paste0(sumstats_name_out, "/", basename(sumstats_name_out), "_coloc_regions.RDS"))
  return("PASS")
}





