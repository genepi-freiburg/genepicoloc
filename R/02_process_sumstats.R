#' Read and format input summary statistics (under development)
#' @param sumstats_file path to sumstats file.
#' @param sumstats sumstats data.frame loaded into R.
#' @param read_method "read.delim" or "data.table". If data.table package is
#' available, "data.table" can speed up reading (usually 2-3 times faster)
#' @param sep_var Separator type (',', '\t', ' ', etc). Should be provided only if 
#' read_method is read.delim (data.table usually determined separator automatically)
#' @param Name update development
#' @param CHR update development
#' @param POS update development
#' @param A1 update development
#' @param A2 update development
#' @param BETA update development
#' @param SE update development
#' @param nlog10p_value update development
#' @param AF update development
#' @param N update development
#' @param other_columns vector with names of other columns
#' @return data frame with formatted columns
#' @examples
#' under development
#' @export
read_sumstats <- function(sumstats, sumstats_file=NULL,
                          Name, Name_new = "Name",
                          A1, A1_new = "A1",
                          BETA, BETA_new = "BETA",
                          SE, SE_new = "SE",
                          nlog10p_value, nlog10p_value_new = "nlog10P",
                          rsID = NULL, rsID_new = "rsID",
                          CHR = NULL, CHR_new = "CHR",
                          POS = NULL, POS_new = "POS",
                          A2 = NULL, A2_new = "A2",
                          AF = NULL, AF_new = "AF",
                          N = NULL, N_new = "N",
                          p_value=NULL,
                          other_columns = NULL) {
  if (!is.null(p_value)) {
    stop(paste("This function now accepts only nlog10p_value, please add it to sumstats first",
               " Please ensure that underflow is handled properly (use handle_overflow() if needed)."))
  }
  if (!is.null(sumstats_file)) {
    stop(paste0("The function supports only 'sumstats' input.",
                "Please read the file into memory first and then pass it here."))
  }
  # if ("data.table" %in% rownames(installed.packages())) {
  #   read_method <- "data.table"
  # } else { 
  #   read_method <- "read.delim"
  #   if (missing(sep_var)) { stop("Please provide separator type (',', '\t', ' ', etc)") }
  # }
  # deprecated: "sep_var"
  ### Checks
  # TODO identify necessary and optional columns
  if (missing(Name)) {stop("Name is required as CHR:POS:REF:ALT, please create it first.")}
  if (missing(A1)) {stop("A1 (effect allele) is required.")}
  if (missing(BETA)) {stop("BETA is required.")}
  if (missing(SE)) {stop("SE is required.")}
  if (missing(nlog10p_value)) {stop("nlog10p_value is required.")}
  ### Read
  # if (!missing(sumstats_file)) {
  #   print("'sumstats_file' provided, sumstats will be read from the file")
  #   if (read_method == "data.table") {
  #     sumstats <- data.table::fread(input = sumstats_file)
  #     sumstats <- as.data.frame(sumstats)
  #   } else if (read_method == "read.delim") {
  #     sumstats <- read.delim(file = sumstats_file, header = T, sep = sep_var)
  #   }
  # } else if (!missing(sumstats)) {
  # print(paste0("'sumstats' object provided, it will be used for formatting. ",
  #              "If you like to read a file from the disk, please use 'sumstats_file' argument instead."))
  sumstats <- as.data.frame(sumstats)
  # }
  ### Select
  all_cols <- c(Name, rsID, CHR, POS,
                A1, A2, BETA, SE,
                nlog10p_value, AF, N, other_columns)
  if (!all(all_cols %in% colnames(sumstats))) {
    stop("Not all provided colnames match colnames in the sumstats")
  }
  sumstats <- sumstats[,all_cols]
  ### Mandatory columns
  colnames(sumstats)[colnames(sumstats) == Name] <- Name_new
  colnames(sumstats)[colnames(sumstats) == A1] <- A1_new
  if (class(sumstats[[BETA]]) != "numeric") {
    stop("Beta column is not of class 'numeric', please check the data")
  }
  colnames(sumstats)[colnames(sumstats) == BETA] <- BETA_new
  if (class(sumstats[[SE]]) != "numeric") {
    warning("SE column is not of class 'numeric', converting to numeric")
    sumstats[[SE]] <- as.numeric(sumstats[[SE]])
  }
  colnames(sumstats)[colnames(sumstats) == SE] <- SE_new
  if (class(sumstats[[nlog10p_value]]) != "numeric") {
    stop(paste0("The nlog10p column is not 'numeric' class, please check sumstats."))
  }
  colnames(sumstats)[colnames(sumstats) == nlog10p_value] <- nlog10p_value_new
  ### Optional columns
  if (!is.null(rsID)) {
    colnames(sumstats)[colnames(sumstats) == rsID] <- rsID_new
  }
  if (!is.null(CHR)) {
    colnames(sumstats)[colnames(sumstats) == CHR] <- CHR_new
  }
  if (!is.null(POS)) {
    if (class(sumstats[[POS]]) != "integer") {
      warning("POS column is not of the 'integer' class, converting to integer")
      sumstats[[POS]] <- as.integer(sumstats[[POS]])
    }
    colnames(sumstats)[colnames(sumstats) == POS] <- POS_new
  }
  if (!is.null(A2)) {
    colnames(sumstats)[colnames(sumstats) == A2] <- A2_new
  }
  if (!is.null(AF)) {
    colnames(sumstats)[colnames(sumstats) == AF] <- AF_new
  }
  if (!is.null(N)) {
    colnames(sumstats)[colnames(sumstats) == N] <- N_new
  }
  return(sumstats)
}

#' Process p-values to handle underflow
#' under development
#' @export
handle_underflow <- function(pvalue_vec,
                             return_nlog10P=F) {
  if (!"Rmpfr" %in% rownames(installed.packages())) {
    stop("Rmpfr is required to run this function")
  }
  message("Converting to numeric using Rmpfr::mpfr() ... ")
  pvalue_vec <- Rmpfr::mpfr(pvalue_vec)
  if (return_nlog10P) {
    message("Converting p-values to negative log10 scale")
    pvalue_vec <- as.numeric(-log10(pvalue_vec))
  }
  return(pvalue_vec)
}

#' Subset chromosomes with significant signals
#' to speed up computation
#' under development
#' @export
subset_chromosomes <- function(sumstats, CHR_name, nlog10p_value_name, nlogP_threshold) {
  CHR_to_keep <- unique(sumstats[[CHR_name]][sumstats[[nlog10p_value_name]] > nlogP_threshold])
  sumstats <- sumstats[sumstats[[CHR_name]] %in% CHR_to_keep,]
  message(paste0("Found significant regions on chromosomes: ", paste(CHR_to_keep, collapse = ", ")))
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
                              p_value_name = NULL,
                              nlog10p_value_name = "nlog10P",
                              CHR_out = "CHR_var",
                              BP_START_var_out = "BP_START_var",
                              BP_STOP_var_out = "BP_STOP_var",
                              nlogP_threshold = 5e-8,
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
    regions_log <- c(regions_log, paste0("Solving region ", region_var))
    message(paste0("Solving region ", region_var))
    # Rmpfr has potential bug with which.min, therefore a fix
    which_max <- which(sumstats[[nlog10p_value_name]] == max(sumstats[[nlog10p_value_name]]))
    min_p_row <- sumstats[which_max,]
    min_p_row[[nlog10p_value_name]] <- as.numeric(min_p_row[[nlog10p_value_name]])
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
  coloc_regions_PASS <- subset(coloc_regions, comment == "PASS", c("CHR_var", "BP_START_var", "BP_STOP_var"))
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






