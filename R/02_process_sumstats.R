#' Create data.frame with input parameters for coloc
#' @param sumstats_path path to folder with indexed sumstats.
#' @param sumstats_pattern usually "gz$".
#' @param grep_invert expression to exclude some found files
#' (e.g., annotation files in the same folder).
#' @param sumstats_function function name to query sumstats.
#' Read sumstats 1 and format columns
#' @return data frame with formatted columns
#' @examples
#' under development
#' @export
read_sumstats_1 <- function(sumstats_file,
                            Name_name = NULL,
                            rsID_name = NULL,
                            CHR_name = NULL,
                            POS_name = NULL,
                            A1_name = NULL,
                            A2_name = NULL,
                            BETA_name = NULL,
                            SE_name = NULL,
                            p_value_name = NULL,
                            AF_name = NULL,
                            N_name = NULL,
                            other_columns = NULL,
                            format_columns = F) {
  if (is.null(CHR_name)) { stop("CHR_name is empty") }
  if (is.null(POS_name)) { stop("POS_name is empty") }
  if (is.null(p_value_name)) { stop("p_value_name is empty") }
  if ("data.table" %in% rownames(installed.packages())) {
    if (format_columns == T) {
      sumstats <- data.table::fread(cmd=sumstats_file)
    }
    if (format_columns == F) {
        sumstats <- data.table::fread(sumstats_file)
    }
    sumstats <- as.data.frame(sumstats)
  } else {
    if (format_columns == T) {
      sumstats <- read.delim(text = system(sumstats_file, intern = T), header = T)
    }
    if (format_columns == F) {
      sumstats <- read.delim(sumstats_file, header = T)
    }
  }
  all_cols <- c(Name_name, rsID_name, CHR_name, POS_name,
                A1_name, A2_name, BETA_name, SE_name,
                p_value_name, AF_name, N_name, other_columns)
  sumstats <- sumstats[,all_cols]
  colnames(sumstats)[colnames(sumstats) == CHR_name] <- "CHR"
  colnames(sumstats)[colnames(sumstats) == POS_name] <- "POS"
  colnames(sumstats)[colnames(sumstats) == p_value_name] <- "P"
  if (!is.null(Name_name)) 
    colnames(sumstats)[colnames(sumstats) == Name_name] <- "Name"
  if (!is.null(rsID_name)) 
    colnames(sumstats)[colnames(sumstats) == rsID_name] <- "rsID"
  if (!is.null(A1_name)) 
    colnames(sumstats)[colnames(sumstats) == A1_name] <- "A1"
  if (!is.null(A2_name)) 
    colnames(sumstats)[colnames(sumstats) == A2_name] <- "A2"
  if (!is.null(BETA_name)) 
    colnames(sumstats)[colnames(sumstats) == BETA_name] <- "BETA"
  if (!is.null(SE_name)) 
    colnames(sumstats)[colnames(sumstats) == SE_name] <- "SE"
  if (!is.null(p_value_name)) 
    colnames(sumstats)[colnames(sumstats) == p_value_name] <- "P"
  if (!is.null(AF_name)) 
    colnames(sumstats)[colnames(sumstats) == AF_name] <- "AF"
  if (!is.null(N_name)) 
    colnames(sumstats)[colnames(sumstats) == N_name] <- "N"
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
                              POS_name = "POS",
                              p_value_name = "P",
                              p_threshold = 5e-8,
                              halfwindow = 500000
                              ) {
  if (is.character(sumstats[[p_value_name]])) {
    warning(paste0(p_value_name, " column is character, converting to numeric"))
    sumstats[[p_value_name]] <- as.numeric(sumstats[[p_value_name]])
  }
  # set up variables
  coloc_regions <- data.frame()
  regions_log <- c()
  # function-specific constants
  region_var <- 1
  comment_var <- "PASS"
  # start iterations
  while(min(sumstats[[p_value_name]], na.rm = T) < p_threshold) {
    regions_log <- c(regions_log, paste0("Solving region ", region_var))
    print(paste0("Solving region ", region_var))
    min_p_row <- sumstats[which.min(sumstats[[p_value_name]]),]
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
    # return results
    colnames(coloc_regions)[colnames(coloc_regions) == "CHR"] <- "CHR_var"
    colnames(coloc_regions)[colnames(coloc_regions) == "BP_START"] <- "BP_START_var"
    colnames(coloc_regions)[colnames(coloc_regions) == "BP_STOP"] <- "BP_STOP_var"
    start_cols <- which(colnames(coloc_regions) %in% c("CHR_var", "BP_START_var", "BP_STOP_var"))
    end_cols <- which(!colnames(coloc_regions) %in% c("CHR_var", "BP_START_var", "BP_STOP_var"))
    coloc_regions <- coloc_regions[,c(start_cols, end_cols)]
    rownames(coloc_regions) <- NULL
  }
  return(list(coloc_regions = coloc_regions, regions_log = regions_log))
}

#' subset_sumstats_1
#' @description under development
#' @export
subset_sumstats_1 <- function(CHR_var, BP_START_var, BP_STOP_var,
                              sumstats,
                              CHR_name = "CHR",
                              POS_name = "POS",
                              dbSNP_file = NULL,
                              CHR_var_col_name = NULL,
                              BP_START_var_col_name = NULL,
                              BP_STOP_var_col_name = NULL,
                              remove_duplicates = T,
                              do_match_rs = F, N_mc.cores = 10, ...) {
  sumstats_filt <- subset(sumstats, sumstats[[CHR_name]] == CHR_var & 
                            sumstats[[POS_name]] >= BP_START_var & 
                            sumstats[[POS_name]] <= BP_STOP_var)
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
    # sumstats_filt_rs_unmatched <- subset(sumstats_filt_rs[!matched_alleles,], !SNP %in% sumstats_filt_rs_matched$SNP)
    # sumstats_filt <- sumstats_filt_rs_matched
  }
  if (remove_duplicates) {
    sumstats_filt <- unique(sumstats_filt)
  }
  return(sumstats_filt)
}


#' process sumstats 1
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
  sumstats_1 <- read_sumstats_1(sumstats_file,
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
                                function(i) do.call(subset_sumstats_1, c(coloc_regions[i,], extra_args)),
                                mc.cores = 8)
  sumstats_1_subset <- unique(do.call(rbind, sumstats_1_subset))
  # save
  writeLines(regions_log, con = paste0(sumstats_name_out, "/", basename(sumstats_name_out), "_log.txt"))
  saveRDS(sumstats_1_subset, paste0(sumstats_name_out, "/", basename(sumstats_name_out), "_subset.RDS"))
  saveRDS(coloc_regions, paste0(sumstats_name_out, "/", basename(sumstats_name_out), "_coloc_regions.RDS"))
  return("PASS")
}




#' match_rs
#' @description under development
#' @export
match_rs <- function(dbSNP_file, sumstats,
                     CHR_var_col_name = "CHR", BP_START_var_col_name = "POS",
                     BP_STOP_var_col_name = "POS",
                     A1_name = "A1", A2_name = "A2", SNP_name = "rsID", ...) {
  # rs_df
  sumstats$comment <- NA
  rsID <- sumstats[[SNP_name]]
  # test if several rsIDs are present
  rsID_split <- strsplit(rsID, ";")[[1]]
  if (length(rsID_split) > 1) {
    rs_df_list <- lapply(rsID_split, function(rsID) {
      CHR_var <- sumstats[[CHR_var_col_name]]
      BP_START_var <- sumstats[[BP_START_var_col_name]]
      BP_STOP_var <- sumstats[[BP_STOP_var_col_name]]
      rs_df <- query_dbSNP(dbSNP_file, CHR_var, BP_START_var, BP_STOP_var, rsID)
      return(rs_df)
    })
    rs_df <- do.call(rbind, rs_df_list)
  } else {
    CHR_var <- sumstats[[CHR_var_col_name]]
    BP_START_var <- sumstats[[BP_START_var_col_name]]
    BP_STOP_var <- sumstats[[BP_STOP_var_col_name]]
    rs_df <- query_dbSNP(dbSNP_file, CHR_var, BP_START_var, BP_STOP_var, rsID)
  }
  if (length(rs_df[[SNP_name]]) == 1) {
    if (is.na(rs_df[[SNP_name]])) {
      sumstats[[SNP_name]] <- NA
      sumstats$Name <- NA
      return(sumstats)
    }
  }
  rs_df <- subset(rs_df, rsID %in% sumstats[[SNP_name]])
  rs_df <- unique(rs_df)
  # sumstats
  sumstats_before_merge <- sumstats
  sumstats <- merge(sumstats, by.x=SNP_name, all.x=T, 
                    rs_df, by.y="rsID")
  sumstats$comment <- "not_inverted"
  matched_alleles <- ((sumstats$REF == sumstats[[A1_name]]) & (sumstats$ALT == sumstats[[A2_name]])) |
    ((sumstats$REF == sumstats[[A2_name]]) & (sumstats$ALT == sumstats[[A1_name]]))
  sumstats_matched <- sumstats[matched_alleles,]
  if (nrow(sumstats_matched) == 0) {
    sumstats[[paste0(A1_name, "_flip")]] <- flip_alleles(sumstats[[A1_name]])
    sumstats[[paste0(A2_name, "_flip")]] <- flip_alleles(sumstats[[A2_name]])
    matched_alleles <- ((sumstats$REF == sumstats[[paste0(A1_name, "_flip")]]) & (sumstats$ALT == sumstats[[paste0(A2_name, "_flip")]])) |
      ((sumstats$REF == sumstats[[paste0(A2_name, "_flip")]]) & (sumstats$ALT == sumstats[[paste0(A1_name, "_flip")]]))
    sumstats[[paste0(A1_name, "_flip")]] <- NULL
    sumstats[[paste0(A2_name, "_flip")]] <- NULL
    sumstats_matched <- sumstats[matched_alleles,]
    if(nrow(sumstats_matched) > 0) {
      sumstats_matched$comment <- "inverted"
    } else {
      sumstats$comment <- "alleles_not_matched"
      sumstats_matched <- sumstats_before_merge
      sumstats_matched$Name <- NA
    }
  }
  sumstats_matched$REF <- NULL
  sumstats_matched$ALT <- NULL
  return(sumstats_matched)
}

