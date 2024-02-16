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
  if (is.numeric(pvalue_vec)) {
    pvalue_vec <- Rmpfr::mpfr(pvalue_vec, precBits=100)
  } else {
    pvalue_vec <- Rmpfr::mpfr(pvalue_vec)
  }
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
  # sort
  sumstats_filt <- sumstats_filt[with(sumstats_filt, order(sumstats_filt[[CHR_name]], sumstats_filt[[POS_name]])), ]
  # return
  return(list(coloc_regions = coloc_regions,
              coloc_regions_PASS = coloc_regions_PASS,
              regions_log = regions_log,
              sumstats_filt = sumstats_filt))
}

#' Save coloc regions
#' @export
save_coloc_regions <- function(coloc_regions_list, sumstats_name, max_row=100000,
                               SKIP_name="Name", CHR_place=3, POS_place=4,
                               bgzip_bin="bgzip", tabix_bin="tabix",
                               regions_log="regions_log", coloc_regions="coloc_regions",
                               coloc_regions_PASS="coloc_regions_PASS",
                               sumstats_filt="sumstats_filt") {
  writeLines(coloc_regions_list[[regions_log]], con = paste0(sumstats_name, "_log.txt"))
  write.table(coloc_regions_list[[coloc_regions]], paste0(sumstats_name, "_", coloc_regions, ".tsv"),
              sep="\t", row.names = F, col.names = T, quote = F)
  write.table(coloc_regions_list[[coloc_regions_PASS]], paste0(sumstats_name, "_", coloc_regions_PASS, ".tsv"),
              sep="\t", row.names = F, col.names = T, quote = F)
  # if (nrow(coloc_regions_list[[sumstats_filt]]) > max_row) {
  write.table(coloc_regions_list[[sumstats_filt]], paste0(sumstats_name, "_subset.tsv"),
              sep="\t", row.names = F, col.names = T, quote = F)
  system(paste0(bgzip_bin, " -f ", sumstats_name, "_subset.tsv"))
  system(paste0("tabix -f -s", CHR_place, " -b", POS_place, " -e", POS_place, " ", sumstats_name, "_subset.tsv.gz -c ", SKIP_name))
  # } else {
  #   saveRDS(coloc_regions_list[[sumstats_filt]], paste0(sumstats_name, "_subset.RDS"))
  # }
}

#' genepi liftOver
#' @export
genepi_liftOver <- function(sumstats, CHR_name, POS_name, A1_name, A2_name,
                            liftOver_bin, liftOver_chain_hg19ToHg38, dbSNP_file,
                            unique_ID_name="unique_ID",
                            mc_cores=4, keep_lower=F, do_soring=T, rm_tmp_liftOver=T) {
  if (missing(liftOver_bin)) {
    stop(paste0("Please provide path to liftOver_bin. ",
                "See https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"))
  }
  if (missing(liftOver_chain_hg19ToHg38)) {
    stop(paste0("Please provide path to liftOver_chain_hg19ToHg38 ",
                "For 'hg19ToHg38' see https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"))
  }
  if (missing(liftOver_chain_hg19ToHg38)) {stop("Please provide path to liftOver_chain_hg19ToHg38")}
  if (missing(dbSNP_file)) {stop("Please provide path to dbSNP_file")}
  if (!"data.table" %in% rownames(installed.packages())) {
    stop("data.table is required to run this function")
  }
  if (!"parallel" %in% rownames(installed.packages())) {
    stop("parallel is required to run this function")
  }
  library(data.table)
  library(parallel)
  if (! "data.table" %in% class(sumstats)) {
    message("sumstats class is data.frame, while data.table is required. Chaning to data.table")
    data.table::setDT(sumstats)
  }
  # set tmp name - sumstats will be saved to the disk for the liftOver, then it will be deleted.
  tmp_name <- paste(sample(letters, 20, replace = T), collapse = "")
  # add index column
  sumstats[[unique_ID_name]] <- paste0("row_", rownames(sumstats))
  # check and add chr prefix, if needed
  if (! any(grepl("chr", sumstats[[CHR_name]]))) {
    message("No 'chr' prefix in the chromosome column found (e.g., '21' should be 'chr21' as required by UCSC liftOver), adding 'chr' prefix. This is a regular message, not a warning.")
    sumstats[[CHR_name]] <- paste0("chr", sumstats[[CHR_name]])
    change_chr_back <- T
  } else { 
    change_chr_back <- F
  }
  # set alleles to upper case
  if ((any(grepl("[[:lower:]]", sumstats[[A1_name]])) | any(grepl("[[:lower:]]", sumstats[[A2_name]]))) & (!keep_lower)) {
    message("lower case alleles coding detected, changing automatically to upper case. Turn off this behaviour using keep_lower=T option, if needed.")
    sumstats[[A1_name]] <- toupper(sumstats[[A1_name]])
    sumstats[[A2_name]] <- toupper(sumstats[[A2_name]])
  }
  # Add necessary columns for liftOver
  sumstats[["score"]] <- "."
  sumstats[["strand"]] <- "+"
  data.table::fwrite(sumstats[,c(CHR_name, POS_name, POS_name, unique_ID_name, "score", "strand"), with=F],
                     paste0(tmp_name, "_liftOver.bed"),
                     sep="\t", quote=F, col.names = F, row.names = F)
  message(paste0("Temporary liftOver files will be written to the disk in the current working directory. ",
                 "After the funcion is finished, these files will be automatically removed. ",
                 "If any issues are encountered, ensure that there is enough space on the disk. ",
                 "Usually sumstats file will have a size of 300-600 Mb"))
  system(paste(liftOver_bin,
               paste0(tmp_name, "_liftOver.bed"),
               liftOver_chain_hg19ToHg38,
               paste0(tmp_name, "_liftOver_newFile.bed"),
               paste0(tmp_name, "_liftOver_unMapped.bed")))
  newFile <- data.table::fread(paste0(tmp_name, "_liftOver_newFile.bed"), header = F)
  unMapped <- tryCatch(read.table(paste0(tmp_name, "_liftOver_unMapped.bed"), header = F), error=function(e) NULL)
  # clean tmp liftOver files
  if (rm_tmp_liftOver) {
    system(paste("rm ", paste0(tmp_name, "_liftOver.bed"),
                 paste0(tmp_name, "_liftOver_newFile.bed"),
                 paste0(tmp_name, "_liftOver_unMapped.bed")))
  }
  if (!is.null(unMapped)) {
    message(paste0("There were ", nrow(unMapped), " unmapped variants (out of ",
                   nrow(sumstats), " variants), they will be discarded. ",
                   "Usually this should be 0.01-0.1%. ",
                   "If it is substantially higher, please rerun with rm_tmp_liftOver=F and check the issue"))
    # Discard unMapped variants
    sumstats <- sumstats[!sumstats[[unique_ID_name]] %in% unMapped[["V4"]],]
  }
  # merge with newFile
  if ((nrow(sumstats) != nrow(newFile)) | any(newFile[["V2"]] != newFile[["V3"]])) {
    warning("Number of rows in sumstats does not match with liftOver results, please check the cause.")
  }
  # remove unused columns and update colnames
  newFile[["V3"]] <- newFile[["V5"]] <- NULL
  colnames(newFile) <- c("CHR_hg38", "POS_hg38", unique_ID_name, "strand_hg38")
  message("Merging with new coordinates and performing flip...")
  sumstats_liftOver <- merge(sumstats, newFile, by=unique_ID_name, sort = F)
  if (nrow(sumstats_liftOver) != nrow(newFile)) {
    warning("Number of rows before and after merge does not match, please check the cause.")
  }
  # flip if necessary
  indeces_to_flip <- sumstats_liftOver$strand_hg38 == "-"
  message(paste0("Flipping strand for ", sum(indeces_to_flip), " variants."))
  sumstats_liftOver[["A1_hg38"]] <- ifelse(indeces_to_flip,
                                           flip_alleles(sumstats_liftOver[[A1_name]]),
                                           sumstats_liftOver[[A1_name]])
  sumstats_liftOver[["A2_hg38"]] <- ifelse(indeces_to_flip,
                                           flip_alleles(sumstats_liftOver[[A2_name]]),
                                           sumstats_liftOver[[A2_name]])
  # Creating new Name column
  ### Major function to append Name with REF:ALT matched by position
  sumstats_liftOver <- Name_by_position(sumstats=sumstats_liftOver,
                                        tmp_name=tmp_name,
                                        dbSNP_file=dbSNP_file)
  # remove "chr"
  if (change_chr_back) {
    sumstats_liftOver[[CHR_name]] <- gsub("chr", "", sumstats_liftOver[[CHR_name]])
  }
  # remove liftOver columns
  for (i in c("score", "strand", "strand_hg38", unique_ID_name)) {
    sumstats_liftOver[[i]] <- NULL }
  sumstats_liftOver[["CHR_hg38"]] <- gsub("chr", "", sumstats_liftOver[["CHR_hg38"]])
  return(sumstats_liftOver)
}


#' Find name by position
#' @export
Name_by_position <- function(sumstats, tmp_name=NULL,
                             CHR_name="CHR_hg38", POS_name="POS_hg38",
                             A1_name="A1_hg38", A2_name="A2_hg38",
                             Name_out="Name_hg38", rs_name="rs",
                             unique_ID_name="unique_ID",
                             tabix_bin, dbSNP_file,
                             do_soring=T, mc_cores=4) {
  if (missing(tabix_bin)) {stop("Please provide path to tabix_bin")}
  if (missing(dbSNP_file)) {stop("Please provide path to dbSNP_file")}
  if (!"data.table" %in% rownames(installed.packages())) {
    stop("data.table is required to run this function")
  }
  if (!"parallel" %in% rownames(installed.packages())) {
    stop("data.table is required to run this function")
  }
  library(data.table); library(parallel)
  if (is.null(tmp_name)) {
    tmp_name <- paste(sample(letters, 20, replace = T), collapse = "")
  }
  # check and add chr prefix, if needed
  if (! any(grepl("chr", sumstats[[CHR_name]]))) {
    message(paste0("No 'chr' prefix in the chromosome column found ",
                   " (e.g., '21' should be 'chr21' as required by UCSC liftOver)",
                   " Adding now 'chr' prefix and after the function is finished,",
                   " chromosomes will be renamed back ('chr' prefix will be removed)."))
    sumstats[[CHR_name]] <- paste0("chr", sumstats[[CHR_name]])
    change_chr_back <- T
  } else { 
    change_chr_back <- F
  }
  message(paste0("Starting name matching by position. ",
                 "Different sumstats can have arbitrary allele order. ",
                 "Especially if they are coming from a metal meta-analysis. ",
                 "This function will try to ensure that REF allele is REF, and ALT is ALT (to harmonize with other sumstats). ",
                 "This procedue is not always possible for ambiguous SNPs or indels.",
                 "If non-standard CHR are present (e.g., chr7_KI270803v1_alt), they will be ignored."))
  CHR_vec <- unique(sumstats[[CHR_name]])
  # extract tabix
  out_list <- mclapply(CHR_vec, function(CHR_var) {
    sumstats_chr <- sumstats[sumstats[[CHR_name]] == CHR_var]
    # add index column
    sumstats_chr[[unique_ID_name]] <- paste0("row_", rownames(sumstats_chr))
    sumstats_chr <- sumstats_chr[order(sumstats_chr[[POS_name]])]
    Name_match <- paste0(tmp_name, "_", CHR_var, "_tabix")
    data.table::fwrite(sumstats_chr[,c(CHR_name, POS_name), with=F],
                       Name_match, sep="\t", col.names = F)
    system(paste0(tabix_bin, " -h ", dbSNP_file, " -R ", Name_match, " -cache 5000 > ", Name_match, "_out"))
    dbSNP_subset <- suppressWarnings(data.table::fread(paste0(Name_match, "_out")))
    stopifnot(all(colnames(dbSNP_subset) == paste0("V", 1:6)))
    # cleaning
    system(paste0("rm ", Name_match, " ", Name_match, "_out"))
    if (nrow(dbSNP_subset) == 0) {
      vec_out <- NA
      names(vec_out) <- CHR_var
      message(paste0("tabix for CHR '", CHR_var, "' returned 0 lines (that is ok for nonstandard chromosomes)."))
      return(c(vec_out))
    }
    # remove extra lines from indels (e.g., indel with 10 BP starting at 5000 will match to POS 5000:5009) 
    dbSNP_subset <- dbSNP_subset[dbSNP_subset[["V2"]] %in% sumstats_chr[[POS_name]],]
    # second allele is comma-separated, replace by 1 allele at each line
    dbSNP_subset[["V5"]] <- gsub(".*:(.*)", "\\1", dbSNP_subset[["V6"]])
    dbSNP_subset[["V1"]] <- NULL
    colnames(dbSNP_subset)[colnames(dbSNP_subset) == "V3"] <- rs_name
    colnames(dbSNP_subset)[colnames(dbSNP_subset) == "V6"] <- Name_out
    sumstats_chr <- merge(sumstats_chr, by.x=POS_name,
                          dbSNP_subset, by.y="V2", all.x=T,
                          allow.cartesian=TRUE)
    index_found <- (sumstats_chr[["V4"]] == sumstats_chr[[A1_name]] & sumstats_chr[["V5"]] == sumstats_chr[[A2_name]]) |
      (sumstats_chr[["V5"]] == sumstats_chr[[A1_name]] & sumstats_chr[["V4"]] == sumstats_chr[[A2_name]])
    sumstats_chr <- sumstats_chr[index_found,]
    sumstats_chr <- sumstats_chr[!duplicated(sumstats_chr[[unique_ID_name]]),]
    sumstats_chr[["V4"]] <- sumstats_chr[["V5"]] <- NULL
    # sort for indexing
    if (do_soring) {
      sumstats_chr <- sumstats_chr[order(sumstats_chr[[POS_name]])]
    }
    message(paste0("Name annotation for CHR ", CHR_var, " is done."))
    return(sumstats_chr)
  }, mc.cores = mc_cores)
  discarded_names <- names(unlist(out_list[!sapply(out_list, is.data.frame)]))
  if (is.null(discarded_names)) {
    message("Tabix successfully processed all chromosomes.")
  } else {
    message(paste0("Tabix did not found any matching variants for the following chromosomes: ",
                   paste(discarded_names, collapse = ", "), ".\n",
                   "These chromosomes will be discarded"))
    out_list <- out_list[sapply(out_list, is.data.frame)]
  }
  # rbind and get stats
  out_dt <- data.table::rbindlist(out_list)
  In_situ_matching_1 <- paste0(out_dt[[CHR_name]], ":",
                               out_dt[[POS_name]], ":",
                               out_dt[[A2_name]], ":",
                               out_dt[[A1_name]])
  In_situ_matching_2 <- paste0(out_dt[[CHR_name]], ":",
                               out_dt[[POS_name]], ":",
                               out_dt[[A1_name]], ":",
                               out_dt[[A2_name]])
  message(paste0("Here are some matching statistics:\n",
                 "Assuming that A2 was REF: ",
                 sum(out_dt[[Name_out]] == In_situ_matching_1),
                 " out of ", length(In_situ_matching_1), " could have been matched.\n",
                 "Assuming that A1 was REF: ",
                 sum(out_dt[[Name_out]] == In_situ_matching_2),
                 " out of ", length(In_situ_matching_2), " could have been matched.\n",
                 "If one of the options gives a good yield (>80-90%), then one can assume that ",
                 "A1/A2 corresponded indeed REF/ALT (or inversely). If the yield is around 50%, then it is ",
                 "necessary to use names matched-by-position which is already appended as a new column."))
  # remove "chr"
  if (change_chr_back) {
    out_dt[[CHR_name]] <- gsub("chr", "", out_dt[[CHR_name]])
  }
  return(out_dt)
}


#' match rs
match_rs <- function(dbSNP_file, sumstats,
                     CHR_var_col_name = "CHR", BP_START_var_col_name = "POS",
                     BP_STOP_var_col_name = "POS",
                     A1_name = "A1", A2_name = "A2", SNP_name = "rsID",
                     window_size = 100, ...) {
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
    BP_START_var <- sumstats[[BP_START_var_col_name]]-window_size
    BP_STOP_var <- sumstats[[BP_STOP_var_col_name]]+window_size
    rs_df <- query_dbSNP(dbSNP_file, CHR_var, BP_START_var, BP_STOP_var, rsID)
  }
  if (length(rs_df[[SNP_name]]) == 1) {
    if (is.na(rs_df[[SNP_name]])) {
      sumstats[[SNP_name]] <- NA
      sumstats[["Name_rs_matching"]] <- NA
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
      sumstats_matched[["Name_rs_matching"]] <- NA
    }
  }
  sumstats_matched$REF <- NULL
  sumstats_matched$ALT <- NULL
  return(sumstats_matched)
}
