#' @title Read sumstats 1
#' @param sumstats_file path to sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @export
query_sumstats_1 <- function(sumstats_file,
                             CHR_var, BP_START_var, BP_STOP_var,
                             ...,
                             read_mode = "tabix") {
  if (read_mode == "tabix") {
    sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file, " ",
                                              CHR_var, ":", BP_START_var, "-",
                                              BP_STOP_var), intern = T), sep = "\t", header = T)
  }
  if (read_mode == "RDS") {
    sumstats <- readRDS(sumstats_file)
  }
  if (read_mode == "read.csv") {
    sumstats <- read.csv(sumstats_file)
  }
  if (read_mode == "get") {
    sumstats <- get(sumstats_file)
  }
  if (nrow(sumstats) == 0) { return(sumstats) }
  sumstats <- subset(sumstats, CHR == CHR_var & POS >= BP_START_var & POS <= BP_STOP_var)
  return(sumstats)
}


#' @title Query Kidney eQTL
#' @description Query GTEx v8 GWAS data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @examples
#' under development
#' @export
query_kidney_eQTL <- function(sumstats_file,
                              CHR_var, BP_START_var, BP_STOP_var,
                              ...) {
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T), header = T)
  if (nrow(sumstats) == 0) { sumstats <- data.frame(); return(sumstats) }
  sumstats <- unique(sumstats)
  # format by phenotype ID
  sumstats_list <- lapply(unique(sumstats$GeneID), function(x) {
    sumstats <- subset(sumstats, GeneID == x)
    sumstats$AF <- NA
    sumstats$N <- NA
    sumstats <- sumstats[,c("Name", "rsID", "CHR", "POS_hg38", "Alt", "Ref", "Beta", "Std", "Pvalue", "AF", "N", "GeneID")]
    colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N", "Phenotype")
    sumstats <- subset(sumstats, (!is.na(BETA)) & (!is.na(SE)))
    sumstats <- subset(sumstats, (! BETA %in% c(Inf, -Inf)) & (! SE %in% c(Inf, -Inf)))
    if (length(unique(sumstats$Phenotype)) > 1) {stop("Phenotype not unique in output query")}
    return(sumstats)
  })
  return(sumstats_list)
}

#' @title query UKB PPP pGWAS
#' @description Query UKB PPP pGWAS data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_UKB_PPP_EUR <- function(sumstats_file,
                              CHR_var, BP_START_var, BP_STOP_var, ...,
                              handle_underflow=F,
                              colClasses_int=c(4L,5L), ncol_sumstats=14) {
  if (CHR_var == "X") {CHR_var <- "23"}
  colClasses <- readtable_colCl(ncol_sumstats, colClasses_int)
  sumstats <- tabix_fun(sumstats_file, CHR_var, BP_START_var, BP_STOP_var, colClasses)
  sumstats$CHROM[sumstats$CHROM == "23"] <- "X"
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$ID <- paste0("chr", sumstats$CHROM, ":", sumstats$GENPOS, ":", sumstats$ALLELE0, ":", sumstats$ALLELE1)
  sumstats$rsID <- NA
  if (handle_underflow) {
    sumstats$P <- 10^(-handle_underflow(sumstats[["LOG10P"]]))
  } else {
    sumstats$P <- 10^(-sumstats$LOG10P)
  }
  sumstats <- sumstats[,c("ID", "rsID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "BETA", "SE", "P", "A1FREQ", "N")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  return(sumstats)
}

#' @title Query finngen GWAS
#' @description Query finngen GWAS data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @examples
#' query_finngen_GWAS(sumstats_file = "finngen_R9_E4_DM2REN.gz", CHR_var = "1", BP_START_var = 100000, BP_STOP_var = 110000)
#' @export
query_finngen_GWAS <- function(sumstats_file,
                               CHR_var, BP_START_var, BP_STOP_var, ...) {
  if (CHR_var == "X") {CHR_var <- "23"}
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T),  header = T, stringsAsFactors = FALSE, 
                         colClasses = c(alt = "character", ref= "character"))
  
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$X.chrom[sumstats$X.chrom == "23"] <- "X"
  sumstats$Name <- paste0("chr", sumstats$X.chrom, ":", sumstats$pos, ":", sumstats$ref, ":", sumstats$alt)
  sumstats <- sumstats[,c("Name", "rsids", "X.chrom", "pos", "alt", "ref", "beta", "sebeta", "pval", "af_alt")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF")
  return(sumstats)
}


### Helper functions ----

#' tabix_fun
#' Helper function to create tabix_cmd
tabix_fun <- function(sumstats_file, CHR_var, BP_START_var, BP_STOP_var,
                      colClasses=NULL, sep = "\t", header = T,
                      show_cmd=F) {
  tabix_cmd <- paste0("tabix -h ", sumstats_file, " ", CHR_var, ":", BP_START_var, "-", BP_STOP_var)
  if (show_cmd) {message(tabix_cmd)}
  tabix_txt <- system(tabix_cmd, intern = T)
  if (!identical(tabix_txt, character(0))) {
    if (!is.null(colClasses)) {
      sumstats <- read.table(text=tabix_txt, sep = sep, header = header, colClasses = colClasses)
    } else {
      sumstats <- read.table(text=tabix_txt, sep = sep, header = header)
    }
  } else {
    sumstats <- data.frame()
  }
  return(sumstats)
}

#' readtable_colCl
#' Helper function to create tabix_cmd
readtable_colCl <- function(ncol_sumstats, colClasses_int) {
  colClasses <- as.character(rep(NA, ncol_sumstats))
  colClasses[colClasses_int] <- "character"
  colClasses
}
