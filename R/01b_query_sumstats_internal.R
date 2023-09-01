#' @title Query GCKD pGWAS
#' @description Internal function, under development
#' @param sumstats_file path to GCKD Olink sumstats.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_GCKD_pGWAS <- function(sumstats_file, CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T), header = T)
  if (nrow(sumstats) == 0) { return(NA) }
  # format
  sumstats$rsID <- NA
  sumstats <- sumstats[,c("SNP", "rsID", "chr", "position", "noncoded_all", "coded_all", "beta", "SE", "pval", "AF_coded_all", "n_total")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  sumstats$BETA <- -sumstats$BETA
  sumstats$AF <- 1 - sumstats$AF
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  # output
  return(sumstats)
}

#' @title Query GCKD mGWAS plasma or urine
#' @description Internal function, under development
#' @param sumstats_file path to sumstats.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_GCKD_mGWAS <- function(sumstats_file,
                             CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T), header = T)
  if (nrow(sumstats) == 0) { return(NA) }
  # format
  sumstats$rsID <- NA
  sumstats$P <- 10^(-sumstats$LOG10P)
  sumstats <- sumstats[,c("ID", "rsID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "BETA", "SE", "P", "A1FREQ", "N")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  return(sumstats)
}

#' @title query TAS
#' @description Internal function, under development
#' @param sumstats_file path to sumstats.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_TAS <- function(sumstats_file,
                      CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var),
                                     intern = T), header = T)
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$rsID <- sumstats$ID
  sumstats$Name <- paste0("chr", sumstats$CHR, ":", sumstats$BP, ":", sumstats$Other_allele, ":", sumstats$Effect_Allele)
  sumstats$AF <- NA
  sumstats$N <- NA
  sumstats <- sumstats[,c("Name", "ID", "CHR", "BP", "Effect_Allele", "Other_allele", "Beta", "SE", "P", "AF", "N")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  # checks
  # stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
  #                 all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  return(sumstats)
}


#' @title query CAD
#' @description Internal function, under development
#' @param sumstats_file path to sumstats.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_CAD <- function(sumstats_file,
                      CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var),
                                     intern = T), header = T)
  if (nrow(sumstats) == 0) { return(sumstats) }
  sumstats$rsID <- sumstats$ID
  sumstats$Name <- paste0("chr", sumstats$CHR, ":", sumstats$BP, ":", sumstats$Other_allele, ":", sumstats$Effect_allele)
  # format
  sumstats$AF <- NA
  sumstats$N <- NA
  sumstats <- sumstats[,c("Name", "ID", "CHR", "BP", "Effect_allele", "Other_allele", "Beta", "SE", "P.value", "AF", "N")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  # checks
  # stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
  #                 all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  return(sumstats)
}


