# TODO
# - Add check for reference allele
# - Split between mandatory and optional columns
#   - mandatory columns for binary traits: ID, BETA, SE 
#   - mandatory columns for quantitative traits: ID, BETA, SE 
# - Add input option for sdY
# - Add to the file with internal examples
# - Add colClasses, e.g., `colClasses = c(alt = "character", ref= "character")`

#' @title Query finngen GWAS
#' @description Query finngen GWAS data to extract a region of interest
#' @param sumstats_file path to finngen sumstats.
#' @param annotation_file annotation file with number of samples for each GWAS.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @export


query_finngen_GWAS <- function(sumstats_file,
                           CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T),  header = T, stringsAsFactors = FALSE, 
                         colClasses = c(alt = "character", ref= "character"))
  
  if (nrow(sumstats) == 0) { return(NA) }
  # format
  sumstats$Name <- paste0("chr", sumstats$X.chrom, ":", sumstats$pos, ":", sumstats$ref, ":", sumstats$alt)
  sumstats <- sumstats[,c("Name", "rsids", "X.chrom", "pos", "alt", "ref", "beta", "sebeta", "pval", "af_alt")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF")
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  return(sumstats)
}


#' @title Query GTEx v7 GWAS
#' @description Query GTEx v7 GWAS data to extract a region of interest
#' @param sumstats_file path to GTEx v7 sumstats.
#' @param annotation_file annotation file with number of samples for each GWAS.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @export

query_GTEXv7_GWAS <- function(sumstats_file,
                               CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T),  header = T,  stringsAsFactors = FALSE, 
                         colClasses = c(alt = "character", ref= "character"))
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T),  header = T,  stringsAsFactors = T)
  
  if (nrow(sumstats) == 0) { return(NA) }
  # format
  sumstats$Name <- paste0("chr", sumstats$chr, ":", sumstats$pos, ":", sumstats$ref, ":", sumstats$alt)
  sumstats$rsids<- NA
  sumstats <- sumstats[,c("Name", "rsids", "chr", "pos", "alt", "ref", "slope", "slope_se", "pval_nominal", "maf")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF")
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  return(sumstats)
}


#' @title Query GTEx v8 GWAS
#' @description Query GTEx v8 GWAS data to extract a region of interest
#' @param sumstats_file path to GTEx v7 sumstats.
#' @param annotation_file annotation file with number of samples for each GWAS.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @export

query_GTEXv8_GWAS <- function(sumstats_file,
                               CHR_var, BP_START_var, BP_STOP_var,
                               phenotype_id_var) {
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " chr",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T),  header = F,  stringsAsFactors = FALSE, 
                         colClasses = c(V5 = "character", V6= "character"))
  
  if (nrow(sumstats) == 0) { return(NA) }
  # phenotype ID
  sumstats <- subset(sumstats, V8 == phenotype_id_var)
  # format
  sumstats$Name <- paste0(sumstats$V1, ":", sumstats$V4, ":", sumstats$V5, ":", sumstats$V6)
  sumstats$rsids<- NA
  sumstats <- sumstats[,c("Name", "rsids", "V19", "V4", "V6", "V5", "V15", "V16", "V14", "V11", "V8")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "Phenotype")
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  return(sumstats)
}
