# TODO
# p_threshold p-value to include sumstats in coloc, default 1e-6
# maf allele frequency to filter out SNPs, default 0.01
# checks: first allele is reference
# checks: duplicated SNPs
# checks: A1 is alternative
# coloc: add sdY

#' @title Query GCKD pGWAS
#' @description Query GCKD pGWAS data to extract a region of interest
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
#' @description Query GCKD mGWAS plasma or urine data to extract a region of interest
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

#' @title Query UKB GWAS
#' @description Query UK Biobank GWAS data to extract a region of interest
#' @param sumstats_file path to UKB sumstats.
#' @param annotation_file annotation file with number of samples for each GWAS.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @export
query_UKB_GWAS <- function(sumstats_file,
                           CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T), sep = "\t", header = T)
  if (nrow(sumstats) == 0) { return(NA) }
  # format
  sumstats$Name <- paste0("chr", sumstats$chrom, ":", sumstats$pos, ":", sumstats$ref, ":", sumstats$alt)
  sumstats <- sumstats[,c("Name", "rsids", "chrom", "pos", "alt", "ref", "beta", "sebeta", "pval", "af")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF")
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  return(sumstats)
}

#' @title query Ferkingstad pGWAS
#' @description Query Ferkingstad pGWAS (FpG) data to extract
#' every pQTL region according to GCKD
#'
#' @param sumstats_file path to sumstats sumstats.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @param assocvariants_annotate_file FpG file with corrected maf and other QC
#' @return data frame with extracted sumstats
#' @export
query_Ferkingstad_pGWAS <- function(sumstats_file,
                                    CHR_var, BP_START_var, BP_STOP_var,
                                    assocvariants_annotate_file) {
  sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file,
                                            " chr", CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T), header = T)
  if (nrow(sumstats) == 0) { return(NA) }
  # FpG specific data processing
  assocvariants_annotated <- read.table(text=system(paste0("tabix -h ", assocvariants_annotate_file,
                                                           " chr", CHR_var, ":",
                                                           BP_START_var, "-", BP_STOP_var), intern = T), header = T)
  sumstats <- merge(sumstats, assocvariants_annotated[,c("Name","effectAlleleFreq")], by="Name")
  # format
  sumstats$Chrom <- gsub("chr", "", sumstats$Chrom)
  sumstats <- sumstats[,c("Name", "rsids", "Chrom", "Pos", "effectAllele", "otherAllele", "Beta", "SE", "Pval", "effectAlleleFreq", "N")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  sumstats$Name <- gsub("(.*:.*):(.*):(.*)", "\\1:\\3:\\2", sumstats$Name)
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  # output
  return(sumstats)
}


#' @title query UKB PPP pGWAS
#' @description Query UKB PPP pGWAS data to extract a region of interest
#' @param sumstats_file path to sumstats.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_UKB_PPP_EUR <- function(sumstats_file,
                              CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T), header = T)
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$ID <- paste0("chr", sumstats$CHROM, ":", sumstats$GENPOS, ":", sumstats$ALLELE0, ":", sumstats$ALLELE1)
  sumstats$rsID <- NA
  sumstats$P <- 10^(-sumstats$LOG10P)
  sumstats <- sumstats[,c("ID", "rsID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "BETA", "SE", "P", "A1FREQ", "N")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  return(sumstats)
}


#' @title query ARIC pGWAS
#' @description Query ARIC pGWAS data to extract a region of interest
#' @param sumstats_file path to sumstats.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_ARIC_pGWAS <- function(sumstats_file,
                             CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var, " | sed 's/\\#CHROM/CHROM/g'"),
                                     intern = T), header = T)
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$rsID <- sumstats$ID
  indeces_to_change <- which(!sumstats$ALT == sumstats$A1)
  REF_to_change <- sumstats[indeces_to_change,]$REF
  ALT_to_change <- sumstats[indeces_to_change,]$ALT
  sumstats[indeces_to_change,]$REF <- ALT_to_change
  sumstats[indeces_to_change,]$ALT <- REF_to_change
  sumstats[indeces_to_change,]$BETA <- -sumstats[indeces_to_change,]$BETA
  sumstats[indeces_to_change,]$A1_FREQ <- 1-sumstats[indeces_to_change,]$A1_FREQ
  sumstats$ID <- paste0("chr", sumstats$CHROM, ":", sumstats$POS, ":", sumstats$REF, ":", sumstats$ALT)
  sumstats <- sumstats[,c("ID", "rsID", "CHROM", "POS", "ALT", "REF", "BETA", "SE", "P", "A1_FREQ", "OBS_CT")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  # output
  return(sumstats)
}

#' @title query CKD Gen
#' @description Query CKD Gen data to extract a region of interest
#' @param sumstats_file path to sumstats.
#' @param CHR_var CHR (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_CKD_pGWAS <- function(sumstats_file,
                            CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var),
                                     intern = T), header = T)
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$rsID <- sumstats$variant_id
  sumstats$Name <- paste0("chr", gsub("_", ":", sumstats$hm_variant_id))
  sumstats <- sumstats[,c("Name", "rsID", "hm_chrom", "hm_pos", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "p_value", "hm_effect_allele_frequency", "n")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  # output
  return(sumstats)
}

#' @title query TAS
#' @description Query Tricuspid Aortic Stenosis data to extract a region of interest
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
#' @description Query Coronary Artery Disease data to extract a region of interest
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


