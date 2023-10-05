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
                             read_mode = "RDS") {
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

#' @title Query UKB GWAS
#' @description Query UK Biobank GWAS data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' @export
query_UKB_GWAS <- function(sumstats_file,
                           CHR_var, BP_START_var, BP_STOP_var, ...) {
  sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T), sep = "\t", header = T)
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$Name <- paste0("chr", sumstats$chrom, ":", sumstats$pos, ":", sumstats$ref, ":", sumstats$alt)
  sumstats <- sumstats[,c("Name", "rsids", "chrom", "pos", "alt", "ref", "beta", "sebeta", "pval", "af")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF")
  return(sumstats)
}

#' @title query Ferkingstad pGWAS
#' @description Query Ferkingstad pGWAS (FpG) data
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @param assocvariants_annotate_file FpG file with corrected maf and other QC
#' @return data frame with extracted sumstats
#' @export
query_Ferkingstad_pGWAS <- function(sumstats_file,
                                    CHR_var, BP_START_var, BP_STOP_var,
                                    assocvariants_annotate_file, ...) {
  sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file,
                                            " chr", CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T), header = T)
  if (nrow(sumstats) == 0) { return(sumstats) }
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
  # output
  return(sumstats)
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
                              CHR_var, BP_START_var, BP_STOP_var, ...) {
  if (CHR_var == "X") {CHR_var <- "23"}
  sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T), header = T)
  sumstats$CHROM[sumstats$CHROM == "23"] <- "X"
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$ID <- paste0("chr", sumstats$CHROM, ":", sumstats$GENPOS, ":", sumstats$ALLELE0, ":", sumstats$ALLELE1)
  sumstats$rsID <- NA
  sumstats$P <- 10^(-sumstats$LOG10P)
  sumstats <- sumstats[,c("ID", "rsID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "BETA", "SE", "P", "A1FREQ", "N")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  return(sumstats)
}


#' @title query ARIC pGWAS
#' @description Query ARIC pGWAS data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_ARIC_pGWAS <- function(sumstats_file,
                             CHR_var, BP_START_var, BP_STOP_var, ...) {
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
  # output
  return(sumstats)
}

#' @title query CKD Gen
#' @description Query CKD Gen data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats
#' @export
query_CKD_pGWAS <- function(sumstats_file,
                            CHR_var, BP_START_var, BP_STOP_var, ...) {
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
  # output
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


#' @title Query GTEx v7 GWAS
#' @description Query GTEx v7 GWAS data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' #' @examples
#' query_GTEXv7_GWAS(sumstats_file = "Whole_Blood.allpairs_tmp_CHR_BP_sorted.txt.gz", CHR_var = "1", BP_START_var = 1000000, 
#' BP_STOP_var = 1010000, phenotype_id_var = "ENSG00000177757.1")
#' @export

query_GTEXv7_GWAS <- function(sumstats_file,
                              CHR_var, BP_START_var, BP_STOP_var,
                              phenotype_id_var, ...) {
  sumstats <- read.delim(text=system(paste0("tabix -h ", sumstats_file, " ",
                                            CHR_var, ":", BP_START_var, "-",
                                            BP_STOP_var), intern = T),  header = T,  stringsAsFactors = FALSE, 
                         colClasses = c(alt = "character", ref= "character"))
  
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format
  sumstats$Name <- paste0("chr", sumstats$chr, ":", sumstats$pos, ":", sumstats$ref, ":", sumstats$alt)
  sumstats$rsids<- NA
  sumstats <- sumstats[,c("Name", "rsids", "chr", "pos", "alt", "ref", "slope", "slope_se", "pval_nominal", 
                          "maf", "ma_samples", "gene_id")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N", "Phenotype")
  return(sumstats)
}


#' @title Query GTEx v8 GWAS
#' @description Query GTEx v8 GWAS data to extract a region of interest
#' @param sumstats_file path tabix-indexed sumstats.
#' @param CHR_var chromosome (as.character "1", "2", ..., "X").
#' @param BP_START_var start of region, integer
#' @param BP_STOP_var end of region, integer
#' @return data frame with extracted sumstats.
#' #' #' @examples
#' query_GTEXv8_GWAS(sumstats_file = "GTEx_V8.gz", CHR_var = "X", BP_START_var = 20000, 
#' BP_STOP_var = 28000, phenotype_id_var = "ENSG00000167393.17")
#' @export
query_GTEXv8_GWAS <- function(sumstats_file,
                              CHR_var, BP_START_var, BP_STOP_var,
                              ...) {
  text_out <- system(paste0("tabix -h ", sumstats_file, " chr",
                            CHR_var, ":", BP_START_var, "-",
                            BP_STOP_var), intern = T)
  if (identical(text_out, character(0))) {
    sumstats <- data.frame()
  } else {
    sumstats <- read.delim(text=text_out, header = F, stringsAsFactors = FALSE,
                           colClasses = c(V5 = "character", V6 = "character"))
  }
  if (nrow(sumstats) == 0) { return(sumstats) }
  # format by phenotype ID
  sumstats_list <- lapply(unique(sumstats$V8), function(x) {
    sumstats <- subset(sumstats, V8 == x)
    sumstats$Name <- paste0(sumstats$V1, ":", sumstats$V4, ":", sumstats$V5, ":", sumstats$V6)
    sumstats$rsids <- NA
    sumstats <- sumstats[,c("Name", "rsids", "V19", "V4", "V6", "V5", "V15", "V16", "V14", "V11", "V12", "V8")]
    colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N", "Phenotype")
    sumstats <- subset(sumstats, (!is.na(BETA)) & (!is.na(SE)) & (!is.na(P)))
    sumstats <- subset(sumstats, (! BETA %in% c(Inf, -Inf)) & (! SE %in% c(Inf, -Inf)))
    if (length(unique(sumstats$Phenotype)) > 1) {stop("Phenotype not unique in output query")}
    return(sumstats)
  })
  return(sumstats_list)
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
  sumstats <- unique(sumstats)
  if (nrow(sumstats) == 0) { return(sumstats) }
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

#' Query dbSNP_file to get REF and ALT
#' @param dbSNP_file path to dbSNP_file
#' @param CHR_var path to dbSNP_file
#' @param BP_START_var path to dbSNP_file
#' @param BP_STOP_var path to dbSNP_file
#' @return data frame with rs to REF-ALT mapping
#' @examples
#' under development
#' @export
query_dbSNP <- function(dbSNP_file,
                        CHR_var, BP_START_var, BP_STOP_var, rsID, ...) {
  system_out <- suppressWarnings(
    system(paste0("tabix -h ", dbSNP_file, " chr", CHR_var, ":",
                  BP_START_var, "-", BP_STOP_var, " | grep -wF ", rsID),
           intern = T)
  )
  if (is.null(attr(system_out,""))) {
    sumstats <- read.table(text=system_out, header = F)
  } else if (attr(system_out, "status") == 1) {
    sumstats <- data.frame(rsID = NA, Name = NA, REF = NA, ALT = NA)
    return(sumstats)
  }
  sumstats <- sumstats[,c(3,6)]
  colnames(sumstats) <- c("rsID", "Name")
  sumstats$REF <- gsub("chr[0-9]+:[0-9]+:(.*):.*", "\\1", sumstats$Name)
  sumstats$ALT <- gsub("chr[0-9]+:[0-9]+:.*:(.*)", "\\1", sumstats$Name)
  stopifnot(all(names(table(sapply(strsplit(sumstats$rsID, "rs"), length))) == "2"))
  return(sumstats)
}


