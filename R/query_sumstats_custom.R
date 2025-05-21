#' tabix_Icelanders_pGWAS
tabix_Icelanders_pGWAS <- function(sumstats_file, coloc_regions_PASS) {
  # add "chr" to coloc_regions_PASS
  coloc_regions_PASS$CHR_var <- paste0("chr", coloc_regions_PASS$CHR_var)
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                             coloc_regions_PASS=coloc_regions_PASS)
  # Icelander specific data processing
  assocvariants_file <- paste0(dirname(sumstats_file), "/assocvariants.annotated.txt.gz")
  if (!file.exists(assocvariants_file)) {
    stop(paste0("Check annotation file for Icelanders_pGWAS: ", assocvariants_file))
  }
  assocvariants <- retrieve_sumstats_tabix(sumstats_file=assocvariants_file,
                                  coloc_regions_PASS=coloc_regions_PASS)
  assocvariants <- assocvariants[,c("Name","effectAlleleFreq")]
  assocvariants_ord <- assocvariants[match(sumstats$Name, assocvariants$Name),]
  if (nrow(assocvariants_ord) != nrow(sumstats) | 
      ! all(assocvariants_ord$Name == sumstats$Name, na.rm=T)) {
    stop("Matching with annotation file for Icelanders_pGWAS failed")
  }
  sumstats$effectAlleleFreq <- assocvariants_ord$effectAlleleFreq
  sumstats <- subset(sumstats, !is.na(effectAlleleFreq))
  # Formatting sumstats
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_Icelanders_pGWAS(sumstats=sumstats)
  }
  return(sumstats)
}

#' format_Icelanders_pGWAS
format_Icelanders_pGWAS <- function(sumstats) {
  # check
  cols_from <- c("Name", "Chrom", "Pos", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "minus_log10_pval", "SE", "N", "ImpMAF", "effectAlleleFreq")
  # format
  sumstats$Chrom <- gsub("chr", "", sumstats$Chrom)
  sumstats$Name <- gsub("(.*:.*):(.*):(.*)", "\\1:\\3:\\2", sumstats$Name)
  # colnames
  sumstats <- match_cols(sumstats=sumstats,
                         Name="Name",
                         rsID="rsids",
                         CHR="Chrom",
                         POS="Pos",
                         A1="effectAllele",
                         A2="otherAllele",
                         BETA="Beta",
                         SE="SE",
                         nlog10P="minus_log10_pval",
                         AF="effectAlleleFreq",
                         N="N")
  return(sumstats)
}

#' tabix_UKB_PPP_EUR
tabix_UKB_PPP_EUR <- function(sumstats_file, coloc_regions_PASS) {
  # change X to 23
  coloc_regions_PASS$CHR_var[coloc_regions_PASS$CHR_var == "X"] <- "23"
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                             coloc_regions_PASS=coloc_regions_PASS)
  sumstats <- format_UKB_PPP_EUR(sumstats=sumstats)
  return(sumstats)
}

#' format_UKB_PPP_EUR
format_UKB_PPP_EUR <- function(sumstats) {
  # check
  cols_from <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")
  # format
  sumstats$rsID <- NA
  sumstats$CHROM[sumstats$CHROM == "23"] <- "X"
  sumstats$Name <- paste0("chr", sumstats$CHROM, ":", sumstats$GENPOS, ":", sumstats$ALLELE0, ":", sumstats$ALLELE1)
  # colnames
  sumstats <- match_cols(sumstats=sumstats,
                         Name="Name",
                         rsID="rsID",
                         CHR="CHROM",
                         POS="GENPOS",
                         A1="ALLELE1",
                         A2="ALLELE0",
                         BETA="BETA",
                         SE="SE",
                         nlog10P="LOG10P",
                         AF="A1FREQ",
                         N="N")
  return(sumstats)
}

#' tabix_UKB_TOPMed
tabix_UKB_TOPMed <- function(sumstats_file, coloc_regions_PASS) {
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                             coloc_regions_PASS=coloc_regions_PASS)
  sumstats <- format_UKB_TOPMed(sumstats=sumstats)
  return(sumstats)
}

#' format_UKB_TOPMed
format_UKB_TOPMed <- function(sumstats) {
  # check
  cols_from <- c("chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "consequence", "pval", "beta", "sebeta", "af", "case_af", "control_af", "tstat")
  # format
  sumstats$Name <- paste0("chr", sumstats$chrom, ":", sumstats$pos, ":", sumstats$ref, ":", sumstats$alt)
  sumstats$N <- NA
  sumstats$nlog10P <- -log10(sumstats$pval)
  # colnames
  sumstats <- match_cols(sumstats=sumstats,
                         Name="Name",
                         rsID="rsids",
                         CHR="chrom",
                         POS="pos",
                         A1="alt",
                         A2="ref",
                         BETA="beta",
                         SE="sebeta",
                         nlog10P="nlog10P",
                         AF="af",
                         N="N")
  return(sumstats)
  
}

#' tabix_FinnGen_r9
#' https://finngen.gitbook.io/documentation/data-description
tabix_FinnGen_r9 <- function(sumstats_file, coloc_regions_PASS) {
  coloc_regions_PASS$CHR_var[coloc_regions_PASS$CHR_var == "X"] <- "23"
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                             coloc_regions_PASS=coloc_regions_PASS)
  sumstats <- format_FinnGen_r9(sumstats=sumstats)
  return(sumstats)
}

#' format_FinnGen_r9
format_FinnGen_r9 <- function(sumstats) {
  # check
  cols_from <- c("#chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls")
  # format
  sumstats[["#chrom"]][sumstats[["#chrom"]] == "23"] <- "X"
  sumstats$Name <- paste0("chr", sumstats[["#chrom"]], ":", sumstats$pos, ":", sumstats$ref, ":", sumstats$alt)
  sumstats$N <- NA
  # colnames
  sumstats <- match_cols(sumstats=sumstats,
                         Name="Name",
                         rsID="rsids",
                         CHR="#chrom",
                         POS="pos",
                         A1="alt",
                         A2="ref",
                         BETA="beta",
                         SE="sebeta",
                         nlog10P="mlogp",
                         AF="af_alt",
                         N="N")
  return(sumstats)
}

#' tabix_eQTLGen
tabix_eQTLGen <- function(sumstats_file, coloc_regions_PASS) {
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                                      coloc_regions_PASS=coloc_regions_PASS)
  attr(sumstats, "sumstats_str") <- "multiple"
  sumstats <- format_eQTLGen(sumstats=sumstats)
  return(sumstats)
}

#' format_eQTLGen
format_eQTLGen <- function(sumstats) {
  # check
  cols_from <- c("Name","rsID","CHR","POS","A1","A2","BETA","SE","P","AF","N","Phenotype")
  if (!identical(cols_from, colnames(sumstats))) stop("Column mismatch when reading eQTLGen")
  # format
  sumstats$nlog10P <- -log10(sumstats$P)
  # colnames
  sumstats <- match_cols(sumstats=sumstats,
                         Name="Name",
                         rsID="rsID",
                         CHR="CHR",
                         POS="POS",
                         A1="A1",
                         A2="A2",
                         BETA="BETA",
                         SE="SE",
                         nlog10P="nlog10P",
                         AF="AF",
                         N="N",
                         Phenotype="Phenotype")
  return(sumstats)
}





