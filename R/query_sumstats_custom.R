#' tabix_GTEXv8
tabix_GTEXv8 <- function(sumstats_file, coloc_regions_PASS) {
  coloc_regions_PASS$CHR_var <- paste0("chr", coloc_regions_PASS$CHR_var)
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                                      coloc_regions_PASS=coloc_regions_PASS,
                                      sumstats_pheno="multiple")
  # GTEXv8 does not have header, so if there are no lines
  # it will return Null data.table (0 rows and 0 cols)
  # which is handled as "tabix_failed"
  # fix "tabix_failed" in case there was no error
  tabix_cmd <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                                       coloc_regions_PASS=coloc_regions_PASS,
                                       file_remove=FALSE,
                                       return_tabix_cmd=TRUE)
  sys_out <- system(tabix_cmd, ignore.stdout=T)
  file.remove(attr(tabix_cmd, "file_regions"))
  if (identical(sys_out, 0L) && nrow(sumstats) == 0) {
    attr(sumstats, "tabix") <- "tabix_ok_no_data"
    sumstats_fixed <- data.table::data.table(matrix(ncol = 19, nrow=0))
    class(sumstats_fixed) <- c("sumstats", class(sumstats_fixed))
    for (i in c("tabix", "sumstats_file", "coloc_regions_PASS", "sumstats_pheno")) {
      attr(sumstats_fixed, i) <- attr(sumstats, i)
    }
    sumstats <- sumstats_fixed
  }
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_GTEXv8(sumstats=sumstats)
  }
  return(sumstats)
}

#' format_GTEXv8
format_GTEXv8 <- function(sumstats) {
  # check
  cols_from <- paste0("V", 1:19)
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
  # format
  sumstats$nlog10P <- -log10(sumstats$V14)
  sumstats$Name <- paste0(sumstats$V1, ":", sumstats$V4, ":", sumstats$V5, ":", sumstats$V6)
  sumstats$rsids <- NA
  # colnames
  sumstats <- match_cols(sumstats=sumstats,
                         Name="Name",
                         rsID="rsids",
                         CHR="V19",
                         POS="V4",
                         A1="V6",
                         A2="V5",
                         BETA="V15",
                         SE="V16",
                         nlog10P="nlog10P",
                         AF="V11",
                         N="V12",
                         Phenotype="V8")
  return(sumstats)
}

#' tabix_Kidney_eQTL
tabix_Kidney_eQTL <- function(sumstats_file, coloc_regions_PASS) {
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                                      coloc_regions_PASS=coloc_regions_PASS,
                                      sumstats_pheno="multiple")
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_Kidney_eQTL(sumstats=sumstats)
  }
  return(sumstats)
}

#' format_Kidney_eQTL
format_Kidney_eQTL <- function(sumstats) {
  # check
  cols_from <- c("CHR","POS_hg38","Name","GeneID","Gene_Symbol","rsID","SNP","Ref","Alt","Beta","Std","Pvalue","SNP_Location","Compartment")
  # handle special case
  if (basename(attr(sumstats, "sumstats_file")) == "Kidney_eQTL_Meta_S686_Significant.q0.01_hg38.txt.gz") {
    cols_from <- c("CHR","POS_hg38","Name","GeneSymbol","GeneID","rsID","SNP_Location","Ref","Alt","Direction","Beta","Std","Pvalue","SNP")
  }
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
  # format
  sumstats$nlog10P <- -log10(sumstats$Pvalue)
  sumstats$AF <- NA
  sumstats$N <- NA
  # colnames
  sumstats <- match_cols(sumstats=sumstats,
                         Name="Name",
                         rsID="rsID",
                         CHR="CHR",
                         POS="POS_hg38",
                         A1="Alt",
                         A2="Ref",
                         BETA="Beta",
                         SE="Std",
                         nlog10P="nlog10P",
                         AF="AF",
                         N="N",
                         Phenotype="GeneID")
  return(sumstats)
}

#' tabix_eQTLGen
tabix_eQTLGen <- function(sumstats_file, coloc_regions_PASS) {
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                                      coloc_regions_PASS=coloc_regions_PASS,
                                      sumstats_pheno="multiple")
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_eQTLGen(sumstats=sumstats)
  }
  return(sumstats)
}

#' format_eQTLGen
format_eQTLGen <- function(sumstats) {
  # check
  cols_from <- c("Name","rsID","CHR","POS","A1","A2","BETA","SE","P","AF","N","Phenotype")
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
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
  cols_from <- c("Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "minus_log10_pval", "SE", "N", "ImpMAF", "effectAlleleFreq")
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
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
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_UKB_PPP_EUR(sumstats=sumstats)
  }
  return(sumstats)
}

#' format_UKB_PPP_EUR
format_UKB_PPP_EUR <- function(sumstats) {
  # check
  cols_from <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
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
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_UKB_TOPMed(sumstats=sumstats)
  }
  return(sumstats)
}

#' format_UKB_TOPMed
format_UKB_TOPMed <- function(sumstats) {
  # check
  cols_from <- c("chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "consequence", "pval", "beta", "sebeta", "af", "case_af", "control_af", "tstat")
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
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


tabix_GCKD_pGWAS <- function(sumstats_file, coloc_regions_PASS) {
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                                      coloc_regions_PASS=coloc_regions_PASS)
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_GCKD_pGWAS(sumstats=sumstats)
    OID <- gsub(".*(OID[0-9]+).*", "\\1", sumstats_file)
    attr(sumstats, "sumstats_file") <- 
      paste0(attr(sumstats, "sumstats_file"), "_", OID)
  }
  return(sumstats)
}

format_GCKD_pGWAS <- function(sumstats) {
  # check
  cols_from <- c("SNP","chr","position","coded_all","noncoded_all","strand_genome","beta","SE","pval","AF_coded_all","HWE_pval","FreqSE","MinFreq","MaxFreq","Direction","HetISq","HetChiSq","HetDf","HetPVal","callrate","n_total","oevar_imp","imputed")
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
  # format
  pval <- sumstats[["pval"]]
  sumstats$nlog10P <- NA
  # handle underflow if p<1e-322
  if (is.character(pval)) {
    l1 <- grepl("e", pval)
    l2 <- !l1
    # Scientific notation cases
    sumstats$nlog10P[l1] <- sapply(pval[l1], function(x) {
      parts <- as.numeric(strsplit(x, "e")[[1]])
      mantissa <- parts[1]
      exponent <- parts[2]
      -log10(mantissa) - exponent
    })
    # Regular decimal cases
    sumstats$nlog10P[l2] <- -log10(as.numeric(pval[l2]))
  } else {
    sumstats$nlog10P <- -log10(sumstats$pval)
  }
  sumstats$rsID <- NA
  sumstats <- match_cols(sumstats=sumstats,
                         Name="SNP",
                         rsID="rsID",
                         CHR="chr",
                         POS="position",
                         A1="coded_all",
                         A2="noncoded_all",
                         BETA="beta",
                         SE="SE",
                         nlog10P="nlog10P",
                         AF="AF_coded_all",
                         N="n_total")
  return(sumstats)
}



tabix_GCKD_mGWAS <- function(sumstats_file, coloc_regions_PASS) {
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                                      coloc_regions_PASS=coloc_regions_PASS)
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_GCKD_mGWAS(sumstats=sumstats)
  }
  return(sumstats)
}

format_GCKD_mGWAS <- function(sumstats) {
  # check
  cols_from <- c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
  # format
  sumstats$rsID <- NA
  sumstats <- match_cols(sumstats=sumstats,
                         Name="ID",
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



#' https://finngen.gitbook.io/documentation/data-description
tabix_FinnGen_r9 <- function(sumstats_file, coloc_regions_PASS) {
  coloc_regions_PASS$CHR_var[coloc_regions_PASS$CHR_var == "X"] <- "23"
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                             coloc_regions_PASS=coloc_regions_PASS)
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_FinnGen_r9(sumstats=sumstats)
  }
  return(sumstats)
}

#' format_FinnGen_r9
format_FinnGen_r9 <- function(sumstats) {
  # check
  cols_from <- c("#chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls")
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
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


tabix_CKDGen_r4 <- function(sumstats_file, coloc_regions_PASS) {
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                                      coloc_regions_PASS=coloc_regions_PASS)
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_CKDGen_r4(sumstats=sumstats)
  }
  return(sumstats)
}

format_CKDGen_r4 <- function(sumstats) {
  # check
  cols_from <- c("CHR_hg38","POS_hg38","A1_hg38","A2_hg38","Name_hg38","RSID","Freq1","Effect","StdErr","P-value","n_total_sum")
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
  # format
  sumstats$nlog10P <- -log10(sumstats$`P-value`)
  sumstats <- match_cols(sumstats=sumstats,
                         Name="Name_hg38",
                         rsID="RSID",
                         CHR="CHR_hg38",
                         POS="POS_hg38",
                         A1="A1_hg38",
                         A2="A2_hg38",
                         BETA="Effect",
                         SE="StdErr",
                         nlog10P="nlog10P",
                         AF="Freq1",
                         N="n_total_sum")
  return(sumstats)
}


tabix_CKDGen_r4 <- function(sumstats_file, coloc_regions_PASS) {
  sumstats <- retrieve_sumstats_tabix(sumstats_file=sumstats_file,
                                      coloc_regions_PASS=coloc_regions_PASS)
  if (attr(sumstats, "tabix") != "tabix_failed") {
    sumstats <- format_CKDGen_r4(sumstats=sumstats)
  }
  return(sumstats)
}

format_CKDGen_r4 <- function(sumstats) {
  # check
  cols_from <- c("CHR_hg38","POS_hg38","A1_hg38","A2_hg38","Name_hg38","RSID","Freq1","Effect","StdErr","P-value","n_total_sum")
  if (!identical(cols_from, colnames(sumstats))) {
    stop("Column mismatch when reading ", attr(sumstats, "sumstats_file"))
  }
  # format
  sumstats$nlog10P <- -log10(sumstats$`P-value`)
  sumstats <- match_cols(sumstats=sumstats,
                         Name="Name_hg38",
                         rsID="RSID",
                         CHR="CHR_hg38",
                         POS="POS_hg38",
                         A1="A1_hg38",
                         A2="A2_hg38",
                         BETA="Effect",
                         SE="StdErr",
                         nlog10P="nlog10P",
                         AF="Freq1",
                         N="n_total_sum")
  return(sumstats)
}


