# TODO
# p_threshold p-value to include sumstats in coloc, default 1e-6
# maf allele frequency to filter out SNPs, default 0.01
# - Add check for reference allele (first allele should be reference)
# - Add check that A1 is alternative allele (first allele should be reference)
# checks: duplicated SNPs
# - Add input option for sdY
# - Split between mandatory and optional columns
#   - mandatory columns for binary traits: ID, BETA, SE 
#   - mandatory columns for quantitative traits: ID, BETA, SE 
# - Add to the file with internal examples
# - Add colClasses, e.g., `colClasses = c(alt = "character", ref= "character")`
# - Check there are no NAs in betas
# - Phenotype (genes) in GTEX: ENSG00000198744.5 is not a real gene. It is ENSG00000198744
# - Check chr X codification
# coloc: add sdY


query_Spanish_GWAS <- function(sumstats_file,
                              CHR_var, BP_START_var, BP_STOP_var) {
  sumstats <- read.csv(sumstats_file)
  sumstats<- subset(sumstats, CHR== CHR_var & pos>= BP_START_var & pos<=BP_STOP_var)
  if (nrow(sumstats) == 0) { return(NA) }
  # format
  sumstats$Name <- paste0("chr", sumstats$CHR, ":", sumstats$pos, ":", sumstats$REF.0, ":", sumstats$ALT.1)
  sumstats$rsids<- NA
  sumstats <- sumstats[,c("Name", "rsids", "CHR", "pos", "ALT.1.", "REF.0.", "beta", "se", "pvalue", "MAF", "N")]
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "AF", "N")
  # checks
  stopifnot(all(gsub(".*:.*:.*:(.*)", "\\1", sumstats$Name) == sumstats$A1 &
                  all(gsub(".*:.*:(.*):.*", "\\1", sumstats$Name) == sumstats$A2)))
  return(sumstats)
}



