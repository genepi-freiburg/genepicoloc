# =============================================================================
# Preprocessing GWAS Summary Statistics for genepicoloc
# =============================================================================
#
# This script demonstrates the complete preprocessing workflow:
# 1. Download summary statistics from public source
# 2. LiftOver coordinates from hg19 to hg38
# 3. Add rsID and standardized variant names (optional)
# 4. Calculate -log10(P) with underflow handling
# 5. Harmonize column names to genepicoloc standard format
#
# Example data: CKDGen Round 4 eGFR (Wuttke et al. 2019, PMID: 31152163)
#
# For detailed explanations, see the companion tutorial:
#   inst/scripts/preprocessing_tutorial.Rmd
#
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup and Dependencies
# -----------------------------------------------------------------------------

library(data.table)
# library(genepicoloc)

# Check external tool dependencies
check_dependencies <- function() {
  # Check tabix and bgzip (from htslib)
  tabix_path <- Sys.which("tabix")
  bgzip_path <- Sys.which("bgzip")

  if (tabix_path == "" || bgzip_path == "") {
    stop("tabix/bgzip not found. Install htslib:\n",
         "  Ubuntu/Debian: sudo apt install tabix\n",
         "  macOS: brew install htslib\n",
         "  conda: conda install -c bioconda htslib")
  }
  message("tabix: ", tabix_path)
  message("bgzip: ", bgzip_path)

  # Check liftOver
  liftOver_path <- Sys.which("liftOver")
  if (liftOver_path == "") {
    message("liftOver not found in PATH (will be downloaded automatically)")
  } else {
    message("liftOver: ", liftOver_path)
  }

  invisible(TRUE)
}

check_dependencies()

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------

# Working directory for downloads and output
# Using a subdirectory in user's home by default
work_dir <- file.path(Sys.getenv("HOME"), "genepicoloc_tutorial")
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
message("Working directory: ", work_dir)

# For testing, use only chromosome 21 (faster, ~113K variants)
# Set to NULL to process all chromosomes
test_chr <- 21

# Tool paths (modify if not in PATH)
liftOver_bin <- Sys.which("liftOver")
if (liftOver_bin == "") liftOver_bin <- file.path(work_dir, "liftOver")

liftOver_chain <- file.path(work_dir, "hg19ToHg38.over.chain.gz")

# Ensembl VCF for rsID lookup (optional)
ensembl_vcf <- file.path(work_dir, "homo_sapiens-chr21.vcf.gz")

# -----------------------------------------------------------------------------
# 2. Download Required Files
# -----------------------------------------------------------------------------

message("\n=== Downloading required files ===\n")

# 2a. Download liftOver binary (if not installed)
if (!file.exists(liftOver_bin)) {
  message("Downloading liftOver...")
  liftOver_url <- "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver"
  download.file(liftOver_url, liftOver_bin, mode = "wb")
  Sys.chmod(liftOver_bin, "755")
  message("  Downloaded: ", liftOver_bin)
}

# 2b. Download chain file
if (!file.exists(liftOver_chain)) {
  message("Downloading hg19 to hg38 chain file...")
  chain_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
  download.file(chain_url, liftOver_chain, mode = "wb")
  message("  Downloaded: ", liftOver_chain)
}

# 2c. Download Ensembl VCF for rsID lookup (optional, chr21 only for testing)
if (!file.exists(ensembl_vcf)) {
  message("Downloading Ensembl variation VCF (chr21)...")
  vcf_url <- "https://ftp.ensembl.org/pub/release-115/variation/vcf/homo_sapiens/homo_sapiens-chr21.vcf.gz"
  download.file(vcf_url, ensembl_vcf, mode = "wb")
  # Also download index
  download.file(paste0(vcf_url, ".csi"), paste0(ensembl_vcf, ".csi"), mode = "wb")
  message("  Downloaded: ", ensembl_vcf)
}

# 2d. Download CKDGen eGFR summary statistics
input_url <- "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.txt.gz"
input_file <- file.path(work_dir, basename(input_url))

if (!file.exists(input_file)) {
  message("Downloading CKDGen eGFR summary statistics...")
  message("  Source: Wuttke et al. 2019 (PMID: 31152163)")
  download.file(input_url, input_file, mode = "wb")
  message("  Downloaded: ", input_file)
}

# -----------------------------------------------------------------------------
# 3. Load and Inspect Data
# -----------------------------------------------------------------------------

message("\n=== Loading summary statistics ===\n")

sumstats <- fread(input_file)
message("Loaded ", format(nrow(sumstats), big.mark = ","), " variants")
message("Columns: ", paste(names(sumstats), collapse = ", "))

# Filter to test chromosome if specified
if (!is.null(test_chr)) {
  sumstats <- sumstats[Chr == test_chr]
  message("Filtered to chr", test_chr, ": ", format(nrow(sumstats), big.mark = ","), " variants")
}

# Show first few rows
message("\nFirst few rows:")
print(head(sumstats, 3))

# -----------------------------------------------------------------------------
# 4. LiftOver from hg19 to hg38
# -----------------------------------------------------------------------------

message("\n=== LiftOver: hg19 -> hg38 ===\n")

sumstats_hg38 <- genepi_liftover(

  sumstats = sumstats,
  CHR_name = "Chr",
  POS_name = "Pos_b37",
  A1_name = "Allele1",
  A2_name = "Allele2",
  liftOver_bin = liftOver_bin,
  liftOver_chain = liftOver_chain
)

message("After liftOver: ", format(nrow(sumstats_hg38), big.mark = ","), " variants")

# -----------------------------------------------------------------------------
# 5. Add Standardized Variant Names
# -----------------------------------------------------------------------------

message("\n=== Variant Name Harmonization ===\n")

# This step is required for creating standardized chr:pos:ref:alt names
# It also adds rsIDs for display

if (file.exists(ensembl_vcf)) {
  sumstats_hg38 <- name_by_position(
    sumstats = sumstats_hg38,
    CHR_name = "CHR_hg38",
    POS_name = "POS_hg38",
    A1_name = "A1_hg38",
    A2_name = "A2_hg38",
    dbSNP_file = ensembl_vcf
  )
  message("Variants with rsID: ", format(nrow(sumstats_hg38), big.mark = ","))
} else {
  # Create Name column manually without rsID lookup
  message("Skipping rsID lookup (Ensembl VCF not found)")
  sumstats_hg38[, Name_hg38 := paste0("chr", CHR_hg38, ":", POS_hg38, ":", A2_hg38, ":", A1_hg38)]
}

# -----------------------------------------------------------------------------
# 6. Calculate -log10(P) - required column
# -----------------------------------------------------------------------------

message("\n=== Calculating -log10(P) ===\n")

sumstats_hg38[, nlog10P := -log10(`P-value`)]

# Handle underflow (Inf values)
if (any(!is.finite(sumstats_hg38$nlog10P))) {
  max_finite <- max(sumstats_hg38$nlog10P[is.finite(sumstats_hg38$nlog10P)])
  n_inf <- sum(!is.finite(sumstats_hg38$nlog10P))
  sumstats_hg38[!is.finite(nlog10P), nlog10P := max_finite + 1]
  message("Fixed ", n_inf, " infinite values (set to ", round(max_finite + 1, 2), ")")
}

message("Max -log10(P): ", round(max(sumstats_hg38$nlog10P), 2))

# -----------------------------------------------------------------------------
# 7. Harmonize Column Names
# -----------------------------------------------------------------------------

message("\n=== Harmonizing column names ===\n")

# Map original columns to genepicoloc standard names
col_mapping <- c(
  "CHR" = "CHR_hg38",
  "POS" = "POS_hg38",
  "A1" = "A1_hg38",
  "A2" = "A2_hg38",
  "Name" = "Name_hg38",
  "AF" = "Freq1",
  "BETA" = "Effect",
  "SE" = "StdErr",
  "P" = "P-value",
  "nlog10P" = "nlog10P",
  "N" = "n_total_sum"
)

# Add rsID
col_mapping <- c(col_mapping, "rsID" = "rs")

# Select and rename columns
existing_cols <- col_mapping[col_mapping %in% names(sumstats_hg38)]
sumstats_final <- sumstats_hg38[, existing_cols, with = FALSE]
setnames(sumstats_final, existing_cols, names(existing_cols))

# Sort by position
setorder(sumstats_final, CHR, POS)

message("Final columns: ", paste(names(sumstats_final), collapse = ", "))
message("Final variant count: ", format(nrow(sumstats_final), big.mark = ","))

# -----------------------------------------------------------------------------
# 8. Save Output (bgzip + tabix indexed)
# -----------------------------------------------------------------------------

message("\n=== Saving output ===\n")

output_base <- file.path(work_dir, "CKDGen_eGFR_chr21_hg38.tsv")
output_file <- paste0(output_base, ".gz")

# Write uncompressed first
fwrite(sumstats_final, output_base, sep = "\t")

# Compress with bgzip (required for tabix indexing)
message("Compressing with bgzip...")
if (file.exists(output_file)) file.remove(output_file)
system2("bgzip", output_base)

# Create tabix index (by CHR and POS columns)
message("Creating tabix index...")
system2("tabix", c("-s", "1", "-b", "2", "-e", "2", "-S", "1", output_file))

message("Output file: ", output_file)
message("Index file:  ", paste0(output_file, ".tbi"))
message("File size: ", round(file.size(output_file) / 1024^2, 1), " MB")

# -----------------------------------------------------------------------------
# 9. Summary
# -----------------------------------------------------------------------------

message("\n", strrep("=", 60))
message("Preprocessing Complete!")
message(strrep("=", 60))
message("\nInput:  ", basename(input_file))
message("Output: ", basename(output_file))
message("Build:  GRCh38/hg38")
message("Variants: ", format(nrow(sumstats_final), big.mark = ","))

message("\nTop 5 signals:")
print(sumstats_final[order(-nlog10P)][1:5, .(CHR, POS, Name, BETA, nlog10P)])

message("\nNext steps:")
message("  1. Use match_cols() to prepare for colocalization")
message("  2. Run get_coloc_regions() to identify significant loci")
message("  3. See vignette('genepicoloc') for full workflow")
