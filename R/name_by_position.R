#' Find SNP Names by Genomic Position
#'
#' Matches genomic variants to dbSNP identifiers and standardized names using 
#' chromosomal position and allele information. This function queries a tabix-indexed
#' dbSNP VCF file to find matching variants and harmonizes allele coding.
#'
#' @param sumstats A data.table containing summary statistics with genomic coordinates
#' @param CHR_name Character. Column name for chromosome (default: "CHR_hg38")
#' @param POS_name Character. Column name for genomic position (default: "POS_hg38")
#' @param A1_name Character. Column name for first allele (default: "A1_hg38")
#' @param A2_name Character. Column name for second allele (default: "A2_hg38")
#' @param tabix_bin Character. Path to tabix binary (default: "/usr/bin/tabix")
#' @param dbSNP_file Character. Path to tabix-indexed dbSNP VCF file
#' @param tmp_name Character. Temporary file prefix (default: NULL, uses tempfile())
#' @param do_sorting Logical. Whether to sort output by position (default: TRUE)
#' @param keep_lower Logical. Whether to keep lowercase alleles (default: FALSE)
#'
#' @return A data.table with additional columns:
#' \itemize{
#'   \item \strong{rs}: dbSNP rs identifiers
#'   \item \strong{Name_hg38}: Standardized variant names (CHR:POS:REF:ALT format)
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Validates input data and converts to appropriate formats
#'   \item Adds 'chr' prefix to chromosomes if missing (required for UCSC format)
#'   \item Converts alleles to uppercase (unless keep_lower=TRUE)
#'   \item Uses tabix to query dbSNP file for each chromosome
#'   \item Matches variants by position and allele combination
#'   \item Returns harmonized variant names and rs identifiers
#' }
#'
#' The function handles ambiguous allele orders by trying both A1/A2 and A2/A1 
#' combinations and provides matching statistics to help assess data quality.
#'
#' @section Required Files:
#' You need a tabix-indexed dbSNP VCF file. To create one:
#' \enumerate{
#'   \item Download dbSNP VCF from NCBI or Ensembl
#'   \item Process to create a 6-column format: CHR, POS, RS_ID, REF, ALT, NAME
#'   \item Compress with bgzip and index with tabix
#' }
#'
#' Expected file format (tab-separated):
#' \preformatted{
#' chr1    10001   rs1570391677   T   A   chr1:10001:T:A
#' chr1    10001   rs1570391677   T   C   chr1:10001:T:C
#' chr1    10002   rs1570391692   A   C   chr1:10002:A:C
#' }
#'
#' @section Dependencies:
#' Requires tabix binary to be installed and accessible. Install via:
#' \itemize{
#'   \item Ubuntu/Debian: \code{sudo apt-get install tabix}
#'   \item CentOS/RHEL: \code{sudo yum install tabix}
#'   \item macOS: \code{brew install htslib}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage with sample data
#' library(data.table)
#' 
#' # Create sample summary statistics
#' sumstats <- data.table(
#'   CHR_hg38 = c("1", "1", "2"),
#'   POS_hg38 = c(10001, 10002, 20001),
#'   A1_hg38 = c("T", "A", "G"),
#'   A2_hg38 = c("A", "C", "C")
#' )
#' 
#' # Find names by position
#' result <- name_by_position(
#'   sumstats = sumstats,
#'   dbSNP_file = "path/to/dbSNP_file.vcf.gz",
#'   tabix_bin = "/usr/bin/tabix"
#' )
#' 
#' # Check matching statistics in console output
#' print(result)
#' }
#'
#' @importFrom data.table fread fwrite setDT rbindlist setnames setorder copy
#' @export
name_by_position <- function(sumstats,
                             CHR_name = "CHR_hg38",
                             POS_name = "POS_hg38",
                             A1_name = "A1_hg38",
                             A2_name = "A2_hg38",
                             tabix_bin = "tabix",
                             dbSNP_file,
                             keep_lower = FALSE) {

  message("=== Variant Name Matching ===")

 # Validate dependencies
  if (Sys.which(tabix_bin) == "" && !file.exists(tabix_bin)) {
    stop("tabix not found. Install htslib or specify path.")
  }

  if (!file.exists(dbSNP_file)) {
    stop("dbSNP file not found: ", dbSNP_file)
  }

  # Work on a copy
  if (!data.table::is.data.table(sumstats)) {
    sumstats <- data.table::as.data.table(sumstats)
  } else {
    sumstats <- data.table::copy(sumstats)
  }

  n_input <- nrow(sumstats)
  message("Input: ", format(n_input, big.mark = ","), " variants")

  # Ensure position is numeric
  if (!is.numeric(sumstats[[POS_name]])) {
    sumstats[, (POS_name) := as.numeric(get(POS_name))]
  }

  # Detect chromosome naming in dbSNP file
  dbsnp_chroms <- system2(tabix_bin, c(dbSNP_file, "-l"), stdout = TRUE, stderr = FALSE)
  dbsnp_has_chr <- any(grepl("^chr", dbsnp_chroms))
  sumstats_has_chr <- any(grepl("^chr", sumstats[[CHR_name]]))

  # Adjust chromosome naming to match dbSNP
  if (dbsnp_has_chr && !sumstats_has_chr) {
    sumstats[, (CHR_name) := paste0("chr", get(CHR_name))]
  } else if (!dbsnp_has_chr && sumstats_has_chr) {
    sumstats[, (CHR_name) := gsub("^chr", "", get(CHR_name))]
  }

  # Handle allele case
  if (!keep_lower) {
    if (any(grepl("[a-z]", sumstats[[A1_name]])) || any(grepl("[a-z]", sumstats[[A2_name]]))) {
      message("Converting alleles to uppercase...")
      sumstats[, (A1_name) := toupper(get(A1_name))]
      sumstats[, (A2_name) := toupper(get(A2_name))]
    }
  }

  # Add row ID for tracking
  sumstats[, .row_id := .I]

  # Process each chromosome
  chromosomes <- unique(sumstats[[CHR_name]])
  message("Processing ", length(chromosomes), " chromosome(s)...")

  results_list <- list()

  for (chr in chromosomes) {
    sumstats_chr <- sumstats[get(CHR_name) == chr]

    # Get position range for tabix query
    pos_min <- min(sumstats_chr[[POS_name]])
    pos_max <- max(sumstats_chr[[POS_name]])
    region <- paste0(chr, ":", pos_min, "-", pos_max)

    # Query dbSNP with tabix
    dbsnp_raw <- tryCatch({
      system2(tabix_bin, c(dbSNP_file, region), stdout = TRUE, stderr = FALSE)
    }, error = function(e) character(0))

    if (length(dbsnp_raw) == 0) {
      message("  Chr ", chr, ": no dbSNP data in region")
      next
    }

    # Parse tabix output (VCF format: CHR, POS, ID, REF, ALT, ...)
    dbsnp <- data.table::fread(text = dbsnp_raw, header = FALSE, sep = "\t", select = 1:5,
                                col.names = c("CHR", "POS", "rsID", "REF", "ALT"))

    # Filter to exact positions in sumstats
    dbsnp <- dbsnp[POS %in% sumstats_chr[[POS_name]]]

    if (nrow(dbsnp) == 0) {
      message("  Chr ", chr, ": no position matches")
      next
    }

    # Create Name column (always with chr prefix)
    chr_prefix <- ifelse(grepl("^chr", dbsnp$CHR), dbsnp$CHR, paste0("chr", dbsnp$CHR))
    dbsnp[, Name_hg38 := paste0(chr_prefix, ":", POS, ":", REF, ":", ALT)]

    # Merge with sumstats
    merged <- merge(sumstats_chr, dbsnp,
                    by.x = POS_name, by.y = "POS",
                    all.x = FALSE, allow.cartesian = TRUE)

    # Match alleles (either orientation)
    matched <- merged[
      (REF == get(A1_name) & ALT == get(A2_name)) |
      (ALT == get(A1_name) & REF == get(A2_name))
    ]

    # Remove duplicates (keep first match per variant)
    matched <- matched[!duplicated(.row_id)]

    # Rename and clean columns
    data.table::setnames(matched, "rsID", "rs")
    matched[, c("CHR", "REF", "ALT") := NULL]

    n_matched <- nrow(matched)
    message("  Chr ", chr, ": matched ", format(n_matched, big.mark = ","), "/",
            format(nrow(sumstats_chr), big.mark = ","), " variants")

    results_list[[chr]] <- matched
  }

  # Combine results
  if (length(results_list) == 0) {
    warning("No variants matched. Check chromosome naming and positions.")
    return(sumstats[0])  # Return empty data.table with same structure
  }

  result <- data.table::rbindlist(results_list, fill = TRUE)

  # Clean up
  result[, .row_id := NULL]

  message("Output: ", format(nrow(result), big.mark = ","), " variants with rsID and Name")

  return(result)
}

#' Genomic Coordinate Lift Over with Variant Harmonization
#'
#' Performs genomic coordinate conversion between genome builds (e.g., hg19 to hg38)
#' using UCSC liftOver tool, with optional variant name harmonization via dbSNP.
#'
#' @param sumstats A data.table containing summary statistics with genomic coordinates
#' @param CHR_name Character. Column name for chromosome
#' @param POS_name Character. Column name for genomic position  
#' @param A1_name Character. Column name for first allele
#' @param A2_name Character. Column name for second allele
#' @param liftOver_bin Character. Path to liftOver binary (default: "liftOver")
#' @param liftOver_chain Character. Path to liftOver chain file
#' @param keep_lower Logical. Whether to keep lowercase alleles (default: FALSE)
#' @param rm_tmp_liftOver Logical. Whether to remove temporary files (default: TRUE)
#'
#' @return A data.table with lifted coordinates in new columns:
#' \itemize{
#'   \item \strong{CHR_hg38}: Chromosome in target build
#'   \item \strong{POS_hg38}: Position in target build
#'   \item \strong{A1_hg38}: First allele (strand-flipped if necessary)
#'   \item \strong{A2_hg38}: Second allele (strand-flipped if necessary)
#' }
#'
#' @note For variant name harmonization (adding Name_hg38 and rsID columns),
#' use \code{\link{name_by_position}} as a separate step after liftOver.
#'
#' @details
#' This function performs comprehensive genomic coordinate conversion:
#' \enumerate{
#'   \item Validates input data and required external tools
#'   \item Converts coordinates using UCSC liftOver
#'   \item Handles strand flipping for negative strand variants
#'   \item Optionally harmonizes variant names using dbSNP
#'   \item Returns cleaned data with new genomic coordinates
#' }
#'
#' Variants that cannot be lifted over are automatically removed with a summary
#' of the number of failed variants.
#'
#' @section Required Tools:
#' 
#' **liftOver Binary:**
#' Download from UCSC Genome Browser:
#' \url{https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver}
#' 
#' **Chain Files:**
#' Download appropriate chain file from:
#' \url{https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/}
#' 
#' Common chain files:
#' \itemize{
#'   \item hg19 to hg38: \code{hg19ToHg38.over.chain.gz}
#'   \item hg18 to hg38: \code{hg18ToHg38.over.chain.gz}
#'   \item hg38 to hg19: \code{hg38ToHg19.over.chain.gz}
#' }
#' 
#' **Installation:**
#' \preformatted{
#' # Download liftOver binary
#' wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
#' chmod +x liftOver
#' 
#' # Download chain file
#' wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
#' }
#'
#' @section dbSNP File Setup:
#' For variant name harmonization, you need a processed dbSNP file. See 
#' \code{\link{name_by_position}} for details on file format and creation.
#'
#' @examples
#' \dontrun{
#' # Basic liftOver from hg19 to hg38
#' library(data.table)
#' 
#' # Sample data in hg19 coordinates
#' sumstats_hg19 <- data.table(
#'   CHR = c("1", "2", "3"),
#'   POS = c(1000000, 2000000, 3000000),
#'   A1 = c("A", "G", "T"),
#'   A2 = c("G", "C", "C"),
#'   PVAL = c(1e-8, 1e-6, 1e-5)
#' )
#' 
#' # Perform liftOver
#' sumstats_hg38 <- genepi_liftover(
#'   sumstats = sumstats_hg19,
#'   CHR_name = "CHR",
#'   POS_name = "POS", 
#'   A1_name = "A1",
#'   A2_name = "A2",
#'   liftOver_bin = "path/to/liftOver",
#'   liftOver_chain = "path/to/hg19ToHg38.over.chain.gz"
#' )
#' 
#' # Check results
#' print(sumstats_hg38)
#' }
#'
#' @importFrom data.table fread fwrite setDT copy set
#' @export
genepi_liftover <- function(sumstats,
                            CHR_name,
                            POS_name,
                            A1_name,
                            A2_name,
                            liftOver_bin = "liftOver",
                            liftOver_chain,
                            keep_lower = FALSE) {

  message("=== Genomic Coordinate Lift Over ===")

  # Validate required tools and files
  if (!file.exists(liftOver_bin) && Sys.which(liftOver_bin) == "") {
    stop("liftOver binary not found: ", liftOver_bin,
         "\nDownload from: https://hgdownload.soe.ucsc.edu/admin/exe/")
  }

  if (!file.exists(liftOver_chain)) {
    stop("Chain file not found: ", liftOver_chain,
         "\nDownload from: https://hgdownload.soe.ucsc.edu/goldenPath/")
  }

  # Work on a copy to avoid modifying original
  if (!data.table::is.data.table(sumstats)) {
    sumstats <- data.table::as.data.table(sumstats)
  } else {
    sumstats <- data.table::copy(sumstats)
  }

  n_input <- nrow(sumstats)
  message("Input: ", format(n_input, big.mark = ","), " variants")

  # Create temp files using proper tempfile()
  bed_file <- tempfile(pattern = "liftover_input_", fileext = ".bed")
  output_file <- tempfile(pattern = "liftover_output_", fileext = ".bed")
  unmapped_file <- tempfile(pattern = "liftover_unmapped_", fileext = ".bed")

  # Ensure cleanup on exit (even if function fails)
  on.exit({
    unlink(c(bed_file, output_file, unmapped_file), force = TRUE)
  }, add = TRUE)

  # Add row index for tracking through liftOver
  row_id_col <- ".liftover_row_id"
  sumstats[, (row_id_col) := .I]

  # Handle chromosome naming (liftOver expects 'chr' prefix)
  chr_had_prefix <- any(grepl("^chr", sumstats[[CHR_name]], ignore.case = TRUE))
  if (!chr_had_prefix) {
    sumstats[, (CHR_name) := paste0("chr", get(CHR_name))]
  }

  # Handle allele case
  if (!keep_lower) {
    if (any(grepl("[a-z]", sumstats[[A1_name]])) || any(grepl("[a-z]", sumstats[[A2_name]]))) {
      message("Converting alleles to uppercase...")
      sumstats[, (A1_name) := toupper(get(A1_name))]
      sumstats[, (A2_name) := toupper(get(A2_name))]
    }
  }

  # Write BED file for liftOver (BED format: chr, start, end, name, score, strand)
  bed_data <- sumstats[, .(
    chr = get(CHR_name),
    start = get(POS_name),
    end = get(POS_name),
    name = get(row_id_col),
    score = ".",
    strand = "+"
  )]
  data.table::fwrite(bed_data, bed_file, sep = "\t", col.names = FALSE)

  # Execute liftOver
  message("Running liftOver...")
  exit_code <- system2(
    liftOver_bin,
    args = c(bed_file, liftOver_chain, output_file, unmapped_file),
    stdout = FALSE, stderr = FALSE
  )

  if (exit_code != 0) {
    stop("liftOver failed with exit code ", exit_code)
  }

  # Read lifted coordinates
  lifted <- data.table::fread(output_file, header = FALSE,
                               col.names = c("CHR_hg38", "POS_hg38", "end", row_id_col, "score", "strand_hg38"))
  lifted[, c("end", "score") := NULL]

  # Read unmapped variants (if any)
  n_unmapped <- 0

  if (file.exists(unmapped_file) && file.size(unmapped_file) > 0) {
    unmapped <- tryCatch({
      data.table::fread(unmapped_file, header = FALSE, fill = TRUE)
    }, error = function(e) NULL)
    if (!is.null(unmapped) && nrow(unmapped) > 0) {
      # Unmapped file has comment lines starting with #
      unmapped <- unmapped[!grepl("^#", V1)]
      n_unmapped <- nrow(unmapped)
    }
  }

  # Report statistics
  n_lifted <- nrow(lifted)
  message("Lifted: ", format(n_lifted, big.mark = ","), " variants (",
          round(100 * n_lifted / n_input, 1), "%)")
  if (n_unmapped > 0) {
    message("Unmapped: ", format(n_unmapped, big.mark = ","), " variants (",
            round(100 * n_unmapped / n_input, 2), "%)")
  }

  # Merge lifted coordinates back to original data
  sumstats_lifted <- merge(sumstats, lifted, by = row_id_col, all = FALSE)

  # Handle strand flipping for negative strand
  n_flipped <- sum(sumstats_lifted$strand_hg38 == "-")
  if (n_flipped > 0) {
    message("Flipping ", format(n_flipped, big.mark = ","), " variants on negative strand...")
    flip_idx <- sumstats_lifted$strand_hg38 == "-"
    sumstats_lifted[flip_idx, A1_hg38 := flip_alleles(get(A1_name))]
    sumstats_lifted[flip_idx, A2_hg38 := flip_alleles(get(A2_name))]
    sumstats_lifted[!flip_idx, A1_hg38 := get(A1_name)]
    sumstats_lifted[!flip_idx, A2_hg38 := get(A2_name)]
  } else {
    sumstats_lifted[, A1_hg38 := get(A1_name)]
    sumstats_lifted[, A2_hg38 := get(A2_name)]
  }

  # Clean up: remove helper columns, restore chr naming
  sumstats_lifted[, c(row_id_col, "strand_hg38") := NULL]

  # Remove 'chr' prefix from output columns
  sumstats_lifted[, CHR_hg38 := gsub("^chr", "", CHR_hg38)]

  # Restore original chr column if we added prefix
  if (!chr_had_prefix) {
    sumstats_lifted[, (CHR_name) := gsub("^chr", "", get(CHR_name))]
  }

  message("Output: ", format(nrow(sumstats_lifted), big.mark = ","), " variants with hg38 coordinates")

  return(sumstats_lifted)
}

#' Create dbSNP Reference File for Position-Based Matching
#'
#' Helper function to create a properly formatted dbSNP file for use with
#' \code{\link{name_by_position}}. This function processes raw dbSNP VCF files
#' into the required 6-column format.
#'
#' @param input_vcf Character. Path to input dbSNP VCF file (can be gzipped)
#' @param output_file Character. Path for output processed file
#' @param compress Logical. Whether to compress output with bgzip (default: TRUE)
#' @param index Logical. Whether to create tabix index (default: TRUE)
#' @param tabix_bin Character. Path to tabix binary (default: "/usr/bin/tabix")
#' @param bgzip_bin Character. Path to bgzip binary (default: "/usr/bin/bgzip")
#'
#' @return Character. Path to the created file
#'
#' @details
#' This function processes a standard dbSNP VCF file to create the required
#' 6-column format for position-based matching:
#' \enumerate{
#'   \item CHR - Chromosome 
#'   \item POS - Position
#'   \item RS_ID - dbSNP rs identifier
#'   \item REF - Reference allele
#'   \item ALT - Alternative allele (one per line)
#'   \item NAME - Standardized name (CHR:POS:REF:ALT)
#' }
#'
#' For multi-allelic variants, each alternative allele gets its own line.
#'
#' @section Data Sources:
#' You can download dbSNP VCF files from:
#' \itemize{
#'   \item NCBI dbSNP: \url{https://ftp.ncbi.nih.gov/snp/}
#'   \item Ensembl: \url{https://ftp.ensembl.org/pub/current_variation/vcf/}
#' }
#'
#' @examples
#' \dontrun{
#' # Process dbSNP VCF file
#' create_dbsnp_file(
#'   input_vcf = "dbSNP_155.vcf.gz",
#'   output_file = "dbSNP_processed.vcf.gz"
#' )
#' }
#'
#' @importFrom data.table fread fwrite
#' @export
create_dbsnp_file <- function(input_vcf,
                              output_file,
                              compress = TRUE,
                              index = TRUE,
                              tabix_bin = "/usr/bin/tabix",
                              bgzip_bin = "/usr/bin/bgzip") {
  
  message("Processing dbSNP VCF file: ", input_vcf)
  
  # Read VCF file
  vcf_data <- data.table::fread(input_vcf, skip = "^#CHROM")
  
  # Process multi-allelic variants
  processed_data <- vcf_data[, {
    alt_alleles <- strsplit(ALT, ",")[[1]]
    lapply(alt_alleles, function(alt) {
      name <- paste(`#CHROM`, POS, REF, alt, sep = ":")
      data.table(
        CHR = `#CHROM`,
        POS = POS,
        RS_ID = ID,
        REF = REF,
        ALT = alt,
        NAME = name
      )
    })
  }, by = seq_len(nrow(vcf_data))]
  
  # Flatten the results
  final_data <- rbindlist(processed_data$V1)
  
  # Write output
  temp_file <- paste0(output_file, ".tmp")
  data.table::fwrite(final_data, temp_file, sep = "\t", col.names = FALSE)
  
  if (compress) {
    message("Compressing with bgzip...")
    system(paste(bgzip_bin, temp_file))
    file.rename(paste0(temp_file, ".gz"), output_file)
  } else {
    file.rename(temp_file, output_file)
  }
  
  if (index && compress) {
    message("Creating tabix index...")
    system(paste(tabix_bin, "-s", "1", "-b", "2", "-e", "2", output_file))
  }
  
  message("Created dbSNP file: ", output_file)
  return(output_file)
}
