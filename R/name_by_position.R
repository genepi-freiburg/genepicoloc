#' Find SNP Names by Genomic Position
#'
#' Matches genomic variants to dbSNP identifiers and standardized names using
#' chromosomal position and allele information. This function queries a tabix-indexed
#' dbSNP VCF file to find matching variants and harmonizes allele coding.
#'
#' @param sumstats A data.frame or data.table containing summary statistics with genomic coordinates
#' @param CHR_name Character. Column name for chromosome (default: "CHR_hg38")
#' @param POS_name Character. Column name for genomic position (default: "POS_hg38")
#' @param A1_name Character. Column name for first allele (default: "A1_hg38")
#' @param A2_name Character. Column name for second allele (default: "A2_hg38")
#' @param tabix_bin Character. Path to tabix binary (default: "tabix")
#' @param dbSNP_file Character. Path to single tabix-indexed dbSNP VCF file (optional if dbSNP_dir provided)
#' @param dbSNP_dir Character. Path to directory containing per-chromosome VCF files
#'   named homo_sapiens-chr{CHR}.vcf.gz (optional if dbSNP_file provided)
#' @param keep_lower Logical. Whether to keep lowercase alleles (default: FALSE)
#' @param chunk_size Integer. Number of variants to process per chunk for memory efficiency.
#'   Set to NULL to disable chunking. Default is 1000000 (1M variants).
#'
#' @return A data.frame with additional columns:
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
#' For large datasets (>1M variants), the function processes data in chunks to
#' reduce memory usage. Each chunk is processed independently and results are
#' combined at the end.
#'
#' The function handles ambiguous allele orders by trying both A1/A2 and A2/A1
#' combinations and provides matching statistics to help assess data quality.
#'
#' @section dbSNP File Options:
#' You can provide either:
#' \itemize{
#'   \item \strong{dbSNP_file}: Single VCF file containing all chromosomes
#'   \item \strong{dbSNP_dir}: Directory with per-chromosome files (homo_sapiens-chr{CHR}.vcf.gz)
#' }
#'
#' Per-chromosome files can be downloaded from Ensembl:
#' \preformatted{
#' for chr in {1..22} X Y MT; do
#'   wget https://ftp.ensembl.org/pub/release-115/variation/vcf/homo_sapiens/homo_sapiens-chr$chr.vcf.gz
#'   wget https://ftp.ensembl.org/pub/release-115/variation/vcf/homo_sapiens/homo_sapiens-chr$chr.vcf.gz.csi
#' done
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
#' # Example with single file
#' result <- name_by_position(
#'   sumstats = sumstats,
#'   dbSNP_file = "path/to/dbSNP_all.vcf.gz"
#' )
#'
#' # Example with per-chromosome directory
#' result <- name_by_position(
#'   sumstats = sumstats,
#'   dbSNP_dir = "/path/to/ensembl_vcf/"
#' )
#'
#' # Process very large file with smaller chunks to reduce memory
#' result <- name_by_position(
#'   sumstats = large_sumstats,
#'   dbSNP_dir = "/path/to/ensembl_vcf/",
#'   chunk_size = 500000
#' )
#' }
#'
#' @importFrom data.table fread
#' @export
name_by_position <- function(sumstats,
                             CHR_name = "CHR_hg38",
                             POS_name = "POS_hg38",
                             A1_name = "A1_hg38",
                             A2_name = "A2_hg38",
                             tabix_bin = "tabix",
                             dbSNP_file = NULL,
                             dbSNP_dir = NULL,
                             keep_lower = FALSE,
                             chunk_size = 1000000) {

  message("=== Variant Name Matching ===")

  # Validate dependencies
  if (Sys.which(tabix_bin) == "" && !file.exists(tabix_bin)) {
    stop("tabix not found. Install htslib or specify path.")
  }

  # Validate dbSNP source - need either file or directory
  use_per_chr <- FALSE
  if (!is.null(dbSNP_dir)) {
    if (!dir.exists(dbSNP_dir)) {
      stop("dbSNP directory not found: ", dbSNP_dir)
    }
    use_per_chr <- TRUE
    message("Using per-chromosome VCF files from: ", dbSNP_dir)
  } else if (!is.null(dbSNP_file)) {
    if (!file.exists(dbSNP_file)) {
      stop("dbSNP file not found: ", dbSNP_file)
    }
    message("Using single VCF file: ", dbSNP_file)
  } else {
    stop("Must provide either dbSNP_file or dbSNP_dir")
  }

  # Work on a copy as data.frame
  df <- as.data.frame(sumstats)
  n_input <- nrow(df)
  message("Input: ", format(n_input, big.mark = ","), " variants")

  # Ensure position is numeric
  if (!is.numeric(df[[POS_name]])) {
    df[[POS_name]] <- as.numeric(df[[POS_name]])
  }

  # Handle allele case
  if (!keep_lower) {
    if (any(grepl("[a-z]", df[[A1_name]])) || any(grepl("[a-z]", df[[A2_name]]))) {
      message("Converting alleles to uppercase...")
      df[[A1_name]] <- toupper(df[[A1_name]])
      df[[A2_name]] <- toupper(df[[A2_name]])
    }
  }

  # Add row ID for tracking
  df$.row_id <- seq_len(nrow(df))

  # Determine if chunking is needed
  use_chunking <- !is.null(chunk_size) && n_input > chunk_size

  if (use_chunking) {
    # Process in chunks
    n_chunks <- ceiling(n_input / chunk_size)
    message("Processing in ", n_chunks, " chunks of ~", format(chunk_size, big.mark = ","), " variants each")

    chunk_indices <- split(seq_len(n_input), ceiling(seq_len(n_input) / chunk_size))
    chunk_results <- list()

    for (i in seq_along(chunk_indices)) {
      message("\n--- Chunk ", i, "/", n_chunks, " ---")
      df_chunk <- df[chunk_indices[[i]], , drop = FALSE]

      chunk_result <- .process_variants_chunk(
        df_chunk = df_chunk,
        CHR_name = CHR_name,
        POS_name = POS_name,
        A1_name = A1_name,
        A2_name = A2_name,
        tabix_bin = tabix_bin,
        dbSNP_file = dbSNP_file,
        dbSNP_dir = dbSNP_dir,
        use_per_chr = use_per_chr
      )

      if (!is.null(chunk_result) && nrow(chunk_result) > 0) {
        chunk_results[[i]] <- chunk_result
      }

      # Clean up to free memory
      rm(df_chunk, chunk_result)
      gc(verbose = FALSE)
    }

    # Combine chunk results
    if (length(chunk_results) == 0) {
      warning("No variants matched. Check chromosome naming and positions.")
      df$.row_id <- NULL
      return(df[0, , drop = FALSE])
    }

    result <- do.call(rbind, chunk_results)

  } else {
    # Process all at once
    result <- .process_variants_chunk(
      df_chunk = df,
      CHR_name = CHR_name,
      POS_name = POS_name,
      A1_name = A1_name,
      A2_name = A2_name,
      tabix_bin = tabix_bin,
      dbSNP_file = dbSNP_file,
      dbSNP_dir = dbSNP_dir,
      use_per_chr = use_per_chr
    )

    if (is.null(result) || nrow(result) == 0) {
      warning("No variants matched. Check chromosome naming and positions.")
      df$.row_id <- NULL
      return(df[0, , drop = FALSE])
    }
  }

  rownames(result) <- NULL
  result$.row_id <- NULL

  message("\nOutput: ", format(nrow(result), big.mark = ","), " variants with rsID and Name")

  return(result)
}

#' Process a chunk of variants for name matching
#'
#' Internal function that processes a subset of variants against dbSNP.
#'
#' @param df_chunk Data frame chunk to process
#' @param CHR_name Column name for chromosome
#' @param POS_name Column name for position
#' @param A1_name Column name for allele 1
#' @param A2_name Column name for allele 2
#' @param tabix_bin Path to tabix binary
#' @param dbSNP_file Path to single dbSNP file (or NULL)
#' @param dbSNP_dir Path to per-chromosome dbSNP directory (or NULL)
#' @param use_per_chr Whether to use per-chromosome files
#'
#' @return Data frame with matched variants, or NULL if no matches
#' @keywords internal
.process_variants_chunk <- function(df_chunk, CHR_name, POS_name, A1_name, A2_name,
                                     tabix_bin, dbSNP_file, dbSNP_dir, use_per_chr) {

  chromosomes <- unique(df_chunk[[CHR_name]])
  message("Processing ", length(chromosomes), " chromosome(s)...")

  results_list <- list()

  for (chr_val in chromosomes) {
    # Subset for this chromosome
    chr_idx <- df_chunk[[CHR_name]] == chr_val
    df_chr <- df_chunk[chr_idx, , drop = FALSE]
    n_chr <- nrow(df_chr)

    # Determine dbSNP file for this chromosome
    if (use_per_chr) {
      chr_for_file <- gsub("^chr", "", chr_val)
      dbSNP_file_chr <- file.path(dbSNP_dir, paste0("homo_sapiens-chr", chr_for_file, ".vcf.gz"))

      if (!file.exists(dbSNP_file_chr)) {
        message("  Chr ", chr_val, ": VCF file not found, skipping (", basename(dbSNP_file_chr), ")")
        next
      }
    } else {
      dbSNP_file_chr <- dbSNP_file
    }

    # Detect chromosome naming in dbSNP file
    dbsnp_chroms <- system2(tabix_bin, c(dbSNP_file_chr, "-l"), stdout = TRUE, stderr = FALSE)
    dbsnp_has_chr <- any(grepl("^chr", dbsnp_chroms))
    sumstats_chr_has_chr <- grepl("^chr", chr_val)

    # Adjust chromosome for tabix query
    chr_query <- chr_val
    if (dbsnp_has_chr && !sumstats_chr_has_chr) {
      chr_query <- paste0("chr", chr_val)
    } else if (!dbsnp_has_chr && sumstats_chr_has_chr) {
      chr_query <- gsub("^chr", "", chr_val)
    }

    # Get position range for tabix query
    pos_min <- min(df_chr[[POS_name]])
    pos_max <- max(df_chr[[POS_name]])
    region <- paste0(chr_query, ":", pos_min, "-", pos_max)

    message("  Chr ", chr_val, ": querying ", format(n_chr, big.mark = ","),
            " variants (", format(pos_min, big.mark = ","), "-", format(pos_max, big.mark = ","), ")...")

    # Query dbSNP with tabix
    dbsnp_raw <- tryCatch({
      system2(tabix_bin, c(dbSNP_file_chr, region), stdout = TRUE, stderr = FALSE)
    }, error = function(e) character(0))

    if (length(dbsnp_raw) == 0) {
      message("  Chr ", chr_val, ": no dbSNP data in region")
      next
    }

    # Parse tabix output
    dbsnp <- data.table::fread(text = dbsnp_raw, header = FALSE, sep = "\t", select = 1:5,
                                col.names = c("CHR", "POS", "rsID", "REF", "ALT"))
    dbsnp <- as.data.frame(dbsnp)

    # Filter to exact positions in sumstats
    dbsnp <- dbsnp[dbsnp$POS %in% df_chr[[POS_name]], , drop = FALSE]

    if (nrow(dbsnp) == 0) {
      message("  Chr ", chr_val, ": no position matches")
      next
    }

    # Expand multi-allelic variants
    if (any(grepl(",", dbsnp$ALT))) {
      expanded_list <- lapply(seq_len(nrow(dbsnp)), function(i) {
        row <- dbsnp[i, , drop = FALSE]
        alts <- unlist(strsplit(row$ALT, ","))
        data.frame(
          CHR = row$CHR,
          POS = row$POS,
          rsID = row$rsID,
          REF = row$REF,
          ALT = alts,
          stringsAsFactors = FALSE
        )
      })
      dbsnp <- do.call(rbind, expanded_list)
    }

    # Create Name column (always with chr prefix)
    chr_prefix <- ifelse(grepl("^chr", dbsnp$CHR), dbsnp$CHR, paste0("chr", dbsnp$CHR))
    dbsnp$Name_hg38 <- paste0(chr_prefix, ":", dbsnp$POS, ":", dbsnp$REF, ":", dbsnp$ALT)

    # Merge with sumstats
    merged <- merge(df_chr, dbsnp, by.x = POS_name, by.y = "POS", all.x = FALSE)

    if (nrow(merged) == 0) {
      message("  Chr ", chr_val, ": no position matches after merge")
      next
    }

    # Match alleles (either orientation)
    match_forward <- merged$REF == merged[[A1_name]] & merged$ALT == merged[[A2_name]]
    match_reverse <- merged$ALT == merged[[A1_name]] & merged$REF == merged[[A2_name]]
    matched <- merged[match_forward | match_reverse, , drop = FALSE]

    # Remove duplicates (keep first match per variant)
    matched <- matched[!duplicated(matched$.row_id), , drop = FALSE]

    # Rename rsID to rs and remove helper columns
    names(matched)[names(matched) == "rsID"] <- "rs"
    matched$CHR <- NULL
    matched$REF <- NULL
    matched$ALT <- NULL

    n_matched <- nrow(matched)
    pct_matched <- round(100 * n_matched / n_chr, 1)
    message("  Chr ", chr_val, ": matched ", format(n_matched, big.mark = ","), "/",
            format(n_chr, big.mark = ","), " variants (", pct_matched, "%)")

    results_list[[chr_val]] <- matched
  }

  # Combine results
  if (length(results_list) == 0) {
    return(NULL)
  }

  do.call(rbind, results_list)
}

#' Genomic Coordinate Lift Over with Variant Harmonization
#'
#' Performs genomic coordinate conversion between genome builds (e.g., hg19 to hg38)
#' using UCSC liftOver tool, with optional variant name harmonization via dbSNP.
#'
#' @param sumstats A data.frame or data.table containing summary statistics with genomic coordinates
#' @param CHR_name Character. Column name for chromosome
#' @param POS_name Character. Column name for genomic position
#' @param A1_name Character. Column name for first allele
#' @param A2_name Character. Column name for second allele
#' @param liftOver_bin Character. Path to liftOver binary (default: "liftOver")
#' @param liftOver_chain Character. Path to liftOver chain file
#' @param keep_lower Logical. Whether to keep lowercase alleles (default: FALSE)
#'
#' @return A data.frame with lifted coordinates in new columns:
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
#' @examples
#' \dontrun{
#' # Basic liftOver from hg19 to hg38
#'
#' # Sample data in hg19 coordinates
#' sumstats_hg19 <- data.frame(
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
#' @importFrom data.table fread fwrite
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

  # Work on a copy as data.frame
  df <- as.data.frame(sumstats)

  n_input <- nrow(df)
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
  df$.liftover_row_id <- seq_len(nrow(df))

  # Handle chromosome naming (liftOver expects 'chr' prefix)
  chr_had_prefix <- any(grepl("^chr", df[[CHR_name]], ignore.case = TRUE))
  if (!chr_had_prefix) {
    df[[CHR_name]] <- paste0("chr", df[[CHR_name]])
  }

  # Handle allele case
  if (!keep_lower) {
    if (any(grepl("[a-z]", df[[A1_name]])) || any(grepl("[a-z]", df[[A2_name]]))) {
      message("Converting alleles to uppercase...")
      df[[A1_name]] <- toupper(df[[A1_name]])
      df[[A2_name]] <- toupper(df[[A2_name]])
    }
  }

  # Write BED file for liftOver (BED format: chr, start, end, name, score, strand)
  bed_data <- data.frame(
    chr = df[[CHR_name]],
    start = df[[POS_name]],
    end = df[[POS_name]],
    name = df$.liftover_row_id,
    score = ".",
    strand = "+",
    stringsAsFactors = FALSE
  )
  write.table(bed_data, bed_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

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
                               col.names = c("CHR_hg38", "POS_hg38", "end", ".liftover_row_id", "score", "strand_hg38"))
  lifted <- as.data.frame(lifted)
  lifted$end <- NULL
  lifted$score <- NULL

  # Read unmapped variants (if any)
  n_unmapped <- 0
  if (file.exists(unmapped_file) && file.size(unmapped_file) > 0) {
    unmapped_lines <- readLines(unmapped_file, warn = FALSE)
    # Filter out comment lines
    unmapped_lines <- unmapped_lines[!grepl("^#", unmapped_lines)]
    n_unmapped <- length(unmapped_lines)
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
  df_lifted <- merge(df, lifted, by = ".liftover_row_id", all = FALSE)

  # Handle strand flipping for negative strand
  flip_idx <- df_lifted$strand_hg38 == "-"
  n_flipped <- sum(flip_idx)

  if (n_flipped > 0) {
    message("Flipping ", format(n_flipped, big.mark = ","), " variants on negative strand...")
    df_lifted$A1_hg38 <- df_lifted[[A1_name]]
    df_lifted$A2_hg38 <- df_lifted[[A2_name]]
    df_lifted$A1_hg38[flip_idx] <- flip_alleles(df_lifted[[A1_name]][flip_idx])
    df_lifted$A2_hg38[flip_idx] <- flip_alleles(df_lifted[[A2_name]][flip_idx])
  } else {
    df_lifted$A1_hg38 <- df_lifted[[A1_name]]
    df_lifted$A2_hg38 <- df_lifted[[A2_name]]
  }

  # Clean up: remove helper columns
  df_lifted$.liftover_row_id <- NULL
  df_lifted$strand_hg38 <- NULL

  # Remove 'chr' prefix from output columns
  df_lifted$CHR_hg38 <- gsub("^chr", "", df_lifted$CHR_hg38)

  # Restore original chr column if we added prefix
  if (!chr_had_prefix) {
    df_lifted[[CHR_name]] <- gsub("^chr", "", df_lifted[[CHR_name]])
  }

  message("Output: ", format(nrow(df_lifted), big.mark = ","), " variants with hg38 coordinates")

  return(df_lifted)
}

#' Flip alleles to complementary bases
#'
#' @param alleles Character vector of alleles (A, T, G, C)
#' @return Character vector of complementary alleles
#' @keywords internal
flip_alleles <- function(alleles) {
  # Create mapping for complement
  mapping <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G",
               "a" = "t", "t" = "a", "g" = "c", "c" = "g")

  # For each allele, flip each base
 sapply(alleles, function(allele) {
    bases <- strsplit(allele, "")[[1]]
    flipped <- mapping[bases]
    # If any base not in mapping, keep original
    flipped[is.na(flipped)] <- bases[is.na(flipped)]
    paste(flipped, collapse = "")
  }, USE.NAMES = FALSE)
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
#' @param tabix_bin Character. Path to tabix binary (default: "tabix")
#' @param bgzip_bin Character. Path to bgzip binary (default: "bgzip")
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
#' @importFrom data.table fread
#' @export
create_dbsnp_file <- function(input_vcf,
                              output_file,
                              compress = TRUE,
                              index = TRUE,
                              tabix_bin = "tabix",
                              bgzip_bin = "bgzip") {

  message("Processing dbSNP VCF file: ", input_vcf)

  # Read VCF file
  vcf_data <- data.table::fread(input_vcf, skip = "#CHROM")
  vcf_data <- as.data.frame(vcf_data)

  # Process multi-allelic variants
  result_list <- lapply(seq_len(nrow(vcf_data)), function(i) {
    row <- vcf_data[i, ]
    alt_alleles <- unlist(strsplit(row$ALT, ","))
    data.frame(
      CHR = row$`#CHROM`,
      POS = row$POS,
      RS_ID = row$ID,
      REF = row$REF,
      ALT = alt_alleles,
      NAME = paste(row$`#CHROM`, row$POS, row$REF, alt_alleles, sep = ":"),
      stringsAsFactors = FALSE
    )
  })

  final_data <- do.call(rbind, result_list)

  # Write output
  temp_file <- paste0(output_file, ".tmp")
  write.table(final_data, temp_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  if (compress) {
    message("Compressing with bgzip...")
    system2(bgzip_bin, temp_file)
    file.rename(paste0(temp_file, ".gz"), output_file)
  } else {
    file.rename(temp_file, output_file)
  }

  if (index && compress) {
    message("Creating tabix index...")
    system2(tabix_bin, c("-s", "1", "-b", "2", "-e", "2", output_file))
  }

  message("Created dbSNP file: ", output_file)
  return(output_file)
}
