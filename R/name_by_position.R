#' Find SNP Names by Genomic Position
#'
#' Matches genomic variants to dbSNP identifiers and standardized names using
#' chromosomal position and allele information. This function queries a tabix-indexed
#' dbSNP VCF file to find matching variants and harmonizes allele coding.
#'
#' @param sumstats A data.frame or data.table containing summary statistics with genomic coordinates.
#'   Can be NULL if input_file is provided.
#' @param input_file Character. Path to input file (supports .gz). Alternative to sumstats parameter
#'   for memory-efficient processing of large files. When provided, chunks are read directly from
#'   disk using bash, avoiding loading the full file into R.
#' @param output_file Character. Path to output file. If provided, results are written to disk
#'   instead of returned as data.frame. Recommended for large files.
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
#'   Set to NULL to disable chunking. Default is 10000 (10K variants) for name_by_position
#'   due to higher RAM usage from dbSNP queries.
#'
#' @return A data.frame with additional columns (or NULL if output_file is provided):
#' \itemize{
#'   \item \strong{rs}: dbSNP rs identifiers
#'   \item \strong{Name_hg38}: Standardized variant names (CHR:POS:REF:ALT format)
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Validates input data and converts to appropriate formats
#'   \item Converts alleles to uppercase (unless keep_lower=TRUE)
#'   \item Uses tabix to query dbSNP file for each chromosome
#'   \item Matches variants by position and allele combination
#'   \item Returns harmonized variant names and rs identifiers
#' }
#'
#' For very large files (>10M variants), use input_file parameter instead of loading
#' data into R first. This enables bash-level chunking where only one chunk is ever
#' in memory at a time.
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
#' # Example with data.frame
#' result <- name_by_position(
#'   sumstats = sumstats,
#'   dbSNP_file = "path/to/dbSNP_all.vcf.gz"
#' )
#'
#' # Example with file input (memory efficient for large files)
#' result <- name_by_position(
#'   input_file = "large_sumstats.tsv.gz",
#'   output_file = "result.tsv.gz",
#'   dbSNP_dir = "/path/to/ensembl_vcf/"
#' )
#' }
#'
#' @importFrom data.table fread fwrite
#' @export
name_by_position <- function(sumstats = NULL,
                             input_file = NULL,
                             output_file = NULL,
                             CHR_name = "CHR_hg38",
                             POS_name = "POS_hg38",
                             A1_name = "A1_hg38",
                             A2_name = "A2_hg38",
                             tabix_bin = "tabix",
                             dbSNP_file = NULL,
                             dbSNP_dir = NULL,
                             keep_lower = FALSE,
                             chunk_size = 10000) {

  message("=== Variant Name Matching ===")

  # Validate input source
  if (is.null(sumstats) && is.null(input_file)) {
    stop("Must provide either sumstats or input_file")
  }
  if (!is.null(sumstats) && !is.null(input_file)) {
    stop("Provide either sumstats or input_file, not both")
  }

  # Validate dependencies
  if (Sys.which(tabix_bin) == "" && !file.exists(tabix_bin)) {
    stop("tabix not found. Install htslib or specify path.")
  }

  # Validate dbSNP source
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

  # ===== FILE-BASED PROCESSING (memory efficient) =====
  if (!is.null(input_file)) {
    if (!file.exists(input_file)) {
      stop("Input file not found: ", input_file)
    }

    message("Using file-based chunking for memory efficiency")

    # Determine if gzipped (.gz or .bgz)
    is_gzipped <- grepl("\\.(gz|bgz)$", input_file)
    cat_cmd <- if (is_gzipped) "zcat" else "cat"

    # Create temp directory for chunks and results
    temp_dir <- tempfile(pattern = "name_by_pos_")
    dir.create(temp_dir)
    chunk_dir <- file.path(temp_dir, "chunks")
    dir.create(chunk_dir)
    on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

    # Get header
    header_cmd <- sprintf("%s '%s' | head -1", cat_cmd, input_file)
    header_line <- system(header_cmd, intern = TRUE)

    # Check for decompression failure (binary data in header)
    if (length(header_line) == 0 || any(charToRaw(header_line) < 0x20 & charToRaw(header_line) != 0x09)) {
      stop("Failed to read header. File may be corrupted or not properly decompressed.\n",
           "Command: ", header_cmd)
    }

    header <- strsplit(header_line, "\t")[[1]]
    if (length(header) < 2) {
      # Try space delimiter as fallback
      header <- strsplit(header_line, "\\s+")[[1]]
      if (length(header) < 2) {
        stop("Header has only ", length(header), " column(s). Could not parse as tab or space-separated.\n",
             "Header line: ", substr(header_line, 1, 100))
      }
      message("Note: Header parsed as space-separated (not tab-separated)")
    }

    # Count lines
    count_cmd <- sprintf("%s '%s' | tail -n +2 | wc -l", cat_cmd, input_file)
    total_lines <- as.integer(system(count_cmd, intern = TRUE))
    message("Input file: ", format(total_lines, big.mark = ","), " variants")

    if (total_lines == 0) {
      warning("Input file has no data rows")
      return(data.frame())
    }

    # Find CHR and POS column indices for sorting
    chr_col <- which(header == CHR_name)
    pos_col <- which(header == POS_name)
    if (length(chr_col) == 0 || length(pos_col) == 0) {
      stop("Could not find CHR column '", CHR_name, "' or POS column '", POS_name, "' in header")
    }

    # Sort by CHR and POS, then split into chunks
    if (is.null(chunk_size)) chunk_size <- total_lines
    message("Sorting and splitting file into chunks...")
    sort_split_cmd <- sprintf(
      "%s '%s' | tail -n +2 | sort -t'\t' -k%d,%dV -k%d,%dn | split -l %d - '%s/chunk_'",
      cat_cmd, input_file, chr_col, chr_col, pos_col, pos_col, chunk_size, chunk_dir
    )
    system(sort_split_cmd)

    chunk_files <- sort(list.files(chunk_dir, full.names = TRUE, pattern = "^chunk_"))
    n_chunks <- length(chunk_files)
    message("Processing ", n_chunks, " chunks of ~", format(chunk_size, big.mark = ","), " variants each")

    total_matched <- 0

    for (i in seq_len(n_chunks)) {
      message("\n--- Chunk ", i, "/", n_chunks, " ---")

      # Read chunk file directly (constant time)
      df_chunk <- data.table::fread(chunk_files[i], header = FALSE, col.names = header)
      df_chunk <- as.data.frame(df_chunk)

      # Add row ID
      df_chunk$.row_id <- seq_len(nrow(df_chunk))

      # Handle allele case
      if (!keep_lower) {
        df_chunk[[A1_name]] <- toupper(df_chunk[[A1_name]])
        df_chunk[[A2_name]] <- toupper(df_chunk[[A2_name]])
      }

      # Process chunk
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
        chunk_result$.row_id <- NULL
        total_matched <- total_matched + nrow(chunk_result)

        # Write chunk result
        result_file <- file.path(temp_dir, paste0("result_", i, ".tsv"))
        data.table::fwrite(chunk_result, result_file, sep = "\t")
      }

      rm(df_chunk, chunk_result)
      gc(verbose = FALSE)
    }

    message("\n--- Summary ---")
    message("Total matched: ", format(total_matched, big.mark = ","), " variants")

    # Combine results
    result_files <- list.files(temp_dir, pattern = "^result_.*\\.tsv$", full.names = TRUE)

    if (length(result_files) == 0) {
      warning("No variants matched. Check chromosome naming and positions.")
      return(data.frame())
    }

    # If output file specified, combine and write
    if (!is.null(output_file)) {
      message("Writing results to: ", output_file)

      # Determine write path (avoid writing to .gz then compressing)
      if (grepl("\\.(gz|bgz)$", output_file)) {
        write_path <- sub("\\.(gz|bgz)$", "", output_file)
        do_compress <- TRUE
      } else {
        write_path <- output_file
        do_compress <- FALSE
      }

      # Write header from first file
      first_result <- data.table::fread(result_files[1], nrows = 0)
      data.table::fwrite(first_result, write_path, sep = "\t")

      # Append all results
      for (rf in result_files) {
        append_cmd <- sprintf("tail -n +2 '%s' >> '%s'", rf, write_path)
        system(append_cmd)
      }

      # Compress if requested
      if (do_compress) {
        system2("gzip", c("-f", write_path))
      }

      message("Output: ", format(total_matched, big.mark = ","), " variants written to file")
      return(invisible(NULL))
    }

    # Otherwise, read and combine all results
    result <- data.table::rbindlist(lapply(result_files, data.table::fread))
    result <- as.data.frame(result)

    message("Output: ", format(nrow(result), big.mark = ","), " variants with rsID and Name")
    return(result)
  }

  # ===== DATA.FRAME-BASED PROCESSING =====
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

      rm(df_chunk, chunk_result)
      gc(verbose = FALSE)
    }

    if (length(chunk_results) == 0) {
      warning("No variants matched. Check chromosome naming and positions.")
      df$.row_id <- NULL
      return(df[0, , drop = FALSE])
    }

    result <- do.call(rbind, chunk_results)

  } else {
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

  # Write to file if requested
  if (!is.null(output_file)) {
    message("Writing results to: ", output_file)
    if (grepl("\\.(gz|bgz)$", output_file)) {
      write_path <- sub("\\.(gz|bgz)$", "", output_file)
      data.table::fwrite(result, write_path, sep = "\t")
      system2("gzip", c("-f", write_path))
    } else {
      data.table::fwrite(result, output_file, sep = "\t")
    }
    message("\nOutput: ", format(nrow(result), big.mark = ","), " variants written to file")
    return(invisible(NULL))
  }

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

    # Expand multi-allelic variants (use data.table for 50x speedup)
    if (any(grepl(",", dbsnp$ALT))) {
      dbsnp <- data.table::as.data.table(dbsnp)
      dbsnp <- dbsnp[, .(ALT = unlist(strsplit(ALT, ","))), by = .(CHR, POS, rsID, REF)]
      dbsnp <- as.data.frame(dbsnp)
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
#' @param sumstats A data.frame or data.table containing summary statistics with genomic coordinates.
#'   Can be NULL if input_file is provided.
#' @param input_file Character. Path to input file (supports .gz). Alternative to sumstats parameter
#'   for memory-efficient processing of large files. When provided, chunks are read directly from
#'   disk using bash, avoiding loading the full file into R.
#' @param output_file Character. Path to output file. If provided, results are written to disk
#'   instead of returned as data.frame. Recommended for large files.
#' @param CHR_name Character. Column name for chromosome
#' @param POS_name Character. Column name for genomic position
#' @param A1_name Character. Column name for first allele
#' @param A2_name Character. Column name for second allele
#' @param liftOver_bin Character. Path to liftOver binary (default: "liftOver")
#' @param liftOver_chain Character. Path to liftOver chain file
#' @param keep_lower Logical. Whether to keep lowercase alleles (default: FALSE)
#' @param chunk_size Integer. Number of variants to process per chunk for memory efficiency.
#'   Set to NULL to disable chunking. Default is 100000 (100K variants).
#'
#' @return A data.frame with lifted coordinates in new columns (or NULL if output_file is provided):
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
#' For very large files (>10M variants), use input_file parameter instead of loading
#' data into R first. This enables bash-level chunking where only one chunk is ever
#' in memory at a time.
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
#' # Basic liftOver from hg19 to hg38 (data.frame mode)
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
#' # File-based mode (memory efficient for large files)
#' genepi_liftover(
#'   input_file = "large_sumstats.tsv.gz",
#'   output_file = "lifted_sumstats.tsv.gz",
#'   CHR_name = "Chr",
#'   POS_name = "Pos_b37",
#'   A1_name = "Allele1",
#'   A2_name = "Allele2",
#'   liftOver_bin = "path/to/liftOver",
#'   liftOver_chain = "path/to/hg19ToHg38.over.chain.gz"
#' )
#' }
#'
#' @importFrom data.table fread fwrite
#' @export
genepi_liftover <- function(sumstats = NULL,
                            input_file = NULL,
                            output_file = NULL,
                            CHR_name,
                            POS_name,
                            A1_name,
                            A2_name,
                            liftOver_bin = "liftOver",
                            liftOver_chain,
                            keep_lower = FALSE,
                            chunk_size = 100000) {

  message("=== Genomic Coordinate Lift Over ===")

  # Validate input source
  if (is.null(sumstats) && is.null(input_file)) {
    stop("Must provide either sumstats or input_file")
  }
  if (!is.null(sumstats) && !is.null(input_file)) {
    stop("Provide either sumstats or input_file, not both")
  }

  # Validate required tools and files
  if (!file.exists(liftOver_bin) && Sys.which(liftOver_bin) == "") {
    stop("liftOver binary not found: ", liftOver_bin,
         "\nDownload from: https://hgdownload.soe.ucsc.edu/admin/exe/")
  }

  if (!file.exists(liftOver_chain)) {
    stop("Chain file not found: ", liftOver_chain,
         "\nDownload from: https://hgdownload.soe.ucsc.edu/goldenPath/")
  }

  # ===== FILE-BASED PROCESSING (memory efficient) =====
  if (!is.null(input_file)) {
    if (!file.exists(input_file)) {
      stop("Input file not found: ", input_file)
    }

    message("Using file-based chunking for memory efficiency")

    # Determine if gzipped (.gz or .bgz)
    is_gzipped <- grepl("\\.(gz|bgz)$", input_file)
    cat_cmd <- if (is_gzipped) "zcat" else "cat"

    # Create temp directory for chunks and results
    temp_dir <- tempfile(pattern = "liftover_")
    dir.create(temp_dir)
    chunk_dir <- file.path(temp_dir, "chunks")
    dir.create(chunk_dir)
    on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

    # Get header
    header_cmd <- sprintf("%s '%s' | head -1", cat_cmd, input_file)
    header_line <- system(header_cmd, intern = TRUE)

    # Check for decompression failure (binary data in header)
    if (length(header_line) == 0 || any(charToRaw(header_line) < 0x20 & charToRaw(header_line) != 0x09)) {
      stop("Failed to read header. File may be corrupted or not properly decompressed.\n",
           "Command: ", header_cmd)
    }

    header <- strsplit(header_line, "\t")[[1]]
    if (length(header) < 2) {
      # Try space delimiter as fallback
      header <- strsplit(header_line, "\\s+")[[1]]
      if (length(header) < 2) {
        stop("Header has only ", length(header), " column(s). Could not parse as tab or space-separated.\n",
             "Header line: ", substr(header_line, 1, 100))
      }
      message("Note: Header parsed as space-separated (not tab-separated)")
    }

    # Count lines
    count_cmd <- sprintf("%s '%s' | tail -n +2 | wc -l", cat_cmd, input_file)
    total_lines <- as.integer(system(count_cmd, intern = TRUE))
    message("Input file: ", format(total_lines, big.mark = ","), " variants")

    if (total_lines == 0) {
      warning("Input file has no data rows")
      return(data.frame())
    }

    # Pre-split file into chunks (one-time cost, true O(1) per chunk read)
    if (is.null(chunk_size)) chunk_size <- total_lines
    message("Splitting file into chunks...")
    split_cmd <- sprintf("%s '%s' | tail -n +2 | split -l %d - '%s/chunk_'",
                         cat_cmd, input_file, chunk_size, chunk_dir)
    system(split_cmd)

    chunk_files <- sort(list.files(chunk_dir, full.names = TRUE, pattern = "^chunk_"))
    n_chunks <- length(chunk_files)
    message("Processing ", n_chunks, " chunks of ~", format(chunk_size, big.mark = ","), " variants each")

    total_lifted <- 0
    total_unmapped <- 0

    for (i in seq_len(n_chunks)) {
      message("\n--- Chunk ", i, "/", n_chunks, " ---")

      # Read chunk file directly (constant time)
      df_chunk <- data.table::fread(chunk_files[i], header = FALSE, col.names = header)
      df_chunk <- as.data.frame(df_chunk)

      # Add row ID
      df_chunk$.liftover_row_id <- seq_len(nrow(df_chunk))

      # Handle chromosome naming (liftOver expects 'chr' prefix)
      chr_had_prefix <- any(grepl("^chr", df_chunk[[CHR_name]], ignore.case = TRUE))
      if (!chr_had_prefix) {
        df_chunk[[CHR_name]] <- paste0("chr", df_chunk[[CHR_name]])
      }

      # Handle allele case
      if (!keep_lower) {
        df_chunk[[A1_name]] <- toupper(df_chunk[[A1_name]])
        df_chunk[[A2_name]] <- toupper(df_chunk[[A2_name]])
      }

      # Process chunk
      chunk_result <- .process_liftover_chunk(
        df_chunk = df_chunk,
        CHR_name = CHR_name,
        POS_name = POS_name,
        A1_name = A1_name,
        A2_name = A2_name,
        liftOver_bin = liftOver_bin,
        liftOver_chain = liftOver_chain
      )

      if (!is.null(chunk_result$result) && nrow(chunk_result$result) > 0) {
        result_df <- chunk_result$result
        result_df$.liftover_row_id <- NULL

        # Remove 'chr' prefix from output columns
        result_df$CHR_hg38 <- gsub("^chr", "", result_df$CHR_hg38)

        # Restore original chr column if we added prefix
        if (!chr_had_prefix) {
          result_df[[CHR_name]] <- gsub("^chr", "", result_df[[CHR_name]])
        }

        total_lifted <- total_lifted + nrow(result_df)

        # Write chunk result
        result_file <- file.path(temp_dir, paste0("result_", i, ".tsv"))
        data.table::fwrite(result_df, result_file, sep = "\t")
      }
      total_unmapped <- total_unmapped + chunk_result$n_unmapped

      rm(df_chunk, chunk_result)
      gc(verbose = FALSE)
    }

    message("\n--- Summary ---")
    message("Lifted: ", format(total_lifted, big.mark = ","), " variants (",
            round(100 * total_lifted / total_lines, 1), "%)")
    if (total_unmapped > 0) {
      message("Unmapped: ", format(total_unmapped, big.mark = ","), " variants (",
              round(100 * total_unmapped / total_lines, 2), "%)")
    }

    # Combine results
    result_files <- list.files(temp_dir, pattern = "^result_.*\\.tsv$", full.names = TRUE)

    if (length(result_files) == 0) {
      warning("No variants lifted. Check input coordinates.")
      return(data.frame())
    }

    # If output file specified, combine and write
    if (!is.null(output_file)) {
      message("Writing results to: ", output_file)

      # Determine write path (avoid writing to .gz then compressing)
      if (grepl("\\.(gz|bgz)$", output_file)) {
        write_path <- sub("\\.(gz|bgz)$", "", output_file)
        do_compress <- TRUE
      } else {
        write_path <- output_file
        do_compress <- FALSE
      }

      # Write header from first file
      first_result <- data.table::fread(result_files[1], nrows = 0)
      data.table::fwrite(first_result, write_path, sep = "\t")

      # Append all results
      for (rf in result_files) {
        append_cmd <- sprintf("tail -n +2 '%s' >> '%s'", rf, write_path)
        system(append_cmd)
      }

      # Compress if requested
      if (do_compress) {
        system2("gzip", c("-f", write_path))
      }

      message("Output: ", format(total_lifted, big.mark = ","), " variants written to file")
      return(invisible(NULL))
    }

    # Otherwise, read and combine all results
    result <- data.table::rbindlist(lapply(result_files, data.table::fread))
    result <- as.data.frame(result)

    message("Output: ", format(nrow(result), big.mark = ","), " variants with hg38 coordinates")
    return(result)
  }

  # ===== DATA.FRAME-BASED PROCESSING =====
  # Work on a copy as data.frame
  df <- as.data.frame(sumstats)
  n_input <- nrow(df)
  message("Input: ", format(n_input, big.mark = ","), " variants")

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

  # Add row index for tracking
  df$.liftover_row_id <- seq_len(nrow(df))

  # Determine if chunking is needed
  use_chunking <- !is.null(chunk_size) && n_input > chunk_size

  if (use_chunking) {
    # Process in chunks
    n_chunks <- ceiling(n_input / chunk_size)
    message("Processing in ", n_chunks, " chunks of ~", format(chunk_size, big.mark = ","), " variants each")

    chunk_results <- list()
    total_unmapped <- 0

    for (i in seq_len(n_chunks)) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, n_input)

      message("\n--- Chunk ", i, "/", n_chunks, " (rows ", format(start_idx, big.mark = ","),
              "-", format(end_idx, big.mark = ","), ") ---")

      df_chunk <- df[start_idx:end_idx, , drop = FALSE]

      chunk_result <- .process_liftover_chunk(
        df_chunk = df_chunk,
        CHR_name = CHR_name,
        POS_name = POS_name,
        A1_name = A1_name,
        A2_name = A2_name,
        liftOver_bin = liftOver_bin,
        liftOver_chain = liftOver_chain
      )

      if (!is.null(chunk_result$result) && nrow(chunk_result$result) > 0) {
        chunk_results[[i]] <- chunk_result$result
      }
      total_unmapped <- total_unmapped + chunk_result$n_unmapped

      # Clean up
      rm(df_chunk, chunk_result)
      gc(verbose = FALSE)
    }

    # Combine results
    if (length(chunk_results) == 0) {
      warning("No variants lifted. Check input coordinates.")
      df$.liftover_row_id <- NULL
      return(df[0, , drop = FALSE])
    }

    df_lifted <- do.call(rbind, chunk_results)

    # Summary statistics
    message("\n--- Summary ---")
    message("Lifted: ", format(nrow(df_lifted), big.mark = ","), " variants (",
            round(100 * nrow(df_lifted) / n_input, 1), "%)")
    if (total_unmapped > 0) {
      message("Unmapped: ", format(total_unmapped, big.mark = ","), " variants (",
              round(100 * total_unmapped / n_input, 2), "%)")
    }

  } else {
    # Process all at once
    result <- .process_liftover_chunk(
      df_chunk = df,
      CHR_name = CHR_name,
      POS_name = POS_name,
      A1_name = A1_name,
      A2_name = A2_name,
      liftOver_bin = liftOver_bin,
      liftOver_chain = liftOver_chain
    )

    if (is.null(result$result) || nrow(result$result) == 0) {
      warning("No variants lifted. Check input coordinates.")
      df$.liftover_row_id <- NULL
      return(df[0, , drop = FALSE])
    }

    df_lifted <- result$result
  }

  rownames(df_lifted) <- NULL

  # Remove helper columns
  df_lifted$.liftover_row_id <- NULL

  # Remove 'chr' prefix from output columns
  df_lifted$CHR_hg38 <- gsub("^chr", "", df_lifted$CHR_hg38)

  # Restore original chr column if we added prefix
  if (!chr_had_prefix) {
    df_lifted[[CHR_name]] <- gsub("^chr", "", df_lifted[[CHR_name]])
  }

  # Write to file if requested
  if (!is.null(output_file)) {
    message("Writing results to: ", output_file)
    if (grepl("\\.(gz|bgz)$", output_file)) {
      write_path <- sub("\\.(gz|bgz)$", "", output_file)
      data.table::fwrite(df_lifted, write_path, sep = "\t")
      system2("gzip", c("-f", write_path))
    } else {
      data.table::fwrite(df_lifted, output_file, sep = "\t")
    }
    message("\nOutput: ", format(nrow(df_lifted), big.mark = ","), " variants written to file")
    return(invisible(NULL))
  }

  message("\nOutput: ", format(nrow(df_lifted), big.mark = ","), " variants with hg38 coordinates")

  return(df_lifted)
}

#' Process a chunk of variants through liftOver
#'
#' Internal function that processes a subset of variants through UCSC liftOver.
#'
#' @param df_chunk Data frame chunk to process
#' @param CHR_name Column name for chromosome
#' @param POS_name Column name for position
#' @param A1_name Column name for allele 1
#' @param A2_name Column name for allele 2
#' @param liftOver_bin Path to liftOver binary
#' @param liftOver_chain Path to chain file
#'
#' @return List with 'result' (data frame) and 'n_unmapped' (integer)
#' @keywords internal
.process_liftover_chunk <- function(df_chunk, CHR_name, POS_name, A1_name, A2_name,
                                     liftOver_bin, liftOver_chain) {

  n_chunk <- nrow(df_chunk)

  # Create temp files
  bed_file <- tempfile(pattern = "liftover_input_", fileext = ".bed")
  output_file <- tempfile(pattern = "liftover_output_", fileext = ".bed")
  unmapped_file <- tempfile(pattern = "liftover_unmapped_", fileext = ".bed")

  on.exit({
    unlink(c(bed_file, output_file, unmapped_file), force = TRUE)
  }, add = TRUE)

  # Write BED file
  bed_data <- data.frame(
    chr = df_chunk[[CHR_name]],
    start = df_chunk[[POS_name]],
    end = df_chunk[[POS_name]],
    name = df_chunk$.liftover_row_id,
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
  if (!file.exists(output_file) || file.size(output_file) == 0) {
    message("No variants lifted in this chunk")
    return(list(result = NULL, n_unmapped = n_chunk))
  }

  lifted <- data.table::fread(output_file, header = FALSE,
                               col.names = c("CHR_hg38", "POS_hg38", "end", ".liftover_row_id", "score", "strand_hg38"))
  lifted <- as.data.frame(lifted)
  lifted$end <- NULL
  lifted$score <- NULL

  # Count unmapped
  n_unmapped <- 0
  if (file.exists(unmapped_file) && file.size(unmapped_file) > 0) {
    unmapped_lines <- readLines(unmapped_file, warn = FALSE)
    unmapped_lines <- unmapped_lines[!grepl("^#", unmapped_lines)]
    n_unmapped <- length(unmapped_lines)
  }

  # Report
  n_lifted <- nrow(lifted)
  message("Lifted: ", format(n_lifted, big.mark = ","), "/", format(n_chunk, big.mark = ","),
          " variants (", round(100 * n_lifted / n_chunk, 1), "%)")

  # Merge
  df_lifted <- merge(df_chunk, lifted, by = ".liftover_row_id", all = FALSE)

  # Handle strand flipping
  flip_idx <- df_lifted$strand_hg38 == "-"
  n_flipped <- sum(flip_idx)

  if (n_flipped > 0) {
    message("Flipping ", format(n_flipped, big.mark = ","), " variants on negative strand")
    df_lifted$A1_hg38 <- df_lifted[[A1_name]]
    df_lifted$A2_hg38 <- df_lifted[[A2_name]]
    df_lifted$A1_hg38[flip_idx] <- flip_alleles(df_lifted[[A1_name]][flip_idx])
    df_lifted$A2_hg38[flip_idx] <- flip_alleles(df_lifted[[A2_name]][flip_idx])
  } else {
    df_lifted$A1_hg38 <- df_lifted[[A1_name]]
    df_lifted$A2_hg38 <- df_lifted[[A2_name]]
  }

  # Remove strand column
  df_lifted$strand_hg38 <- NULL

  return(list(result = df_lifted, n_unmapped = n_unmapped))
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
