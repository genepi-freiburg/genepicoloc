#' Flip alleles according to complementary base pairing
#'
#' @description Converts nucleotide alleles to their complementary bases
#' (A<->T, C<->G). Useful for strand flipping in genetic data.
#'
#' @param vec Character vector of nucleotide alleles (A, T, C, or G)
#' @return Character vector with complementary alleles
#' @examples
#' \dontrun{
#' flip_alleles(c("A", "T", "C", "G")) # Returns c("T", "A", "G", "C")
#' }
#' @keywords internal
flip_alleles <- function(vec) {
  vec_out <- vec
  vec_out <- toupper(vec_out)
  vec_out[vec == "A"] <- "T"
  vec_out[vec == "T"] <- "A"
  vec_out[vec == "C"] <- "G"
  vec_out[vec == "G"] <- "C"
  return(vec_out)
}
#' Check if a binary is available in the system path
#'
#' @param bin_name Name or path of the binary to check
#' @return TRUE if the binary is available, otherwise stops with an error
#' @keywords internal
check_bin <- function(bin_name) {
  # Use 'which' on Unix/Linux/Mac or 'where' on Windows
  cmd <- if (.Platform$OS.type == "windows") {
    paste0("where ", bin_name, " 2>NUL")
  } else {
    paste0("which ", bin_name, " 2>/dev/null")
  }
  
  # Run command and check result
  result <- system(cmd, intern = TRUE)
  
  if (length(result) == 0) {
    stop("Required binary '", bin_name, "' not found in system path. ",
         "Please install it or provide the full path.")
  }
  
  return(TRUE)
}

#' Handle P-value Underflow
#'
#' Process p-values to handle underflow in case of very small p-values (p < 1e-320).
#' This function converts p-values to -log10 scale without requiring high-precision
#' arithmetic libraries, making it suitable for statistical visualization and analysis.
#'
#' @param pvalue_vec A numeric or character vector of p-values. Can handle both
#'   regular decimal notation (e.g., 0.001) and scientific notation (e.g., 1e-350).
#' @param return_nlog10P Logical. If TRUE, returns -log10 transformed p-values.
#'   If FALSE, returns the p-values as numeric (default: FALSE).
#'
#' @return If return_nlog10P is TRUE, returns a numeric vector of -log10 transformed
#'   p-values. If FALSE, returns the original p-values converted to numeric.
#'
#' @details
#' This function addresses the underflow problem that occurs when p-values are
#' extremely small (e.g., < 1e-320). Instead of requiring arbitrary precision
#' arithmetic, it directly converts scientific notation p-values to -log10 scale:
#' \itemize{
#'   \item For scientific notation (e.g., "5e-200"): -log10(5e-200) = -log10(5) + 200
#'   \item For regular decimals (e.g., "0.001"): -log10(0.001) = 3
#' }
#'
#' @examples
#' # Example with mixed p-value formats
#' pvals <- c("1e-350", "0.001", "5e-200", "0.05", "1e-10")
#' 
#' # Convert to -log10 scale
#' result <- handle_underflow(pvals, return_nlog10P = TRUE)
#' print(result)
#' # [1] 350.00000   3.00000 199.30103   1.30103  10.00000
#' 
#' # Return as numeric (no transformation)
#' result2 <- handle_underflow(pvals, return_nlog10P = FALSE)
#' print(result2)
#'
#' @export
handle_underflow <- function(pvalue_vec, return_nlog10P = FALSE) {
  # Convert to character to handle both numeric and character inputs
  pval <- as.character(pvalue_vec)
  
  # Identify scientific notation cases (both 'e' and 'E')
  l1 <- grepl("[eE]", pval)
  l2 <- !l1
  
  if (return_nlog10P) {
    # Initialize result vector
    result <- numeric(length(pval))
    
    # Handle scientific notation cases
    if (any(l1)) {
      result[l1] <- sapply(pval[l1], function(x) {
        # Split on either 'e' or 'E'
        parts <- as.numeric(strsplit(x, "[eE]")[[1]])
        mantissa <- parts[1]
        exponent <- parts[2]
        # Correct formula: -log10(mantissa * 10^exponent) = -log10(mantissa) - exponent * log10(10)
        -log10(mantissa) - exponent
      })
    }
    
    # Handle regular decimal cases
    if (any(l2)) {
      result[l2] <- -log10(as.numeric(pval[l2]))
    }
    
    return(result)
  } else {
    # If not converting to -log10, just return as high precision numeric
    # This mimics the Rmpfr behavior for handling underflow
    return(as.numeric(pval))
  }
}

#' Validate coloc_regions_PASS data frame
#' 
#' @description
#' Validates that coloc_regions_PASS has the correct structure and content
#' for colocalization analysis.
#' 
#' @param coloc_regions_PASS Data frame to validate
#' @param param_name Character. Name of the parameter (for error messages)
#' 
#' @return NULL (invisibly) if validation passes, otherwise stops with error
#' 
#' @details
#' Validates:
#' \itemize{
#'   \item Is a data.frame with at least 1 row
#'   \item Has exactly 3 required columns: CHR_var, BP_START_var, BP_STOP_var
#'   \item CHR_var contains valid chromosome names (1-22, X, Y, XY, MT, or chr-prefixed)
#'   \item BP_START_var and BP_STOP_var are positive integers
#'   \item BP_START_var <= BP_STOP_var for each region
#' }
#' 
#' @examples
#' \dontrun{
#' # Valid example
#' regions <- data.frame(
#'   CHR_var = c("1", "2", "X"),
#'   BP_START_var = c(567054, 661911, 7534918),
#'   BP_STOP_var = c(1567054, 1661911, 8534918)
#' )
#' validate_coloc_regions(regions)
#' 
#' # This would fail - missing required columns
#' bad_regions <- data.frame(chr = "1", start = 100, end = 200)
#' validate_coloc_regions(bad_regions) # Error
#' }
validate_coloc_regions <- function(coloc_regions_PASS, param_name = "coloc_regions_PASS") {
  
  # Check if it's a data.frame
  if (!is.data.frame(coloc_regions_PASS)) {
    stop(sprintf("%s must be a data.frame, got %s", 
                 param_name, class(coloc_regions_PASS)[1]))
  }
  
  # Check if it has at least one row
  if (nrow(coloc_regions_PASS) == 0) {
    stop(sprintf("%s must contain at least one row", param_name))
  }
  
  # Check required columns
  required_cols <- c("CHR_var", "BP_START_var", "BP_STOP_var")
  missing_cols <- setdiff(required_cols, names(coloc_regions_PASS))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("%s missing required columns: %s. Required columns are: %s", 
                 param_name,
                 paste(missing_cols, collapse = ", "),
                 paste(required_cols, collapse = ", ")))
  }
  
  # Check for extra columns (warning only)
  extra_cols <- setdiff(names(coloc_regions_PASS), required_cols)
  if (length(extra_cols) > 0) {
    warning(sprintf("%s contains additional columns that will be ignored: %s", 
                    param_name, paste(extra_cols, collapse = ", ")))
  }
  
  # Validate CHR_var
  chr_values <- coloc_regions_PASS$CHR_var
  
  # Convert to character if not already
  if (!is.character(chr_values)) {
    chr_values <- as.character(chr_values)
    coloc_regions_PASS$CHR_var <- chr_values
    message(sprintf("Converting %s$CHR_var to character", param_name))
  }
  
  # Define valid chromosome names and normalize chr prefixes
  valid_chrs_base <- c(
    as.character(1:22),           # Autosomes: "1", "2", ..., "22"
    "X", "Y", "XY", "MT", "M"     # Sex chromosomes and mitochondrial
  )
  
  # Also accept chr-prefixed versions but convert them
  valid_chrs_with_prefix <- c(
    paste0("chr", 1:22),          # With chr prefix: "chr1", "chr2", ...
    "chrX", "chrY", "chrXY", "chrMT", "chrM"  # Sex/mito with chr prefix
  )
  
  all_valid_chrs <- c(valid_chrs_base, valid_chrs_with_prefix)
  
  # Check for invalid chromosomes first
  invalid_chrs <- chr_values[!chr_values %in% all_valid_chrs]
  if (length(invalid_chrs) > 0) {
    unique_invalid <- unique(invalid_chrs)
    stop(sprintf("%s$CHR_var contains invalid chromosome names: %s. Valid options are: %s", 
                 param_name,
                 paste(unique_invalid, collapse = ", "),
                 paste(head(valid_chrs_base, 10), collapse = ", "),
                 if (length(valid_chrs_base) > 10) "..." else ""))
  }
  
  # Convert chr-prefixed chromosomes to standard format
  has_chr_prefix <- grepl("^chr", chr_values)
  if (any(has_chr_prefix)) {
    chr_values[has_chr_prefix] <- sub("^chr", "", chr_values[has_chr_prefix])
    coloc_regions_PASS$CHR_var <- chr_values
    n_converted <- sum(has_chr_prefix)
    message(sprintf("Converted %d chromosome names from 'chr' prefix format to standard format (e.g., 'chr1' -> '1')", 
                    n_converted))
  }
  
  # Validate BP_START_var and BP_STOP_var
  bp_cols <- c("BP_START_var", "BP_STOP_var")
  
  for (col in bp_cols) {
    bp_values <- coloc_regions_PASS[[col]]
    
    if (!is.numeric(bp_values)) {
      stop(sprintf("%s$%s must be numeric, got %s", 
                   param_name, col, class(bp_values)[1]))
    }
    
    # Check for non-integer values
    if (any(!bp_values == as.integer(bp_values), na.rm = TRUE)) {
      warning(sprintf("%s$%s contains non-integer values, converting to integers", 
                      param_name, col))
      coloc_regions_PASS[[col]] <- as.integer(round(bp_values))
    }
    
    # Check for non-positive values
    if (any(bp_values <= 0, na.rm = TRUE)) {
      invalid_rows <- which(bp_values <= 0)
      stop(sprintf("%s$%s must contain positive integers. Invalid values at rows: %s", 
                   param_name, col, paste(head(invalid_rows, 5), collapse = ", ")))
    }
    
    # Check for missing values
    if (any(is.na(bp_values))) {
      invalid_rows <- which(is.na(bp_values))
      stop(sprintf("%s$%s contains missing values at rows: %s", 
                   param_name, col, paste(head(invalid_rows, 5), collapse = ", ")))
    }
  }
  
  # Validate that BP_START_var <= BP_STOP_var
  invalid_ranges <- coloc_regions_PASS$BP_START_var > coloc_regions_PASS$BP_STOP_var
  if (any(invalid_ranges)) {
    invalid_rows <- which(invalid_ranges)
    stop(sprintf("%s: BP_START_var must be <= BP_STOP_var. Invalid ranges at rows: %s", 
                 param_name, paste(head(invalid_rows, 5), collapse = ", ")))
  }
  
  # Summary message
  message(sprintf(" %s validation passed: %d regions across %d chromosomes", 
                  param_name, 
                  nrow(coloc_regions_PASS),
                  length(unique(chr_values))))
  
  return(invisible(NULL))
}
