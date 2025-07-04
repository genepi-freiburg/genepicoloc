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
