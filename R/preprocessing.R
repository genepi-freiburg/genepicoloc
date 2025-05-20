# preprocessing ----

#' Extract GWAS summary statistics from tabix-indexed file
#'
#' @param sumstats_file Path to tabix-indexed summary statistics file
#' @param coloc_regions_PASS Data frame with regions to extract
#' @param verbose Whether to print progress messages
#' @param test_mode If TRUE, only use the first region
#' @param file_remove Whether to remove temporary files
#' @return A sumstats object
sumstats_tabix <- function(sumstats_file, coloc_regions_PASS,
                           verbose = FALSE, test_mode = FALSE, file_remove = TRUE) {
  # Check required columns
  req_cols <- c("CHR_var", "BP_START_var", "BP_STOP_var")
  if (!all(req_cols %in% colnames(coloc_regions_PASS))) {
    stop("coloc_regions_PASS must contain: ", paste(req_cols, collapse = ", "))
  }
  
  # Create temporary file with regions
  file_regions <- tempfile()
  if (test_mode) coloc_regions_PASS <- coloc_regions_PASS[1, ]
  
  # Sort and write regions
  sorted <- coloc_regions_PASS[order(coloc_regions_PASS$CHR_var, 
                                     coloc_regions_PASS$BP_START_var), ]
  write.table(sorted, file_regions, sep = "\t", row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  
  # Run tabix
  tabix_cmd <- paste0("tabix -h ", sumstats_file, " -R ", file_regions)
  if (verbose) message("Running tabix: ", tabix_cmd)
  text_in <- suppressWarnings(system(tabix_cmd, intern = TRUE, ignore.stderr = TRUE))
  
  # Process results
  if (is.null(attr(text_in, "status"))) {
    sumstats <- data.table::fread(text = text_in, showProgress = FALSE)
    tabix_attr <- if (nrow(sumstats) == 0) "tabix_ok_no_data" else "tabix_ok"
  } else {
    tabix_attr <- "tabix_failed"
    sumstats <- data.table::data.table()
  }
  
  # Set attributes
  attr(sumstats, "tabix") <- tabix_attr
  attr(sumstats, "sumstats_file") <- sumstats_file
  attr(sumstats, "QC") <- "not_checked"
  class(sumstats) <- unique(c("sumstats", class(sumstats)))
  
  # Clean up
  if (file_remove) file.remove(file_regions)
  
  return(sumstats)
}


#' Calculate and set max_nlog10P attribute for a sumstats object
#'
#' @param sumstats A sumstats object
#' @return The sumstats object with updated max_nlog10P attribute
sumstats_nlog10P <- function(sumstats) {
  if (!inherits(sumstats, "sumstats")) 
    stop("Expected a sumstats object")
  
  # Calculate max_nlog10P
  max_nlog10P <- if (nrow(sumstats) > 0 && "nlog10P" %in% names(sumstats)) {
    max(sumstats[["nlog10P"]], na.rm = TRUE)
  } else {
    NA
  }
  
  attr(sumstats, "max_nlog10P") <- max_nlog10P
  return(sumstats)
}

#' Add core attributes to a sumstats object
#'
#' @param sumstats A sumstats object
#' @param sumstats_function Name of the function that created the data
#' @param sumstats_type Type of summary statistics ('quantitative' or 'binary')
#' @param sumstats_sdY Standard deviation for quantitative trait
#' @return sumstats object with additional attributes
sumstats_attr <- function(sumstats, sumstats_function,
                          sumstats_type, sumstats_sdY) {
  if (!inherits(sumstats, "sumstats")) 
    stop("Expected a sumstats object")
  
  # Add attributes
  attr(sumstats, "sumstats_function") <- sumstats_function
  attr(sumstats, "sumstats_type") <- sumstats_type
  attr(sumstats, "sumstats_sdY") <- sumstats_sdY
  
  return(sumstats)
}
#' Check and clean a sumstats object
#'
#' @param sumstats A sumstats object
#' @param verbose Whether to print progress messages
#' @return Validated and cleaned sumstats object
sumstats_check <- function(sumstats, verbose = FALSE) {
  # Check object type and required attributes
  if (!inherits(sumstats, "sumstats")) 
    stop("Expected a sumstats object")
  
  # Check required attributes
  req_attrs <- c("sumstats_type", "sumstats_sdY", "sumstats_file", 
                 "max_nlog10P", "tabix")
  missing <- req_attrs[sapply(req_attrs, function(a) is.null(attr(sumstats, a)))]
  if (length(missing) > 0)
    stop("Missing attributes: ", paste(missing, collapse = ", "))
  
  # Check column names
  if (!all(colnames(sumstats) == get_cols_to()))
    stop("Column names don't match expected format")
  
  # Start QC tracking
  QC <- ""
  
  # Remove duplicates
  dups <- duplicated(sumstats$Name)
  if (any(dups)) {
    if (verbose) message("Removing ", sum(dups), " duplicated variants")
    QC <- paste0(QC, "_dupvar")
    sumstats <- sumstats[!dups, ]
  }
  
  # Check numeric columns
  for (col in c("BETA", "SE")) {
    # Remove infinite values
    inf_vals <- is.infinite(sumstats[[col]])
    if (any(inf_vals)) {
      if (verbose) message("Removing ", sum(inf_vals), " infinite values in ", col)
      QC <- paste0(QC, "_INFin", col)
      sumstats <- sumstats[!inf_vals, ]
    }
    
    # Remove NAs
    na_vals <- is.na(sumstats[[col]])
    if (any(na_vals)) {
      if (verbose) message("Removing ", sum(na_vals), " NA values in ", col)
      QC <- paste0(QC, "_NAin", col)
      sumstats <- sumstats[!na_vals, ]
    }
    
    # Remove zeros
    zero_vals <- sumstats[[col]] == 0
    if (any(zero_vals)) {
      if (verbose) message("Removing ", sum(zero_vals), " zero values in ", col)
      QC <- paste0(QC, "_0in", col)
      sumstats <- sumstats[!zero_vals, ]
    }
  }
  
  # Set QC status
  attr(sumstats, "QC") <- if (QC == "") "ok" else substring(QC, 2)
  
  if (verbose) message("QC completed")
  return(sumstats)
}

#' Get standard column names for GWAS summary statistics
#'
#' @description Returns the standard column names used for summary statistics
#' in the colocalization analysis. All summary statistics objects must have 
#' these exact column names in this exact order.
#'
#' @return Character vector of standard column names
#' @export
get_cols_to <- function() {
  c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "nlog10P", "AF", "N") 
}

#' Map non-standard column names to standard format
#'
#' @description Converts a GWAS summary statistics data frame with non-standard
#' column names to the standardized format required for colocalization analysis.
#'
#' @param sumstats Data frame or data.table containing GWAS summary statistics
#' @param Name Column name for variant identifier
#' @param rsID Column name for rs ID
#' @param CHR Column name for chromosome
#' @param POS Column name for position
#' @param A1 Column name for effect allele
#' @param A2 Column name for other allele
#' @param BETA Column name for effect size
#' @param SE Column name for standard error
#' @param nlog10P Column name for -log10(p-value)
#' @param AF Column name for allele frequency
#' @param N Column name for sample size
#' @return Data frame with standardized column names and order
#' @export
#' @examples
#' \dontrun{
#' # Convert column names in a GWAS summary statistics data frame
#' standardized_gwas <- match_cols(
#'   my_gwas,
#'   Name = "variant_id",
#'   rsID = "rs_id",
#'   CHR = "chromosome", 
#'   POS = "position",
#'   A1 = "effect_allele",
#'   A2 = "other_allele",
#'   BETA = "effect_size",
#'   SE = "standard_error",
#'   nlog10P = "neg_log10_pvalue",
#'   AF = "allele_freq",
#'   N = "sample_size"
#' )
#' }
match_cols <- function(sumstats, Name, rsID, CHR, POS, A1, A2, 
                       BETA, SE, nlog10P, AF, N) {
  if (!is.data.table(sumstats)) stop("data.table object expected")
  
  # Check if all specified columns exist in the data
  input_cols <- c(Name, rsID, CHR, POS, A1, A2, BETA, SE, nlog10P, AF, N)
  if (!all(input_cols %in% colnames(sumstats))) {
    missing_cols <- input_cols[!input_cols %in% colnames(sumstats)]
    stop("Missing columns in input data: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  # Create mapping between input and standard column names
  col_mapping <- c(
    Name = Name,
    rsID = rsID,
    CHR = CHR,
    POS = POS,
    A1 = A1,
    A2 = A2,
    BETA = BETA,
    SE = SE,
    nlog10P = nlog10P,
    AF = AF,
    N = N
  )
  
  # Rename columns using the mapping
  for (std_name in names(col_mapping)) {
    input_name <- col_mapping[std_name]
    colnames(sumstats)[colnames(sumstats) == input_name] <- std_name
  }
  
  # Select and reorder columns to match standard format
  std_cols <- get_cols_to()
  sumstats <- sumstats[, ..std_cols]
  
  return(sumstats)
}


