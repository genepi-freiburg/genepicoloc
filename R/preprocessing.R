# preprocessing ----
retrieve_sumstats_wrapper <- function(sumstats_function,
                                      sumstats_file,
                                      coloc_regions_PASS) {
  sumstats <- do.call(what = sumstats_function, 
                      args = list(sumstats_file = sumstats_file,
                                  coloc_regions_PASS = coloc_regions_PASS))
  attr(sumstats, "sumstats_function") <- sumstats_function
  return(sumstats)
}


#' Extract GWAS summary statistics from tabix-indexed file
#'
#' @param sumstats_file Path to tabix-indexed summary statistics file
#' @param coloc_regions_PASS Data frame with regions to extract
#' @param verbose Whether to print progress messages
#' @param test_mode If TRUE, only use the first region
#' @param file_remove Whether to remove temporary files
#' @return A sumstats object
retrieve_sumstats_tabix <- function(sumstats_file, coloc_regions_PASS,
                                    verbose = FALSE, test_mode = FALSE,
                                    file_remove = TRUE) {
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
  
  # Run tabix with fread
  tabix_cmd <- paste0("tabix -h ", shQuote(sumstats_file), " -R ",
                      shQuote(file_regions), " 2>/dev/null")
  if (verbose) message("Running tabix: ", tabix_cmd)
  
  sumstats <- suppressWarnings(
    data.table::fread(cmd = tabix_cmd, showProgress = FALSE)
  )
  
  # Determine status
  if (is.null(sumstats) || (nrow(sumstats) == 0 && ncol(sumstats) == 0)) {
    tabix_attr <- "tabix_failed"
    sumstats <- data.table::data.table()
  } else {
    tabix_attr <- if (nrow(sumstats) == 0) "tabix_ok_no_data" else "tabix_ok"
  }
  
  # Set attributes
  attr(sumstats, "tabix") <- tabix_attr
  attr(sumstats, "sumstats_file") <- sumstats_file
  # attr(sumstats, "sumstats_str") <- sumstats_str
  class(sumstats) <- unique(c("sumstats", class(sumstats)))
  
  # Clean up
  if (file_remove) file.remove(file_regions)
  
  return(sumstats)
}

#' Calculate and set max_nlog10P attribute for a sumstats object
#'
#' @param sumstats A sumstats object
#' @return The sumstats object with updated max_nlog10P attribute
set_max_nlog10P <- function(sumstats) {
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


#' Perform QC on a sumstats object
#'
#' @param sumstats A validated sumstats object
#' @param verbose Whether to print progress messages
#' @return QC'd sumstats object with QC attribute set
perform_sumstats_qc <- function(sumstats, verbose = FALSE) {
  # Skip QC if tabix failed
  if (attr(sumstats, "tabix") == "tabix_failed") {
    attr(sumstats, "QC") <- "skipped_tabix_failed"
    if (verbose) message("Skipping QC - tabix failed")
    return(sumstats)
  }
  
  QC <- ""
  
  # Remove duplicated variants
  dups <- duplicated(sumstats$Name)
  if (any(dups)) {
    if (verbose) message("Removing ", sum(dups), " duplicated variants")
    QC <- paste0(QC, "_dupvar")
    sumstats <- sumstats[!dups, ]
  }
  
  # Clean numeric columns
  for (col in c("BETA", "SE")) {
    inf_vals <- is.infinite(sumstats[[col]])
    if (any(inf_vals)) {
      if (verbose) message("Removing ", sum(inf_vals), " infinite values in ", col)
      QC <- paste0(QC, "_INFin", col)
      sumstats <- sumstats[!inf_vals, ]
    }
    
    na_vals <- is.na(sumstats[[col]])
    if (any(na_vals)) {
      if (verbose) message("Removing ", sum(na_vals), " NA values in ", col)
      QC <- paste0(QC, "_NAin", col)
      sumstats <- sumstats[!na_vals, ]
    }
    
    zero_vals <- sumstats[[col]] == 0
    if (any(zero_vals)) {
      if (verbose) message("Removing ", sum(zero_vals), " zero values in ", col)
      QC <- paste0(QC, "_0in", col)
      sumstats <- sumstats[!zero_vals, ]
    }
  }
  
  attr(sumstats, "QC") <- if (QC == "") "ok" else substring(QC, 2)
  
  if (verbose) message("QC completed")
  return(sumstats)
}

#' Check required attributes of a sumstats object
#'
#' @param sumstats A sumstats object
#' @param sumstats_function Name of the function that created the data
#' @param sumstats_type Type of summary statistics ('quantitative' or 'binary')
#' @param sumstats_sdY Standard deviation for quantitative trait
#' @return The input object if valid; otherwise, stops with an error
check_sumstats_attributes <- function(sumstats, sumstats_function,
                                      sumstats_type, sumstats_sdY) {
  # Check object type
  if (!inherits(sumstats, "sumstats")) 
    stop("Expected a sumstats object")
  
  # Add attributes if provided
    attr(sumstats, "sumstats_type") <- sumstats_type
  if (!missing(sumstats_sdY)) 
    attr(sumstats, "sumstats_sdY") <- sumstats_sdY
  
  # Required attributes
  req_attrs <- c("sumstats_type", "sumstats_sdY", "sumstats_file", 
                 "sumstats_function", "tabix")
  missing <- req_attrs[sapply(req_attrs, function(a) is.null(attr(sumstats, a)))]
  if (length(missing) > 0)
    stop("Missing attributes: ", paste(missing, collapse = ", "))
  
  # Check expected columns
  if (!all(colnames(sumstats) == get_cols_to()))
    stop("Column names don't match expected format")
  
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
                       BETA, SE, nlog10P, AF, N, Phenotype=NULL) {
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
  if (!is.null(Phenotype)) col_mapping <- c(col_mapping, Phenotype = Phenotype)
  
  # Rename columns using the mapping
  for (std_name in names(col_mapping)) {
    input_name <- col_mapping[std_name]
    colnames(sumstats)[colnames(sumstats) == input_name] <- std_name
  }
  
  # Select and reorder columns to match standard format
  std_cols <- get_cols_to()
  if (!is.null(Phenotype)) std_cols <- c(std_cols, Phenotype = Phenotype)
  sumstats <- sumstats[, ..std_cols]
  
  return(sumstats)
}


