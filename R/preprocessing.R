# Group 1: Data Retrieval and Processing Functions ----

#' Retrieve raw summary statistics using specified function
#' 
#' @description
#' A wrapper function that calls the appropriate retrieval function for
#' different summary statistics formats. This allows flexibility in handling
#' various data sources (tabix-indexed files, API calls, etc.).
#' 
#' @param sumstats_function Character string. Name of the function to use for
#'   retrieving summary statistics (default: "retrieve_sumstats_tabix").
#' @param sumstats_file Character string. Path to the summary statistics file
#'   or identifier for API-based retrieval.
#' @param coloc_regions_PASS A data.frame containing genomic regions to extract.
#'   Must contain columns CHR_var, BP_START_var, and BP_STOP_var.
#' 
#' @return A sumstats object with an additional "sumstats_function" attribute
#'   indicating which function was used for retrieval.
#' 
#' @details
#' This function provides a unified interface for different data retrieval
#' methods. Common retrieval functions include:
#' \itemize{
#'   \item retrieve_sumstats_tabix: For tabix-indexed files
#'   \item query_eQTL_Catalogue: For eQTL Catalogue API
#'   \item tabix_GTEXv8: For GTEx v8 format
#'   \item tabix_UKB_PPP_EUR: For UK Biobank format
#' }
#' 
#' @seealso 
#' \code{\link{retrieve_sumstats_tabix}} for the default retrieval method
#' 
#' @export
retrieve_sumstats_raw <- function(sumstats_function = "retrieve_sumstats_tabix",
                                  sumstats_file,
                                  coloc_regions_PASS) {
  
  # Input validation
  if (!is.character(sumstats_function) || length(sumstats_function) != 1) {
    stop("sumstats_function must be a single character string")
  }
  
  if (!is.data.frame(coloc_regions_PASS)) {
    stop("coloc_regions_PASS must be a data.frame")
  }
  
  # Check if the function exists
  if (!exists(sumstats_function, mode = "function")) {
    stop("Function '", sumstats_function, "' not found")
  }
  
  # Call the specified function
  sumstats <- do.call(
    what = sumstats_function, 
    args = list(
      sumstats_file = sumstats_file,
      coloc_regions_PASS = coloc_regions_PASS
    )
  )
  
  # Add attribute to track which function was used
  attr(sumstats, "sumstats_function") <- sumstats_function
  
  return(sumstats)
}


#' Extract GWAS summary statistics from tabix-indexed file
#'
#' @description
#' Retrieves summary statistics from a tabix-indexed file for specified
#' genomic regions. This is the primary method for accessing large GWAS
#' files efficiently.
#'
#' @param sumstats_file Character string. Path to a tabix-indexed (.gz + .tbi)
#'   summary statistics file.
#' @param coloc_regions_PASS A data.frame with columns CHR_var, BP_START_var,
#'   and BP_STOP_var defining genomic regions to extract.
#' @param sumstats_pheno Character string. Either "single" for single phenotype
#'   or "multiple" for multi-phenotype data (default: "single").
#' @param verbose Logical. Whether to print progress messages (default: FALSE).
#' @param test_mode Logical. If TRUE, only process the first region (default: FALSE).
#' @param file_remove Logical. Whether to remove temporary files (default: TRUE).
#' @param return_tabix_cmd Logical. If TRUE, return the tabix command instead of
#'   executing it (default: FALSE).
#'
#' @return A sumstats object (data.table) with attributes:
#'   \itemize{
#'     \item tabix: Status ("tabix_ok", "tabix_ok_no_data", or "tabix_failed")
#'     \item sumstats_file: Path to the source file
#'     \item coloc_regions_PASS: Regions that were queried
#'     \item sumstats_pheno: Phenotype type
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Creates a temporary file with region coordinates
#'   \item Uses tabix to extract data for these regions
#'   \item Parses the output using data.table::fread
#'   \item Sets appropriate status attributes
#' }
#'
#' Tabix status codes:
#' \itemize{
#'   \item "tabix_ok": Data successfully retrieved
#'   \item "tabix_ok_no_data": Query successful but no data in regions
#'   \item "tabix_failed": Query failed (file not found, not indexed, etc.)
#' }
#'
#' @importFrom data.table fread data.table
#' @export
retrieve_sumstats_tabix <- function(sumstats_file, 
                                    coloc_regions_PASS,
                                    sumstats_pheno = "single",
                                    verbose = FALSE, 
                                    test_mode = FALSE,
                                    file_remove = TRUE,
                                    return_tabix_cmd = FALSE) {
  
  # Input validation
  if (!file.exists(sumstats_file)) {
    warning("Summary statistics file not found: ", sumstats_file)
    # Return empty sumstats with failed status
    sumstats <- data.table::data.table()
    attr(sumstats, "tabix") <- "tabix_failed"
    attr(sumstats, "sumstats_file") <- sumstats_file
    attr(sumstats, "coloc_regions_PASS") <- coloc_regions_PASS
    attr(sumstats, "sumstats_pheno") <- sumstats_pheno
    class(sumstats) <- unique(c("sumstats", class(sumstats)))
    return(sumstats)
  }
  
  # Check for tabix index
  if (!file.exists(paste0(sumstats_file, ".tbi"))) {
    warning("Tabix index not found for: ", sumstats_file)
  }
  
  # Validate required columns
  req_cols <- c("CHR_var", "BP_START_var", "BP_STOP_var")
  if (!all(req_cols %in% colnames(coloc_regions_PASS))) {
    stop("coloc_regions_PASS must contain: ", paste(req_cols, collapse = ", "))
  }
  
  # Handle test mode
  if (test_mode) {
    if (verbose) message("Test mode: using only first region")
    coloc_regions_PASS <- coloc_regions_PASS[1, , drop = FALSE]
  }
  
  # Create temporary file for regions
  file_regions <- tempfile(pattern = "regions_", fileext = ".bed")
  
  # Set up cleanup
  if (file_remove && !return_tabix_cmd) {
    on.exit({
      if (file.exists(file_regions)) file.remove(file_regions)
    }, add = TRUE)
  }
  
  # Sort regions by chromosome and start position
  sorted <- coloc_regions_PASS[order(coloc_regions_PASS$CHR_var, 
                                     coloc_regions_PASS$BP_START_var), ]
  
  # Write regions to temporary file
  write.table(sorted, file_regions, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE)
  
  # Construct tabix command
  tabix_cmd <- paste0("tabix -h ", shQuote(sumstats_file), 
                      " -R ", shQuote(file_regions), 
                      " 2>/dev/null")
  
  # Return command if requested
  if (return_tabix_cmd) {
    attr(tabix_cmd, "file_regions") <- file_regions
    return(tabix_cmd)
  }
  
  if (verbose) message("Running tabix: ", tabix_cmd)
  
  # Execute tabix and read results
  sumstats <- tryCatch({
    suppressWarnings(
      data.table::fread(cmd = tabix_cmd, showProgress = FALSE)
    )
  }, error = function(e) {
    if (verbose) message("Error reading tabix output: ", e$message)
    NULL
  })
  
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
  attr(sumstats, "coloc_regions_PASS") <- coloc_regions_PASS
  attr(sumstats, "sumstats_pheno") <- sumstats_pheno
  
  # Set class
  class(sumstats) <- unique(c("sumstats", class(sumstats)))
  
  if (verbose) {
    message("Tabix status: ", tabix_attr)
    if (tabix_attr == "tabix_ok") {
      message("Retrieved ", nrow(sumstats), " variants")
    }
  }
  
  return(sumstats)
}


#' Calculate and set maximum -log10(p-value) attribute
#'
#' @description
#' Calculates the maximum -log10(p-value) in a sumstats object and stores
#' it as an attribute. This is used to quickly assess whether a region
#' contains significant signals.
#'
#' @param sumstats A sumstats object containing a nlog10P column.
#'
#' @return The input sumstats object with an added/updated "max_nlog10P" attribute.
#'
#' @details
#' The maximum value is calculated excluding NA values. If the object is empty
#' or the nlog10P column is missing, the attribute is set to NA.
#'
#' @export
set_max_nlog10P <- function(sumstats) {
  
  # Input validation
  if (!inherits(sumstats, "sumstats")) {
    stop("Expected a sumstats object")
  }
  
  # Calculate max_nlog10P
  if (nrow(sumstats) > 0 && "nlog10P" %in% names(sumstats)) {
    # Remove NA values before calculating max
    nlog10P_values <- sumstats[["nlog10P"]]
    nlog10P_values <- nlog10P_values[!is.na(nlog10P_values)]
    
    if (length(nlog10P_values) > 0) {
      max_nlog10P <- max(nlog10P_values)
    } else {
      max_nlog10P <- NA_real_
    }
  } else {
    max_nlog10P <- NA_real_
  }
  
  # Set attribute
  attr(sumstats, "max_nlog10P") <- max_nlog10P
  
  return(sumstats)
}

# Group 2: QC and formatting functions ----
#' Perform quality control on summary statistics
#'
#' @description
#' Applies quality control filters to remove problematic variants including
#' duplicates, infinite values, NAs, and invalid allele frequencies.
#'
#' @param sumstats A sumstats object to be quality controlled.
#' @param verbose Logical. Whether to print progress messages (default: FALSE).
#'
#' @return A sumstats object with problematic variants removed and a "QC" 
#'   attribute describing what was filtered.
#'
#' @details
#' QC steps include:
#' \enumerate{
#'   \item Skip if tabix failed or no data
#'   \item Remove duplicate variants by Name
#'   \item Remove variants with infinite BETA, SE, or AF
#'   \item Remove variants with NA in BETA, SE, or AF
#'   \item Remove variants with zero BETA, SE, or AF
#'   \item Remove variants with AF = 1 (fixed alleles)
#' }
#'
#' The QC attribute contains codes for what was filtered:
#' \itemize{
#'   \item "ok": No issues found
#'   \item "dupvar": Duplicated variants removed
#'   \item "INFin{column}": Infinite values removed
#'   \item "NAin{column}": NA values removed
#'   \item "0in{column}": Zero values removed
#'   \item "1inAF": Fixed alleles removed
#' }
#'
#' @export
perform_sumstats_qc <- function(sumstats, verbose = FALSE) {
  
  # Input validation
  if (!inherits(sumstats, "sumstats")) {
    stop("Expected a sumstats object")
  }
  
  # Skip QC if tabix failed
  tabix_status <- attr(sumstats, "tabix")
  if (!is.null(tabix_status)) {
    if (tabix_status == "tabix_failed") {
      attr(sumstats, "QC") <- "skipped_tabix_failed"
      return(sumstats)
    }
    if (tabix_status == "tabix_ok_no_data") {
      attr(sumstats, "QC") <- "no_data"
      return(sumstats)
    }
  }
  
  # Track QC steps
  QC_codes <- character()
  initial_rows <- nrow(sumstats)
  
  # Remove duplicated variants
  if ("Name" %in% colnames(sumstats)) {
    dups <- duplicated(sumstats$Name)
    if (any(dups)) {
      if (verbose) message("Removing ", sum(dups), " duplicated variants")
      QC_codes <- c(QC_codes, "dupvar")
      sumstats <- sumstats[!dups, ]
    }
  }
  
  # QC numeric columns
  numeric_cols <- c("BETA", "SE", "AF")
  
  for (col in numeric_cols) {
    if (!col %in% colnames(sumstats)) {
      if (verbose) message("Column ", col, " not found, skipping QC")
      next
    }
    
    # Check for infinite values
    inf_vals <- is.infinite(sumstats[[col]])
    if (any(inf_vals)) {
      if (verbose) message("Removing ", sum(inf_vals), " infinite values in ", col)
      QC_codes <- c(QC_codes, paste0("INFin", col))
      sumstats <- sumstats[!inf_vals, ]
    }
    
    # Check for NA values
    na_vals <- is.na(sumstats[[col]])
    if (any(na_vals)) {
      if (verbose) message("Removing ", sum(na_vals), " NA values in ", col)
      QC_codes <- c(QC_codes, paste0("NAin", col))
      sumstats <- sumstats[!na_vals, ]
    }
    
    # Check for zero values (problematic for SE and AF)
    zero_vals <- sumstats[[col]] == 0
    if (any(zero_vals)) {
      if (verbose) message("Removing ", sum(zero_vals), " zero values in ", col)
      QC_codes <- c(QC_codes, paste0("0in", col))
      sumstats <- sumstats[!zero_vals, ]
    }
  }
  
  # Special check for AF = 1 (fixed alleles)
  if ("AF" %in% colnames(sumstats)) {
    one_vals <- sumstats[["AF"]] == 1
    if (any(one_vals)) {
      if (verbose) message("Removing ", sum(one_vals), " fixed alleles (AF=1)")
      QC_codes <- c(QC_codes, "1inAF")
      sumstats <- sumstats[!one_vals, ]
    }
  }
  
  # Set QC attribute
  if (length(QC_codes) == 0) {
    attr(sumstats, "QC") <- "ok"
  } else {
    attr(sumstats, "QC") <- paste(QC_codes, collapse = "_")
  }
  
  # Report total filtered
  if (verbose) {
    filtered_rows <- initial_rows - nrow(sumstats)
    if (filtered_rows > 0) {
      message("QC completed: removed ", filtered_rows, " of ", initial_rows, 
              " variants (", round(filtered_rows/initial_rows * 100, 1), "%)")
    } else {
      message("QC completed: no variants removed")
    }
  }
  
  return(sumstats)
}


#' Validate and set attributes for sumstats object
#'
#' @description
#' Ensures a sumstats object has all required attributes and validates
#' that column names match the expected standardized format.
#'
#' @param sumstats A sumstats object to validate.
#' @param sumstats_function Character string. Name of the function that created the data.
#' @param sumstats_type Character string. Type of summary statistics ('quant' or 'cc').
#' @param sumstats_sdY Numeric. Standard deviation for quantitative traits.
#'
#' @return The validated sumstats object with all required attributes set.
#'
#' @details
#' Required attributes:
#' \itemize{
#'   \item sumstats_function: Function used to retrieve data
#'   \item sumstats_file: Path to source file
#'   \item tabix: Tabix status
#'   \item coloc_regions_PASS: Regions queried
#'   \item sumstats_type: Trait type (quant/cc)
#'   \item sumstats_sdY: Standard deviation (for quant)
#'   \item sumstats_pheno: Phenotype type (single/multiple)
#' }
#'
#' The function also verifies that column names match the standard format
#' returned by \code{get_cols_to()}.
#'
#' @seealso 
#' \code{\link{get_cols_to}} for standard column names
#'
#' @export
check_sumstats_attributes <- function(sumstats, 
                                      sumstats_function,
                                      sumstats_type, 
                                      sumstats_sdY) {
  
  # Check object type
  if (!inherits(sumstats, "sumstats")) {
    stop("Expected a sumstats object")
  }
  
  # Validate sumstats_type
  if (!missing(sumstats_type)) {
    if (!sumstats_type %in% c("quant", "cc")) {
      stop("sumstats_type must be 'quant' or 'cc'")
    }
    attr(sumstats, "sumstats_type") <- sumstats_type
  }
  
  # Set sumstats_sdY if provided
  if (!missing(sumstats_sdY)) {
    if (!is.na(sumstats_sdY) && (!is.numeric(sumstats_sdY) || sumstats_sdY <= 0)) {
      warning("sumstats_sdY should be a positive number or NA")
    }
    attr(sumstats, "sumstats_sdY") <- sumstats_sdY
  }
  
  # List of required attributes
  req_attrs <- c("sumstats_function", "sumstats_file",
                 "tabix", "coloc_regions_PASS",
                 "sumstats_type", "sumstats_sdY",
                 "sumstats_pheno")
  
  # Check for missing attributes
  missing_attrs <- req_attrs[sapply(req_attrs, function(a) is.null(attr(sumstats, a)))]
  
  if (length(missing_attrs) > 0) {
    stop("Missing required attributes: ", paste(missing_attrs, collapse = ", "))
  }
  
  # Validate column names if data is not empty
  if (nrow(sumstats) > 0) {
    expected_cols <- get_cols_to()
    actual_cols <- colnames(sumstats)
    
    # Check if all expected columns are present
    missing_cols <- setdiff(expected_cols, actual_cols)
    if (length(missing_cols) > 0) {
      stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # Check column order (warning only)
    if (!identical(actual_cols[1:length(expected_cols)], expected_cols)) {
      warning("Column order doesn't match expected format. ",
              "Expected: ", paste(expected_cols, collapse = ", "))
    }
  }
  
  return(sumstats)
}


#' Format summary statistics with QC and validation
#'
#' @description
#' A wrapper function that performs quality control and attribute validation
#' on summary statistics. This is typically called after retrieving raw data.
#'
#' @param sumstats A sumstats object to format.
#' @param sumstats_type Character string. Type of trait ('quant' or 'cc').
#' @param sumstats_sdY Numeric. Standard deviation for quantitative traits.
#' @param verbose Logical. Whether to print QC messages (default: FALSE).
#'
#' @return A formatted sumstats object that has passed QC and validation.
#'
#' @details
#' This function combines two steps:
#' \enumerate{
#'   \item Quality control via \code{perform_sumstats_qc}
#'   \item Attribute validation via \code{check_sumstats_attributes}
#' }
#'
#' Note: set_max_nlog10P is not called here as it's used separately
#' downstream in the pipeline.
#'
#' @seealso 
#' \code{\link{perform_sumstats_qc}} for QC details
#' \code{\link{check_sumstats_attributes}} for validation details
#'
#' @export
format_sumstats <- function(sumstats,
                            sumstats_type,
                            sumstats_sdY,
                            verbose = FALSE) {
  
  # Step 1: Perform quality control
  sumstats <- perform_sumstats_qc(
    sumstats = sumstats, 
    verbose = verbose
  )
  
  # Step 2: Check and set attributes
  sumstats <- check_sumstats_attributes(
    sumstats = sumstats,
    sumstats_function = attr(sumstats, "sumstats_function"),
    sumstats_type = sumstats_type,
    sumstats_sdY = sumstats_sdY
  )
  
  return(sumstats)
}


# Group 3: column standardization functions ----
#' Get standard column names for GWAS summary statistics
#'
#' @description
#' Returns the standardized column names required for the colocalization
#' pipeline. All summary statistics must be converted to this format.
#'
#' @return Character vector of standard column names in the required order:
#'   \itemize{
#'     \item Name: Variant identifier (e.g., "chr1:12345:A:G")
#'     \item rsID: RS identifier (e.g., "rs12345")
#'     \item CHR: Chromosome (numeric or "X", "Y")
#'     \item POS: Base pair position
#'     \item A1: Effect allele
#'     \item A2: Other allele
#'     \item BETA: Effect size
#'     \item SE: Standard error
#'     \item nlog10P: -log10(p-value)
#'     \item AF: Allele frequency of A1
#'     \item N: Sample size
#'   }
#'
#' @details
#' This standardized format ensures compatibility across different
#' data sources and analysis functions. The order is important for
#' some legacy functions.
#'
#' @export
get_cols_to <- function() {
  c("Name", "rsID", "CHR", "POS", "A1", "A2", 
    "BETA", "SE", "nlog10P", "AF", "N") 
}


#' Map non-standard column names to standard format
#'
#' @description
#' Converts a GWAS summary statistics data.table with non-standard column
#' names to the standardized format required for colocalization analysis.
#'
#' @param sumstats A data.table containing GWAS summary statistics.
#' @param Name Character string. Column name for variant identifier.
#' @param rsID Character string. Column name for RS ID.
#' @param CHR Character string. Column name for chromosome.
#' @param POS Character string. Column name for position.
#' @param A1 Character string. Column name for effect allele.
#' @param A2 Character string. Column name for other allele.
#' @param BETA Character string. Column name for effect size.
#' @param SE Character string. Column name for standard error.
#' @param nlog10P Character string. Column name for -log10(p-value).
#' @param AF Character string. Column name for allele frequency.
#' @param N Character string. Column name for sample size.
#' @param Phenotype Character string. Optional column name for phenotype
#'   (used for multi-phenotype datasets).
#'
#' @return A data.table with standardized column names in the correct order.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Validates that all specified columns exist
#'   \item Renames columns to standard names
#'   \item Reorders columns to match the standard format
#'   \item Optionally includes a Phenotype column for multi-trait data
#' }
#'
#' @importFrom data.table is.data.table
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert UKBB format to standard
#' standardized <- match_cols(
#'   ukbb_data,
#'   Name = "variant",
#'   rsID = "rsid",
#'   CHR = "chromosome",
#'   POS = "position",
#'   A1 = "alt",
#'   A2 = "ref",
#'   BETA = "beta",
#'   SE = "se",
#'   nlog10P = "neglog10_pval",
#'   AF = "minor_AF",
#'   N = "n_complete_samples"
#' )
#' }
match_cols <- function(sumstats, Name, rsID, CHR, POS, A1, A2, 
                       BETA, SE, nlog10P, AF, N, Phenotype = NULL) {
  
  # Input validation
  if (!data.table::is.data.table(sumstats)) {
    stop("Input must be a data.table object")
  }
  
  # Create list of input columns
  input_cols <- c(Name, rsID, CHR, POS, A1, A2, BETA, SE, nlog10P, AF, N)
  
  # Check for missing columns
  missing_cols <- setdiff(input_cols, colnames(sumstats))
  if (length(missing_cols) > 0) {
    stop("Missing columns in input data: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  # Handle optional Phenotype column
  if (!is.null(Phenotype)) {
    if (!Phenotype %in% colnames(sumstats)) {
      stop("Specified Phenotype column '", Phenotype, "' not found")
    }
    input_cols <- c(input_cols, Phenotype)
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
  
  if (!is.null(Phenotype)) {
    col_mapping <- c(col_mapping, Phenotype = Phenotype)
  }
  
  # Make a copy to avoid modifying the original
  sumstats_copy <- copy(sumstats)
  
  # Rename columns using setnames for efficiency
  for (std_name in names(col_mapping)) {
    old_name <- col_mapping[std_name]
    if (old_name != std_name) {  # Only rename if different
      data.table::setnames(sumstats_copy, old = old_name, new = std_name)
    }
  }
  
  # Get standard column order
  std_cols <- get_cols_to()
  if (!is.null(Phenotype)) {
    std_cols <- c(std_cols, "Phenotype")
  }
  
  # Select and reorder columns
  sumstats_copy <- sumstats_copy[, ..std_cols]
  
  return(sumstats_copy)
}
