# Group 1: Ontology Lookup Service (OLS) ----
#' Search for EFO terms using the Ontology Lookup Service
#' 
#' @description
#' Queries the EBI Ontology Lookup Service (OLS) API to search for terms
#' in the Experimental Factor Ontology (EFO). This is useful for finding
#' standardized trait identifiers for GWAS catalog queries.
#' 
#' @param query Character string. Search term(s) to query in the EFO.
#' @param rows Integer. Maximum number of results to return (default: 50).
#' 
#' @return A data.frame with columns:
#'   \itemize{
#'     \item label: Human-readable term label
#'     \item iri: Full IRI (Internationalized Resource Identifier)
#'     \item short_form: EFO identifier (e.g., "EFO_0000270")
#'   }
#'   Returns NULL if no results found or if the query fails.
#' 
#' @details
#' The function searches across all EFO fields including labels, synonyms,
#' and definitions. Only terms with EFO identifiers (starting with "EFO_")
#' are returned.
#' 
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' 
#' @examples
#' \dontrun{
#' # Search for kidney-related terms
#' kidney_terms <- search_efo_terms("kidney disease")
#' 
#' # Search for diabetes
#' diabetes_terms <- search_efo_terms("type 2 diabetes")
#' }
#' 
#' @export
search_efo_terms <- function(query, rows = 50) {
  
  # Input validation
  if (!is.character(query) || length(query) != 1) {
    stop("query must be a single character string")
  }
  
  if (!is.numeric(rows) || rows < 1 || rows > 1000) {
    stop("rows must be a number between 1 and 1000")
  }
  
  # API endpoint
  base_url <- "https://www.ebi.ac.uk/ols/api/search"
  
  # Make API request
  response <- tryCatch({
    httr::GET(
      base_url, 
      query = list(
        q = query,
        ontology = "efo",
        rows = rows
      ),
      httr::timeout(30)
    )
  }, error = function(e) {
    warning("Failed to query OLS API: ", e$message)
    return(NULL)
  })
  
  # Check response
  if (is.null(response)) return(NULL)
  
  if (httr::status_code(response) == 200) {
    # Parse JSON response
    content <- jsonlite::fromJSON(
      rawToChar(response$content),
      simplifyVector = TRUE
    )
    
    # Extract results if available
    if (!is.null(content$response$docs) && length(content$response$docs) > 0) {
      results <- data.frame(
        label = content$response$docs$label,
        iri = content$response$docs$iri,
        short_form = content$response$docs$short_form,
        stringsAsFactors = FALSE
      )
      
      # Filter for EFO terms only
      efo_results <- results[grepl("^EFO_", results$short_form), ]
      
      if (nrow(efo_results) == 0) {
        message("No EFO terms found for query: ", query)
        return(NULL)
      }
      
      return(efo_results)
    }
  } else {
    warning("OLS API returned status code: ", httr::status_code(response))
  }
  
  return(NULL)
}
# Group 2: GWAS Catalog ----
#' Get GWAS associations for an EFO term
#' 
#' @description
#' Retrieves all GWAS catalog associations for a given EFO term. This includes
#' variant information, effect sizes, p-values, and mapped genes.
#' 
#' @param efo_id Character string. EFO identifier (e.g., "EFO_0000270").
#' @param efo_id_label Character string. Human-readable label for the EFO term.
#' 
#' @return A data.frame with columns:
#'   \itemize{
#'     \item variant_risk_allele: Risk allele identifier (rsID or chr:pos)
#'     \item location: Genomic location (chr:position format)
#'     \item pvalue: Association p-value
#'     \item raf: Risk allele frequency (or "NR" if not reported)
#'     \item beta: Effect size with direction
#'     \item ci: Confidence interval (or "-" if not available)
#'     \item mapped_gene: Genes mapped to the variant
#'     \item efo_id: EFO identifier
#'     \item efo_name: EFO term label
#'   }
#'   Returns a minimal data.frame with just efo_id and efo_name if no associations found.
#' 
#' @details
#' The function processes the GWAS catalog API response to extract key
#' association information. P-values are formatted in scientific notation,
#' and beta values include direction signs.
#' 
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' 
#' @seealso 
#' \code{\link{search_efo_terms}} to find EFO identifiers
#' \code{\link{get_snp_details}} for variant location lookup
#' 
#' @export
get_gwas_associations <- function(efo_id, efo_id_label) {
  
  # Input validation
  if (!grepl("^EFO_[0-9]+$", efo_id)) {
    warning("Invalid EFO ID format: ", efo_id)
  }
  
  message("Querying GWAS catalog for ", efo_id, " (", efo_id_label, ")")
  
  # API endpoint
  base_url <- paste0("https://www.ebi.ac.uk/gwas/rest/api/efoTraits/", 
                     efo_id, "/associations")
  
  # Make API request
  response <- tryCatch({
    httr::GET(base_url, httr::timeout(60))
  }, error = function(e) {
    warning("Failed to query GWAS catalog: ", e$message)
    return(NULL)
  })
  
  # Return minimal data.frame if request failed
  if (is.null(response)) {
    return(data.frame(
      efo_id = efo_id,
      efo_name = efo_id_label,
      stringsAsFactors = FALSE
    ))
  }
  
  if (httr::status_code(response) == 200) {
    # Parse response
    content <- jsonlite::fromJSON(
      rawToChar(response$content),
      simplifyVector = FALSE
    )
    
    # Check if associations exist
    if (!is.null(content$`_embedded`$associations) && 
        length(content$`_embedded`$associations) > 0) {
      
      associations <- content$`_embedded`$associations
      
      # Process each association
      results <- lapply(seq_along(associations$loci), function(i) {
        tryCatch({
          # Extract variant information
          locus <- associations$loci[[i]]
          risk_alleles <- locus$strongestRiskAlleles[[1]]
          risk_allele <- risk_alleles$riskAlleleName[1]
          
          # Get genomic location
          if (grepl("^rs", risk_allele)) {
            # For rsIDs, query for position
            location <- get_snp_details(risk_allele)
          } else {
            # For chr:pos format, clean up
            location <- gsub("-\\?", "", risk_allele)
            location <- gsub("chr", "", location)
          }
          
          # Format p-value
          pvalue <- tryCatch({
            as.numeric(sprintf("%.0e", 
                               associations$pvalueMantissa[i] * 
                                 10^associations$pvalueExponent[i]))
          }, error = function(e) NA)
          
          # Format beta with direction
          beta_value <- if (!is.na(associations$betaNum[i])) {
            sign_char <- ifelse(associations$betaDirection[i] == "decrease", "-", "+")
            as.numeric(paste0(sign_char, sprintf("%.3f", associations$betaNum[i])))
          } else {
            NA
          }
          
          # Extract mapped genes
          mapped_genes <- if (length(locus$authorReportedGenes) > 0 && 
                              length(locus$authorReportedGenes[[1]]) > 0) {
            paste(unique(unlist(locus$authorReportedGenes)), collapse = ", ")
          } else {
            "-"
          }
          
          # Create result row
          data.frame(
            variant_risk_allele = risk_allele,
            location = ifelse(is.na(location), "-", location),
            pvalue = pvalue,
            raf = ifelse(is.na(associations$riskFrequency[i]), "NR", 
                         sprintf("%.3f", associations$riskFrequency[i])),
            beta = beta_value,
            ci = ifelse(is.na(associations$range[i]), "-", associations$range[i]),
            mapped_gene = mapped_genes,
            efo_id = efo_id,
            efo_name = efo_id_label,
            stringsAsFactors = FALSE
          )
        }, error = function(e) {
          warning("Error processing association ", i, ": ", e$message)
          NULL
        })
      })
      
      # Remove NULL results and combine
      results <- results[!sapply(results, is.null)]
      if (length(results) > 0) {
        return(do.call(rbind, results))
      }
    }
  } else {
    warning("GWAS catalog API returned status: ", httr::status_code(response))
  }
  
  # Return minimal data.frame if no associations found
  return(data.frame(
    efo_id = efo_id,
    efo_name = efo_id_label,
    stringsAsFactors = FALSE
  ))
}


#' Get genomic coordinates for an rsID
#' 
#' @description
#' Queries the GWAS catalog API to retrieve chromosomal coordinates
#' for a given rsID.
#' 
#' @param rsid Character string. The rsID to look up (e.g., "rs12345").
#' 
#' @return Character string in "chr:position" format, or NA if lookup fails.
#' 
#' @details
#' This helper function is used internally by get_gwas_associations to
#' convert rsIDs to genomic coordinates. It handles compound rsIDs
#' (e.g., "rs12345-A") by extracting just the rsID portion.
#' 
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' 
#' @export
get_snp_details <- function(rsid) {
  
  # Extract just the rsID part (handle compound IDs like "rs12345-A")
  rsid_clean <- strsplit(rsid, "-")[[1]][1]
  
  # Validate rsID format
  if (!grepl("^rs[0-9]+$", rsid_clean)) {
    warning("Invalid rsID format: ", rsid)
    return(NA)
  }
  
  # API endpoint
  snp_url <- paste0("https://www.ebi.ac.uk/gwas/rest/api/",
                    "singleNucleotidePolymorphisms/", rsid_clean)
  
  # Query API
  response <- tryCatch({
    httr::GET(snp_url, httr::timeout(10))
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(response) && httr::status_code(response) == 200) {
    content <- tryCatch({
      jsonlite::fromJSON(rawToChar(response$content))
    }, error = function(e) NULL)
    
    if (!is.null(content) && !is.null(content$locations)) {
      # Extract chromosome and position
      chr <- content$locations$chromosomeName[1]
      pos <- content$locations$chromosomePosition[1]
      
      if (!is.na(chr) && !is.na(pos)) {
        return(paste(chr, pos, sep = ":"))
      }
    }
  }
  
  return(NA)
}

# Group 3: GWAS Catalog Download Functions ----
#' Download GWAS Catalog study files
#' 
#' @description
#' Downloads all files for a GWAS Catalog study (GCST identifier) from the
#' EBI FTP server. Automatically detects and downloads harmonized files if available.
#' 
#' @param url_in Character string. URL to the study on GWAS catalog
#'   (e.g., "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCSTXXXXXX").
#' @param dest_folder Character string. Local directory to save files.
#' 
#' @details
#' The function:
#' \enumerate{
#'   \item Converts HTTP URLs to FTP protocol
#'   \item Lists all files in the study directory
#'   \item Checks for harmonized subdirectory and uses it if available
#'   \item Downloads files that don't already exist locally
#' }
#' 
#' Files are downloaded using wget for reliability. Existing files are
#' not re-downloaded.
#' 
#' @importFrom RCurl getURL
#' 
#' @examples
#' \dontrun{
#' # Download a GWAS study
#' download_GCST(
#'   url_in = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST123456",
#'   dest_folder = "data/gwas/GCST123456"
#' )
#' }
#' 
#' @export
download_GCST <- function(url_in, dest_folder) {
  
  # Input validation
  if (missing(dest_folder) || !is.character(dest_folder)) {
    stop("dest_folder must be provided as a character string")
  }
  
  if (!is.character(url_in) || length(url_in) != 1) {
    stop("url_in must be a single character string")
  }
  
  # Create destination directory
  if (!dir.exists(dest_folder)) {
    dir.create(dest_folder, recursive = TRUE)
  }
  
  # Convert HTTP to FTP and ensure trailing slash
  url_in <- gsub("^http://", "ftp://", url_in)
  url_in <- gsub("/*$", "/", url_in)
  
  tryCatch({
    # List files in directory
    message("Listing files from: ", url_in)
    files <- RCurl::getURL(url_in, ftp.use.epsv = FALSE, dirlistonly = TRUE)
    file_list <- strsplit(files, "\n|\r")[[1]]
    file_list <- file_list[file_list != ""]  # Remove empty entries
    
    # Check for harmonised subdirectory
    if ("harmonised" %in% file_list) {
      message("Found harmonised directory, using harmonised files")
      url_in <- paste0(url_in, "harmonised/")
      files <- RCurl::getURL(url_in, ftp.use.epsv = FALSE, dirlistonly = TRUE)
      file_list <- strsplit(files, "\n|\r")[[1]]
      file_list <- file_list[file_list != ""]
    }
    
    # Download each file
    message("Found ", length(file_list), " files to process")
    
    for (dest_file in file_list) {
      local_path <- file.path(dest_folder, dest_file)
      
      if (!file.exists(local_path)) {
        message("Downloading: ", dest_file)
        download_url <- paste0(url_in, dest_file)
        
        # Use wget for reliable download
        cmd <- sprintf("wget -q '%s' -O '%s'", download_url, local_path)
        result <- system(cmd)
        
        if (result != 0) {
          warning("Failed to download: ", dest_file)
          if (file.exists(local_path)) file.remove(local_path)
        }
      } else {
        message("File already exists, skipping: ", dest_file)
      }
    }
    
    message("Download complete")
    
  }, error = function(e) {
    stop("Error accessing ", url_in, ": ", e$message)
  })
}


#' Process GWAS Catalog summary statistics file
#' 
#' @description
#' Reads and standardizes GWAS Catalog summary statistics files, handling
#' various column naming conventions and calculating missing values where possible.
#' 
#' @param file_in Character string. Path to the GWAS catalog summary statistics file.
#' @param test_mode Logical. If TRUE, only read first 2 lines for testing (default: FALSE).
#' 
#' @return A data.table with standardized columns:
#'   chromosome, base_pair_location, effect_allele, other_allele,
#'   beta/effect_size, standard_error, p_value, effect_allele_frequency, n
#' 
#' @details
#' The function handles:
#' \itemize{
#'   \item Multiple p-value column names (pval, p_value)
#'   \item Odds ratios: converts to log scale and calculates SE from CI
#'   \item Missing sample sizes: attempts to read from metadata YAML
#'   \item Special cases like GCST files with boolean CI columns
#' }
#' 
#' @importFrom data.table fread
#' 
#' @export
process_GCST <- function(file_in, test_mode = FALSE) {
  
  # Input validation
  if (!file.exists(file_in)) {
    stop("File not found: ", file_in)
  }
  
  # Read data
  if (test_mode) {
    message("Test mode: reading first 2 lines only")
    sumstats <- data.table::fread(cmd = paste0("zcat ", shQuote(file_in), " | head -2"))
  } else {
    sumstats <- data.table::fread(file_in, showProgress = FALSE)
  }
  
  if (nrow(sumstats) == 0) {
    stop("No data found in file: ", file_in)
  }
  
  # Find p-value column
  p_value_col <- NULL
  for (col in c("pval", "p_value", "p-value", "P-value")) {
    if (col %in% colnames(sumstats)) {
      p_value_col <- col
      break
    }
  }
  if (is.null(p_value_col)) {
    stop("P-value column not found. Available columns: ", 
         paste(colnames(sumstats), collapse = ", "))
  }
  
  # Find effect size column
  es_col <- "beta"
  if (!es_col %in% colnames(sumstats)) {
    if (all(c("odds_ratio", "ci_lower", "ci_upper") %in% colnames(sumstats))) {
      message("Converting odds ratios to log scale")
      
      # Calculate beta as log(OR)
      sumstats[[es_col]] <- log(sumstats$odds_ratio)
      
      # Calculate standard error from confidence intervals
      # Special case: when CI columns are boolean (TRUE)
      if (is.logical(sumstats$ci_lower) && all(sumstats$ci_lower == TRUE) && 
          all(sumstats$ci_upper == TRUE) && p_value_col %in% colnames(sumstats)) {
        message("CI columns are boolean, calculating SE from p-values")
        z_score <- qnorm(sumstats[[p_value_col]]/2, lower.tail = FALSE)
        sumstats$standard_error <- abs(sumstats[[es_col]]) / z_score
      } else {
        # Standard calculation from CI
        sumstats$standard_error <- (log(sumstats$ci_upper) - log(sumstats$ci_lower)) / (2 * 1.96)
      }
    } else {
      stop("Effect size column not found. Available columns: ", 
           paste(colnames(sumstats), collapse = ", "))
    }
  }
  
  # Find or calculate sample size
  N_col <- "n"
  if (!N_col %in% colnames(sumstats)) {
    if (all(c("num_cases", "num_controls") %in% colnames(sumstats))) {
      message("Calculating sample size from cases and controls")
      sumstats[[N_col]] <- sumstats$num_cases + sumstats$num_controls
      
    } else if (file.exists(paste0(file_in, "-meta.yaml"))) {
      # Try to extract from metadata
      message("Attempting to extract sample size from metadata YAML")
      yaml_file <- paste0(file_in, "-meta.yaml")
      
      tryCatch({
        ss_txt <- system(paste0("grep sample_size ", shQuote(yaml_file)), 
                         intern = TRUE, ignore.stderr = TRUE)
        if (length(ss_txt) > 0) {
          # Extract numeric values
          n_values <- as.numeric(gsub("[^0-9]", "", ss_txt))
          n_values <- n_values[!is.na(n_values)]
          
          if (length(n_values) > 0) {
            n_tmp <- sum(unique(n_values))
            sumstats[[N_col]] <- n_tmp
            message("Sample size from YAML: ", n_tmp, " (please verify)")
          } else {
            stop("Could not parse sample size from YAML")
          }
        } else {
          stop("No sample_size found in YAML")
        }
      }, error = function(e) {
        stop("Sample size column not found and could not read from YAML: ", e$message)
      })
    } else {
      stop("Sample size information not found")
    }
  }
  
  # Standardize and validate output columns
  required_cols <- c("chromosome", "base_pair_location", "effect_allele", 
                     "other_allele", es_col, "standard_error", 
                     p_value_col, "effect_allele_frequency", N_col)
  
  missing_cols <- setdiff(required_cols, colnames(sumstats))
  if (length(missing_cols) > 0) {
    stop("Missing required columns after processing: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  # Select and reorder columns
  sumstats <- sumstats[, ..required_cols]
  
  # Standardize column names
  data.table::setnames(sumstats, 
                       old = required_cols,
                       new = c("CHR", "POS", "A1", "A2", "BETA", "SE", 
                               "P", "AF", "N"))
  
  return(sumstats)
}

# eQTL Catalogue ----
#' Query eQTL Catalogue for associations in a genomic region
#' 
#' @description
#' Retrieves eQTL associations from the eQTL Catalogue API for a specific
#' dataset and genomic region. Handles pagination automatically to retrieve
#' all results.
#' 
#' @param sumstats_file Character string. Dataset identifier from eQTL Catalogue.
#' @param CHR_var Character or numeric. Chromosome.
#' @param BP_START_var Numeric. Start position of the region.
#' @param BP_STOP_var Numeric. End position of the region.
#' @param start Integer. Starting index for pagination (default: 0).
#' @param size Integer. Number of results per page (default: 1000).
#' 
#' @return A list of data.frames, one per molecular trait (gene), with
#'   standardized column names suitable for colocalization analysis.
#' 
#' @details
#' The function:
#' \enumerate{
#'   \item Queries the eQTL Catalogue API with pagination
#'   \item Handles rate limiting with built-in delays
#'   \item Formats results to standard column names
#'   \item Splits results by molecular trait (gene)
#' }
#' 
#' @importFrom httr GET accept_json content status_code
#' @importFrom jsonlite fromJSON
#' 
#' @export
query_eQTL_Catalogue <- function(sumstats_file,
                                 CHR_var, BP_START_var, BP_STOP_var,
                                 start = 0, size = 1000) {
  
  # Input validation
  if (!is.numeric(size) || size < 1 || size > 1000) {
    stop("size must be between 1 and 1000")
  }
  
  # Format position string
  pos <- paste0(CHR_var, ":", BP_START_var, "-", BP_STOP_var)
  message("Querying eQTL Catalog for dataset: ", sumstats_file, ", region: ", pos)
  
  # Initialize results
  all_results <- list()
  
  # Pagination loop
  repeat {
    # Construct API URL
    sumstats_URL <- paste0(
      "https://www.ebi.ac.uk/eqtl/api/v2/datasets/",
      sumstats_file, "/associations?",
      "size=", size, 
      "&start=", start,
      "&pos=", pos
    )
    
    message("Fetching records: ", start, "-", start + size, "... ", appendLF = FALSE)
    
    # Make API request
    response <- tryCatch({
      httr::GET(sumstats_URL, httr::accept_json(), httr::timeout(60))
    }, error = function(e) {
      warning("API request failed: ", e$message)
      return(NULL)
    })
    
    # Rate limiting - wait before next request
    Sys.sleep(2)
    
    # Check response
    if (is.null(response)) {
      message("Failed")
      break
    }
    
    if (httr::status_code(response) == 200) {
      message("Success")
      
      # Parse JSON response
      content_text <- httr::content(response, "text", encoding = "UTF-8")
      sumstats_chunk <- jsonlite::fromJSON(content_text, simplifyVector = TRUE)
      
      # Check if we got data
      if (is.null(sumstats_chunk) || length(sumstats_chunk) == 0) {
        message("No more data available")
        break
      }
      
      # Accumulate results
      if (start == 0) {
        sumstats <- sumstats_chunk
      } else {
        sumstats <- rbind(sumstats, sumstats_chunk)
      }
      
      # Check if we got fewer results than requested (last page)
      if (nrow(sumstats_chunk) < size) {
        message("Reached end of results")
        break
      }
      
      # Move to next page
      start <- start + size
      
    } else {
      # Handle different status codes
      if (start == 0) {
        message("No data available for this region")
        return(data.frame())
      } else {
        message("Reached end of available data")
        break
      }
    }
  }
  
  # Format results if we have any
  if (exists("sumstats") && nrow(sumstats) > 0) {
    message("Retrieved ", nrow(sumstats), " associations")
    sumstats <- format_eQTL_Catalogue(sumstats)
  } else {
    message("No associations found")
    return(list())
  }
  
  return(sumstats)
}


#' Format eQTL Catalogue results to standard format
#' 
#' @description
#' Converts eQTL Catalogue API results to the standardized format required
#' for colocalization analysis, splitting by molecular trait.
#' 
#' @param sumstats Data.frame from eQTL Catalogue API query.
#' 
#' @return A named list of data.frames, one per molecular trait (gene),
#'   with standardized column names.
#' 
#' @details
#' The function:
#' \itemize{
#'   \item Creates variant names in chr:pos:ref:alt format
#'   \item Renames columns to standard names
#'   \item Splits data by molecular trait (gene)
#'   \item Returns a list suitable for multi-phenotype analysis
#' }
#' 
#' @export
format_eQTL_Catalogue <- function(sumstats) {
  
  # Check for required columns
  required_cols <- c("chromosome", "position", "ref", "alt", "rsid", 
                     "beta", "se", "nlog10p", "maf", "an", 
                     "molecular_trait_id")
  
  missing_cols <- setdiff(required_cols, colnames(sumstats))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create standardized variant name
  sumstats$Name <- paste0("chr", sumstats$chromosome, ":", 
                          sumstats$position, ":",
                          sumstats$ref, ":", sumstats$alt)
  
  # Select and rename columns
  cols_to_keep <- c("Name", "rsid", "chromosome", "position", "alt", "ref", 
                    "beta", "se", "nlog10p", "maf", "an", "molecular_trait_id")
  sumstats <- sumstats[, cols_to_keep]
  
  # Rename to standard column names
  colnames(sumstats) <- c("Name", "rsID", "CHR", "POS", "A1", "A2", 
                          "BETA", "SE", "nlog10P", "AF", "N", "Phenotype")
  
  # Split by molecular trait (gene)
  phenotypes <- unique(sumstats$Phenotype)
  
  if (length(phenotypes) == 0) {
    warning("No phenotypes found in data")
    return(list())
  }
  
  # Create named list of data.frames
  sumstats_list <- lapply(phenotypes, function(pheno) {
    subset(sumstats, Phenotype == pheno)
  })
  names(sumstats_list) <- phenotypes
  
  message("Formatted data for ", length(phenotypes), " molecular traits")
  
  return(sumstats_list)
}


#' Get eQTL Catalogue dataset information
#' 
#' @description
#' Retrieves metadata for all datasets available in the eQTL Catalogue,
#' including study information, tissue types, and analysis methods.
#' 
#' @return A data.frame containing dataset metadata with columns including:
#'   dataset_id, study_id, tissue_id, condition, quant_method, etc.
#' 
#' @details
#' This function is useful for:
#' \itemize{
#'   \item Finding available datasets for a tissue of interest
#'   \item Filtering datasets by quantification method (ge, tx, etc.)
#'   \item Getting study metadata for result annotation
#' }
#' 
#' @importFrom httr GET accept_json content status_code
#' @importFrom jsonlite fromJSON
#' 
#' @examples
#' \dontrun{
#' # Get all available datasets
#' datasets <- get_datasets_eQTL_Catalogue()
#' 
#' # Filter for gene expression in specific tissue
#' blood_ge <- subset(datasets, quant_method == "ge" & 
#'                    grepl("blood", tissue_label, ignore.case = TRUE))
#' }
#' 
#' @export
get_datasets_eQTL_Catalogue <- function() {
  
  message("Retrieving eQTL Catalogue dataset information...")
  
  # API endpoint
  datasets_url <- "https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size=1000&start=0"
  
  # Make API request
  response <- tryCatch({
    httr::GET(datasets_url, httr::accept_json(), httr::timeout(30))
  }, error = function(e) {
    stop("Failed to query eQTL Catalogue: ", e$message)
  })
  
  # Check response
  if (httr::status_code(response) != 200) {
    stop("eQTL Catalogue API returned status: ", httr::status_code(response))
  }
  
  # Parse response
  content_text <- httr::content(response, "text", encoding = "UTF-8")
  datasets <- jsonlite::fromJSON(content_text, simplifyVector = TRUE)
  
  if (is.null(datasets) || nrow(datasets) == 0) {
    warning("No datasets retrieved from eQTL Catalogue")
    return(data.frame())
  }
  
  message("Retrieved information for ", nrow(datasets), " datasets")
  
  return(datasets)
}


#' Annotate colocalization results with eQTL Catalogue metadata
#' 
#' @description
#' Adds gene annotations and dataset metadata to colocalization results
#' from eQTL Catalogue analyses.
#' 
#' @param coloc_out Data.frame. Colocalization results with sumstats_2_file column.
#' @param annotation_df Data.frame. Gene annotations (e.g., ENSG to gene symbol mapping).
#' @param datasets_eQTL_Catalogue Data.frame. Dataset metadata from get_datasets_eQTL_Catalogue().
#' 
#' @return The input coloc_out data.frame with additional annotation columns:
#'   \itemize{
#'     \item dataset_id: Extracted dataset identifier
#'     \item gene_id: Extracted gene identifier (ENSG)
#'     \item Columns from annotation_df (e.g., gene_symbol)
#'     \item Columns from datasets_eQTL_Catalogue (e.g., tissue_label, study_id)
#'   }
#' 
#' @details
#' The function parses the sumstats_2_file name to extract dataset and gene
#' identifiers, then merges relevant metadata. File names are expected to
#' follow the pattern: "dataset_id_gene_id".
#' 
#' @export
annotate_eQTL_Catalog <- function(coloc_out, annotation_df,
                                  datasets_eQTL_Catalogue) {
  
  # Input validation
  if (!is.data.frame(coloc_out) || nrow(coloc_out) == 0) {
    warning("Empty colocalization results provided")
    return(coloc_out)
  }
  
  if (!"sumstats_2_file" %in% colnames(coloc_out)) {
    stop("sumstats_2_file column not found in coloc_out")
  }
  
  # Extract dataset and gene IDs from filename
  # Expected format: "dataset_id_gene_id"
  coloc_out$dataset_id <- gsub("(.*)_.*", "\\1", coloc_out$sumstats_2_file)
  coloc_out$gene_id <- gsub(".*_(.*)", "\\1", coloc_out$sumstats_2_file)
  
  # Add gene annotations if provided
  if (!is.null(annotation_df) && nrow(annotation_df) > 0) {
    # Determine merge key (assuming gene_id or similar)
    merge_key <- intersect(c("gene_id", "ensembl_gene_id", "ENSG"), 
                           colnames(annotation_df))[1]
    
    if (!is.na(merge_key)) {
      coloc_out <- merge(coloc_out, annotation_df, 
                         by.x = "gene_id", by.y = merge_key,
                         all.x = TRUE)
      
      # Reorder to keep original columns first
      orig_cols <- setdiff(colnames(coloc_out), colnames(annotation_df))
      new_cols <- setdiff(colnames(annotation_df), merge_key)
      coloc_out <- coloc_out[, c(orig_cols, new_cols)]
    }
  }
  
  # Add dataset metadata if provided
  if (!is.null(datasets_eQTL_Catalogue) && nrow(datasets_eQTL_Catalogue) > 0) {
    coloc_out <- merge(coloc_out, datasets_eQTL_Catalogue,
                       by = "dataset_id", all.x = TRUE)
    
    # Reorder to keep original columns first
    orig_cols <- setdiff(colnames(coloc_out), colnames(datasets_eQTL_Catalogue))
    new_cols <- setdiff(colnames(datasets_eQTL_Catalogue), "dataset_id")
    coloc_out <- coloc_out[, c(orig_cols, new_cols)]
  }
  
  return(coloc_out)
}


#' Create argument list for eQTL Catalogue colocalization analysis
#' 
#' @description
#' Generates a properly formatted argument list for running colocalization
#' analysis with eQTL Catalogue datasets.
#' 
#' @param QTS Character vector. Study IDs to include (default: "QTS000015").
#'   Use get_datasets_eQTL_Catalogue() to see available studies.
#' @param quant_method Character string. Quantification method to use
#'   (default: "ge" for gene expression). Options include:
#'   \itemize{
#'     \item "ge": Gene expression
#'     \item "tx": Transcript expression
#'     \item "txrev": Transcript usage
#'     \item "microarray": Microarray gene expression
#'   }
#' 
#' @return A list with components:
#'   \itemize{
#'     \item sumstats_2_files: Dataset identifiers to analyze
#'     \item sumstats_2_function: Function name for data retrieval
#'     \item sumstats_2_type: Summary statistics type ("quant")
#'     \item sumstats_2_sdY: Standard deviation (NA for eQTL data)
#'   }
#' 
#' @details
#' This function simplifies the setup for eQTL colocalization by:
#' \enumerate{
#'   \item Fetching current dataset metadata
#'   \item Filtering by study and quantification method
#'   \item Formatting arguments for genepicoloc_wrapper
#' }
#' 
#' @examples
#' \dontrun{
#' # Get arguments for default study
#' eqtl_args <- make_eQTL_Catalogue_args()
#' 
#' # Get arguments for multiple studies
#' eqtl_args <- make_eQTL_Catalogue_args(
#'   QTS = c("QTS000001", "QTS000002"),
#'   quant_method = "ge"
#' )
#' 
#' # Use in colocalization analysis
#' args_df <- data.frame(
#'   sumstats_2_study = eqtl_args$sumstats_2_files,
#'   sumstats_2_file = eqtl_args$sumstats_2_files,
#'   sumstats_2_function = eqtl_args$sumstats_2_function,
#'   sumstats_2_type = eqtl_args$sumstats_2_type,
#'   sumstats_2_sdY = eqtl_args$sumstats_2_sdY
#' )
#' }
#' 
#' @export
make_eQTL_Catalogue_args <- function(QTS = "QTS000015", quant_method = "ge") {
  
  # Input validation
  if (!is.character(QTS) || length(QTS) == 0) {
    stop("QTS must be a character vector of study IDs")
  }
  
  valid_methods <- c("ge", "tx", "txrev", "microarray", "leafcutter", "aptamer")
  if (!quant_method %in% valid_methods) {
    stop("quant_method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Get dataset metadata
  message("Fetching eQTL Catalogue datasets...")
  datasets_eQTL_Catalogue <- get_datasets_eQTL_Catalogue()
  
  # Filter datasets
  datasets_filtered <- subset(
    datasets_eQTL_Catalogue,
    quant_method == quant_method & study_id %in% QTS
  )
  
  if (nrow(datasets_filtered) == 0) {
    warning("No datasets found for studies: ", paste(QTS, collapse = ", "),
            " with method: ", quant_method)
    
    # Show available studies for this method
    available_studies <- unique(
      datasets_eQTL_Catalogue$study_id[datasets_eQTL_Catalogue$quant_method == quant_method]
    )
    message("Available studies for ", quant_method, ": ", 
            paste(available_studies, collapse = ", "))
    
    return(list(
      sumstats_2_files = character(),
      sumstats_2_function = "query_eQTL_Catalogue",
      sumstats_2_type = "quant",
      sumstats_2_sdY = NA
    ))
  }
  
  message("Found ", nrow(datasets_filtered), " datasets:")
  message(paste("  -", datasets_filtered$dataset_id, 
                "(", datasets_filtered$tissue_label, ")", 
                collapse = "\n"))
  
  # Create argument list
  eQTL_Catalogue_args <- list(
    sumstats_2_files = datasets_filtered$dataset_id,
    sumstats_2_function = "query_eQTL_Catalogue",
    sumstats_2_type = "quant",
    sumstats_2_sdY = NA
  )
  
  return(eQTL_Catalogue_args)
}