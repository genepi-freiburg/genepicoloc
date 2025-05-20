# main wrapper ----
#' Wrapper for parallel colocalization analysis
#' 
#' @param dir_out Output directory for results
#' @param sumstats_1 Primary summary statistics object
#' @param coloc_regions_PASS Regions to analyze
#' @param args_2_df Data frame with secondary dataset parameters
#' @param mc_cores Number of cores for parallel processing
#' @param use_pbmcapply Whether to use pbmcapply for progress tracking
#' @param test_mode Run in test mode (limited execution)
#' @param verbose Print detailed messages
#' @param debug_mode Run sequentially for debugging
#' @return Invisible result of parallel execution
#' @importFrom parallel mcmapply
#' @importFrom pbmcapply pbmcmapply
genepicoloc_wrapper <- function(dir_out, sumstats_1, coloc_regions_PASS, args_2_df,
                                mc_cores = 10, use_pbmcapply = FALSE, test_mode = FALSE, 
                                verbose = FALSE, debug_mode = FALSE) {
  
  # Setup common arguments for all iterations
  shared_args <- list(
    dir_out = dir_out,
    sumstats_1 = sumstats_1,
    coloc_regions_PASS = coloc_regions_PASS,
    test_mode = test_mode,
    verbose = verbose
  )
  
  # Choose parallel function based on parameters
  if (debug_mode) {
    parallel_func <- "mapply"  # Sequential for debugging
  } else if (use_pbmcapply) {
    parallel_func <- "pbmcmapply"  # With progress bar
  } else {
    parallel_func <- "mcmapply"  # Standard parallel
  }
  
  # Prepare arguments from the data frame
  mapply_args <- list(
    FUN = process_sumstats_2,
    study = args_2_df$study,
    sumstats_2_file = args_2_df$sumstats_2_files,
    sumstats_2_function = args_2_df$sumstats_2_function,
    sumstats_2_type = args_2_df$sumstats_2_type,
    sumstats_2_sdY = args_2_df$sumstats_2_sdY,
    MoreArgs = shared_args
  )
  
  # Execute parallel mapping function
  invisible(do.call(parallel_func, mapply_args))
}

# colocalization functions ----
#' Process colocalization analysis for a secondary trait
#'
#' @description Processes a secondary trait summary statistics file, performs
#' colocalization analysis against a primary trait across multiple genomic regions,
#' and saves the results.
#'
#' @param study Study identifier
#' @param sumstats_2_file Path to secondary summary statistics file
#' @param sumstats_2_function Function used to process secondary data
#' @param sumstats_2_type Type of secondary summary statistics ('quant' or 'cc')
#' @param sumstats_2_sdY Standard deviation for secondary trait (for quantitative traits)
#' @param sumstats_1 Primary summary statistics object
#' @param coloc_regions_PASS Data frame with regions to analyze (must contain CHR_var, BP_START_var, BP_STOP_var)
#' @param dir_out Base output directory
#' @param test_mode If TRUE, limit analysis to first two regions (default: FALSE)
#' @param verbose If TRUE, print detailed progress messages (default: FALSE)
#' @param write_excel If TRUE, create Excel files with filtered results (default: FALSE)
#' @param PP_H4_threshold Minimum PP.H4.abf value for filtered results (default: 0.5)
#' @return Invisibly returns path to the output file
#' @importFrom data.table rbindlist
#' @export
process_sumstats_2 <- function(study,
                               sumstats_2_file,
                               sumstats_2_function,
                               sumstats_2_type,
                               sumstats_2_sdY,
                               sumstats_1,
                               coloc_regions_PASS,
                               dir_out,
                               test_mode = FALSE,
                               verbose = FALSE,
                               write_excel = FALSE,
                               PP_H4_threshold = 0.5) {
  
  # Set up output paths and limit regions for test mode
  if (test_mode) coloc_regions_PASS <- coloc_regions_PASS[1:2, ]
  
  study_dir <- file.path(dir_out, study)
  if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE)
  
  file_base <- tools::file_path_sans_ext(basename(sumstats_2_file))
  file_out <- file.path(study_dir, paste0(file_base, "_gc.txt.gz"))
  if (test_mode) file_out <- paste0(file_out, "_test_mode")
  
  # Process secondary summary statistics
  if (verbose) message("Processing secondary summary statistics for ", study)
  
  sumstats_2 <- sumstats_tabix(
    sumstats_file = sumstats_2_file,
    coloc_regions_PASS = coloc_regions_PASS
  )
  
  sumstats_2 <- sumstats_nlog10P(sumstats_2)
  
  sumstats_2 <- sumstats_attr(
    sumstats = sumstats_2,
    sumstats_function = sumstats_2_function,
    sumstats_type = sumstats_2_type,
    sumstats_sdY = sumstats_2_sdY
  )
  
  sumstats_2 <- sumstats_check(sumstats = sumstats_2)
  
  # Run colocalization on each region
  if (verbose) message("Running colocalization across ", nrow(coloc_regions_PASS), " regions")
  
  coloc_results <- mapply(
    run_region,
    CHR_var = coloc_regions_PASS$CHR_var,
    BP_START_var = coloc_regions_PASS$BP_START_var,
    BP_STOP_var = coloc_regions_PASS$BP_STOP_var,
    MoreArgs = list(sumstats_1 = sumstats_1, sumstats_2 = sumstats_2),
    SIMPLIFY = FALSE
  )
  
  # Combine results
  coloc_df <- data.table::rbindlist(coloc_results)
  
  # Write results using the appropriate function
  save_coloc_results(
    coloc_df = coloc_df, 
    file_out = file_out,
    write_excel = write_excel,
    PP_H4_threshold = PP_H4_threshold,
    verbose = verbose
  )
  
  if (verbose) message("Colocalization analysis completed for ", study)
  
  return(invisible(file_out))
}

#' Run colocalization for a specific genomic region
#' 
#' @param sumstats_1 Primary summary statistics object
#' @param sumstats_2 Secondary summary statistics object
#' @param CHR_var Chromosome
#' @param BP_START_var Start position
#' @param BP_STOP_var End position
#' @param min_nlog10P Minimum -log10(p-value) threshold
#' @return Data frame with colocalization results
run_region <- function(sumstats_1, sumstats_2,
                       CHR_var, BP_START_var, BP_STOP_var,
                       min_nlog10P = 5) {
  # Check input objects
  if (!inherits(sumstats_1, "sumstats")) 
    stop("Expected a sumstats object for sumstats_1")
  
  if (!inherits(sumstats_2, "sumstats")) 
    stop("Expected a sumstats object for sumstats_2")
  
  # Extract region-specific data
  region_filter <- function(data) {
    subset(data, CHR == CHR_var & POS > BP_START_var & POS < BP_STOP_var)
  }
  
  sumstats_1_sub <- region_filter(sumstats_1)
  sumstats_2_sub <- region_filter(sumstats_2)
  
  # Get max significance values
  sumstats_1_sub <- sumstats_nlog10P(sumstats_1_sub)
  sumstats_2_sub <- sumstats_nlog10P(sumstats_2_sub)
  sumstats_1_max_nlog10P <- attr(sumstats_1_sub, "max_nlog10P")
  sumstats_2_max_nlog10P <- attr(sumstats_2_sub, "max_nlog10P")
  
  # Create output template with metadata
  result <- create_result_template(
    CHR_var = CHR_var, 
    BP_START_var = BP_START_var,
    BP_STOP_var = BP_STOP_var,
    sumstats_1_file = attr(sumstats_1, "sumstats_file"),
    sumstats_1_max_nlog10P = sumstats_1_max_nlog10P,
    sumstats_2_file = attr(sumstats_2, "sumstats_file"),
    sumstats_2_max_nlog10P = sumstats_2_max_nlog10P,
    tabix = paste0("s1_", attr(sumstats_1, "tabix"), "_s2_", attr(sumstats_2, "tabix")),
    sumstats_1_QC = attr(sumstats_1, "QC"),
    sumstats_2_QC = attr(sumstats_2, "QC")
  )
  
  # Check if primary dataset has significant SNPs
  if (is.na(sumstats_1_max_nlog10P) || sumstats_1_max_nlog10P < min_nlog10P) {
    stop(if (is.na(sumstats_1_max_nlog10P)) 
      "No data found in sumstats_1" 
      else "No significant SNPs in sumstats_1")
  }
  
  # Check if secondary dataset has significant SNPs
  if (is.na(sumstats_2_max_nlog10P) || sumstats_2_max_nlog10P < min_nlog10P) {
    result$coloc <- ifelse(is.na(sumstats_2_max_nlog10P), "no_data_sumstats_2", "no_signif_sumstats_2")
    return(result)
  }
  
  # Check if SNPs overlap between datasets
  if (!any(sumstats_1_sub$Name %in% sumstats_2_sub$Name)) {
    result$coloc <- "no_SNP_intersect"
    return(result)
  }
  
  # Run colocalization
  coloc_output <- run_coloc(
    sumstats_1_df = sumstats_1_sub,
    sumstats_2_df = sumstats_2_sub
  )
  
  # Format and return results
  result$coloc <- "coloc_done"
  result <- format_coloc_output(coloc_output, result)
  
  return(result)
}

#' Run colocalization analysis using coloc.abf
#'
#' @param sumstats_1 Primary summary statistics object
#' @param sumstats_2 Secondary summary statistics object
#' @param silent Whether to suppress coloc.abf output messages
#' @return Results from coloc.abf function
#' @importFrom coloc coloc.abf
#' @export
run_coloc <- function(sumstats_1, sumstats_2, silent = TRUE) {
  # Check input objects
  if (!inherits(sumstats_1, "sumstats")) 
    stop("Expected a sumstats object for sumstats_1")
  
  if (!inherits(sumstats_2, "sumstats")) 
    stop("Expected a sumstats object for sumstats_2")
  
  # Extract attributes needed for colocalization
  sumstats_1_type <- attr(sumstats_1, "sumstats_type")
  sumstats_1_sdY <- attr(sumstats_1, "sumstats_sdY")
  sumstats_2_type <- attr(sumstats_2, "sumstats_type")
  sumstats_2_sdY <- attr(sumstats_2, "sumstats_sdY")
  
  # Format both datasets for coloc.abf input
  dataset_1 <- prepare_coloc_dataset(sumstats_1, sumstats_1_type, sumstats_1_sdY)
  dataset_2 <- prepare_coloc_dataset(sumstats_2, sumstats_2_type, sumstats_2_sdY)
  
  # Run coloc.abf with or without output suppression
  if (silent) {
    # Temporarily redirect output to a file
    tmp_file <- tempfile()
    sink(tmp_file, type = "output")
    on.exit({
      sink()  # Restore output
      unlink(tmp_file)  # Remove temporary file
    })
    
    # Run coloc with warnings suppressed
    coloc_output <- suppressWarnings(coloc.abf(dataset1 = dataset_1, dataset2 = dataset_2))
  } else {
    # Run coloc normally
    coloc_output <- coloc.abf(dataset1 = dataset_1, dataset2 = dataset_2)
  }
  
  return(coloc_output)
}

#' Format sumstats data.frame for the coloc.abf input
#'
#' @param sumstats Summary statistics object
#' @param sumstats_type Type of summary statistics ('quant' or 'cc')
#' @param sumstats_sdY Standard deviation for quantitative trait
#' @return List formatted for coloc.abf input
#' @export
prepare_coloc_dataset <- function(sumstats, sumstats_type, sumstats_sdY) {
  # Create base dataset with required fields
  dataset <- list(
    beta = sumstats$BETA,
    varbeta = (sumstats$SE)^2,
    snp = sumstats$Name,
    type = sumstats_type
  )
  
  # Add additional fields based on trait type
  if (sumstats_type == "quant") {
    if (!is.na(sumstats_sdY)) {
      dataset$sdY <- sumstats_sdY
    } else {
      # Fallback if sdY not provided
      dataset$MAF <- sumstats$AF
      dataset$N <- sumstats$N
    }
  } else if (sumstats_type == "cc") {
    # For case-control studies
    # Note: May need to add fields here depending on coloc.abf requirements
    dataset$s <- NULL  # Uncomment and add proportion of cases if available
    dataset$N <- sumstats$N  # Total sample size
  }
  
  return(dataset)
}


#' Format the output of colocalization analysis
#'
#' @param coloc_output Output from coloc::coloc.abf function
#' @param coloc_template Results template from create_result_template
#' @param n_top_snps Number of SNPs with highest PP.H4 to include
#' @return Formatted results of colocalization analysis
#' @export
format_coloc_output <- function(coloc_output, coloc_template, n_top_snps = 5) {
  # Extract top SNPs by posterior probability
  top_snps <- coloc_output$results[, c("snp", "SNP.PP.H4")]
  top_snps <- head(top_snps[order(top_snps$SNP.PP.H4, decreasing = TRUE), ], n_top_snps)
  
  # Format SNP names and probabilities
  snp_names <- paste(top_snps$snp, collapse = ", ")
  snp_probs <- top_snps$SNP.PP.H4
  formatted_probs <- ifelse(
    snp_probs < 0.01,
    format(snp_probs, scientific = TRUE, digits = 2),
    sprintf("%.2f", round(snp_probs, 2))
  )
  formatted_probs <- paste(formatted_probs, collapse = ", ")
  
  # Add results to template
  coloc_template$Top_coloc_SNP <- snp_names
  coloc_template$Top_coloc_SNP.PP.H4 <- formatted_probs
  coloc_template$priors <- paste(coloc_output$priors, collapse = ", ")
  coloc_template$nsnps <- coloc_output$summary["nsnps"]
  
  # Add posterior probabilities for each hypothesis
  for (h in c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")) {
    coloc_template[[h]] <- coloc_output$summary[h]
  }
  
  return(coloc_template)
}

#' Create a template for colocalization results
#'
#' @param CHR_var Chromosome
#' @param BP_START_var Start position
#' @param BP_STOP_var End position
#' @param sumstats_1_file Primary summary statistics file
#' @param sumstats_1_max_nlog10P Maximum -log10(p-value) in primary data
#' @param sumstats_2_file Secondary summary statistics file
#' @param sumstats_2_max_nlog10P Maximum -log10(p-value) in secondary data
#' @param nsnps Number of SNPs in analysis
#' @param tabix Tabix status
#' @param sumstats_1_QC QC status of primary data
#' @param sumstats_2_QC QC status of secondary data
#' @param coloc Colocalization status
#' @return Data frame template for colocalization results
#' @export
create_result_template <- function(CHR_var, BP_START_var, BP_STOP_var,
                                   sumstats_1_file, sumstats_1_max_nlog10P,
                                   sumstats_2_file, sumstats_2_max_nlog10P, 
                                   nsnps = NA, tabix = NA, 
                                   sumstats_1_QC = NA, sumstats_2_QC = NA, coloc = NA) {
  
  # Create results data frame with NA values for fields to be filled later
  data.frame(
    # Region information
    CHR_var = CHR_var,
    BP_START_var = BP_START_var,
    BP_STOP_var = BP_STOP_var,
    
    # Summary statistics metadata
    sumstats_1_file = sumstats_1_file,
    sumstats_1_max_nlog10P = sumstats_1_max_nlog10P,
    sumstats_2_file = sumstats_2_file,
    sumstats_2_max_nlog10P = sumstats_2_max_nlog10P,
    
    # Analysis metadata
    nsnps = nsnps,
    tabix = tabix,
    sumstats_1_QC = sumstats_1_QC,
    sumstats_2_QC = sumstats_2_QC,
    coloc = coloc,
    
    # Colocalization results (to be filled later)
    PP.H0.abf = NA, 
    PP.H1.abf = NA, 
    PP.H2.abf = NA,
    PP.H3.abf = NA, 
    PP.H4.abf = NA,
    Top_coloc_SNP = NA, 
    Top_coloc_SNP.PP.H4 = NA, 
    priors = NA,
    
    stringsAsFactors = FALSE
  )
}

#' Save colocalization results to files
#'
#' @description Saves colocalization results to text and optionally RDS/Excel files.
#'
#' @param coloc_df Data frame containing colocalization results
#' @param file_out Path for the main compressed text output file
#' @param write_excel Whether to create Excel files (default: FALSE)
#' @param PP_H4_threshold Threshold for filtering results by PP.H4.abf (default: 0.5)
#' @param verbose Whether to print progress messages (default: FALSE)
#' @return Invisibly returns path to main output file
#' @importFrom data.table fwrite
#' @importFrom writexl write_xlsx
#' @export
save_coloc_results <- function(coloc_df,
                               file_out,
                               write_excel = FALSE,
                               PP_H4_threshold = 0.5,
                               verbose = FALSE) {
  
  # Check inputs
  if (is.null(coloc_df) || nrow(coloc_df) == 0) {
    warning("Empty colocalization results provided")
    return(invisible(NULL))
  }
  
  # Extract directory and base name for additional files
  study_dir <- dirname(file_out)
  base_name <- tools::file_path_sans_ext(basename(file_out))
  
  # Save RDS version
  rds_file <- file.path(study_dir, paste0(base_name, ".RDS"))
  saveRDS(coloc_df, rds_file)
  if (verbose) message("Unfiltered results saved to: ", rds_file)
  
  # Optionally create Excel files
  if (write_excel) {
    # Filter results by PP.H4.abf threshold
    coloc_df_filt <- coloc_df[!is.na(coloc_df$PP.H4.abf) & coloc_df$PP.H4.abf >= PP_H4_threshold, ]
    
    # Create Excel file
    xlsx_file <- file.path(study_dir, paste0(base_name, "_filt.xlsx"))
    writexl::write_xlsx(coloc_df_filt, xlsx_file)
    if (verbose) message("Filtered results saved to: ", xlsx_file)
  }
  
  return(invisible(file_out))
}