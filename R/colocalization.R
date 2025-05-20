# coloc ----
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

#' Process colocalization for a secondary dataset across multiple regions
#' 
#' @param study Study identifier
#' @param sumstats_2_file Path to secondary summary statistics file
#' @param sumstats_2_function Function used to process secondary data
#' @param sumstats_2_type Type of secondary summary statistics
#' @param sumstats_2_sdY Standard deviation for secondary trait
#' @param sumstats_1 Primary summary statistics object
#' @param coloc_regions_PASS Regions to analyze
#' @param dir_out Base output directory
#' @param test_mode Run in test mode (limited regions)
#' @param verbose Print detailed messages
#' @return NULL (results are written to file)
process_sumstats_2 <- function(study,
                               sumstats_2_file,
                               sumstats_2_function,
                               sumstats_2_type,
                               sumstats_2_sdY,
                               sumstats_1,
                               coloc_regions_PASS,
                               dir_out,
                               test_mode = FALSE,
                               verbose = FALSE) {
  
  # Set up output paths and limit regions for test mode
  if (test_mode) coloc_regions_PASS <- coloc_regions_PASS[1:2, ]
  
  study_dir <- file.path(dir_out, study)
  if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE)
  
  file_out <- file.path(study_dir, paste0(basename(sumstats_2_file), "_gc.txt.gz"))
  if (test_mode) file_out <- paste0(file_out, "_test_mode")
  
  # Process secondary summary statistics
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
  coloc_results <- mapply(
    run_region,
    CHR_var = coloc_regions_PASS$CHR_var,
    BP_START_var = coloc_regions_PASS$BP_START_var,
    BP_STOP_var = coloc_regions_PASS$BP_STOP_var,
    MoreArgs = list(sumstats_1 = sumstats_1, sumstats_2 = sumstats_2),
    SIMPLIFY = FALSE
  )
  
  # Combine results and write to file
  coloc_df <- data.table::rbindlist(coloc_results)
  data.table::fwrite(coloc_df, file_out)
  
  if (verbose) message("Results written to ", file_out)
  
  return(NULL)
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
