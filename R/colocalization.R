# main wrapper ----
#' Wrapper for parallel colocalization analysis
#' 
#' @param dir_out Output directory for results
#' @param sumstats_1 Primary summary statistics object
#' @param coloc_regions_PASS Regions to analyze
#' @param args_df Data frame with secondary dataset parameters
#' @param mc_cores Number of cores for parallel processing
#' @param use_pbmcapply Whether to use pbmcapply for progress tracking
#' @param test_mode Run in test mode (limited execution)
#' @param verbose Print detailed messages
#' @param debug_mode Run sequentially for debugging
#' @return Invisible result of parallel execution
#' @importFrom parallel mcmapply
#' @importFrom pbmcapply pbmcmapply
genepicoloc_wrapper <- function(dir_out,
                                sumstats_1_form,
                                args_df,
                                mc_cores = 10,
                                use_pbmcapply = FALSE,
                                test_mode = FALSE, 
                                verbose = FALSE,
                                debug_mode = FALSE) {
  # Limit regions for test mode
  if (test_mode) {
    args_df <- args_df[, .SD[1], by = sumstats_2_study]
    coloc_regions_PASS <- attr(sumstats_1_form, "coloc_regions_PASS")[1,]
    dir_out <- paste0(dir_out, "_test_mode")
  }
  
  message("Starting genepicoloc: ", nrow(args_df), " sumstats to be processed.")

  # Create folder structure
  if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
  if (verbose) message("Output will be written to ", dir_out)
  
  # Setup common arguments for all iterations
  shared_args <- list(
    dir_out = dir_out,
    sumstats_1 = sumstats_1_form,
    verbose = verbose
  )
  
  # Prepare arguments from the data frame
  mapply_args <- list(
    FUN = process_sumstats_form,
    sumstats_2_study = args_df$sumstats_2_study,
    sumstats_2_file = args_df$sumstats_2_file,
    sumstats_2_function = args_df$sumstats_2_function,
    sumstats_2_type = args_df$sumstats_2_type,
    sumstats_2_sdY = args_df$sumstats_2_sdY,
    MoreArgs = shared_args
  )
  
  # Choose parallel function based on parameters
  if (debug_mode) {
    parallel_func <- "mapply"  # Sequential for debugging
  } else if (use_pbmcapply) {
    parallel_func <- "pbmcmapply"  # With progress bar
    mapply_args <- c(mapply_args, mc.cores=mc_cores)
  } else {
    parallel_func <- "mcmapply"  # Standard parallel
    mapply_args <- c(mapply_args, mc.cores=mc_cores)
  }
  
  # save args_df
  args_df_file <- paste0(dir_out, "/args_df.csv.gz")
  if (verbose) message("Saving args_df_file to ", args_df_file)
  data.table::fwrite(args_df, args_df_file)

  # Execute parallel mapping function
  invisible(do.call(parallel_func, mapply_args))
}

# colocalization functions ----
format_sumstats_1 <- function(coloc_regions_PASS,
                              sumstats_1_function,
                              sumstats_1_file,
                              sumstats_1_type,
                              sumstats_1_sdY) {
  sumstats_1_raw <- retrieve_sumstats_raw(
    sumstats_function = sumstats_1_function,
    sumstats_file = sumstats_1_file,
    coloc_regions_PASS = coloc_regions_PASS
  )
  sumstats_1_form <- format_sumstats(
    sumstats = sumstats_1_raw,
    sumstats_type = sumstats_1_type,
    sumstats_sdY = sumstats_1_sdY
  )
}

#' Process colocalization analysis for a secondary trait
#'
#' @description Processes a secondary trait summary statistics file, performs
#' colocalization analysis against a primary trait across multiple genomic regions,
#' and saves the results.
#'
#' @param sumstats_2_study Study identifier
#' @param sumstats_2_file Path to secondary summary statistics file
#' @param sumstats_2_function Function used to process secondary data
#' @param sumstats_2_type Type of secondary summary statistics ('quant' or 'cc')
#' @param sumstats_2_sdY Standard deviation for secondary trait (for quantitative traits)
#' @param sumstats_1 Primary summary statistics object
#' @param coloc_regions_PASS Data frame with regions to analyze (must contain CHR_var, BP_START_var, BP_STOP_var)
#' @param dir_out Base output directory
#' @param verbose If TRUE, print detailed progress messages (default: FALSE)
#' @param write_excel If TRUE, create Excel files with filtered results (default: FALSE)
#' @param PP_H4_threshold Minimum PP.H4.abf value for filtered results (default: 0.5)
#' @return Invisibly returns path to the output file
#' @importFrom data.table rbindlist
#' @export
format_sumstats_2 <- function(coloc_regions_PASS,
                              sumstats_2_function,
                              sumstats_2_file,
                              sumstats_2_type,
                              sumstats_2_sdY,
                              sumstats_2_study) {
  sumstats_2_raw <- retrieve_sumstats_raw(
    sumstats_function = sumstats_2_function,
    sumstats_file = sumstats_2_file,
    coloc_regions_PASS = coloc_regions_PASS
  )

  # check that input is valid
  if (is.null(attr(sumstats_2_raw, "sumstats_pheno"))) {
    stop("sumstats_pheno attribute in sumstats_2_raw not found")
  }
  if (!attr(sumstats_2_raw, "sumstats_pheno") %in% c("single", "multiple")) {
    stop("sumstats_pheno should be either single or multiple")
  }
  
  if (attr(sumstats_2_raw, "sumstats_pheno") == "multiple") {
    if (! "Phenotype" %in% colnames(sumstats_2_raw)) {
      stop("'Phenotype' column in sumstats_2_raw not found, ",
           "required for 'multiple' sumstats_pheno attribute")
    }
    # handle case with tabix_ok_no_data - tmp
    # TODO implement proper handling using sumstats_2_raw class
    Phenotypes <- unique(sumstats_2_raw$Phenotype)
    if (length(Phenotypes) == 0) {
      attr(sumstats_2_raw, "sumstats_pheno") <- "single"
      sumstats_2_raw$Phenotype <- NULL
      sumstats_2_form <- list(format_sumstats(
        sumstats = sumstats_2_raw,
        sumstats_type = sumstats_2_type,
        sumstats_sdY = sumstats_2_sdY))
    } else {
      sumstats_2_form <- 
        lapply(Phenotypes, function(x) {
          sumstats_2_pheno <- subset(sumstats_2_raw, Phenotype == x)
          attr(sumstats_2_pheno, "sumstats_file") <- 
            paste0(attr(sumstats_2_pheno, "sumstats_file"), "_", x)
          attr(sumstats_2_pheno, "sumstats_pheno") <- "single"
          sumstats_2_pheno$Phenotype <- NULL
          sumstats_2_pheno <- format_sumstats(
            sumstats = sumstats_2_pheno,
            sumstats_type = sumstats_2_type,
            sumstats_sdY = sumstats_2_sdY)
        })
    }
    attr(sumstats_2_form, "sumstats_pheno") <- "multiple"
    attr(sumstats_2_form, "sumstats_file") <- sumstats_2_file
  } else {
    sumstats_2_form <- format_sumstats(
      sumstats = sumstats_2_raw,
      sumstats_type = sumstats_2_type,
      sumstats_sdY = sumstats_2_sdY)
    attr(sumstats_2_form, "sumstats_pheno") <- "single"
  }
  class(sumstats_2_form) <- c("sumstats_2_form", class(sumstats_2_form))
  return(sumstats_2_form)
}


process_sumstats_form <- function(dir_out,
                                  sumstats_1_form,
                                  sumstats_2_function,
                                  sumstats_2_file,
                                  sumstats_2_type,
                                  sumstats_2_sdY,
                                  sumstats_2_study,
                                  write_excel = FALSE,
                                  PP_H4_threshold = 0.5,
                                  verbose = TRUE) {
  
  # check input, consider creating class "sumstats_1_form"
  if (is.null(attr(sumstats_1_form, "coloc_regions_PASS"))) {
    stop("'coloc_regions_PASS' attribute not found in sumstats_1_form object")
  } else {
    coloc_regions_PASS <- attr(sumstats_1_form, "coloc_regions_PASS")
  }
  
  sumstats_2_form <- format_sumstats_2(
    coloc_regions_PASS = coloc_regions_PASS,
    sumstats_2_function = sumstats_2_function,
    sumstats_2_file = sumstats_2_file,
    sumstats_2_type = sumstats_2_type,
    sumstats_2_sdY = sumstats_2_sdY,
    sumstats_2_study = sumstats_2_study
  )
  
  # TODO create a helper to check for sumstats_2_form
  # check that is has sumstats_pheno attribute
  # and the values should be only single or multiple
  if (!inherits(sumstats_2_form, "sumstats_2_form")) {
    stop("sumstats_2_form should have 'sumstats_2_form' class")
  }
  
  # Run colocalization on each region
  if (verbose) message("Running colocalization across ", nrow(coloc_regions_PASS), " regions")

  file_out <- set_output_path(
    sumstats_2_file = attr(sumstats_2_form, "sumstats_file"),
    dir_out = dir_out,
    sumstats_2_study = sumstats_2_study
  )
  if (verbose) message("Output file: ", file_out, "\nStarting coloc.abf")
  
  # now there are two options, if sumstats_pheno == "single",
  # i.e., there is only one phenotype
  if (attr(sumstats_2_form, "sumstats_pheno") == "single") {
    coloc_results <- run_single_phenotype(
      sumstats_1 = sumstats_1_form,
      sumstats_2 = sumstats_2_form,
      coloc_regions_PASS = coloc_regions_PASS)
    # Combine all regions (there is just a single phenotype)
    coloc_results <- data.table::rbindlist(coloc_results)
  } else { # if there are multiple phenotypes
    coloc_results <- lapply(sumstats_2_form, function(sumstats_2) {
      coloc_results <- run_single_phenotype(
        sumstats_1 = sumstats_1_form,
        sumstats_2 = sumstats_2,
        coloc_regions_PASS = coloc_regions_PASS)
      # Combine all regions across each phenotype
      coloc_results <- data.table::rbindlist(coloc_results)
    })
    # Combine all phenotypes
    coloc_results <- data.table::rbindlist(coloc_results)
  }
  
  
  # Write results using the appropriate function
  save_coloc_results(
    coloc_results = coloc_results, 
    file_out = file_out,
    write_excel = write_excel,
    PP_H4_threshold = PP_H4_threshold
  )
  
  if (verbose) message("Colocalization analysis completed for ", sumstats_2_study)
  
  return(invisible(file_out))
  
}


# Set up output paths
set_output_path <- function(dir_out,
                            sumstats_2_file,
                            sumstats_2_study) {
  
  study_dir <- file.path(dir_out, sumstats_2_study)
  if (!dir.exists(study_dir)) dir.create(study_dir, recursive = TRUE)
  
  file_base <- basename(sumstats_2_file)
  file_out <- file.path(study_dir, file_base)
  
  return(file_out)
}

#
run_single_phenotype <- function(sumstats_1, sumstats_2, coloc_regions_PASS) {
  MoreArgs <- list(sumstats_1 = sumstats_1,
                   sumstats_2 = sumstats_2)
  coloc_results <- with(coloc_regions_PASS, 
                        mapply(
                          FUN = run_region,
                          CHR_var = CHR_var,
                          BP_START_var = BP_START_var,
                          BP_STOP_var = BP_STOP_var,
                          MoreArgs = MoreArgs,
                          SIMPLIFY = FALSE
                        ))
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

  # Create output template with metadata
  result <- create_result_template(
    CHR_var = CHR_var, 
    BP_START_var = BP_START_var,
    BP_STOP_var = BP_STOP_var,
    sumstats_1_file = attr(sumstats_1, "sumstats_file"),
    sumstats_2_file = attr(sumstats_2, "sumstats_file"),
    sumstats_1_tabix = attr(sumstats_1, "tabix"),
    sumstats_2_tabix = attr(sumstats_2, "tabix"),
    sumstats_1_QC = attr(sumstats_1, "QC"),
    sumstats_2_QC = attr(sumstats_2, "QC")
  )
  
  if (attr(sumstats_2, "tabix") == "tabix_failed" || 
      attr(sumstats_2, "tabix") == "tabix_ok_no_data") {
    return(result)
  }
  
  # Extract region-specific data
  region_filter <- function(data) {
    subset(data, CHR == CHR_var & POS > BP_START_var & POS < BP_STOP_var)
  }
  sumstats_1_sub <- region_filter(sumstats_1)
  if (nrow(sumstats_1_sub) == 0) stop ("No data found in sumstats_1")
  sumstats_2_sub <- region_filter(sumstats_2)
  
  # Get max significance values
  sumstats_1_sub <- set_max_nlog10P(sumstats_1_sub)
  sumstats_2_sub <- set_max_nlog10P(sumstats_2_sub)
  sumstats_1_max_nlog10P <- attr(sumstats_1_sub, "max_nlog10P")
  sumstats_2_max_nlog10P <- attr(sumstats_2_sub, "max_nlog10P")
  result$sumstats_1_max_nlog10P <- sumstats_1_max_nlog10P
  result$sumstats_2_max_nlog10P <- sumstats_2_max_nlog10P
  
  # Check if primary dataset has significant SNPs
  if (is.na(sumstats_1_max_nlog10P) || sumstats_1_max_nlog10P < min_nlog10P) {
    stop("No significant SNPs in sumstats_1")
  }
  
  # Check if secondary dataset has significant SNPs
  if (is.na(sumstats_2_max_nlog10P) || sumstats_2_max_nlog10P < min_nlog10P) {
    result$coloc <- "no_signif_sumstats_2"
    return(result)
  }
  
  # Check if SNPs overlap between datasets
  if (!any(sumstats_1_sub$Name %in% sumstats_2_sub$Name)) {
    result$coloc <- "no_SNP_intersect"
    return(result)
  }
  
  # beta effect directionality module - simple implementation first
  directionality_df <- get_directionality(
    sumstats_1 = sumstats_1_sub,
    sumstats_2 = sumstats_2_sub
  )
  for (i in colnames(directionality_df)) result[[i]] <- directionality_df[[i]]

  # Run colocalization
  coloc_output <- run_coloc(
    sumstats_1 = sumstats_1_sub,
    sumstats_2 = sumstats_2_sub
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

get_directionality <- function(sumstats_1, sumstats_2) {
  common_names <- intersect(sumstats_1$Name, sumstats_2$Name)
  sumstats_1 <- subset(sumstats_1, Name %in% common_names)
  sumstats_2 <- subset(sumstats_2, Name %in% common_names)
  
  index_1 <- sumstats_1[which.max(sumstats_1$nlog10P),]
  index_2 <- subset(sumstats_2, Name  == index_1$Name)
  
  directionality_df <- data.frame(
    sumstats_1_ind_Name = index_1$Name,
    sumstats_1_ind_A1 = index_1$A1,
    sumstats_1_ind_A2 = index_1$A2,
    sumstats_1_ind_nlog10P = index_1$nlog10P,
    sumstats_1_ind_BETA = index_1$BETA,
    sumstats_1_ind_BETA_sign = sign(index_1$BETA),
    sumstats_2_ind_A1 = index_2$A1,
    sumstats_2_ind_A2 = index_2$A2,
    sumstats_2_ind_nlog10P = index_2$nlog10P,
    sumstats_2_ind_BETA = index_2$BETA,
    sumstats_2_ind_BETA_sign = sign(index_2$BETA)
  )
  directionality_df <- detect_directionality(ind_df=directionality_df)
  return(directionality_df)
}

detect_directionality <- function(ind_df) {
  # case 1, same A1
  if (ind_df$sumstats_1_ind_A1 == ind_df$sumstats_2_ind_A1) {
    directionality <- ifelse(
      test = ind_df$sumstats_1_ind_BETA_sign == ind_df$sumstats_2_ind_BETA_sign,
      yes = 1,
      no = -1
    )
  } else {
    # try to compare with A2
    if (ind_df$sumstats_1_ind_A1 == ind_df$sumstats_2_ind_A2) {
      # if yes, then compare opposite directionality
      directionality <- ifelse(
        test = ind_df$sumstats_1_ind_BETA_sign == ind_df$sumstats_2_ind_BETA_sign,
        yes = -1,
        no = 1
      )
    } else {
      # directionality undefined - need to handle this special case
      directionality <- NA
    }
  }
  ind_df$directionality <- directionality
  return(ind_df)
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
  } 
  # TODO: Add support for p-values depending on coloc.abf requirements
  # https://github.com/chr1swallace/coloc/blob/main/R/claudia.R#L222
  
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
create_result_template <- function(CHR_var,
                                   BP_START_var,
                                   BP_STOP_var,
                                   sumstats_1_file,
                                   sumstats_2_file,
                                   sumstats_1_tabix, 
                                   sumstats_2_tabix, 
                                   sumstats_1_max_nlog10P = NA,
                                   sumstats_2_max_nlog10P = NA, 
                                   nsnps = NA,
                                   sumstats_1_QC = NA,
                                   sumstats_2_QC = NA,
                                   coloc = NA) {
  
  # Create results data frame with NA values for fields to be filled later
  data.frame(
    # Region information
    CHR_var = CHR_var,
    BP_START_var = BP_START_var,
    BP_STOP_var = BP_STOP_var,
    
    # Summary statistics metadata
    sumstats_1_file = sumstats_1_file,
    sumstats_1_tabix = sumstats_1_tabix,
    sumstats_1_max_nlog10P = sumstats_1_max_nlog10P,
    sumstats_2_file = sumstats_2_file,
    sumstats_2_max_nlog10P = sumstats_2_max_nlog10P,
    sumstats_2_tabix = sumstats_2_tabix,
    
    # Analysis metadata
    nsnps = nsnps,
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
    
    # directionality module
    sumstats_1_ind_Name = NA,
    sumstats_1_ind_A1 = NA,
    sumstats_1_ind_A2 = NA,
    sumstats_1_ind_nlog10P = NA,
    sumstats_1_ind_BETA = NA,
    sumstats_1_ind_BETA_sign = NA,
    sumstats_2_ind_A1 = NA,
    sumstats_2_ind_A2 = NA,
    sumstats_2_ind_nlog10P = NA,
    sumstats_2_ind_BETA = NA,
    sumstats_2_ind_BETA_sign = NA,
    directionality = NA,
    
    # stringsAsFactors
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
save_coloc_results <- function(coloc_results,
                               file_out,
                               write_excel = FALSE,
                               PP_H4_threshold = 0.5,
                               verbose = TRUE,
                               sans_ext = TRUE) {
  
  # Check inputs
  if (is.null(coloc_results) || nrow(coloc_results) == 0) {
    warning("Empty colocalization results provided")
    return(invisible(NULL))
  }
  
  # Save RDS version
  if (sans_ext) file_out <- tools::file_path_sans_ext(file_out)
  csv_file <- paste0(file_out, "_ge.csv.gz")
  data.table::fwrite(coloc_results, csv_file)
  if (verbose) message("Unfiltered results saved to: ", csv_file)
  
  # Optionally create Excel files
  if (write_excel) {
    # Filter results by PP.H4.abf threshold
    coloc_df_filt <- coloc_results[!is.na(coloc_results$PP.H4.abf) & coloc_results$PP.H4.abf >= PP_H4_threshold, ]
    
    # Create Excel file
    xlsx_file <- paste0(file_out, "_filt.xlsx")
    writexl::write_xlsx(coloc_df_filt, xlsx_file)
    if (verbose) message("Filtered results saved to: ", xlsx_file)
  }
  
  return(invisible(file_out))
}