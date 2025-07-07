#' Wrapper for parallel colocalization analysis
#' 
#' @description
#' Orchestrates parallel colocalization analysis across multiple secondary datasets
#' against a primary dataset. Supports different parallelization strategies and 
#' includes options for debugging and testing.
#' 
#' @param dir_out Character string. Output directory path where results will be saved.
#'   Directory will be created if it doesn't exist.
#' @param sumstats_1_form A formatted sumstats object containing the primary dataset.
#'   Must have a "coloc_regions_PASS" attribute containing regions to analyze.
#' @param args_df A data.table/data.frame with columns:
#'   \itemize{
#'     \item sumstats_2_study: Study identifier for secondary dataset
#'     \item sumstats_2_file: Path to secondary summary statistics file
#'     \item sumstats_2_function: Function name to process secondary data
#'     \item sumstats_2_type: Type of summary statistics ('quant' or 'cc')
#'     \item sumstats_2_sdY: Standard deviation for quantitative traits
#'   }
#' @param mc_cores Integer. Number of cores for parallel processing (default: 10).
#'   Ignored when debug_mode = TRUE.
#' @param test_mode Logical. Run in test mode with limited execution (default: FALSE).
#'   Processes only the first region and first entry per study.
#' @param verbose Logical. Print detailed progress messages (default: TRUE).
#' @param debug_mode Logical. Run sequentially for debugging (default: FALSE).
#'   Overrides mc_cores setting.
#' @param progress_interval Numeric. How often to check progress in seconds (default: 5).
#'   Only used when verbose = TRUE.
#' @param max_regions_per_job Integer. Maximum number of regions per job (default: 100).
#'   Jobs will be split if coloc_regions_PASS exceeds this number.
#' @param collect_output Logical. Collect and consolidate output from job subfolders (default: TRUE).
#'   Creates consolidated files in main directory and archives job folders.
#' 
#' @return Invisibly returns NULL. Results are written to disk in the specified
#'   output directory.
#' 
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates output directory structure
#'   \item Saves input arguments to args_df.csv.gz for reproducibility
#'   \item Runs colocalization analysis in parallel (or sequentially if debugging)
#'   \item Each secondary dataset is processed by \code{\link{process_sumstats_form}}
#'   \item If collect_output = TRUE, consolidates results from job subdirectories
#' }
#' 
#' In test mode, the output directory name is appended with "_test_mode" and
#' only the first region and first dataset per study are processed.
#' 
#' When coloc_regions_PASS has more than max_regions_per_job rows, the analysis
#' is automatically split into multiple jobs to manage computational load.
#' 
#' @importFrom parallel mcmapply
#' @importFrom data.table fwrite fread rbindlist .SD
#' 
#' @examples
#' \dontrun{
#' # Prepare arguments data frame
#' args_df <- data.frame(
#'   sumstats_2_study = c("Study1", "Study2"),
#'   sumstats_2_file = c("path/to/study1.gz", "path/to/study2.gz"),
#'   sumstats_2_function = "retrieve_sumstats_tabix",
#'   sumstats_2_type = "quant",
#'   sumstats_2_sdY = NA
#' )
#' 
#' # Run colocalization
#' genepicoloc_wrapper(
#'   dir_out = "results/coloc",
#'   sumstats_1_form = my_sumstats,
#'   args_df = args_df,
#'   mc_cores = 4,
#'   verbose = TRUE
#' )
#' }
#' 
#' @export
genepicoloc_wrapper <- function(dir_out,
                                sumstats_1_form,
                                args_df,
                                mc_cores = 10,
                                test_mode = FALSE, 
                                verbose = TRUE,
                                debug_mode = FALSE,
                                progress_interval = 5,
                                max_regions_per_job = 100,
                                collect_output = TRUE) {
  
  # Input validation
  if (!inherits(sumstats_1_form, "sumstats")) {
    stop("sumstats_1_form must be a sumstats object")
  }
  if (is.null(attr(sumstats_1_form, "coloc_regions_PASS"))) {
    stop("sumstats_1_form must have a 'coloc_regions_PASS' attribute")
  }
  
  # Get coloc_regions_PASS
  coloc_regions_PASS <- attr(sumstats_1_form, "coloc_regions_PASS")
  n_regions <- nrow(coloc_regions_PASS)
  
  # Configure test mode
  if (test_mode) {
    # Process only first dataset per study
    args_df <- data.table::data.table(args_df[, .SD[1], by = sumstats_2_study])
    # Process only first region
    coloc_regions_PASS <- coloc_regions_PASS[1,, drop = FALSE]
    attr(sumstats_1_form, "coloc_regions_PASS") <- coloc_regions_PASS
    # Modify output directory name
    dir_out <- paste0(dir_out, "_test_mode")
    n_regions <- 1  # Update region count for test mode
  }
  
  # Calculate number of jobs needed
  n_jobs <- ceiling(n_regions / max_regions_per_job)
  
  if (verbose) {
    message(sprintf("coloc_regions_PASS has %d rows. Splitting into %d jobs (max %d regions per job).",
                    n_regions, n_jobs, max_regions_per_job))
  }
  
  # Run jobs sequentially
  all_results <- list()
  
  for (job_idx in 1:n_jobs) {
    # Calculate region indices for this job
    start_idx <- (job_idx - 1) * max_regions_per_job + 1
    end_idx <- min(job_idx * max_regions_per_job, n_regions)
    
    if (verbose) {
      message(sprintf("\n=== Running job %d/%d (regions %d-%d) ===", 
                      job_idx, n_jobs, start_idx, end_idx))
    }
    
    # Create a copy of sumstats_1_form with subset of regions
    sumstats_1_form_subset <- sumstats_1_form
    attr(sumstats_1_form_subset, "coloc_regions_PASS") <- 
      coloc_regions_PASS[start_idx:end_idx,, drop = FALSE]
    
    # Create job-specific output directory
    dir_out_job <- file.path(dir_out, sprintf("job_%03d", job_idx))
    
    # Run the analysis for this subset
    job_results <- genepicoloc_wrapper_single_job(
      dir_out = dir_out_job,
      sumstats_1_form = sumstats_1_form_subset,
      args_df = args_df,
      mc_cores = mc_cores,
      verbose = verbose,
      debug_mode = debug_mode,
      progress_interval = progress_interval,
      is_subjob = TRUE,
      job_info = list(job_idx = job_idx, n_jobs = n_jobs, 
                      start_idx = start_idx, end_idx = end_idx)
    )
    
    all_results[[job_idx]] <- job_results
  }
  
  if (verbose) {
    message(sprintf("\nAll %d jobs completed successfully!", n_jobs))
    message(sprintf("Results are organized in subdirectories under: %s", dir_out))
  }
  
  # Save job splitting information
  job_info_file <- file.path(dir_out, "job_splitting_info.txt")
  writeLines(c(
    sprintf("Total regions: %d", n_regions),
    sprintf("Regions per job: %d", max_regions_per_job),
    sprintf("Number of jobs: %d", n_jobs),
    sprintf("Job directories: job_001 to job_%03d", n_jobs)
  ), job_info_file)
  
  # Collect and consolidate output if requested
  if (collect_output) {
    if (verbose) message("\n=== Collecting and consolidating output files ===")
    collect_job_outputs(dir_out, n_jobs, verbose)
  }
  
  return(invisible(all_results))
}

#' Collect and consolidate output files from job subdirectories
#' 
#' @param dir_out Main output directory containing job subdirectories
#' @param n_jobs Number of job subdirectories to process
#' @param verbose Logical. Print progress messages
#' 
#' @importFrom data.table fread fwrite rbindlist
collect_job_outputs <- function(dir_out, n_jobs, verbose = TRUE) {
  
  # Get all job directories
  job_dirs <- file.path(dir_out, sprintf("job_%03d", 1:n_jobs))
  
  # Check that all job directories exist
  if (!all(dir.exists(job_dirs))) {
    stop("Not all job directories exist. Please check the output.")
  }
  
  # 1. Find all unique output files (excluding args_df.csv.gz and job_info.txt)
  all_files <- list()
  for (job_dir in job_dirs) {
    files <- list.files(job_dir, pattern = "\\.csv\\.gz$", 
                        recursive = TRUE, full.names = FALSE)
    # Exclude args_df.csv.gz from the main level
    files <- files[!files %in% c("args_df.csv.gz")]
    all_files[[job_dir]] <- files
  }
  
  # Get unique file paths (relative to job directory)
  unique_files <- unique(unlist(all_files))
  
  if (verbose) {
    message(sprintf("Found %d unique output files to consolidate", length(unique_files)))
  }
  
  # 2. Consolidate each unique file
  for (file_rel_path in unique_files) {
    if (verbose) message(sprintf("Consolidating: %s", file_rel_path))
    
    # Read data from all job directories
    data_list <- list()
    for (job_idx in 1:n_jobs) {
      job_file <- file.path(dir_out, sprintf("job_%03d", job_idx), file_rel_path)
      if (file.exists(job_file)) {
        tryCatch({
          data_list[[job_idx]] <- data.table::fread(job_file)
        }, error = function(e) {
          warning(sprintf("Failed to read %s: %s", job_file, e$message))
        })
      }
    }
    
    # Combine data if any was read successfully
    if (length(data_list) > 0) {
      combined_data <- data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
      
      # Create output directory if needed
      output_file <- file.path(dir_out, file_rel_path)
      output_dir <- dirname(output_file)
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      
      # Write combined data
      data.table::fwrite(combined_data, output_file)
      
      if (verbose) {
        message(sprintf("  Combined %d files -> %d rows", 
                        length(data_list), nrow(combined_data)))
      }
    }
  }
  
  # 3. Consolidate args_df.csv.gz files with job column
  if (verbose) message("Consolidating args_df.csv.gz files...")
  
  args_list <- list()
  for (job_idx in 1:n_jobs) {
    args_file <- file.path(dir_out, sprintf("job_%03d", job_idx), "args_df.csv.gz")
    if (file.exists(args_file)) {
      args_data <- data.table::fread(args_file)
      args_data$job <- job_idx
      args_list[[job_idx]] <- args_data
    }
  }
  
  if (length(args_list) > 0) {
    combined_args <- data.table::rbindlist(args_list, use.names = TRUE)
    # Write to main directory
    data.table::fwrite(combined_args, file.path(dir_out, "args_df_combined.csv.gz"))
    if (verbose) {
      message(sprintf("  Combined %d args_df files", length(args_list)))
    }
  }
  
  # 4. Archive job directories
  if (verbose) message("Archiving job directories...")
  
  for (job_idx in 1:n_jobs) {
    job_dir <- sprintf("job_%03d", job_idx)
    job_path <- file.path(dir_out, job_dir)
    archive_name <- file.path(dir_out, paste0(job_dir, ".tar.gz"))
    
    # Create tar.gz archive
    # Using system command for better compression
    tar_cmd <- sprintf("tar -czf '%s' -C '%s' '%s'", 
                       archive_name, dir_out, job_dir)
    
    tryCatch({
      system(tar_cmd, intern = FALSE)
      # Remove directory after successful archiving
      if (file.exists(archive_name)) {
        unlink(job_path, recursive = TRUE)
        if (verbose) message(sprintf("  Archived and removed: %s", job_dir))
      }
    }, error = function(e) {
      warning(sprintf("Failed to archive %s: %s", job_dir, e$message))
    })
  }
  
  if (verbose) {
    message("\nOutput collection completed!")
    message(sprintf("Consolidated files are in: %s", dir_out))
    message("Job directories have been archived as .tar.gz files")
  }
}

# Internal function for running a single job (original logic)
genepicoloc_wrapper_single_job <- function(dir_out,
                                           sumstats_1_form,
                                           args_df,
                                           mc_cores = 10,
                                           verbose = TRUE,
                                           debug_mode = FALSE,
                                           progress_interval = 5,
                                           is_subjob = FALSE,
                                           job_info = NULL) {
  
  message("Starting genepicoloc: ", nrow(args_df), " sumstats to be processed.")
  
  # Create output directory
  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }
  if (verbose) message("Output will be written to ", dir_out)
  
  # Save arguments for reproducibility
  args_df_file <- file.path(dir_out, "args_df.csv.gz")
  if (verbose) message("Saving args_df to ", args_df_file)
  data.table::fwrite(args_df, args_df_file)
  
  # Save job info if this is a subjob
  if (is_subjob && !is.null(job_info)) {
    job_info_file <- file.path(dir_out, "job_info.txt")
    writeLines(c(
      sprintf("Job %d of %d", job_info$job_idx, job_info$n_jobs),
      sprintf("Regions %d to %d", job_info$start_idx, job_info$end_idx)
    ), job_info_file)
  }
  
  # Setup common arguments for all iterations
  shared_args <- list(
    dir_out = dir_out,
    sumstats_1_form = sumstats_1_form,
    verbose = verbose
  )
  
  # Prepare arguments for mapply
  mapply_args <- list(
    FUN = process_sumstats_form,
    sumstats_2_study = args_df$sumstats_2_study,
    sumstats_2_file = args_df$sumstats_2_file,
    sumstats_2_function = args_df$sumstats_2_function,
    sumstats_2_type = args_df$sumstats_2_type,
    sumstats_2_sdY = args_df$sumstats_2_sdY,
    MoreArgs = shared_args,
    SIMPLIFY = FALSE
  )
  
  # Choose execution strategy
  if (debug_mode) {
    # Sequential execution for debugging
    parallel_func <- mapply
  } else {
    # Standard parallel execution
    parallel_func <- mcmapply
    mapply_args$mc.cores <- mc_cores
  }
  
  # Print study summary before starting
  if (verbose) {
    total_datasets <- nrow(args_df)
    study_counts <- table(args_df$sumstats_2_study)
    message("Study breakdown:")
    for (study in names(study_counts)) {
      message(sprintf("  - %s: %d datasets", study, study_counts[study]))
    }
    message("Starting analysis...")
    
    # Simple progress estimation
    start_time <- Sys.time()
    message(sprintf("Processing %d datasets across %d studies...", 
                    total_datasets, length(study_counts)))
    
    # Start progress monitoring by watching output files
    # Function to monitor output files
    monitor_output <- function() {
      last_count <- 0
      while (TRUE) {
        
        # Count completed output files (assuming they have a specific pattern)
        tryCatch({
          output_files <- list.files(dir_out, pattern = "\\.csv\\.gz$", 
                                     recursive = TRUE, full.names = FALSE)
          # Exclude args_df.csv.gz from count
          output_files <- output_files[!grepl("^args_df\\.csv\\.gz$", basename(output_files))]
          current_count <- length(output_files)
          
          if (current_count > last_count && current_count > 0) {
            elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
            rate <- current_count / elapsed
            eta <- if (rate > 0) round((total_datasets - current_count) / rate, 1) else NA
            
            progress_msg <- sprintf("Progress: %d/%d files completed (%.1f%%) | %.1f min elapsed", 
                                    current_count, total_datasets, 
                                    (current_count/total_datasets)*100, elapsed)
            if (!is.na(eta) && eta > 0) {
              progress_msg <- paste0(progress_msg, sprintf(" | ETA: %.1f min", eta))
            }
            message(progress_msg)
            last_count <- current_count
          }
          
          # Stop monitoring if we've reached the expected number
          if (current_count >= total_datasets) break
        }, error = function(e) {
          # Continue monitoring if file listing fails
        })
        
        Sys.sleep(progress_interval)  # Check at specified interval
        
      }
    }
    
    # Start monitoring in background (Unix systems only)
    if (.Platform$OS.type == "unix") {
      monitor_pid <- parallel::mcparallel(monitor_output())
    }
  }
  
  # Execute analysis
  results <- do.call(parallel_func, mapply_args)
  
  if (verbose) {
    Sys.sleep(progress_interval)  # Check at specified interval
    end_time <- Sys.time()
    duration <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 1)
    message(sprintf("Analysis completed successfully in %.1f minutes!", duration))
  }
  
  # After the Sys.sleep and before the end
  if (exists("monitor_pid") && !is.null(monitor_pid) && .Platform$OS.type == "unix") {
    tryCatch({
      parallel::mccollect(monitor_pid, wait = FALSE)
    }, error = function(e) {
      # Ignore if already collected
    })
  }
  
  # Return invisibly
  invisible(results)
}

# Group 1: Format and Process Functions ----

#' Format primary summary statistics for colocalization analysis
#' 
#' @description
#' Retrieves and formats primary GWAS summary statistics for use in colocalization
#' analysis. This function serves as a wrapper that first retrieves raw data from
#' a specified source and then formats it according to the standardized format
#' required by the colocalization pipeline.
#' 
#' @param coloc_regions_PASS A data.frame containing genomic regions to analyze.
#'   Must contain columns CHR_var, BP_START_var, and BP_STOP_var defining the
#'   chromosome and position boundaries for each region.
#' @param sumstats_1_function Character string. Name of the function to use for
#'   retrieving summary statistics (e.g., "retrieve_sumstats_tabix", "tabix_UKB_PPP_EUR").
#' @param sumstats_1_file Character string. Path to the primary summary statistics file.
#'   Should be a tabix-indexed file or compatible with the specified retrieval function.
#' @param sumstats_1_type Character string. Type of trait for primary summary statistics.
#'   Must be either "quant" (quantitative) or "cc" (case-control).
#' @param sumstats_1_sdY Numeric. Standard deviation of the trait for quantitative
#'   phenotypes. Use NA for case-control studies or when not available.
#' 
#' @return A formatted sumstats object of class "sumstats" with standardized column
#'   names and attributes including:
#'   \itemize{
#'     \item Standardized columns: Name, rsID, CHR, POS, A1, A2, BETA, SE, nlog10P, AF, N
#'     \item Attributes: sumstats_type, sumstats_sdY, coloc_regions_PASS, and others
#'   }
#' 
#' @details
#' This function performs a two-step process:
#' \enumerate{
#'   \item Retrieves raw summary statistics using the specified function
#'   \item Formats the data including QC and standardization
#' }
#' 
#' The formatting step includes:
#' \itemize{
#'   \item Quality control (removing duplicates, infinite values, NAs)
#'   \item Column name standardization
#'   \item Attribute validation
#' }
#' 
#' @seealso 
#' \code{\link{retrieve_sumstats_raw}} for data retrieval
#' \code{\link{format_sumstats}} for data formatting
#' \code{\link{format_sumstats_2}} for formatting secondary datasets
#' 
#' @examples
#' \dontrun{
#' # Define regions to analyze
#' regions <- data.frame(
#'   CHR_var = c(1, 2),
#'   BP_START_var = c(1000000, 2000000),
#'   BP_STOP_var = c(2000000, 3000000)
#' )
#' 
#' # Format primary GWAS data
#' sumstats_1 <- format_sumstats_1(
#'   coloc_regions_PASS = regions,
#'   sumstats_1_function = "retrieve_sumstats_tabix",
#'   sumstats_1_file = "path/to/gwas.gz",
#'   sumstats_1_type = "quant",
#'   sumstats_1_sdY = 1.5
#' )
#' }
#' 
#' @export
format_sumstats_1 <- function(coloc_regions_PASS,
                              sumstats_1_function,
                              sumstats_1_file,
                              sumstats_1_type,
                              sumstats_1_sdY) {
  
  # Input validation
  if (!is.data.frame(coloc_regions_PASS)) {
    stop("coloc_regions_PASS must be a data.frame")
  }
  
  required_cols <- c("CHR_var", "BP_START_var", "BP_STOP_var")
  missing_cols <- setdiff(required_cols, colnames(coloc_regions_PASS))
  if (length(missing_cols) > 0) {
    stop("coloc_regions_PASS missing required columns: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  if (!sumstats_1_type %in% c("quant", "cc")) {
    stop("sumstats_1_type must be either 'quant' or 'cc'")
  }
  
  if (!file.exists(sumstats_1_file)) {
    stop("sumstats_1_file does not exist: ", sumstats_1_file)
  }
  
  # Retrieve raw summary statistics
  sumstats_1_raw <- retrieve_sumstats_raw(
    sumstats_function = sumstats_1_function,
    sumstats_file = sumstats_1_file,
    coloc_regions_PASS = coloc_regions_PASS
  )
  
  # Check retrieval status
  if (attr(sumstats_1_raw, "tabix") == "tabix_failed") {
    stop("Failed to retrieve summary statistics from: ", sumstats_1_file)
  }
  
  # Format summary statistics
  sumstats_1_form <- format_sumstats(
    sumstats = sumstats_1_raw,
    sumstats_type = sumstats_1_type,
    sumstats_sdY = sumstats_1_sdY
  )
  
  # Add coloc_regions_PASS as an attribute for downstream use
  attr(sumstats_1_form, "coloc_regions_PASS") <- coloc_regions_PASS
  
  return(sumstats_1_form)
}


#' Format secondary summary statistics for colocalization analysis
#' 
#' @description
#' Retrieves and formats secondary GWAS summary statistics, handling both single
#' and multi-phenotype datasets (e.g., eQTL data with multiple genes).
#' 
#' @param coloc_regions_PASS A data.frame containing genomic regions to analyze.
#'   Must contain columns CHR_var, BP_START_var, and BP_STOP_var.
#' @param sumstats_2_function Character string. Name of the function to use for
#'   retrieving summary statistics (e.g., "retrieve_sumstats_tabix").
#' @param sumstats_2_file Character string. Path to the secondary summary statistics file.
#' @param sumstats_2_type Character string. Type of trait for secondary summary statistics.
#'   Must be either "quant" (quantitative) or "cc" (case-control).
#' @param sumstats_2_sdY Numeric. Standard deviation of the trait for quantitative
#'   phenotypes. Use NA for case-control studies.
#' @param sumstats_2_study Character string. Study identifier for tracking purposes.
#' 
#' @return An object of class "sumstats_2_form" which is either:
#'   \itemize{
#'     \item A single sumstats object (if sumstats_pheno = "single")
#'     \item A list of sumstats objects (if sumstats_pheno = "multiple")
#'   }
#'   With attributes: sumstats_pheno, sumstats_file
#' 
#' @details
#' This function handles the complexity of datasets that may contain multiple
#' phenotypes (e.g., gene expression for multiple genes) by splitting them
#' appropriately and formatting each separately.
#' 
#' @seealso 
#' \code{\link{format_sumstats_1}} for formatting primary datasets
#' \code{\link{process_sumstats_form}} which calls this function
#' 
#' @export
format_sumstats_2 <- function(coloc_regions_PASS,
                              sumstats_2_function,
                              sumstats_2_file,
                              sumstats_2_type,
                              sumstats_2_sdY,
                              sumstats_2_study) {
  
  # Input validation
  if (!sumstats_2_type %in% c("quant", "cc")) {
    stop("sumstats_2_type must be either 'quant' or 'cc'")
  }
  
  # Retrieve raw summary statistics
  sumstats_2_raw <- retrieve_sumstats_raw(
    sumstats_function = sumstats_2_function,
    sumstats_file = sumstats_2_file,
    coloc_regions_PASS = coloc_regions_PASS
  )
  
  # Validate required attributes
  if (is.null(attr(sumstats_2_raw, "sumstats_pheno"))) {
    stop("sumstats_pheno attribute in sumstats_2_raw not found")
  }
  
  if (!attr(sumstats_2_raw, "sumstats_pheno") %in% c("single", "multiple")) {
    stop("sumstats_pheno should be either 'single' or 'multiple'")
  }
  
  # Process based on phenotype type
  if (attr(sumstats_2_raw, "sumstats_pheno") == "multiple") {
    # Multi-phenotype processing (e.g., eQTL data)
    if (!"Phenotype" %in% colnames(sumstats_2_raw)) {
      stop("'Phenotype' column in sumstats_2_raw not found, ",
           "required for 'multiple' sumstats_pheno attribute")
    }
    
    # Get unique phenotypes
    Phenotypes <- unique(sumstats_2_raw$Phenotype)
    
    # Handle empty data case
    if (length(Phenotypes) == 0) {
      # Convert to single phenotype with empty data
      attr(sumstats_2_raw, "sumstats_pheno") <- "single"
      sumstats_2_raw$Phenotype <- NULL
      sumstats_2_form <- format_sumstats(
        sumstats = sumstats_2_raw,
        sumstats_type = sumstats_2_type,
        sumstats_sdY = sumstats_2_sdY
      )
      # Wrap in list for consistency
      sumstats_2_form <- list(sumstats_2_form)
    } else {
      # Process each phenotype separately
      sumstats_2_form <- lapply(Phenotypes, function(pheno) {
        # Subset to specific phenotype
        sumstats_2_pheno <- subset(sumstats_2_raw, Phenotype == pheno)
        
        # Update file attribute to include phenotype
        attr(sumstats_2_pheno, "sumstats_file") <- 
          paste0(attr(sumstats_2_raw, "sumstats_file"), "_", pheno)
        
        # Set as single phenotype for formatting
        attr(sumstats_2_pheno, "sumstats_pheno") <- "single"
        
        # Remove Phenotype column after subsetting
        sumstats_2_pheno$Phenotype <- NULL
        
        # Format the subset
        sumstats_2_formatted <- format_sumstats(
          sumstats = sumstats_2_pheno,
          sumstats_type = sumstats_2_type,
          sumstats_sdY = sumstats_2_sdY
        )
        
        # Get the unique chromosomes in the formatted sumstats
        unique_variant_chrs <- unique(sumstats_2_formatted[["CHR"]])
        
        # Get the chromosome and position ranges for variants
        variant_chr_ranges <- data.frame(
          CHR_var = unique_variant_chrs,
          BP_START_var = sapply(unique_variant_chrs, function(chr) min(sumstats_2_formatted[CHR == chr, POS])),
          BP_STOP_var = sapply(unique_variant_chrs, function(chr) max(sumstats_2_formatted[CHR == chr, POS]))
        )
        
        # Subset the coloc_regions_PASS attribute
        coloc_regions_subset <- as.data.frame(attr(sumstats_2_formatted, "coloc_regions_PASS"))
        coloc_regions_subset <- coloc_regions_subset[
          coloc_regions_subset$CHR_var %in% unique_variant_chrs & 
            (
              (coloc_regions_subset$BP_START_var <= variant_chr_ranges$BP_STOP_var[match(coloc_regions_subset$CHR_var, variant_chr_ranges$CHR_var)]) & 
                (coloc_regions_subset$BP_STOP_var >= variant_chr_ranges$BP_START_var[match(coloc_regions_subset$CHR_var, variant_chr_ranges$CHR_var)])
            ),
        ]        
        
        # Update the attribute
        attr(sumstats_2_formatted, "coloc_regions_PASS") <- coloc_regions_subset
        
        return(sumstats_2_formatted)
      })
      
      names(sumstats_2_form) <- Phenotypes
    }
    
    # Set attributes for multiple phenotype object
    attr(sumstats_2_form, "sumstats_pheno") <- "multiple"
    attr(sumstats_2_form, "sumstats_file") <- sumstats_2_file
    
  } else {
    # Single phenotype processing
    sumstats_2_form <- format_sumstats(
      sumstats = sumstats_2_raw,
      sumstats_type = sumstats_2_type,
      sumstats_sdY = sumstats_2_sdY
    )
    attr(sumstats_2_form, "sumstats_pheno") <- "single"
  }
  
  # Set class
  class(sumstats_2_form) <- c("sumstats_2_form", class(sumstats_2_form))
  
  return(sumstats_2_form)
}


#' Process secondary summary statistics and run colocalization
#'
#' @description
#' Main processing function that formats secondary summary statistics and runs
#' colocalization analysis against a primary dataset across all specified regions.
#' This is typically called by \code{genepicoloc_wrapper} for each secondary dataset.
#'
#' @param dir_out Character string. Base output directory for results.
#' @param sumstats_1_form A formatted sumstats object containing the primary dataset.
#'   Must have a "coloc_regions_PASS" attribute.
#' @param sumstats_2_function Character string. Function name for retrieving secondary data.
#' @param sumstats_2_file Character string. Path to secondary summary statistics file.
#' @param sumstats_2_type Character string. Type of secondary trait ('quant' or 'cc').
#' @param sumstats_2_sdY Numeric. Standard deviation for quantitative traits (NA for cc).
#' @param sumstats_2_study Character string. Study identifier for output organization.
#' @param write_excel Logical. Whether to create Excel files with filtered results 
#'   (default: FALSE).
#' @param PP_H4_threshold Numeric. Minimum PP.H4.abf value for filtered results 
#'   (default: 0.5).
#' @param verbose Logical. Print detailed progress messages (default: TRUE).
#'
#' @return Invisibly returns the path to the output file.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Formats the secondary summary statistics
#'   \item Runs colocalization for each region
#'   \item Handles both single and multi-phenotype datasets
#'   \item Saves results to compressed CSV and optionally Excel files
#' }
#'
#' @importFrom data.table rbindlist
#' @export
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
  
  # Extract regions from primary dataset
  if (is.null(attr(sumstats_1_form, "coloc_regions_PASS"))) {
    stop("'coloc_regions_PASS' attribute not found in sumstats_1_form object")
  }
  coloc_regions_PASS <- attr(sumstats_1_form, "coloc_regions_PASS")
  
  # Format secondary summary statistics
  sumstats_2_form <- format_sumstats_2(
    coloc_regions_PASS = coloc_regions_PASS,
    sumstats_2_function = sumstats_2_function,
    sumstats_2_file = sumstats_2_file,
    sumstats_2_type = sumstats_2_type,
    sumstats_2_sdY = sumstats_2_sdY,
    sumstats_2_study = sumstats_2_study
  )
  
  # Validate formatted object
  if (!inherits(sumstats_2_form, "sumstats_2_form")) {
    stop("sumstats_2_form should have 'sumstats_2_form' class")
  }

  # Set up output path
  file_out <- set_output_path(
    sumstats_2_file = attr(sumstats_2_form, "sumstats_file"),
    dir_out = dir_out,
    sumstats_2_study = sumstats_2_study
  )

  # Run colocalization based on phenotype type
  if (attr(sumstats_2_form, "sumstats_pheno") == "single") {
    # Single phenotype analysis
    coloc_results <- run_single_phenotype(
      sumstats_1 = sumstats_1_form,
      sumstats_2 = sumstats_2_form,
      coloc_regions_PASS = coloc_regions_PASS
    )
    # Combine results from all regions
    coloc_results <- data.table::rbindlist(coloc_results, fill = TRUE)
    
  } else {
    # Multiple phenotype analysis
    coloc_results <- lapply(sumstats_2_form, function(sumstats_2) {
      region_results <- run_single_phenotype(
        sumstats_1 = sumstats_1_form,
        sumstats_2 = sumstats_2,
        coloc_regions_PASS = attr(sumstats_2, "coloc_regions_PASS")
      )
      # Combine regions for this phenotype
      data.table::rbindlist(region_results, fill = TRUE)
    })
    # Combine all phenotypes
    coloc_results <- data.table::rbindlist(coloc_results, fill = TRUE)
  }
  
  # Save results
  save_coloc_results(
    coloc_results = coloc_results, 
    file_out = file_out,
    write_excel = write_excel,
    PP_H4_threshold = PP_H4_threshold,
    verbose = verbose
  )

  return(invisible(file_out))
}


#' Set up output paths for colocalization results
#'
#' @description
#' Creates a hierarchical directory structure for organizing colocalization
#' results by study and generates the output file path.
#'
#' @param dir_out Character string. Base output directory.
#' @param sumstats_2_file Character string. Path to secondary summary statistics file.
#' @param sumstats_2_study Character string. Study identifier for directory organization.
#'
#' @return Character string. Full path for the output file.
#'
#' @details
#' Creates directory structure: dir_out/sumstats_2_study/basename(sumstats_2_file)
#'
#' @export
set_output_path <- function(dir_out,
                            sumstats_2_file,
                            sumstats_2_study) {
  
  # Create study-specific subdirectory
  study_dir <- file.path(dir_out, sumstats_2_study)
  if (!dir.exists(study_dir)) {
    dir.create(study_dir, recursive = TRUE)
  }
  
  # Generate output filename based on input filename
  file_base <- basename(sumstats_2_file)
  file_out <- file.path(study_dir, file_base)
  
  return(file_out)
}


# Group 2: Core colocalization functions ----
#' Run colocalization analysis for a single phenotype across all regions
#' 
#' @description
#' Executes colocalization analysis for a single phenotype (or single dataset)
#' across all specified genomic regions. This function applies \code{run_region}
#' to each region in parallel using mapply.
#' 
#' @param sumstats_1 A sumstats object containing the primary dataset.
#' @param sumstats_2 A sumstats object containing the secondary dataset.
#' @param coloc_regions_PASS A data.frame with columns CHR_var, BP_START_var, 
#'   and BP_STOP_var defining the genomic regions to analyze.
#' 
#' @return A list where each element contains the colocalization results for
#'   one region. Each element is a data.frame with colocalization statistics.
#' 
#' @details
#' This function serves as a wrapper that applies region-based analysis across
#' all regions specified in coloc_regions_PASS. Results are returned as a list
#' that can be combined using \code{rbindlist}.
#' 
#' @seealso 
#' \code{\link{run_region}} for single region analysis
#' \code{\link{process_sumstats_form}} which calls this function
#' 
#' @export
run_single_phenotype <- function(sumstats_1, sumstats_2, coloc_regions_PASS) {
  
  # Validate inputs
  if (!inherits(sumstats_1, "sumstats")) {
    stop("sumstats_1 must be a sumstats object")
  }
  
  if (!inherits(sumstats_2, "sumstats")) {
    stop("sumstats_2 must be a sumstats object")
  }
  
  required_cols <- c("CHR_var", "BP_START_var", "BP_STOP_var")
  if (!all(required_cols %in% colnames(coloc_regions_PASS))) {
    stop("coloc_regions_PASS must contain columns: ", 
         paste(required_cols, collapse = ", "))
  }
  
  # Prepare arguments for mapply
  MoreArgs <- list(
    sumstats_1 = sumstats_1,
    sumstats_2 = sumstats_2
  )
  
  # Run analysis for each region
  coloc_results <- with(coloc_regions_PASS, 
                        mapply(
                          FUN = run_region,
                          CHR_var = CHR_var,
                          BP_START_var = BP_START_var,
                          BP_STOP_var = BP_STOP_var,
                          MoreArgs = MoreArgs,
                          SIMPLIFY = FALSE
                        )
  )
  
  return(coloc_results)
}


#' Run colocalization analysis for a specific genomic region
#' 
#' @description
#' Performs colocalization analysis between two traits within a specific genomic
#' region. This is the core function that executes the coloc.abf algorithm after
#' preparing and validating the data.
#' 
#' @param sumstats_1 A sumstats object containing the primary dataset.
#' @param sumstats_2 A sumstats object containing the secondary dataset.
#' @param CHR_var Character or numeric. Chromosome identifier.
#' @param BP_START_var Numeric. Start position of the region (base pairs).
#' @param BP_STOP_var Numeric. End position of the region (base pairs).
#' @param min_nlog10P Numeric. Minimum -log10(p-value) threshold for considering
#'   a region to have significant signal (default: 5, equivalent to p < 1e-5).
#' 
#' @return A data.frame containing one row with colocalization results including:
#'   \itemize{
#'     \item Region coordinates and metadata
#'     \item Maximum significance values for both datasets
#'     \item Colocalization status (coloc_done, no_signif_sumstats_2, no_SNP_intersect)
#'     \item Posterior probabilities (PP.H0-H4.abf) if analysis was performed
#'     \item Directionality analysis results
#'     \item Top colocalized SNPs and their probabilities
#'   }
#' 
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates a result template with metadata
#'   \item Extracts region-specific data from both datasets
#'   \item Checks for sufficient signal in both datasets
#'   \item Verifies SNP overlap between datasets
#'   \item Analyzes effect directionality at the lead SNP
#'   \item Runs coloc.abf if all conditions are met
#'   \item Formats and returns the results
#' }
#' 
#' The analysis requires:
#' \itemize{
#'   \item Significant signal (> min_nlog10P) in the primary dataset
#'   \item At least one overlapping SNP between datasets
#'   \item Valid summary statistics in the region
#' }
#' 
#' @importFrom coloc coloc.abf
#' @export
run_region <- function(sumstats_1, sumstats_2,
                       CHR_var, BP_START_var, BP_STOP_var,
                       min_nlog10P = 5) {
  
  # Input validation
  if (!inherits(sumstats_1, "sumstats")) {
    stop("Expected a sumstats object for sumstats_1")
  }
  
  if (!inherits(sumstats_2, "sumstats")) {
    stop("Expected a sumstats object for sumstats_2")
  }
  
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
  
  # Add gene annotations to the result
  # Note: Gene map can be provided as parameter or will use packaged data
  tryCatch({
    result <- add_gene_annotations(result, gene_map = NULL)
  }, error = function(e) { e$message })
  
  # Early return if secondary data retrieval failed
  if (attr(sumstats_2, "tabix") %in% c("tabix_failed", "tabix_ok_no_data")) {
    result$coloc <- paste0("skipped_", attr(sumstats_2, "tabix"))
    return(result)
  }
  
  # Define region filter function
  region_filter <- function(data) {
    subset(data, CHR == CHR_var & POS > BP_START_var & POS < BP_STOP_var)
  }
  
  # Extract region-specific data
  sumstats_1_sub <- region_filter(sumstats_1)
  
  # Check primary dataset has data
  if (nrow(sumstats_1_sub) == 0) {
    stop("No data found in sumstats_1 for region ", 
         CHR_var, ":", BP_START_var, "-", BP_STOP_var)
  }
  
  sumstats_2_sub <- region_filter(sumstats_2)
  
  # Calculate maximum significance values
  sumstats_1_sub <- set_max_nlog10P(sumstats_1_sub)
  sumstats_2_sub <- set_max_nlog10P(sumstats_2_sub)
  
  sumstats_1_max_nlog10P <- attr(sumstats_1_sub, "max_nlog10P")
  sumstats_2_max_nlog10P <- attr(sumstats_2_sub, "max_nlog10P")
  
  # Store max values in result
  result$sumstats_1_max_nlog10P <- sumstats_1_max_nlog10P
  result$sumstats_2_max_nlog10P <- sumstats_2_max_nlog10P
  
  # Check primary dataset has significant signal
  if (is.na(sumstats_1_max_nlog10P) || sumstats_1_max_nlog10P < min_nlog10P) {
    stop("No significant SNPs in sumstats_1 (max -log10P = ", 
         round(sumstats_1_max_nlog10P, 2), " < ", min_nlog10P, ")")
  }
  
  # Check secondary dataset has significant signal
  if (is.na(sumstats_2_max_nlog10P) || sumstats_2_max_nlog10P < min_nlog10P) {
    result$coloc <- "no_signif_sumstats_2"
    return(result)
  }
  
  # Check for SNP overlap
  common_snps <- intersect(sumstats_1_sub$Name, sumstats_2_sub$Name)
  if (length(common_snps) == 0) {
    result$coloc <- "no_SNP_intersect"
    result$nsnps <- 0
    return(result)
  }
  
  result$nsnps <- length(common_snps)
  
  # Analyze effect directionality at lead SNP
  directionality_df <- get_directionality(
    sumstats_1 = sumstats_1_sub,
    sumstats_2 = sumstats_2_sub
  )
  
  # Add directionality results to output
  for (col in colnames(directionality_df)) {
    result[[col]] <- directionality_df[[col]]
  }
  
  # Run colocalization analysis
  coloc_output <- run_coloc(
    sumstats_1 = sumstats_1_sub,
    sumstats_2 = sumstats_2_sub
  )
  
  # Update status and format results
  result$coloc <- "coloc_done"
  result <- format_coloc_output(coloc_output, result)
  
  return(result)
}


#' Run colocalization analysis using coloc.abf
#'
#' @description
#' Executes the Bayesian colocalization analysis using the coloc.abf function
#' from the coloc package. This function prepares the input data and handles
#' output suppression if requested.
#'
#' @param sumstats_1 A sumstats object containing the primary dataset with
#'   required attributes sumstats_type and sumstats_sdY.
#' @param sumstats_2 A sumstats object containing the secondary dataset with
#'   required attributes sumstats_type and sumstats_sdY.
#' @param silent Logical. Whether to suppress coloc.abf output messages 
#'   (default: TRUE).
#'
#' @return A list containing coloc.abf results including:
#'   \itemize{
#'     \item summary: Named vector with nsnps and posterior probabilities
#'     \item results: Data frame with per-SNP results
#'     \item priors: Prior probabilities used
#'   }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Extracts required attributes from both datasets
#'   \item Prepares data in the format required by coloc.abf
#'   \item Runs the colocalization analysis
#'   \item Optionally suppresses console output for cleaner logs
#' }
#'
#' The coloc.abf function tests five hypotheses:
#' \itemize{
#'   \item H0: No association with either trait
#'   \item H1: Association with trait 1 only
#'   \item H2: Association with trait 2 only
#'   \item H3: Association with both traits, different variants
#'   \item H4: Association with both traits, same variant
#' }
#'
#' @importFrom coloc coloc.abf
#' @seealso 
#' \code{\link{prepare_coloc_dataset}} for data preparation
#' \code{\link[coloc]{coloc.abf}} for the underlying method
#' 
#' @export
run_coloc <- function(sumstats_1, sumstats_2, silent = TRUE) {
  
  # Input validation
  if (!inherits(sumstats_1, "sumstats")) {
    stop("Expected a sumstats object for sumstats_1")
  }
  
  if (!inherits(sumstats_2, "sumstats")) {
    stop("Expected a sumstats object for sumstats_2")
  }
  
  # Extract required attributes
  required_attrs <- c("sumstats_type", "sumstats_sdY")
  
  for (attr_name in required_attrs) {
    if (is.null(attr(sumstats_1, attr_name))) {
      stop("Missing required attribute '", attr_name, "' in sumstats_1")
    }
    if (is.null(attr(sumstats_2, attr_name))) {
      stop("Missing required attribute '", attr_name, "' in sumstats_2")
    }
  }
  
  sumstats_1_type <- attr(sumstats_1, "sumstats_type")
  sumstats_1_sdY <- attr(sumstats_1, "sumstats_sdY")
  sumstats_2_type <- attr(sumstats_2, "sumstats_type")
  sumstats_2_sdY <- attr(sumstats_2, "sumstats_sdY")
  
  # Prepare datasets for coloc.abf
  dataset_1 <- prepare_coloc_dataset(sumstats_1, sumstats_1_type, sumstats_1_sdY)
  dataset_2 <- prepare_coloc_dataset(sumstats_2, sumstats_2_type, sumstats_2_sdY)
  
  # Run colocalization
  if (silent) {
    # Suppress output by redirecting to temporary file
    tmp_file <- tempfile()
    sink(tmp_file, type = "output")
    on.exit({
      sink()  # Restore output
      unlink(tmp_file)  # Clean up
    }, add = TRUE)
    
    # Run with warnings suppressed
    coloc_output <- suppressWarnings(
      coloc::coloc.abf(dataset1 = dataset_1, dataset2 = dataset_2)
    )
  } else {
    # Run normally with output
    coloc_output <- coloc::coloc.abf(dataset1 = dataset_1, dataset2 = dataset_2)
  }
  
  return(coloc_output)
}


# Group 3: Helper Functions ----
#' Analyze effect directionality between datasets at lead SNP
#' 
#' @description
#' Identifies the most significant SNP in the primary dataset and analyzes
#' whether the effects in both datasets have the same direction. This helps
#' determine if the colocalized signal represents the same or opposite effects.
#' 
#' @param sumstats_1 A sumstats object containing the primary dataset.
#' @param sumstats_2 A sumstats object containing the secondary dataset.
#' 
#' @return A data.frame with one row containing:
#'   \itemize{
#'     \item sumstats_1_ind_Name: SNP identifier of lead variant
#'     \item sumstats_1_ind_A1/A2: Alleles in primary dataset
#'     \item sumstats_1_ind_nlog10P: Significance of lead variant
#'     \item sumstats_1_ind_BETA: Effect size in primary dataset
#'     \item sumstats_1_ind_BETA_sign: Sign of effect (+1 or -1)
#'     \item sumstats_2_ind_A1/A2: Alleles in secondary dataset
#'     \item sumstats_2_ind_nlog10P: Significance in secondary dataset
#'     \item sumstats_2_ind_BETA: Effect size in secondary dataset
#'     \item sumstats_2_ind_BETA_sign: Sign of effect (+1 or -1)
#'     \item directionality: 1 if same direction, -1 if opposite, NA if undetermined
#'   }
#' 
#' @details
#' The function:
#' \enumerate{
#'   \item Finds the intersection of SNPs between datasets
#'   \item Identifies the most significant SNP in the primary dataset
#'   \item Extracts information for this SNP from both datasets
#'   \item Determines effect directionality accounting for allele coding
#' }
#' 
#' Directionality is determined by comparing effect signs while accounting
#' for potential allele flips between datasets.
#' 
#' @seealso 
#' \code{\link{detect_directionality}} for the directionality logic
#' \code{\link{run_region}} which calls this function
#' 
#' @export
get_directionality <- function(sumstats_1, sumstats_2) {
  
  # Input validation
  if (!inherits(sumstats_1, "sumstats") || !inherits(sumstats_2, "sumstats")) {
    stop("Both inputs must be sumstats objects")
  }
  
  # Find common SNPs
  common_names <- intersect(sumstats_1$Name, sumstats_2$Name)
  
  if (length(common_names) == 0) {
    stop("No common SNPs found between datasets")
  }
  
  # Subset to common SNPs
  sumstats_1_common <- subset(sumstats_1, Name %in% common_names)
  sumstats_2_common <- subset(sumstats_2, Name %in% common_names)
  
  # Find most significant SNP in primary dataset
  max_idx <- which.max(sumstats_1_common$nlog10P)
  if (length(max_idx) == 0) {
    stop("Could not identify lead SNP")
  }
  
  # Extract lead SNP information
  lead_snp_1 <- sumstats_1_common[max_idx, ]
  
  # Find same SNP in secondary dataset
  lead_snp_2 <- subset(sumstats_2_common, Name == lead_snp_1$Name)
  
  if (nrow(lead_snp_2) == 0) {
    stop("Lead SNP not found in secondary dataset")
  }
  
  # Create output data frame
  directionality_df <- data.frame(
    sumstats_1_ind_Name = lead_snp_1$Name,
    sumstats_1_ind_A1 = lead_snp_1$A1,
    sumstats_1_ind_A2 = lead_snp_1$A2,
    sumstats_1_ind_nlog10P = lead_snp_1$nlog10P,
    sumstats_1_ind_BETA = lead_snp_1$BETA,
    sumstats_1_ind_BETA_sign = sign(lead_snp_1$BETA),
    sumstats_2_ind_A1 = lead_snp_2$A1,
    sumstats_2_ind_A2 = lead_snp_2$A2,
    sumstats_2_ind_nlog10P = lead_snp_2$nlog10P,
    sumstats_2_ind_BETA = lead_snp_2$BETA,
    sumstats_2_ind_BETA_sign = sign(lead_snp_2$BETA),
    stringsAsFactors = FALSE
  )
  
  # Determine directionality
  directionality_df <- detect_directionality(ind_df = directionality_df)
  
  return(directionality_df)
}


#' Detect effect directionality accounting for allele coding
#' 
#' @description
#' Determines whether effects in two datasets have the same or opposite
#' direction, accounting for potential differences in allele coding between
#' datasets.
#' 
#' @param ind_df A data.frame with one row containing allele and effect
#'   information for both datasets. Must include columns:
#'   sumstats_1_ind_A1, sumstats_2_ind_A1, sumstats_2_ind_A2,
#'   sumstats_1_ind_BETA_sign, sumstats_2_ind_BETA_sign
#' 
#' @return The input data.frame with an additional 'directionality' column:
#'   \itemize{
#'     \item 1: Effects are in the same direction
#'     \item -1: Effects are in opposite directions
#'     \item NA: Directionality cannot be determined (alleles don't match)
#'   }
#' 
#' @details
#' The function handles three cases:
#' \enumerate{
#'   \item A1 alleles match: Direct comparison of effect signs
#'   \item A1 in dataset 1 matches A2 in dataset 2: Flipped comparison
#'   \item No allele match: Returns NA
#' }
#' 
#' This is crucial for interpreting colocalization results, as opposite
#' directions might indicate different causal variants or pleiotropic effects.
#' 
#' @examples
#' \dontrun{
#' # Example with matching A1 alleles
#' df <- data.frame(
#'   sumstats_1_ind_A1 = "A",
#'   sumstats_2_ind_A1 = "A", 
#'   sumstats_2_ind_A2 = "G",
#'   sumstats_1_ind_BETA_sign = 1,
#'   sumstats_2_ind_BETA_sign = -1
#' )
#' detect_directionality(df)  # Returns df with directionality = -1
#' }
#' 
#' @export
detect_directionality <- function(ind_df) {
  
  # Input validation
  required_cols <- c("sumstats_1_ind_A1", "sumstats_2_ind_A1", 
                     "sumstats_2_ind_A2", "sumstats_1_ind_BETA_sign", 
                     "sumstats_2_ind_BETA_sign")
  
  missing_cols <- setdiff(required_cols, colnames(ind_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (nrow(ind_df) != 1) {
    stop("ind_df must have exactly one row")
  }
  
  # Extract values for clarity
  a1_dataset1 <- ind_df$sumstats_1_ind_A1
  a1_dataset2 <- ind_df$sumstats_2_ind_A1
  a2_dataset2 <- ind_df$sumstats_2_ind_A2
  sign1 <- ind_df$sumstats_1_ind_BETA_sign
  sign2 <- ind_df$sumstats_2_ind_BETA_sign
  
  # Determine directionality based on allele matching
  if (a1_dataset1 == a1_dataset2) {
    # Case 1: Same A1 allele - direct comparison
    directionality <- ifelse(sign1 == sign2, 1, -1)
    
  } else if (a1_dataset1 == a2_dataset2) {
    # Case 2: A1 in dataset1 matches A2 in dataset2 - flipped comparison
    # If signs are same, effects are actually opposite (due to allele flip)
    directionality <- ifelse(sign1 == sign2, -1, 1)
    
  } else {
    # Case 3: Alleles don't match - cannot determine directionality
    # This might happen with tri-allelic sites or data errors
    directionality <- NA
    warning("Alleles don't match between datasets for SNP: ", 
            ind_df$sumstats_1_ind_Name, 
            ". Dataset1 A1=", a1_dataset1, 
            ", Dataset2 alleles=", a1_dataset2, "/", a2_dataset2)
  }
  
  # Add to output
  ind_df$directionality <- directionality
  
  return(ind_df)
}


#' Prepare summary statistics for coloc.abf input
#'
#' @description
#' Formats a sumstats object into the specific list structure required by
#' the coloc.abf function. Handles both quantitative and case-control traits.
#'
#' @param sumstats A sumstats object containing standardized summary statistics.
#' @param sumstats_type Character string. Type of trait ('quant' or 'cc').
#' @param sumstats_sdY Numeric. Standard deviation of quantitative trait.
#'   Required for quant traits if MAF/N not available.
#'
#' @return A list formatted for coloc.abf with elements:
#'   \itemize{
#'     \item beta: Effect sizes
#'     \item varbeta: Variance of effect sizes (SE squared)
#'     \item snp: SNP identifiers
#'     \item type: Trait type ('quant' or 'cc')
#'     \item For quant traits: Either sdY or (MAF and N)
#'     \item For cc traits: Additional parameters may be added
#'   }
#'
#' @details
#' The coloc package requires specific input formats:
#' \itemize{
#'   \item For quantitative traits: Preferred to provide sdY (trait SD)
#'   \item If sdY unavailable: Falls back to MAF and N
#'   \item For case-control: Standard errors must account for case/control ratio
#' }
#'
#' @references
#' Giambartolomei et al. (2014) "Bayesian Test for Colocalisation between
#' Pairs of Genetic Association Studies Using Summary Statistics." PLOS Genetics.
#'
#' @seealso 
#' \code{\link[coloc]{coloc.abf}} for input requirements
#' \code{\link{run_coloc}} which calls this function
#'
#' @export
prepare_coloc_dataset <- function(sumstats, sumstats_type, sumstats_sdY) {
  
  # Input validation
  if (!inherits(sumstats, "sumstats")) {
    stop("Input must be a sumstats object")
  }
  
  if (!sumstats_type %in% c("quant", "cc")) {
    stop("sumstats_type must be 'quant' or 'cc'")
  }
  
  # Check required columns
  required_cols <- c("BETA", "SE", "Name")
  missing_cols <- setdiff(required_cols, colnames(sumstats))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create base dataset with essential fields
  dataset <- list(
    beta = sumstats$BETA,
    varbeta = sumstats$SE^2,
    snp = sumstats$Name,
    type = sumstats_type
  )
  
  # Add trait-specific parameters
  if (sumstats_type == "quant") {
    # For quantitative traits
    if (!is.na(sumstats_sdY) && is.numeric(sumstats_sdY) && sumstats_sdY > 0) {
      # Preferred: Use provided standard deviation
      dataset$sdY <- sumstats_sdY
    } else {
      # Fallback: Use MAF and N if available
      if (all(c("AF", "N") %in% colnames(sumstats))) {
        # Check for valid values
        if (any(is.na(sumstats$AF)) || any(is.na(sumstats$N))) {
          warning("Missing AF or N values; coloc.abf may not run properly")
        }
        dataset$MAF <- sumstats$AF
        dataset$N <- sumstats$N
      } else {
        warning("No sdY provided and AF/N columns missing for quantitative trait. ",
                "Coloc.abf may not run properly.")
      }
    }
  }
  # For case-control traits
  # Note: coloc.abf can handle cc traits with just beta/varbeta/type
  # but additional parameters like s (case proportion) can improve accuracy
  # TODO: Add case/control counts if available
  # if (all(c("N_cases", "N_controls") %in% colnames(sumstats))) {
  #   dataset$s <- sumstats$N_cases / (sumstats$N_cases + sumstats$N_controls)
  # }
  
  # Validate output
  if (any(is.infinite(dataset$varbeta))) {
    warning("Infinite variance values detected; these SNPs may cause issues")
  }
  
  if (any(dataset$varbeta == 0)) {
    warning("Zero variance values detected; these SNPs may cause issues")
  }
  
  return(dataset)
}
# Group 4: Output Functions ----
#' Format colocalization analysis output
#'
#' @description
#' Extracts and formats key results from the coloc.abf output, including
#' posterior probabilities and top colocalized SNPs, into a structured format
#' for reporting.
#'
#' @param coloc_output A list returned by coloc.abf containing:
#'   \itemize{
#'     \item summary: Named vector with posterior probabilities
#'     \item results: Data frame with per-SNP results
#'     \item priors: Prior probabilities used
#'   }
#' @param coloc_template A data.frame created by create_result_template
#'   containing region metadata and placeholder fields.
#' @param n_top_snps Integer. Number of top SNPs by PP.H4 to include 
#'   (default: 5).
#'
#' @return The coloc_template data.frame with added fields:
#'   \itemize{
#'     \item Top_coloc_SNP: Comma-separated list of top SNP identifiers
#'     \item Top_coloc_SNP.PP.H4: Formatted posterior probabilities
#'     \item priors: Prior probabilities used (as string)
#'     \item nsnps: Number of SNPs analyzed
#'     \item PP.H0-H4.abf: Posterior probabilities for each hypothesis
#'   }
#'
#' @details
#' The function formats probabilities using scientific notation for very
#' small values (< 0.01) and fixed decimal places otherwise. This ensures
#' readable output while preserving precision where needed.
#'
#' @seealso 
#' \code{\link{run_region}} which calls this function
#' \code{\link{create_result_template}} for the template structure
#'
#' @export
format_coloc_output <- function(coloc_output, coloc_template, n_top_snps = 5) {
  
  # Input validation
  if (!is.list(coloc_output) || !all(c("summary", "results", "priors") %in% names(coloc_output))) {
    stop("coloc_output must be a list with 'summary', 'results', and 'priors' elements")
  }
  
  if (!is.data.frame(coloc_template)) {
    stop("coloc_template must be a data.frame")
  }
  
  # Ensure n_top_snps is valid
  n_top_snps <- min(n_top_snps, nrow(coloc_output$results))
  
  if (n_top_snps > 0) {
    # Extract top SNPs by posterior probability
    results_df <- coloc_output$results
    
    # Sort by PP.H4 descending
    results_df <- results_df[order(results_df$SNP.PP.H4, decreasing = TRUE), ]
    
    # Get top SNPs
    top_snps <- head(results_df[, c("snp", "SNP.PP.H4")], n_top_snps)
    
    # Format SNP names
    snp_names <- paste(top_snps$snp, collapse = ", ")
    
    # Format probabilities with appropriate precision
    snp_probs <- top_snps$SNP.PP.H4
    formatted_probs <- sapply(snp_probs, function(p) {
      if (is.na(p)) {
        return("NA")
      } else if (p < 0.01) {
        # Use scientific notation for small probabilities
        return(format(p, scientific = TRUE, digits = 2))
      } else {
        # Use fixed notation with 2 decimal places
        return(sprintf("%.2f", p))
      }
    })
    formatted_probs_str <- paste(formatted_probs, collapse = ", ")
    
    # Add to template
    coloc_template$Top_coloc_SNP <- snp_names
    coloc_template$Top_coloc_SNP.PP.H4 <- formatted_probs_str
  } else {
    coloc_template$Top_coloc_SNP <- "None"
    coloc_template$Top_coloc_SNP.PP.H4 <- "NA"
  }
  
  # Format priors as string
  prior_values <- unlist(coloc_output$priors)
  prior_names <- names(prior_values)
  if (is.null(prior_names)) {
    prior_names <- paste0("p", seq_along(prior_values))
  }
  priors_formatted <- paste(paste0(prior_names, "=", format(prior_values, scientific = TRUE, digits = 2)), 
                            collapse = ", ")
  coloc_template$priors <- priors_formatted
  
  # Add number of SNPs
  coloc_template$nsnps <- coloc_output$summary["nsnps"]
  
  # Add posterior probabilities for each hypothesis
  pp_cols <- c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
  for (h in pp_cols) {
    if (h %in% names(coloc_output$summary)) {
      coloc_template[[h]] <- as.numeric(coloc_output$summary[h])
    } else {
      warning("Missing ", h, " in coloc output")
      coloc_template[[h]] <- NA
    }
  }
  
  return(coloc_template)
}


#' Create a template for colocalization results with gene annotation
#'
#' @description
#' Creates a standardized data.frame template for storing colocalization
#' results. This ensures consistent output format across all analyses.
#'
#' @param CHR_var Character or numeric. Chromosome identifier.
#' @param BP_START_var Numeric. Start position of the region.
#' @param BP_STOP_var Numeric. End position of the region.
#' @param sumstats_1_file Character. Path to primary summary statistics file.
#' @param sumstats_2_file Character. Path to secondary summary statistics file.
#' @param sumstats_1_tabix Character. Tabix status for primary data.
#' @param sumstats_2_tabix Character. Tabix status for secondary data.
#' @param sumstats_1_max_nlog10P Numeric. Maximum -log10(p) in primary data.
#' @param sumstats_2_max_nlog10P Numeric. Maximum -log10(p) in secondary data.
#' @param nsnps Numeric. Number of SNPs in the analysis.
#' @param sumstats_1_QC Character. QC status of primary data.
#' @param sumstats_2_QC Character. QC status of secondary data.
#' @param coloc Character. Colocalization analysis status.
#'
#' @return A data.frame with one row containing all metadata fields and
#'   placeholder NA values for results to be filled later.
#'
#' @details
#' The template includes fields for:
#' \itemize{
#'   \item Region coordinates and metadata
#'   \item Gene annotation (nearest genes, region center annotation)
#'   \item File paths and data quality indicators
#'   \item Colocalization results (PP.H0-H4)
#'   \item Directionality analysis results
#'   \item Top SNP information
#' }
#'
#' This standardized format facilitates downstream analysis and reporting.
#'
#' @seealso 
#' \code{\link{run_region}} which creates and populates this template
#' \code{\link{format_coloc_output}} which adds results to the template
#'
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
  
  # Create results data frame with all fields
  result_df <- data.frame(
    # Region information
    CHR_var = CHR_var,
    BP_START_var = BP_START_var,
    BP_STOP_var = BP_STOP_var,
    
    # Gene annotation fields (to be filled by gene annotation step)
    region_center_pos = NA_real_,
    region_center_gene = NA_character_,
    nearest_genes_3 = NA_character_,
    nearest_genes_5 = NA_character_,
    nearest_genes_10 = NA_character_,
    
    # Summary statistics file information
    sumstats_1_file = sumstats_1_file,
    sumstats_1_tabix = sumstats_1_tabix,
    sumstats_1_max_nlog10P = sumstats_1_max_nlog10P,
    sumstats_1_QC = sumstats_1_QC,
    
    sumstats_2_file = sumstats_2_file,
    sumstats_2_tabix = sumstats_2_tabix,
    sumstats_2_max_nlog10P = sumstats_2_max_nlog10P,
    sumstats_2_QC = sumstats_2_QC,
    
    # Analysis metadata
    nsnps = nsnps,
    coloc = coloc,
    
    # Colocalization results (to be filled later)
    PP.H0.abf = NA_real_, 
    PP.H1.abf = NA_real_, 
    PP.H2.abf = NA_real_,
    PP.H3.abf = NA_real_, 
    PP.H4.abf = NA_real_,
    Top_coloc_SNP = NA_character_, 
    Top_coloc_SNP.PP.H4 = NA_character_, 
    priors = NA_character_,
    
    # Directionality analysis results
    sumstats_1_ind_Name = NA_character_,
    sumstats_1_ind_A1 = NA_character_,
    sumstats_1_ind_A2 = NA_character_,
    sumstats_1_ind_nlog10P = NA_real_,
    sumstats_1_ind_BETA = NA_real_,
    sumstats_1_ind_BETA_sign = NA_integer_,
    sumstats_2_ind_A1 = NA_character_,
    sumstats_2_ind_A2 = NA_character_,
    sumstats_2_ind_nlog10P = NA_real_,
    sumstats_2_ind_BETA = NA_real_,
    sumstats_2_ind_BETA_sign = NA_integer_,
    directionality = NA_integer_,
    
    # Ensure proper data types
    stringsAsFactors = FALSE
  )
  
  # Ensure single row
  if (nrow(result_df) != 1) {
    stop("Template creation resulted in multiple rows")
  }
  
  return(result_df)
}

#' Add gene annotations to colocalization results
#'
#' @description
#' Annotates genomic regions with nearest genes using the gene annotation
#' system. Adds multiple gene annotation fields to provide context for
#' colocalization results.
#'
#' @param result A data.frame with region information (CHR_var, BP_START_var, BP_STOP_var).
#' @param gene_map Data.frame with gene annotations from setup_gene_annotation()$gene_map.
#'   If NULL, attempts to use packaged gene data.
#' @param verbose Logical. Print gene annotation details (default: FALSE).
#'
#' @return The input data.frame with added gene annotation columns:
#'   \itemize{
#'     \item region_center_pos: Center position of the region
#'     \item region_center_gene: Gene annotation at region center
#'     \item nearest_genes_3: 3 nearest genes (simple format)
#'     \item nearest_genes_5: 5 nearest genes (simple format)
#'     \item nearest_genes_10: 10 nearest genes (simple format)
#'   }
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Calculates the center position of the genomic region
#'   \item Annotates the center position with the nearest gene
#'   \item Finds multiple nearest genes for broader context
#'   \item Uses simple gene name format for easy interpretation
#' }
#'
#' @seealso 
#' \code{\link{annotate_position}} for the underlying annotation function
#' \code{\link{setup_gene_annotation}} for gene map preparation
#'
#' @examples
#' \dontrun{
#' # Prepare gene annotation
#' gene_setup <- setup_gene_annotation()
#' 
#' # Create a result template
#' result <- create_result_template(
#'   CHR_var = "1", BP_START_var = 1000000, BP_STOP_var = 2000000,
#'   sumstats_1_file = "file1", sumstats_2_file = "file2",
#'   sumstats_1_tabix = "ok", sumstats_2_tabix = "ok"
#' )
#' 
#' # Add gene annotations
#' result_annotated <- add_gene_annotations(result, gene_setup$gene_map)
#' }
#'
#' @export
add_gene_annotations <- function(result, gene_map = NULL, verbose = FALSE) {
  
  # Input validation
  if (!is.data.frame(result)) {
    stop("result must be a data.frame")
  }
  
  required_cols <- c("CHR_var", "BP_START_var", "BP_STOP_var")
  missing_cols <- setdiff(required_cols, colnames(result))
  if (length(missing_cols) > 0) {
    stop("result must contain columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Set up gene annotation if not provided
  if (is.null(gene_map)) {
    if (verbose) message("Loading packaged gene annotation data...")
    tryCatch({
      gene_setup <- setup_gene_annotation(verbose = FALSE)
      gene_map <- gene_setup$gene_map
    }, error = function(e) {
      warning("Failed to load gene annotation data: ", e$message)
      return(result)  # Return unchanged if gene annotation fails
    })
  }
  
  # Validate gene_map
  if (!is.data.frame(gene_map)) {
    warning("Invalid gene_map provided; skipping gene annotation")
    return(result)
  }
  
  # Process each row (though typically there's only one)
  for (i in seq_len(nrow(result))) {
    chr <- as.character(result$CHR_var[i])
    start_pos <- as.numeric(result$BP_START_var[i])
    stop_pos <- as.numeric(result$BP_STOP_var[i])
    
    # Calculate region center
    center_pos <- round((start_pos + stop_pos) / 2)
    result$region_center_pos[i] <- center_pos
    
    if (verbose) {
      message("Annotating region ", chr, ":", start_pos, "-", stop_pos, 
              " (center: ", center_pos, ")")
    }
    
    # Annotate center position
    tryCatch({
      # Get detailed annotation for center
      center_annotation <- annotate_position(
        chr = chr, 
        pos = center_pos, 
        gene_map = gene_map,
        n_nearest = 1,
        output_format = "default"
      )
      result$region_center_gene[i] <- center_annotation
      
      # Get multiple nearest genes in simple format
      nearest_3 <- annotate_position(
        chr = chr, 
        pos = center_pos, 
        gene_map = gene_map,
        n_nearest = 3,
        output_format = "simple"
      )
      result$nearest_genes_3[i] <- nearest_3
      
      nearest_5 <- annotate_position(
        chr = chr, 
        pos = center_pos, 
        gene_map = gene_map,
        n_nearest = 5,
        output_format = "simple"
      )
      result$nearest_genes_5[i] <- nearest_5
      
      nearest_10 <- annotate_position(
        chr = chr, 
        pos = center_pos, 
        gene_map = gene_map,
        n_nearest = 10,
        output_format = "simple"
      )
      result$nearest_genes_10[i] <- nearest_10
      
      if (verbose) {
        message("  Center gene: ", center_annotation)
        message("  5 nearest: ", nearest_5)
      }
      
    }, error = function(e) {
      warning("Gene annotation failed for region ", chr, ":", start_pos, "-", stop_pos, 
              ": ", e$message)
      # Leave as NA if annotation fails
    })
  }
  
  return(result)
}

#' Save colocalization results to files
#'
#' @description
#' Saves colocalization results to compressed CSV format and optionally creates
#' filtered Excel files for easy viewing of significant results.
#'
#' @param coloc_results A data.frame or data.table containing colocalization
#'   results from one or more regions/phenotypes.
#' @param file_out Character string. Base path for output files.
#' @param write_excel Logical. Whether to create Excel file with filtered
#'   results (default: FALSE).
#' @param PP_H4_threshold Numeric. Minimum PP.H4.abf value for filtering
#'   results in Excel output (default: 0.5).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#' @param sans_ext Logical. Whether to remove extension from file_out before
#'   adding new extensions (default: TRUE).
#'
#' @return Invisibly returns the base path used for output files.
#'
#' @details
#' Output files created:
#' \itemize{
#'   \item {basename}_ge.csv.gz: Complete results in compressed CSV format
#'   \item {basename}_filt.xlsx: Filtered results (if write_excel = TRUE)
#' }
#'
#' The Excel file only includes results where PP.H4.abf >= PP_H4_threshold,
#' making it easier to focus on likely colocalized signals.
#'
#' @importFrom data.table fwrite
#' @importFrom writexl write_xlsx
#' 
#' @seealso 
#' \code{\link{process_sumstats_form}} which calls this function
#'
#' @examples
#' \dontrun{
#' # Save results with Excel output for strong colocalizations
#' save_coloc_results(
#'   coloc_results = my_results,
#'   file_out = "output/my_analysis",
#'   write_excel = TRUE,
#'   PP_H4_threshold = 0.8
#' )
#' }
#'
#' @export
save_coloc_results <- function(coloc_results,
                               file_out,
                               write_excel = FALSE,
                               PP_H4_threshold = 0.5,
                               verbose = TRUE,
                               sans_ext = TRUE) {
  
  # Input validation
  if (is.null(coloc_results) || nrow(coloc_results) == 0) {
    warning("Empty colocalization results provided")
    return(invisible(NULL))
  }
  
  if (!is.data.frame(coloc_results)) {
    stop("coloc_results must be a data.frame or data.table")
  }
  
  if (!is.numeric(PP_H4_threshold) || PP_H4_threshold < 0 || PP_H4_threshold > 1) {
    stop("PP_H4_threshold must be a number between 0 and 1")
  }
  
  # Prepare base filename
  if (sans_ext) {
    file_out <- tools::file_path_sans_ext(file_out)
  }
  
  # Save complete results as compressed CSV
  csv_file <- paste0(file_out, "_ge.csv.gz")
  
  tryCatch({
    data.table::fwrite(coloc_results, csv_file, compress = "gzip")
  }, error = function(e) {
    stop("Failed to write CSV file: ", e$message)
  })
  
  # Optionally create Excel file with filtered results
  if (write_excel) {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      warning("Package 'writexl' is required for Excel output but not installed")
      return(invisible(file_out))
    }
    
    # Filter for significant colocalizations
    if ("PP.H4.abf" %in% colnames(coloc_results)) {
      coloc_df_filt <- coloc_results[!is.na(coloc_results$PP.H4.abf) & 
                                       coloc_results$PP.H4.abf >= PP_H4_threshold, ]
      
      if (nrow(coloc_df_filt) > 0) {
        # Sort by PP.H4.abf descending for easier review
        coloc_df_filt <- coloc_df_filt[order(coloc_df_filt$PP.H4.abf, 
                                             decreasing = TRUE), ]
        
        xlsx_file <- paste0(file_out, "_filt.xlsx")
        
        tryCatch({
          writexl::write_xlsx(coloc_df_filt, xlsx_file)
        }, error = function(e) {
          warning("Failed to write Excel file: ", e$message)
        })
      } else {
      }
    } else {
      warning("PP.H4.abf column not found; cannot create filtered Excel file")
    }
  }
  
  return(invisible(file_out))
}