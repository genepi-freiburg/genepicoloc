# Group 0: Pipeline orchestration ----
#' Three-level pipeline hierarchy:
#' 1. genepicoloc_wrapper() - Splits regions into jobs, manages parallelization
#' 2. genepicoloc_job() - Processes one job (region subset, all datasets) 
#' 3. genepicoloc_run() - Analyzes one secondary dataset against primary

#' Wrapper for parallel colocalization analysis
#' 
#' @description
#' Orchestrates parallel colocalization analysis across multiple secondary datasets
#' against a primary dataset. Automatically splits large region sets into manageable
#' jobs and handles parallelization.
#' 
#' @param dir_out Character. Output directory path for results.
#' @param sumstats_1_args Named list containing primary dataset parameters:
#'   \itemize{
#'     \item coloc_regions_PASS: Data frame with CHR_var, BP_START_var, BP_STOP_var
#'     \item sumstats_1_function: Function name for data retrieval
#'     \item sumstats_1_file: Path to primary summary statistics
#'     \item sumstats_1_type: 'quant' or 'cc'
#'     \item sumstats_1_sdY: Standard deviation (NA for cc)
#'   }
#' @param args_df Data frame with secondary dataset specifications:
#'   \itemize{
#'     \item sumstats_2_study: Study identifier
#'     \item sumstats_2_file: Path to summary statistics
#'     \item sumstats_2_function: Retrieval function name
#'     \item sumstats_2_type: 'quant' or 'cc'
#'     \item sumstats_2_sdY: Standard deviation (NA for cc)
#'   }
#' @param mc_cores Integer. Cores for parallel processing (default: 10).
#' @param verbose Logical. Print progress messages (default: TRUE).
#' @param debug_mode Logical. Run sequentially for debugging (default: FALSE).
#' @param max_regions_per_job Integer. Maximum regions per job (default: 10).
#'   For eQTL studies with multiple phenotypes, recommended value is not more
#'   than 10. For other studies with a single phenotype it can be increased to 100.
#' @param save_sumstats Logical. Whether to save filtered summary statistics 
#'   (default: FALSE). Set to TRUE only if you need the filtered summary statistics
#'   for downstream analyses like fine-mapping. WARNING: This substantially increases
#'   disk usage (~20-30GB per 1,00 regions and 10,000 datasets with typical settings).
#'   When FALSE, only colocalization results are saved (~2-3GB).
#' @param p_filt Numeric. Maximum p-value threshold for variants to include in 
#'   saved summary statistics (default: 1, includes all variants). Only applies 
#'   when save_sumstats = TRUE. Set to a lower value (e.g., 0.05) to reduce 
#'   storage by excluding variants with high p-values. Note: For fine-mapping 
#'   and conditional analyses, keep at 1 to preserve complete LD structure. 
#' @param p_min_save Numeric. Minimum p-value threshold for saving summary 
#'   statistics sumstats_2 (default: 5e-8). Only regions with at least one variant 
#'   with p-value below this threshold will be saved when save_sumstats is TRUE
#' @param batch_size Integer or NULL. Number of secondary datasets to process per 
#'   subjob (default: NULL, which sets it to mc_cores * 10). Controls Level 2 
#'   parallelization batching.
#'
#' @return NULL (invisibly). This function is called for its side effects of 
#'   writing analysis results to disk in the specified output directory.
#'
#' @details
#' The function implements a two-level parallelization strategy:
#' \enumerate{
#'   \item Level 1: Splits regions into jobs (max_regions_per_job each)
#'   \item Level 2: Processes secondary datasets in parallel within each job
#' }
#' 
#' Output structure:
#' \itemize{
#'   \item job_XXXX/: Subdirectories for each region job
#'   \item sumstats_1.RDS: Primary dataset for job's regions
#'   \item sumstats.tar: Archive of secondary datasets
#'   \item coloc.tar: Archive of colocalization results
#'   \item Metadata files in parent directory
#' }
#' 
#' @examples
#' \dontrun{
#' sumstats_1_args <- list(
#'   coloc_regions_PASS = coloc_regions_PASS,
#'   sumstats_1_function = "retrieve_sumstats_tabix",
#'   sumstats_1_file = "sumstats_1.gz",
#'   sumstats_1_type = "quant",
#'   sumstats_1_sdY = 1
#' )
#' 
#' args_df <- data.frame(
#'   sumstats_2_study = c("Study1", "Study2"),
#'   sumstats_2_file = c("study1_sumstats_2.gz", "study2_sumstats_2.gz"),
#'   sumstats_2_function = "retrieve_sumstats_tabix",
#'   sumstats_2_type = "quant",
#'   sumstats_2_sdY = NA
#' )
#' 
#' genepicoloc_wrapper(
#'   dir_out = "results",
#'   sumstats_1_args = sumstats_1_args,
#'   args_df = args_df,
#'   mc_cores = 4
#' )
#' }
#' 
#' @export
genepicoloc_wrapper <- function(dir_out,
                                sumstats_1_args,
                                args_df,
                                mc_cores = 2,
                                verbose = TRUE,
                                debug_mode = FALSE,
                                max_regions_per_job = 10,
                                save_sumstats = FALSE,
                                p_filt=1,
                                p_min_save = 5e-8,
                                batch_size = NULL) {
  
  # Warning for save_sumstats = TRUE
  if (save_sumstats) {
    warning_msg <- paste0(
      "\n",
      "================================================================================\n",
      "WARNING: save_sumstats = TRUE will generate substantial disk usage!\n",
      "--------------------------------------------------------------------------------\n",
      "To reduce storage, consider:\n",
      "  1. Set save_sumstats = FALSE (saves ~90% space)\n",
      "  2. Use stricter p_min_save threshold (e.g., 1e-10)\n",
      "  3. Process fewer regions or traits at once\n",
      "================================================================================\n"
    )
    
    if (interactive()) {
      # In interactive mode, ask for confirmation
      cat(warning_msg)
      response <- readline(prompt = "Continue with save_sumstats = TRUE? (y/n): ")
      if (!tolower(response) %in% c("y", "yes")) {
        stop("Analysis cancelled by user. Set save_sumstats = FALSE to save disk space.")
      }
    } else {
      # In non-interactive mode, just show warning
      warning(warning_msg, immediate. = TRUE)
    }
  } else {
    message("save_sumstats is FALSE (default), ",
            "saving only colocalization results. ",
            "This option is sufficient for most use cases.")
  }
  
  # Input validation for sumstats_1_args
  required_args <- c("coloc_regions_PASS", "sumstats_1_function", 
                     "sumstats_1_file", "sumstats_1_type", "sumstats_1_sdY")
  missing_args <- setdiff(required_args, names(sumstats_1_args))
  if (length(missing_args) > 0) {
    stop("sumstats_1_args missing required elements: ", 
         paste(missing_args, collapse = ", "))
  }
  
  # Extract coloc_regions_PASS from the args
  coloc_regions_PASS <- sumstats_1_args$coloc_regions_PASS
  
  # Jobs metadata
  study_counts <- table(args_df$sumstats_2_study)
  if (verbose) {
    message("Total number of regions: ", nrow(coloc_regions_PASS))
    message("Total number of argument: ", nrow(args_df))
    message("Study breakdown:")
    for (study in names(study_counts)) {
      message(sprintf("  - %s: %d datasets", study, study_counts[study]))
    }
  }
  
  ### STEP 1: Process region-based jobs ###
  # Determine parallelization for level 1 - regions
  n_regions <- nrow(coloc_regions_PASS)
  n_jobs_1 <- ceiling(n_regions / max_regions_per_job)
  if (verbose) {
    message(sprintf("Splitting %d regions into %d jobs (max %d regions per job).",
                    n_regions, n_jobs_1, max_regions_per_job))
  }
  
  job_times <- sapply(1:n_jobs_1, function(job_idx) {
    
    job_start_time <- Sys.time()
    
    # Calculate region indices for this job
    start_idx <- (job_idx - 1) * max_regions_per_job + 1
    end_idx <- min(job_idx * max_regions_per_job, n_regions)
    
    if (verbose) {
      message(sprintf("\n=== Running job %d/%d (regions %d-%d) ===", 
                      job_idx, n_jobs_1, start_idx, end_idx))
    }
    
    # Create job-specific output directory
    dir_out <- file.path(dir_out, "jobs", sprintf("job_%04d", job_idx))
    if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
    
    # Format sumstats_1 for this job's regions only
    sumstats_1_form <- format_sumstats_1(
      coloc_regions_PASS = coloc_regions_PASS[start_idx:end_idx,],
      sumstats_1_function = sumstats_1_args$sumstats_1_function,
      sumstats_1_file = sumstats_1_args$sumstats_1_file,
      sumstats_1_type = sumstats_1_args$sumstats_1_type,
      sumstats_1_sdY = sumstats_1_args$sumstats_1_sdY
    )
    
    if (verbose) {
      mb_size <- format(object.size(sumstats_1_form), units = "MB")
      if (mb_size == "0 Mb") mb_size <- format(object.size(sumstats_1_form), units = "KB")
      message("Output will be written to ", dir_out)
      message("sumstats_1: Queried ", format(nrow(sumstats_1_form), big.mark = ","),
              " rows, object size = ", mb_size)
    }
    
    ### Saving sumstats_1 and creating tar archive for sumstats_2
    if (save_sumstats) {
      saveRDS(sumstats_1_form, file.path(dir_out, "sumstats_1.RDS"))
      tar_sumstats_file <- file.path(dir_out, "sumstats.tar")
      exit_code_sumstats <- system(sprintf("tar -cf '%s' -T /dev/null", tar_sumstats_file), ignore.stderr = T)
      if (exit_code_sumstats != 0) {
        stop(sprintf("Failed to create empty sumstats tar archive: %s (exit code: %d)", 
                     tar_sumstats_file, exit_code_sumstats))
      }
    } else {
      tar_sumstats_file <- NULL
    }
    
    # Create empty tar archive
    tar_coloc_file <- file.path(dir_out, "coloc.tar") # save all in one, not by study
    exit_code_coloc <- system(sprintf("tar -cf '%s' -T /dev/null", tar_coloc_file))
    # Check if creation succeeded
    if (exit_code_coloc != 0) {
      stop(sprintf("Failed to create empty coloc tar archive: %s (exit code: %d)", 
                   tar_coloc_file, exit_code_coloc))
    }
    
    ### STEP 2: Process arguments-based subjobs ###
    n_args <- nrow(args_df)
    # Set batch size if not specified
    if (is.null(batch_size)) batch_size <- mc_cores * 10
    
    n_jobs_2 <- ceiling(n_args / batch_size) # max jobs based on mc_cores
    if (verbose) {
      message(sprintf("Splitting %d arguments into %d subjobs (batch size: %d, using %d cores).",
                      n_args, n_jobs_2, batch_size, mc_cores))
    }
    
    sapply(1:n_jobs_2, function(job_idx) {
      # Calculate indices for this subjob
      start_idx <- (job_idx - 1) * batch_size + 1
      end_idx <- min(job_idx * batch_size, n_args)
      n_expected <- end_idx - start_idx + 1
      
      # Run the analysis for this subset
      result <- genepicoloc_job(
        sumstats_1_form = sumstats_1_form,
        args_df = args_df[start_idx:end_idx,],
        tar_sumstats_file = tar_sumstats_file,
        tar_coloc_file = tar_coloc_file,
        mc_cores = mc_cores,
        verbose = verbose,
        debug_mode = debug_mode,
        save_sumstats = save_sumstats,
        p_filt = p_filt,
        p_min_save = p_min_save
      )
      
      # Validate output
      if (result$n_coloc != n_expected) {
        warning(sprintf("Subjob %d: File count mismatch. Expected %d files, got %d coloc files", 
                        job_idx, n_expected, result$n_coloc))
      }
      
      # Progress reporting
      if (verbose) {
        completed_args <- min(job_idx * batch_size, n_args)
        pct_complete <- (completed_args / n_args) * 100
        elapsed <- round(as.numeric(difftime(Sys.time(), job_start_time, units = "mins")), 1)
        
        # Estimate remaining time
        if (completed_args > 0) {
          rate <- completed_args / elapsed
          remaining <- n_args - completed_args
          eta <- if (rate > 0) round(remaining / rate, 1) else NA
          
          progress_msg <- sprintf("Running subjob %d/%d (args %d-%d): %d/%d datasets completed (%.1f%%) | %.1f min elapsed", 
                                  job_idx, n_jobs_2, start_idx, end_idx,
                                  completed_args, n_args, pct_complete, elapsed)
          
          if (!is.na(eta) && eta > 0) {
            progress_msg <- paste0(progress_msg, sprintf(" | ETA: %.1f min", eta))
          }
          
          message(progress_msg)
        }
        
      }
      
    })
    
    # Final summary
    total_time <- round(as.numeric(difftime(Sys.time(), job_start_time, units = "mins")), 1)
    
    # Return the job directory path
    return(invisible(total_time))
  })
  
  # Save arguments and job splitting information
  save_job_metadata(
    dir_out = dir_out,
    args_df = args_df,
    n_regions = n_regions,
    max_regions_per_job = max_regions_per_job,
    n_jobs_1 = n_jobs_1,
    mc_cores = mc_cores,
    study_counts = study_counts,
    job_times = unlist(job_times),
    verbose = verbose
  )
  
  if (verbose) {
    message("\n=== Consolidating results by study ===")
  }
  consolidate_coloc_results(dir_out = dir_out, verbose = verbose)

  return(invisible(NULL))
}

#' Process a single job of colocalization analyses
#' 
#' @description
#' Processes a subset of secondary datasets against a primary dataset for a 
#' specific set of regions. Manages temporary file creation, parallel processing
#' of datasets, and archiving of results.
#' 
#' @param sumstats_1_form Formatted primary dataset for this job's regions.
#' @param args_df Data frame subset with secondary datasets to process.
#' @param tar_sumstats_file Path to tar archive for storing sumstats.
#' @param tar_coloc_file Path to tar archive for storing coloc results.
#' @param mc_cores Integer. Cores for parallel processing (default: 10).
#' @param verbose Logical. Print progress messages (default: TRUE).
#' @param debug_mode Logical. Run sequentially (default: FALSE).
#' @param save_sumstats Logical. Whether to save filtered summary statistics 
#'   (default: FALSE). Set to TRUE only if you need the filtered summary statistics
#'   for downstream analyses like fine-mapping. WARNING: This substantially increases
#'   disk usage (~20-30GB per 1,00 regions and 10,000 datasets with typical settings).
#'   When FALSE, only colocalization results are saved (~2-3GB).
#' @param p_filt Numeric. Maximum p-value threshold for variants to include in 
#'   saved summary statistics (default: 1, includes all variants). Only applies 
#'   when save_sumstats = TRUE. Set to a lower value (e.g., 0.05) to reduce 
#'   storage by excluding variants with high p-values. Note: For fine-mapping 
#'   and conditional analyses, keep at 1 to preserve complete LD structure. 
#' @param p_min_save Numeric. Minimum p-value threshold for saving summary 
#'   statistics sumstats_2 (default: 5e-8). Only regions with at least one variant 
#'   with p-value below this threshold will be saved when save_sumstats is TRUE
#' @return List with two elements:
#'   \itemize{
#'     \item n_sumstats: Number of sumstats files processed
#'     \item n_coloc: Number of coloc result files created
#'   }
#' 
#' @details
#' This function:
#' \enumerate{
#'   \item Creates temporary directories for intermediate files
#'   \item Processes each secondary dataset via genepicoloc_run
#'   \item Archives results into tar files
#'   \item Cleans up temporary directories
#' }
#' 
#' The function uses mapply/mcmapply for parallel processing based on
#' debug_mode setting. Results are appended to existing tar archives,
#' allowing incremental processing across multiple calls.
#' 
#' @examples
#' \dontrun{
#' result <- genepicoloc_job(
#'   sumstats_1_form = formatted_primary_data,
#'   args_df = args_subset,
#'   tar_sumstats_file = "output/sumstats.tar",
#'   tar_coloc_file = "output/coloc.tar",
#'   mc_cores = 4
#' )
#' 
#' # Check processing status
#' if (result$n_sumstats != nrow(args_subset)) {
#'   warning("Some datasets failed to process")
#' }
#' }
#' 
#' @export
genepicoloc_job <- function(sumstats_1_form,
                            args_df,
                            tar_sumstats_file,
                            tar_coloc_file,
                            mc_cores = 10,
                            verbose = TRUE,
                            debug_mode = FALSE,
                            save_sumstats = FALSE,
                            p_filt = 0.05,
                            p_min_save = 5e-8) {
  
  # Create temp directories
  temp_coloc_dir <- tempfile(pattern = "coloc_", tmpdir = tempdir())
  dir.create(temp_coloc_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (save_sumstats) {
    temp_sumstats_dir <- tempfile(pattern = "sumstats_", tmpdir = tempdir())
    dir.create(temp_sumstats_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    temp_sumstats_dir <- NULL  # Pass NULL when not saving
  }
  
  # Register cleanup immediately after creation
  on.exit({
    if (save_sumstats && dir.exists(temp_sumstats_dir)) {
      unlink(temp_sumstats_dir, recursive = TRUE)
    }
    if (dir.exists(temp_coloc_dir)) unlink(temp_coloc_dir, recursive = TRUE)
  }, add = TRUE)
  
  # Prepare arguments for mapply
  mapply_args <- list(
    FUN = genepicoloc_run,
    sumstats_2_study = args_df$sumstats_2_study,
    sumstats_2_file = args_df$sumstats_2_file,
    sumstats_2_function = args_df$sumstats_2_function,
    sumstats_2_type = args_df$sumstats_2_type,
    sumstats_2_sdY = args_df$sumstats_2_sdY,
    MoreArgs = list(
      temp_sumstats_dir = temp_sumstats_dir,
      temp_coloc_dir = temp_coloc_dir,
      sumstats_1_form = sumstats_1_form,
      verbose = verbose,
      save_sumstats = save_sumstats,
      p_filt = p_filt,
      p_min_save = p_min_save
    ),
    SIMPLIFY = FALSE
  )
  
  # Execute analysis
  if (debug_mode) {
    results <- do.call(mapply, mapply_args)
  } else {
    results <- do.call(parallel::mcmapply, c(mapply_args, list(mc.cores = mc_cores)))
  }
  
  # Archive the results - only handle sumstats if saving
  if (save_sumstats) {
    tar_cmd_sumstats <- sprintf(
      "find '%s' -maxdepth 2 -type f -printf '%%P\\0' | tar -rf '%s' --null -C '%s' -T -",
      temp_sumstats_dir, tar_sumstats_file, temp_sumstats_dir
    )
    exit_code_sumstats <- system(tar_cmd_sumstats)
    if (exit_code_sumstats != 0) {
      stop(sprintf("Failed to add sumstats to tar archive (exit code: %d)", exit_code_sumstats))
    }
  }
  
  # Always archive coloc results
  tar_cmd_coloc <- sprintf(
    "find '%s' -maxdepth 2 -type f -printf '%%P\\0' | tar -rf '%s' --null -C '%s' -T -",
    temp_coloc_dir, tar_coloc_file, temp_coloc_dir
  )
  exit_code_coloc <- system(tar_cmd_coloc)
  if (exit_code_coloc != 0) {
    stop(sprintf("Failed to add coloc results to tar archive (exit code: %d)", exit_code_coloc))
  }
  
  # Return number of processed files for validation
  return(list(
    n_sumstats = if (save_sumstats) length(list.files(temp_sumstats_dir)) else 0,
    n_coloc = length(list.files(temp_coloc_dir, recursive = T))
  ))
  
}


#' Process secondary summary statistics and run colocalization
#'
#' @description
#' Formats secondary summary statistics and runs colocalization analysis against
#' a primary dataset across all specified regions. Handles both single and 
#' multi-phenotype datasets.
#'
#' @param temp_sumstats_dir Path to temporary directory for sumstats output.
#' @param temp_coloc_dir Path to temporary directory for coloc results.
#' @param sumstats_1_form Formatted primary dataset with coloc_regions_PASS attribute.
#' @param sumstats_2_function Function name for retrieving secondary data.
#' @param sumstats_2_file Path to secondary summary statistics file.
#' @param sumstats_2_type Type of secondary trait ('quant' or 'cc').
#' @param sumstats_2_sdY Standard deviation for quantitative traits (NA for cc).
#' @param sumstats_2_study Study identifier for tracking.
#' @param p_filt P-value threshold for filtering (default: 0.05).
#' @param p_min_save P-value threshold to save regions (default: 5e-8).
#' @param verbose Logical. Print progress messages (default: TRUE).
#' @param save_sumstats Logical. Whether to save filtered summary statistics 
#'   (default: FALSE). Set to TRUE only if you need the filtered summary statistics
#'   for downstream analyses like fine-mapping. WARNING: This substantially increases
#'   disk usage (~20-30GB per 1,00 regions and 10,000 datasets with typical settings).
#'   When FALSE, only colocalization results are saved (~2-3GB).
#' @param p_filt Numeric. Maximum p-value threshold for variants to include in 
#'   saved summary statistics (default: 1, includes all variants). Only applies 
#'   when save_sumstats = TRUE. Set to a lower value (e.g., 0.05) to reduce 
#'   storage by excluding variants with high p-values. Note: For fine-mapping 
#'   and conditional analyses, keep at 1 to preserve complete LD structure. 
#' @param p_min_save Numeric. Minimum p-value threshold for saving summary 
#'   statistics sumstats_2 (default: 5e-8). Only regions with at least one variant 
#'   with p-value below this threshold will be saved when save_sumstats is TRUE
#' @return Invisibly returns the sumstats_2_file path.
#'
#' @details
#' Processing steps:
#' \enumerate{
#'   \item Formats secondary summary statistics via format_sumstats_2
#'   \item For single phenotype: runs colocalization directly
#'   \item For multiple phenotypes: processes each phenotype separately
#'   \item Saves filtered sumstats (significant regions only)
#'   \item Saves colocalization results as RDS files
#' }
#'
#' The function automatically detects multi-phenotype datasets (e.g., eQTL)
#' and processes each phenotype independently, combining results afterward.
#'
#' @examples
#' \dontrun{
#' genepicoloc_run(
#'   temp_sumstats_dir = "/tmp/sumstats",
#'   temp_coloc_dir = "/tmp/coloc",
#'   sumstats_1_form = primary_data,
#'   sumstats_2_function = "retrieve_sumstats_tabix",
#'   sumstats_2_file = "secondary.gz",
#'   sumstats_2_type = "quant",
#'   sumstats_2_sdY = 1.2,
#'   sumstats_2_study = "StudyA"
#' )
#' }
#'
#' @export
genepicoloc_run <- function(temp_sumstats_dir,
                            temp_coloc_dir,
                            sumstats_1_form,
                            sumstats_2_function,
                            sumstats_2_file,
                            sumstats_2_type,
                            sumstats_2_sdY,
                            sumstats_2_study,
                            verbose = TRUE,
                            save_sumstats = FALSE,
                            p_filt = 0.05,
                            p_min_save = 5e-8) {
  
  # Extract regions from primary dataset
  # NOTE: This will be a subset of the original regions when processing jobs
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
  
  # Calculate filtering and min_save p-value thresholds
  nlog10P_filt <- -log10(p_filt)
  nlog10P_min_save <- -log10(p_min_save)
  
  # Run colocalization based on phenotype type
  if (attr(sumstats_2_form, "sumstats_pheno") == "single") {
    ### Single phenotype analysis
    # Save sumstats (only significant regions)
    if (save_sumstats) {
      if (p_filt < 1) sumstats_2_form <- subset(sumstats_2_form, nlog10P > nlog10P_filt)
      sumstats_2_filt <- filter_significant_regions(sumstats_2_form, nlog10P_min_save)
      if (nrow(sumstats_2_filt) > 0) {
        if (!dir.exists(file.path(temp_sumstats_dir, sumstats_2_study))) {
          dir.create(file.path(temp_sumstats_dir, sumstats_2_study))
        }
        saveRDS(sumstats_2_filt,
                file.path(temp_sumstats_dir, sumstats_2_study,
                          paste0(basename(attr(sumstats_2_form, "sumstats_file")), ".RDS")))
      }
    }
    # Run colocalization
    coloc_results <- run_single_phenotype(
      sumstats_1 = sumstats_1_form,
      sumstats_2 = sumstats_2_form,
      coloc_regions_PASS = coloc_regions_PASS
    )
    # Combine results from all regions
    coloc_results <- data.table::rbindlist(coloc_results, fill = TRUE)
    
  } else if (attr(sumstats_2_form, "sumstats_pheno") == "multiple") {
    ### Multiple phenotype analysis
    # Check if the list is empty or all elements are empty
    if (length(sumstats_2_form) == 0 || 
        all(sapply(sumstats_2_form, function(x) nrow(x) == 0))) {
      # Handle empty multi-phenotype case
      coloc_results <- data.table::data.table()
    } else {
      coloc_results <- lapply(sumstats_2_form, function(sumstats_2) {
        if (nrow(sumstats_2) == 0) return(NULL) # Skip empty phenotypes
        # Save sumstats (only significant regions)
        if (save_sumstats) {
          if (p_filt < 1) sumstats_2 <- subset(sumstats_2, nlog10P > nlog10P_filt)
          sumstats_2_filt <- filter_significant_regions(sumstats_2, nlog10P_min_save)
          if (nrow(sumstats_2_filt) > 0) {
            if (!dir.exists(file.path(temp_sumstats_dir, sumstats_2_study))) {
              dir.create(file.path(temp_sumstats_dir, sumstats_2_study))
            }
            saveRDS(sumstats_2_filt, 
                    file.path(temp_sumstats_dir, sumstats_2_study,
                              paste0(basename(attr(sumstats_2, "sumstats_file")), ".RDS")))
          }
        }
        # Run colocalization
        region_results <- run_single_phenotype(
          sumstats_1 = sumstats_1_form,
          sumstats_2 = sumstats_2,
          coloc_regions_PASS = attr(sumstats_2, "coloc_regions_PASS")
        )
        # Combine regions for this phenotype
        return(data.table::rbindlist(region_results, fill = TRUE))
      })
      
      # Remove NULL results (from empty phenotypes)
      coloc_results <- coloc_results[!sapply(coloc_results, is.null)]
      
      # Combine results from all phenotypes, from all regions
      if (length(coloc_results) > 0) {
        coloc_results <- data.table::rbindlist(coloc_results, fill = TRUE)
      } else {
        coloc_results <- data.table::data.table()
      }
    }
  } else {
    stop("Attribute 'sumstats_pheno' of sumstats_2_form can only be 'single' or 'multiple'")
  }
  
  # Save colocalization results (even if empty, for completeness)
  if (!dir.exists(file.path(temp_coloc_dir, sumstats_2_study))) {
    dir.create(file.path(temp_coloc_dir, sumstats_2_study))
  }
  saveRDS(coloc_results, file.path(temp_coloc_dir, sumstats_2_study,
                                   paste0(basename(attr(sumstats_2_form, "sumstats_file")), ".RDS")))
  
  invisible(sumstats_2_file)
  
}

  
# Group 1: Format and process functions ----

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
#'     \item A list of sumstats objects (if sumstats_pheno = "multiple"), where
#'           each element has its coloc_regions_PASS attribute filtered to regions
#'           containing data for that specific phenotype
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
#' \code{\link{genepicoloc_run}} which calls this function
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
#' \code{\link{genepicoloc_run}} which calls this function
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
  
  # Run analysis for each region
  coloc_results <- mapply(
    FUN = run_region,
    CHR_var = coloc_regions_PASS$CHR_var,
    BP_START_var = coloc_regions_PASS$BP_START_var,
    BP_STOP_var = coloc_regions_PASS$BP_STOP_var,
    MoreArgs = list(sumstats_1 = sumstats_1, sumstats_2 = sumstats_2),
    SIMPLIFY = FALSE
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
  
  # Extract region-specific data
  sumstats_1_sub <- subset(sumstats_1, 
                           CHR == CHR_var & POS > BP_START_var & POS < BP_STOP_var)
  
  # Check primary dataset has data
  if (nrow(sumstats_1_sub) == 0) {
    stop("No data found in sumstats_1 for region ", 
         CHR_var, ":", BP_START_var, "-", BP_STOP_var)
  }
  
  sumstats_2_sub <- subset(sumstats_2, 
                           CHR == CHR_var & POS > BP_START_var & POS < BP_STOP_var)
  
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
  
  # Add directionality results to output except for "sumstats_2_ind_A2"
  for (col in setdiff(colnames(directionality_df), "sumstats_2_ind_A2")) {
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
#'     \item sumstats_2_ind_A1/A2: Alleles in secondary dataset
#'     \item sumstats_2_ind_nlog10P: Significance in secondary dataset
#'     \item sumstats_2_ind_BETA: Effect size in secondary dataset
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
    sumstats_1_ind_nlog10P = lead_snp_1$nlog10P,
    sumstats_1_ind_BETA = lead_snp_1$BETA,
    sumstats_2_ind_A1 = lead_snp_2$A1,
    sumstats_2_ind_A2 = lead_snp_2$A2,
    sumstats_2_ind_nlog10P = lead_snp_2$nlog10P,
    sumstats_2_ind_BETA = lead_snp_2$BETA,
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
#'   sumstats_1_ind_A1, sumstats_2_ind_A1, sumstats_2_ind_A2
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
#'   sumstats_2_ind_A2 = "G"
#' )
#' detect_directionality(df)  # Returns df with directionality = -1
#' }
#' 
#' @export
detect_directionality <- function(ind_df) {
  
  # Input validation
  required_cols <- c("sumstats_1_ind_A1", "sumstats_2_ind_A1", 
                     "sumstats_2_ind_A2")
  
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
  sign1 <- sign(ind_df$sumstats_1_ind_BETA)
  sign2 <- sign(ind_df$sumstats_2_ind_BETA)
  
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

#' Filter summary statistics to regions with significant variants
#' 
#' @description 
#' Filters a formatted summary statistics object to retain only genomic regions 
#' that contain at least one variant exceeding a significance threshold. This is 
#' useful for reducing computational burden by excluding regions unlikely to show 
#' colocalization signals. The function handles overlapping regions correctly by 
#' retaining all variants within any significant region.
#' 
#' @param sumstats_form A formatted summary statistics object of class 'sumstats_form' 
#'   containing columns CHR, POS, and nlog10P. Must have a 'coloc_regions_PASS' 
#'   attribute defining the genomic regions to test.
#' @param nlog10P_min_save Numeric. Minimum -log10(p-value) threshold for considering 
#'   a variant as significant. Default is 6 (corresponding to p < 1e-6). Regions 
#'   containing at least one variant exceeding this threshold are retained.
#' 
#' @return A filtered summary statistics object containing only variants within 
#'   significant regions. Returns an empty data frame with preserved structure if 
#'   no regions meet the significance threshold. The 'coloc_regions_PASS' attribute 
#'   is updated to reflect only the retained regions.
#' 
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Identifies regions containing at least one variant with nlog10P > threshold
#'   \item Collects all variants within these significant regions (handles overlaps)
#'   \item Updates the coloc_regions_PASS attribute to match filtered regions
#' }
#' 
#' Regions may overlap, so the function ensures that all variants belonging to any 
#' significant region are retained, even if they also belong to non-significant regions.
#' 
#' @examples
#' \dontrun{
#' # Filter to keep only regions with suggestive GWAS signals (p < 1e-6)
#' sumstats_filtered <- filter_significant_regions(sumstats_form, nlog10P_min_save = 6)
#' 
#' # Filter to keep only regions with genome-wide significant signals (p < 5e-8)
#' sumstats_filtered <- filter_significant_regions(sumstats_form, nlog10P_min_save = 7.3)
#' }
filter_significant_regions <- function(sumstats_form, nlog10P_min_save = 6) {
  
  # Get regions
  coloc_regions_PASS <- attr(sumstats_form, "coloc_regions_PASS")
  
  # Convert CHR to character if needed for consistency
  if (is.numeric(sumstats_form$CHR)) {
    sumstats_form$CHR <- as.character(sumstats_form$CHR)
  }
  
  # Determine regions with at least one significant SNP
  signif_regions <- sapply(1:nrow(coloc_regions_PASS), function(x) {
    CHR_var <- coloc_regions_PASS[x,]$CHR_var
    BP_START_var <- coloc_regions_PASS[x,]$BP_START_var
    BP_STOP_var <- coloc_regions_PASS[x,]$BP_STOP_var
    any(sumstats_form$CHR == CHR_var & 
          sumstats_form$POS >= BP_START_var & 
          sumstats_form$POS <= BP_STOP_var & 
          sumstats_form$nlog10P > nlog10P_min_save)
  })
  
  # Keep only significant regions
  coloc_regions_PASS_sign <- coloc_regions_PASS[signif_regions,]
  
  # Handle case with no significant regions
  if (nrow(coloc_regions_PASS_sign) == 0) {
    sumstats_filtered <- sumstats_form[0,]  # Fixed: was referencing undefined variable
    attr(sumstats_filtered, "coloc_regions_PASS") <- coloc_regions_PASS_sign
    return(sumstats_filtered)
  }
  
  # Regions might intersect
  # Keep all SNPs that belong to regions with at least one significant SNP
  # Use lapply to ensure we always get a list, then unlist to vector
  idx_list <- lapply(1:nrow(coloc_regions_PASS_sign), function(x) {
    CHR_var <- coloc_regions_PASS_sign[x,]$CHR_var
    BP_START_var <- coloc_regions_PASS_sign[x,]$BP_START_var
    BP_STOP_var <- coloc_regions_PASS_sign[x,]$BP_STOP_var
    which(
      sumstats_form$CHR == CHR_var & 
        sumstats_form$POS >= BP_START_var & 
        sumstats_form$POS <= BP_STOP_var)
  })
  idx_to_keep <- as.numeric(unique(unlist(idx_list)))
  
  # Handle empty index case
  if (length(idx_to_keep) == 0) {
    sumstats_filtered <- sumstats_form[0,]
  } else {
    sumstats_filtered <- sumstats_form[idx_to_keep,]
  }

  # Update the regions attribute
  attr(sumstats_filtered, "coloc_regions_PASS") <- coloc_regions_PASS_sign
  
  return(sumstats_filtered)
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
    nearest_gene_1 = NA_character_,
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
    sumstats_1_ind_nlog10P = NA_real_,
    sumstats_1_ind_BETA = NA_real_,
    sumstats_2_ind_A1 = NA_character_,
    sumstats_2_ind_nlog10P = NA_real_,
    sumstats_2_ind_BETA = NA_real_,
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
#' @param nearest Numeric vector. Number of nearest genes to annotate for each region
#'   (default: c(1, 10)). Creates one column per value.
#'
#' @return The input data.frame with added gene annotation columns:
#'   \itemize{
#'     \item region_center_pos: Center position of the region
#'     \item nearest_genes_N: N nearest genes for each value in 'nearest' parameter
#'   }
#'   For example, with nearest = c(1, 10), creates nearest_genes_1 and nearest_genes_10.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Calculates the center position of the genomic region
#'   \item Finds the specified number of nearest genes for each value in 'nearest'
#'   \item Uses appropriate format: detailed for single gene, simple for multiple
#' }
#'
#' @seealso 
#' \code{\link{annotate_position}} for the underlying annotation function
#' \code{\link{setup_gene_annotation}} for gene map preparation
#'
#' @examples
#' \dontrun{
#' # Default: 1 and 10 nearest genes
#' result_annotated <- add_gene_annotations(result, gene_setup$gene_map)
#' 
#' # Custom: 1, 5, and 20 nearest genes
#' result_annotated <- add_gene_annotations(result, gene_setup$gene_map, 
#'                                          nearest = c(1, 5, 20))
#' }
#'
#' @export
add_gene_annotations <- function(result, gene_map = NULL, verbose = FALSE,
                                 nearest = c(1, 10)) {
  
  # Input validation
  if (!is.data.frame(result)) {
    stop("result must be a data.frame")
  }
  
  required_cols <- c("CHR_var", "BP_START_var", "BP_STOP_var")
  missing_cols <- setdiff(required_cols, colnames(result))
  if (length(missing_cols) > 0) {
    stop("result must contain columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Validate nearest parameter
  if (!is.numeric(nearest) || length(nearest) == 0) {
    stop("nearest must be a numeric vector with at least one value")
  }
  nearest <- sort(unique(as.integer(nearest)))
  if (any(nearest < 1)) {
    stop("nearest values must be >= 1")
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
  
  # Process each row
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
    
    # Annotate for each specified number of nearest genes
    tryCatch({
      for (n in nearest) {
        # Use default format for single gene, simple format for multiple
        output_format <- if (n == 1) "default" else "simple"
        
        annotation <- annotate_position(
          chr = chr, 
          pos = center_pos, 
          gene_map = gene_map,
          n_nearest = n,
          output_format = output_format
        )
        
        # Create column name
        col_name <- if (n == 1) "nearest_gene_1" else paste0("nearest_genes_", n)
        result[[col_name]][i] <- annotation
        
        if (verbose && n %in% c(1, max(nearest))) {
          message("  ", n, " nearest: ", annotation)
        }
      }
      
    }, error = function(e) {
      warning("Gene annotation failed for region ", chr, ":", start_pos, "-", stop_pos, 
              ": ", e$message)
      # Create NA columns for failed annotations
      for (n in nearest) {
        col_name <- paste0("nearest_genes_", n)
        result[[col_name]][i] <- NA_character_
      }
    })
  }
  
  return(result)
}

#' Save job metadata and summary information
#' 
#' @description Helper function to save comprehensive job metadata after completion
#' 
#' @param dir_out Output directory path
#' @param args_df Arguments data frame
#' @param n_regions Total number of regions
#' @param max_regions_per_job Maximum regions per job
#' @param n_jobs_1 Number of jobs created
#' @param mc_cores Number of cores used
#' @param study_counts Table of study counts
#' @param job_times Vector of job processing times
#' @param verbose Whether to print messages
#' 
#' @return Invisibly returns list of saved file paths
save_job_metadata <- function(dir_out, args_df, n_regions, max_regions_per_job, 
                              n_jobs_1, mc_cores, study_counts, job_times, 
                              verbose = TRUE) {
  
  # Define output files
  args_df_file <- file.path(dir_out, "args_df.csv.gz")
  job_info_file <- file.path(dir_out, "job_splitting_info.txt")
  job_summary_file <- file.path(dir_out, "job_summary.csv")
  
  # Console output
  if (verbose) {
    message("\n=== Final Summary ===")
    message(sprintf("Jobs completed successfully: %d", n_jobs_1))
    message(sprintf("Total processing time: %.1f minutes", sum(job_times)))
    message("\nSaving metadata:")
    message("  - Arguments: ", args_df_file)
    message("  - Job info: ", job_info_file)
    message("  - Job summary: ", job_summary_file)
  }
  
  # Save arguments
  data.table::fwrite(args_df, args_df_file, compress = "gzip")
  
  # Create detailed job information text file
  job_info <- c(
    "=== JOB SPLITTING INFORMATION ===",
    sprintf("Finished at: %s", Sys.time()),
    "",
    "=== REGIONS ===",
    sprintf("Total regions: %d", n_regions),
    sprintf("Regions per job: %d (max)", max_regions_per_job),
    sprintf("Number of jobs: %d", n_jobs_1),
    "",
    "=== DATASETS ===",
    sprintf("Total arguments (sumstats_2): %d", nrow(args_df)),
    "",
    "=== STUDY BREAKDOWN ===",
    unlist(lapply(names(study_counts), function(study) {
      sprintf("  - %s: %d datasets", study, study_counts[study])
    })),
    "",
    "=== PARALLELIZATION ===",
    sprintf("Cores per job: %d", mc_cores),
    sprintf("Subjobs per job: %d", ceiling(nrow(args_df) / mc_cores)),
    sprintf("Datasets per subjob: %d (max)", mc_cores),
    "",
    "=== JOB DETAILS ==="
  )
  
  # Add per-job details
  for (i in 1:n_jobs_1) {
    start_idx <- (i - 1) * max_regions_per_job + 1
    end_idx <- min(i * max_regions_per_job, n_regions)
    n_regions_in_job <- end_idx - start_idx + 1
    
    job_info <- c(job_info,
                  sprintf("Job %04d: regions %d-%d (n=%d), %d datasets, time=%.1f min",
                          i, start_idx, end_idx, n_regions_in_job, nrow(args_df), job_times[i])
    )
  }
  
  job_info <- c(job_info,
                "",
                "=== OUTPUT STRUCTURE ===",
                sprintf("Job directories: jobs/job_0001 to jobs/job_%04d", n_jobs_1),
                "Each job contains:",
                "  - sumstats_1.RDS: Primary dataset for job's regions",
                "  - sumstats.tar: All secondary datasets (filtered)",
                "  - coloc.tar: Colocalization results",
                "",
                "=== TOTAL OPERATIONS ===",
                sprintf("Total colocalizations: %d regions x %d datasets = %d analyses",
                        n_regions, nrow(args_df), n_regions * nrow(args_df))
  )
  
  # Write to file
  writeLines(job_info, job_info_file)
  
  # Display the job info to console if verbose
  if (verbose) {
    message("\n", paste(job_info, collapse = "\n"))
  }
  
  # Create CSV with job details for easier parsing
  job_summary_df <- data.frame(
    job_id = sprintf("job_%04d", 1:n_jobs_1),
    start_region = sapply(1:n_jobs_1, function(i) (i-1) * max_regions_per_job + 1),
    end_region = sapply(1:n_jobs_1, function(i) min(i * max_regions_per_job, n_regions)),
    n_regions = sapply(1:n_jobs_1, function(i) {
      end_idx <- min(i * max_regions_per_job, n_regions)
      start_idx <- (i - 1) * max_regions_per_job + 1
      end_idx - start_idx + 1
    }),
    n_datasets = nrow(args_df),
    n_analyses = sapply(1:n_jobs_1, function(i) {
      end_idx <- min(i * max_regions_per_job, n_regions)
      start_idx <- (i - 1) * max_regions_per_job + 1
      (end_idx - start_idx + 1) * nrow(args_df)
    }),
    processing_time_min = job_times
  )
  
  data.table::fwrite(job_summary_df, job_summary_file)
  
  invisible(list(
    args_df = args_df_file,
    job_info = job_info_file,
    job_summary = job_summary_file
  ))
}

#' Consolidate colocalization results by study using job metadata
#' 
#' @description
#' Extracts and combines colocalization results from job-specific tar archives
#' into study-specific RDS files, using the saved job metadata.
#' 
#' @param dir_out Base directory containing job subdirectories and metadata files
#' @param output_subdir Subdirectory name for consolidated results (default: "consolidated")
#' @param verbose Print progress messages (default: TRUE)
#' 
#' @return Invisibly returns list of created output files
#' 
#' @export
consolidate_coloc_results <- function(dir_out, output_subdir = "unfilt", verbose = TRUE) {
  
  # Read metadata files
  args_df_file <- file.path(dir_out, "args_df.csv.gz")
  job_summary_file <- file.path(dir_out, "job_summary.csv")
  
  if (!file.exists(args_df_file)) {
    stop("args_df.csv.gz not found in ", dir_out)
  }
  if (!file.exists(job_summary_file)) {
    stop("job_summary.csv not found in ", dir_out)
  }
  
  args_df <- data.table::fread(args_df_file)
  job_summary <- data.table::fread(job_summary_file)
  
  # Get unique studies
  studies <- unique(args_df$sumstats_2_study)
  if (verbose) {
    message("Found ", length(studies), " studies to process")
    message("Found ", nrow(job_summary), " jobs to consolidate")
  }
  
  # Create output directory
  dir_output <- file.path(dir_out, output_subdir)
  if (!dir.exists(dir_output)) dir.create(dir_output, recursive = TRUE)
  
  # Process each study
  output_files <- character()
  
  for (study in studies) {
    if (verbose) message("\nProcessing study: ", study)
    
    # Get expected files for this study
    study_files <- args_df[args_df$sumstats_2_study == study, basename(sumstats_2_file)]
    expected_count <- length(study_files) * nrow(job_summary)
    
    # Temporary directory for extraction
    temp_dir <- tempfile(pattern = paste0("coloc_", study, "_"))
    dir.create(temp_dir, recursive = TRUE)
    on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
    
    # Extract and collect results from each job
    all_results <- list()
    file_count <- 0
    
    for (i in 1:nrow(job_summary)) {
      job_id <- job_summary$job_id[i]
      coloc_tar <- file.path(dir_out, "jobs", job_id, "coloc.tar")
      
      if (!file.exists(coloc_tar)) {
        warning("coloc.tar not found for ", job_id)
        next
      }
      
      # Extract files for this study
      extract_cmd <- sprintf(
        "tar -xf '%s' -C '%s' '%s/' 2>/dev/null",
        coloc_tar, temp_dir, study
      )
      
      exit_code <- system(extract_cmd, ignore.stdout = TRUE)
      
      # Read extracted RDS files
      study_dir <- file.path(temp_dir, study)
      if (dir.exists(study_dir)) {
        rds_files <- list.files(study_dir, pattern = "\\.RDS$", full.names = TRUE)
        
        for (rds_file in rds_files) {
          result <- readRDS(rds_file)

          all_results[[length(all_results) + 1]] <- result
          file_count <- file_count + 1
        }
        
        # Clean up extracted files for this job
        unlink(study_dir, recursive = TRUE)
      }
    }
    
    # Combine all results
    if (length(all_results) > 0) {
      combined_results <- data.table::rbindlist(all_results, fill = TRUE)

      # Save consolidated results
      output_file <- file.path(dir_output, paste0(study, "_unfilt.RDS"))
      saveRDS(combined_results, output_file)
      output_files <- c(output_files, output_file)
      
      if (verbose) {
        message(sprintf("  - Processed %d/%d expected files", 
                        file_count, expected_count))
        message(sprintf("  - Combined %d rows into %s", 
                        nrow(combined_results), basename(output_file)))
      }
      
      # Validation warning
      if (file_count != expected_count) {
        warning(sprintf("Study %s: Expected %d files (%.0f files x %.0f jobs), found %d", 
                        study, expected_count, length(study_files), nrow(job_summary), file_count))
      }
      
    } else {
      warning("No results found for study: ", study)
    }
  }
  
  if (verbose) {
    message("\n=== Consolidation Summary ===")
    message("Created ", length(output_files), " study files in ", dir_output)
    message("\nYou can load results with:")
    message("  results <- readRDS('", output_files[1], "')")
  }
  
  invisible(output_files)
}
