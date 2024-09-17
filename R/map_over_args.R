#' map_over_args
#' @importFrom parallel mclapply
map_over_args <- function(args_df, mc_cores=10, dry_run=F,
                          debug_mode=F, verbose=T, do_rbind=T, save_tmp=F) {
  if (debug_mode) { mc_cores <- 1; cat("\n**\nDebugging mode\n**\n\n") }
  # eQTL_Catalogue - limit to 1 process to avoid too many API requests
  if (any(args_df$sumstats_2_function == "query_eQTL_Catalogue")) {
    message("eQTL_Catalogue - limit to 1 process to avoid too many API requests")
    mc_cores <- 1
  }
  if (dry_run) {
    max_row <- ifelse (nrow(args_df) < 10, nrow(args_df), 10)
    args_df <- args_df[1:max_row,]
  }
  # map over args_df
  vec_to_map <- 1:nrow(args_df)
  if (verbose) {
    when_message <- sapply(split(vec_to_map, sort(vec_to_map%%100)), function(x) x[1])
    gc_timer(time_start=Sys.time(), nrow_args_df=nrow(args_df),
             mc_cores=mc_cores, sumstats_2_file=args_df$sumstats_2_file[1])
  }
  coloc_out <- parallel::mclapply(vec_to_map, function(i) {
    if (verbose) if (i %in% when_message) cat("|")
    c_out <- do.call(coloc_wrapper, c(args_df[i,]))
    if (save_tmp) saveRDS(c_out, paste0(tempfile(), "_", basename(args_df[i,]$sumstats_2_file)))
    return(c_out)
  }, mc.cores = mc_cores)
  if (verbose) {
    cat(" Done.\n")
    time_stop <- Sys.time()
    message(paste0("Timestamp: ", time_stop))
    message("Time elapsed: ", round(as.numeric(time_stop-time_start,units="secs")), " seconds")
  }
  # rbind list to data.frame
  if (do_rbind) coloc_out <- do.call(rbind, coloc_out)
  return(coloc_out)
}

# helpers
#' gc_timer
gc_timer <- function(gc_speed=6, time_start, nrow_args_df, mc_cores, sumstats_2_file) {
  if (grepl("sQTLs", sumstats_2_file)) {
    gc_speed <- 0.05
    message("sQTLs processing takes substantially more time than other datasets (approximately x10), and the progress bar is not linear.")}
  if (grepl("QTS", sumstats_2_file)) {gc_speed <- 0.6}
  cat("-----------------------------\n")
  message(paste0("Timestamp: ", time_start))
  message(paste0("Total number of jobs: ", nrow_args_df))
  message(paste0("Number of requested cores: ", mc_cores))
  message(paste0("Approximate speed: ", gc_speed, " jobs per core per second"))
  message(paste0("Estimated time: ", round(nrow_args_df/mc_cores/gc_speed), " seconds"))
  progress_bar <- ifelse(nrow_args_df > 100, 100, nrow_args_df)
  cat(paste0(paste0(rep(" ", progress_bar), collapse = ""), "|100%\r"))
}

