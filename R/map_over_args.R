#' map_over_args
#' @importFrom parallel mclapply
map_over_args <- function(args_df, mc_cores=10, dry_run=F,
                          debug_mode=F, verbose=T, do_rbind=T) {
  if (debug_mode) { mc_cores <- 1; cat("\n**\nDebugging mode\n**\n\n") }
  if (dry_run) {
    max_row <- ifelse (nrow(args_df) < 10, nrow(args_df), 10)
    args_df <- args_df[1:max_row,]
  }
  # map over args_df
  vec_to_map <- 1:nrow(args_df)
  if (verbose) {
    when_message <- sapply(split(vec_to_map, sort(vec_to_map%%100)), function(x) x[1])
    time_start <- Sys.time()
    gc_speed <- 6
    cat("-----------------------------\n")
    message(paste0("Timestamp: ", time_start))
    message(paste0("Total number of jobs: ", nrow(args_df)))
    message(paste0("Number of requested cores: ", mc_cores))
    message(paste0("Approximate speed: ", gc_speed, " jobs per core per second"))
    message(paste0("Estimated time: ", round(nrow(args_df)/mc_cores/gc_speed), " seconds"))
    message(paste0("Starting iterations, the % order can be custom due to parallelization:"))
  }
  coloc_out <- parallel::mclapply(vec_to_map, function(i) {
    if (verbose) if (i %in% when_message) {
      cat(paste0(names(when_message[when_message == i]), "%..")) }
    do.call(coloc_wrapper, c(args_df[i,]))
  }, mc.cores = mc_cores)
  if (verbose) {
    if (length(when_message) < 100) cat("......100%..")
    cat(" Done.\n")
    time_stop <- Sys.time()
    message(paste0("Timestamp: ", time_stop))
    message("Time elapsed: ", round(as.numeric(time_stop-time_start,units="secs")), " seconds")
  }
  # rbind list to data.frame
  if (do_rbind) coloc_out <- do.call(rbind, coloc_out)
  return(coloc_out)
}
