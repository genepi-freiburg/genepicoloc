#' create_args_df
#' create a data.frame with arguments to be passed to coloc_wrapper
create_args_df <- function(CHR_var, BP_START_var, BP_STOP_var,
                           sumstats_1_file, sumstats_1_function,
                           sumstats_1_type, sumstats_1_sdY,
                           sumstats_2_files, sumstats_2_function,
                           sumstats_2_type, sumstats_2_sdY) {
  region_args <- data.frame(CHR_var = CHR_var,
                            BP_START_var = BP_START_var,
                            BP_STOP_var = BP_STOP_var)
  sumstats_1_args <- data.frame(sumstats_1_file=sumstats_1_file,
                                sumstats_1_function=sumstats_1_function,
                                sumstats_1_type=sumstats_1_type,
                                sumstats_1_sdY=sumstats_1_sdY)
  sumstats_2_args <- data.frame(sumstats_2_file=sumstats_2_files,
                                sumstats_2_function=sumstats_2_function,
                                sumstats_2_type=sumstats_2_type,
                                sumstats_2_sdY=sumstats_2_sdY)
  args_df <- Reduce(merge, list(region_args, sumstats_1_args, sumstats_2_args))
  return(args_df)
}



