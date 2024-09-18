#' make_eQTL_Catalogue_args
#' @export
make_eQTL_Catalogue_args <- function(QTS="QTS000015") {
  datasets_eQTL_Catalogue <- get_datasets_eQTL_Catalogue()
  datasets_eQTL_Catalogue <- subset(datasets_eQTL_Catalogue,
                                    quant_method == "ge" & study_id %in% QTS)
  eQTL_Catalogue <- list()
  eQTL_Catalogue$sumstats_2_files <- datasets_eQTL_Catalogue$dataset_id
  eQTL_Catalogue$sumstats_2_function <- "query_eQTL_Catalogue"
  eQTL_Catalogue$sumstats_2_type <- "quant"
  eQTL_Catalogue$sumstats_2_sdY <- NA
  eQTL_Catalogue
}

#' create_args_df
#' create a data.frame with arguments to be passed to coloc_wrapper
#' @export
create_args_df <- function(CHR_var, BP_START_var, BP_STOP_var,
                           sumstats_1_file, sumstats_1_function,
                           sumstats_1_type, sumstats_1_sdY,
                           sumstats_2_files, sumstats_2_function,
                           sumstats_2_type, sumstats_2_sdY) {
  region_args <- data.frame(CHR_var = CHR_var,
                            BP_START_var = BP_START_var,
                            BP_STOP_var = BP_STOP_var,
                            stringsAsFactors = F)
  sumstats_1_args <- data.frame(sumstats_1_file=sumstats_1_file,
                                sumstats_1_function=sumstats_1_function,
                                sumstats_1_type=sumstats_1_type,
                                sumstats_1_sdY=sumstats_1_sdY,
                                stringsAsFactors = F)
  sumstats_2_args <- data.frame(sumstats_2_file=sumstats_2_files,
                                sumstats_2_function=sumstats_2_function,
                                sumstats_2_type=sumstats_2_type,
                                sumstats_2_sdY=sumstats_2_sdY,
                                stringsAsFactors = F)
  args_df <- Reduce(merge, list(region_args, sumstats_1_args, sumstats_2_args))
  return(args_df)
}



