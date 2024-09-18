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
