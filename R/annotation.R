create_sumstats_1_args <- function(sumstats_1_file,
                                   sumstats_1_function="query_sumstats_1",
                                   sumstats_1_type,
                                   sumstats_1_sdY) {
  sumstats_1_args <- list()
  sumstats_1_args$sumstats_1_file <- sumstats_1_file
  sumstats_1_args$sumstats_1_function <- sumstats_1_function
  sumstats_1_args$sumstats_1_type <- sumstats_1_type
  sumstats_1_args$sumstats_1_sdY <- sumstats_1_sdY
  return(sumstats_1_args)
}

#' annotate_eQTL_Catalog
annotate_eQTL_Catalog <- function(coloc_out, annotation_df,
                                  datasets_eQTL_Catalogue) {
  coloc_out$dataset_id <- gsub("(.*)_.*", "\\1", coloc_out$sumstats_2_file)
  coloc_out$gene_id <- gsub(".*_(.*)", "\\1", coloc_out$sumstats_2_file)
  coloc_out <- merge(coloc_out, annotation_df, all.x=T)[, union(names(coloc_out), names(annotation_df))]
  coloc_out <- merge(coloc_out, datasets_eQTL_Catalogue, all.x=T)[, union(names(coloc_out), names(datasets_eQTL_Catalogue))]
}
