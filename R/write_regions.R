#' Save coloc regions
#' @export
write_regions <- function(coloc_regions_list, sumstats_name,
                          selected_studies = NULL,
                          SKIP_name="Name", CHR_place=3, POS_place=4,
                          bgzip_bin="bgzip", tabix_bin="tabix",
                          regions_log="regions_log", coloc_regions="coloc_regions",
                          coloc_regions_PASS="coloc_regions_PASS",
                          sumstats_filt="sumstats_filt") {
  which_na <- is.na(coloc_regions_list)
  if(!which_na[regions_log]) {
    writeLines(coloc_regions_list[[regions_log]], con = paste0(sumstats_name, "_log.txt"))
    message(paste0("Written: ", paste0(sumstats_name, "_log.txt")))
  }
  if(!which_na[coloc_regions]){
    write.table(coloc_regions_list[[coloc_regions]], paste0(sumstats_name, "_", coloc_regions, ".tsv"),
                sep="\t", row.names = F, col.names = T, quote = F)
    message(paste0("Written: ", paste0(sumstats_name, "_", coloc_regions, ".tsv")))
  }
  if(!which_na[coloc_regions_PASS]){
    write.table(coloc_regions_list[[coloc_regions_PASS]], paste0(sumstats_name, "_", coloc_regions_PASS, ".tsv"),
                sep="\t", row.names = F, col.names = T, quote = F)
    message(paste0("Written: ", paste0(sumstats_name, "_", coloc_regions_PASS, ".tsv")))
  }
  if(!which_na[sumstats_filt]){
    write.table(coloc_regions_list[[sumstats_filt]], paste0(sumstats_name, "_subset.tsv"),
                sep="\t", row.names = F, col.names = T, quote = F)
    system(paste0(bgzip_bin, " -f ", sumstats_name, "_subset.tsv"))
    system(paste0(tabix_bin, " -f -s", CHR_place, " -b", POS_place, " -e", POS_place, " ", sumstats_name, "_subset.tsv.gz -c ", SKIP_name))
    message(paste0("Written: ", paste0(sumstats_name, "_subset.tsv.gz and tbi")))
  }
  sumstats_1_files <- paste0(sumstats_name, "_subset.tsv.gz")
  return(sumstats_1_files)
}

