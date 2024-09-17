#' write_coloc_out
#' @importFrom writexl write_xlsx
#' @importFrom data.table rbindlist
#' @export
write_coloc_out <- function(coloc_out, study, output_folder,
                            remove_dirname = T,
                            PP.H4.abf_filt=0.5,
                            PP.H3.abf_filt=NULL) {
  if (remove_dirname) {
    coloc_out[["sumstats_1_file"]] <- basename(coloc_out[["sumstats_1_file"]])
    coloc_out[["sumstats_2_file"]] <- basename(coloc_out[["sumstats_2_file"]])
  }
  saveRDS(coloc_out, paste0(output_folder, "/", study, "_unfilt.RDS"))
  message(paste0("Unfiltered file written: ", output_folder, "/", study, "_unfilt.RDS"))
  writexl::write_xlsx(coloc_out, paste0(output_folder, "/", study, "_unfilt.xlsx"))
  message(paste0("Unfiltered file written: ", output_folder, "/", study, "_unfilt.xlsx"))
  coloc_out <- coloc_out_filter(coloc_out,
                                PP.H4.abf_filt=PP.H4.abf_filt,
                                PP.H3.abf_filt=PP.H3.abf_filt)
  writexl::write_xlsx(coloc_out, paste0(output_folder, "/", study, "_filt.xlsx"))
  message(paste0("Filtered file written: ", output_folder, "/", study, "_filt.xlsx"))
}

#' coloc_out_summary
#' @importFrom data.table rbindlist
coloc_out_summary <- function(coloc_out_all, output_folder,
                              remove_dirname = T,
                              PP.H4.abf_filt=0.5,
                              PP.H3.abf_filt=NULL) {
  coloc_out_all <- data.table::rbindlist(coloc_out_all, fill=TRUE, idcol = "Dataset")
  coloc_out_all <- coloc_out_filter(coloc_out_all,
                                    PP.H4.abf_filt=PP.H4.abf_filt,
                                    PP.H3.abf_filt=PP.H3.abf_filt)
  writexl::write_xlsx(coloc_out_all, paste0(output_folder, "/summary.xlsx"))
  message(paste0("Summary file written: ", output_folder, "/summary.xlsx"))
}

# helpers
#' coloc_out_filter
coloc_out_filter <- function(coloc_out, PP.H4.abf_filt=0.5, PP.H3.abf_filt=NULL) {
  coloc_out <- subset(coloc_out, !is.na(PP.H4.abf))
  coloc_out <- subset(coloc_out, PP.H4.abf >= PP.H4.abf_filt)
  if (!is.null(PP.H3.abf_filt)) coloc_out <- subset(coloc_out, PP.H3.abf >= PP.H3.abf_filt)
  coloc_out
}
