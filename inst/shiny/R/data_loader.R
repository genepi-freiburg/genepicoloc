# Data loading and trait name extraction

#' Load colocalization data for a study
#'
#' Reads annot_filt.RDS and normalizes trait names into a single trait_name column.
#'
#' @param config App configuration from load_app_config()
#' @param study_name Study name (e.g., "eGFR")
#' @return data.table with normalized trait_name column
load_coloc_data <- function(config, study_name) {

  f <- config$studies[[study_name]]$file
  if (is.null(f) || !file.exists(f)) {
    warning("File not found for study: ", study_name)
    return(data.table())
  }

  dt <- as.data.table(readRDS(f))
  if (nrow(dt) == 0) return(dt)

  # Normalize trait names from study-specific columns into a single column
  dt[, trait_name := extract_trait_name(.SD), by = seq_len(nrow(dt))]

  # Create trait_id for node identification
  dt[, trait_id := paste0(source_study, ":", trait_name)]

  # Ensure region ID
  dt[, region_id := paste0("chr", CHR_var, ":", BP_START_var, "-", BP_STOP_var)]

  dt
}


#' Extract human-readable trait name from study-specific columns
#'
#' Handles all known study column naming conventions.
#'
#' @param row A single-row data.table
#' @return Character: trait name
extract_trait_name <- function(row) {
  study <- row$source_study[1]

  name <- switch(study,
    "MVP_R4_EUR" = , "MVP_R4" = , "MVP_R4_AFR" = , "MVP_R4_AMR" = ,
    "MVP_R4_EAS" = , "MVP_R4_META" = {
      col <- intersect(c("MVP_R4_EUR_Analyzed.variable", "MVP_R4_Analyzed.variable"), names(row))
      if (length(col) > 0) gsub("\\.", " ", row[[col[1]]]) else basename(row$sumstats_2_file)
    },

    "FinnGen_r9" = {
      if ("FinnGen_r9_phenotype" %in% names(row)) row$FinnGen_r9_phenotype
      else basename(row$sumstats_2_file)
    },

    "UKB_PPP_EUR" = {
      if ("UKB_PPP_EUR_olink_target_fullname" %in% names(row)) row$UKB_PPP_EUR_olink_target_fullname
      else basename(row$sumstats_2_file)
    },

    "Icelanders_pGWAS" = {
      if ("Icelanders_pGWAS_Protein..short.name." %in% names(row)) row$`Icelanders_pGWAS_Protein..short.name.`
      else basename(row$sumstats_2_file)
    },

    "eQTLGen" = {
      if ("eQTLGen_gene_name" %in% names(row)) row$eQTLGen_gene_name
      else basename(row$sumstats_2_file)
    },

    "Kidney_eQTL" = {
      if ("Kidney_eQTL_gene_name" %in% names(row)) row$Kidney_eQTL_gene_name
      else basename(row$sumstats_2_file)
    },

    "CKDGen_r4" = {
      if ("CKDGen_r4_Name" %in% names(row)) row$CKDGen_r4_Name
      else basename(row$sumstats_2_file)
    },

    "CKDGen_r5" = {
      if ("CKDGen_r5_ckdgen_r5_name" %in% names(row)) row$CKDGen_r5_ckdgen_r5_name
      else basename(row$sumstats_2_file)
    },

    "GCKD_mGWAS_plasma" = {
      if ("GCKD_mGWAS_plasma_BIOCHEMICAL" %in% names(row)) row$GCKD_mGWAS_plasma_BIOCHEMICAL
      else basename(row$sumstats_2_file)
    },

    "GCKD_mGWAS_urine" = {
      if ("GCKD_mGWAS_urine_BIOCHEMICAL" %in% names(row)) row$GCKD_mGWAS_urine_BIOCHEMICAL
      else basename(row$sumstats_2_file)
    },

    "UKB_TOPMed" = {
      if ("UKB_TOPMed_phenostring" %in% names(row)) row$UKB_TOPMed_phenostring
      else basename(row$sumstats_2_file)
    },

    "UKB_kidney_vol" = {
      if ("UKB_kidney_vol_phenotype" %in% names(row)) row$UKB_kidney_vol_phenotype
      else gsub("model1_qnorm_(.*)_chr.*", "\\1", basename(row$sumstats_2_file))
    },

    # Default: use filename
    basename(row$sumstats_2_file)
  )

  if (is.null(name) || is.na(name) || name == "") basename(row$sumstats_2_file)
  else as.character(name)
}
