# Trait name extraction from study-specific columns
#
# Maps study-specific column names to human-readable trait names.
# This is the single place to add support for new studies.

# Column mapping: study -> column name containing the trait name
# Add new studies here
trait_column_map <- list(
  "eQTLGen"           = "eQTLGen_gene_name",
  "Kidney_eQTL"       = "Kidney_eQTL_gene_name",
  "GTEXv8_eQTL"       = "GTEXv8_eQTL_gene_name",
  "UKB_PPP_EUR"       = "UKB_PPP_EUR_olink_target_fullname",
  "Icelanders_pGWAS"  = "Icelanders_pGWAS_Protein..short.name.",
  "GCKD_mGWAS_plasma" = "GCKD_mGWAS_plasma_BIOCHEMICAL",
  "GCKD_mGWAS_urine"  = "GCKD_mGWAS_urine_BIOCHEMICAL",
  "FinnGen_r9"        = "FinnGen_r9_phenotype",
  "UKB_TOPMed"        = "UKB_TOPMed_phenostring",
  "MVP_R4"            = "MVP_R4_Title.of.analysis",
  "MVP_R4_EUR"        = "MVP_R4_EUR_Title.of.analysis",
  "MVP_R4_AFR"        = "MVP_R4_AFR_Title.of.analysis",
  "MVP_R4_AMR"        = "MVP_R4_AMR_Title.of.analysis",
  "MVP_R4_EAS"        = "MVP_R4_EAS_Title.of.analysis",
  "MVP_R4_META"       = "MVP_R4_META_Title.of.analysis",
  "CKDGen_r4"         = "CKDGen_r4_Name",
  "CKDGen_r5"         = "CKDGen_r5_ckdgen_r5_name",
  "pho_ca"            = "pho_ca_pho_ca_name"
)

# Display name column mapping (for MVP: use Analyzed.variable for human-readable names)
display_column_map <- list(
  "MVP_R4"     = "MVP_R4_Analyzed.variable",
  "MVP_R4_EUR" = "MVP_R4_EUR_Analyzed.variable",
  "MVP_R4_AFR" = "MVP_R4_AFR_Analyzed.variable",
  "MVP_R4_AMR" = "MVP_R4_AMR_Analyzed.variable",
  "MVP_R4_EAS" = "MVP_R4_EAS_Analyzed.variable",
  "MVP_R4_META"= "MVP_R4_META_Analyzed.variable"
)

#' Extract trait name from a row (for file matching / network node IDs)
#'
#' @param row A single-row data.table/data.frame
#' @return Character: trait name (sanitized for filenames)
get_trait_name <- function(row) {
  source <- row$source_study[1]

  # Special case: eQTLGen appends ENSG ID
  if (source == "eQTLGen") {
    col <- "eQTLGen_gene_name"
    if (col %in% names(row) && !is.na(row[[col]][1])) {
      if ("sumstats_2_file" %in% names(row)) {
        ensg <- gsub(".*_(ENSG[0-9]+).*", "\\1", basename(row$sumstats_2_file[1]))
        return(sanitize_trait_name(paste0(row[[col]][1], "_", ensg)))
      }
      return(sanitize_trait_name(row[[col]][1]))
    }
    return("Unknown")
  }

  # Special case: UKB_kidney_vol - extract from filename
  if (source == "UKB_kidney_vol") {
    if ("UKB_kidney_vol_phenotype" %in% names(row) && !is.na(row$UKB_kidney_vol_phenotype[1])) {
      return(sanitize_trait_name(row$UKB_kidney_vol_phenotype[1]))
    }
    return(sanitize_trait_name(gsub("model1_qnorm_(.*)_chr.*", "\\1", basename(row$sumstats_2_file[1]))))
  }

  # Generic: look up column from map
  col <- trait_column_map[[source]]
  if (!is.null(col) && col %in% names(row) && !is.na(row[[col]][1])) {
    return(sanitize_trait_name(row[[col]][1]))
  }

  # Fallback: use filename
  if ("sumstats_2_file" %in% names(row)) {
    return(sanitize_trait_name(basename(row$sumstats_2_file[1])))
  }

  "Unknown"
}

#' Get human-readable display name for a trait
#'
#' @param row A single-row data.table/data.frame
#' @return Character: display name (human-readable, may contain spaces)
get_trait_display_name <- function(row) {
  source <- row$source_study[1]

  # Check display column map (e.g., MVP full phecode name)
  col <- display_column_map[[source]]
  if (!is.null(col) && col %in% names(row) && !is.na(row[[col]][1])) {
    name <- row[[col]][1]
    name <- gsub("\\.", " ", name)
    name <- gsub("  +", " ", trimws(name))
    return(name)
  }

  # Fall back to trait_name (but without filename sanitization)
  col <- trait_column_map[[source]]
  if (!is.null(col) && col %in% names(row) && !is.na(row[[col]][1])) {
    return(as.character(row[[col]][1]))
  }

  # UKB_kidney_vol
  if (source == "UKB_kidney_vol") {
    return(gsub("model1_qnorm_(.*)_chr.*", "\\1", basename(row$sumstats_2_file[1])))
  }

  get_trait_name(row)
}

#' Sanitize trait name for filename safety
#' @keywords internal
sanitize_trait_name <- function(name) {
  if (is.null(name) || is.na(name) || name == "") return("Unknown")
  name <- gsub("[^A-Za-z0-9_.-]", "_", name)
  name <- gsub("_+", "_", name)
  name <- gsub("^_|_$", "", name)
  name
}

#' Get the trait column name for a given study
#' Used in Trait View to populate dropdown choices
#'
#' @param study Source study name
#' @return Column name string, or NULL
get_trait_column <- function(study) {
  # For display purposes, prefer display_column_map
  col <- display_column_map[[study]]
  if (!is.null(col)) return(col)
  trait_column_map[[study]]
}
