# Study Metadata and Reference Information
# Contains display names, descriptions, PMIDs, and ancestry information

# Display names for studies (user-friendly)
study_display_names <- c(
  "GTEXv8_eQTL" = "GTEx V8 (eQTL)",
  "Kidney_eQTL" = "Kidney eQTL",
  "eQTLGen" = "eQTLGen Phase I",
  "UKB_PPP_EUR" = "UK Biobank Proteomics (PPP)",
  "Icelanders_pGWAS" = "Icelanders pGWAS",
  "GCKD_mGWAS_plasma" = "GCKD mGWAS (Plasma)",
  "GCKD_mGWAS_urine" = "GCKD mGWAS (Urine)",
  "UKB_TOPMed" = "UK Biobank TOPMed",
  "FinnGen_r9" = "FinnGen R9",
  "MVP_R4" = "MVP R4 (Million Veteran)",
  "CKDGen_r5" = "CKDGen R5",
  "pho_ca" = "Phosphate & Calcium"
)

# Detailed study information
study_info <- list(
  "GTEXv8_eQTL" = list(
    name = "GTEx V8 (eQTL)",
    description = "GTEx V8 - 49 tissues, 39,867 unique transcripts",
    url = "https://www.gtexportal.org/home/",
    pmid = NULL,
    ancestry = "Multi-ancestry",
    notes = "Contains chrX"
  ),

  "Kidney_eQTL" = list(
    name = "Kidney eQTL",
    description = "Kidney eQTL - 5 summary statistics, 14,694 unique transcripts",
    url = "https://susztaklab.com/Kidney_eQTL/",
    pmid = NULL,
    ancestry = "European",
    notes = "Does not have chrX"
  ),

  "eQTLGen" = list(
    name = "eQTLGen Phase I",
    description = "eQTLGen Phase I - 19,250 unique transcripts",
    url = "https://eqtlgen.org/phase1.html",
    pmid = NULL,
    ancestry = "European",
    notes = "Contains chrX"
  ),

  "UKB_PPP_EUR" = list(
    name = "UK Biobank Proteomics (PPP)",
    description = "UK Biobank Proteomics - 2,940 unique protein analytes",
    url = "https://www.nature.com/articles/s41586-023-06592-6",
    pmid = "37794186",
    ancestry = "European",
    notes = "Contains chrX"
  ),

  "Icelanders_pGWAS" = list(
    name = "Icelanders pGWAS",
    description = "Icelanders pGWAS - 4,909 unique protein analytes",
    url = "https://pubmed.ncbi.nlm.nih.gov/34857953/",
    pmid = "34857953",
    ancestry = "European (Icelandic)",
    notes = "Contains chrX"
  ),

  "GCKD_mGWAS_plasma" = list(
    name = "GCKD mGWAS (Plasma)",
    description = "GCKD metabolomics - >1,000 metabolites in plasma",
    url = "https://www.nature.com/articles/s41588-023-01409-8",
    pmid = "37231185",
    ancestry = "European",
    notes = "German Chronic Kidney Disease study"
  ),

  "GCKD_mGWAS_urine" = list(
    name = "GCKD mGWAS (Urine)",
    description = "GCKD metabolomics - >1,000 metabolites in urine",
    url = "https://www.nature.com/articles/s41588-023-01409-8",
    pmid = "37231185",
    ancestry = "European",
    notes = "German Chronic Kidney Disease study"
  ),

  "UKB_TOPMed" = list(
    name = "UK Biobank TOPMed",
    description = "UK Biobank TOPMed - 1,419 clinical outcomes (diseases)",
    url = "https://pheweb.org/UKB-TOPMed/about",
    pmid = NULL,
    ancestry = "Multi-ancestry",
    notes = "UK Biobank with TOPMed imputation"
  ),

  "FinnGen_r9" = list(
    name = "FinnGen R9",
    description = "FinnGen Release 9 - 2,272 clinical outcomes (diseases)",
    url = "https://r9.finngen.fi/",
    pmid = NULL,
    ancestry = "European (Finnish)",
    notes = "Finnish biobank"
  ),

  "MVP_R4" = list(
    name = "MVP R4 (Million Veteran)",
    description = "Million Veteran Program - 2,068 quantitative and case-control traits",
    url = "https://www.science.org/doi/10.1126/science.adj1182",
    pmid = "37733863",
    ancestry = "Multi-ancestry",
    notes = "5,000 traits (PheWAS) from VA Million Veteran Program"
  ),

  "CKDGen_r5" = list(
    name = "CKDGen R5",
    description = "CKDGen Consortium - Kidney function and disease traits",
    url = "https://ckdgen.imbi.uni-freiburg.de/",
    pmid = "Multiple (see website)",
    ancestry = "Multi-ancestry",
    notes = "UACR, MA, gout, urate, eGFR, BUN and other kidney traits"
  ),

  "pho_ca" = list(
    name = "Phosphate & Calcium",
    description = "Mineral metabolism traits",
    url = NULL,
    pmid = NULL,
    ancestry = "Multi-ancestry",
    notes = "CKDGen mineral metabolism traits"
  )
)

# Input sumstats metadata (for Home page study cards)
# Display names for input GWAS traits
input_sumstats_display_names <- c(
  "eGFR" = "eGFR (creatinine)",
  "BUN" = "Blood Urea Nitrogen",
  "UACR" = "Urinary Albumin-to-Creatinine Ratio",
  "urate" = "Serum Urate",
  "gout" = "Gout",
  "MA" = "Microalbuminuria"
)

input_sumstats_descriptions <- c(
  "eGFR" = "Estimated glomerular filtration rate - primary kidney function marker (Wuttke et al. 2019)",
  "BUN" = "Blood urea nitrogen - kidney filtration marker (Wuttke et al. 2019)",
  "UACR" = "Urinary albumin-to-creatinine ratio - kidney damage marker (Teumer et al. 2019)",
  "urate" = "Serum urate levels - associated with gout and CKD (Tin et al. 2019)",
  "gout" = "Gout - inflammatory arthritis caused by urate crystal deposition (Tin et al. 2019)",
  "MA" = "Microalbuminuria - early marker of kidney damage (Teumer et al. 2019)"
)

# Build input_sumstats_info dynamically from discovered studies
# Regional plot availability is auto-detected from the data directory
input_sumstats_info <- lapply(names(DEFAULT_AVAILABLE_STUDIES), function(study_name) {
  # Check category layout first, then clean, then legacy
  cats <- attr(DEFAULT_AVAILABLE_STUDIES, "categories")
  cat_name <- if (!is.null(cats)) cats[[study_name]] else NULL
  regional_dir <- NULL
  if (!is.null(cat_name)) {
    candidate <- file.path(DATA_PATH, cat_name, "regional", study_name)
    if (dir.exists(candidate)) regional_dir <- candidate
  }
  if (is.null(regional_dir)) {
    candidate <- file.path(DATA_PATH, "regional", study_name)
    if (dir.exists(candidate)) regional_dir <- candidate
  }
  if (is.null(regional_dir)) {
    candidate <- file.path(DATA_PATH, "regional_plots", study_name)
    if (dir.exists(candidate)) regional_dir <- candidate
  }
  has_regional <- !is.null(regional_dir) && length(list.files(regional_dir)) > 0
  list(
    display_name = if (study_name %in% names(input_sumstats_display_names))
      input_sumstats_display_names[[study_name]] else study_name,
    ancestry = "Multi-ancestry",
    description = if (study_name %in% names(input_sumstats_descriptions))
      input_sumstats_descriptions[[study_name]] else NULL,
    regional_plot = if (has_regional) "Available" else "Not available"
  )
})
names(input_sumstats_info) <- names(DEFAULT_AVAILABLE_STUDIES)

# Helper function to get display name
get_study_display_name <- function(study_code) {
  if (study_code %in% names(study_display_names)) {
    return(study_display_names[study_code])
  }
  return(study_code)  # Return original if not found
}

# Helper function to format study info for modal
format_study_info <- function(study_code) {
  if (!study_code %in% names(study_info)) {
    return(paste("No information available for", study_code))
  }

  info <- study_info[[study_code]]

  html_text <- paste0(
    "<h4>", info$name, "</h4>",
    "<p><strong>Description:</strong> ", info$description, "</p>",
    "<p><strong>Ancestry:</strong> ", info$ancestry, "</p>"
  )

  if (!is.null(info$url)) {
    html_text <- paste0(html_text,
                       "<p><strong>Website:</strong> <a href='", info$url,
                       "' target='_blank'>", info$url, "</a></p>")
  }

  if (!is.null(info$pmid) && info$pmid != "Multiple (see website)") {
    html_text <- paste0(html_text,
                       "<p><strong>PMID:</strong> <a href='https://pubmed.ncbi.nlm.nih.gov/",
                       info$pmid, "/' target='_blank'>", info$pmid, "</a></p>")
  } else if (!is.null(info$pmid)) {
    html_text <- paste0(html_text, "<p><strong>PMID:</strong> ", info$pmid, "</p>")
  }

  if (!is.null(info$notes)) {
    html_text <- paste0(html_text, "<p><em>Note: ", info$notes, "</em></p>")
  }

  return(HTML(html_text))
}
