# Convergence view: region-centric multi-omics category overview
# Groups colocalizations by high-level category (diseases, PheWAS, imaging,
# proteins, transcripts, metabolites) for drill-down from the Convergence tab.

# =============================================================================
# Category mapping
# =============================================================================
# Group source_study values into 6 high-level categories for the overview cards.
# Order defines card display order.
TRAIT_CATEGORIES <- list(
  kidney_disease = list(
    label = "Kidney diseases & traits",
    icon = "\U0001F9EC",
    color = "#F39B7F",
    studies = c("CKDGen_r4", "CKDGen_r5")
  ),
  phewas = list(
    label = "PheWAS (diseases & traits)",
    icon = "\U0001F4CA",
    color = "#2E5C8A",
    studies = c("MVP_R4", "MVP_R4_AFR", "MVP_R4_AMR", "MVP_R4_EAS",
                "MVP_R4_EUR", "MVP_R4_META", "FinnGen_r9", "UKB_TOPMed",
                "PanUKB")
  ),
  imaging = list(
    label = "Kidney imaging",
    icon = "\U0001F9F2",
    color = "#7570B3",
    studies = c("UKB_kidney_vol")
  ),
  proteins = list(
    label = "Proteins (pQTL)",
    icon = "\U0001F9EC",
    color = "#66B266",
    studies = c("UKB_PPP_EUR", "UKB_PPP_AFR", "Icelanders_pGWAS",
                "GCKD_pGWAS")
  ),
  transcripts = list(
    label = "Transcripts (eQTL)",
    icon = "\U0001F9EC",
    color = "#9B59B6",
    studies = c("eQTLGen", "Kidney_eQTL", "GTEXv8_eQTL")
  ),
  metabolites = list(
    label = "Metabolites",
    icon = "\U0001F9EA",
    color = "#E67E22",
    studies = c("GCKD_mGWAS_plasma", "GCKD_mGWAS_urine", "GCKD_uMet")
  )
)

# Given a source_study, return its category id or "other"
category_for_study <- function(study) {
  for (cat_id in names(TRAIT_CATEGORIES)) {
    if (study %in% TRAIT_CATEGORIES[[cat_id]]$studies) return(cat_id)
  }
  "other"
}

# Count traits per category in a filtered data.table
# Returns a named integer vector with one entry per category
count_by_category <- function(dt) {
  zero <- setNames(rep(0L, length(TRAIT_CATEGORIES)), names(TRAIT_CATEGORIES))
  if (is.null(dt) || nrow(dt) == 0) return(zero)
  cats <- vapply(dt$source_study, category_for_study, character(1))
  counts <- table(factor(cats, levels = names(TRAIT_CATEGORIES)))
  setNames(as.integer(counts), names(counts))
}
