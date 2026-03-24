# Configuration file for GWAS Colocalization Network Viewer

# Define color palette for studies
study_colors <- c(
  # Phenotypes
  "CKDGen_r4" = "#F39B7F",
  "FinnGen_r9" = "#7BA3D9",
  "MVP_R4" = "#2E5C8A",
  "MVP_R4_EUR" = "#2E5C8A",
  "UKB_TOPMed" = "#5D9ED3",

  # pQTL
  "UKB_PPP_EUR" = "#66B266",
  "Icelanders_pGWAS" = "#4D9966",

  # eQTL
  "eQTLGen" = "#9B59B6",
  "Kidney_eQTL" = "#B97FC9",

  # mQTL
  "GCKD_mGWAS_plasma" = "#E67E22",
  "GCKD_mGWAS_urine" = "#D35400",

  # Imaging
  "UKB_kidney_vol" = "#7570B3"
)

# Study categories for grouped display
study_categories <- list(
  "Phenotypes" = c("CKDGen_r4", "FinnGen_r9", "MVP_R4", "MVP_R4_EUR", "UKB_TOPMed"),
  "pQTL" = c("UKB_PPP_EUR", "Icelanders_pGWAS"),
  "eQTL" = c("eQTLGen", "Kidney_eQTL"),
  "mQTL" = c("GCKD_mGWAS_plasma", "GCKD_mGWAS_urine"),
  "Imaging" = c("UKB_kidney_vol")
)

# Data path: use option, env var, or default
DATA_PATH <- getOption("genepicoloc.data_path",
  Sys.getenv("GENEPICOLOC_DATA", "data"))

# Auto-discover available studies from RDS files
# Supports two layouts:
#   1. Flat: data/<trait>_annot_filt.RDS (atlas layout)
#   2. Nested: data/<folder>/annot/annot_filt.RDS (CKDGen layout)
discover_studies <- function(data_path) {
  studies <- list()

  # Try flat layout first
  flat_files <- list.files(data_path, pattern = "_annot_filt\\.RDS$", full.names = TRUE)
  if (length(flat_files) > 0) {
    for (f in flat_files) {
      name <- gsub("_annot_filt\\.RDS$", "", basename(f))
      studies[[name]] <- f
    }
    return(studies)
  }

  # Try nested layout
  dirs <- list.dirs(data_path, recursive = FALSE)
  for (d in dirs) {
    f <- file.path(d, "annot", "annot_filt.RDS")
    if (file.exists(f)) {
      name <- basename(d)
      studies[[name]] <- f
    }
  }

  studies
}

DEFAULT_AVAILABLE_STUDIES <- discover_studies(DATA_PATH)

# Default study base path (for CKDGen nested layout compatibility)
DEFAULT_STUDY_BASE_PATH <- DATA_PATH

# Feature toggles
# INCLUDE_CKDGEN_R5 removed - not relevant for public package
SHOW_EXPORT_CUSTOMIZATION <- FALSE
