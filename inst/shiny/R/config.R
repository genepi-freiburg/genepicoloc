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
# In container, /app/data is the mount point; locally, use "data" relative to app dir
DATA_PATH <- getOption("genepicoloc.data_path",
  Sys.getenv("GENEPICOLOC_DATA",
    if (dir.exists("/app/data")) "/app/data" else "data"))

# Unified study registry
# -------------------------------------------------------------------------
# Every study gets the same structure regardless of ancestry count:
#   study_id -> list(
#     category      = "CKDGen_r4",
#     ancestries    = c("EUR"),           # always a character vector
#     coloc_files   = list(EUR = path),   # named by ancestry
#     regional_dirs = list(EUR = path),   # named by ancestry (NA if missing)
#     real_ids      = c("eGFR")           # filesystem trait ids that were merged
#   )
#
# Multi-ancestry studies (e.g. MVPkid_BUN_BSP_Mean_INT_{AFR,EUR,...}) are
# collapsed into one entry with multiple ancestries. Single-ancestry studies
# get ancestries = c("EUR") (or the detected ancestry if suffix is present).
#
# Supports layouts (checked in order):
#   1. Category: data/<category>/coloc/<trait>.RDS
#   2. Flat:     data/coloc/<trait>.RDS
#   3. Legacy:   data/<trait>_annot_filt.RDS
#   4. Nested:   data/<folder>/annot/annot_filt.RDS

discover_all_studies <- function(data_path) {
  # Step 1: collect all individual coloc files with their categories
  raw_studies <- list()   # name -> file path
  raw_categories <- list()  # name -> category

  subdirs <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)

  # Multi-category: data/<category>/coloc/*.RDS
  for (d in subdirs) {
    coloc_dir <- file.path(d, "coloc")
    if (dir.exists(coloc_dir)) {
      cat_name <- basename(d)
      for (f in list.files(coloc_dir, pattern = "\\.RDS$", full.names = TRUE)) {
        name <- tools::file_path_sans_ext(basename(f))
        raw_studies[[name]] <- f
        raw_categories[[name]] <- cat_name
      }
    }
  }

  # Single-category flat: data/coloc/*.RDS
  if (length(raw_studies) == 0) {
    coloc_dir <- file.path(data_path, "coloc")
    if (dir.exists(coloc_dir)) {
      for (f in list.files(coloc_dir, pattern = "\\.RDS$", full.names = TRUE)) {
        name <- tools::file_path_sans_ext(basename(f))
        raw_studies[[name]] <- f
      }
    }
  }

  # Legacy flat: data/<trait>_annot_filt.RDS
  if (length(raw_studies) == 0) {
    for (f in list.files(data_path, pattern = "_annot_filt\\.RDS$", full.names = TRUE)) {
      name <- gsub("_annot_filt\\.RDS$", "", basename(f))
      raw_studies[[name]] <- f
    }
  }

  # Nested: data/<folder>/annot/annot_filt.RDS
  if (length(raw_studies) == 0) {
    for (d in subdirs) {
      f <- file.path(d, "annot", "annot_filt.RDS")
      if (file.exists(f)) {
        name <- basename(d)
        raw_studies[[name]] <- f
      }
    }
  }

  if (length(raw_studies) == 0) return(list())

  # Step 2: group by ancestry suffix
  ids <- names(raw_studies)
  anc_pat <- "_(AFR|AMR|CSA|EAS|EUR|FIN|META|MID|NFE|SAS)$"
  stems <- sub(anc_pat, "", ids)
  has_anc <- stems != ids

  # Find multi-ancestry groups (2+ ancestries sharing a stem + category)
  multi_stems <- character(0)
  if (any(has_anc)) {
    for (stem in unique(stems[has_anc])) {
      idx <- which(stems == stem & has_anc)
      if (length(idx) < 2) next
      # All siblings must share the same category
      cat_names <- unlist(lapply(ids[idx], function(i) raw_categories[[i]]))
      if (length(unique(cat_names)) == 1) {
        multi_stems <- c(multi_stems, stem)
      }
    }
  }

  registry <- list()

  # Step 3a: build multi-ancestry entries
  consumed <- character(0)  # track ids consumed by multi-ancestry grouping
  for (stem in multi_stems) {
    idx <- which(stems == stem & has_anc)
    matched_ids <- ids[idx]
    ancs <- sub(paste0("^", stem, "_"), "", matched_ids)
    category <- raw_categories[[matched_ids[1]]]

    coloc_files <- setNames(as.list(unlist(raw_studies[matched_ids])), ancs)
    regional_dirs <- setNames(lapply(matched_ids, function(i) {
      p <- file.path(data_path, category, "regional", i)
      if (dir.exists(p)) p else NA_character_
    }), ancs)

    registry[[stem]] <- list(
      category = category,
      ancestries = ancs,
      coloc_files = coloc_files,
      regional_dirs = regional_dirs,
      real_ids = matched_ids
    )
    consumed <- c(consumed, matched_ids)
  }

  # Step 3b: build single-ancestry entries for everything else
  remaining <- setdiff(ids, consumed)
  for (id in remaining) {
    category <- raw_categories[[id]]  # may be NULL for legacy layouts
    # Detect ancestry from suffix if present (single-ancestry with suffix)
    anc <- if (has_anc[which(ids == id)]) sub(paste0(".*_"), "", id) else "EUR"

    coloc_files <- setNames(list(raw_studies[[id]]), anc)
    regional_dirs <- setNames(list({
      if (!is.null(category)) {
        p <- file.path(data_path, category, "regional", id)
        if (dir.exists(p)) p else NA_character_
      } else NA_character_
    }), anc)

    registry[[id]] <- list(
      category = category,
      ancestries = anc,
      coloc_files = coloc_files,
      regional_dirs = regional_dirs,
      real_ids = id
    )
  }

  registry
}

DEFAULT_STUDY_REGISTRY <- discover_all_studies(DATA_PATH)

# Backward-compat shims (used during transition, will be removed)
DEFAULT_AVAILABLE_STUDIES <- {
  s <- list()
  cats <- list()
  for (nm in names(DEFAULT_STUDY_REGISTRY)) {
    entry <- DEFAULT_STUDY_REGISTRY[[nm]]
    if (length(entry$ancestries) > 1) {
      # Multi-ancestry: expose individual per-ancestry ids
      for (anc in entry$ancestries) {
        real_id <- entry$real_ids[grep(paste0("_", anc, "$"), entry$real_ids)]
        if (length(real_id) == 1) {
          s[[real_id]] <- entry$coloc_files[[anc]]
          if (!is.null(entry$category)) cats[[real_id]] <- entry$category
        }
      }
    } else {
      # Single ancestry
      real_id <- entry$real_ids[1]
      s[[real_id]] <- entry$coloc_files[[1]]
      if (!is.null(entry$category)) cats[[real_id]] <- entry$category
    }
  }
  attr(s, "categories") <- cats
  s
}

DEFAULT_VIRTUAL_STUDIES <- {
  vs <- list()
  for (nm in names(DEFAULT_STUDY_REGISTRY)) {
    entry <- DEFAULT_STUDY_REGISTRY[[nm]]
    if (length(entry$ancestries) > 1) vs[[nm]] <- entry
  }
  vs
}

# Okabe-Ito colorblind-safe palette for ancestry overlays.
ANCESTRY_COLORS <- c(
  AFR  = "#D55E00",  # vermillion
  AMR  = "#009E73",  # bluish green
  EAS  = "#F0E442",  # yellow
  EUR  = "#0072B2",  # blue
  META = "#CC79A7"   # reddish purple
)

# Gradient for the Manhattan "n ancestries" coloring (1 -> 4).
ANCESTRY_COVERAGE_COLORS <- c(
  "1" = "#fdd49e",
  "2" = "#fdbb84",
  "3" = "#e34a33",
  "4" = "#b30000"
)

# Default study base path (for CKDGen nested layout compatibility)
DEFAULT_STUDY_BASE_PATH <- DATA_PATH

# Feature toggles
# INCLUDE_CKDGEN_R5 removed - not relevant for public package
SHOW_EXPORT_CUSTOMIZATION <- FALSE
