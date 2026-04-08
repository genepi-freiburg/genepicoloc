# Loader and helpers for the bundled Shiny annotation database.
#
# Single source of truth for trait/study labels. Built from the internal
# annotation_db.RDS via inst/scripts/build_shiny_annotation_db.R.

# Resolve the bundled file. The canonical location is
# inst/shiny/extdata/shiny_annotation_db.RDS so both the Dockerfile COPY
# and the deploy.sh dev-mode bind-mount (both rooted at inst/shiny/) can
# reach it.
.find_annotation_db <- function() {
  # Container / production: Dockerfile copies extdata/ to /app/extdata/
  if (file.exists("/app/extdata/shiny_annotation_db.RDS")) {
    return("/app/extdata/shiny_annotation_db.RDS")
  }
  # Dev mode (running from inst/shiny/)
  dev_candidates <- c(
    "extdata/shiny_annotation_db.RDS",       # cwd = inst/shiny
    "../extdata/shiny_annotation_db.RDS",    # cwd = inst/shiny/R
    "inst/shiny/extdata/shiny_annotation_db.RDS"  # cwd = repo root
  )
  for (p in dev_candidates) {
    if (file.exists(p)) return(normalizePath(p))
  }
  NULL
}

ANNOTATION_DB <- local({
  path <- .find_annotation_db()
  if (is.null(path)) {
    warning("shiny_annotation_db.RDS not found - trait labels will fall back ",
            "to raw IDs. Run inst/scripts/build_shiny_annotation_db.R.")
    return(list())
  }
  readRDS(path)
})

# Extract the GCKD urine metabolite key (M4_u_NNNN) from a uMet trait_id
# such as "uMet_gckd_EUR_TopMed_2023-05-19_M4_u_1107_504_M4_u_1107.regenie".
umet_key <- function(trait_id) {
  m <- regmatches(trait_id, regexpr("M4_u_[0-9]+", trait_id))
  if (length(m) == 0 || identical(m, character(0))) return(NA_character_)
  m
}

# Friendly label for a uMet trait_id (e.g. "allantoin (M4_u_1107)").
# Falls back to the raw key if the metabolite isn't in the bundled DB.
umet_label <- function(trait_id) {
  key <- umet_key(trait_id)
  if (is.na(key)) return(trait_id)
  ann <- ANNOTATION_DB[["GCKD_mGWAS_urine"]]
  if (is.null(ann)) return(key)
  hit <- ann$BIOCHEMICAL[ann$merge_column == key]
  if (length(hit) == 0 || is.na(hit[1])) return(key)
  paste0(hit[1], " (", key, ")")
}

# Generic friendly label dispatcher for any trait_id. Today only uMet has a
# rewrite; everything else returns the input unchanged. Add new study types
# here as needed.
trait_label <- function(trait_id) {
  if (length(trait_id) != 1 || is.na(trait_id) || !nzchar(trait_id)) {
    return(trait_id)
  }
  if (grepl("^uMet_", trait_id)) return(umet_label(trait_id))
  trait_id
}
