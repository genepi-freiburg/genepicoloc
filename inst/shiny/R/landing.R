# Landing page category cards.
#
# Pure UI builders - no Shiny reactives are read here. All the data
# this file needs is already available at module load time:
#   DEFAULT_AVAILABLE_STUDIES, DEFAULT_VIRTUAL_STUDIES (from config.R)
#   DATA_PATH                                           (from config.R)
#   umet_label                                          (from annotation_db.R)
#
# The server wires clicks to current_study() via
# Shiny.setInputValue('selected_study', ...) which is handled in app.R.

# Category definitions for the landing page.
# `section` groups cards into rows on the landing page:
#   "featured"   - CKDGen r4 (full atlas) and MVP kidney (multi-ancestry)
#   "additional" - MRI volumes and urine metabolomics (supporting data)
atlas_categories <- list(
  list(
    id = "kidney_disease",
    section = "featured",
    icon = "\U0001F9EC",
    title = "CKDGen Round 4",
    description = paste(
      "Kidney function and disease phenotypes (eGFR, BUN, UACR, urate,",
      "gout, microalbuminuria) from the CKDGen Round 4 meta-analysis,",
      "colocalized against 11 molecular and clinical datasets."
    ),
    traits = list(
      list(id = "eGFR",  label = "eGFR (creatinine)", desc = "Estimated glomerular filtration rate"),
      list(id = "BUN",   label = "Blood Urea Nitrogen", desc = "Kidney filtration marker"),
      list(id = "UACR",  label = "UACR",               desc = "Urinary albumin-to-creatinine ratio"),
      list(id = "urate", label = "Serum Urate",        desc = "Urate levels"),
      list(id = "gout",  label = "Gout",               desc = "Inflammatory arthritis"),
      list(id = "MA",    label = "Microalbuminuria",   desc = "Early kidney damage marker")
    )
  ),
  list(
    id = "mvp_kidney",
    section = "featured",
    icon = "\U0001F30D",
    title = "MVP Kidney (multi-ancestry)",
    description = paste(
      "Million Veteran Program kidney traits colocalized across five",
      "ancestries (AFR, AMR, EAS, EUR, META) against ancestry-matched",
      "MVP_R4 PheWAS. Each trait loads all available ancestries at once",
      "and overlays them in the Manhattan and regional plots."
    ),
    traits = list(),
    autopopulate_virtual = TRUE
  ),
  list(
    id = "kidney_mri",
    section = "additional",
    icon = "\U0001F9F2",
    title = "Kidney MRI Volumes",
    description = "UK Biobank kidney MRI - structural imaging phenotypes (BSA-adjusted).",
    traits = list(
      list(id = "MRI_tkv",     label = "Total Kidney Volume", desc = "BSA-adjusted"),
      list(id = "MRI_cortex",  label = "Cortex Volume",       desc = "BSA-adjusted"),
      list(id = "MRI_medulla", label = "Medulla Volume",      desc = "BSA-adjusted"),
      list(id = "MRI_hilus",   label = "Hilus Volume",        desc = "BSA-adjusted")
    )
  ),
  list(
    id = "metabolomics",
    section = "additional",
    icon = "\U0001F9EA",
    title = "Urine Metabolomics",
    description = paste(
      "GCKD urine metabolome GWAS (Schlosser et al., Nat Genet 2023) -",
      "1,409 urine metabolites measured by Metabolon, colocalized against",
      "Tier 1 datasets."
    ),
    traits = list(),
    autopopulate = "GCKD_uMet"
  )
)

# Is this trait id resolvable? Accepts real study ids (present in
# DEFAULT_AVAILABLE_STUDIES) OR virtual multi-ancestry ids.
is_known_study <- function(id) {
  id %in% names(DEFAULT_AVAILABLE_STUDIES) ||
    id %in% names(DEFAULT_VIRTUAL_STUDIES)
}

# Build the traits list for one landing card.
#
# - If `cat$autopopulate_virtual` is TRUE, emit one entry per virtual
#   study id (MVP kidney uses this).
# - If `cat$autopopulate = "<atlas_category>"`, expand to all studies
#   auto-discovered under that atlas folder, labeled via the bundled
#   annotation DB. Used for urine metabolomics where hand-listing
#   1,409 metabolites is impractical.
# - Otherwise return cat$traits as-is.
expand_card_traits <- function(cat) {
  if (length(cat$traits) > 0) return(cat$traits)

  if (isTRUE(cat$autopopulate_virtual)) {
    if (length(DEFAULT_VIRTUAL_STUDIES) == 0) return(list())
    return(lapply(names(DEFAULT_VIRTUAL_STUDIES), function(vid) {
      v <- DEFAULT_VIRTUAL_STUDIES[[vid]]
      lbl <- sub("^MVPkid_", "", vid)
      list(id = vid, label = lbl,
           desc = paste0(length(v$ancestries), " ancestries: ",
                         paste(v$ancestries, collapse = ", ")))
    }))
  }

  if (is.null(cat$autopopulate)) return(cat$traits)

  cats_attr <- attr(DEFAULT_AVAILABLE_STUDIES, "categories")
  if (is.null(cats_attr)) return(list())
  ids <- names(cats_attr)[unlist(cats_attr) == cat$autopopulate]

  # Drop trait_ids whose coloc RDS has zero rows. Metabolites with no
  # high-confidence colocs would render as empty Region Views, so we
  # hide them from the landing card entirely.
  ids <- Filter(function(id) {
    path <- DEFAULT_AVAILABLE_STUDIES[[id]]
    if (is.null(path) || !file.exists(path)) return(FALSE)
    tryCatch(nrow(readRDS(path)) > 0, error = function(e) FALSE)
  }, ids)

  label_fn <- if (cat$autopopulate == "GCKD_uMet") umet_label else identity
  lapply(ids, function(id) {
    list(id = id, label = label_fn(id), desc = id)
  })
}

# Render one category card (used by the landing renderUI in app.R).
build_card <- function(cat) {
  cat$traits <- expand_card_traits(cat)
  available <- sapply(cat$traits, function(t) {
    !isTRUE(t$disabled) && is_known_study(t$id)
  })
  n_available <- sum(available)
  stats_text <- if (n_available > 0) {
    paste0(n_available, " traits available")
  } else {
    "Coming soon"
  }

  shiny::div(class = "category-card",
    shiny::div(class = "card-icon", cat$icon),
    shiny::h3(cat$title),
    shiny::p(class = "card-description", cat$description),
    shiny::p(class = "card-stats", stats_text),
    shiny::div(class = "trait-list",
      lapply(cat$traits, function(trait) {
        is_virtual <- trait$id %in% names(DEFAULT_VIRTUAL_STUDIES)
        is_disabled <- isTRUE(trait$disabled) || !is_known_study(trait$id)
        cats <- attr(DEFAULT_AVAILABLE_STUDIES, "categories")
        cat_name <- if (!is.null(cats)) cats[[trait$id]] else NULL
        has_regional <- if (is_virtual) {
          vdirs <- DEFAULT_VIRTUAL_STUDIES[[trait$id]]$regional_dirs
          any(!is.na(unlist(vdirs)))
        } else {
          !is_disabled && (
            (!is.null(cat_name) && dir.exists(file.path(DATA_PATH, cat_name, "regional", trait$id))) ||
            dir.exists(file.path(DATA_PATH, "regional", trait$id))
          )
        }
        css_class <- paste("trait-btn",
          if (is_disabled) "disabled",
          if (has_regional) "has-regional"
        )
        if (is_disabled) {
          shiny::tags$span(class = css_class, title = trait$desc, trait$label)
        } else {
          shiny::tags$span(
            class = css_class,
            title = trait$desc,
            onclick = paste0("Shiny.setInputValue('selected_study', '",
                             trait$id, "', {priority: 'event'});"),
            trait$label
          )
        }
      })
    )
  )
}
