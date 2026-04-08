# Bundle key helpers - pure, no Shiny reactive dependencies.
#
# Bundle keys uniquely identify a per-ancestry trait inside a region RDS
# bundle (e.g. "MVP_R4_EUR__MVP_R4.1000G_AGR.A1C_Max_INT.EUR.GIA.dbGaP.txt.gz").
# Consensus keys collapse the per-ancestry variants of the same trait
# into one stable lookup (e.g. "MVP_R4__MVP_R4.1000G_AGR.A1C_Max_INT.GIA.dbGaP.txt.gz").
#
# These functions never touch Shiny inputs or reactives - callers pass
# the virtual-study flag explicitly.

# Ancestry codes recognized by the virtual multi-ancestry codepath.
# New codes (e.g. SAS) can be added here and in ANCESTRY_COLORS in
# config.R.
ANC_CODES <- c("AFR", "AMR", "EAS", "EUR", "META")

# Build the region-bundle lookup key for a trait row. For a virtual
# multi-ancestry study the region bundle is per-ancestry, so the key
# must be reassembled as "<source_study>_<ancestry>__<basename>" even
# though coloc_data() unifies source_study across ancestries (e.g.
# "MVP_R4").
#
# Vector-safe: accepts scalar or data.table-column inputs and returns
# a character vector of the same length.
#
# @param source_study Unified source study column (character scalar or vector)
# @param sumstats_2_file Full file path to the sumstats_2 file (character)
# @param ancestry Ancestry code per row (character, NA for non-virtual)
# @param is_virtual TRUE when the current session is rendering a virtual
#   multi-ancestry study. Only then does the ancestry suffix get inserted.
make_bundle_key <- function(source_study, sumstats_2_file,
                            ancestry = NA_character_, is_virtual = FALSE) {
  base_key <- paste0(source_study, "__", basename(sumstats_2_file))
  if (!isTRUE(is_virtual) || is.null(ancestry)) {
    return(base_key)
  }
  anc_key <- paste0(source_study, "_", ancestry, "__", basename(sumstats_2_file))
  use_anc <- !is.na(ancestry) & nzchar(ancestry)
  ifelse(use_anc, anc_key, base_key)
}

# Consensus trait key: collapses all per-ancestry variants of the same
# trait into one stable key by stripping the "_<ANC>" study suffix and
# the ".<ANC>." filename segment.
make_consensus_trait_key <- function(bundle_key) {
  anc_alt <- paste(ANC_CODES, collapse = "|")
  k <- sub(paste0("_(", anc_alt, ")__"), "__", bundle_key)
  sub(paste0("\\.(", anc_alt, ")\\."), ".", k)
}

# Given a consensus trait key and a per-ancestry region bundle, try to
# locate the matching ancestry-specific bundle key inside that bundle.
# Returns a single key (character) or NULL if nothing matches.
resolve_consensus_in_bundle <- function(consensus_key, bundle_keys, anc) {
  parts <- strsplit(consensus_key, "__", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NULL)
  prefix <- parts[1]; fn <- parts[2]

  # Fast path: rebuild the expected per-ancestry candidate keys.
  cand <- c(
    paste0(prefix, "_", anc, "__",
           sub("\\.GIA", paste0(".", anc, ".GIA"), fn)),
    paste0(prefix, "_", anc, "__",
           sub("^(MVP_R4\\.[^.]+\\.[^.]+)\\.",
               paste0("\\1.", anc, "."), fn))
  )
  hit <- intersect(cand, bundle_keys)
  if (length(hit) > 0) return(hit[1])

  # Fallback: any bundle key that starts with "<prefix>_<anc>__" and
  # reduces (via make_consensus_trait_key) to the same consensus key.
  pref <- paste0(prefix, "_", anc, "__")
  cands <- bundle_keys[startsWith(bundle_keys, pref)]
  if (length(cands) == 0) return(NULL)
  matches <- cands[vapply(cands,
                          function(k) identical(make_consensus_trait_key(k),
                                                consensus_key),
                          logical(1))]
  if (length(matches) > 0) matches[1] else NULL
}
