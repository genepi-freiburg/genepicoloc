# Bundled region RDS loaders - pure data-layer code (no Shiny
# reactives). Callers from the server pass the data path and the
# virtual study metadata in directly.
#
# The slim regional RDS format is documented in docs/multi-ancestry.md
# and CLAUDE.md; briefly:
#
#   <data_path>/<region_id>.RDS  where region_id = "chr<chr>_<start>_<stop>"
#   contains list(base = data.table, traits = named list of data.tables)
#
# Caches
# ------
# Two process-wide caches speed up successive reads from the same
# region. They're plain environments (shared across Shiny sessions)
# because the underlying files are immutable and the cache keys are
# derived from canonicalized file paths. Process lifetime is fine
# because the app is bounced on deploy and data changes happen at
# deploy time.

.region_bundle_cache <- new.env(parent = emptyenv())
.multi_bundle_cache  <- new.env(parent = emptyenv())

# Load bundled region data (one RDS per region).
# Returns list(base, traits) or NULL if the RDS is missing.
load_region_bundle <- function(data_path, region_str) {
  region_id <- parse_region_id(region_str)
  if (is.null(region_id)) return(NULL)

  cache_key <- paste(data_path, region_id, sep = "||")
  cached <- .region_bundle_cache[[cache_key]]
  if (!is.null(cached)) return(cached)

  rds_file <- file.path(data_path, paste0(region_id, ".RDS"))
  if (!file.exists(rds_file)) return(NULL)

  region_data <- readRDS(rds_file)
  if (!is.null(region_data$base)) {
    region_data$base <- reconstruct_name(region_data$base)
  }
  region_data$traits <- lapply(region_data$traits, reconstruct_name)

  assign(cache_key, region_data, envir = .region_bundle_cache)
  region_data
}

# Load the base (sumstats_1) data.table from a bundled region RDS.
load_base_sumstats <- function(data_path, region_str) {
  region_data <- load_region_bundle(data_path, region_str)
  if (is.null(region_data)) return(NULL)
  region_data$base
}

# Look up a single trait's sumstats from a bundled region RDS by the
# trait's bundle key (e.g. "MVP_R4__study_file.txt.gz").
load_trait_sumstats <- function(data_path, region_str, trait_selector) {
  region_data <- load_region_bundle(data_path, region_str)
  if (is.null(region_data)) return(NULL)
  region_data$traits[[trait_selector]]
}

# For a virtual multi-ancestry study, find all per-ancestry region
# RDS files that overlap the consensus cluster containing region_str.
# The consensus coordinates are re-derived from the filesystem by
# scanning each ancestry's regional_dirs entry, so the overlay is
# always consistent with what's on disk.
#
# @param region_str A representative per-ancestry region key (e.g.
#   "16:19330554-21586583"). The loader seeds the consensus window
#   from this region and grows it to cover any overlapping file.
# @param virtual_info The DEFAULT_VIRTUAL_STUDIES entry for the
#   current study, i.e. list(category, ancestries, coloc_files,
#   regional_dirs, real_ids).
# @return A named list(ancestry = region_data) where each entry has
#   $base, $traits, and $.consensus = list(chr, start, stop).
load_multi_region_bundles <- function(region_str, virtual_info) {
  if (is.null(virtual_info)) return(NULL)

  rk <- parse_region_key(region_str)
  if (is.null(rk)) return(NULL)
  chr <- rk$chr; seed_start <- rk$start; seed_stop <- rk$stop

  study_tag <- paste(unlist(virtual_info$coloc_files), collapse = "|")
  cache_key <- paste(study_tag, chr, seed_start, seed_stop, sep = "||")
  cached <- .multi_bundle_cache[[cache_key]]
  if (!is.null(cached)) return(cached)

  # Collect all candidate files across ancestries on this chromosome.
  candidates <- list()
  for (anc in names(virtual_info$regional_dirs)) {
    d <- virtual_info$regional_dirs[[anc]]
    if (is.na(d) || !dir.exists(d)) next
    fs <- list.files(d, pattern = "\\.RDS$", full.names = FALSE)
    for (f in fs) {
      coord <- parse_regional_filename(f)
      if (is.null(coord)) next
      if (sub("^chr", "", coord$chr) != chr) next
      candidates[[length(candidates) + 1]] <- list(
        ancestry = anc, file = file.path(d, f),
        start = coord$start, stop = coord$stop)
    }
  }
  if (length(candidates) == 0) return(NULL)
  cand_dt <- data.table::rbindlist(candidates)

  # Grow the consensus window iteratively.
  cur_start <- seed_start; cur_stop <- seed_stop
  repeat {
    hit <- cand_dt[start <= cur_stop & stop >= cur_start]
    if (nrow(hit) == 0) break
    new_start <- min(hit$start); new_stop <- max(hit$stop)
    if (new_start == cur_start && new_stop == cur_stop) break
    cur_start <- new_start; cur_stop <- new_stop
  }
  hit <- cand_dt[start <= cur_stop & stop >= cur_start]
  if (nrow(hit) == 0) return(NULL)

  # Within each ancestry keep the widest overlapping file.
  hit[, width := stop - start]
  data.table::setorder(hit, ancestry, -width)
  hit <- hit[, .SD[1], by = ancestry]

  # Read and reconstruct each bundle.
  result <- list()
  for (i in seq_len(nrow(hit))) {
    rd <- tryCatch(readRDS(hit$file[i]), error = function(e) NULL)
    if (is.null(rd)) next
    if (!is.null(rd$base)) rd$base <- reconstruct_name(rd$base)
    rd$traits <- lapply(rd$traits, reconstruct_name)
    rd$.consensus <- list(chr = chr, start = cur_start, stop = cur_stop)
    result[[hit$ancestry[i]]] <- rd
  }

  assign(cache_key, result, envir = .multi_bundle_cache)
  result
}
