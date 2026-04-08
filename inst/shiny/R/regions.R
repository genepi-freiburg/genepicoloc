# Region parsing helpers - pure, no Shiny dependencies.
#
# Two string conventions are in play:
#
#   "region_str" (a.k.a. region_key in UI code):
#     the selectize value and the display identifier;
#     format "<chr>:<start>-<stop>" with chr optionally "chr"-prefixed
#     (e.g. "16:19330554-21586583" or "chr16:19330554-21586583")
#
#   "region_id":
#     the filesystem-friendly form of the region used as the basename
#     of the bundled regional RDS file;
#     format "chr<chr>_<start>_<stop>"
#     (e.g. "chr16_19330554_21586583")
#
# parse_region_key()       region_str -> list(chr, start, stop)
# parse_region_id()        region_str -> "chr<chr>_<start>_<stop>"
# parse_regional_filename() "chr<chr>_<start>_<stop>.RDS" -> list(chr, start, stop)

# Parse a region string "chr:start-end" into numeric components.
# Returns NULL when the string is missing or malformed. The `chr`
# component has the "chr" prefix stripped so callers can match against
# data.tables that store bare chromosome codes.
parse_region_key <- function(region_str) {
  if (is.null(region_str) || is.na(region_str) || !nzchar(region_str)) return(NULL)
  parts <- strsplit(region_str, ":", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(NULL)
  chr <- sub("^chr", "", parts[1])
  pos <- strsplit(parts[2], "-", fixed = TRUE)[[1]]
  if (length(pos) < 2) return(NULL)
  start <- suppressWarnings(as.numeric(pos[1]))
  stop  <- suppressWarnings(as.numeric(pos[2]))
  if (is.na(start) || is.na(stop)) return(NULL)
  list(chr = chr, start = start, stop = stop)
}

# Parse a region string into the filesystem-friendly form used as the
# basename of bundled regional RDS files (without the `.RDS` extension).
# Returns NULL for malformed input.
parse_region_id <- function(region_str) {
  rk <- parse_region_key(region_str)
  if (is.null(rk)) return(NULL)
  paste0("chr", rk$chr, "_", rk$start, "_", rk$stop)
}

# Parse a regional RDS filename back into numeric coordinates.
# Returns list(chr, start, stop) or NULL for non-matching filenames.
parse_regional_filename <- function(fname) {
  m <- regmatches(fname,
    regexec("^chr([0-9XY]+)_([0-9]+)_([0-9]+)\\.RDS$", fname))[[1]]
  if (length(m) != 4) return(NULL)
  list(chr = m[2], start = as.numeric(m[3]), stop = as.numeric(m[4]))
}

# Reconstruct the `Name` column (chr:pos:ref:alt) from the slim regional
# RDS format if it is missing. Idempotent - returns the input unchanged
# when `Name` is already present.
reconstruct_name <- function(dt) {
  if (!"Name" %in% names(dt) && all(c("CHR", "POS", "A1", "A2") %in% names(dt))) {
    dt$Name <- paste0("chr", dt$CHR, ":", dt$POS, ":", dt$A1, ":", dt$A2)
  }
  dt
}
