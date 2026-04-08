#!/usr/bin/env Rscript
# Build the slim annotation DB bundled with the Shiny app.
#
# Reads the internal annotation database (genepicoloc_internal/annotations/
# annotation_db.RDS), strips per-study to a hardcoded column whitelist that
# excludes cluster paths and irrelevant QC metrics, and writes
# inst/shiny/extdata/shiny_annotation_db.RDS.
#
# This is the single source of truth for trait/study labels in the Shiny app.
# It lives under inst/shiny/ (not inst/extdata/) so the Dockerfile build
# context and the deploy.sh bind-mount can reach it. Re-run after the
# internal DB is regenerated. Output is committed to git.
#
# Usage:
#   Rscript inst/scripts/build_shiny_annotation_db.R [internal_db_path]

args <- commandArgs(trailingOnly = TRUE)
internal_db_path <- if (length(args) >= 1) {
  args[1]
} else {
  "~/Work/bioinfo/genepicoloc/genepicoloc_internal/annotations/annotation_db.RDS"
}
internal_db_path <- path.expand(internal_db_path)

if (!file.exists(internal_db_path)) {
  stop("Internal annotation DB not found: ", internal_db_path)
}

# Resolve this script's directory (inst/scripts/) so the output always lands
# in inst/extdata/, regardless of the caller's working directory.
get_script_dir <- function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cargs, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  # Sourced interactively
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(normalizePath(sys.frame(1)$ofile)))
  }
  getwd()
}
script_dir <- get_script_dir()
out_path <- file.path(script_dir, "..", "shiny", "extdata",
                      "shiny_annotation_db.RDS")
out_path <- normalizePath(out_path, mustWork = FALSE)

message("Reading internal DB: ", internal_db_path)
db_in <- readRDS(internal_db_path)

# Hardcoded column whitelist per study. Keys not in this list are dropped
# (e.g. gencode_v26, which is duplicated by the bundled gene annotation TSV).
keep_cols <- list(
  GCKD_mGWAS_plasma = c("merge_column", "BIOCHEMICAL"),
  GCKD_mGWAS_urine  = c("merge_column", "BIOCHEMICAL"),
  UKB_TOPMed = c("merge_column", "phenostring", "category",
                 "num_cases", "num_controls", "num_samples"),
  FinnGen_r9 = c("merge_column", "phenotype", "category",
                 "num_cases", "num_controls"),
  CKDGen_r4 = c("merge_column", "Name"),
  MVP_R4 = c("merge_column", "Title.of.analysis", "Analyzed.variable",
             "Phenotypic.trait.type", "Sample.size",
             "sumstats_2_type", "ancestry"),
  CKDGen_r5 = c("merge_column", "ckdgen_r5_name"),
  pho_ca = c("merge_column", "pho_ca_name"),
  UACR = c("merge_column", "Name"),
  olink_protein_map = c("UKBPPP_ProteinID", "Assay", "HGNC.symbol",
                        "UniProt", "olink_target_fullname", "Panel"),
  Icelanders_pGWAS_annotation = c("SeqId", "Protein..short.name.",
                                  "Protein..full.name.", "Gene",
                                  "UniProt", "Type")
  # gencode_v26: intentionally excluded - the Shiny app uses
  # inst/extdata/gene_annotation_hg38_full.tsv for gene annotation.
)

db_out <- list()
for (study in names(keep_cols)) {
  if (!study %in% names(db_in)) {
    warning("Study not in internal DB, skipping: ", study)
    next
  }
  x <- db_in[[study]]
  cols <- keep_cols[[study]]
  missing <- setdiff(cols, colnames(x))
  if (length(missing) > 0) {
    stop("Study ", study, " is missing expected columns: ",
         paste(missing, collapse = ", "))
  }
  db_out[[study]] <- x[, cols, drop = FALSE]
  message(sprintf("  %-30s %5d rows x %2d cols",
                  study, nrow(db_out[[study]]), ncol(db_out[[study]])))
}

# Audit: confirm no path-like strings remain in any column
path_pattern <- "^/(g|data|home|dsk|srv)/"
for (study in names(db_out)) {
  for (col in colnames(db_out[[study]])) {
    vals <- db_out[[study]][[col]]
    if (is.character(vals) && any(grepl(path_pattern, vals))) {
      stop("Path leak in ", study, "$", col)
    }
  }
}
message("Audit passed: no cluster paths in output.")

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(db_out, out_path, compress = "xz")
size_kb <- round(file.info(out_path)$size / 1024, 1)
message("Wrote ", out_path, " (", size_kb, " KB, ",
        length(db_out), " studies)")
