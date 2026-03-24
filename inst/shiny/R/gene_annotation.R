# Gene annotation loading and display

#' Load gene annotation data
#'
#' @param path Path to gene annotation file (TSV or TSV.GZ)
#' @return data.table with gene annotations, or NULL if not found
load_gene_annotation <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    message("Gene annotation file not found: ", path)
    return(NULL)
  }

  tryCatch({
    dt <- data.table::fread(path)
    message("Loaded gene annotation: ", nrow(dt), " genes")
    dt
  }, error = function(e) {
    message("Error loading gene annotation: ", e$message)
    NULL
  })
}
