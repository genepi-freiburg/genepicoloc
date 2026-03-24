#!/usr/bin/env Rscript
# Launch genepicoloc Shiny app
#
# Usage:
#   Rscript run.R                                          # uses data/ subfolder
#   Rscript run.R ~/Work/bioinfo/genepicoloc/data/atlas   # custom data path

library(shiny)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  data_path <- normalizePath(args[1], mustWork = TRUE)
  options(genepicoloc.data_path = data_path)
  cat("Using data from:", data_path, "\n")
}

# Find app directory (same dir as this script)
script_dir <- tryCatch({
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
}, error = function(e) ".")

cat("Starting genepicoloc Shiny app...\n")
runApp(script_dir, port = 3838, host = "127.0.0.1", launch.browser = TRUE)
