# YAML-based study configuration system ----

#' Load study configurations from a YAML file
#'
#' @description Reads a YAML configuration file defining datasets for
#' colocalization analysis. Each study entry specifies file paths, query
#' functions, and metadata needed to build an `args_df` for
#' [genepicoloc_wrapper()].
#'
#' @param config_file Character. Path to a YAML configuration file.
#'   See `system.file("extdata", "example_studies.yaml", package = "genepicoloc")`
#'   for an example.
#' @param internal Logical. If `TRUE`, include studies marked with
#'   `internal: true` (e.g., cluster-only datasets). Default: `FALSE`.
#' @param validate Logical. If `TRUE`, check that file paths exist and
#'   `expected_count` constraints are met. Default: `TRUE`.
#'
#' @return A named list of study configurations, each containing:
#'   \itemize{
#'     \item sumstats_2_file: Character vector of file paths
#'     \item sumstats_2_function: Query function name
#'     \item sumstats_2_type: Character vector ("quant" or "cc" per file)
#'     \item sumstats_2_sdY: Numeric vector (sdY per file, NA for cc)
#'     \item description: Human-readable study description
#'   }
#'
#' @examples
#' \dontrun{
#' config <- load_studies_config("my_studies.yaml")
#' names(config)
#' }
#'
#' @export
load_studies_config <- function(config_file, internal = FALSE, validate = TRUE) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for config file support. ",
         "Install it with: install.packages('yaml')")
  }

  if (!file.exists(config_file)) {
    stop("Config file not found: ", config_file)
  }

  config <- yaml::read_yaml(config_file)

  if (!is.list(config) || length(config) == 0) {
    stop("Config file must contain at least one study entry")
  }

  # Filter out internal studies unless requested
  if (!internal) {
    config <- config[!vapply(config, function(x) {
      isTRUE(x[["internal"]])
    }, logical(1))]
  }

  if (length(config) == 0) {
    stop("No studies available (all are marked internal). ",
         "Use internal=TRUE to include them.")
  }

  # Process each study entry
  result <- lapply(names(config), function(study_name) {
    process_study_entry(study_name, config[[study_name]], validate = validate)
  })
  names(result) <- names(config)

  result
}


#' Process a single study entry from YAML config
#'
#' @description Resolves file paths and metadata for one study definition.
#' Supports three file source modes: explicit file list, directory + pattern,
#' or a metadata TSV file.
#'
#' @param study_name Character. Study identifier.
#' @param entry Named list. Study configuration from YAML.
#' @param validate Logical. Whether to validate file existence and counts.
#'
#' @return Named list with sumstats_2_file, sumstats_2_function,
#'   sumstats_2_type, sumstats_2_sdY, description.
#'
#' @keywords internal
process_study_entry <- function(study_name, entry, validate = TRUE) {

  # Required field
  if (is.null(entry[["query_function"]])) {
    stop("Study '", study_name, "': 'query_function' is required")
  }
  if (is.null(entry[["type"]])) {
    stop("Study '", study_name, "': 'type' is required (quant or cc)")
  }

  query_function <- entry[["query_function"]]
  default_type <- entry[["type"]]
  default_sdY <- if (is.null(entry[["sdY"]])) NA_real_ else as.numeric(entry[["sdY"]])
  description <- if (is.null(entry[["description"]])) "" else entry[["description"]]

  # Resolve files from one of three modes
  if (!is.null(entry[["files"]])) {
    # Mode 1: Explicit file list
    files <- unlist(entry[["files"]])

  } else if (!is.null(entry[["path"]])) {
    # Mode 2: Directory + pattern
    dir_path <- entry[["path"]]
    pattern <- if (is.null(entry[["pattern"]])) "*.gz" else entry[["pattern"]]

    if (validate && !dir.exists(dir_path)) {
      stop("Study '", study_name, "': directory not found: ", dir_path)
    }

    files <- list.files(dir_path, pattern = utils::glob2rx(pattern),
                        full.names = TRUE)

    # Apply exclude_pattern if set
    if (!is.null(entry[["exclude_pattern"]])) {
      files <- files[!grepl(entry[["exclude_pattern"]], basename(files))]
    }

  } else if (!is.null(entry[["metadata_file"]])) {
    # Mode 3: Metadata TSV (MVP-style, per-file metadata)
    meta_file <- entry[["metadata_file"]]
    if (validate && !file.exists(meta_file)) {
      stop("Study '", study_name, "': metadata file not found: ", meta_file)
    }

    meta <- data.table::fread(meta_file)

    # Expect columns: file, type, sdY (at minimum 'file')
    if (!"file" %in% names(meta)) {
      stop("Study '", study_name, "': metadata file must have a 'file' column")
    }

    files <- meta$file

    # Per-file type and sdY from metadata
    file_types <- if ("type" %in% names(meta)) meta$type else rep(default_type, length(files))
    file_sdY <- if ("sdY" %in% names(meta)) as.numeric(meta$sdY) else rep(default_sdY, length(files))

    return(list(
      sumstats_2_file = files,
      sumstats_2_function = query_function,
      sumstats_2_type = file_types,
      sumstats_2_sdY = file_sdY,
      description = description
    ))

  } else {
    stop("Study '", study_name, "': must specify 'files', 'path' + 'pattern', ",
         "or 'metadata_file'")
  }

  # Validate files exist
  if (validate && length(files) > 0) {
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) {
      warning("Study '", study_name, "': ", length(missing), " of ",
              length(files), " files not found. First missing: ", missing[1])
    }
  }

  if (length(files) == 0) {
    warning("Study '", study_name, "': no files found")
  }

  # Check expected_count
  if (!is.null(entry[["expected_count"]])) {
    expected <- as.integer(entry[["expected_count"]])
    if (length(files) != expected) {
      msg <- sprintf("Study '%s': expected %d files but found %d",
                     study_name, expected, length(files))
      if (validate) stop(msg) else warning(msg)
    }
  }

  # Build type and sdY vectors
  file_types <- rep(default_type, length(files))
  file_sdY <- rep(default_sdY, length(files))

  # Apply type_rules for per-file overrides
  if (!is.null(entry[["type_rules"]])) {
    for (rule_type in names(entry[["type_rules"]])) {
      rule_pattern <- entry[["type_rules"]][[rule_type]]
      matches <- grepl(rule_pattern, basename(files))
      file_types[matches] <- rule_type
      # cc studies don't have sdY
      if (rule_type == "cc") {
        file_sdY[matches] <- NA_real_
      }
    }
  }

  list(
    sumstats_2_file = files,
    sumstats_2_function = query_function,
    sumstats_2_type = file_types,
    sumstats_2_sdY = file_sdY,
    description = description
  )
}


#' Create args_df from a YAML configuration file
#'
#' @description Main user-facing function that reads a YAML study configuration
#' and builds an `args_df` data.table compatible with [genepicoloc_wrapper()].
#'
#' @param config_file Character. Path to a YAML configuration file.
#' @param selected_studies Character vector. If provided, only include these
#'   studies. Default: `NULL` (all studies).
#' @param list_studies_only Logical. If `TRUE`, return a data.table summarizing
#'   available studies instead of the full args_df. Default: `FALSE`.
#' @param test_mode Logical. If `TRUE`, limit to max 15 files per study for
#'   quick testing. Default: `FALSE`.
#' @param internal Logical. If `TRUE`, include internal/cluster-only studies.
#'   Default: `FALSE`.
#' @param validate Logical. If `TRUE`, validate file paths and counts.
#'   Default: `TRUE`.
#'
#' @return A data.table with columns: sumstats_2_study, sumstats_2_file,
#'   sumstats_2_function, sumstats_2_type, sumstats_2_sdY. Or if
#'   `list_studies_only = TRUE`, a summary data.table with study names,
#'   descriptions, and file counts.
#'
#' @examples
#' \dontrun{
#' # List available studies
#' list_available_studies("my_studies.yaml")
#'
#' # Build args_df for all studies
#' args_df <- create_args_df_from_config("my_studies.yaml")
#'
#' # Build args_df for selected studies only
#' args_df <- create_args_df_from_config("my_studies.yaml",
#'   selected_studies = c("GTEx_eQTL", "eQTLGen"))
#'
#' # Quick test with limited files
#' args_df <- create_args_df_from_config("my_studies.yaml", test_mode = TRUE)
#' }
#'
#' @export
create_args_df_from_config <- function(config_file,
                                       selected_studies = NULL,
                                       list_studies_only = FALSE,
                                       test_mode = FALSE,
                                       internal = FALSE,
                                       validate = TRUE) {

  config <- load_studies_config(config_file, internal = internal,
                                validate = validate)

  # Filter to selected studies
  if (!is.null(selected_studies)) {
    missing_studies <- setdiff(selected_studies, names(config))
    if (length(missing_studies) > 0) {
      stop("Studies not found in config: ",
           paste(missing_studies, collapse = ", "))
    }
    config <- config[selected_studies]
  }

  # List-only mode: return summary
  if (list_studies_only) {
    summary_dt <- data.frame(
      study = names(config),
      description = vapply(config, function(x) x$description, character(1)),
      n_files = vapply(config, function(x) length(x$sumstats_2_file), integer(1)),
      query_function = vapply(config, function(x) x$sumstats_2_function, character(1)),
      type = vapply(config, function(x) {
        types <- unique(x$sumstats_2_type)
        paste(types, collapse = "/")
      }, character(1)),
      stringsAsFactors = FALSE
    )
    return(summary_dt)
  }

  # Build args_df rows for each study
  rows <- lapply(names(config), function(study_name) {
    entry <- config[[study_name]]
    files <- entry$sumstats_2_file
    types <- entry$sumstats_2_type
    sdYs <- entry$sumstats_2_sdY

    # Test mode: limit files
    if (test_mode && length(files) > 15) {
      files <- files[1:15]
      types <- types[1:15]
      sdYs <- sdYs[1:15]
    }

    data.frame(
      sumstats_2_study = study_name,
      sumstats_2_file = files,
      sumstats_2_function = entry$sumstats_2_function,
      sumstats_2_type = types,
      sumstats_2_sdY = sdYs,
      stringsAsFactors = FALSE
    )
  })

  args_df <- do.call(rbind, rows)

  message(sprintf("Created args_df: %d datasets across %d studies",
                  nrow(args_df), length(config)))
  for (study in names(config)) {
    n <- sum(args_df$sumstats_2_study == study)
    message(sprintf("  - %s: %d datasets", study, n))
  }

  args_df
}


#' List available studies from a YAML configuration file
#'
#' @description Convenience wrapper that displays the studies defined in a
#' YAML configuration file with their descriptions and file counts.
#'
#' @param config_file Character. Path to a YAML configuration file.
#' @param internal Logical. If `TRUE`, include internal/cluster-only studies.
#'   Default: `FALSE`.
#'
#' @return A data.table with columns: study, description, n_files,
#'   query_function, type.
#'
#' @examples
#' \dontrun{
#' list_available_studies("my_studies.yaml")
#' }
#'
#' @export
list_available_studies <- function(config_file, internal = FALSE) {
  create_args_df_from_config(config_file, list_studies_only = TRUE,
                             internal = internal, validate = FALSE)
}
