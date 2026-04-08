# Gene annotation and gene track visualization

# Load gene annotation for gene tracks (Gencode v49 + HGNC + HPO + NCBI + Reactome)
# Files in R/ are auto-sourced BEFORE app.R, so library(data.table) from app.R
# is not yet attached. Use data.table::fread explicitly here.
# Paths are absolute (DATA_PATH) to avoid wd-dependence under shiny runApp.
gene_annotation <- local({
  gene_paths <- c(
    file.path(DATA_PATH, "gene_annotation.tsv"),
    file.path(DATA_PATH, "gene_annotation_hg38_full.tsv"),
    system.file("extdata", "genes_chr.txt.gz", package = "genepicoloc")
  )
  gene_file <- gene_paths[file.exists(gene_paths)][1]
  if (is.na(gene_file)) {
    warning("gene_annotation not found. Checked: ",
            paste(gene_paths, collapse = ", "))
    return(NULL)
  }
  message("Loaded gene_annotation from ", gene_file)
  data.table::fread(gene_file)
})

# Gene track plot function
# Creates a plotly gene track showing genes in a genomic region
plot_gene_track <- function(region_chr, region_start, region_end, max_genes = 30, source = NULL) {
  if (is.null(gene_annotation)) {
    return(plotly_empty() %>% layout(title = "Gene annotation not available"))
  }

  # Filter genes in region
  chr_clean <- gsub("^chr", "", as.character(region_chr))
  genes_in_region <- gene_annotation[
    chr == chr_clean &
    start <= region_end &
    end >= region_start
  ]

  if (nrow(genes_in_region) == 0) {
    return(plotly_empty() %>%
             layout(title = "No genes in region",
                    xaxis = list(range = c(region_start, region_end) / 1e6)))
  }

  # Limit number of genes for readability
  if (nrow(genes_in_region) > max_genes) {
    # Prioritize genes that overlap more with the region
    genes_in_region[, overlap := pmin(end, region_end) - pmax(start, region_start)]
    setorder(genes_in_region, -overlap)
    genes_in_region <- genes_in_region[1:max_genes]
    genes_in_region[, overlap := NULL]
  }

  # Assign rows to avoid overlapping labels (simple greedy algorithm)
  setorder(genes_in_region, start)
  genes_in_region[, row := 1L]

  if (nrow(genes_in_region) > 1) {
    row_ends <- c(genes_in_region$end[1])
    for (i in 2:nrow(genes_in_region)) {
      gene_start <- genes_in_region$start[i]
      # Find first row where gene fits (with padding for label)
      placed <- FALSE
      for (r in seq_along(row_ends)) {
        if (gene_start > row_ends[r] + (region_end - region_start) * 0.02) {  # 2% padding
          genes_in_region[i, row := r]
          row_ends[r] <- genes_in_region$end[i]
          placed <- TRUE
          break
        }
      }
      if (!placed) {
        # Create new row
        genes_in_region[i, row := length(row_ends) + 1L]
        row_ends <- c(row_ends, genes_in_region$end[i])
      }
    }
  }

  n_rows <- max(genes_in_region$row)

  # Create plot (source for click event capture)
  p <- plot_ly(source = source)

  # Add gene rectangles and arrows
  for (i in seq_len(nrow(genes_in_region))) {
    g <- genes_in_region[i]
    y_pos <- -g$row + 0.5

    # Gene body (rectangle)
    gene_start_mb <- max(g$start, region_start) / 1e6
    gene_end_mb <- min(g$end, region_end) / 1e6

    # Build hover text with basic info only (click for full details)
    hover_parts <- paste0("<b>", g$gene_name, "</b>",
                          "<br>chr", g$chr, ":", format(g$start, big.mark = ","),
                          "-", format(g$end, big.mark = ","))

    # Add gene as a line with markers for direction
    p <- p %>% add_trace(
      x = c(gene_start_mb, gene_end_mb),
      y = c(y_pos, y_pos),
      type = "scatter",
      mode = "lines",
      line = list(color = "#2c3e50", width = 6),
      hoverinfo = "text",
      text = hover_parts,
      customdata = g$gene_name,  # For click handler
      showlegend = FALSE
    )

    # Add direction arrow
    arrow_x <- if (g$strand == "+") gene_end_mb else gene_start_mb
    p <- p %>% add_trace(
      x = arrow_x,
      y = y_pos,
      type = "scatter",
      mode = "markers",
      marker = list(
        symbol = if (g$strand == "+") "triangle-right" else "triangle-left",
        size = 8,
        color = "#2c3e50"
      ),
      hoverinfo = "skip",
      showlegend = FALSE
    )

    # Add gene label
    label_x <- (gene_start_mb + gene_end_mb) / 2
    p <- p %>% add_annotations(
      x = label_x,
      y = y_pos + 0.3,
      text = paste0("<i>", g$gene_name, "</i>"),
      showarrow = FALSE,
      font = list(size = 9, color = "#2c3e50"),
      xanchor = "center",
      yanchor = "bottom"
    )
  }

  # Layout
  p <- p %>% layout(
    xaxis = list(
      title = "",
      range = c(region_start / 1e6, region_end / 1e6),
      showgrid = FALSE,
      zeroline = FALSE
    ),
    yaxis = list(
      title = "",
      showticklabels = FALSE,
      showgrid = FALSE,
      zeroline = FALSE,
      range = c(-n_rows - 0.5, 1)
    ),
    margin = list(t = 5, b = 5, l = 60, r = 20),
    showlegend = FALSE,
    hovermode = "closest"
  ) %>%
    config(displayModeBar = FALSE)

  # Register click event if source specified
  if (!is.null(source)) {
    p <- p %>% event_register("plotly_click")
  }

  return(p)
}

# Set study paths and available studies
STUDY_BASE_PATH <- DEFAULT_STUDY_BASE_PATH
available_studies <- DEFAULT_AVAILABLE_STUDIES

# Helper: resolve RDS file path (supports both flat and nested layouts).
# Virtual multi-ancestry studies don't have a single file; return the first
# available ancestry file as a sentinel so file.exists() checks in callers
# (e.g. observeEvent(selected_study)) succeed. The real multi-file load
# happens in coloc_data() via is_virtual_study().
resolve_annot_path <- function(study_name) {
  if (exists("DEFAULT_VIRTUAL_STUDIES") &&
      study_name %in% names(DEFAULT_VIRTUAL_STUDIES)) {
    paths <- unlist(DEFAULT_VIRTUAL_STUDIES[[study_name]]$coloc_files)
    existing <- paths[file.exists(paths)]
    if (length(existing) > 0) return(existing[1])
    return(NULL)
  }
  entry <- available_studies[[study_name]]
  if (is.null(entry)) return(NULL)
  # If entry is a full path to RDS (flat layout)
  if (grepl("\\.RDS$", entry) && file.exists(entry)) return(entry)
  # Nested layout: folder/annot/annot_filt.RDS
  file.path(STUDY_BASE_PATH, entry, "annot", "annot_filt.RDS")
}

# Is a given id a virtual multi-ancestry study?
is_virtual_study <- function(study_name) {
  exists("DEFAULT_VIRTUAL_STUDIES") &&
    !is.null(study_name) &&
    study_name %in% names(DEFAULT_VIRTUAL_STUDIES)
}

