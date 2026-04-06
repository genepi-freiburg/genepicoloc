# Convergence view: gene-centric region view
# Shows trait of interest + gene track + gene/trait coloc matrix + stacked RAPs

# =============================================================================
# Gene-attribution helpers
# =============================================================================

# Extract gene symbol from a coloc row if this is a cis molecular coloc
# Returns NA for trans signals or non-molecular studies (PheWAS/imaging/etc.)
# Input: a data.table row (or 1-row data.table)
get_cis_gene <- function(row) {
  study <- row$source_study[1]

  # Each molecular study has its own gene and cis_trans columns
  gene_col   <- NULL
  cistrans_col <- NULL

  if (study == "eQTLGen") {
    gene_col <- "eQTLGen_gene_name"
    cistrans_col <- "eQTLGen_cis_trans"
  } else if (study == "Kidney_eQTL") {
    gene_col <- "Kidney_eQTL_gene_name"
    cistrans_col <- "Kidney_eQTL_cis_trans"
  } else if (study == "GTEXv8_eQTL") {
    gene_col <- "GTEXv8_eQTL_gene_name"
    cistrans_col <- "GTEXv8_eQTL_cis_trans"
  } else if (study == "UKB_PPP_EUR") {
    gene_col <- "UKB_PPP_EUR_HGNC.symbol"
    cistrans_col <- "UKB_PPP_EUR_cis_trans"
  } else if (study == "Icelanders_pGWAS") {
    gene_col <- "Icelanders_pGWAS_Gene"
    cistrans_col <- "Icelanders_pGWAS_cis_trans"
  }

  if (is.null(gene_col) || !gene_col %in% names(row)) return(NA_character_)

  gene <- row[[gene_col]][1]
  if (is.null(gene) || is.na(gene)) return(NA_character_)

  # Only return gene if cis (or if cis_trans column missing, assume cis)
  if (!is.null(cistrans_col) && cistrans_col %in% names(row)) {
    ct <- row[[cistrans_col]][1]
    if (is.na(ct) || ct != "cis") return(NA_character_)
  }

  # Clean: take first gene if multiple (UKB_PPP has some combinations)
  gene <- strsplit(gene, "[,;/]")[[1]][1]
  trimws(gene)
}

# Apply get_cis_gene to an entire data.table, vectorized
annotate_with_cis_gene <- function(dt) {
  if (nrow(dt) == 0) {
    dt$cis_gene <- character(0)
    return(dt)
  }
  dt$cis_gene <- vapply(seq_len(nrow(dt)),
                         function(i) get_cis_gene(dt[i]),
                         character(1))
  dt
}

# Get trait label for display (reuses get_trait_display_name if available)
get_trait_label <- function(row, fallback_fn = NULL) {
  if (!is.null(fallback_fn)) {
    res <- tryCatch(fallback_fn(row), error = function(e) NULL)
    if (!is.null(res) && !is.na(res) && nchar(res) > 0) return(res)
  }
  # Fallback: basename of sumstats_2_file
  if ("sumstats_2_file" %in% names(row)) {
    return(sub("\\.[^.]+(\\.gz)?$", "", basename(row$sumstats_2_file[1])))
  }
  "Unknown"
}

# =============================================================================
# Gene track data
# =============================================================================

# Get protein-coding genes in a region from gene_annotation
get_genes_in_region <- function(gene_annotation, chr, start, end, pad = 0) {
  if (is.null(gene_annotation)) return(NULL)
  chr_str <- sub("^chr", "", as.character(chr))
  # Gene annotation may or may not have "chr" prefix
  gene_chr <- sub("^chr", "", as.character(gene_annotation$chr))
  in_region <- gene_chr == chr_str &
               gene_annotation$end >= (start - pad) &
               gene_annotation$start <= (end + pad)
  genes <- gene_annotation[in_region, ]
  if (nrow(genes) == 0) return(NULL)
  genes
}

# =============================================================================
# Convergence matrix (gene x trait)
# =============================================================================

# Build a long-form table: one row per (gene, trait, PP.H4)
# Separates gene-attributed (cis) vs locus-level rows
build_convergence_table <- function(filtered_region_data, trait_name_fn = NULL) {
  dt <- filtered_region_data
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  dt <- annotate_with_cis_gene(dt)

  # Get trait label per row
  dt$trait_label <- vapply(seq_len(nrow(dt)),
                            function(i) get_trait_label(dt[i], trait_name_fn),
                            character(1))

  # Level: "gene" if cis_gene present, else "locus"
  dt$evidence_level <- ifelse(!is.na(dt$cis_gene) & nchar(dt$cis_gene) > 0,
                               "gene", "locus")

  dt[, .(source_study, trait_label, cis_gene, evidence_level,
         PP.H4.abf, sumstats_2_max_nlog10P, directionality,
         sumstats_2_file)]
}

# =============================================================================
# Plots
# =============================================================================

# Gene track is provided by plot_gene_track() in gene_track.R (reused).

# NOTE: plot_gene_track_compact below is unused but kept for future refinements
# with gene highlighting based on cis colocs. Currently the app uses
# plot_gene_track() from gene_track.R.
plot_gene_track_compact <- function(genes, region_start, region_end,
                                     highlight_genes = character(0)) {
  if (is.null(genes) || nrow(genes) == 0) {
    return(plot_ly(
      x = c(region_start, region_end),
      y = c(0, 1),
      type = "scatter", mode = "markers",
      marker = list(opacity = 0, size = 1),
      hoverinfo = "none"
    ) %>% layout(
      xaxis = list(range = c(region_start, region_end), title = "",
                   showgrid = FALSE, zeroline = FALSE),
      yaxis = list(visible = FALSE),
      showlegend = FALSE,
      margin = list(t = 5, b = 20, l = 5, r = 5),
      annotations = list(list(text = "No genes in annotation",
                               x = (region_start + region_end) / 2, y = 0.5,
                               showarrow = FALSE))
    ))
  }

  # Simple layout: alternate genes to two rows to avoid overlap
  genes <- genes[order(genes$start), ]
  genes$row <- rep(c(0, 1), length.out = nrow(genes))

  # Color by whether highlighted
  is_highlight <- genes$gene_name %in% highlight_genes
  genes$color <- ifelse(is_highlight, "#e74c3c", "#95a5a6")

  # Draw genes as a single scatter trace with wide rectangles per gene
  # Use error bars as rectangles: x centered at midpoint, y at row, with
  # error_x spanning full gene length. Simpler: use bar-like markers.
  # Best approach: one scatter trace with markers for midpoints + line segments.
  # Easiest: use 'bar' trace with base + width... but bar is vertical.
  #
  # Alternative: multiple scatter traces (one per gene as line) - works reliably.
  hover_text <- paste0(
    "<b>", genes$gene_name, "</b><br>",
    "chr", genes$chr, ":", format(genes$start, big.mark = ","),
    "-", format(genes$end, big.mark = ","), "<br>",
    ifelse(!is.na(genes$full_name), genes$full_name, "")
  )

  # Build one trace per row (two rows, to stack non-overlapping)
  p <- plot_ly()
  for (r in c(0, 1)) {
    row_genes <- genes[genes$row == r, ]
    if (nrow(row_genes) == 0) next

    # Each gene as a horizontal segment: x = c(start, end), y = c(r, r)
    # Use NA separators to draw multiple segments in one trace
    xs <- as.numeric(rbind(row_genes$start, row_genes$end, NA))
    ys <- rep(r, length(xs))
    # Color - use gray for mixed; individual coloring would need per-trace split
    p <- p %>% add_trace(
      x = xs, y = ys,
      type = "scatter", mode = "lines",
      line = list(color = "#95a5a6", width = 10),
      hoverinfo = "none", showlegend = FALSE
    )

    # Highlight segments - separate trace
    hl_genes <- row_genes[row_genes$gene_name %in% highlight_genes, ]
    if (nrow(hl_genes) > 0) {
      xs_hl <- as.numeric(rbind(hl_genes$start, hl_genes$end, NA))
      ys_hl <- rep(r, length(xs_hl))
      p <- p %>% add_trace(
        x = xs_hl, y = ys_hl,
        type = "scatter", mode = "lines",
        line = list(color = "#e74c3c", width = 12),
        hoverinfo = "none", showlegend = FALSE
      )
    }

    # Midpoint markers for hover
    mid_x <- (row_genes$start + row_genes$end) / 2
    mid_y <- rep(r, nrow(row_genes))
    row_hover <- paste0(
      "<b>", row_genes$gene_name, "</b><br>",
      "chr", row_genes$chr, ":",
      format(row_genes$start, big.mark = ","), "-",
      format(row_genes$end, big.mark = ",")
    )
    p <- p %>% add_trace(
      x = mid_x, y = mid_y,
      type = "scatter", mode = "markers",
      marker = list(size = 1, opacity = 0),
      text = row_hover, hoverinfo = "text",
      showlegend = FALSE
    )
  }

  # Gene name labels: highlighted + top by width
  label_mask <- genes$gene_name %in% highlight_genes
  if (sum(label_mask) < 10) {
    widths <- genes$end - genes$start
    top_idx <- order(-widths)[1:min(10, nrow(genes))]
    label_mask[top_idx] <- TRUE
  }

  annotations <- lapply(which(label_mask), function(i) {
    list(
      x = (genes$start[i] + genes$end[i]) / 2,
      y = genes$row[i],
      text = genes$gene_name[i],
      showarrow = FALSE,
      font = list(size = 9,
                  color = if (is_highlight[i]) "#e74c3c" else "#333"),
      yshift = if (genes$row[i] == 0) -14 else 14
    )
  })

  p %>% layout(
    annotations = annotations,
    xaxis = list(range = c(region_start, region_end),
                 title = paste0("chr", sub("^chr", "", genes$chr[1]),
                                " position (bp)"),
                 showgrid = FALSE,
                 tickformat = ",d"),
    yaxis = list(range = c(-1, 2), visible = FALSE, fixedrange = TRUE),
    showlegend = FALSE,
    margin = list(t = 5, b = 35, l = 5, r = 5),
    hovermode = "closest"
  ) %>% config(displayModeBar = FALSE)
}

# Convergence heatmap: rows = genes, columns = individual traits
# Only shows gene-attributed (cis) evidence
plot_convergence_heatmap <- function(convergence_table) {
  gene_rows <- convergence_table[evidence_level == "gene"]
  if (is.null(gene_rows) || nrow(gene_rows) == 0) {
    return(plotly_empty() %>% layout(
      title = list(text = "No cis molecular colocalizations in this region",
                    font = list(size = 12))
    ))
  }

  # Use trait label + study for unique column ID
  gene_rows[, col_id := paste0(trait_label, " (", source_study, ")")]

  # Pivot to matrix: rows = genes, cols = traits
  genes <- sort(unique(gene_rows$cis_gene))
  traits <- unique(gene_rows$col_id)

  # Order traits: group by study, then by best PP.H4
  study_order <- gene_rows[, .(max_pp = max(PP.H4.abf)), by = source_study]
  study_order <- study_order[order(-max_pp)]
  gene_rows[, source_study := factor(source_study, levels = study_order$source_study)]
  setorder(gene_rows, source_study, -PP.H4.abf)
  traits <- unique(gene_rows$col_id)

  # Build matrix
  mat <- matrix(NA_real_, nrow = length(genes), ncol = length(traits),
                 dimnames = list(genes, traits))
  for (i in seq_len(nrow(gene_rows))) {
    g <- gene_rows$cis_gene[i]
    t <- gene_rows$col_id[i]
    # Take max PP.H4 if multiple
    cur <- mat[g, t]
    new <- gene_rows$PP.H4.abf[i]
    if (is.na(cur) || new > cur) mat[g, t] <- new
  }

  # Hover text matrix
  hover_mat <- matrix("", nrow = length(genes), ncol = length(traits),
                       dimnames = list(genes, traits))
  for (i in seq_len(nrow(gene_rows))) {
    g <- gene_rows$cis_gene[i]
    t <- gene_rows$col_id[i]
    hover_mat[g, t] <- paste0(
      "<b>", g, "</b><br>",
      t, "<br>",
      "PP.H4: ", round(gene_rows$PP.H4.abf[i], 3)
    )
  }

  plot_ly(
    x = traits,
    y = genes,
    z = mat,
    type = "heatmap",
    colorscale = list(c(0, "#f8f9fa"), c(0.5, "#fde68a"), c(1, "#27ae60")),
    zmin = 0.8, zmax = 1,
    text = hover_mat,
    hoverinfo = "text",
    showscale = TRUE,
    colorbar = list(title = "PP.H4", thickness = 10, len = 0.6)
  ) %>% layout(
    xaxis = list(title = "", tickangle = -45, tickfont = list(size = 9)),
    yaxis = list(title = "", tickfont = list(size = 11, color = "#333")),
    margin = list(t = 10, b = 120, l = 80, r = 10)
  ) %>% config(displayModeBar = FALSE)
}

# Locus-level traits bar chart (one row, colored by study)
plot_locus_level_bar <- function(convergence_table) {
  locus_rows <- convergence_table[evidence_level == "locus"]
  if (is.null(locus_rows) || nrow(locus_rows) == 0) {
    return(plotly_empty() %>% layout(
      title = list(text = "No locus-level traits in this region",
                    font = list(size = 12))
    ))
  }

  # Group by study, count traits
  summary <- locus_rows[, .(n_traits = .N, max_pp = max(PP.H4.abf)),
                        by = source_study]
  setorder(summary, -n_traits)

  # Colors from study_colors if available, else fallback
  cols <- sapply(summary$source_study, function(s) {
    if (exists("study_colors") && s %in% names(study_colors)) {
      study_colors[[s]]
    } else {
      "#7f7f7f"
    }
  })

  plot_ly(
    x = summary$source_study,
    y = summary$n_traits,
    type = "bar",
    marker = list(color = cols),
    text = paste0(summary$n_traits, " traits<br>max PP.H4: ", round(summary$max_pp, 2)),
    hoverinfo = "text"
  ) %>% layout(
    xaxis = list(title = "", tickangle = -30),
    yaxis = list(title = "# colocalized traits"),
    margin = list(t = 10, b = 60, l = 50, r = 10)
  ) %>% config(displayModeBar = FALSE)
}

