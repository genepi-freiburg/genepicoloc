# Regional association plot functions (LocusZoom-style)

# Interactive regional association plot using plotly
# Creates an interactive Manhattan-style plot with hover tooltips.
# Only variants with -log10(P) > 5 (p < 1e-5) have hover annotations for
# performance.
#
# If `sumstats_dt` is a **named list** of data.tables it is treated as an
# overlay request (one trace per element, named by ancestry). Marker
# colors cycle through ANCESTRY_COLORS.
plot_regional_association_interactive <- function(sumstats_dt, title = "", highlight_snp = NULL,
                                                   y_col = "nlog10P", color = "#3498db",
                                                   x_range = NULL) {
  # Overlay path: named list of data.tables (one per ancestry).
  if (is.list(sumstats_dt) && !is.data.frame(sumstats_dt)) {
    return(plot_regional_overlay(sumstats_dt, title = title, y_col = y_col,
                                 x_range = x_range))
  }

  if (is.null(sumstats_dt) || nrow(sumstats_dt) == 0) {
    # Return empty plotly with message
    return(plotly_empty() %>%
             layout(title = list(text = "No data available", x = 0.5, y = 0.5)))
  }

  # Ensure data.table
  if (!is.data.table(sumstats_dt)) {
    sumstats_dt <- as.data.table(sumstats_dt)
  }

  # Separate significant (p < 1e-5, i.e., -log10(P) > 5) and non-significant variants
  is_significant <- sumstats_dt[[y_col]] > 5
  is_lead <- if (!is.null(highlight_snp)) sumstats_dt$Name == highlight_snp else rep(FALSE, nrow(sumstats_dt))

  # Create hover text only for significant variants (for performance)
  hover_text_sig <- paste0(
    "<b>", sumstats_dt$Name[is_significant], "</b>",
    ifelse("rsID" %in% names(sumstats_dt) & !is.na(sumstats_dt$rsID[is_significant]),
           paste0("<br>rsID: ", sumstats_dt$rsID[is_significant]), ""),
    "<br>Position: ", format(sumstats_dt$POS[is_significant], big.mark = ","),
    "<br>-log10(P): ", round(sumstats_dt[[y_col]][is_significant], 2),
    ifelse("BETA" %in% names(sumstats_dt) & !is.na(sumstats_dt$BETA[is_significant]),
           paste0("<br>Beta: ", round(sumstats_dt$BETA[is_significant], 4)), ""),
    ifelse("SE" %in% names(sumstats_dt) & !is.na(sumstats_dt$SE[is_significant]),
           paste0("<br>SE: ", round(sumstats_dt$SE[is_significant], 4)), ""),
    ifelse("AF" %in% names(sumstats_dt) & !is.na(sumstats_dt$AF[is_significant]),
           paste0("<br>AF: ", round(sumstats_dt$AF[is_significant], 3)), "")
  )

  # Create hover text for lead SNP
  hover_text_lead <- if (any(is_lead)) {
    paste0(
      "<b>", sumstats_dt$Name[is_lead], "</b>",
      ifelse("rsID" %in% names(sumstats_dt) & !is.na(sumstats_dt$rsID[is_lead]),
             paste0("<br>rsID: ", sumstats_dt$rsID[is_lead]), ""),
      "<br>Position: ", format(sumstats_dt$POS[is_lead], big.mark = ","),
      "<br>-log10(P): ", round(sumstats_dt[[y_col]][is_lead], 2),
      ifelse("BETA" %in% names(sumstats_dt) & !is.na(sumstats_dt$BETA[is_lead]),
             paste0("<br>Beta: ", round(sumstats_dt$BETA[is_lead], 4)), ""),
      ifelse("SE" %in% names(sumstats_dt) & !is.na(sumstats_dt$SE[is_lead]),
             paste0("<br>SE: ", round(sumstats_dt$SE[is_lead], 4)), ""),
      ifelse("AF" %in% names(sumstats_dt) & !is.na(sumstats_dt$AF[is_lead]),
             paste0("<br>AF: ", round(sumstats_dt$AF[is_lead], 3)), "")
    )
  } else NULL

  # Non-significant, non-lead variants (no hover)
  is_nonsig_nonlead <- !is_significant & !is_lead
  # Significant, non-lead variants (with hover)
  is_sig_nonlead <- is_significant & !is_lead

  # Create the plot - start with non-significant variants (no hover for performance)
  p <- plot_ly()

  # Add non-significant variants (no hover text for performance)
  if (any(is_nonsig_nonlead)) {
    p <- p %>%
      add_trace(
        data = sumstats_dt[is_nonsig_nonlead],
        x = ~POS,
        y = ~get(y_col),
        type = "scatter",
        mode = "markers",
        marker = list(
          color = color,
          size = 5,
          opacity = 0.4
        ),
        hoverinfo = "none",
        name = "Variants (p > 1e-5)"
      )
  }

  # Add significant variants (with hover)
  if (any(is_sig_nonlead)) {
    p <- p %>%
      add_trace(
        data = sumstats_dt[is_sig_nonlead],
        x = ~POS,
        y = ~get(y_col),
        type = "scatter",
        mode = "markers",
        marker = list(
          color = color,
          size = 7,
          opacity = 0.8
        ),
        text = hover_text_sig[!is_lead[is_significant]],
        hoverinfo = "text",
        name = "Variants (p < 1e-5)"
      )
  }

  # Add lead SNP as separate trace if present
  if (any(is_lead)) {
    p <- p %>%
      add_trace(
        data = sumstats_dt[is_lead],
        x = ~POS,
        y = ~get(y_col),
        type = "scatter",
        mode = "markers",
        marker = list(
          color = "purple",
          size = 12,
          symbol = "diamond"
        ),
        text = hover_text_lead,
        hoverinfo = "text",
        name = "Lead SNP"
      )
  }

  # Add significance lines and layout
  chr_label <- if ("CHR" %in% names(sumstats_dt)) sumstats_dt$CHR[1] else "?"

  # Set x-axis range (use provided range for synchronization, or auto)
  xaxis_config <- list(
    title = paste0("Position (Chromosome ", chr_label, ")"),
    tickformat = ",d"
  )
  if (!is.null(x_range)) {
    xaxis_config$range <- x_range
  }

  p <- p %>%
    layout(
      title = list(text = title, font = list(size = 14)),
      xaxis = xaxis_config,
      yaxis = list(title = "-log<sub>10</sub>(P)"),
      shapes = list(
        # Genome-wide significance line (7.3 = -log10(5e-8))
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
             y0 = 7.3, y1 = 7.3, line = list(color = "red", dash = "dash", width = 1)),
        # Suggestive significance line
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
             y0 = 5, y1 = 5, line = list(color = "blue", dash = "dot", width = 0.8))
      ),
      showlegend = FALSE,
      hovermode = "closest",
      margin = list(t = 40, b = 60, l = 60, r = 20)
    ) %>%
    config(displayModeBar = TRUE, displaylogo = FALSE,
           modeBarButtonsToRemove = c("lasso2d", "select2d"))

  return(p)
}

# Multi-ancestry overlay regional plot. Takes a named list of data.tables
# (one per ancestry, names are the ancestry codes) and produces a single
# plotly scatter with one trace per ancestry, colored via ANCESTRY_COLORS.
# If the list carries an attribute `coloc_ancestries`, ancestries NOT in
# that vector are rendered with reduced opacity and marked "(context)"
# in the legend - they show raw sumstats for comparison even though the
# trait did not colocalize in that ancestry.
plot_regional_overlay <- function(bundles, title = "", y_col = "nlog10P",
                                  x_range = NULL) {
  coloc_ancs <- attr(bundles, "coloc_ancestries")
  # Drop empty/missing entries.
  bundles <- bundles[vapply(bundles, function(x) {
    !is.null(x) && is.data.frame(x) && nrow(x) > 0
  }, logical(1))]
  if (length(bundles) == 0) {
    return(plotly_empty() %>%
             layout(title = list(text = "No data available", x = 0.5, y = 0.5)))
  }

  p <- plot_ly()
  chr_label <- NA

  # Stable ancestry order for consistent legend.
  anc_order <- intersect(names(ANCESTRY_COLORS), names(bundles))
  anc_order <- c(anc_order, setdiff(names(bundles), anc_order))

  for (anc in anc_order) {
    dt <- bundles[[anc]]
    if (!is.data.table(dt)) dt <- as.data.table(dt)
    col <- if (anc %in% names(ANCESTRY_COLORS)) ANCESTRY_COLORS[[anc]] else "#7f7f7f"
    if (is.na(chr_label) && "CHR" %in% names(dt)) chr_label <- dt$CHR[1]

    # Is this ancestry among the colocalized ones? If coloc_ancs is NULL
    # (attribute not set), treat all traces as colocalized.
    is_coloc <- is.null(coloc_ancs) || anc %in% coloc_ancs
    legend_name <- if (is_coloc) anc else paste0(anc, " (no coloc.)")
    bg_alpha  <- if (is_coloc) 0.35 else 0.12
    sig_alpha <- if (is_coloc) 0.85 else 0.40

    yvec <- dt[[y_col]]
    is_sig <- yvec > 5

    # Hover text only for significant variants (performance).
    hover_sig <- if (any(is_sig)) {
      paste0(
        "<b>", anc, "</b><br>",
        if ("Name" %in% names(dt)) dt$Name[is_sig] else "",
        "<br>Position: ", format(dt$POS[is_sig], big.mark = ","),
        "<br>-log10(P): ", round(yvec[is_sig], 2),
        if ("BETA" %in% names(dt))
          paste0("<br>Beta: ", round(dt$BETA[is_sig], 4)) else "",
        if ("AF" %in% names(dt))
          paste0("<br>AF: ", round(dt$AF[is_sig], 3)) else ""
      )
    } else character(0)

    # Non-significant background trace (no hover, small, transparent).
    if (any(!is_sig)) {
      p <- p %>% add_trace(
        data = dt[!is_sig],
        x = ~POS, y = ~get(y_col),
        type = "scatter", mode = "markers",
        marker = list(color = col, size = 4, opacity = bg_alpha),
        hoverinfo = "none",
        name = paste0(anc, " (p > 1e-5)"),
        legendgroup = anc,
        showlegend = FALSE
      )
    }

    # Significant variants with hover.
    if (any(is_sig)) {
      p <- p %>% add_trace(
        data = dt[is_sig],
        x = ~POS, y = ~get(y_col),
        type = "scatter", mode = "markers",
        marker = list(color = col, size = 7, opacity = sig_alpha),
        text = hover_sig, hoverinfo = "text",
        name = legend_name,
        legendgroup = anc,
        showlegend = TRUE
      )
    } else {
      # Add an invisible trace so the legend still lists this ancestry.
      p <- p %>% add_trace(
        x = numeric(0), y = numeric(0),
        type = "scatter", mode = "markers",
        marker = list(color = col, size = 7, opacity = sig_alpha),
        name = legend_name,
        legendgroup = anc,
        showlegend = TRUE
      )
    }
  }

  xaxis_config <- list(
    title = paste0("Position (Chromosome ", ifelse(is.na(chr_label), "?", chr_label), ")"),
    tickformat = ",d"
  )
  if (!is.null(x_range)) xaxis_config$range <- x_range

  p %>% layout(
    title = list(text = title, font = list(size = 14)),
    xaxis = xaxis_config,
    yaxis = list(title = "-log<sub>10</sub>(P)"),
    shapes = list(
      list(type = "line", x0 = 0, x1 = 1, xref = "paper",
           y0 = 7.3, y1 = 7.3, line = list(color = "red", dash = "dash", width = 1)),
      list(type = "line", x0 = 0, x1 = 1, xref = "paper",
           y0 = 5, y1 = 5, line = list(color = "blue", dash = "dot", width = 0.8))
    ),
    showlegend = TRUE,
    legend = list(orientation = "h", x = 0, y = 1.05),
    hovermode = "closest",
    margin = list(t = 40, b = 60, l = 60, r = 20)
  ) %>%
  config(displayModeBar = TRUE, displaylogo = FALSE,
         modeBarButtonsToRemove = c("lasso2d", "select2d"))
}

# Load gene annotation for gene tracks (Gencode v49 + HGNC + HPO + NCBI + Reactome)
# gene_annotation is loaded in gene_track.R

# Gene track plot function
