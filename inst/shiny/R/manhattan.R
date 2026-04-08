# Mini-Manhattan: pure table builder + pure plotly renderer.
#
# The mini Manhattan is the top-of-page "pick a region" plot. It reads
# two inputs from the server:
#
#   dt       - coloc_data(): one row per coloc hit
#   regions  - regions_data(): one row per (CHR_var, BP_START_var,
#              BP_STOP_var) with nearest_gene_1 / Prioritized_Gene /
#              nearest_genes_10 metadata
#
# and one optional argument:
#
#   is_multi - TRUE when the current study is a virtual multi-ancestry
#              one (coloc_data carries an `ancestry` column). Forward
#              consensus clustering and ancestry-coverage coloring
#              activate only in this mode.
#
# manhattan_consensus_table() returns a data.table with all columns the
# renderer needs (x_pos, chr_color, region_key, hover, gene, is_selected).
# The `selected_region_key` argument is the currently highlighted region
# key (may be NA/NULL); the function just tags a boolean column, no
# reactive reads.

# Build the full `region_stats` table used by the Manhattan renderer.
manhattan_consensus_table <- function(dt, regions, is_multi,
                                      selected_region_key = NA_character_) {
  if (isTRUE(is_multi)) {
    # Per (CHR, START, STOP) aggregation (per-ancestry regions may
    # repeat if multiple coloc hits fall into them).
    per_region <- dt[, .(max_nlog10P = max(sumstats_1_max_nlog10P, na.rm = TRUE),
                         n_coloc = .N,
                         ancs = list(sort(unique(ancestry)))),
                     by = .(CHR_var, BP_START_var, BP_STOP_var)]

    # Cluster overlapping regions within each chromosome (single pass,
    # sort by start, merge if next start <= running end).
    per_region <- per_region[order(CHR_var, BP_START_var)]
    per_region[, cluster_id := {
      cid <- integer(.N); running_end <- -Inf; cur <- 0L
      for (i in seq_len(.N)) {
        if (BP_START_var[i] > running_end) {
          cur <- cur + 1L
          running_end <- BP_STOP_var[i]
        } else if (BP_STOP_var[i] > running_end) {
          running_end <- BP_STOP_var[i]
        }
        cid[i] <- cur
      }
      cid
    }, by = CHR_var]

    # Representative per-ancestry region key per cluster = the one with
    # the most colocs (fallback: the widest). Click routing uses this
    # key so downstream filtered_region_data() still finds matching
    # rows without knowing about consensus coordinates.
    per_region[, width := BP_STOP_var - BP_START_var]
    per_region[, rep_row := order(-n_coloc, -width)[1] == seq_len(.N),
               by = .(CHR_var, cluster_id)]
    rep_keys <- per_region[rep_row == TRUE,
      .(CHR_var, cluster_id,
        rep_start = BP_START_var, rep_stop = BP_STOP_var)]

    # Collapse to consensus region per cluster.
    region_stats <- per_region[, .(
      BP_START_var = min(BP_START_var),
      BP_STOP_var  = max(BP_STOP_var),
      max_nlog10P  = max(max_nlog10P, na.rm = TRUE),
      n_coloc      = sum(n_coloc),
      ancs         = list(sort(unique(unlist(ancs))))
    ), by = .(CHR_var, cluster_id)]
    region_stats <- merge(region_stats, rep_keys,
                          by = c("CHR_var", "cluster_id"), all.x = TRUE)
    region_stats[, n_ancestries := lengths(ancs)]
    region_stats[, anc_label := vapply(ancs, paste, character(1), collapse = ", ")]
    region_stats[, cluster_id := NULL]
  } else {
    region_stats <- dt[, .(max_nlog10P = max(sumstats_1_max_nlog10P, na.rm = TRUE),
                           n_coloc = .N),
                       by = .(CHR_var, BP_START_var, BP_STOP_var)]
  }

  # Numeric chr for x-axis, cumulative layout positions
  region_stats[, CHR_num := as.integer(gsub("X", "23", gsub("Y", "24", CHR_var)))]
  region_stats[, mid_pos := (BP_START_var + BP_STOP_var) / 2]
  chr_offsets <- region_stats[, .(max_pos = max(mid_pos)), by = CHR_num]
  data.table::setorder(chr_offsets, CHR_num)
  chr_offsets[, offset := cumsum(c(0, head(max_pos, -1))) +
                          (seq_len(.N) - 1) * 5e7]
  region_stats <- merge(region_stats, chr_offsets[, .(CHR_num, offset)],
                        by = "CHR_num")
  region_stats[, x_pos := mid_pos + offset]

  # Gene labels from the per-region metadata
  if (isTRUE(is_multi)) {
    # Consensus region coordinates don't match the per-ancestry
    # BP_START_var in `regions`, so pick the per-ancestry row whose
    # midpoint is closest to each consensus region's midpoint on the
    # same chromosome. Also carry nearest_genes_10 through so the
    # tooltip can show the full gene neighborhood.
    cols <- c("CHR_var", "BP_START_var", "BP_STOP_var",
              "nearest_gene_1", "Prioritized_Gene")
    if ("nearest_genes_10" %in% names(regions)) cols <- c(cols, "nearest_genes_10")
    reg_meta <- unique(regions[, ..cols])
    reg_meta[, r_mid := (BP_START_var + BP_STOP_var) / 2]
    has_n10 <- "nearest_genes_10" %in% names(reg_meta)
    gene_for_cluster <- region_stats[, {
      cands <- reg_meta[CHR_var == .BY$CHR_var]
      if (nrow(cands) == 0) {
        list(nearest_gene_1 = NA_character_,
             Prioritized_Gene = NA_character_,
             nearest_genes_10 = NA_character_)
      } else {
        j <- which.min(abs(cands$r_mid - mid_pos))
        list(nearest_gene_1 = cands$nearest_gene_1[j],
             Prioritized_Gene = cands$Prioritized_Gene[j],
             nearest_genes_10 = if (has_n10) cands$nearest_genes_10[j] else NA_character_)
      }
    }, by = .(CHR_var, BP_START_var, BP_STOP_var)]
    region_stats <- merge(region_stats, gene_for_cluster,
                          by = c("CHR_var", "BP_START_var", "BP_STOP_var"),
                          all.x = TRUE)
  } else {
    region_stats <- merge(region_stats,
      unique(regions[, .(CHR_var, BP_START_var, nearest_gene_1, Prioritized_Gene)]),
      by = c("CHR_var", "BP_START_var"), all.x = TRUE)
  }
  region_stats[, gene := ifelse(!is.na(Prioritized_Gene), Prioritized_Gene, nearest_gene_1)]

  # Region click key: for multi-ancestry virtual studies this points at
  # the representative per-ancestry region so the downstream filter
  # recognizes it. For single-ancestry studies it's just the row's
  # coordinates.
  if (isTRUE(is_multi)) {
    region_stats[, region_key := paste0(CHR_var, ":", rep_start, "-", rep_stop)]
  } else {
    region_stats[, region_key := paste0(CHR_var, ":", BP_START_var, "-", BP_STOP_var)]
  }

  # Dot color: ancestry-coverage gradient in multi mode, alternating
  # chromosome colors otherwise.
  if (isTRUE(is_multi)) {
    region_stats[, chr_color := ANCESTRY_COVERAGE_COLORS[as.character(n_ancestries)]]
  } else {
    region_stats[, chr_color := ifelse(CHR_num %% 2 == 0, "#3498db", "#2c3e50")]
  }

  # Selection highlight
  sel <- selected_region_key
  if (is.null(sel) || length(sel) == 0) sel <- NA_character_
  region_stats[, is_selected := !is.na(sel) & region_key == sel]

  # Hover text
  if (isTRUE(is_multi)) {
    has_n10 <- "nearest_genes_10" %in% names(region_stats)
    region_stats[, hover := paste0(
      "<b>", gene, "</b><br>",
      if (has_n10) paste0("Nearby: ",
                          ifelse(is.na(nearest_genes_10), "", nearest_genes_10),
                          "<br>") else "",
      "chr", CHR_var, ":", format(BP_START_var, big.mark = ","),
      "-", format(BP_STOP_var, big.mark = ","), "<br>",
      "-log10(P): ", round(max_nlog10P, 1), "<br>",
      n_coloc, " colocalizations<br>",
      n_ancestries, "/4 ancestries: ", anc_label
    )]
  } else {
    region_stats[, hover := paste0(
      "<b>", gene, "</b><br>",
      "chr", CHR_var, ":", format(BP_START_var, big.mark = ","), "<br>",
      "-log10(P): ", round(max_nlog10P, 1), "<br>",
      n_coloc, " colocalizations"
    )]
  }

  region_stats
}

# Render the Manhattan plot from a prebuilt region_stats table. Returns
# a plotly object.
manhattan_plot <- function(region_stats, plot_title = "") {
  chr_ticks <- region_stats[, .(mid_x = mean(x_pos)), by = CHR_num]
  data.table::setorder(chr_ticks, CHR_num)

  p <- plotly::plot_ly() %>%
    plotly::add_trace(
      data = region_stats[is_selected == FALSE],
      x = ~x_pos, y = ~max_nlog10P,
      type = "scatter", mode = "markers",
      marker = list(color = ~chr_color,
                    size = ~pmin(5 + n_coloc / 10, 15),
                    opacity = 0.7),
      text = ~hover, hoverinfo = "text",
      customdata = ~region_key,
      name = "Regions"
    )

  if (any(region_stats$is_selected)) {
    p <- p %>% plotly::add_trace(
      data = region_stats[is_selected == TRUE],
      x = ~x_pos, y = ~max_nlog10P,
      type = "scatter", mode = "markers",
      marker = list(color = "red", size = 14, symbol = "diamond", opacity = 1),
      text = ~hover, hoverinfo = "text",
      customdata = ~region_key,
      name = "Selected"
    )
  }

  # Label top 20 regions by coloc count, PLUS the currently selected
  # region (even if not in top 20). Selected gets bold red label.
  top_labeled <- region_stats[order(-n_coloc, -max_nlog10P)][
    1:min(20, nrow(region_stats))]
  selected_row <- region_stats[is_selected == TRUE]
  if (nrow(selected_row) > 0 &&
      !selected_row$region_key[1] %in% top_labeled$region_key) {
    top_labeled <- rbind(top_labeled, selected_row)
  }
  top_labeled[, gene_short := gsub("\\(.*\\)", "", gene)]
  top_labeled[, gene_short := trimws(gsub("INTERGENIC: ", "", gene_short))]

  gene_annotations <- lapply(seq_len(nrow(top_labeled)), function(i) {
    is_sel <- isTRUE(top_labeled$is_selected[i])
    list(
      x = top_labeled$x_pos[i],
      y = top_labeled$max_nlog10P[i],
      text = if (is_sel) paste0("<b>", top_labeled$gene_short[i], "</b>")
             else top_labeled$gene_short[i],
      showarrow = FALSE,
      font = list(size = if (is_sel) 15 else 13,
                  color = if (is_sel) "#e74c3c" else "#333"),
      yshift = 14
    )
  })

  p %>% plotly::layout(
    title = list(text = plot_title,
                 font = list(size = 18, color = "#2c3e50"),
                 x = 0.5, xanchor = "center", y = 0.98),
    annotations = gene_annotations,
    xaxis = list(
      title = "", showgrid = FALSE,
      tickvals = chr_ticks$mid_x,
      ticktext = chr_ticks$CHR_num,
      tickfont = list(size = 13),
      range = c(min(region_stats$x_pos) * 0.98,
                max(region_stats$x_pos) * 1.02)
    ),
    yaxis = list(title = list(text = "-log10(P)", font = list(size = 14)),
                 tickfont = list(size = 12)),
    margin = list(t = 32, b = 25, l = 40, r = 5),
    showlegend = FALSE,
    hovermode = "closest"
  ) %>%
    plotly::config(displayModeBar = FALSE) %>%
    plotly::event_register("plotly_click")
}
