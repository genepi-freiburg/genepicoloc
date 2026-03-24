# Regional association plots (placeholder for future implementation)
# Will be populated when regional sumstats data is available

#' Plot regional association (LocusZoom-style)
#'
#' @param sumstats data.frame with CHR, POS, nlog10P, Name, rsID
#' @param title Plot title
#' @param lead_snp Name of lead SNP to highlight
#' @return plotly object
plot_regional_association <- function(sumstats, title = "", lead_snp = NULL) {
  # Placeholder - returns empty plot with message
  plot_ly() %>%
    layout(
      title = title,
      annotations = list(
        text = "Regional plots require sumstats data (save_sumstats = TRUE)",
        x = 0.5, y = 0.5, xref = "paper", yref = "paper",
        showarrow = FALSE, font = list(size = 14, color = "gray")
      )
    )
}
