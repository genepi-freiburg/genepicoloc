# Network visualization functions

#' Build region-centered network
#'
#' Central star node = genomic region, satellite nodes = colocalized traits.
#'
#' @param dt Filtered data.table with trait_name, PP.H4.abf, source_study, directionality
#' @param study_colors Named vector of hex colors per study
#' @return visNetwork object
build_region_network <- function(dt, study_colors) {

  if (nrow(dt) == 0) return(visNetwork(data.frame(), data.frame()))

  # Central node: region
  gene <- dt$nearest_gene_1[1]
  gene_clean <- gsub("\\(within\\)|INTERGENIC: ", "", gene)

  center_node <- data.frame(
    id = "center",
    label = gene_clean,
    shape = "star",
    size = 30,
    color = "#f1c40f",
    font.size = 16,
    title = paste0(
      "<b>", gene, "</b><br>",
      "chr", dt$CHR_var[1], ":", dt$BP_START_var[1], "-", dt$BP_STOP_var[1], "<br>",
      nrow(dt), " colocalizations"
    ),
    group = "Region",
    stringsAsFactors = FALSE
  )

  # Deduplicate: keep highest PP.H4 per trait_id
  dt <- dt[order(-PP.H4.abf)]
  dt <- dt[!duplicated(trait_id)]

  # Trait nodes
  trait_nodes <- data.frame(
    id = dt$trait_id,
    label = substr(dt$trait_name, 1, 30),  # truncate long names
    shape = "dot",
    size = 10 + 20 * dt$PP.H4.abf,
    color = vapply(dt$source_study, function(s) {
      if (s %in% names(study_colors)) study_colors[s] else "#999999"
    }, character(1)),
    font.size = 10,
    title = paste0(
      "<b>", dt$trait_name, "</b><br>",
      "Study: ", dt$source_study, "<br>",
      "PP.H4: ", round(dt$PP.H4.abf, 3), "<br>",
      "nsnps: ", dt$nsnps
    ),
    group = dt$source_study,
    stringsAsFactors = FALSE
  )

  nodes <- rbind(center_node, trait_nodes)

  # Edges with directional coloring
  edge_color <- ifelse(
    is.na(dt$directionality) | dt$directionality == 0, "#848484",
    ifelse(dt$directionality > 0, "#e67e22", "#3498db")
  )

  edges <- data.frame(
    from = "center",
    to = dt$trait_id,
    color = edge_color,
    width = 3,
    title = paste0("PP.H4: ", round(dt$PP.H4.abf, 3)),
    stringsAsFactors = FALSE
  )

  visNetwork(nodes, edges) %>%
    visPhysics(
      enabled = nrow(nodes) > 20,
      stabilization = list(enabled = TRUE, iterations = 200)
    ) %>%
    visInteraction(
      hover = TRUE,
      tooltipDelay = 100,
      navigationButtons = TRUE
    ) %>%
    visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1),
      nodesIdSelection = FALSE
    ) %>%
    visEvents(select = "function(nodes) {
      Shiny.setInputValue('region_network_selected', nodes.nodes[0]);
    }")
}


#' Build trait-centered network
#'
#' Central diamond node = selected trait, satellite nodes = genomic regions.
#'
#' @param dt data.table filtered to one trait across all regions
#' @param trait_name Selected trait name
#' @param study_colors Named vector of hex colors
#' @return visNetwork object
build_trait_network <- function(dt, trait_name, study_colors) {

  if (nrow(dt) == 0) return(visNetwork(data.frame(), data.frame()))

  # Central node: trait
  center_node <- data.frame(
    id = "center",
    label = substr(trait_name, 1, 40),
    shape = "diamond",
    size = 30,
    color = "#9b59b6",
    font.size = 14,
    title = paste0("<b>", trait_name, "</b><br>", nrow(dt), " regions"),
    stringsAsFactors = FALSE
  )

  # Region nodes
  gene_labels <- gsub("\\(within\\)|INTERGENIC: ", "", dt$nearest_gene_1)
  region_nodes <- data.frame(
    id = dt$region_id,
    label = gene_labels,
    shape = "dot",
    size = 10 + 20 * dt$PP.H4.abf,
    color = "#4DBBD5",
    font.size = 10,
    title = paste0(
      "<b>", dt$nearest_gene_1, "</b><br>",
      dt$region_id, "<br>",
      "PP.H4: ", round(dt$PP.H4.abf, 3)
    ),
    stringsAsFactors = FALSE
  )

  # Deduplicate regions (same region might appear if multiple sub-regions)
  region_nodes <- region_nodes[!duplicated(region_nodes$id), ]

  nodes <- rbind(center_node, region_nodes)

  edges <- data.frame(
    from = "center",
    to = dt$region_id,
    color = ifelse(is.na(dt$directionality) | dt$directionality == 0, "#848484",
                   ifelse(dt$directionality > 0, "#e67e22", "#3498db")),
    width = 3,
    title = paste0("PP.H4: ", round(dt$PP.H4.abf, 3)),
    stringsAsFactors = FALSE
  )

  visNetwork(nodes, edges) %>%
    visPhysics(
      enabled = nrow(nodes) > 20,
      stabilization = list(enabled = TRUE, iterations = 200)
    ) %>%
    visInteraction(hover = TRUE, tooltipDelay = 100, navigationButtons = TRUE) %>%
    visOptions(highlightNearest = list(enabled = TRUE, degree = 1)) %>%
    visEvents(select = "function(nodes) {
      Shiny.setInputValue('trait_network_selected', nodes.nodes[0]);
    }")
}
