# GWAS Colocalization Network Viewer with Prioritized Genes
# Shiny application for visualizing GWAS colocalization results
#
# Loading order (managed by Shiny):
#   1. R/*.R   - auto-sourced by Shiny BEFORE app.R (shiny.autoload.r=TRUE).
#                Use data.table::fread etc. since packages aren't attached yet.
#   2. app.R   - this file. library() calls below attach packages for runtime use.
#
# Do NOT source("R/*.R") here - auto-loader handles it (and duplicate sourcing
# into different environments creates scoping bugs, see rstudio/shiny#2547).

library(shiny)
library(visNetwork)
library(data.table)
library(plotly)

  # UI ----
  ui <- fluidPage(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "landing.css"),
      # JavaScript for PNG export
      tags$script(HTML("
        // Region View network export
        Shiny.addCustomMessageHandler('exportNetwork', function(message) {
          var container = document.getElementById('network');
          if (container) {
            var canvas = container.querySelector('canvas');
            if (canvas) {
              var link = document.createElement('a');
              link.download = message.filename;
              link.href = canvas.toDataURL('image/png');
              link.click();
            } else {
              alert('Network canvas not found. Please wait for the network to render.');
            }
          } else {
            alert('Network container not found.');
          }
        });
      "))
    ),
    
    tabsetPanel(
      id = "main_tabs",
      tabPanel("Atlas",
               # Landing page with category cards
               div(class = "landing-hero",
                   h1("Kidney Genomics Colocalization Atlas"),
                   p(class = "subtitle", "Explore colocalization results across kidney disease, imaging, and metabolomics datasets")
               ),

               # Category cards
               uiOutput("landing_categories"),

               div(class = "landing-legend",
                   p(tags$span(style = "border-left: 3px solid #27ae60; padding-left: 6px; margin-right: 12px;",
                               "Regional plots available"),
                     tags$span(style = "opacity: 0.5; margin-right: 12px;", "Coming soon"))
               )
      ),
      tabPanel("Region View",
               # Add notification area for study loading
               fluidRow(
                 column(12,
                        uiOutput("study_notification")
                 )
               ),

               # Top row: Search (left) + mini Manhattan (right)
               fluidRow(
                 column(
                   width = 2,
                   div(style = "padding: 8px;",
                     selectizeInput("selected_region",
                       "Search by gene:",
                       choices = NULL,
                       selected = NULL,
                       options = list(
                         placeholder = "Type gene name...",
                         closeAfterSelect = TRUE
                       ),
                       width = "100%")
                   )
                 ),
                 column(10, uiOutput("mini_manhattan_ui"))
               ),

               fluidRow(
                 # Left sidebar panel
                 column(
                   width = 2,
                   wellPanel(
                     style = "padding: 10px;",
                     # Current study indicator
                     uiOutput("current_study_info"),

                     # Back to study selection
                     actionButton("back_to_home", "Different study",
                                  class = "btn-xs btn-info btn-block",
                                  icon = icon("arrow-left"),
                                  style = "margin-bottom: 8px;"),

                     hr(style = "margin: 8px 0;"),

                     # Colocalized traits - category cards (moved from main panel)
                     h6("Colocalized traits", style = "margin: 0 0 4px 0; color: #555;"),
                     div(style = "font-size: 10px; color: #999; margin-bottom: 4px;",
                       "Click to explore."),
                     uiOutput("conv_category_cards"),

                     hr(style = "margin: 8px 0;"),

                     # Filter options
                     # Filters (collapsible, collapsed by default to save space)
                     tags$details(
                       style = "margin-bottom: 6px;",
                       tags$summary(style = "cursor: pointer; font-size: 11px; color: #555;",
                                    "Filters"),
                       div(style = "margin-top: 6px;",
                         sliderInput("min_pph4", "Min PP.H4:",
                                     min = 0.8, max = 1, value = 0.8, step = 0.01),
                         sliderInput("min_nlog10P", "Min -log10(P):",
                                     min = 7.308, max = 100, value = 7.308, step = 0.5),
                         sliderInput("max_traits", "Max traits (network):",
                                     min = 10, max = 500, value = 50, step = 10)
                       )
                     ),

                     # Include Studies (collapsed by default, dynamic)
                     tags$details(
                       style = "margin-bottom: 6px;",
                       tags$summary(
                         style = "cursor: pointer; font-size: 11px; color: #555;",
                         "Include studies"
                       ),
                       div(style = "margin-top: 6px;",
                         fluidRow(
                           column(6, actionButton("select_all", "All",
                                                  class = "btn-xs btn-block",
                                                  style = "margin-bottom: 3px;")),
                           column(6, actionButton("select_none", "None",
                                                  class = "btn-xs btn-block",
                                                  style = "margin-bottom: 3px;"))
                         ),
                         uiOutput("study_selector")
                       )
                     ),

                     # Display / download
                     tags$details(
                       style = "margin-bottom: 6px;",
                       tags$summary(style = "cursor: pointer; font-size: 11px; color: #555;",
                                    "More"),
                       div(style = "margin-top: 6px;",
                         checkboxInput("physics", "Enable physics simulation", value = FALSE),
                         downloadButton("download_data", "Download filtered",
                                        class = "btn-xs btn-block")
                       )
                     ),

                     # Hidden dummy input to preserve reactivity for legacy coord selector
                     tags$div(style = "display: none;",
                       selectizeInput("selected_region_coord", NULL,
                                      choices = NULL, selected = NULL)),

                     # Customize for Export panel (conditional based on config)
                     if (SHOW_EXPORT_CUSTOMIZATION) {
                       tagList(
                         hr(),
                         tags$details(open = "open",
                           tags$summary(style = "cursor: pointer; font-weight: bold;",
                                       "Customize for Export"),
                           div(style = "margin-top: 8px;",
                             # Node Selection
                             h6("Node Selection:", style = "margin-bottom: 5px;"),
                             fluidRow(
                               column(4, actionButton("nodes_all", "All", class = "btn-xs btn-block")),
                               column(4, actionButton("nodes_none", "None", class = "btn-xs btn-block")),
                               column(4, actionButton("nodes_top10", "Top 10", class = "btn-xs btn-block"))
                             ),
                             br(),
                             selectizeInput("visible_nodes",
                                           label = NULL,
                                           choices = NULL,
                                           multiple = TRUE,
                                           options = list(
                                             plugins = list('remove_button'),
                                             placeholder = 'Select nodes to display...'
                                           )),

                             hr(style = "margin: 10px 0;"),

                             # Study colors
                             h6("Study Colors:", style = "margin-bottom: 5px;"),
                             uiOutput("study_color_pickers"),

                             hr(style = "margin: 10px 0;"),

                             # Edit selected node label
                             h6("Edit Node Label:", style = "margin-bottom: 5px;"),
                             textInput("custom_node_label", label = NULL,
                                      placeholder = "Click a node to edit its label"),
                             actionButton("apply_label", "Apply Label",
                                         class = "btn-xs btn-success btn-block")
                           )
                         )
                       )
                     }
                   )  # End wellPanel
                 ),  # End left column

                 # Main panel (network + regional plot tabs)
                 column(
                   width = 10,
                   tabsetPanel(
                     id = "region_view_tabs",
                     type = "tabs",

                     # Tab 0: Convergence (region-centric multi-omics view)
                     tabPanel(
                       "Convergence",
                       # Panel 1: trait of interest regional plot
                       h6("Trait of interest", style = "margin: 8px 0 2px 0; color: #555;"),
                       plotlyOutput("conv_base_plot", height = "200px"),
                       # Panel 1b: Selected trait RAP (appears below base when user clicks a trait)
                       uiOutput("conv_trait_plot_header"),
                       conditionalPanel(
                         condition = "output.conv_trait_has_data",
                         plotlyOutput("conv_trait_plot", height = "200px")
                       ),
                       # Panel 2: gene track + gene info panel (click a gene)
                       h6("Genes in region", style = "margin: 8px 0 2px 0; color: #555;"),
                       div(style = "font-size: 10px; color: #999; margin-bottom: 2px;",
                         "Click a gene for details."),
                       plotlyOutput("conv_gene_track", height = "100px"),
                       uiOutput("gene_info_panel"),
                       # Panel 3: trait tiles for selected category
                       uiOutput("conv_drilldown_header"),
                       uiOutput("conv_trait_tiles")
                     ),

                     # Tab 1: Network visualization
                     tabPanel(
                       "Network",
                       # Export buttons - compact, right-aligned
                       div(style = "text-align: right; margin: 4px 0;",
                         actionButton("export_png", "PNG",
                                     icon = icon("image"),
                                     class = "btn-xs btn-default",
                                     style = "margin-right: 3px;"),
                         downloadButton("export_html", "HTML",
                                       class = "btn-xs btn-default")
                       ),
                       # Network visualization
                       visNetworkOutput("network", height = "650px"),
                       h5("Selected Node Details", style = "margin: 4px 0;"),
                       verbatimTextOutput("node_info")
                     ),  # End Network tab

                     # Regional Plot tab removed - functionality moved into
                     # the Convergence tab (click a trait tile, click a gene
                     # on the gene track for details).
                   )  # End tabsetPanel
                 )  # End main column
               )  # End fluidRow
      ),
      tabPanel("Documentation",
               tags$iframe(
                 src = "documentation.html",
                 style = "width: 100%; height: calc(100vh - 80px); border: none;"
               )
      )
    )
  )
  
  server <- function(input, output, session) {

    # === Landing Page ===

    # Category definitions for the landing page
    atlas_categories <- list(
      list(
        id = "kidney_disease",
        icon = "\U0001F9EC",
        title = "Kidney Disease Traits",
        description = "CKDGen Round 4 GWAS - kidney function and disease phenotypes colocalized against 11 molecular and clinical datasets",
        traits = list(
          list(id = "eGFR", label = "eGFR (creatinine)", desc = "Estimated glomerular filtration rate"),
          list(id = "BUN", label = "Blood Urea Nitrogen", desc = "Kidney filtration marker"),
          list(id = "UACR", label = "UACR", desc = "Urinary albumin-to-creatinine ratio"),
          list(id = "urate", label = "Serum Urate", desc = "Urate levels"),
          list(id = "gout", label = "Gout", desc = "Inflammatory arthritis"),
          list(id = "MA", label = "Microalbuminuria", desc = "Early kidney damage marker")
        )
      ),
      list(
        id = "kidney_mri",
        icon = "\U0001F9F2",
        title = "Kidney MRI Volumes",
        description = "UK Biobank kidney MRI - structural imaging phenotypes (BSA-adjusted)",
        traits = list(
          list(id = "MRI_tkv", label = "Total Kidney Volume", desc = "BSA-adjusted"),
          list(id = "MRI_cortex", label = "Cortex Volume", desc = "BSA-adjusted"),
          list(id = "MRI_medulla", label = "Medulla Volume", desc = "BSA-adjusted"),
          list(id = "MRI_hilus", label = "Hilus Volume", desc = "BSA-adjusted")
        )
      ),
      list(
        id = "metabolomics",
        icon = "\U0001F9EA",
        title = "Urine Metabolomics",
        description = paste(
          "GCKD urine metabolome GWAS (Schlosser et al., Nat Genet 2023) -",
          "1,409 urine metabolites measured by Metabolon, colocalized against",
          "Tier 1 datasets."
        ),
        # Traits are populated dynamically below from auto-discovered uMet
        # files so we don't have to enumerate 1,409 metabolites by hand.
        traits = list(),
        autopopulate = "GCKD_uMet"
      ),
      list(
        id = "mvp_kidney",
        icon = "\U0001F30D",
        title = "MVP Kidney (multi-ancestry)",
        description = paste(
          "Million Veteran Program kidney traits colocalized across five",
          "ancestries (AFR, AMR, EAS, EUR, META) against ancestry-matched",
          "MVP_R4 PheWAS. Each trait loads all available ancestries at once",
          "and overlays them in the Manhattan and regional plots."
        ),
        # Traits populated from virtual multi-ancestry studies.
        traits = list(),
        autopopulate_virtual = TRUE
      )
    )

    # Build trait list for one landing card. If the card declares
    # `autopopulate = "<atlas_category>"`, expand its (empty) traits list to
    # all studies discovered under that atlas folder, labeled via the bundled
    # annotation DB. Used for urine metabolomics where hand-listing 1,409
    # metabolites is impractical.
    expand_card_traits <- function(cat) {
      if (length(cat$traits) > 0) return(cat$traits)

      # Virtual multi-ancestry card: one entry per virtual study id.
      if (isTRUE(cat$autopopulate_virtual)) {
        if (length(DEFAULT_VIRTUAL_STUDIES) == 0) return(list())
        return(lapply(names(DEFAULT_VIRTUAL_STUDIES), function(vid) {
          v <- DEFAULT_VIRTUAL_STUDIES[[vid]]
          # Strip MVPkid_ prefix for display; keep the virtual id as $id.
          lbl <- sub("^MVPkid_", "", vid)
          list(id = vid, label = lbl,
               desc = paste0(length(v$ancestries), " ancestries: ",
                             paste(v$ancestries, collapse = ", ")))
        }))
      }

      if (is.null(cat$autopopulate)) return(cat$traits)

      cats_attr <- attr(DEFAULT_AVAILABLE_STUDIES, "categories")
      if (is.null(cats_attr)) return(list())
      ids <- names(cats_attr)[unlist(cats_attr) == cat$autopopulate]

      # Drop trait_ids whose coloc RDS has zero rows. Metabolites with no
      # high-confidence colocs would render as empty Region Views, so we
      # hide them from the landing card entirely. This is a performance-
      # cheap filter (readRDS on ~322 KB slim files) done once at UI build.
      ids <- Filter(function(id) {
        path <- DEFAULT_AVAILABLE_STUDIES[[id]]
        if (is.null(path) || !file.exists(path)) return(FALSE)
        tryCatch(nrow(readRDS(path)) > 0, error = function(e) FALSE)
      }, ids)

      label_fn <- if (cat$autopopulate == "GCKD_uMet") umet_label else identity
      lapply(ids, function(id) {
        list(id = id, label = label_fn(id), desc = id)
      })
    }

    # Is this trait id resolvable? Accepts real study ids (present in
    # DEFAULT_AVAILABLE_STUDIES) OR virtual multi-ancestry ids.
    is_known_study <- function(id) {
      id %in% names(DEFAULT_AVAILABLE_STUDIES) ||
        id %in% names(DEFAULT_VIRTUAL_STUDIES)
    }

    # Build the region-bundle lookup key for a trait row. For a virtual
    # multi-ancestry study the region bundle is per-ancestry, so the key
    # must be reassembled as "<source_study>_<ancestry>__<basename>" even
    # though coloc_data() has unified source_study (e.g. "MVP_R4").
    make_bundle_key <- function(source_study, sumstats_2_file, ancestry = NULL) {
      # Vector-safe: works for scalar or data.table-column inputs.
      if (is.null(ancestry)) {
        return(paste0(source_study, "__", basename(sumstats_2_file)))
      }
      is_virt <- !is.null(current_study()) && is_virtual_study(current_study())
      base_key <- paste0(source_study, "__", basename(sumstats_2_file))
      anc_key  <- paste0(source_study, "_", ancestry, "__", basename(sumstats_2_file))
      use_anc  <- is_virt & !is.na(ancestry) & nzchar(ancestry)
      ifelse(use_anc, anc_key, base_key)
    }

    # Consensus trait key for a virtual multi-ancestry study: collapses all
    # per-ancestry variants of the same trait into one stable key.
    # Turns "MVP_R4_EUR__MVP_R4.1000G_AGR.A1C_Max_INT.EUR.GIA.dbGaP.txt.gz"
    # into "MVP_R4__MVP_R4.1000G_AGR.A1C_Max_INT.GIA.dbGaP.txt.gz".
    ANC_CODES <- c("AFR","AMR","EAS","EUR","META")
    make_consensus_trait_key <- function(bundle_key) {
      anc_alt <- paste(ANC_CODES, collapse = "|")
      # Strip "_<ANC>" that immediately precedes "__" in the study prefix.
      k <- sub(paste0("_(", anc_alt, ")__"), "__", bundle_key)
      # Also drop ".<ANC>." segment from the filename portion so keys match
      # across ancestries whose filenames differ only by that token.
      k <- sub(paste0("\\.(", anc_alt, ")\\."), ".", k)
      k
    }

    # Given a consensus trait key and a per-ancestry region bundle, try to
    # locate the matching ancestry-specific bundle key inside that bundle.
    resolve_consensus_in_bundle <- function(consensus_key, bundle_keys, anc) {
      # Rebuild expected ancestry-specific candidate keys:
      #  1. Insert "_<anc>" after the study prefix
      #  2. Insert ".<anc>" before ".GIA" (common MVP naming)
      #  3. Or more generally: if consensus key split at "__" -> prefix,fn,
      #     candidates are <prefix>_<anc>__<fn with .<anc>. spliced in>.
      parts <- strsplit(consensus_key, "__", fixed = TRUE)[[1]]
      if (length(parts) != 2) return(NULL)
      prefix <- parts[1]; fn <- parts[2]

      cand <- c(
        paste0(prefix, "_", anc, "__", sub("\\.GIA", paste0(".", anc, ".GIA"), fn)),
        paste0(prefix, "_", anc, "__", sub("^(MVP_R4\\.[^.]+\\.[^.]+)\\.", paste0("\\1.", anc, "."), fn))
      )
      hit <- intersect(cand, bundle_keys)
      if (length(hit) > 0) return(hit[1])

      # Fallback: fuzzy match - any bundle key that starts with
      # "<prefix>_<anc>__" and reduces (via make_consensus_trait_key) to
      # the same consensus_key.
      pref <- paste0(prefix, "_", anc, "__")
      cands <- bundle_keys[startsWith(bundle_keys, pref)]
      if (length(cands) == 0) return(NULL)
      matches <- cands[vapply(cands, function(k) identical(make_consensus_trait_key(k), consensus_key), logical(1))]
      if (length(matches) > 0) matches[1] else NULL
    }

    output$landing_categories <- renderUI({
      div(class = "category-grid",
        lapply(atlas_categories, function(cat) {
          cat$traits <- expand_card_traits(cat)
          available <- sapply(cat$traits, function(t) {
            !isTRUE(t$disabled) && is_known_study(t$id)
          })
          n_available <- sum(available)
          n_total <- length(cat$traits)
          stats_text <- if (n_available > 0) {
            paste0(n_available, " traits available")
          } else {
            "Coming soon"
          }

          div(class = "category-card",
            div(class = "card-icon", cat$icon),
            h3(cat$title),
            p(class = "card-description", cat$description),
            p(class = "card-stats", stats_text),
            div(class = "trait-list",
              lapply(cat$traits, function(trait) {
                is_virtual <- trait$id %in% names(DEFAULT_VIRTUAL_STUDIES)
                is_disabled <- isTRUE(trait$disabled) || !is_known_study(trait$id)
                # Check regional data across all possible layouts
                cats <- attr(DEFAULT_AVAILABLE_STUDIES, "categories")
                cat_name <- if (!is.null(cats)) cats[[trait$id]] else NULL
                has_regional <- if (is_virtual) {
                  # Virtual study: any ancestry with a regional dir counts.
                  vdirs <- DEFAULT_VIRTUAL_STUDIES[[trait$id]]$regional_dirs
                  any(!is.na(unlist(vdirs)))
                } else {
                  !is_disabled && (
                    (!is.null(cat_name) && dir.exists(file.path(DATA_PATH, cat_name, "regional", trait$id))) ||
                    dir.exists(file.path(DATA_PATH, "regional", trait$id)) ||
                    dir.exists(file.path(DATA_PATH, "regional_plots", trait$id))
                  )
                }
                css_class <- paste("trait-btn",
                  if (is_disabled) "disabled",
                  if (has_regional) "has-regional"
                )
                if (is_disabled) {
                  tags$span(class = css_class, title = trait$desc, trait$label)
                } else {
                  tags$span(
                    class = css_class,
                    title = trait$desc,
                    onclick = paste0("Shiny.setInputValue('selected_study', '", trait$id, "', {priority: 'event'});"),
                    trait$label
                  )
                }
              })
            )
          )
        })
      )
    })

    # === Reactive Values ===

    # Reactive value to store the current loaded study
    current_study <- reactiveVal(NULL)

    # === Event Handlers ===
    
    # Handle study selection from home page
    observeEvent(input$selected_study, {
      req(input$selected_study)
      
      # Resolve file path (supports flat and nested layouts)
      file_path <- resolve_annot_path(input$selected_study)
      
      # Check if file exists
      if (file.exists(file_path)) {
        # Store current study (now storing the name, not the folder)
        current_study(input$selected_study)

        # Switch to Region View tab
        updateTabsetPanel(session, "main_tabs", selected = "Region View")
        
        # Show loading notification
        showNotification(
          paste("Loading study:", input$selected_study),
          type = "message",
          duration = 3
        )
      } else {
        showNotification(
          paste("Error: Could not find data file for", input$selected_study),
          type = "error",
          duration = 5
        )
      }
    })

    # Handle back to home button
    observeEvent(input$back_to_home, {
      updateTabsetPanel(session, "main_tabs", selected = "Atlas")
    })

    # Select all studies (dynamic)
    observeEvent(input$select_all, {
      lapply(available_study_names(), function(study) {
        updateCheckboxInput(session, paste0("study_", study), value = TRUE)
      })
    })

    # Select no studies (dynamic)
    observeEvent(input$select_none, {
      lapply(available_study_names(), function(study) {
        updateCheckboxInput(session, paste0("study_", study), value = FALSE)
      })
    })

    # Select studies by group
    observeEvent(input$select_phenotypes, {
      lapply(study_categories$Phenotypes, function(study) {
        if (study %in% names(study_colors)) {
          updateCheckboxInput(session, paste0("study_", study), value = TRUE)
        }
      })
    })

    observeEvent(input$select_pqtl, {
      lapply(study_categories$pQTL, function(study) {
        if (study %in% names(study_colors)) {
          updateCheckboxInput(session, paste0("study_", study), value = TRUE)
        }
      })
    })

    observeEvent(input$select_eqtl, {
      lapply(study_categories$eQTL, function(study) {
        if (study %in% names(study_colors)) {
          updateCheckboxInput(session, paste0("study_", study), value = TRUE)
        }
      })
    })

    observeEvent(input$select_mqtl, {
      lapply(study_categories$mQTL, function(study) {
        if (study %in% names(study_colors)) {
          updateCheckboxInput(session, paste0("study_", study), value = TRUE)
        }
      })
    })

    # === Node Selection for Export ===

    # Store available nodes (all nodes from network_data)
    available_nodes <- reactiveVal(NULL)

    # Update node selector choices when network data changes
    observe({
      req(network_data())
      net_data <- network_data()

      if (nrow(net_data$nodes) > 0) {
        # Create node choices: id -> label (group)
        node_choices <- setNames(
          net_data$nodes$id,
          paste0(net_data$nodes$label, " (", net_data$nodes$group, ")")
        )
        # Store for quick selection buttons
        available_nodes(net_data$nodes)

        # Update selectize with all nodes selected by default
        updateSelectizeInput(session, "visible_nodes",
                            choices = node_choices,
                            selected = net_data$nodes$id)

        # Auto-enable physics simulation for large networks (>150 nodes)
        if (nrow(net_data$nodes) > 150 && !input$physics) {
          updateCheckboxInput(session, "physics", value = TRUE)
        }
      }
    })

    # Node selection: All
    observeEvent(input$nodes_all, {
      req(available_nodes())
      updateSelectizeInput(session, "visible_nodes",
                          selected = available_nodes()$id)
    })

    # Node selection: None
    observeEvent(input$nodes_none, {
      updateSelectizeInput(session, "visible_nodes", selected = character(0))
    })

    # Node selection: Top 10 by PP.H4
    observeEvent(input$nodes_top10, {
      req(filtered_data(), available_nodes())
      dt <- filtered_data()

      # Get top 10 by PP.H4 (excluding the central region node)
      top_10 <- head(dt[order(-PP.H4.abf)], 10)

      # Build matching node IDs
      top_node_ids <- paste0(top_10$source_study, "_", seq_len(nrow(top_10)))

      # Always include the region node
      selected_ids <- c("region", top_node_ids)

      updateSelectizeInput(session, "visible_nodes", selected = selected_ids)
    })

    # === Customize for Export (Colors & Labels) ===

    # Reactive values for custom colors and labels
    custom_colors <- reactiveVal(study_colors)  # Start with default colors
    custom_labels <- reactiveVal(list())  # node_id -> custom_label

    # Render color pickers for each study in the current data
    output$study_color_pickers <- renderUI({
      req(filtered_data())
      dt <- filtered_data()

      # Get unique studies in current data
      studies_in_data <- unique(dt$source_study)
      studies_in_data <- studies_in_data[studies_in_data %in% names(study_colors)]

      if (length(studies_in_data) == 0) {
        return(p("No studies to customize", style = "color: #999; font-size: 11px;"))
      }

      # Create color picker for each study - use study_colors directly for initial values
      tagList(
        lapply(studies_in_data, function(study) {
          display_name <- get_study_display_name(study)
          # Use study_colors for initial value (unname to get plain hex string)
          initial_color <- unname(study_colors[study])
          div(style = "display: flex; align-items: center; margin-bottom: 3px;",
            textInput(
              paste0("color_", study),
              label = NULL,
              value = initial_color
            ),
            span(display_name, style = "font-size: 11px; margin-left: 5px;")
          )
        })
      )
    })

    # Observe color picker changes and update custom_colors
    observe({
      req(filtered_data())
      dt <- filtered_data()
      studies_in_data <- unique(dt$source_study)
      studies_in_data <- studies_in_data[studies_in_data %in% names(study_colors)]

      current_colors <- custom_colors()

      for (study in studies_in_data) {
        color_input <- input[[paste0("color_", study)]]
        if (!is.null(color_input)) {
          current_colors[[study]] <- color_input
        }
      }

      custom_colors(current_colors)
    })

    # When a node is clicked, populate the label input with current label
    observe({
      req(input$selected_node)
      net_data <- network_data()

      node <- net_data$nodes[net_data$nodes$id == input$selected_node, ]
      if (nrow(node) > 0) {
        # Check if there's a custom label, otherwise use the original
        labels <- custom_labels()
        current_label <- if (input$selected_node %in% names(labels)) {
          labels[[input$selected_node]]
        } else {
          node$label
        }
        updateTextInput(session, "custom_node_label", value = current_label)
      }
    })

    # Apply custom label when button is clicked
    observeEvent(input$apply_label, {
      req(input$selected_node, input$custom_node_label)

      labels <- custom_labels()
      labels[[input$selected_node]] <- input$custom_node_label
      custom_labels(labels)

      # Show confirmation
      showNotification(
        paste0("Label updated for: ", input$selected_node),
        type = "message",
        duration = 2
      )
    })

    # PNG Export via JavaScript
    observeEvent(input$export_png, {
      # Get current region name for filename
      region <- current_region()
      gene <- if (!is.null(filtered_data()) && nrow(filtered_data()) > 0) {
        fd <- filtered_data()
        if (!is.na(fd$Prioritized_Gene[1])) fd$Prioritized_Gene[1] else "network"
      } else {
        "network"
      }
      filename <- paste0(gene, "_", gsub(":", "-", region), ".png")

      # Trigger JavaScript export
      session$sendCustomMessage("exportNetwork", list(filename = filename))
    })

    # HTML Export (interactive) - Region View
    output$export_html <- downloadHandler(
      filename = function() {
        region <- current_region()
        gene <- if (!is.null(filtered_data()) && nrow(filtered_data()) > 0) {
          fd <- filtered_data()
          if (!is.na(fd$Prioritized_Gene[1])) fd$Prioritized_Gene[1] else "network"
        } else {
          "network"
        }
        paste0(gene, "_", gsub(":", "-", region), ".html")
      },
      content = function(file) {
        net_data <- network_data()

        # Filter nodes based on selection
        visible <- input$visible_nodes
        if (!is.null(visible) && length(visible) > 0) {
          nodes_filtered <- net_data$nodes[net_data$nodes$id %in% visible, ]
          edges_filtered <- net_data$edges[
            net_data$edges$from %in% visible & net_data$edges$to %in% visible,
          ]
        } else {
          nodes_filtered <- net_data$nodes
          edges_filtered <- net_data$edges
        }

        # Create visNetwork and save
        network <- visNetwork(nodes_filtered, edges_filtered) %>%
          visNodes(
            font = list(color = "#000000", size = 14, face = "arial")
          ) %>%
          visOptions(
            highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)
          ) %>%
          visPhysics(enabled = FALSE) %>%
          visLayout(randomSeed = 123) %>%
          visInteraction(
            dragNodes = TRUE,
            dragView = TRUE,
            zoomView = TRUE
          )

        visSave(network, file)
      }
    )

    # Update region selectors
    # Map from bare region_key -> gene-prefixed selectize value
    region_to_value_map <- reactiveVal(character(0))

    observe({
      req(regions_data())
      regions <- regions_data()
      message("[gene_search] regions: ", nrow(regions),
              " gene_annotation null?: ", is.null(gene_annotation),
              " nrow: ", if (is.null(gene_annotation)) 0 else nrow(gene_annotation))

      # No regions for the current study (e.g. a uMet metabolite with zero
      # colocs passing filters): clear the selectizes and bail. Without this
      # guard the fallback below builds a degenerate length-1 choice with an
      # NA name and updateSelectizeInput crashes inside data.frame().
      if (nrow(regions) == 0) {
        region_to_value_map(character(0))
        updateSelectizeInput(session, "selected_region",
                             choices = character(0), server = TRUE)
        updateSelectizeInput(session, "selected_region_coord",
                             choices = character(0), server = TRUE)
        return()
      }

      # Build gene -> region mapping by overlapping gene annotation with
      # region coordinates. Any gene overlapping a region becomes searchable
      # and jumps to that region.
      gene_choices <- character(0)
      if (!is.null(gene_annotation) && nrow(gene_annotation) > 0) {
        regions_for_search <- as.data.table(regions)
        regions_for_search[, chr_clean := sub("^chr", "", as.character(CHR_var))]
        regions_for_search[, region_key := paste0(
          CHR_var, ":", BP_START_var, "-", BP_STOP_var
        )]
        # Preserve original region coords so non-equi join columns don't
        # overwrite them (data.table uses BP_START_var/BP_STOP_var as
        # comparison RHS columns in the join).
        regions_for_search[, orig_start := BP_START_var]
        regions_for_search[, orig_stop := BP_STOP_var]
        ga <- as.data.table(gene_annotation)
        ga[, chr_clean := sub("^chr", "", as.character(chr))]

        # Non-equi join: find genes overlapping each region
        matched <- regions_for_search[
          ga,
          on = .(chr_clean, BP_START_var <= end, BP_STOP_var >= start),
          nomatch = 0,
          .(gene_name = i.gene_name, region_key, CHR_var,
            BP_START_var = orig_start, BP_STOP_var = orig_stop)
        ]

        message("[gene_search] matched: ", nrow(matched),
                " UMOD: ", nrow(matched[gene_name == "UMOD"]))

        if (nrow(matched) > 0) {
          # If a gene maps to multiple regions, disambiguate with suffix
          matched[, gene_occurrences := .N, by = gene_name]
          matched[, occurrence_idx := seq_len(.N), by = gene_name]
          matched[, label := ifelse(
            gene_occurrences > 1,
            paste0(gene_name, " (chr", CHR_var, ":",
                   format(BP_START_var, big.mark = ","), ")"),
            gene_name
          )]
          setorder(matched, gene_name, occurrence_idx)
          matched <- matched[!duplicated(label)]
          # IMPORTANT: selectize deduplicates by VALUE, so we need a unique
          # value per gene. Encode gene name into the value, then parse back
          # in an observer to set the real region.
          matched[, value := paste0(gene_name, "::", region_key)]
          gene_choices <- setNames(matched$value, matched$label)
          # Cache a region_key -> gene_prefixed_value map so the Manhattan
          # click handler can resolve a bare region_key to a valid selectize
          # choice (pick the first gene in that region).
          first_value_per_region <- matched[, .(first_value = value[1]),
                                             by = region_key]
          region_to_value_map(setNames(first_value_per_region$first_value,
                                         first_value_per_region$region_key))
        }
      }

      # Fallback: if gene annotation is missing, use nearest_gene_1 as before
      if (length(gene_choices) == 0) {
        gene_display <- ifelse(!is.na(regions$Prioritized_Gene),
                              regions$Prioritized_Gene,
                              regions$nearest_gene_1)
        gene_clean <- gsub("\\(.*\\)", "", gene_display)
        gene_clean <- trimws(gsub("INTERGENIC:\\s*", "", gene_clean))
        gene_clean <- sub("\\s*\\([^)]*\\)$", "", gene_clean)
        gene_choices <- setNames(
          paste0(regions$CHR_var, ":", regions$BP_START_var, "-", regions$BP_STOP_var),
          gene_clean
        )
      }

      message("[gene_search] final gene_choices length: ", length(gene_choices))
      updateSelectizeInput(session, "selected_region", choices = gene_choices,
                           server = TRUE)
      # Legacy coord selector - same choices
      updateSelectizeInput(session, "selected_region_coord",
                           choices = gene_choices, server = TRUE)
    })
    
    # Synchronize selectors (legacy coord selector is hidden but kept in sync)
    .is_valid_sel <- function(x) {
      !is.null(x) && length(x) > 0 && !is.na(x[1]) && nzchar(x[1])
    }
    observeEvent(input$selected_region, {
      if (.is_valid_sel(input$selected_region)) {
        updateSelectizeInput(session, "selected_region_coord",
                             selected = input$selected_region)
      }
    }, ignoreNULL = FALSE)

    observeEvent(input$selected_region_coord, {
      if (.is_valid_sel(input$selected_region_coord)) {
        updateSelectizeInput(session, "selected_region",
                             selected = input$selected_region_coord)
      }
    }, ignoreNULL = FALSE)
    
    # === Mini Manhattan Plot ===

    output$mini_manhattan_ui <- renderUI({
      req(regions_data(), coloc_data())
      plotlyOutput("mini_manhattan", height = "180px")
    })

    output$mini_manhattan <- renderPlotly({
      req(regions_data(), coloc_data())
      regions <- regions_data()
      dt <- coloc_data()

      # Multi-ancestry virtual study: cluster overlapping per-ancestry
      # regions into consensus regions so the Manhattan shows one dot per
      # locus, colored by the number of ancestries that hit it. The region
      # lookup (for clicks + drilldown) still uses the underlying per-ancestry
      # region keys to match `filtered_data()`.
      is_multi <- !is.null(current_study()) && is_virtual_study(current_study()) &&
                  "ancestry" %in% names(dt)

      if (is_multi) {
        # Per (CHR, START, STOP) aggregation first (per-ancestry regions may
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
        # the most colocs (fallback: the widest). This key is what the
        # click handler forwards to selectize/filtered_data(), so all the
        # existing region-match plumbing stays unchanged.
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
        # Get max nlog10P per region from sumstats_1
        region_stats <- dt[, .(max_nlog10P = max(sumstats_1_max_nlog10P, na.rm = TRUE),
                               n_coloc = .N),
                           by = .(CHR_var, BP_START_var, BP_STOP_var)]
      }

      # Numeric chr for x-axis
      region_stats[, CHR_num := as.integer(gsub("X", "23", gsub("Y", "24", CHR_var)))]
      region_stats[, mid_pos := (BP_START_var + BP_STOP_var) / 2]

      # Cumulative position for Manhattan layout
      chr_offsets <- region_stats[, .(max_pos = max(mid_pos)), by = CHR_num]
      setorder(chr_offsets, CHR_num)
      chr_offsets[, offset := cumsum(c(0, head(max_pos, -1))) + (seq_len(.N) - 1) * 5e7]
      region_stats <- merge(region_stats, chr_offsets[, .(CHR_num, offset)], by = "CHR_num")
      region_stats[, x_pos := mid_pos + offset]

      # Gene labels
      if (is_multi) {
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
      if (is_multi) {
        # region_key points at the representative per-ancestry region so
        # clicks forward to a key the downstream filter recognizes.
        region_stats[, region_key := paste0(CHR_var, ":", rep_start, "-", rep_stop)]
      } else {
        region_stats[, region_key := paste0(CHR_var, ":", BP_START_var, "-", BP_STOP_var)]
      }

      # Dot color: alternating chr (default) or ancestry-coverage gradient
      # for the multi-ancestry virtual study.
      if (is_multi) {
        region_stats[, chr_color := ANCESTRY_COVERAGE_COLORS[as.character(n_ancestries)]]
      } else {
        region_stats[, chr_color := ifelse(CHR_num %% 2 == 0, "#3498db", "#2c3e50")]
      }

      # Highlight currently selected region (strips gene:: prefix internally)
      sel <- current_region()
      if (is.null(sel)) sel <- NA_character_
      region_stats[, is_selected := !is.na(sel) & region_key == sel]

      # Hover text
      if (is_multi) {
        has_n10 <- "nearest_genes_10" %in% names(region_stats)
        region_stats[, hover := paste0(
          "<b>", gene, "</b><br>",
          if (has_n10) paste0("Nearby: ", ifelse(is.na(nearest_genes_10), "", nearest_genes_10), "<br>") else "",
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

      # Chr tick labels
      chr_ticks <- region_stats[, .(mid_x = mean(x_pos)), by = CHR_num]
      setorder(chr_ticks, CHR_num)

      p <- plot_ly() %>%
        add_trace(
          data = region_stats[is_selected == FALSE],
          x = ~x_pos, y = ~max_nlog10P,
          type = "scatter", mode = "markers",
          marker = list(color = ~chr_color, size = ~pmin(5 + n_coloc / 10, 15), opacity = 0.7),
          text = ~hover, hoverinfo = "text",
          customdata = ~region_key,
          name = "Regions"
        )

      if (any(region_stats$is_selected)) {
        p <- p %>% add_trace(
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
      top_labeled <- region_stats[order(-n_coloc, -max_nlog10P)][1:min(20, nrow(region_stats))]
      selected_row <- region_stats[is_selected == TRUE]
      if (nrow(selected_row) > 0 && !selected_row$region_key[1] %in% top_labeled$region_key) {
        top_labeled <- rbind(top_labeled, selected_row)
      }
      # Clean gene name: remove "(within)" etc
      top_labeled[, gene_short := gsub("\\(.*\\)", "", gene)]
      top_labeled[, gene_short := trimws(gsub("INTERGENIC: ", "", gene_short))]

      gene_annotations <- lapply(seq_len(nrow(top_labeled)), function(i) {
        is_sel <- isTRUE(top_labeled$is_selected[i])
        list(
          x = top_labeled$x_pos[i],
          y = top_labeled$max_nlog10P[i],
          text = if (is_sel) paste0("<b>", top_labeled$gene_short[i], "</b>") else top_labeled$gene_short[i],
          showarrow = FALSE,
          font = list(size = if (is_sel) 15 else 13,
                      color = if (is_sel) "#e74c3c" else "#333"),
          yshift = 14
        )
      })

      p %>% layout(
        annotations = gene_annotations,
        xaxis = list(
          title = "", showgrid = FALSE,
          tickvals = chr_ticks$mid_x,
          ticktext = chr_ticks$CHR_num,
          tickfont = list(size = 13),
          range = c(min(region_stats$x_pos) * 0.98, max(region_stats$x_pos) * 1.02)
        ),
        yaxis = list(title = list(text = "-log10(P)", font = list(size = 14)),
                     tickfont = list(size = 12)),
        margin = list(t = 5, b = 25, l = 40, r = 5),
        showlegend = FALSE,
        hovermode = "closest"
      ) %>% config(displayModeBar = FALSE) %>%
        event_register("plotly_click")
    })

    # Handle Manhattan click -> select region
    # The selectize values use "gene::region" format; look up the first gene
    # in the clicked region from region_to_value_map.
    observeEvent(event_data("plotly_click", source = "A"), {
      click <- event_data("plotly_click", source = "A")
      if (!is.null(click) && !is.null(click$customdata)) {
        region_key <- click$customdata
        # Resolve to gene-prefixed value that selectize knows about
        map <- region_to_value_map()
        selected_val <- if (length(map) > 0 && region_key %in% names(map)) {
          map[[region_key]]
        } else {
          region_key  # fallback - strip_gene_prefix handles bare keys
        }
        updateSelectizeInput(session, "selected_region", selected = selected_val)
      }
    })

    # === Reactive Data Processing ===

    # Get current selected region (guard against NULL, NA, empty string)
    # Values from the gene-search selectize are prefixed with "gene::region"
    # so multiple genes in the same region can be distinguished. Strip the
    # gene prefix here to get the plain region key.
    strip_gene_prefix <- function(x) {
      if (is.null(x) || length(x) == 0 || is.na(x)) return(x)
      sub("^[^:]*::", "", x)
    }
    current_region <- reactive({
      if (.is_valid_sel(input$selected_region)) {
        strip_gene_prefix(input$selected_region)
      } else if (.is_valid_sel(input$selected_region_coord)) {
        strip_gene_prefix(input$selected_region_coord)
      } else {
        NULL
      }
    })
    
    # Load and process data
    coloc_data <- reactive({
      # First priority: auto-loaded study from home page
      if (!is.null(current_study())) {
        study <- current_study()

        # --- Virtual multi-ancestry study: rbind all ancestries ---
        if (is_virtual_study(study)) {
          vinfo <- DEFAULT_VIRTUAL_STUDIES[[study]]
          parts <- lapply(vinfo$ancestries, function(anc) {
            f <- vinfo$coloc_files[[anc]]
            if (is.null(f) || !file.exists(f)) return(NULL)
            dt <- readRDS(f)
            if (!is.data.table(dt)) dt <- as.data.table(dt)
            if (nrow(dt) == 0) return(NULL)
            dt[, ancestry := anc]
            # Strip ancestry suffix from per-study metadata columns so the
            # rest of the app sees a single "MVP_R4_*" schema instead of
            # MVP_R4_{AFR,AMR,...}_*.
            suffix <- paste0("^MVP_R4_", anc, "_")
            renamed <- grep(suffix, names(dt), value = TRUE)
            if (length(renamed) > 0) {
              setnames(dt, renamed, sub(suffix, "MVP_R4_", renamed))
            }
            # Unify source_study so it is stable across ancestries
            if ("source_study" %in% names(dt)) {
              dt[, source_study := "MVP_R4"]
            }
            dt
          })
          parts <- parts[!vapply(parts, is.null, logical(1))]
          if (length(parts) == 0) return(NULL)
          data <- data.table::rbindlist(parts, fill = TRUE)

          # Add common fallback columns (same as real-study branch below).
          if (!"Prioritized_Gene" %in% names(data)) {
            data[, Prioritized_Gene := nearest_gene_1]
          }
          if (!"clump_index_Name" %in% names(data)) {
            data[, clump_index_Name := NA_character_]
          }
          if (!"clump_index_P" %in% names(data)) {
            data[, clump_index_P := NA_real_]
          }
          if (!"region_center_pos" %in% names(data)) {
            data[, region_center_pos := as.numeric(BP_START_var + BP_STOP_var) / 2]
          }

          showNotification(
            paste0("Loaded ", study, " (",
                   paste(sort(unique(data$ancestry)), collapse = ", "), ")"),
            type = "message", duration = 3
          )
          return(data)
        }

        # --- Regular single-study path ---
        file_path <- resolve_annot_path(study)
        if (file.exists(file_path)) {
          data <- readRDS(file_path)
          # Ensure data is a data.table
          if (!is.data.table(data)) {
            data <- as.data.table(data)
          }

          # Add missing columns with fallbacks
          if (!"Prioritized_Gene" %in% names(data)) {
            data[, Prioritized_Gene := nearest_gene_1]
          }
          if (!"clump_index_Name" %in% names(data)) {
            data[, clump_index_Name := NA_character_]
          }
          if (!"clump_index_P" %in% names(data)) {
            data[, clump_index_P := NA_real_]
          }
          if (!"region_center_pos" %in% names(data)) {
            data[, region_center_pos := as.numeric(BP_START_var + BP_STOP_var) / 2]
          }

          # Filter out disabled studies

          showNotification(
            paste("Successfully loaded:", study),
            type = "message",
            duration = 3
          )
          return(data)
        }
      }

      return(NULL)
    })
    
    # Get unique regions
    regions_data <- reactive({
      req(coloc_data())
      dt <- coloc_data()

      # Select whichever metadata columns exist (atlas data may lack
      # Prioritized_Gene, clump_index_Name, etc.)
      meta_cols <- c("CHR_var", "BP_START_var", "BP_STOP_var")
      optional_cols <- c("region_center_pos", "nearest_gene_1",
                         "nearest_genes_10", "Prioritized_Gene",
                         "clump_index_Name")
      sel_cols <- c(meta_cols, intersect(optional_cols, names(dt)))

      regions <- unique(dt[, ..sel_cols])

      # Ensure columns the rest of the app expects are present (NA if missing)
      for (col in optional_cols) {
        if (!col %in% names(regions)) regions[, (col) := NA_character_]
      }
      regions[, region_id := paste0("chr", CHR_var, ":", BP_START_var, "-",
                                     BP_STOP_var)]

      # Count colocalizations per region
      region_counts <- dt[, .N, by = .(CHR_var, BP_START_var, BP_STOP_var)]
      regions <- merge(regions, region_counts,
                       by = c("CHR_var", "BP_START_var", "BP_STOP_var"))

      # Sort numerically by chromosome then position
      regions[, CHR_num := as.integer(gsub("X", "23", gsub("Y", "24", CHR_var)))]
      setorder(regions, CHR_num, BP_START_var)
      regions[, CHR_num := NULL]
      regions
    })
    
    # Filter data for selected region - full (no max_traits limit)
    # Used by convergence view which needs to see all cis molecular colocs
    filtered_region_data <- reactive({
      req(coloc_data(), current_region())

      region_parts <- strsplit(current_region(), ":")[[1]]
      chr <- region_parts[1]
      pos_parts <- strsplit(region_parts[2], "-")[[1]]
      start_pos <- as.numeric(pos_parts[1])
      end_pos <- as.numeric(pos_parts[2])

      # Virtual multi-ancestry: widen the match to the consensus window so
      # colocs from all ancestries that hit the same locus are included.
      if (!is.null(current_study()) && is_virtual_study(current_study())) {
        bundles <- load_multi_region_bundles(current_region())
        cons <- if (!is.null(bundles) && length(bundles) > 0) bundles[[1]]$.consensus else NULL
        if (!is.null(cons)) {
          dt <- coloc_data()[as.character(CHR_var) == as.character(cons$chr) &
                              BP_START_var <= cons$stop &
                              BP_STOP_var  >= cons$start]
        } else {
          dt <- coloc_data()[CHR_var == chr & BP_START_var == start_pos & BP_STOP_var == end_pos]
        }
      } else {
        dt <- coloc_data()[CHR_var == chr & BP_START_var == start_pos & BP_STOP_var == end_pos]
      }
      dt <- dt[PP.H4.abf >= input$min_pph4]
      dt <- dt[sumstats_2_max_nlog10P >= input$min_nlog10P]
      dt <- dt[source_study %in% selected_studies()]
      dt
    })

    # Filter data limited to max_traits (for network rendering performance)
    filtered_data <- reactive({
      dt <- filtered_region_data()
      if (nrow(dt) > input$max_traits) {
        dt <- dt[order(-PP.H4.abf)][1:input$max_traits]
      }
      dt
    })
    
    # Create network data
    network_data <- reactive({
      req(filtered_data())
      dt <- filtered_data()
      
      if (nrow(dt) == 0) {
        return(list(nodes = data.frame(), edges = data.frame()))
      }
      
      # Create central node (the region)
      # Select columns, including clump index if available
      cols_to_select <- c("Prioritized_Gene", "nearest_gene_1", "CHR_var", "BP_START_var", "BP_STOP_var")
      if ("clump_index_Name" %in% names(dt)) cols_to_select <- c(cols_to_select, "clump_index_Name")
      if ("clump_index_P" %in% names(dt)) cols_to_select <- c(cols_to_select, "clump_index_P")

      region_info <- unique(dt[, ..cols_to_select])

      # Use Prioritized_Gene if available, otherwise use nearest_gene_1
      gene_label <- if (!is.na(region_info$Prioritized_Gene[1]) && region_info$Prioritized_Gene[1] != "") {
        region_info$Prioritized_Gene[1]
      } else {
        region_info$nearest_gene_1[1]
      }

      # Build tooltip with both genes and clump index info
      tooltip_text <- paste0("Region: chr", region_info$CHR_var[1], ":",
                            region_info$BP_START_var[1], "-", region_info$BP_STOP_var[1], "<br>",
                            "Prioritized Gene: ", ifelse(!is.na(region_info$Prioritized_Gene[1]),
                                                         region_info$Prioritized_Gene[1], "N/A"), "<br>",
                            "Nearest Gene: ", region_info$nearest_gene_1[1])

      # Add clump index information if available
      if ("clump_index_Name" %in% names(region_info)) {
        tooltip_text <- paste0(tooltip_text, "<br>Clump Index Name: ",
                              ifelse(!is.na(region_info$clump_index_Name[1]),
                                    region_info$clump_index_Name[1], "N/A"))
      }
      if ("clump_index_P" %in% names(region_info)) {
        # Format p-value in scientific notation with 2 decimal places
        formatted_p <- ifelse(!is.na(region_info$clump_index_P[1]),
                             sprintf("%.2e", region_info$clump_index_P[1]),
                             "N/A")
        tooltip_text <- paste0(tooltip_text, "<br>Clump Index P: ", formatted_p)
      }

      central_node <- data.frame(
        id = "region",
        label = gene_label,
        title = tooltip_text,
        group = "region",
        shape = "star",
        size = 30,
        color = "#FFD700",
        stringsAsFactors = FALSE
      )
      
      # Create trait nodes - make sure column structure matches
      trait_nodes_list <- dt[, {
        # Extract trait name based on source study
        trait_name <- switch(source_study[1],
                             "CKDGen_r4" = if ("CKDGen_r4_Name" %in% names(.SD)) CKDGen_r4_Name else NA,
                             "CKDGen_r5" = CKDGen_r5_ckdgen_r5_name,
                             "FinnGen_r9" = FinnGen_r9_phenotype,
                             "GCKD_mGWAS_plasma" = GCKD_mGWAS_plasma_BIOCHEMICAL,
                             "GCKD_mGWAS_urine" = GCKD_mGWAS_urine_BIOCHEMICAL,
                             "Icelanders_pGWAS" = Icelanders_pGWAS_Protein..short.name.,
                             "MVP_R4" = MVP_R4_Analyzed.variable,
                             "MVP_R4_EUR" = if ("MVP_R4_EUR_Analyzed.variable" %in% names(.SD)) MVP_R4_EUR_Analyzed.variable else NA,
                             "UKB_PPP_EUR" = UKB_PPP_EUR_olink_target_fullname,
                             "UKB_TOPMed" = UKB_TOPMed_phenostring,
                             "UKB_kidney_vol" = gsub("model1_qnorm_(.*)_chr.*", "\\1", basename(sumstats_2_file)),
                             "Kidney_eQTL" = Kidney_eQTL_gene_name,
                             "pho_ca" = pho_ca_pho_ca_name,
                             "eQTLGen" = eQTLGen_gene_name,
                             "GTEXv8_eQTL" = GTEXv8_eQTL_gene_name,
                             NA
        )

        # Create node ID
        node_id <- paste0(source_study, "_", .I)

        list(
          id = node_id,
          label = ifelse(is.na(trait_name), source_study,
                         substr(trait_name, 1, 30)),  # Truncate long names
          title = paste0(
            "Study: ", source_study, "<br>",
            "Trait: ", ifelse(is.na(trait_name), "N/A", trait_name), "<br>",
            "PP.H4: ", round(PP.H4.abf, 3), "<br>",
            "Significance: ", round(sumstats_2_max_nlog10P, 2), " (-log10 P)"
          ),
          group = source_study,
          shape = "dot",
          size = 10 + 20 * PP.H4.abf  # Size by PP.H4
        )
      }, by = .I]

      # Convert to data.frame
      trait_nodes <- data.frame(
        id = trait_nodes_list$id,
        label = trait_nodes_list$label,
        title = trait_nodes_list$title,
        group = trait_nodes_list$group,
        shape = trait_nodes_list$shape,
        size = trait_nodes_list$size,
        stringsAsFactors = FALSE
      )

      # Add colors: use color picker input if available, otherwise default
      # Note: #000000 is textInput's default, treat it as "not set"
      trait_nodes$color <- unname(sapply(trait_nodes$group, function(study) {
        color_input <- input[[paste0("color_", study)]]
        if (!is.null(color_input) && color_input != "#000000") {
          color_input
        } else {
          study_colors[study]
        }
      }))

      nodes <- rbind(central_node, trait_nodes)

      # Apply custom labels if any
      labels <- custom_labels()
      if (length(labels) > 0) {
        for (node_id in names(labels)) {
          if (node_id %in% nodes$id) {
            nodes$label[nodes$id == node_id] <- labels[[node_id]]
          }
        }
      }

      # Create edges with directionality colors and tooltips
      # Get current study name for tooltip
      study_name <- if (!is.null(current_study())) current_study() else "Region"

      edge_list <- dt[, {
        # Determine edge color and label based on directionality (colorblind-friendly)
        dir_label <- "N/A"
        edge_color <- "#848484"  # Default gray

        if ("directionality" %in% names(dt)) {
          dir_val <- as.character(directionality)
          if (!is.na(directionality) && dir_val != "") {
            # Check if positive
            if (grepl("positive|concordant", dir_val, ignore.case = TRUE) ||
                (suppressWarnings(!is.na(as.numeric(dir_val))) && as.numeric(dir_val) > 0)) {
              edge_color <- "#e67e22"  # Orange for positive/concordant
              dir_label <- "Positive"
            }
            # Check if negative
            else if (grepl("negative|discordant", dir_val, ignore.case = TRUE) ||
                     (suppressWarnings(!is.na(as.numeric(dir_val))) && as.numeric(dir_val) < 0)) {
              edge_color <- "#3498db"  # Blue for negative/discordant
              dir_label <- "Negative"
            }
          }
        }

        # Build edge tooltip
        edge_title <- paste0(
          "PP.H4: ", round(PP.H4.abf, 3), "<br>",
          "Directionality: ", dir_label, "<br>"
        )

        # Add index SNP info if available
        if ("sumstats_1_ind_Name" %in% names(dt)) {
          edge_title <- paste0(edge_title, "Index SNP: ",
                              ifelse(is.na(sumstats_1_ind_Name), "N/A", sumstats_1_ind_Name), "<br>")
        }
        if ("sumstats_1_ind_nlog10P" %in% names(dt)) {
          edge_title <- paste0(edge_title, "P-value in ", study_name, ": ",
                              round(sumstats_1_ind_nlog10P, 2), " (-log10)", "<br>")
        }
        if ("sumstats_2_ind_nlog10P" %in% names(dt)) {
          edge_title <- paste0(edge_title, "P-value in sumstats_2: ",
                              round(sumstats_2_ind_nlog10P, 2), " (-log10)")
        }

        list(
          edge_color = edge_color,
          edge_title = edge_title
        )
      }, by = .I]

      edges <- data.frame(
        from = rep("region", nrow(trait_nodes)),
        to = trait_nodes$id,
        width = 4,  # Thicker edges
        color = edge_list$edge_color,
        title = edge_list$edge_title,
        stringsAsFactors = FALSE
      )

      list(nodes = nodes, edges = edges)
    })

    # === Outputs ===
    
    # Display current study info
    output$current_study_info <- renderUI({
      if (!is.null(current_study())) {
        div(
          style = "background-color: #e8f4f8; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
          p(HTML(paste0("<strong>Currently Loaded:</strong> ",
                        trait_label(current_study()))),
            style = "margin: 0;")
        )
      } else {
        div(
          style = "background-color: #f5f5f5; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
          p("No study auto-loaded. Upload your own data below.", style = "margin: 0; color: #666;")
        )
      }
    })
    
    # Create combined study selector with colors, checkboxes, and info buttons
    # Dynamically discover studies from current coloc data
    available_study_names <- reactive({
      req(coloc_data())
      sort(unique(coloc_data()$source_study))
    })

    output$study_selector <- renderUI({
      req(available_study_names())
      counts <- study_counts_for_region()
      studies <- available_study_names()

      study_items <- lapply(studies, function(study) {
        display_name <- get_study_display_name(study)
        color <- if (study %in% names(study_colors)) study_colors[[study]] else "#7f7f7f"

        # Count for current region
        study_count <- if (!is.null(counts) && study %in% names(counts)) counts[study] else 0
        count_style <- if (study_count == 0) {
          "font-size: 10px; color: #999; margin-left: 3px;"
        } else {
          "font-size: 10px; color: #27ae60; font-weight: bold; margin-left: 3px;"
        }

        div(style = "margin: 0; padding: 0; display: flex; align-items: center;",
            tags$div(style = "margin: 0; padding: 0;",
              checkboxInput(paste0("study_", study),
                           label = NULL, value = TRUE, width = "20px")),
            tags$span(style = paste0("display: inline-block; width: 12px; height: 12px; ",
                                    "background-color: ", color,
                                    "; margin-left: -5px; margin-right: 6px; border: 1px solid #ddd;")),
            tags$span(style = "font-size: 11px;", display_name),
            tags$span(style = count_style, paste0("(", study_count, ")"))
        )
      })

      do.call(tagList, study_items)
    })

    # Create reactive for selected studies (combines individual checkboxes)
    # Uses dynamic study list from current coloc data
    selected_studies <- reactive({
      studies <- available_study_names()
      if (length(studies) == 0) return(character(0))
      selected <- sapply(studies, function(study) {
        val <- input[[paste0("study_", study)]]
        if (is.null(val)) TRUE else val  # Default to TRUE
      })
      studies[selected]
    })

    # Compute coloc counts per study for current region (before study filtering)
    study_counts_for_region <- reactive({
      req(coloc_data(), current_region())

      # Parse region
      region_parts <- strsplit(current_region(), ":")[[1]]
      chr <- region_parts[1]
      pos_parts <- strsplit(region_parts[2], "-")[[1]]
      start_pos <- as.numeric(pos_parts[1])
      end_pos <- as.numeric(pos_parts[2])

      # Filter data for region only (apply PP.H4 and nlog10P filters, but NOT study filter)
      dt <- coloc_data()[CHR_var == chr & BP_START_var == start_pos & BP_STOP_var == end_pos]
      dt <- dt[PP.H4.abf >= input$min_pph4]
      dt <- dt[sumstats_2_max_nlog10P >= input$min_nlog10P]

      # Count by source_study
      counts <- dt[, .N, by = source_study]
      setNames(counts$N, counts$source_study)
    })

    # Create observers for info buttons
    lapply(names(study_colors), function(study) {
      observeEvent(input[[paste0("info_", study)]], {
        showModal(modalDialog(
          title = get_study_display_name(study),
          format_study_info(study),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      })
    })
    
    # === Convergence view ===

    # Reactive: current region coords (chr, start, end)
    conv_region_coords <- reactive({
      req(current_region())
      parts <- strsplit(current_region(), ":")[[1]]
      if (length(parts) < 2) return(NULL)
      chr <- sub("^chr", "", parts[1])
      pos <- strsplit(parts[2], "-")[[1]]
      list(chr = chr,
           start = as.numeric(pos[1]),
           end   = as.numeric(pos[2]))
    })

    # Panel 1: base trait regional plot (reuses existing plotting logic)
    output$conv_base_plot <- renderPlotly({
      req(regional_base_sumstats())
      # Use full region bounds as x-range so the base plot and the selected
      # trait plot (and gene track) share the same x-axis for comparison
      coords <- conv_region_coords()
      xr <- if (!is.null(coords)) c(coords$start, coords$end) else NULL
      plot_regional_association_interactive(
        regional_base_sumstats(),
        title = paste0(current_study(), " - ", current_region()),
        highlight_snp = NULL,
        color = "#2c3e50",
        x_range = xr
      )
    })

    # Panel 2: gene track - reuse existing plot_gene_track from gene_track.R
    output$conv_gene_track <- renderPlotly({
      coords <- conv_region_coords()
      req(coords)
      plot_gene_track(coords$chr, coords$start, coords$end,
                      source = "conv_gene_track")
    })

    # Handle click on convergence gene track - populates selected_gene
    # (the same reactiveVal used by gene_info_panel)
    observeEvent(event_data("plotly_click", source = "conv_gene_track"), {
      click_data <- event_data("plotly_click", source = "conv_gene_track")
      if (!is.null(click_data) && !is.null(click_data$customdata)) {
        clicked_gene <- click_data$customdata
        if (!is.null(gene_annotation)) {
          gene_info <- gene_annotation[gene_name == clicked_gene]
          if (nrow(gene_info) > 0) {
            selected_gene(gene_info[1])
          }
        }
      }
    })

    # Clear gene info panel when navigating to a new region
    observeEvent(current_region(), {
      selected_gene(NULL)
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    # === Panel 3: category overview cards ===

    # Which category is currently selected for drill-down (NULL = none)
    conv_selected_cat <- reactiveVal(NULL)

    # Reset selection when region changes
    observeEvent(current_region(), {
      conv_selected_cat(NULL)
    })

    # Category counts for current region
    conv_cat_counts <- reactive({
      dt <- filtered_region_data()
      count_by_category(dt)
    })

    output$conv_category_cards <- renderUI({
      counts <- conv_cat_counts()
      selected <- conv_selected_cat()

      cards <- lapply(names(TRAIT_CATEGORIES), function(cat_id) {
        meta <- TRAIT_CATEGORIES[[cat_id]]
        n <- tryCatch(counts[[cat_id]], error = function(e) 0L)
        if (is.null(n) || length(n) == 0 || is.na(n)) n <- 0L
        is_empty <- n == 0
        is_active <- !is.null(selected) && selected == cat_id

        # Build inline style based on state
        border_style <- if (is_active) {
          paste0("border: 2px solid ", meta$color, ";")
        } else {
          "border: 1px solid #ddd;"
        }
        opacity_style <- if (is_empty) "opacity: 0.45;" else ""
        cursor_style <- if (is_empty) "cursor: default;" else "cursor: pointer;"
        bg_style <- if (is_active) {
          paste0("background: ", meta$color, "15;")  # 15 = 8% alpha
        } else {
          "background: white;"
        }

        actionButton(
          inputId = paste0("conv_cat_btn_", cat_id),
          label = HTML(paste0(
            "<div style='display: flex; align-items: center; gap: 10px; text-align: left;'>",
            "<div style='width: 4px; height: 32px; background: ", meta$color,
              "; border-radius: 2px;'></div>",
            "<div style='font-size: 1.4em;'>", meta$icon, "</div>",
            "<div style='flex: 1; min-width: 0;'>",
              "<div style='font-size: 12px; font-weight: 600; color: #2c3e50; line-height: 1.2;'>",
                meta$label, "</div>",
              "<div style='font-size: 16px; font-weight: 700; color: ",
                if (is_empty) "#aaa" else "#2c3e50", ";'>", n, "</div>",
            "</div>",
            "</div>"
          )),
          style = paste0(
            "width: 100%; padding: 10px; margin-bottom: 6px; ",
            "border-radius: 8px; ",
            border_style, opacity_style, cursor_style, bg_style,
            "text-align: left; box-shadow: 0 1px 3px rgba(0,0,0,0.04);"
          ),
          disabled = is_empty
        )
      })

      do.call(tagList, cards)
    })

    # Handle category button clicks (one observer per category)
    lapply(names(TRAIT_CATEGORIES), function(cat_id) {
      btn_id <- paste0("conv_cat_btn_", cat_id)
      observeEvent(input[[btn_id]], {
        counts <- conv_cat_counts()
        n <- tryCatch(counts[[cat_id]], error = function(e) 0L)
        if (is.null(n) || length(n) == 0 || is.na(n) || n == 0) return()
        # Toggle: click active category to close
        if (identical(conv_selected_cat(), cat_id)) {
          conv_selected_cat(NULL)
        } else {
          conv_selected_cat(cat_id)
        }
      }, ignoreInit = TRUE)
    })

    # Drill-down header (shows selected category and clear button)
    output$conv_drilldown_header <- renderUI({
      sel <- conv_selected_cat()
      if (is.null(sel)) {
        return(div(style = "font-size: 11px; color: #999; margin: 12px 0 4px 0; font-style: italic;",
                   "Select a category above to see traits."))
      }
      meta <- TRAIT_CATEGORIES[[sel]]
      counts <- conv_cat_counts()
      n_total <- tryCatch(counts[[sel]], error = function(e) 0L)
      if (is.null(n_total) || length(n_total) == 0 || is.na(n_total)) n_total <- 0L
      n_shown <- min(n_total, 50L)
      h6(style = "margin: 12px 0 4px 0;",
         meta$icon, " ", meta$label,
         tags$span(style = "font-size: 11px; color: #888; margin-left: 8px;",
                   sprintf("(showing %d of %d)", n_shown, n_total)))
    })

    # Data for drill-down: filter by selected category, limit to 50 by signal strength.
    # For virtual multi-ancestry studies we additionally collapse rows per
    # consensus trait (one row per trait across ancestries) and attach:
    #   - `ancestries`   : comma-separated list of ancestries that hit
    #   - `n_ancestries` : count
    #   - `consensus_key`: canonical trait key for overlay plot lookup
    # The max sumstats_2_max_nlog10P is kept for sorting. Columns expected by
    # downstream code (source_study, sumstats_2_file, ancestry) are taken
    # from the ancestry with the strongest signal so display helpers keep
    # working without branching.
    conv_drilldown_data <- reactive({
      sel <- conv_selected_cat()
      if (is.null(sel)) return(NULL)
      dt <- filtered_region_data()
      studies_in_cat <- TRAIT_CATEGORIES[[sel]]$studies
      dt <- dt[source_study %in% studies_in_cat]
      if (nrow(dt) == 0) return(NULL)

      if (!is.null(current_study()) && is_virtual_study(current_study()) &&
          "ancestry" %in% names(dt)) {
        dt[, .ck := make_consensus_trait_key(
          make_bundle_key(source_study, sumstats_2_file, ancestry))]
        # Rank rows by signal so the "first row per consensus" is the
        # strongest ancestry signal.
        if ("sumstats_2_max_nlog10P" %in% names(dt)) {
          setorder(dt, -sumstats_2_max_nlog10P)
        }
        dt <- dt[, {
          # Capture group-level ancestry BEFORE any nested data.table op;
          # once we enter `first[, ... ]` the inner .SD refers to `first`
          # (one row) rather than the outer group.
          group_ancs <- sort(unique(ancestry))
          group_ck   <- .ck[1]
          first <- .SD[1]  # strongest row (group is pre-sorted by signal)
          first[, ancestries    := paste(group_ancs, collapse = ",")]
          first[, n_ancestries  := length(group_ancs)]
          first[, consensus_key := group_ck]
          first
        }, by = .ck]
        dt[, .ck := NULL]
        if ("sumstats_2_max_nlog10P" %in% names(dt)) {
          setorder(dt, -sumstats_2_max_nlog10P)
        }
        return(dt)
      }

      # Sort by sumstats_2 signal strength (full list; tiles cap to 50)
      if ("sumstats_2_max_nlog10P" %in% names(dt)) {
        setorder(dt, -sumstats_2_max_nlog10P)
      }
      dt
    })

    # Top 50 subset for tile rendering performance
    conv_drilldown_data_top <- reactive({
      dt <- conv_drilldown_data()
      if (is.null(dt) || nrow(dt) <= 50) return(dt)
      dt[1:50]
    })

    # Drill-down tiles: top 50 traits by -log10P (full list via dropdown)
    output$conv_trait_tiles <- renderUI({
      dt <- conv_drilldown_data_top()
      if (is.null(dt) || nrow(dt) == 0) {
        return(div(style = "font-size: 11px; color: #999; font-style: italic; margin: 12px 0;",
                   "Select a category from the left panel to see traits."))
      }

      sel_cat <- conv_selected_cat()
      cat_meta <- TRAIT_CATEGORIES[[sel_cat]]

      # Build metadata for each tile and the id->bundle_key map. For
      # virtual studies the drilldown already carries a `consensus_key`
      # and `ancestries` column; prefer those so clicking a tile selects
      # the consensus trait (overlay) rather than one ancestry.
      is_virt_dd <- is_virtual_study(current_study()) && "consensus_key" %in% names(dt)
      tiles_list <- dt[, {
        study <- source_study[1]
        anc <- if ("ancestry" %in% names(.SD)) ancestry[1] else NA_character_
        trait_name <- tryCatch(get_trait_display_name(.SD),
                               error = function(e) basename(sumstats_2_file[1]))
        if (is.null(trait_name) || is.na(trait_name) || trait_name == "") {
          trait_name <- basename(sumstats_2_file[1])
        }
        study_color <- if (exists("study_colors") && study %in% names(study_colors)) {
          study_colors[[study]]
        } else cat_meta$color
        nlp <- if ("sumstats_2_max_nlog10P" %in% names(.SD)) sumstats_2_max_nlog10P[1] else NA
        bk <- if (is_virt_dd) consensus_key[1] else make_bundle_key(study, sumstats_2_file[1], anc)
        anc_badge <- if (is_virt_dd && "ancestries" %in% names(.SD)) ancestries[1] else NA_character_
        list(
          tile_id = paste0("conv_t_", .I),
          bundle_key = bk,
          trait_name = trait_name,
          study = study,
          study_color = study_color,
          nlp = nlp,
          anc_badge = anc_badge
        )
      }, by = .I]

      # Update node/tile -> bundle_key map
      conv_node_keys(setNames(tiles_list$bundle_key, tiles_list$tile_id))

      # Build tile elements
      tiles <- lapply(seq_len(nrow(tiles_list)), function(i) {
        tl <- tiles_list[i]
        is_selected <- !is.null(conv_selected_trait()) &&
                        conv_selected_trait() == tl$bundle_key
        border_color <- if (is_selected) "#3498db" else "#e0e0e0"
        border_width <- if (is_selected) "2px" else "1px"
        bg <- if (is_selected) "#f0f8ff" else "white"
        p_text <- if (is.na(tl$nlp)) "-" else sprintf("%.1f", tl$nlp)

        div(
          class = "conv-tile",
          onclick = paste0("Shiny.setInputValue('conv_trait_click','",
                            tl$tile_id, "', {priority:'event'});"),
          style = paste0(
            "cursor: pointer; padding: 6px 8px; margin: 3px; ",
            "border: ", border_width, " solid ", border_color, "; ",
            "border-radius: 6px; background: ", bg, "; ",
            "transition: all 0.15s; ",
            "display: flex; flex-direction: column; gap: 2px; ",
            "min-height: 52px;"
          ),
          title = paste0(tl$trait_name, " (", tl$study,
                         ", -log10P=", p_text, ")"),
          # Trait name (truncated with ellipsis via CSS)
          div(style = "font-size: 11px; font-weight: 600; color: #2c3e50; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;",
              tl$trait_name),
          # Study badge + p-value on one line.
          # For virtual studies, replace the study badge with ancestry
          # badges so the tile immediately shows how broadly the trait
          # colocalized.
          div(style = "display: flex; justify-content: space-between; align-items: center; gap: 4px;",
              if (!is.null(tl$anc_badge) && !is.na(tl$anc_badge) && nzchar(tl$anc_badge)) {
                tags$span(
                  style = paste0(
                    "display: inline-block; padding: 1px 5px; ",
                    "border-radius: 8px; font-size: 9px; font-weight: 600; ",
                    "color: white; background: ", tl$study_color, "; ",
                    "white-space: nowrap; overflow: hidden; text-overflow: ellipsis; ",
                    "max-width: 70%;"
                  ),
                  tl$anc_badge
                )
              } else {
                tags$span(
                  style = paste0(
                    "display: inline-block; padding: 1px 5px; ",
                    "border-radius: 8px; font-size: 9px; font-weight: 600; ",
                    "color: white; background: ", tl$study_color, "; ",
                    "white-space: nowrap; overflow: hidden; text-overflow: ellipsis; ",
                    "max-width: 70%;"
                  ),
                  tl$study
                )
              },
              tags$span(
                style = "font-size: 10px; color: #555; font-weight: 600;",
                paste0("-log10P=", p_text)
              )
          )
        )
      })

      # Wrap tiles in a responsive grid (5 columns on wide screens)
      div(
        style = "display: grid; grid-template-columns: repeat(auto-fill, minmax(160px, 1fr)); gap: 4px; margin-top: 8px;",
        tiles
      )
    })

    # Track which trait is selected for regional plot display
    conv_selected_trait <- reactiveVal(NULL)
    # Map node id -> bundle key (populated by drilldown renderer)
    conv_node_keys <- reactiveVal(character(0))

    observeEvent(input$conv_trait_click, {
      node_id <- input$conv_trait_click
      keys <- conv_node_keys()
      bundle_key <- keys[[node_id]]
      if (!is.null(bundle_key)) {
        conv_selected_trait(bundle_key)
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # When region changes: pick the first non-empty category and its top trait
    observeEvent(current_region(), {
      counts <- conv_cat_counts()
      # Find first category with > 0 entries
      non_empty <- names(counts)[counts > 0]
      if (length(non_empty) == 0) {
        conv_selected_cat(NULL)
        conv_selected_trait(NULL)
        return()
      }
      conv_selected_cat(non_empty[1])
      # The drill-down data reactive will emit; we set the trait in a
      # dedicated observer below (driven by conv_drilldown_data)
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    # When drill-down data updates (region or category changed), auto-select
    # the top trait by signal strength
    observeEvent(conv_drilldown_data(), {
      dt <- conv_drilldown_data()
      if (is.null(dt) || nrow(dt) == 0) {
        conv_selected_trait(NULL)
        return()
      }
      # Top row (already sorted by -sumstats_2_max_nlog10P in conv_drilldown_data).
      # Prefer the consensus key if the drilldown is virtual-collapsed.
      bundle_key <- if ("consensus_key" %in% names(dt) && !is.na(dt$consensus_key[1])) {
        dt$consensus_key[1]
      } else {
        anc <- if ("ancestry" %in% names(dt)) dt$ancestry[1] else NA_character_
        make_bundle_key(dt$source_study[1], dt$sumstats_2_file[1], anc)
      }
      conv_selected_trait(bundle_key)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    # Load selected trait sumstats. For virtual multi-ancestry studies
    # conv_selected_trait() is a consensus key - resolve it against each
    # ancestry's loaded bundle and return a named list (one per ancestry)
    # so plot_regional_association_interactive() overlays the traces.
    # Attach an attribute `coloc_ancestries` listing which ancestries
    # actually colocalized this trait (rest shown as raw-sumstats context).
    conv_trait_sumstats <- reactive({
      req(regional_data_path(), regional_layout(), current_region(), conv_selected_trait())
      sel <- conv_selected_trait()

      if (!is.null(current_study()) && is_virtual_study(current_study())) {
        bundles <- load_multi_region_bundles(current_region())
        if (is.null(bundles) || length(bundles) == 0) return(NULL)
        out <- lapply(names(bundles), function(anc) {
          bk <- names(bundles[[anc]]$traits)
          ak <- resolve_consensus_in_bundle(sel, bk, anc)
          if (is.null(ak)) return(NULL)
          bundles[[anc]]$traits[[ak]]
        })
        names(out) <- names(bundles)
        out <- out[!vapply(out, is.null, logical(1))]
        if (length(out) == 0) return(NULL)

        # Tag which ancestries actually colocalized (from the drilldown
        # data for the current consensus_key). Used by the overlay to
        # fade non-colocalized traces.
        dd <- conv_drilldown_data()
        coloc_ancs <- character(0)
        if (!is.null(dd) && "consensus_key" %in% names(dd)) {
          row <- dd[consensus_key == sel]
          if (nrow(row) > 0) {
            coloc_ancs <- strsplit(row$ancestries[1], ",", fixed = TRUE)[[1]]
          }
        }
        attr(out, "coloc_ancestries") <- coloc_ancs
        return(out)
      }

      load_trait_sumstats(regional_data_path(), current_region(),
                          sel, regional_layout())
    })

    # Display trait name for header
    conv_selected_trait_label <- reactive({
      sel <- conv_selected_trait()
      if (is.null(sel)) return(NULL)
      # Virtual studies: look up the consensus row in the drilldown data
      # and render "<name> [ANC1,ANC2,...]".
      if (is_virtual_study(current_study())) {
        dd <- conv_drilldown_data()
        if (!is.null(dd) && "consensus_key" %in% names(dd)) {
          row <- dd[consensus_key == sel]
          if (nrow(row) > 0) {
            nm <- tryCatch(get_trait_display_name(row[1]),
                           error = function(e) NULL)
            if (!is.null(nm) && !is.na(nm) && nm != "") {
              return(paste0(nm, " [", row$ancestries[1], "]"))
            }
          }
        }
      }
      # Parse "study__file"
      parts <- strsplit(sel, "__", fixed = TRUE)[[1]]
      if (length(parts) < 2) return(sel)
      study <- parts[1]
      file_part <- paste(parts[-1], collapse = "__")
      # Try to find display name from filtered_region_data
      dt <- filtered_region_data()
      if (!is.null(dt) && nrow(dt) > 0) {
        row <- dt[source_study == study & basename(sumstats_2_file) == file_part]
        if (nrow(row) > 0) {
          nm <- tryCatch(get_trait_display_name(row[1]),
                         error = function(e) NULL)
          if (!is.null(nm) && !is.na(nm) && nm != "") {
            return(paste0(nm, " (", study, ")"))
          }
        }
      }
      sel
    })

    # Signal for conditionalPanel: is a trait selected?
    output$conv_trait_has_data <- reactive({
      !is.null(conv_selected_trait())
    })
    outputOptions(output, "conv_trait_has_data", suspendWhenHidden = FALSE)

    # Regional plot for selected trait
    output$conv_trait_plot <- renderPlotly({
      req(conv_selected_trait())
      dt <- conv_trait_sumstats()
      # dt can be a single data.table or a named list of them (multi-ancestry
      # overlay). Treat both uniformly for the empty check.
      is_empty <- is.null(dt) || (is.data.frame(dt) && nrow(dt) == 0) ||
                  (is.list(dt) && !is.data.frame(dt) && length(dt) == 0)
      if (is_empty) {
        return(plotly_empty() %>% layout(
          title = list(text = "No regional plot data available for this trait",
                        font = list(size = 11))))
      }
      label <- conv_selected_trait_label()
      coords <- conv_region_coords()
      # For virtual studies x_range should span the consensus window.
      xr <- if (!is.null(current_study()) && is_virtual_study(current_study())) {
        bundles <- load_multi_region_bundles(current_region())
        cons <- if (!is.null(bundles) && length(bundles) > 0) bundles[[1]]$.consensus else NULL
        if (!is.null(cons)) c(cons$start, cons$stop) else NULL
      } else if (!is.null(coords)) {
        c(coords$start, coords$end)
      } else NULL
      plot_regional_association_interactive(
        dt,
        title = label,
        color = "#E67E22",
        x_range = xr
      )
    })

    output$conv_trait_plot_header <- renderUI({
      req(conv_selected_trait())
      dt <- conv_drilldown_data()
      if (is.null(dt) || nrow(dt) == 0) {
        return(h6(style = "margin: 12px 0 2px 0; color: #555;",
                  icon("chart-line"), " Selected trait"))
      }

      # Build "label -> bundle_key" choices for all traits in the current
      # category (same order as the tiles: by -log10P desc). Virtual
      # multi-ancestry studies use consensus keys with ancestry badges.
      is_virt_dd <- is_virtual_study(current_study()) && "consensus_key" %in% names(dt)
      choices_list <- lapply(seq_len(nrow(dt)), function(i) {
        row <- dt[i]
        study <- row$source_study
        anc <- if ("ancestry" %in% names(row)) row$ancestry else NA_character_
        trait_name <- tryCatch(get_trait_display_name(row),
                               error = function(e) basename(row$sumstats_2_file))
        if (is.null(trait_name) || is.na(trait_name) || trait_name == "") {
          trait_name <- basename(row$sumstats_2_file)
        }
        nlp <- if ("sumstats_2_max_nlog10P" %in% names(row)) row$sumstats_2_max_nlog10P else NA
        if (is_virt_dd) {
          ancs <- row$ancestries
          label <- paste0(trait_name, " [", ancs, "], -log10P=",
                          if (is.na(nlp)) "-" else sprintf("%.1f", nlp))
          key <- row$consensus_key
        } else {
          disp_study <- if (!is.na(anc) && nzchar(anc)) paste0(study, " ", anc) else study
          label <- paste0(trait_name, " (", disp_study, ", -log10P=",
                          if (is.na(nlp)) "-" else sprintf("%.1f", nlp), ")")
          key <- make_bundle_key(study, row$sumstats_2_file, anc)
        }
        list(label = label, key = key)
      })
      labels <- vapply(choices_list, `[[`, character(1), "label")
      keys   <- vapply(choices_list, `[[`, character(1), "key")
      choices <- setNames(keys, labels)

      div(style = "display: flex; align-items: center; gap: 8px; margin: 12px 0 2px 0;",
        h6(style = "margin: 0; color: #555; white-space: nowrap;",
           icon("chart-line"), " Selected trait:"),
        div(style = "flex: 1;",
          selectizeInput("conv_trait_dropdown", label = NULL,
                         choices = choices,
                         selected = conv_selected_trait(),
                         width = "100%",
                         options = list(closeAfterSelect = TRUE))
        )
      )
    })

    # Dropdown change -> update selected trait (and tile highlight follows)
    observeEvent(input$conv_trait_dropdown, {
      val <- input$conv_trait_dropdown
      if (is.null(val) || length(val) == 0 || is.na(val) || val == "") return()
      if (!identical(conv_selected_trait(), val)) {
        conv_selected_trait(val)
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # Render network
    output$network <- renderVisNetwork({
      net_data <- network_data()

      if (nrow(net_data$nodes) == 0) {
        visNetwork(data.frame(id = 1, label = "No data"), data.frame()) %>%
          visNodes(shape = "text", font = list(size = 20))
      } else {
        # Filter nodes based on selection (if any selected)
        visible <- input$visible_nodes
        if (!is.null(visible) && length(visible) > 0) {
          # Filter nodes to only show selected ones
          nodes_filtered <- net_data$nodes[net_data$nodes$id %in% visible, ]
          # Filter edges to only include those connecting visible nodes
          edges_filtered <- net_data$edges[
            net_data$edges$from %in% visible & net_data$edges$to %in% visible,
          ]
        } else {
          # No selection means show all (during initial load)
          nodes_filtered <- net_data$nodes
          edges_filtered <- net_data$edges
        }

        # Show notification for large networks (dismissed when stabilized)
        n_nodes <- nrow(nodes_filtered)
        if (n_nodes > 100) {
          showNotification(
            paste0("Rendering network (", n_nodes, " nodes)... please wait"),
            id = "region_network_loading",
            duration = NULL,
            type = "message"
          )
        }

        visNetwork(nodes_filtered, edges_filtered) %>%
          visNodes(
            font = list(color = "#000000", size = 14, face = "arial")
          ) %>%
          visOptions(
            highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)
          ) %>%
          visPhysics(enabled = input$physics,
                     stabilization = list(enabled = TRUE, iterations = 200)) %>%
          visLayout(randomSeed = 123) %>%
          visInteraction(
            dragNodes = TRUE,
            dragView = TRUE,
            zoomView = TRUE
          ) %>%
          visEvents(
            select = "function(nodes) {
              Shiny.setInputValue('selected_node', nodes.nodes[0]);
            }",
            stabilizationIterationsDone = "function() {
              Shiny.setInputValue('region_network_stabilized', Math.random());
            }"
          )
      }
    })

    # Dismiss loading notification when region network stabilizes
    observeEvent(input$region_network_stabilized, {
      removeNotification(id = "region_network_loading")
    })
    
    # Show node info when selected
    output$node_info <- renderPrint({
      req(input$selected_node)
      net_data <- network_data()
      
      node <- net_data$nodes[net_data$nodes$id == input$selected_node, ]
      if (nrow(node) > 0) {
        cat("Node ID:", node$id, "\n")
        cat("Label:", node$label, "\n")
        cat("Group:", node$group, "\n")
        # Parse and display HTML title content
        title_text <- gsub("<br>", "\n", node$title)
        cat("\nDetails:\n", title_text)
      }
    })
    # Download handler
    output$download_data <- downloadHandler(
      filename = function() {
        paste0("coloc_", gsub(":", "_", current_region()), "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        fwrite(filtered_data(), file)
      }
    )

    # === Regional Plot Logic (Region View) ===

    # Path to regional plot data
    # Supports three layouts:
    #   1. Category: {DATA_PATH}/<category>/regional/<study>/   (multi-category atlas)
    #   2. Bundled:  {DATA_PATH}/regional/<study>/              (single-category)
    #   3. Legacy:   {DATA_PATH}/regional_plots/<study>/        (old layout)
    regional_data_path <- reactive({
      req(current_study())
      study <- current_study()

      # Virtual multi-ancestry study: return the first ancestry regional dir
      # that actually contains the selected representative region's RDS
      # file (so `load_region_bundle` picks the right one). Fall back to the
      # first existing dir. Step 2 will generalize this to a list for
      # overlay rendering.
      if (is_virtual_study(study)) {
        vdirs <- DEFAULT_VIRTUAL_STUDIES[[study]]$regional_dirs
        vdirs <- vdirs[!is.na(unlist(vdirs))]
        if (length(vdirs) == 0) return(NULL)

        sel <- current_region()
        if (!is.null(sel) && !is.na(sel)) {
          rid <- parse_region_id(sel)
          if (!is.null(rid)) {
            for (d in unlist(vdirs)) {
              if (file.exists(file.path(d, paste0(rid, ".RDS")))) return(d)
            }
          }
        }
        return(unlist(vdirs)[1])
      }

      # Look up category from discover_studies() attribute
      cats <- attr(available_studies, "categories")
      if (!is.null(cats) && !is.null(cats[[study]])) {
        cat_path <- file.path(DATA_PATH, cats[[study]], "regional", study)
        if (dir.exists(cat_path)) return(cat_path)
      }

      # Fall back to single-category bundled layout
      bundled_path <- file.path(DATA_PATH, "regional", study)
      if (dir.exists(bundled_path)) return(bundled_path)

      # Fall back to legacy layout
      legacy_path <- file.path(DATA_PATH, "regional_plots", study)
      if (dir.exists(legacy_path)) return(legacy_path)

      NULL
    })

    # Detect layout: bundled (one RDS per region) vs legacy (directory per region)
    regional_layout <- reactive({
      req(regional_data_path())
      # Bundled layout has .RDS files directly in the trait folder
      rds_files <- list.files(regional_data_path(), pattern = "\\.RDS$", full.names = FALSE)
      if (length(rds_files) > 0) "bundled" else "legacy"
    })

    # Get available traits for current region
    regional_traits_available <- reactive({
      req(regional_data_path(), regional_layout(), current_region())

      if (regional_layout() == "bundled") {
        is_virt <- !is.null(current_study()) && is_virtual_study(current_study()) &&
                   "ancestry" %in% names(coloc_data())

        # --- Virtual multi-ancestry: consensus-keyed trait dropdown ---
        if (is_virt) {
          bundles <- load_multi_region_bundles(current_region())
          if (is.null(bundles) || length(bundles) == 0) return(NULL)
          cons <- bundles[[1]]$.consensus
          if (is.null(cons)) return(NULL)

          # Pull colocs for ANY per-ancestry region that overlaps the
          # consensus window (not filtered_data(), which is locked to the
          # representative per-ancestry region and would hide other
          # ancestries' colocs).
          all_dt <- coloc_data()
          dt <- all_dt[as.character(CHR_var) == as.character(cons$chr) &
                       BP_START_var <= cons$stop &
                       BP_STOP_var  >= cons$start]
          if (nrow(dt) == 0) return(NULL)
          dt <- dt[PP.H4.abf >= input$min_pph4]
          dt <- dt[sumstats_2_max_nlog10P >= input$min_nlog10P]
          dt <- dt[source_study %in% selected_studies()]
          if (nrow(dt) == 0) return(NULL)

          # Union of bundle keys across ancestries (each mapped to its
          # consensus key) tells us which consensus traits are available.
          all_keys <- unlist(lapply(bundles, function(b) names(b$traits)),
                             use.names = FALSE)
          if (length(all_keys) == 0) return(NULL)
          consensus_by_key <- make_consensus_trait_key(all_keys)

          # Build a consensus row per (source_study, sumstats_2_file)
          # collapsing across ancestries and recording which ancestries
          # colocalized.
          trait_info <- unique(dt[, .(
            source_study, ancestry,
            trait_display = tryCatch(get_trait_display_name(.SD), error = function(e) "Unknown"),
            ancestry_key  = make_bundle_key(source_study, sumstats_2_file, ancestry),
            consensus_key = make_consensus_trait_key(
              make_bundle_key(source_study, sumstats_2_file, ancestry))
          ), by = seq_len(nrow(dt))])

          # Only keep rows whose per-ancestry key exists in at least one
          # loaded bundle.
          trait_info <- trait_info[ancestry_key %in% all_keys]
          if (nrow(trait_info) == 0) return(NULL)

          consensus <- trait_info[, .(
            trait_display = trait_display[1],
            ancestries = paste(sort(unique(ancestry)), collapse = ","),
            n_anc = length(unique(ancestry))
          ), by = consensus_key]
          # Sort: more ancestries first, then alphabetical by display name.
          setorder(consensus, -n_anc, trait_display)

          choices <- setNames(
            consensus$consensus_key,
            paste0(consensus$trait_display, " [", consensus$ancestries, "]")
          )
          return(choices)
        }

        # --- Single-study bundled layout (original path) ---
        req(filtered_data())
        dt <- filtered_data()
        if (nrow(dt) == 0) return(NULL)

        region_data <- load_region_bundle(regional_data_path(), current_region())
        if (is.null(region_data) || length(region_data$traits) == 0) return(NULL)

        bundle_keys <- names(region_data$traits)

        trait_info <- unique(dt[, .(
          source_study,
          trait_display = tryCatch(get_trait_display_name(.SD), error = function(e) "Unknown"),
          bundle_key = make_bundle_key(source_study, sumstats_2_file)
        ), by = seq_len(nrow(dt))])

        trait_info <- trait_info[bundle_key %in% bundle_keys]
        if (nrow(trait_info) == 0) return(NULL)

        trait_info[, trait_label := paste0(trait_display, " (", source_study, ")")]
        setNames(trait_info$bundle_key, trait_info$trait_label)
      } else {
        # Legacy layout: use get_trait_name for filename matching
        req(filtered_data())
        dt <- filtered_data()
        if (nrow(dt) == 0) return(NULL)
        trait_info <- unique(dt[, .(
          source_study,
          trait_name = tryCatch(get_trait_name(.SD), error = function(e) "Unknown"),
          trait_display = tryCatch(get_trait_display_name(.SD), error = function(e) "Unknown")
        ), by = seq_len(nrow(dt))])
        trait_info[, trait_label := paste0(trait_display, " (", source_study, ")")]

        setNames(
          paste0(trait_info$source_study, "__", trait_info$trait_name),
          trait_info$trait_label
        )
      }
    })

    # Helper to extract trait name from filtered_data row
    # MUST match file naming in get_trait_name_for_regional() from ckdgen_genepicoloc_functions.R
    get_trait_name <- function(row) {
      source <- row$source_study[1]
      if (is.null(source) || is.na(source)) return("Unknown")

      trait_name <- switch(source,
        "eQTLGen" = {
          # eQTLGen uses gene_name_ENSG format
          if ("eQTLGen_gene_name" %in% names(row) && !is.na(row$eQTLGen_gene_name[1])) {
            # Try to extract ENSG from sumstats_2_file if available
            if ("sumstats_2_file" %in% names(row)) {
              ensg <- gsub(".*_(ENSG[0-9]+).*", "\\1", basename(row$sumstats_2_file[1]))
              paste0(row$eQTLGen_gene_name[1], "_", ensg)
            } else {
              row$eQTLGen_gene_name[1]
            }
          } else {
            "Unknown"
          }
        },
        "UKB_PPP_EUR" = {
          if ("UKB_PPP_EUR_olink_target_fullname" %in% names(row)) {
            row$UKB_PPP_EUR_olink_target_fullname[1]
          } else "Unknown"
        },
        "Icelanders_pGWAS" = {
          if ("Icelanders_pGWAS_Protein..short.name." %in% names(row)) {
            row$`Icelanders_pGWAS_Protein..short.name.`[1]
          } else "Unknown"
        },
        "GCKD_mGWAS_plasma" = {
          if ("GCKD_mGWAS_plasma_BIOCHEMICAL" %in% names(row)) {
            row$GCKD_mGWAS_plasma_BIOCHEMICAL[1]
          } else "Unknown"
        },
        "GCKD_mGWAS_urine" = {
          if ("GCKD_mGWAS_urine_BIOCHEMICAL" %in% names(row)) {
            row$GCKD_mGWAS_urine_BIOCHEMICAL[1]
          } else "Unknown"
        },
        "FinnGen_r9" = {
          if ("FinnGen_r9_phenotype" %in% names(row)) {
            row$FinnGen_r9_phenotype[1]
          } else "Unknown"
        },
        "UKB_TOPMed" = {
          if ("UKB_TOPMed_phenostring" %in% names(row)) {
            row$UKB_TOPMed_phenostring[1]
          } else "Unknown"
        },
        "MVP_R4" = {
          if ("MVP_R4_Title.of.analysis" %in% names(row)) {
            row$`MVP_R4_Title.of.analysis`[1]
          } else "Unknown"
        },
        "CKDGen_r4" = {
          if ("CKDGen_r4_Name" %in% names(row)) {
            row$CKDGen_r4_Name[1]
          } else "Unknown"
        },
        "CKDGen_r5" = {
          if ("CKDGen_r5_ckdgen_r5_name" %in% names(row)) {
            row$CKDGen_r5_ckdgen_r5_name[1]
          } else "Unknown"
        },
        "pho_ca" = {
          if ("pho_ca_pho_ca_name" %in% names(row)) {
            row$pho_ca_pho_ca_name[1]
          } else "Unknown"
        },
        "MVP_R4_EUR" = , "MVP_R4_AFR" = , "MVP_R4_AMR" = , "MVP_R4_EAS" = , "MVP_R4_META" = {
          col <- paste0(source, "_Title.of.analysis")
          if (col %in% names(row)) row[[col]][1] else "Unknown"
        },
        "UKB_kidney_vol" = {
          gsub("model1_qnorm_(.*)_chr.*", "\\1", basename(row$sumstats_2_file[1]))
        },
        "Kidney_eQTL" = {
          if ("Kidney_eQTL_gene_name" %in% names(row)) {
            row$Kidney_eQTL_gene_name[1]
          } else "Unknown"
        },
        "Unknown"
      )

      # Clean trait name for filename safety (match extraction function)
      if (!is.null(trait_name) && !is.na(trait_name) && trait_name != "Unknown") {
        trait_name <- gsub("[^A-Za-z0-9_.-]", "_", trait_name)
        trait_name <- gsub("_+", "_", trait_name)
        trait_name <- gsub("^_|_$", "", trait_name)
      }

      return(trait_name)
    }

    # Helper to get human-readable display name (for dropdown labels)
    # For most studies, same as get_trait_name; for MVP_R4, use full phecode description
    get_trait_display_name <- function(row) {
      source <- row$source_study[1]
      if (is.null(source) || is.na(source)) return("Unknown")

      # MVP_R4: Use Analyzed.variable (full phecode name) instead of Title.of.analysis (phecode)
      if (source == "MVP_R4") {
        if ("MVP_R4_Analyzed.variable" %in% names(row) && !is.na(row$`MVP_R4_Analyzed.variable`[1])) {
          display_name <- row$`MVP_R4_Analyzed.variable`[1]
          # Clean up: replace dots with spaces for readability
          display_name <- gsub("\\.", " ", display_name)
          display_name <- gsub("  +", " ", display_name)
          display_name <- trimws(display_name)
          return(display_name)
        }
      }

      # For all other studies, use the same as get_trait_name (already clean)
      get_trait_name(row)
    }

    # Update trait selector when region changes
    observe({
      traits <- tryCatch(regional_traits_available(), error = function(e) NULL)

      if (is.null(traits) || length(traits) == 0) {
        updateSelectInput(session, "regional_trait_selector",
                         choices = c("No colocalized traits" = ""),
                         selected = "")
      } else {
        updateSelectInput(session, "regional_trait_selector",
                         choices = traits,
                         selected = traits[1])
      }
    })

    # Helper: parse region string "chr:start-end" into folder name "chr_start_end"
    parse_region_id <- function(region_str) {
      region_parts <- strsplit(region_str, ":")[[1]]
      if (length(region_parts) < 2) return(NULL)
      chr <- region_parts[1]
      if (!grepl("^chr", chr)) chr <- paste0("chr", chr)
      pos_parts <- strsplit(region_parts[2], "-")[[1]]
      paste0(chr, "_", pos_parts[1], "_", pos_parts[2])
    }

    # Helper: reconstruct Name column if missing (slim format)
    reconstruct_name <- function(dt) {
      if (!"Name" %in% names(dt) && all(c("CHR", "POS", "A1", "A2") %in% names(dt))) {
        dt$Name <- paste0("chr", dt$CHR, ":", dt$POS, ":", dt$A1, ":", dt$A2)
      }
      dt
    }

    # Load bundled region data (one RDS per region)
    # Returns list(base = dt, traits = list("study__trait" = dt, ...)) or NULL
    # Cached per region to avoid re-reading
    loaded_region_cache <- reactiveVal(list(id = NULL, data = NULL))

    load_region_bundle <- function(data_path, region_str) {
      region_id <- parse_region_id(region_str)
      if (is.null(region_id)) return(NULL)

      # Check cache
      cache <- loaded_region_cache()
      cache_key <- paste(data_path, region_id)
      if (identical(cache$id, cache_key)) return(cache$data)

      rds_file <- file.path(data_path, paste0(region_id, ".RDS"))
      if (!file.exists(rds_file)) return(NULL)

      region_data <- readRDS(rds_file)

      # Reconstruct Name in base and all traits
      if (!is.null(region_data$base)) {
        region_data$base <- reconstruct_name(region_data$base)
      }
      region_data$traits <- lapply(region_data$traits, reconstruct_name)

      loaded_region_cache(list(id = cache_key, data = region_data))
      region_data
    }

    # Load base sumstats: supports bundled and legacy layouts
    load_base_sumstats <- function(data_path, region_str, layout) {
      if (layout == "bundled") {
        region_data <- load_region_bundle(data_path, region_str)
        if (is.null(region_data)) return(NULL)
        return(region_data$base)
      }
      # Legacy layout
      region_id <- parse_region_id(region_str)
      sumstats_file <- file.path(data_path, "regions", region_id, "sumstats_1.RDS")
      if (!file.exists(sumstats_file)) return(NULL)
      reconstruct_name(readRDS(sumstats_file))
    }

    # Parse a regional bundle filename (chr16_19330554_21586583.RDS) back
    # into numeric coordinates. Returns list(chr, start, stop) or NULL.
    parse_regional_filename <- function(fname) {
      m <- regmatches(fname,
        regexec("^chr([0-9XY]+)_([0-9]+)_([0-9]+)\\.RDS$", fname))[[1]]
      if (length(m) != 4) return(NULL)
      list(chr = m[2], start = as.numeric(m[3]), stop = as.numeric(m[4]))
    }

    # For a virtual multi-ancestry study, find all per-ancestry region RDS
    # files that overlap the consensus cluster containing region_str. The
    # consensus coordinates are derived by scanning all ancestries' region
    # files once and merging intervals, seeded from region_str. Returns a
    # named list: ancestry -> region_data (list(base, traits)), cached to
    # avoid re-reading on every render.
    loaded_multi_cache <- reactiveVal(list(id = NULL, data = NULL))

    load_multi_region_bundles <- function(region_str) {
      study <- current_study()
      if (is.null(study) || !is_virtual_study(study)) return(NULL)
      vinfo <- DEFAULT_VIRTUAL_STUDIES[[study]]

      # Seed chromosome + coordinates from region_str (representative key).
      rp <- strsplit(region_str, ":")[[1]]
      if (length(rp) != 2) return(NULL)
      chr <- sub("^chr", "", rp[1])
      pp <- as.numeric(strsplit(rp[2], "-")[[1]])
      seed_start <- pp[1]; seed_stop <- pp[2]

      # Check cache
      cache_key <- paste(study, chr, seed_start, seed_stop)
      cache <- loaded_multi_cache()
      if (identical(cache$id, cache_key)) return(cache$data)

      # Collect all candidate files across ancestries on this chromosome.
      candidates <- list()  # rows: ancestry, file, start, stop
      for (anc in names(vinfo$regional_dirs)) {
        d <- vinfo$regional_dirs[[anc]]
        if (is.na(d) || !dir.exists(d)) next
        fs <- list.files(d, pattern = "\\.RDS$", full.names = FALSE)
        for (f in fs) {
          coord <- parse_regional_filename(f)
          if (is.null(coord)) next
          if (sub("^chr", "", coord$chr) != chr) next
          candidates[[length(candidates) + 1]] <- list(
            ancestry = anc, file = file.path(d, f),
            start = coord$start, stop = coord$stop)
        }
      }
      if (length(candidates) == 0) return(NULL)
      cand_dt <- data.table::rbindlist(candidates)

      # Grow the consensus window iteratively: start from the seed and
      # include any candidate overlapping the running window; repeat until
      # stable.
      cur_start <- seed_start; cur_stop <- seed_stop
      repeat {
        hit <- cand_dt[start <= cur_stop & stop >= cur_start]
        if (nrow(hit) == 0) break
        new_start <- min(hit$start); new_stop <- max(hit$stop)
        if (new_start == cur_start && new_stop == cur_stop) break
        cur_start <- new_start; cur_stop <- new_stop
      }
      hit <- cand_dt[start <= cur_stop & stop >= cur_start]
      if (nrow(hit) == 0) return(NULL)

      # Within each ancestry keep the widest overlapping file (usually
      # there's only one; if multiple we prefer the largest window).
      hit[, width := stop - start]
      setorder(hit, ancestry, -width)
      hit <- hit[, .SD[1], by = ancestry]

      # Read and reconstruct each bundle.
      result <- list()
      for (i in seq_len(nrow(hit))) {
        rd <- tryCatch(readRDS(hit$file[i]), error = function(e) NULL)
        if (is.null(rd)) next
        if (!is.null(rd$base)) rd$base <- reconstruct_name(rd$base)
        rd$traits <- lapply(rd$traits, reconstruct_name)
        rd$.consensus <- list(chr = chr, start = cur_start, stop = cur_stop)
        result[[hit$ancestry[i]]] <- rd
      }

      loaded_multi_cache(list(id = cache_key, data = result))
      result
    }

    # Load base study sumstats for current region. For virtual studies,
    # returns a named list(ancestry = dt) for overlay rendering; otherwise
    # returns a single data.table.
    regional_base_sumstats <- reactive({
      req(regional_data_path(), regional_layout(), current_region())
      study <- current_study()
      if (!is.null(study) && is_virtual_study(study)) {
        bundles <- load_multi_region_bundles(current_region())
        if (is.null(bundles) || length(bundles) == 0) return(NULL)
        return(lapply(bundles, `[[`, "base"))
      }
      load_base_sumstats(regional_data_path(), current_region(), regional_layout())
    })

    # Load trait sumstats: supports bundled and legacy layouts
    load_trait_sumstats <- function(data_path, region_str, trait_selector, layout) {
      if (layout == "bundled") {
        region_data <- load_region_bundle(data_path, region_str)
        if (is.null(region_data)) return(NULL)
        return(region_data$traits[[trait_selector]])
      }
      # Legacy layout
      region_id <- parse_region_id(region_str)
      trait_file <- file.path(data_path, "regions", region_id, "traits",
                              paste0(trait_selector, ".RDS"))
      if (!file.exists(trait_file)) return(NULL)

      dt <- reconstruct_name(readRDS(trait_file))
      # Filter to region bounds (legacy files may span multiple regions)
      region_parts <- strsplit(region_str, ":")[[1]]
      pos_parts <- strsplit(region_parts[2], "-")[[1]]
      if (is.data.frame(dt) && nrow(dt) > 0 && "POS" %in% names(dt)) {
        dt <- dt[dt$POS >= as.numeric(pos_parts[1]) & dt$POS <= as.numeric(pos_parts[2]), ]
      }
      dt
    }

    # Load trait sumstats for selected colocalized trait. For virtual
    # multi-ancestry studies the selector is a consensus key; resolve it
    # against each per-ancestry region bundle and return a named list.
    regional_trait_sumstats <- reactive({
      req(regional_data_path(), regional_layout(), input$regional_trait_selector, current_region())
      sel <- input$regional_trait_selector
      if (is.null(sel) || sel == "") return(NULL)

      if (!is.null(current_study()) && is_virtual_study(current_study())) {
        bundles <- load_multi_region_bundles(current_region())
        if (is.null(bundles) || length(bundles) == 0) return(NULL)
        out <- lapply(names(bundles), function(anc) {
          bk <- names(bundles[[anc]]$traits)
          ak <- resolve_consensus_in_bundle(sel, bk, anc)
          if (is.null(ak)) return(NULL)
          bundles[[anc]]$traits[[ak]]
        })
        names(out) <- names(bundles)
        out <- out[!vapply(out, is.null, logical(1))]
        if (length(out) == 0) return(NULL)
        return(out)
      }

      load_trait_sumstats(regional_data_path(), current_region(),
                          sel, regional_layout())
    })

    # Get lead SNP for current region
    regional_lead_snp <- reactive({
      req(filtered_data())
      dt <- filtered_data()

      if (nrow(dt) > 0 && "sumstats_1_ind_Name" %in% names(dt)) {
        dt$sumstats_1_ind_Name[1]
      } else {
        NULL
      }
    })

    # Render base study plot title
    output$regional_plot_title_base <- renderText({
      req(current_study())
      paste0("Base Study: ", current_study())
    })

    # Render trait plot title
    output$regional_plot_title_trait <- renderText({
      req(input$regional_trait_selector)
      if (input$regional_trait_selector == "") return("Select a trait")
      # Get display name from available traits (the label, not the value)
      traits <- regional_traits_available()
      if (!is.null(traits) && input$regional_trait_selector %in% traits) {
        names(traits)[which(traits == input$regional_trait_selector)]
      } else {
        gsub("__", " - ", input$regional_trait_selector)
      }
    })

    # Render coloc info
    output$regional_coloc_info <- renderText({
      req(filtered_data(), input$regional_trait_selector)
      dt <- filtered_data()

      if (input$regional_trait_selector == "" || nrow(dt) == 0) {
        return("")
      }

      # Find PP.H4 for selected trait
      parts <- strsplit(input$regional_trait_selector, "__")[[1]]
      if (length(parts) == 2) {
        trait_rows <- dt[source_study == parts[1]]
        if (nrow(trait_rows) > 0) {
          max_pph4 <- max(trait_rows$PP.H4.abf, na.rm = TRUE)
          return(paste0("PP.H4: ", round(max_pph4, 3)))
        }
      }
      ""
    })

    # Calculate shared X-axis range for Region View plots
    regional_x_range <- reactive({
      # Virtual multi-ancestry: use the consensus window from the loaded
      # multi bundles (widest overlap across ancestries).
      req(current_region())
      if (!is.null(current_study()) && is_virtual_study(current_study())) {
        bundles <- load_multi_region_bundles(current_region())
        cons <- if (!is.null(bundles) && length(bundles) > 0)
          bundles[[1]]$.consensus else NULL
        if (!is.null(cons)) return(c(cons$start, cons$stop))
      }

      # Use region bounds as x_range (not data extent)
      region_parts <- strsplit(current_region(), ":")[[1]]
      if (length(region_parts) >= 2) {
        pos_parts <- as.numeric(strsplit(region_parts[2], "-")[[1]])
        return(pos_parts)
      }

      # Fallback: use data extent
      base_sumstats <- regional_base_sumstats()
      trait_sumstats <- regional_trait_sumstats()

      # Handle overlay list payloads.
      gather_pos <- function(x) {
        if (is.null(x)) return(numeric(0))
        if (is.data.frame(x)) return(x$POS)
        if (is.list(x)) return(unlist(lapply(x, function(y) if (is.data.frame(y)) y$POS else numeric(0))))
        numeric(0)
      }
      all_pos <- c(gather_pos(base_sumstats), gather_pos(trait_sumstats))

      if (length(all_pos) > 0) {
        # Add small padding (1% on each side)
        range_span <- max(all_pos) - min(all_pos)
        padding <- range_span * 0.01
        return(c(min(all_pos) - padding, max(all_pos) + padding))
      }
      NULL
    })

    # Render base study regional plot (interactive)
    output$regional_plot_base <- renderPlotly({
      sumstats <- regional_base_sumstats()
      lead_snp <- if (input$regional_highlight_lead) regional_lead_snp() else NULL
      x_range <- regional_x_range()

      plot_regional_association_interactive(
        sumstats_dt = sumstats,
        title = "",
        highlight_snp = lead_snp,
        color = "#2ecc71",  # Green for base study
        x_range = x_range
      )
    })

    # Render trait regional plot (interactive)
    output$regional_plot_trait <- renderPlotly({
      sumstats <- regional_trait_sumstats()
      lead_snp <- if (input$regional_highlight_lead) regional_lead_snp() else NULL
      x_range <- regional_x_range()

      # Get color based on study type
      trait_sel <- input$regional_trait_selector
      study_color <- "#3498db"  # Default blue
      if (!is.null(trait_sel) && trait_sel != "") {
        study <- strsplit(trait_sel, "__")[[1]][1]
        if (study %in% names(study_colors)) {
          study_color <- study_colors[[study]]
        }
      }

      plot_regional_association_interactive(
        sumstats_dt = sumstats,
        title = "",
        highlight_snp = lead_snp,
        color = study_color,
        x_range = x_range
      )
    })

    # Render gene track (Region View)
    output$regional_gene_track <- renderPlotly({
      req(current_region())

      # Parse region coordinates from current_region()
      region_sel <- current_region()
      region_parts <- strsplit(region_sel, ":")[[1]]
      chr <- region_parts[1]
      pos_parts <- strsplit(region_parts[2], "-")[[1]]
      start_pos <- as.numeric(pos_parts[1])
      end_pos <- as.numeric(pos_parts[2])

      plot_gene_track(chr, start_pos, end_pos, source = "regional_gene_track")
    })

    # Store selected gene for info panel
    selected_gene <- reactiveVal(NULL)

    # Handle click on gene track (Region View)
    observeEvent(event_data("plotly_click", source = "regional_gene_track"), {
      click_data <- event_data("plotly_click", source = "regional_gene_track")
      if (!is.null(click_data) && !is.null(click_data$customdata)) {
        clicked_gene <- click_data$customdata
        # Look up gene in annotation
        if (!is.null(gene_annotation)) {
          gene_info <- gene_annotation[gene_name == clicked_gene]
          if (nrow(gene_info) > 0) {
            selected_gene(gene_info[1])
          }
        }
      }
    })

    # Render gene info panel
    output$gene_info_panel <- renderUI({
      gene <- selected_gene()
      if (is.null(gene)) return(NULL)

      # Helper for source tooltip (ⓘ with hover description, not clickable)
      source_tip <- function(source_name, desc, source_url) {
        tags$span(
          title = paste0(source_name, " - ", desc),
          style = "color: #999; margin-left: 4px; cursor: help;",
          "ⓘ"
        )
      }

      # Build info sections
      info_items <- list()

      # Gene name and full name + close button (right-aligned)
      info_items[[length(info_items) + 1]] <- div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        tags$h4(style = "margin: 0;",
          gene$gene_name,
          if (!is.null(gene$full_name) && gene$full_name != "")
            tags$span(style = "color: #666; margin-left: 10px; font-weight: normal;", gene$full_name)
        ),
        actionButton("gene_info_close", "Close",
                     icon = icon("times"),
                     class = "btn-xs btn-default")
      )

      # Coordinates (from Gencode)
      info_items[[length(info_items) + 1]] <- tags$p(
        tags$strong("Location: "),
        paste0("chr", gene$chr, ":", format(gene$start, big.mark = ","), "-", format(gene$end, big.mark = ","),
               " (", gene$strand, " strand)"),
        source_tip("Gencode v49", "comprehensive gene annotation", "https://www.gencodegenes.org/")
      )

      # Gene family (from HGNC)
      if (!is.null(gene$gene_family) && gene$gene_family != "") {
        info_items[[length(info_items) + 1]] <- tags$p(
          tags$strong("Gene Family: "), gene$gene_family,
          source_tip("HGNC", "HUGO Gene Nomenclature Committee", "https://www.genenames.org/")
        )
      }

      # NCBI description
      if (!is.null(gene$ncbi_description) && gene$ncbi_description != "") {
        info_items[[length(info_items) + 1]] <- tags$p(
          tags$strong("Description: "), gene$ncbi_description,
          source_tip("NCBI Gene", "gene-specific information", "https://www.ncbi.nlm.nih.gov/gene/")
        )
      }

      # Pathways (from Reactome, top 5)
      if (!is.null(gene$pathways) && gene$pathways != "") {
        info_items[[length(info_items) + 1]] <- tags$p(
          tags$strong("Pathways (top 5): "), gene$pathways,
          source_tip("Reactome", "curated biological pathways", "https://reactome.org/")
        )
      }

      # Phenotypes (from HPO, top 5)
      if (!is.null(gene$phenotypes) && gene$phenotypes != "") {
        info_items[[length(info_items) + 1]] <- tags$p(
          tags$strong("HPO Phenotypes (top 5): "), gene$phenotypes,
          source_tip("HPO", "Human Phenotype Ontology", "https://hpo.jax.org/")
        )
      }

      # External links
      info_items[[length(info_items) + 1]] <- tags$p(
        tags$a(href = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene$gene_name),
               target = "_blank", "GeneCards"),
        tags$span(" | ", style = "color: #999;"),
        tags$a(href = paste0("https://www.genenames.org/tools/search/#!/?query=", gene$gene_name),
               target = "_blank", "HGNC"),
        tags$span(" | ", style = "color: #999;"),
        tags$a(href = paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", gene$gene_name),
               target = "_blank", "NCBI"),
        tags$span(" | ", style = "color: #999;"),
        tags$a(href = paste0("https://hpo.jax.org/browse/search?q=", gene$gene_name, "&navFilter=all"),
               target = "_blank", "HPO")
      )

      # Wrap in a styled panel
      wellPanel(
        style = "background-color: #f8f9fa; border: 1px solid #dee2e6; margin-top: 10px;",
        do.call(tagList, info_items)
      )
    })

    # Close gene info panel
    observeEvent(input$gene_info_close, {
      selected_gene(NULL)
    })

  }

# Run the application
shinyApp(ui = ui, server = server)
