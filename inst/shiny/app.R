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

        // Trait View network export
        Shiny.addCustomMessageHandler('exportTraitNetwork', function(message) {
          var container = document.getElementById('trait_network');
          if (container) {
            var canvas = container.querySelector('canvas');
            if (canvas) {
              var link = document.createElement('a');
              link.download = message.filename;
              link.href = canvas.toDataURL('image/png');
              link.click();
            } else {
              alert('Trait network canvas not found. Please wait for the network to render.');
            }
          } else {
            alert('Trait network container not found.');
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
      tabPanel("Trait View",
               titlePanel("Trait-Centered Network View"),

               fluidRow(
                 # Left sidebar panel
                 column(
                   width = 3,
                   wellPanel(
                     # Study type selector
                     selectInput("trait_view_type",
                                 "Select Study Type:",
                                 choices = c("Phenotypes" = "Phenotypes",
                                           "pQTL" = "pQTL",
                                           "eQTL" = "eQTL",
                                           "mQTL" = "mQTL"),
                                 selected = "Phenotypes"),

                     hr(),

                     # Phenotypes study and trait selectors (conditional)
                     conditionalPanel(
                       condition = "input.trait_view_type == 'Phenotypes'",
                       selectInput("selected_pheno_study",
                                   "Select Study:",
                                   choices = c("MVP Million Veteran Program (MVP_R4)" = "MVP_R4",
                                             "FinnGen (FinnGen_r9)" = "FinnGen_r9",
                                             "UK Biobank TOPMed (UKB_TOPMed)" = "UKB_TOPMed"),
                                   selected = "MVP_R4"),

                       selectInput("selected_trait_pheno",
                                   "Select Trait/Phenotype:",
                                   choices = NULL,
                                   selected = NULL)
                     ),

                     # pQTL study and protein selectors (conditional)
                     conditionalPanel(
                       condition = "input.trait_view_type == 'pQTL'",
                       selectInput("selected_pqtl_study",
                                   "Select Study:",
                                   choices = c("UKB Plasma Proteomics (UKB_PPP_EUR)" = "UKB_PPP_EUR",
                                             "Icelanders Proteomics (Icelanders_pGWAS)" = "Icelanders_pGWAS"),
                                   selected = "UKB_PPP_EUR"),

                       selectInput("selected_protein_pqtl",
                                   "Select Protein:",
                                   choices = NULL,
                                   selected = NULL)
                     ),

                     # eQTL study and gene selectors (conditional)
                     conditionalPanel(
                       condition = "input.trait_view_type == 'eQTL'",
                       selectInput("selected_eqtl_study",
                                   "Select Study:",
                                   choices = c("eQTLGen Phase I (eQTLGen)" = "eQTLGen"
                                             # "Kidney eQTL (Kidney_eQTL)" = "Kidney_eQTL",  # Uncomment when data available
                                             # "GTEx V8 (GTEXv8_eQTL)" = "GTEXv8_eQTL"      # Uncomment when data available
                                             ),
                                   selected = "eQTLGen"),

                       selectInput("selected_gene_eqtl",
                                   "Select Gene:",
                                   choices = NULL,
                                   selected = NULL)
                     ),

                     # mQTL study and metabolite selectors (conditional)
                     conditionalPanel(
                       condition = "input.trait_view_type == 'mQTL'",
                       selectInput("selected_mqtl_study",
                                   "Select Study:",
                                   choices = c("GCKD mGWAS Plasma (GCKD_mGWAS_plasma)" = "GCKD_mGWAS_plasma",
                                             "GCKD mGWAS Urine (GCKD_mGWAS_urine)" = "GCKD_mGWAS_urine"),
                                   selected = "GCKD_mGWAS_plasma"),

                       selectInput("selected_metabolite_mqtl",
                                   "Select Metabolite:",
                                   choices = NULL,
                                   selected = NULL)
                     ),

                     hr(),

                     # PP.H4 filter (shared)
                     sliderInput("min_pph4_trait",
                                 "Minimum PP.H4:",
                                 min = 0.8, max = 1, value = 0.8, step = 0.01),

                     # p-value filter (shared)
                     sliderInput("min_nlog10P_trait",
                                 "Minimum -log10(p-value):",
                                 min = 7.308, max = 100, value = 7.308, step = 0.5),

                     hr(),

                     # Edge color legend
                     h5("Edge Colors (Directionality)"),
                     HTML('
                       <div style="margin: 2px 0;">
                         <span style="display: inline-block; width: 30px; height: 4px;
                           background-color: #e67e22; margin-right: 5px; vertical-align: middle;"></span>
                         <span style="font-size: 11px;">Concordant</span>
                       </div>
                       <div style="margin: 2px 0;">
                         <span style="display: inline-block; width: 30px; height: 4px;
                           background-color: #3498db; margin-right: 5px; vertical-align: middle;"></span>
                         <span style="font-size: 11px;">Discordant</span>
                       </div>
                       <div style="margin: 2px 0;">
                         <span style="display: inline-block; width: 30px; height: 4px;
                           background-color: #848484; margin-right: 5px; vertical-align: middle;"></span>
                         <span style="font-size: 11px;">Unknown/N/A</span>
                       </div>
                     '),

                     hr(),

                     # Display options
                     checkboxInput("trait_physics", "Enable physics simulation", value = FALSE),

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
                               column(4, actionButton("trait_nodes_all", "All", class = "btn-xs btn-block")),
                               column(4, actionButton("trait_nodes_none", "None", class = "btn-xs btn-block")),
                               column(4, actionButton("trait_nodes_top10", "Top 10", class = "btn-xs btn-block"))
                             ),
                             br(),
                             selectizeInput("trait_visible_nodes",
                                           label = NULL,
                                           choices = NULL,
                                           multiple = TRUE,
                                           options = list(
                                             plugins = list('remove_button'),
                                             placeholder = 'Select nodes to display...'
                                           )),

                             hr(style = "margin: 10px 0;"),

                             # Node colors
                             h6("Node Colors:", style = "margin-bottom: 5px;"),
                             div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                               textInput("trait_central_color", label = NULL,
                                 value = "#9b59b6"),
                               span("Central (Trait/Protein)", style = "font-size: 11px; margin-left: 5px;")
                             ),
                             div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                               textInput("trait_region_color", label = NULL,
                                 value = "#3498db"),
                               span("Regions (Genes)", style = "font-size: 11px; margin-left: 5px;")
                             ),

                             hr(style = "margin: 10px 0;"),

                             # Edit selected node label
                             h6("Edit Node Label:", style = "margin-bottom: 5px;"),
                             textInput("trait_custom_node_label", label = NULL,
                                      placeholder = "Click a node to edit its label"),
                             actionButton("trait_apply_label", "Apply Label",
                                         class = "btn-xs btn-success btn-block")
                           )
                         )
                       )
                     }
                   )  # End wellPanel
                 ),  # End left column

                 # Main panel (network + regional plot tabs)
                 column(
                   width = 9,
                   tabsetPanel(
                     id = "trait_view_tabs",
                     type = "tabs",

                     # Tab 1: Network visualization
                     tabPanel(
                       "Network",
                       br(),
                       # Export buttons row
                       fluidRow(
                         column(12,
                           div(style = "text-align: right; margin-bottom: 5px;",
                             actionButton("trait_export_png", "Export PNG",
                                         icon = icon("image"),
                                         class = "btn-sm btn-default",
                                         style = "margin-right: 5px;"),
                             downloadButton("trait_export_html", "Export HTML",
                                           class = "btn-sm btn-primary")
                           )
                         )
                       ),
                       # Network visualization
                       visNetworkOutput("trait_network", height = "600px"),
                       br(),
                       p("Center node (purple diamond) = Selected trait/protein/transcript/metabolite. Surrounding nodes (blue) = Colocalizing genomic regions/genes from the selected CKDGen trait.",
                         style = "font-size: 12px; color: #7f8c8d;")
                     ),  # End Network tab

                     # Tab 2: Regional Association Plot
                     tabPanel(
                       "Regional Plot",
                       br(),
                       fluidRow(
                         column(4,
                           selectInput("trait_regional_region_selector",
                                      "Select Colocalized Region:",
                                      choices = NULL)
                         ),
                         column(4,
                           checkboxInput("trait_regional_highlight_lead",
                                        "Highlight Lead SNP",
                                        value = TRUE)
                         ),
                         column(4,
                           div(style = "margin-top: 25px;",
                             textOutput("trait_regional_coloc_info")
                           )
                         )
                       ),
                       # Two stacked regional plots (interactive with plotly)
                       h5(textOutput("trait_regional_plot_title_base"), style = "margin-top: 10px;"),
                       plotlyOutput("trait_regional_plot_base", height = "300px"),
                       tags$p(style = "color: #888; font-size: 12px; margin: 5px 0;",
                              "Displayed loci are defined as 1 Mb windows around clumping index variants and may span multiple association peaks.",
                              tags$span(
                                title = "Displayed loci capture the local genetic architecture. Look for association peaks at similar positions in both traits - overlapping signals support true colocalization. See Documentation for additional interpretation guidance.",
                                style = "cursor: help; margin-left: 4px;",
                                "ⓘ"
                              )),
                       hr(),
                       h5(textOutput("trait_regional_plot_title_trait")),
                       plotlyOutput("trait_regional_plot_trait", height = "300px"),
                       hr(),
                       h5("Genes in Region", tags$span(style = "color: #666; font-weight: normal;", " - click gene for details")),
                       plotlyOutput("trait_regional_gene_track", height = "150px"),
                       uiOutput("trait_gene_info_panel")
                     )  # End Regional Plot tab
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
      )
    )

    # Build trait list for one landing card. If the card declares
    # `autopopulate = "<atlas_category>"`, expand its (empty) traits list to
    # all studies discovered under that atlas folder, labeled via the bundled
    # annotation DB. Used for urine metabolomics where hand-listing 1,409
    # metabolites is impractical.
    expand_card_traits <- function(cat) {
      if (length(cat$traits) > 0 || is.null(cat$autopopulate)) {
        return(cat$traits)
      }
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

    output$landing_categories <- renderUI({
      div(class = "category-grid",
        lapply(atlas_categories, function(cat) {
          cat$traits <- expand_card_traits(cat)
          available <- sapply(cat$traits, function(t) {
            !isTRUE(t$disabled) && t$id %in% names(DEFAULT_AVAILABLE_STUDIES)
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
                is_disabled <- isTRUE(trait$disabled) || !trait$id %in% names(DEFAULT_AVAILABLE_STUDIES)
                # Check regional data across all possible layouts
                cats <- attr(DEFAULT_AVAILABLE_STUDIES, "categories")
                cat_name <- if (!is.null(cats)) cats[[trait$id]] else NULL
                has_regional <- !is_disabled && (
                  (!is.null(cat_name) && dir.exists(file.path(DATA_PATH, cat_name, "regional", trait$id))) ||
                  dir.exists(file.path(DATA_PATH, "regional", trait$id)) ||
                  dir.exists(file.path(DATA_PATH, "regional_plots", trait$id))
                )
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

      # Get max nlog10P per region from sumstats_1
      region_stats <- dt[, .(max_nlog10P = max(sumstats_1_max_nlog10P, na.rm = TRUE),
                             n_coloc = .N),
                         by = .(CHR_var, BP_START_var, BP_STOP_var)]

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
      region_stats <- merge(region_stats,
        unique(regions[, .(CHR_var, BP_START_var, nearest_gene_1, Prioritized_Gene)]),
        by = c("CHR_var", "BP_START_var"), all.x = TRUE)
      region_stats[, gene := ifelse(!is.na(Prioritized_Gene), Prioritized_Gene, nearest_gene_1)]
      region_stats[, region_key := paste0(CHR_var, ":", BP_START_var, "-", BP_STOP_var)]

      # Alternating chr colors
      region_stats[, chr_color := ifelse(CHR_num %% 2 == 0, "#3498db", "#2c3e50")]

      # Highlight currently selected region (strips gene:: prefix internally)
      sel <- current_region()
      if (is.null(sel)) sel <- NA_character_
      region_stats[, is_selected := !is.na(sel) & region_key == sel]

      # Hover text
      region_stats[, hover := paste0(
        "<b>", gene, "</b><br>",
        "chr", CHR_var, ":", format(BP_START_var, big.mark = ","), "<br>",
        "-log10(P): ", round(max_nlog10P, 1), "<br>",
        n_coloc, " colocalizations"
      )]

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
        file_path <- resolve_annot_path(current_study())
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
            paste("Successfully loaded:", current_study()),
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
                         "Prioritized_Gene", "clump_index_Name")
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

      dt <- coloc_data()[CHR_var == chr & BP_START_var == start_pos & BP_STOP_var == end_pos]
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

    # === Unified Trait View Reactives ===

    # Get trait data based on selected study type
    trait_data <- reactive({
      req(coloc_data(), input$trait_view_type)
      dt <- coloc_data()

      if (input$trait_view_type == "Phenotypes") {
        req(input$selected_pheno_study)
        dt <- dt[source_study == input$selected_pheno_study]
      } else if (input$trait_view_type == "pQTL") {
        req(input$selected_pqtl_study)
        dt <- dt[source_study == input$selected_pqtl_study]
      } else if (input$trait_view_type == "eQTL") {
        req(input$selected_eqtl_study)
        dt <- dt[source_study == input$selected_eqtl_study]
      } else if (input$trait_view_type == "mQTL") {
        req(input$selected_mqtl_study)
        dt <- dt[source_study == input$selected_mqtl_study]
      }

      dt
    })

    # Update Phenotypes trait selector
    observe({
      req(trait_data(), input$trait_view_type, input$selected_pheno_study)
      if (input$trait_view_type == "Phenotypes") {
        dt <- trait_data()

        trait_col <- switch(input$selected_pheno_study,
                           "MVP_R4" = "MVP_R4_Analyzed.variable",
                           "MVP_R4_EUR" = "MVP_R4_EUR_Analyzed.variable",
                           "FinnGen_r9" = "FinnGen_r9_phenotype",
                           "UKB_TOPMed" = "UKB_TOPMed_phenostring",
                           "CKDGen_r4" = "CKDGen_r4_Name",
                           NULL)

        if (!is.null(trait_col) && nrow(dt) > 0 && trait_col %in% names(dt)) {
          # Count colocs per trait
          trait_counts <- dt[, .N, by = trait_col]
          setnames(trait_counts, trait_col, "trait")
          trait_counts <- trait_counts[!is.na(trait)]
          trait_counts <- trait_counts[order(trait)]  # Sort alphabetically

          # Create named choices with counts
          choices <- setNames(
            trait_counts$trait,
            paste0(trait_counts$trait, " (", trait_counts$N, ")")
          )
          updateSelectInput(session, "selected_trait_pheno", choices = choices)
        }
      }
    })

    # Update pQTL protein selector
    observe({
      req(trait_data(), input$trait_view_type, input$selected_pqtl_study)
      if (input$trait_view_type == "pQTL") {
        dt <- trait_data()

        protein_col <- switch(input$selected_pqtl_study,
                             "UKB_PPP_EUR" = "UKB_PPP_EUR_olink_target_fullname",
                             "Icelanders_pGWAS" = "Icelanders_pGWAS_Protein..short.name.",
                             NULL)

        if (!is.null(protein_col) && nrow(dt) > 0 && protein_col %in% names(dt)) {
          # Count colocs per protein
          protein_counts <- dt[, .N, by = protein_col]
          setnames(protein_counts, protein_col, "protein")
          protein_counts <- protein_counts[!is.na(protein)]
          protein_counts <- protein_counts[order(protein)]  # Sort alphabetically

          # Create named choices with counts
          choices <- setNames(
            protein_counts$protein,
            paste0(protein_counts$protein, " (", protein_counts$N, ")")
          )
          updateSelectInput(session, "selected_protein_pqtl", choices = choices)
        }
      }
    })

    # Update eQTL gene selector
    observe({
      req(trait_data(), input$trait_view_type, input$selected_eqtl_study)
      if (input$trait_view_type == "eQTL") {
        dt <- trait_data()

        gene_col <- switch(input$selected_eqtl_study,
                          "eQTLGen" = "eQTLGen_gene_name",
                          # "Kidney_eQTL" = "Kidney_eQTL_gene_name",    # Uncomment when data available
                          # "GTEXv8_eQTL" = "GTEXv8_eQTL_gene_name",    # Uncomment when data available
                          NULL)

        if (!is.null(gene_col) && nrow(dt) > 0 && gene_col %in% names(dt)) {
          # Count colocs per gene
          gene_counts <- dt[, .N, by = gene_col]
          setnames(gene_counts, gene_col, "gene")
          gene_counts <- gene_counts[!is.na(gene)]
          gene_counts <- gene_counts[order(gene)]  # Sort alphabetically

          # Create named choices with counts
          choices <- setNames(
            gene_counts$gene,
            paste0(gene_counts$gene, " (", gene_counts$N, ")")
          )
          updateSelectInput(session, "selected_gene_eqtl", choices = choices)
        }
      }
    })

    # Update mQTL metabolite selector
    observe({
      req(trait_data(), input$trait_view_type, input$selected_mqtl_study)
      if (input$trait_view_type == "mQTL") {
        dt <- trait_data()

        metabolite_col <- switch(input$selected_mqtl_study,
                                "GCKD_mGWAS_plasma" = "GCKD_mGWAS_plasma_BIOCHEMICAL",
                                "GCKD_mGWAS_urine" = "GCKD_mGWAS_urine_BIOCHEMICAL",
                                NULL)

        if (!is.null(metabolite_col) && nrow(dt) > 0 && metabolite_col %in% names(dt)) {
          # Count colocs per metabolite
          metabolite_counts <- dt[, .N, by = metabolite_col]
          setnames(metabolite_counts, metabolite_col, "metabolite")
          metabolite_counts <- metabolite_counts[!is.na(metabolite)]
          metabolite_counts <- metabolite_counts[order(metabolite)]  # Sort alphabetically

          # Create named choices with counts
          choices <- setNames(
            metabolite_counts$metabolite,
            paste0(metabolite_counts$metabolite, " (", metabolite_counts$N, ")")
          )
          updateSelectInput(session, "selected_metabolite_mqtl", choices = choices)
        }
      }
    })

    # Filter by selected trait/protein/gene/metabolite, PP.H4, and p-value
    trait_filtered_data <- reactive({
      req(trait_data(), input$trait_view_type)
      dt <- trait_data()

      # Filter by trait/protein/gene/metabolite based on study type
      if (input$trait_view_type == "Phenotypes") {
        req(input$selected_trait_pheno, input$selected_pheno_study)

        trait_col <- switch(input$selected_pheno_study,
                           "MVP_R4" = "MVP_R4_Analyzed.variable",
                           "MVP_R4_EUR" = "MVP_R4_EUR_Analyzed.variable",
                           "FinnGen_r9" = "FinnGen_r9_phenotype",
                           "UKB_TOPMed" = "UKB_TOPMed_phenostring",
                           "CKDGen_r4" = "CKDGen_r4_Name",
                           NULL)

        if (!is.null(trait_col) && trait_col %in% names(dt)) {
          dt <- dt[get(trait_col) == input$selected_trait_pheno]
        }
      } else if (input$trait_view_type == "pQTL") {
        req(input$selected_protein_pqtl, input$selected_pqtl_study)

        protein_col <- switch(input$selected_pqtl_study,
                             "UKB_PPP_EUR" = "UKB_PPP_EUR_olink_target_fullname",
                             "Icelanders_pGWAS" = "Icelanders_pGWAS_Protein..short.name.",
                             NULL)

        if (!is.null(protein_col) && protein_col %in% names(dt)) {
          dt <- dt[get(protein_col) == input$selected_protein_pqtl]
        }
      } else if (input$trait_view_type == "eQTL") {
        req(input$selected_gene_eqtl, input$selected_eqtl_study)

        gene_col <- switch(input$selected_eqtl_study,
                          "eQTLGen" = "eQTLGen_gene_name",
                          # "Kidney_eQTL" = "Kidney_eQTL_gene_name",    # Uncomment when data available
                          # "GTEXv8_eQTL" = "GTEXv8_eQTL_gene_name",    # Uncomment when data available
                          NULL)

        if (!is.null(gene_col) && gene_col %in% names(dt)) {
          dt <- dt[get(gene_col) == input$selected_gene_eqtl]
        }
      } else if (input$trait_view_type == "mQTL") {
        req(input$selected_metabolite_mqtl, input$selected_mqtl_study)

        metabolite_col <- switch(input$selected_mqtl_study,
                                "GCKD_mGWAS_plasma" = "GCKD_mGWAS_plasma_BIOCHEMICAL",
                                "GCKD_mGWAS_urine" = "GCKD_mGWAS_urine_BIOCHEMICAL",
                                NULL)

        if (!is.null(metabolite_col) && metabolite_col %in% names(dt)) {
          dt <- dt[get(metabolite_col) == input$selected_metabolite_mqtl]
        }
      }

      # Apply filters (shared for all types)
      dt <- dt[PP.H4.abf >= input$min_pph4_trait]

      if ("sumstats_2_max_nlog10P" %in% names(dt)) {
        dt <- dt[sumstats_2_max_nlog10P >= input$min_nlog10P_trait]
      }

      dt
    })

    # Create trait-centered network (unified)
    trait_network_data <- reactive({
      req(trait_filtered_data(), input$trait_view_type)
      dt <- trait_filtered_data()

      if (nrow(dt) == 0) {
        return(list(nodes = data.frame(), edges = data.frame()))
      }

      # Determine label and study info based on type
      if (input$trait_view_type == "Phenotypes") {
        req(input$selected_trait_pheno, input$selected_pheno_study)
        node_label <- input$selected_trait_pheno
        study_info <- input$selected_pheno_study
        node_type <- "Trait"
      } else if (input$trait_view_type == "pQTL") {
        req(input$selected_protein_pqtl, input$selected_pqtl_study)
        node_label <- input$selected_protein_pqtl
        study_info <- input$selected_pqtl_study
        node_type <- "Protein"
      } else if (input$trait_view_type == "eQTL") {
        req(input$selected_gene_eqtl, input$selected_eqtl_study)
        node_label <- input$selected_gene_eqtl
        study_info <- input$selected_eqtl_study
        node_type <- "Gene"
      } else {  # mQTL
        req(input$selected_metabolite_mqtl, input$selected_mqtl_study)
        node_label <- input$selected_metabolite_mqtl
        study_info <- input$selected_mqtl_study
        node_type <- "Metabolite"
      }

      # Get custom colors from inputs (with defaults)
      # Note: #000000 is textInput's default, treat it as "not set"
      central_color <- if (!is.null(input$trait_central_color) && input$trait_central_color != "#000000") {
        input$trait_central_color
      } else {
        "#9b59b6"
      }
      region_color <- if (!is.null(input$trait_region_color) && input$trait_region_color != "#000000") {
        input$trait_region_color
      } else {
        "#3498db"
      }

      # Central node: The trait/protein/gene/metabolite (purple diamond)
      central_node <- data.frame(
        id = "trait",
        label = node_label,
        title = paste0(node_type, ": ", node_label, "<br>",
                      "Study: ", study_info, "<br>",
                      "Genomic regions: ", length(unique(paste(dt$CHR_var, dt$BP_START_var, dt$BP_STOP_var)))),
        shape = "diamond",
        size = 40,
        color = central_color,
        stringsAsFactors = FALSE
      )

      # Get unique regions
      regions <- unique(dt[, .(CHR_var, BP_START_var, BP_STOP_var, Prioritized_Gene, nearest_gene_1)])

      # Create region nodes
      region_nodes_list <- lapply(1:nrow(regions), function(i) {
        gene_label <- ifelse(!is.na(regions$Prioritized_Gene[i]) && regions$Prioritized_Gene[i] != "",
                            regions$Prioritized_Gene[i],
                            regions$nearest_gene_1[i])

        # Get PP.H4 for this region (take max if multiple)
        region_dt <- dt[CHR_var == regions$CHR_var[i] &
                       BP_START_var == regions$BP_START_var[i] &
                       BP_STOP_var == regions$BP_STOP_var[i]]
        max_pph4 <- max(region_dt$PP.H4.abf)

        data.frame(
          id = paste0("region_", i),
          label = gene_label,
          title = paste0("Gene: ", gene_label, "<br>",
                        "Region: chr", regions$CHR_var[i], ":",
                        regions$BP_START_var[i], "-", regions$BP_STOP_var[i], "<br>",
                        "Max PP.H4: ", round(max_pph4, 3)),
          shape = "dot",
          size = 10 + 30 * max_pph4,  # Size by PP.H4
          color = region_color,
          stringsAsFactors = FALSE
        )
      })

      region_nodes <- do.call(rbind, region_nodes_list)
      nodes <- rbind(central_node, region_nodes)

      # Apply custom labels if any
      labels <- trait_custom_labels()
      if (length(labels) > 0) {
        for (node_id in names(labels)) {
          if (node_id %in% nodes$id) {
            nodes$label[nodes$id == node_id] <- labels[[node_id]]
          }
        }
      }

      # Create edges with directionality
      edge_list <- lapply(1:nrow(regions), function(i) {
        # Get first matching row for this region
        region_dt <- dt[CHR_var == regions$CHR_var[i] &
                       BP_START_var == regions$BP_START_var[i] &
                       BP_STOP_var == regions$BP_STOP_var[i]]
        row <- region_dt[1]  # Take first row

        # Determine color and label based on directionality
        edge_color <- "#848484"  # Default gray
        dir_label <- "N/A"

        if ("directionality" %in% names(row) && !is.na(row$directionality)) {
          dir_val <- as.character(row$directionality)
          if (dir_val != "") {
            # Check if positive
            if (grepl("positive|concordant", dir_val, ignore.case = TRUE) ||
                (suppressWarnings(!is.na(as.numeric(dir_val))) && as.numeric(dir_val) > 0)) {
              edge_color <- "#e67e22"  # Orange
              dir_label <- "Positive"
            }
            # Check if negative
            else if (grepl("negative|discordant", dir_val, ignore.case = TRUE) ||
                     (suppressWarnings(!is.na(as.numeric(dir_val))) && as.numeric(dir_val) < 0)) {
              edge_color <- "#3498db"  # Blue
              dir_label <- "Negative"
            }
          }
        }

        data.frame(
          from = "trait",
          to = paste0("region_", i),
          width = 4,
          color = edge_color,
          title = paste0("PP.H4: ", round(row$PP.H4.abf, 3), "<br>",
                        "Directionality: ", dir_label),
          stringsAsFactors = FALSE
        )
      })

      edges <- do.call(rbind, edge_list)

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

    # Data for drill-down: filter by selected category, limit to 50 by signal strength
    conv_drilldown_data <- reactive({
      sel <- conv_selected_cat()
      if (is.null(sel)) return(NULL)
      dt <- filtered_region_data()
      studies_in_cat <- TRAIT_CATEGORIES[[sel]]$studies
      dt <- dt[source_study %in% studies_in_cat]
      if (nrow(dt) == 0) return(NULL)
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

      # Build metadata for each tile and the id->bundle_key map
      tiles_list <- dt[, {
        study <- source_study[1]
        trait_name <- tryCatch(get_trait_display_name(.SD),
                               error = function(e) basename(sumstats_2_file[1]))
        if (is.null(trait_name) || is.na(trait_name) || trait_name == "") {
          trait_name <- basename(sumstats_2_file[1])
        }
        study_color <- if (exists("study_colors") && study %in% names(study_colors)) {
          study_colors[[study]]
        } else cat_meta$color
        nlp <- if ("sumstats_2_max_nlog10P" %in% names(.SD)) sumstats_2_max_nlog10P[1] else NA
        list(
          tile_id = paste0("conv_t_", .I),
          bundle_key = paste0(study, "__", basename(sumstats_2_file[1])),
          trait_name = trait_name,
          study = study,
          study_color = study_color,
          nlp = nlp
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
          # Study badge + p-value on one line
          div(style = "display: flex; justify-content: space-between; align-items: center; gap: 4px;",
              tags$span(
                style = paste0(
                  "display: inline-block; padding: 1px 5px; ",
                  "border-radius: 8px; font-size: 9px; font-weight: 600; ",
                  "color: white; background: ", tl$study_color, "; ",
                  "white-space: nowrap; overflow: hidden; text-overflow: ellipsis; ",
                  "max-width: 70%;"
                ),
                tl$study
              ),
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
      # Top row (already sorted by -sumstats_2_max_nlog10P in conv_drilldown_data)
      study <- dt$source_study[1]
      bundle_key <- paste0(study, "__", basename(dt$sumstats_2_file[1]))
      conv_selected_trait(bundle_key)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    # Load selected trait sumstats
    conv_trait_sumstats <- reactive({
      req(regional_data_path(), regional_layout(), current_region(), conv_selected_trait())
      load_trait_sumstats(regional_data_path(), current_region(),
                          conv_selected_trait(), regional_layout())
    })

    # Display trait name for header
    conv_selected_trait_label <- reactive({
      sel <- conv_selected_trait()
      if (is.null(sel)) return(NULL)
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
      if (is.null(dt) || nrow(dt) == 0) {
        return(plotly_empty() %>% layout(
          title = list(text = "No regional plot data available for this trait",
                        font = list(size = 11))))
      }
      label <- conv_selected_trait_label()
      coords <- conv_region_coords()
      xr <- if (!is.null(coords)) c(coords$start, coords$end) else NULL
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
      # category (same order as the tiles: by -log10P desc)
      choices_list <- lapply(seq_len(nrow(dt)), function(i) {
        row <- dt[i]
        study <- row$source_study
        trait_name <- tryCatch(get_trait_display_name(row),
                               error = function(e) basename(row$sumstats_2_file))
        if (is.null(trait_name) || is.na(trait_name) || trait_name == "") {
          trait_name <- basename(row$sumstats_2_file)
        }
        nlp <- if ("sumstats_2_max_nlog10P" %in% names(row)) row$sumstats_2_max_nlog10P else NA
        label <- paste0(trait_name, " (", study, ", -log10P=",
                        if (is.na(nlp)) "-" else sprintf("%.1f", nlp), ")")
        key <- paste0(study, "__", basename(row$sumstats_2_file))
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
      req(filtered_data(), regional_data_path(), regional_layout(), current_region())
      dt <- filtered_data()
      if (nrow(dt) == 0) return(NULL)

      if (regional_layout() == "bundled") {
        # Bundled layout: match annotation rows to bundle keys
        region_data <- load_region_bundle(regional_data_path(), current_region())
        if (is.null(region_data) || length(region_data$traits) == 0) return(NULL)

        bundle_keys <- names(region_data$traits)

        # Build key from annotation: source_study__basename(sumstats_2_file)
        trait_info <- unique(dt[, .(
          source_study,
          trait_display = tryCatch(get_trait_display_name(.SD), error = function(e) "Unknown"),
          bundle_key = paste0(source_study, "__", basename(sumstats_2_file))
        ), by = seq_len(nrow(dt))])

        # Keep only traits that exist in the bundle
        trait_info <- trait_info[bundle_key %in% bundle_keys]
        if (nrow(trait_info) == 0) return(NULL)

        trait_info[, trait_label := paste0(trait_display, " (", source_study, ")")]
        setNames(trait_info$bundle_key, trait_info$trait_label)
      } else {
        # Legacy layout: use get_trait_name for filename matching
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

    # Load base study sumstats for current region
    regional_base_sumstats <- reactive({
      req(regional_data_path(), regional_layout(), current_region())
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

    # Load trait sumstats for selected colocalized trait
    regional_trait_sumstats <- reactive({
      req(regional_data_path(), regional_layout(), input$regional_trait_selector, current_region())
      if (input$regional_trait_selector == "") return(NULL)
      load_trait_sumstats(regional_data_path(), current_region(),
                          input$regional_trait_selector, regional_layout())
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
      # Use region bounds as x_range (not data extent)
      req(current_region())
      region_parts <- strsplit(current_region(), ":")[[1]]
      if (length(region_parts) >= 2) {
        pos_parts <- as.numeric(strsplit(region_parts[2], "-")[[1]])
        return(pos_parts)
      }

      # Fallback: use data extent
      base_sumstats <- regional_base_sumstats()
      trait_sumstats <- regional_trait_sumstats()

      all_pos <- c()
      if (!is.null(base_sumstats) && nrow(base_sumstats) > 0) {
        all_pos <- c(all_pos, base_sumstats$POS)
      }
      if (!is.null(trait_sumstats) && nrow(trait_sumstats) > 0) {
        all_pos <- c(all_pos, trait_sumstats$POS)
      }

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

    # === Trait View Outputs ===

    # === Trait View Node Selection for Export ===

    # Store available trait nodes
    trait_available_nodes <- reactiveVal(NULL)

    # Update trait node selector choices when network data changes
    observe({
      req(trait_network_data())
      net_data <- trait_network_data()

      if (nrow(net_data$nodes) > 0) {
        # Create node choices: id -> label (group)
        node_choices <- setNames(
          net_data$nodes$id,
          paste0(net_data$nodes$label, " (", net_data$nodes$group, ")")
        )
        # Store for quick selection buttons
        trait_available_nodes(net_data$nodes)

        # Update selectize with all nodes selected by default
        updateSelectizeInput(session, "trait_visible_nodes",
                            choices = node_choices,
                            selected = net_data$nodes$id)

        # Auto-enable physics simulation for large networks (>150 nodes)
        if (nrow(net_data$nodes) > 150 && !input$trait_physics) {
          updateCheckboxInput(session, "trait_physics", value = TRUE)
        }
      }
    })

    # Trait node selection: All
    observeEvent(input$trait_nodes_all, {
      req(trait_available_nodes())
      updateSelectizeInput(session, "trait_visible_nodes",
                          selected = trait_available_nodes()$id)
    })

    # Trait node selection: None
    observeEvent(input$trait_nodes_none, {
      updateSelectizeInput(session, "trait_visible_nodes", selected = character(0))
    })

    # Trait node selection: Top 10 by PP.H4
    observeEvent(input$trait_nodes_top10, {
      req(trait_filtered_data(), trait_available_nodes())
      dt <- trait_filtered_data()

      # Get top 10 by PP.H4
      top_10 <- head(dt[order(-PP.H4.abf)], 10)

      # Build matching node IDs (region nodes have format "region_X")
      top_node_ids <- paste0("region_", seq_len(nrow(top_10)))

      # Always include the central trait node
      selected_ids <- c("trait", top_node_ids)

      updateSelectizeInput(session, "trait_visible_nodes", selected = selected_ids)
    })

    # === Trait View Customize for Export ===

    # Reactive values for trait view custom labels
    trait_custom_labels <- reactiveVal(list())  # node_id -> custom_label

    # When a trait node is clicked, populate the label input with current label
    observe({
      req(input$trait_selected_node)
      net_data <- trait_network_data()

      node <- net_data$nodes[net_data$nodes$id == input$trait_selected_node, ]
      if (nrow(node) > 0) {
        # Check if there's a custom label, otherwise use the original
        labels <- trait_custom_labels()
        current_label <- if (input$trait_selected_node %in% names(labels)) {
          labels[[input$trait_selected_node]]
        } else {
          node$label
        }
        updateTextInput(session, "trait_custom_node_label", value = current_label)
      }
    })

    # Apply custom label when button is clicked
    observeEvent(input$trait_apply_label, {
      req(input$trait_selected_node, input$trait_custom_node_label)

      labels <- trait_custom_labels()
      labels[[input$trait_selected_node]] <- input$trait_custom_node_label
      trait_custom_labels(labels)

      # Show confirmation
      showNotification(
        paste0("Label updated for: ", input$trait_selected_node),
        type = "message",
        duration = 2
      )
    })

    # Trait PNG Export via JavaScript
    observeEvent(input$trait_export_png, {
      # Get current trait name for filename
      trait_name <- switch(input$trait_view_type,
        "Phenotypes" = input$selected_trait_pheno,
        "pQTL" = input$selected_protein_pqtl,
        "eQTL" = input$selected_gene_eqtl,
        "mQTL" = input$selected_metabolite_mqtl,
        "trait"
      )
      filename <- paste0(gsub("[^a-zA-Z0-9]", "_", trait_name), "_network.png")

      # Trigger JavaScript export (reuse same handler, different network ID)
      session$sendCustomMessage("exportTraitNetwork", list(filename = filename))
    })

    # Trait HTML Export (interactive)
    output$trait_export_html <- downloadHandler(
      filename = function() {
        trait_name <- switch(input$trait_view_type,
          "Phenotypes" = input$selected_trait_pheno,
          "pQTL" = input$selected_protein_pqtl,
          "eQTL" = input$selected_gene_eqtl,
          "mQTL" = input$selected_metabolite_mqtl,
          "trait"
        )
        if (is.null(trait_name) || trait_name == "") trait_name <- "trait"
        paste0(gsub("[^a-zA-Z0-9]", "_", trait_name), "_network.html")
      },
      content = function(file) {
        net <- trait_network_data()

        # Check if we have data
        if (is.null(net) || nrow(net$nodes) == 0) {
          # Create empty network
          network <- visNetwork(data.frame(id = 1, label = "No data"), data.frame())
          visSave(network, file)
          return()
        }

        # Filter nodes based on selection
        visible <- input$trait_visible_nodes
        if (!is.null(visible) && length(visible) > 0) {
          nodes_filtered <- net$nodes[net$nodes$id %in% visible, ]
          edges_filtered <- net$edges[
            net$edges$from %in% visible & net$edges$to %in% visible,
          ]
        } else {
          nodes_filtered <- net$nodes
          edges_filtered <- net$edges
        }

        # Create visNetwork and save
        network <- visNetwork(nodes_filtered, edges_filtered) %>%
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

    # Render trait network
    output$trait_network <- renderVisNetwork({
      net <- trait_network_data()

      if (nrow(net$nodes) == 0) {
        visNetwork(data.frame(id = 1, label = "No data"), data.frame()) %>%
          visNodes(shape = "text", font = list(size = 20))
      } else {
        # Filter nodes based on selection (if any selected)
        visible <- input$trait_visible_nodes
        if (!is.null(visible) && length(visible) > 0) {
          # Filter nodes to only show selected ones
          nodes_filtered <- net$nodes[net$nodes$id %in% visible, ]
          # Filter edges to only include those connecting visible nodes
          edges_filtered <- net$edges[
            net$edges$from %in% visible & net$edges$to %in% visible,
          ]
        } else {
          # No selection means show all (during initial load)
          nodes_filtered <- net$nodes
          edges_filtered <- net$edges
        }

        # Show notification for large networks (dismissed when stabilized)
        n_nodes <- nrow(nodes_filtered)
        if (n_nodes > 100) {
          showNotification(
            paste0("Rendering network (", n_nodes, " nodes)... please wait"),
            id = "trait_network_loading",
            duration = NULL,
            type = "message"
          )
        }

        visNetwork(nodes_filtered, edges_filtered) %>%
          visOptions(
            highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)
          ) %>%
          visPhysics(enabled = input$trait_physics,
                     stabilization = list(enabled = TRUE, iterations = 200)) %>%
          visLayout(randomSeed = 123) %>%
          visInteraction(
            dragNodes = TRUE,
            dragView = TRUE,
            zoomView = TRUE
          ) %>%
          visEvents(
            select = "function(nodes) {
              Shiny.setInputValue('trait_selected_node', nodes.nodes[0]);
            }",
            stabilizationIterationsDone = "function() {
              Shiny.setInputValue('trait_network_stabilized', Math.random());
            }"
          )
      }
    })

    # Dismiss loading notification when trait network stabilizes
    observeEvent(input$trait_network_stabilized, {
      removeNotification(id = "trait_network_loading")
    })

    # Trait/Protein/Gene/Metabolite information summary (unified)
    output$trait_info <- renderPrint({
      req(trait_filtered_data(), input$trait_view_type)
      dt <- trait_filtered_data()

      if (input$trait_view_type == "Phenotypes") {
        req(input$selected_trait_pheno, input$selected_pheno_study)
        cat("Study:", input$selected_pheno_study, "\n")
        cat("Trait:", input$selected_trait_pheno, "\n\n")
      } else if (input$trait_view_type == "pQTL") {
        req(input$selected_protein_pqtl, input$selected_pqtl_study)
        cat("Study:", input$selected_pqtl_study, "\n")
        cat("Protein:", input$selected_protein_pqtl, "\n\n")
      } else if (input$trait_view_type == "eQTL") {
        req(input$selected_gene_eqtl, input$selected_eqtl_study)
        cat("Study:", input$selected_eqtl_study, "\n")
        cat("Gene:", input$selected_gene_eqtl, "\n\n")
      } else {  # mQTL
        req(input$selected_metabolite_mqtl, input$selected_mqtl_study)
        cat("Study:", input$selected_mqtl_study, "\n")
        cat("Metabolite:", input$selected_metabolite_mqtl, "\n\n")
      }

      cat("Genomic regions:", length(unique(paste(dt$CHR_var, dt$BP_START_var, dt$BP_STOP_var))), "\n")
      cat("Total colocalizations:", nrow(dt), "\n\n")
      cat("Mean PP.H4:", round(mean(dt$PP.H4.abf), 3), "\n")
      cat("Median PP.H4:", round(median(dt$PP.H4.abf), 3), "\n")
      cat("Max PP.H4:", round(max(dt$PP.H4.abf), 3), "\n")
    })

    # === Trait View Regional Plot Logic ===

    # Get available regions for current trait (from trait_filtered_data)
    trait_regional_regions_available <- reactive({
      req(trait_filtered_data())
      dt <- trait_filtered_data()

      if (nrow(dt) == 0) return(NULL)

      # Get unique regions with their gene labels
      regions <- unique(dt[, .(
        CHR_var, BP_START_var, BP_STOP_var,
        gene = ifelse(!is.na(Prioritized_Gene) & Prioritized_Gene != "",
                     Prioritized_Gene, nearest_gene_1),
        max_pph4 = max(PP.H4.abf)
      ), by = .(CHR_var, BP_START_var, BP_STOP_var)])

      # Build choices: region_id -> display label
      regions[, region_id := paste0("chr", CHR_var, ":", BP_START_var, "-", BP_STOP_var)]
      regions[, region_label := paste0(gene, " (", region_id, ") PP.H4=", round(max_pph4, 2))]

      # Sort by PP.H4
      regions <- regions[order(-max_pph4)]

      setNames(regions$region_id, regions$region_label)
    })

    # Update region selector when trait changes
    observe({
      regions <- trait_regional_regions_available()

      if (is.null(regions) || length(regions) == 0) {
        updateSelectInput(session, "trait_regional_region_selector",
                         choices = c("No colocalized regions" = ""),
                         selected = "")
      } else {
        updateSelectInput(session, "trait_regional_region_selector",
                         choices = regions,
                         selected = regions[1])
      }
    })

    # Get current trait name for file lookup
    trait_view_current_trait_name <- reactive({
      req(input$trait_view_type)

      trait_name <- switch(input$trait_view_type,
        "Phenotypes" = input$selected_trait_pheno,
        "pQTL" = input$selected_protein_pqtl,
        "eQTL" = input$selected_gene_eqtl,
        "mQTL" = input$selected_metabolite_mqtl,
        NULL
      )

      trait_name
    })

    # Get current study type for file lookup
    trait_view_current_study_type <- reactive({
      req(input$trait_view_type)

      study_type <- switch(input$trait_view_type,
        "Phenotypes" = input$selected_pheno_study,
        "pQTL" = input$selected_pqtl_study,
        "eQTL" = "eQTLGen",
        "mQTL" = input$selected_mqtl_study,
        NULL
      )

      study_type
    })

    # Load base study sumstats for selected region in Trait View
    trait_regional_base_sumstats <- reactive({
      req(regional_data_path(), regional_layout(), input$trait_regional_region_selector)
      if (input$trait_regional_region_selector == "") return(NULL)
      load_base_sumstats(regional_data_path(), input$trait_regional_region_selector, regional_layout())
    })

    # Helper function to get trait file name from trait_filtered_data
    # This maps from the UI trait name (e.g. MVP_R4_Analyzed.variable) to the file name column
    get_trait_file_name <- function(dt, study, ui_trait_name, region_chr, region_start, region_end) {
      # Map study to UI column and file name column
      col_map <- list(
        "MVP_R4" = list(ui = "MVP_R4_Analyzed.variable", file = "MVP_R4_Title.of.analysis"),
        "FinnGen_r9" = list(ui = "FinnGen_r9_phenotype", file = "FinnGen_r9_phenotype"),
        "UKB_TOPMed" = list(ui = "UKB_TOPMed_phenostring", file = "UKB_TOPMed_phenostring"),
        "UKB_PPP_EUR" = list(ui = "UKB_PPP_EUR_olink_target_fullname", file = "UKB_PPP_EUR_olink_target_fullname"),
        "Icelanders_pGWAS" = list(ui = "Icelanders_pGWAS_Protein..short.name.", file = "Icelanders_pGWAS_Protein..short.name."),
        "eQTLGen" = list(ui = "eQTLGen_gene_name", file = NULL),  # eQTLGen uses special ENSG format
        "GCKD_mGWAS_plasma" = list(ui = "GCKD_mGWAS_plasma_BIOCHEMICAL", file = "GCKD_mGWAS_plasma_BIOCHEMICAL"),
        "GCKD_mGWAS_urine" = list(ui = "GCKD_mGWAS_urine_BIOCHEMICAL", file = "GCKD_mGWAS_urine_BIOCHEMICAL")
      )

      if (!study %in% names(col_map)) return(ui_trait_name)

      ui_col <- col_map[[study]]$ui
      file_col <- col_map[[study]]$file

      # Filter to matching region and trait
      matching <- dt[CHR_var == region_chr &
                     BP_START_var == region_start &
                     BP_STOP_var == region_end &
                     source_study == study]

      if (ui_col %in% names(matching)) {
        matching <- matching[get(ui_col) == ui_trait_name]
      }

      if (nrow(matching) == 0) return(ui_trait_name)

      # Get file name from file column
      if (!is.null(file_col) && file_col %in% names(matching)) {
        file_name <- matching[[file_col]][1]
        if (!is.na(file_name) && file_name != "") {
          # Clean for filename safety (same as get_trait_name_for_regional)
          file_name <- gsub("[^A-Za-z0-9_.-]", "_", file_name)
          file_name <- gsub("_+", "_", file_name)
          file_name <- gsub("^_|_$", "", file_name)
          return(file_name)
        }
      }

      # For eQTLGen, construct from gene name and ENSG
      if (study == "eQTLGen" && "eQTLGen_gene_name" %in% names(matching) &&
          "sumstats_2_file" %in% names(matching)) {
        gene_name <- matching$eQTLGen_gene_name[1]
        ensg <- gsub(".*_(ENSG[0-9]+).*", "\\1", basename(matching$sumstats_2_file[1]))
        if (!is.na(gene_name) && gene_name != "" && grepl("^ENSG", ensg)) {
          return(paste0(gene_name, "_", ensg))
        }
      }

      return(ui_trait_name)
    }

    # Load trait sumstats for current trait in Trait View
    trait_regional_trait_sumstats <- reactive({
      req(regional_data_path(), trait_view_current_study_type(), trait_view_current_trait_name(),
          input$trait_regional_region_selector, trait_filtered_data())

      if (input$trait_regional_region_selector == "") return(NULL)

      # Parse region for folder path
      region_parts <- strsplit(input$trait_regional_region_selector, ":")[[1]]
      chr <- region_parts[1]
      # Add "chr" prefix if not present (for folder path)
      chr_folder <- chr
      if (!grepl("^chr", chr_folder)) chr_folder <- paste0("chr", chr_folder)
      pos_parts <- strsplit(region_parts[2], "-")[[1]]
      start_pos <- pos_parts[1]
      end_pos <- pos_parts[2]
      region_folder <- paste0(chr_folder, "_", start_pos, "_", end_pos)

      # Get file name from data (maps UI name to file name)
      study <- trait_view_current_study_type()
      ui_trait <- trait_view_current_trait_name()
      dt <- trait_filtered_data()

      # Strip "chr" prefix for data lookup (data uses "12", not "chr12")
      chr_for_lookup <- gsub("^chr", "", chr)

      file_trait_name <- get_trait_file_name(dt, study, ui_trait, chr_for_lookup,
                                             as.numeric(start_pos), as.numeric(end_pos))

      trait_key <- paste0(study, "__", file_trait_name)
      load_trait_sumstats(regional_data_path(), input$trait_regional_region_selector,
                          trait_key, regional_layout())
    })

    # Get lead SNP for selected region in Trait View
    trait_regional_lead_snp <- reactive({
      req(trait_filtered_data(), input$trait_regional_region_selector)

      if (input$trait_regional_region_selector == "") return(NULL)

      dt <- trait_filtered_data()

      # Parse region
      region_parts <- strsplit(input$trait_regional_region_selector, ":")[[1]]
      chr <- region_parts[1]
      pos_parts <- strsplit(region_parts[2], "-")[[1]]
      start_pos <- as.numeric(pos_parts[1])
      end_pos <- as.numeric(pos_parts[2])

      # Filter to selected region
      region_dt <- dt[CHR_var == chr & BP_START_var == start_pos & BP_STOP_var == end_pos]

      if (nrow(region_dt) > 0 && "sumstats_1_ind_Name" %in% names(region_dt)) {
        region_dt$sumstats_1_ind_Name[1]
      } else {
        NULL
      }
    })

    # Render base study plot title (Trait View)
    output$trait_regional_plot_title_base <- renderText({
      req(current_study())
      paste0("Base Study: ", current_study())
    })

    # Render trait plot title (Trait View)
    output$trait_regional_plot_title_trait <- renderText({
      req(trait_view_current_study_type(), trait_view_current_trait_name())
      paste0(trait_view_current_study_type(), ": ", trait_view_current_trait_name())
    })

    # Render coloc info (Trait View)
    output$trait_regional_coloc_info <- renderText({
      req(trait_filtered_data(), input$trait_regional_region_selector)

      if (input$trait_regional_region_selector == "") return("")

      dt <- trait_filtered_data()

      # Parse region
      region_parts <- strsplit(input$trait_regional_region_selector, ":")[[1]]
      chr <- region_parts[1]
      pos_parts <- strsplit(region_parts[2], "-")[[1]]
      start_pos <- as.numeric(pos_parts[1])
      end_pos <- as.numeric(pos_parts[2])

      # Filter to selected region
      region_dt <- dt[CHR_var == chr & BP_START_var == start_pos & BP_STOP_var == end_pos]

      if (nrow(region_dt) > 0) {
        max_pph4 <- max(region_dt$PP.H4.abf, na.rm = TRUE)
        paste0("PP.H4: ", round(max_pph4, 3))
      } else {
        ""
      }
    })

    # Calculate shared X-axis range for Trait View plots
    trait_regional_x_range <- reactive({
      base_sumstats <- trait_regional_base_sumstats()
      trait_sumstats <- trait_regional_trait_sumstats()

      # Collect all positions from both datasets
      all_pos <- c()
      if (!is.null(base_sumstats) && nrow(base_sumstats) > 0) {
        all_pos <- c(all_pos, base_sumstats$POS)
      }
      if (!is.null(trait_sumstats) && nrow(trait_sumstats) > 0) {
        all_pos <- c(all_pos, trait_sumstats$POS)
      }

      if (length(all_pos) > 0) {
        # Add small padding (1% on each side)
        range_span <- max(all_pos) - min(all_pos)
        padding <- range_span * 0.01
        return(c(min(all_pos) - padding, max(all_pos) + padding))
      }
      NULL
    })

    # Render base study regional plot (Trait View - interactive)
    output$trait_regional_plot_base <- renderPlotly({
      sumstats <- trait_regional_base_sumstats()
      lead_snp <- if (input$trait_regional_highlight_lead) trait_regional_lead_snp() else NULL
      x_range <- trait_regional_x_range()

      plot_regional_association_interactive(
        sumstats_dt = sumstats,
        title = "",
        highlight_snp = lead_snp,
        color = "#2ecc71",  # Green for base study
        x_range = x_range
      )
    })

    # Render trait regional plot (Trait View - interactive)
    output$trait_regional_plot_trait <- renderPlotly({
      sumstats <- trait_regional_trait_sumstats()
      lead_snp <- if (input$trait_regional_highlight_lead) trait_regional_lead_snp() else NULL
      x_range <- trait_regional_x_range()

      # Get color based on study type
      study_type <- trait_view_current_study_type()
      study_color <- "#3498db"  # Default blue
      if (!is.null(study_type) && study_type %in% names(study_colors)) {
        study_color <- study_colors[[study_type]]
      }

      plot_regional_association_interactive(
        sumstats_dt = sumstats,
        title = "",
        highlight_snp = lead_snp,
        color = study_color,
        x_range = x_range
      )
    })

    # Render gene track (Trait View)
    output$trait_regional_gene_track <- renderPlotly({
      req(input$trait_regional_region_selector)

      # Parse region coordinates
      region_sel <- input$trait_regional_region_selector
      region_parts <- strsplit(region_sel, ":")[[1]]
      chr <- region_parts[1]
      pos_parts <- strsplit(region_parts[2], "-")[[1]]
      start_pos <- as.numeric(pos_parts[1])
      end_pos <- as.numeric(pos_parts[2])

      plot_gene_track(chr, start_pos, end_pos, source = "trait_regional_gene_track")
    })

    # Store selected gene for Trait View info panel
    trait_selected_gene <- reactiveVal(NULL)

    # Handle click on gene track (Trait View)
    observeEvent(event_data("plotly_click", source = "trait_regional_gene_track"), {
      click_data <- event_data("plotly_click", source = "trait_regional_gene_track")
      if (!is.null(click_data) && !is.null(click_data$customdata)) {
        clicked_gene <- click_data$customdata
        if (!is.null(gene_annotation)) {
          gene_info <- gene_annotation[gene_name == clicked_gene]
          if (nrow(gene_info) > 0) {
            trait_selected_gene(gene_info[1])
          }
        }
      }
    })

    # Render gene info panel (Trait View)
    output$trait_gene_info_panel <- renderUI({
      gene <- trait_selected_gene()
      if (is.null(gene)) return(NULL)

      # Helper for source tooltip (ⓘ with hover description, not clickable)
      source_tip <- function(source_name, desc, source_url) {
        tags$span(
          title = paste0(source_name, " - ", desc),
          style = "color: #999; margin-left: 4px; cursor: help;",
          "ⓘ"
        )
      }

      info_items <- list()

      info_items[[length(info_items) + 1]] <- tags$h4(
        gene$gene_name,
        if (!is.null(gene$full_name) && gene$full_name != "")
          tags$span(style = "color: #666; margin-left: 10px; font-weight: normal;", gene$full_name)
      )

      info_items[[length(info_items) + 1]] <- tags$p(
        tags$strong("Location: "),
        paste0("chr", gene$chr, ":", format(gene$start, big.mark = ","), "-", format(gene$end, big.mark = ","),
               " (", gene$strand, " strand)"),
        source_tip("Gencode v49", "comprehensive gene annotation", "https://www.gencodegenes.org/")
      )

      if (!is.null(gene$gene_family) && gene$gene_family != "") {
        info_items[[length(info_items) + 1]] <- tags$p(
          tags$strong("Gene Family: "), gene$gene_family,
          source_tip("HGNC", "HUGO Gene Nomenclature Committee", "https://www.genenames.org/")
        )
      }

      if (!is.null(gene$ncbi_description) && gene$ncbi_description != "") {
        info_items[[length(info_items) + 1]] <- tags$p(
          tags$strong("Description: "), gene$ncbi_description,
          source_tip("NCBI Gene", "gene-specific information", "https://www.ncbi.nlm.nih.gov/gene/")
        )
      }

      if (!is.null(gene$pathways) && gene$pathways != "") {
        info_items[[length(info_items) + 1]] <- tags$p(
          tags$strong("Pathways (top 5): "), gene$pathways,
          source_tip("Reactome", "curated biological pathways", "https://reactome.org/")
        )
      }

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

      wellPanel(
        style = "background-color: #f8f9fa; border: 1px solid #dee2e6; margin-top: 10px;",
        do.call(tagList, info_items)
      )
    })

    # Workflow diagram for Documentation tab
    output$workflow_diagram <- renderPlot({
      par(mar = c(1, 1, 3, 1), bg = "white")
      plot(NULL, xlim = c(0, 10), ylim = c(0, 11), 
           xlab = "", ylab = "", axes = FALSE)
      
      title("From GWAS Discovery to Biological Insight", 
            cex.main = 1.5, font.main = 2)
      
      # Box colors
      col_gwas <- "#3498db"
      col_challenge <- "#e74c3c"
      col_resources <- "#2ecc71"
      col_solution <- "#f39c12"
      
      # Level 1: GWAS Studies
      rect(1, 9.2, 3, 10.3, col = col_gwas, border = "black", lwd = 2)
      text(2, 9.9, "Multiple GWAS", cex = 0.9, font = 2, col = "white")
      text(2, 9.5, "Meta-analyses", cex = 0.8, col = "white")
      
      # Arrow down
      arrows(2, 9.2, 2, 8.3, lwd = 2, length = 0.15)
      
      # Level 2: Loci discovered
      rect(0.8, 7.3, 3.2, 8.2, col = col_challenge, border = "black", lwd = 2)
      text(2, 7.9, "Hundreds-Thousands", cex = 0.85, font = 2, col = "white")
      text(2, 7.5, "of Significant Loci", cex = 0.85, font = 2, col = "white")
      
      # Challenge boxes
      rect(4, 8.4, 5.8, 9.2, col = "#ecf0f1", border = col_challenge, lwd = 2)
      text(4.9, 8.8, "Many are intergenic", cex = 0.7)
      
      rect(4, 7.3, 5.8, 8.1, col = "#ecf0f1", border = col_challenge, lwd = 2)
      text(4.9, 7.7, "Gene function", cex = 0.7)
      text(4.9, 7.5, "often unknown", cex = 0.7)
      
      # Arrow to question
      arrows(2, 7.3, 2, 6.4, lwd = 2, length = 0.15)
      
      # Question box
      rect(0.8, 5.4, 3.2, 6.3, col = col_challenge, border = "black", lwd = 2)
      text(2, 6, "What biological/", cex = 0.85, font = 2, col = "white")
      text(2, 5.6, "clinical pathways?", cex = 0.85, font = 2, col = "white")
      
      # Resources column (right side)
      text(7.5, 10.2, "Integration Resources", cex = 1.1, font = 2)
      
      resource_y <- seq(9.5, 6, length.out = 7)
      resources <- c("PheWAS", "Phenotypes", "pGWAS", "Proteomics", 
                     "eQTL (blood)", "Tissue eQTL", "Other GWAS")
      
      for (i in seq_along(resources)) {
        rect(6.5, resource_y[i] - 0.25, 8.5, resource_y[i] + 0.25, 
             col = col_resources, border = "black", lwd = 1)
        text(7.5, resource_y[i], resources[i], cex = 0.7, col = "white", font = 2)
      }
      
      # Arrows connecting
      arrows(3.2, 5.85, 4.2, 4.8, lwd = 2.5, length = 0.15, col = col_solution)
      arrows(6.5, 7.75, 5.8, 4.8, lwd = 2.5, length = 0.15, col = col_solution)
      
      # Solution box
      rect(3.5, 3.8, 6.5, 5.3, col = col_solution, border = "black", lwd = 3)
      text(5, 4.95, "Colocalization", cex = 1.1, font = 2)
      text(5, 4.55, "Analysis", cex = 1.1, font = 2)
      text(5, 4.1, "Network Viewer", cex = 0.9, font = 1)
      
      # Bottom: Benefits
      rect(0.5, 0.5, 9.5, 2.8, col = "#ecf0f1", border = col_solution, lwd = 2)
      text(5, 2.4, "Benefits of This Tool:", cex = 1, font = 2)
      text(5, 1.9, "• Identify multi-trait associations across all resources simultaneously", 
           cex = 0.75, adj = 0.5)
      text(5, 1.4, "• Discover biological pathways and cross-phenotype relationships", 
           cex = 0.75, adj = 0.5)
      text(5, 0.9, "• Generate hypotheses about gene function and clinical connections", 
           cex = 0.75, adj = 0.5)
      
    }, height = 450)
    
  }

# Run the application
shinyApp(ui = ui, server = server)
