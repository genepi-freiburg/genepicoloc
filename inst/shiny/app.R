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
library(data.table)
library(plotly)
library(DT)

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
               # Top row: Search + Chromosome (left) + mini Manhattan (right).
               # Bottom border separates the study/region selection area from
               # the lower region-drilldown area so new users have a clear
               # "top = pick a region, bottom = explore it" mental model.
               div(style = "border-bottom: 1px solid #dce1e7; padding-bottom: 8px; margin-bottom: 8px;",
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
                         width = "100%"),
                       # Chromosome zoom: restrict the mini Manhattan to one
                       # chromosome. Choices are populated dynamically from
                       # coloc_data() so MVP chr16-only test bundles only
                       # show chr16.
                       selectInput("manhattan_chr",
                         "Chromosome:",
                         choices = c("All" = "all"),
                         selected = "all",
                         width = "100%")
                     )
                   ),
                   column(10, uiOutput("mini_manhattan_ui"))
                 )
               ),

               fluidRow(
                 # Left side: compact study header, category strip, trait
                 # list (DT). Width 5 leaves width 7 for the RAPs. The
                 # old left sidebar and all its collapsible filter blocks
                 # (PP.H4, nlog10P, study selector) are gone - column
                 # filters in the DT cover the same use cases.
                 column(
                   width = 5,
                   div(style = "padding: 4px 6px;",
                     # Hidden dummy input to preserve reactivity for legacy coord selector
                     tags$div(style = "display: none;",
                       selectizeInput("selected_region_coord", NULL,
                                      choices = NULL, selected = NULL)),

                     # Horizontal category strip (6 compact buttons)
                     uiOutput("conv_category_cards"),

                     # Drilldown header + CSV download on the same row
                     div(style = "display: flex; align-items: center; gap: 8px; margin: 2px 0 4px 0;",
                       div(style = "flex: 1; min-width: 0;",
                         uiOutput("conv_drilldown_header")),
                       downloadButton("download_data", "CSV",
                                      class = "btn-xs btn-default")
                     ),
                     DT::dataTableOutput("conv_trait_list", height = "100%")
                   )
                 ),

                 # Right side: region-centric multi-omics view.
                 # Left border acts as a vertical separator between the
                 # left (category + trait list) and right (RAP + genes)
                 # halves of the lower area.
                 column(
                   width = 7,
                   style = "border-left: 1px solid #dce1e7;",
                   div(style = "padding-left: 10px;",
                     # Trait of interest RAP (title is the Manhattan;
                     # no per-panel header to save vertical space)
                     plotlyOutput("conv_base_plot", height = "200px"),
                     # Selected trait RAP - shown only when a trait is
                     # picked in the middle-panel list. Title suppressed
                     # so the plot itself is the indicator.
                     conditionalPanel(
                       condition = "output.conv_trait_has_data",
                       plotlyOutput("conv_trait_plot", height = "200px")
                     ),
                     # Gene track + gene info panel (click a gene)
                     plotlyOutput("conv_gene_track", height = "100px"),
                     uiOutput("gene_info_panel")
                   )
                 )  # End main column
               )  # End fluidRow
      ),
      tabPanel("Documentation",
               div(
                 style = "max-width: 900px; margin: 20px auto; padding: 0 16px;",
                 includeMarkdown("www/docs.md")
               )
      )
    )
  )
  
  server <- function(input, output, session) {

    # -----------------------------------------------------------------
    # Server structure (top -> bottom)
    #
    #   Landing page render                   output$landing_categories
    #   Reactive values (all reactiveVals)    current_study, selected_gene, ...
    #   Event handlers
    #     - study selection card click        observeEvent(selected_study)
    #     - gene / region selectize sync      observe({ ... })
    #     - Manhattan chr selector
    #     - Manhattan click
    #   Core data reactives                   coloc_data, regions_data, filtered_*
    #   Mini Manhattan output                 output$mini_manhattan*
    #   Convergence reactives                 conv_*
    #   Regional bundle reactives             regional_data_path, regional_base_sumstats
    #   Trait name helpers                    get_trait_name, get_trait_display_name
    #   Gene info panel                       output$gene_info_panel
    #   Download handler                      output$download_data
    #
    # Pure helpers (categories, bundle keys, region parsing, bundles,
    # Manhattan builder, plots, gene track, annotations) live in
    # inst/shiny/R/*.R and are auto-sourced before this file runs.
    # -----------------------------------------------------------------

    # === Landing Page ===

    output$landing_categories <- renderUI({
      featured   <- Filter(function(c) isTRUE(c$section == "featured"),   atlas_categories)
      additional <- Filter(function(c) isTRUE(c$section == "additional"), atlas_categories)
      # Categories without an explicit section fall into "featured" so
      # adding a new card without touching the helper still works.
      extras     <- Filter(function(c) is.null(c$section), atlas_categories)
      featured   <- c(featured, extras)

      tagList(
        div(class = "landing-section-header",
            h2("Featured atlases", style = "margin: 4px 0 10px 0;")),
        div(class = "category-grid",
          lapply(featured, build_card)),
        if (length(additional) > 0) tagList(
          div(class = "landing-section-header",
              h2("Additional phenotypes",
                 style = "margin: 28px 0 10px 0; color: #555; font-size: 18px;")),
          div(class = "category-grid",
            lapply(additional, build_card))
        )
      )
    })

    # === Reactive values (all declared up front) ===
    #
    # current_study         - user's currently loaded study (card click)
    # selected_gene         - clicked gene on the gene track
    # region_to_value_map   - bare region_key -> gene-prefixed selectize id
    # conv_selected_cat     - convergence drilldown: selected category
    # conv_selected_trait   - convergence drilldown: selected trait key
    current_study       <- reactiveVal(NULL)
    selected_gene       <- reactiveVal(NULL)
    region_to_value_map <- reactiveVal(character(0))
    conv_selected_cat   <- reactiveVal(NULL)
    conv_selected_trait <- reactiveVal(NULL)

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

        # Reset the chromosome zoom filter so the new study's Manhattan
        # shows "All" by default. (The observer that populates the
        # choices fires on coloc_data() change and will restore "All"
        # when the new list is built, but we also clear it eagerly so
        # a stale selection doesn't flash before that happens.)
        updateSelectInput(session, "manhattan_chr", selected = "all")

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

    # Update region selectors.
    # region_to_value_map (declared at the top) holds the mapping from
    # bare region keys to gene-prefixed selectize values.
    observe({
      req(regions_data())
      regions <- regions_data()

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

    # Populate the chromosome selector from whichever chromosomes the
    # current study actually has colocs on. Numeric sort, X last.
    observe({
      dt <- coloc_data()
      if (is.null(dt) || nrow(dt) == 0) {
        updateSelectInput(session, "manhattan_chr",
                          choices = c("All" = "all"), selected = "all")
        return()
      }
      chrs <- unique(as.character(dt$CHR_var))
      chr_num <- suppressWarnings(as.integer(chrs))
      chrs <- chrs[order(is.na(chr_num), chr_num, chrs)]
      choices <- c("All" = "all", setNames(chrs, paste0("chr", chrs)))
      current <- isolate(input$manhattan_chr)
      sel <- if (!is.null(current) && current %in% choices) current else "all"
      updateSelectInput(session, "manhattan_chr",
                        choices = choices, selected = sel)
    })

    output$mini_manhattan_ui <- renderUI({
      req(regions_data(), coloc_data())
      plotlyOutput("mini_manhattan", height = "200px")
    })

    # Mini Manhattan: the heavy lifting (consensus clustering, gene
    # labels, plotly assembly) lives in R/manhattan.R as two pure
    # helpers. This server reactive only wires inputs to those helpers,
    # handles the chromosome-zoom filter, and computes the title.
    output$mini_manhattan <- renderPlotly({
      req(regions_data(), coloc_data())
      regions <- regions_data()
      dt <- coloc_data()

      # Chromosome filter (zoom): when the user picks a specific chromosome
      # from the dropdown, restrict both `dt` and `regions` to that chr.
      chr_sel <- input$manhattan_chr
      if (!is.null(chr_sel) && chr_sel != "all") {
        dt      <- dt[as.character(CHR_var) == chr_sel]
        regions <- regions[as.character(CHR_var) == chr_sel]
        if (nrow(dt) == 0) {
          return(plotly_empty() %>% layout(
            title = list(text = paste0("No colocs on chr", chr_sel),
                         x = 0.5, y = 0.5)))
        }
      }

      is_multi <- !is.null(current_virtual_info()) && "ancestry" %in% names(dt)
      region_stats <- manhattan_consensus_table(
        dt, regions, is_multi,
        selected_region_key = current_region()
      )

      # Build the Manhattan title: "<study> - N regions".
      study_label <- if (!is.null(current_study())) {
        tryCatch(trait_label(current_study()),
                 error = function(e) current_study())
      } else {
        "Atlas"
      }
      plot_title <- paste0(
        study_label, " - ",
        format(nrow(region_stats), big.mark = ","),
        if (nrow(region_stats) == 1) " region" else " regions"
      )

      manhattan_plot(region_stats, plot_title = plot_title)
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
    
    # Load and process data - unified path for all studies.
    # Every study gets an `ancestry` column (single-ancestry = "EUR").
    coloc_data <- reactive({
      if (is.null(current_study())) return(NULL)
      study <- current_study()
      entry <- DEFAULT_STUDY_REGISTRY[[study]]
      if (is.null(entry)) return(NULL)

      # Load all ancestry files and tag with ancestry column
      parts <- lapply(entry$ancestries, function(anc) {
        f <- entry$coloc_files[[anc]]
        if (is.null(f) || !file.exists(f)) return(NULL)
        dt <- readRDS(f)
        if (!data.table::is.data.table(dt)) dt <- data.table::as.data.table(dt)
        if (nrow(dt) == 0) return(NULL)
        dt[, ancestry := anc]
        # Strip per-ancestry prefixes from metadata columns (e.g.
        # MVP_R4_AFR_* -> MVP_R4_*) so downstream sees a unified schema.
        anc_prefix <- paste0("_", anc, "_")
        renamed <- grep(anc_prefix, names(dt), value = TRUE, fixed = TRUE)
        if (length(renamed) > 0) {
          setnames(dt, renamed, sub(anc_prefix, "_", renamed, fixed = TRUE))
        }
        # Unify source_study (strip ancestry suffix if present)
        if ("source_study" %in% names(dt)) {
          dt[, source_study := sub(paste0("_", anc, "$"), "", source_study)]
        }
        dt
      })
      parts <- parts[!vapply(parts, is.null, logical(1))]
      if (length(parts) == 0) return(NULL)
      data <- data.table::rbindlist(parts, fill = TRUE)

      # Add common fallback columns
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

      ancs <- sort(unique(data$ancestry))
      showNotification(
        paste0("Loaded ", study, " (", paste(ancs, collapse = ", "), ")"),
        type = "message", duration = 3
      )
      data
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
    
    # Study metadata for the current study. Always returns the registry
    # entry (never NULL for a valid study). During transition, callers
    # that check `!is.null(current_virtual_info())` to detect multi-
    # ancestry will still work because this returns NULL for unknown studies.
    current_study_info <- reactive({
      study <- current_study()
      if (is.null(study)) return(NULL)
      DEFAULT_STUDY_REGISTRY[[study]]
    })

    # Backward compat: returns NULL for single-ancestry studies
    # (will be removed once all callers migrate to current_study_info)
    current_virtual_info <- reactive({
      info <- current_study_info()
      if (is.null(info) || length(info$ancestries) <= 1) return(NULL)
      info
    })

    # Filter data for selected region - full (no max_traits limit)
    # Used by convergence view which needs to see all cis molecular colocs
    filtered_region_data <- reactive({
      req(coloc_data(), current_region())

      rk <- parse_region_key(current_region())
      if (is.null(rk)) return(coloc_data()[0])
      chr <- rk$chr; start_pos <- rk$start; end_pos <- rk$stop

      # Virtual multi-ancestry: widen the match to the consensus window so
      # colocs from all ancestries that hit the same locus are included.
      vinfo <- current_virtual_info()
      if (!is.null(vinfo)) {
        bundles <- load_multi_region_bundles(current_region(), vinfo)
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
      # PP.H4 and nlog10P filters used to live here, driven by sidebar
      # sliders. The sliders are gone - the DT trait list has per-column
      # range filters on both, and the atlas is already pre-filtered to
      # PP.H4 >= 0.8 at extraction time, so the sidebar filters were
      # redundant. source_study filter is gone because the Include
      # Studies sidebar block is gone (it duplicated the category
      # cards).
      dt
    })

    # Historically this capped rows at `input$max_traits` so the Network
    # tab could render large regions without choking. With the Network
    # view gone, the uncapped `filtered_region_data()` is exactly what
    # every downstream consumer wants, so this is just an alias kept for
    # backwards compatibility with the existing callsites.
    filtered_data <- reactive({
      filtered_region_data()
    })

    # === Convergence view ===

    # Reactive: current region coords (chr, start, end)
    conv_region_coords <- reactive({
      req(current_region())
      rk <- parse_region_key(current_region())
      if (is.null(rk)) return(NULL)
      list(chr = rk$chr, start = rk$start, end = rk$stop)
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

    # === Panel 3: category overview cards ===

    # Category counts for current region
    conv_cat_counts <- reactive({
      dt <- filtered_region_data()
      count_by_category(dt)
    })

    # Horizontal strip of compact category buttons. Each button shows
    # icon + short label + count. The full long label is available as
    # the native title attribute for hover (0 perf cost, no delay for
    # assistive tech users).
    output$conv_category_cards <- renderUI({
      counts <- conv_cat_counts()
      selected <- conv_selected_cat()

      cards <- lapply(names(TRAIT_CATEGORIES), function(cat_id) {
        meta <- TRAIT_CATEGORIES[[cat_id]]
        n <- tryCatch(counts[[cat_id]], error = function(e) 0L)
        if (is.null(n) || length(n) == 0 || is.na(n)) n <- 0L
        is_empty <- n == 0
        is_active <- !is.null(selected) && selected == cat_id

        border_style <- if (is_active) {
          paste0("border: 2px solid ", meta$color, ";")
        } else {
          "border: 1px solid #ddd;"
        }
        bg_style <- if (is_active) {
          paste0("background: ", meta$color, "15;")
        } else {
          "background: white;"
        }
        opacity_style <- if (is_empty) "opacity: 0.45;" else ""
        cursor_style <- if (is_empty) "cursor: default;" else "cursor: pointer;"

        short <- if (!is.null(meta$short)) meta$short else meta$label

        tags$button(
          id = paste0("conv_cat_btn_", cat_id),
          class = "action-button",
          type = "button",
          title = meta$label,
          disabled = if (is_empty) NA else NULL,
          style = paste0(
            "flex: 1 1 0; min-width: 0; padding: 3px 6px; ",
            "border-radius: 6px; ",
            border_style, bg_style, opacity_style, cursor_style,
            "display: flex; align-items: center; justify-content: center; ",
            "gap: 4px; text-align: center; font-family: inherit; ",
            "line-height: 1.15; white-space: nowrap; overflow: hidden;"
          ),
          tags$span(style = "font-size: 14px;", HTML(meta$icon)),
          tags$span(style = "font-size: 11px; font-weight: 600; color: #2c3e50; overflow: hidden; text-overflow: ellipsis;",
                    short),
          tags$span(style = paste0("font-size: 12px; font-weight: 700; color: ",
                                   if (is_empty) "#aaa" else "#2c3e50", ";"),
                    n)
        )
      })

      # 3x2 grid so each button gets ~double the width of a single row
      # of 6 - fits the full short label ("Metabolites", "PheWAS")
      # without truncating.
      div(style = paste0(
        "display: grid; grid-template-columns: repeat(3, 1fr); ",
        "gap: 4px; margin: 4px 0 6px 0;"),
        cards)
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

    # Drill-down header: shows the currently selected category above the
    # DT trait list. "Clear" is implicit - the user clicks a different
    # category card in the left sidebar.
    output$conv_drilldown_header <- renderUI({
      sel <- conv_selected_cat()
      if (is.null(sel)) {
        return(div(style = "font-size: 11px; color: #999; margin: 4px 0 6px 0; font-style: italic;",
                   "Select a category in the left sidebar."))
      }
      meta <- TRAIT_CATEGORIES[[sel]]
      counts <- conv_cat_counts()
      n_total <- tryCatch(counts[[sel]], error = function(e) 0L)
      if (is.null(n_total) || length(n_total) == 0 || is.na(n_total)) n_total <- 0L
      div(style = "margin: 2px 0 6px 0;",
          h6(style = "margin: 0; color: #2c3e50;",
             meta$icon, " ", meta$label,
             tags$span(style = "font-size: 11px; color: #888; margin-left: 6px;",
                       sprintf("(%d)", n_total))))
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
          make_bundle_key(source_study, sumstats_2_file, ancestry,
                          is_virtual = TRUE))]
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

    # Flat data.frame shape for the trait list DT + a parallel
    # bundle-key vector (stable across sort/filter operations on the
    # client, because DT's server-side selection maps back to original
    # row indices). Rebuilt whenever the drilldown data changes.
    conv_trait_table <- reactive({
      dt <- conv_drilldown_data()
      if (is.null(dt) || nrow(dt) == 0) {
        return(list(df = NULL, keys = character(0)))
      }
      is_virt_dd <- is_virtual_study(current_study()) &&
                    "consensus_key" %in% names(dt)

      trait_name_fn <- function(row) {
        nm <- tryCatch(get_trait_display_name(row),
                       error = function(e) basename(row$sumstats_2_file))
        if (is.null(nm) || is.na(nm) || !nzchar(nm)) basename(row$sumstats_2_file) else nm
      }
      names_vec <- vapply(seq_len(nrow(dt)),
                          function(i) trait_name_fn(dt[i]),
                          character(1))

      badge_vec <- if (is_virt_dd) {
        as.character(dt$ancestries)
      } else {
        as.character(dt$source_study)
      }

      nlp_vec <- if ("sumstats_2_max_nlog10P" %in% names(dt)) {
        round(as.numeric(dt$sumstats_2_max_nlog10P), 1)
      } else rep(NA_real_, nrow(dt))

      pph4_vec <- if ("PP.H4.abf" %in% names(dt)) {
        round(as.numeric(dt$PP.H4.abf), 3)
      } else rep(NA_real_, nrow(dt))

      keys <- if (is_virt_dd) {
        as.character(dt$consensus_key)
      } else {
        vapply(seq_len(nrow(dt)), function(i)
          make_bundle_key(dt$source_study[i], dt$sumstats_2_file[i]),
          character(1))
      }

      list(
        df = data.frame(
          Trait = names_vec,
          Ancestries = badge_vec,
          nlog10P = nlp_vec,
          PP.H4 = pph4_vec,
          check.names = FALSE,
          stringsAsFactors = FALSE
        ),
        keys = keys
      )
    })

    # Render the DT trait list. Single-row selection feeds
    # conv_selected_trait via the rows_selected observer below.
    output$conv_trait_list <- DT::renderDataTable({
      tbl <- conv_trait_table()
      if (is.null(tbl$df) || nrow(tbl$df) == 0) {
        return(DT::datatable(
          data.frame(Trait = character(0), Ancestries = character(0),
                     nlog10P = numeric(0), PP.H4 = numeric(0),
                     check.names = FALSE),
          selection = "none",
          options = list(dom = 't',
                         language = list(emptyTable = "Pick a category above to see traits."))
        ))
      }

      # Preselect the current conv_selected_trait if it maps back to a
      # row in this table (otherwise DT defaults to nothing).
      sel_key <- conv_selected_trait()
      sel_row <- if (!is.null(sel_key)) which(tbl$keys == sel_key) else integer(0)
      sel_row <- sel_row[1]

      DT::datatable(
        tbl$df,
        selection = list(mode = "single",
                         selected = if (length(sel_row) && !is.na(sel_row)) sel_row else NULL,
                         target = "row"),
        rownames = FALSE,
        filter = list(position = "top", clear = FALSE, plain = TRUE),
        class = "compact stripe hover",
        options = list(
          dom = 't',             # table only (per-column filters above)
          pageLength = -1,       # show all rows
          scrollY = "540px",
          scrollCollapse = TRUE,
          paging = FALSE,
          order = list(list(2, 'desc')),  # default sort by nlog10P
          autoWidth = FALSE,
          columnDefs = list(
            list(className = "dt-right", targets = c(2, 3)),
            list(width = "50%", targets = 0),
            list(width = "22%", targets = 1),
            list(width = "14%", targets = 2),
            list(width = "14%", targets = 3)
          )
        )
      )
    }, server = TRUE)

    # DT row-selection -> consensus/bundle key.
    observeEvent(input$conv_trait_list_rows_selected, {
      row <- input$conv_trait_list_rows_selected
      if (is.null(row) || length(row) == 0) return()
      keys <- conv_trait_table()$keys
      if (length(keys) >= row[1]) {
        conv_selected_trait(keys[row[1]])
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    # On region change: reset per-region state (gene info panel,
    # currently-selected category and trait), then auto-select the
    # first non-empty category. The trait inside that category is
    # picked by a separate observer that reacts to conv_drilldown_data.
    observeEvent(current_region(), {
      selected_gene(NULL)
      conv_selected_trait(NULL)

      counts <- conv_cat_counts()
      non_empty <- names(counts)[counts > 0]
      if (length(non_empty) == 0) {
        conv_selected_cat(NULL)
        return()
      }
      conv_selected_cat(non_empty[1])
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
        make_bundle_key(dt$source_study[1], dt$sumstats_2_file[1])
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
      req(regional_data_path(), current_region(), conv_selected_trait())
      sel <- conv_selected_trait()

      vinfo <- current_virtual_info()
      if (!is.null(vinfo)) {
        bundles <- load_multi_region_bundles(current_region(), vinfo)
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

      load_trait_sumstats(regional_data_path(), current_region(), sel)
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
      vinfo <- current_virtual_info()
      xr <- if (!is.null(vinfo)) {
        bundles <- load_multi_region_bundles(current_region(), vinfo)
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

    # Download handler
    output$download_data <- downloadHandler(
      filename = function() {
        paste0("coloc_", gsub(":", "_", current_region()), "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        fwrite(filtered_data(), file)
      }
    )

    # === Regional bundle data ===

    # Path to the bundled regional RDS directory for the current study.
    # For virtual multi-ancestry studies: returns the ancestry regional
    # dir that contains the current representative region, fallback to
    # the first available ancestry dir. For single-ancestry studies:
    # returns the multi-category path first, fallback to single-category.
    regional_data_path <- reactive({
      req(current_study())
      study <- current_study()

      # Virtual multi-ancestry study: return the first ancestry regional
      # dir that actually contains the selected representative region's
      # RDS file (so `load_region_bundle` picks the right one). Fall
      # back to the first existing dir.
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

      # Single-ancestry study. Category layout first (the current atlas
      # layout), then single-category fallback. The "legacy" regional_plots
      # layout is gone.
      # TODO: collapse once every study is routed through the virtual
      # multi-ancestry codepath (single-ancestry studies with
      # ancestries = c("EUR") or similar).
      cats <- attr(available_studies, "categories")
      if (!is.null(cats) && !is.null(cats[[study]])) {
        cat_path <- file.path(DATA_PATH, cats[[study]], "regional", study)
        if (dir.exists(cat_path)) return(cat_path)
      }
      bundled_path <- file.path(DATA_PATH, "regional", study)
      if (dir.exists(bundled_path)) return(bundled_path)

      NULL
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

    # Load base study sumstats for current region. For virtual studies,
    # returns a named list(ancestry = dt) for overlay rendering; otherwise
    # returns a single data.table.
    regional_base_sumstats <- reactive({
      req(regional_data_path(), current_region())
      vinfo <- current_virtual_info()
      if (!is.null(vinfo)) {
        bundles <- load_multi_region_bundles(current_region(), vinfo)
        if (is.null(bundles) || length(bundles) == 0) return(NULL)
        return(lapply(bundles, `[[`, "base"))
      }
      load_base_sumstats(regional_data_path(), current_region())
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
