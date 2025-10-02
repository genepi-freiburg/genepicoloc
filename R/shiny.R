# shiny visualization for genepicoloc

#' Launch GWAS Colocalization Network Viewer
#' 
#' @param port Port number (default: 3838)
#' @param study_base_path Base path for study data
#' @param available_studies Named list of available studies or path to file containing the list
#' @param enable_llm Enable AI/LLM features (default: FALSE)
#' @export
launch_coloc_viewer <- function(port = 3838,
                                STUDY_BASE_PATH,
                                available_studies) {

  # visualization
  library(shiny)
  library(visNetwork)
  
  # data processing
  library(jsonlite)
  library(data.table)
  library(DT)
  
  # disable AI features for initial deployment
  ENABLE_LLM <- FALSE
  if (ENABLE_LLM) {
    rag_scripts <- ""
    invisible(sapply(list.files(rag_scripts, full.names = T), source))
    force_recreate <- TRUE
    library(ellmer)
    library(ragnar)
    ### Initialize RAG system once
    dir_rag <- ""
    ragnar_model <- "all-minilm:latest_t20"
    ollama_model <- "llama3.2:3b_t20" # "gemma3:1b_t10"
    STORE_PATH <- setup_genomics_rag(ragnar_model = ragnar_model,
                                     dir_rag = dir_rag,
                                     force_recreate = force_recreate)
  } else {
    setup_genomics_rag <- function(...) { return(NULL) }
    interpret_with_ragnar <- function(...) {
      return("AI features disabled for testing")
    }
  }
  
  # Define color palette for studies
  study_colors <- c(
    "CKDGen_r5" = "#1f77b4",
    "FinnGen_r9" = "#ff7f0e", 
    "GCKD_mGWAS_plasma" = "#2ca02c",
    "GCKD_mGWAS_urine" = "#d62728",
    "Icelanders_pGWAS" = "#9467bd",
    "MVP_R4" = "#8c564b",
    "UKB_PPP_EUR" = "#e377c2",
    "UKB_TOPMed" = "#7f7f7f",
    "pho_ca" = "#bcbd22",
    "eQTLGen" = "#17becf",
    "Kidney_eQTL" = "#8dd3c7",
    "GTEXv8_eQTL" = "#8dd0d7"
  )
  
  # UI ----
  ui <- fluidPage(
    tags$head(
      tags$style(HTML("
    .study-card {
      border: 1px solid #ddd;
      border-radius: 8px;
      padding: 15px;
      margin: 10px;
      cursor: pointer;
      transition: all 0.3s;
      background-color: white;
      box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    .study-card:hover {
      box-shadow: 0 4px 8px rgba(0,0,0,0.2);
      transform: translateY(-2px);
    }
    .study-card h4 {
      margin-top: 0;
      color: #2c3e50;
    }
    .study-card p {
      color: #7f8c8d;
      margin-bottom: 10px;
    }
    .category-header {
      color: #34495e;
      border-bottom: 2px solid #3498db;
      padding-bottom: 10px;
      margin-top: 30px;
      margin-bottom: 20px;
    }
  "))
    ),
    
    tabsetPanel(
      id = "main_tabs",
      tabPanel("Home", 
               fluidRow(
                 column(12,
                        h1("GWAS Colocalization Network Viewer", style = "text-align: center; margin-bottom: 30px;"),
                        h3("Select a Study to Analyze", style = "text-align: center; color: #7f8c8d; margin-bottom: 40px;")
                 )
               ),
               
               fluidRow(
                 column(12,
                        h3("Available Studies", class = "category-header"),
                        p("Click on any study below to load and analyze its colocalization data:", 
                          style = "text-align: center; margin-bottom: 30px;")
                 )
               ),
               
               # Create study cards in a grid layout
               fluidRow(
                 lapply(names(available_studies), function(study_name) {
                   column(4,
                          div(class = "study-card",
                              onclick = paste0("Shiny.setInputValue('selected_study', '", study_name, "', {priority: 'event'});"),
                              h4(study_name),
                              p(style = "font-size: 12px; font-style: italic;", 
                                paste("Folder:", available_studies[[study_name]])),
                              tags$button(class = "btn btn-primary btn-sm", "Load Study")
                          )
                   )
                 })
               ),
               
               hr(),
               
               fluidRow(
                 column(12, style = "text-align: center; margin-top: 30px;",
                        h4("Or upload your own data:"),
                        actionButton("go_to_analysis_manual", "Go to Manual Upload", 
                                     class = "btn-success btn-lg", icon = icon("upload"))
                 )
               )
      ),
      tabPanel("Analysis",
               titlePanel("GWAS Colocalization Network Viewer"),
               
               # Add notification area for study loading
               fluidRow(
                 column(12,
                        uiOutput("study_notification")
                 )
               ),
               
               fluidRow(
                 # Left sidebar panel
                 column(
                   width = 3,
                   wellPanel(
                     # Current study indicator
                     uiOutput("current_study_info"),
                     
                     # File upload - KEEP ONLY THIS ONE
                     fileInput("upload_file", 
                               "Or upload different data (RDS/TSV/CSV):",
                               accept = c(".tsv", ".csv", ".txt", ".gz", ".RDS")),
                     
                     # Or use demo data
                     actionButton("use_demo", "Load Demo Data", 
                                  class = "btn-warning", width = "100%"),
                     
                     # Back to study selection
                     actionButton("back_to_home", "Select Different Study", 
                                  class = "btn-info", width = "100%", icon = icon("arrow-left")),
                     
                     hr(),
                     
                     # Region selector
                     selectInput("selected_region",
                                 "Select Region (by nearest gene):",
                                 choices = NULL,
                                 selected = NULL),
                     
                     # Alternative region selector by coordinates
                     selectInput("selected_region_coord",
                                 "Or select by coordinates:",
                                 choices = NULL,
                                 selected = NULL),
                     
                     hr(),
                     
                     # Filter options
                     sliderInput("min_pph4",
                                 "Minimum PP.H4:",
                                 min = 0.5, max = 1, value = 0.5, step = 0.01),
                     
                     sliderInput("min_nlog10P",
                                 "Minimum -log10(p-value):",
                                 min = 6, max = 100, value = -log10(1e-7), step = 0.5),
                     
                     checkboxGroupInput("filter_studies",
                                        "Include Studies:",
                                        choices = names(study_colors),
                                        selected = names(study_colors)),
                     
                     # Select all/none buttons
                     fluidRow(
                       column(6, actionButton("select_all", "Select All", class = "btn-sm")),
                       column(6, actionButton("select_none", "Select None", class = "btn-sm"))
                     ),
                     
                     hr(),
                     
                     # Color legend
                     h5("Study Colors"),
                     htmlOutput("color_legend"),
                     
                     hr(),
                     
                     # Display options
                     checkboxInput("show_labels", "Show node labels", value = TRUE),
                     checkboxInput("physics", "Enable physics simulation", value = FALSE),
                     
                     hr(),
                     
                     # Region info
                     h4("Region Information"),
                     verbatimTextOutput("region_info"),
                     
                     hr(),
                     
                     # Download button
                     downloadButton("download_data", "Download filtered data")
                   )  # End wellPanel
                 ),  # End left column
                 
                 # Main panel
                 column(
                   width = 6,
                   tabsetPanel(
                     tabPanel("Network", 
                              visNetworkOutput("network", height = "600px"),
                              br(),
                              h4("Selected Node Details"),
                              verbatimTextOutput("node_info")),
                     
                     tabPanel("Data Table",
                              br(),
                              DTOutput("data_table")),
                     
                     tabPanel("Summary",
                              br(),
                              h4("Colocalization Summary for Selected Region"),
                              plotOutput("summary_plot", height = "400px"),
                              br(),
                              verbatimTextOutput("summary_stats"))
                   )  # End tabsetPanel
                 ),  # End main column
                 
                 # Right LLM panel
                 column(
                   width = 3,
                   wellPanel(
                     h4("AI Assistant"),
                     actionButton("llm_interpret", "AI Interpret Region", 
                                  class = "btn-success", width = "100%",
                                  icon = icon("brain")),
                     br(),
                     br(),
                     textAreaInput("llm_query", 
                                   "Ask a follow-up question:", 
                                   placeholder = "e.g., What diseases is this gene associated with?",
                                   rows = 3),
                     actionButton("llm_submit", "Ask", class = "btn-primary", width = "100%"),
                     hr(),
                     h5("Response:"),
                     uiOutput("llm_response")
                   )  # End wellPanel
                 )  # End right column
               )  # End fluidRow
      ),
      tabPanel("Documentation",
               fluidRow(
                 column(10, offset = 1,
                        h2("GWAS Colocalization Network Viewer Documentation"),
                        br(),
                        
                        # Add the visualization section
                        h3("The Challenge: From GWAS Hits to Biological Insights"),
                        plotOutput("workflow_diagram", height = "400px"),
                        br(),
                        
                        h3("Getting Started"),
                        p("This application visualizes GWAS colocalization results across multiple studies."),
                        
                        h4("1. Study Selection"),
                        p("On the home page, select a pre-configured study by clicking its card, or click 'Go to Manual Upload' to upload your own data."),
                        
                        h4("2. Data Format"),
                        p("The app accepts RDS files containing colocalization results with the following required columns:"),
                        tags$ul(
                          tags$li("CHR_var, BP_START_var, BP_STOP_var - genomic coordinates"),
                          tags$li("nearest_gene_1 - nearest gene symbol"),
                          tags$li("PP.H4.abf - posterior probability of colocalization"),
                          tags$li("source_study - study identifier"),
                          tags$li("sumstats_2_max_nlog10P - significance measure")
                        ),
                        
                        h4("3. Network Visualization"),
                        p("The network displays colocalizations with:"),
                        tags$ul(
                          tags$li("Central star node - the genomic region"),
                          tags$li("Surrounding nodes - colocalized traits/studies"),
                          tags$li("Node size - proportional to PP.H4 value"),
                          tags$li("Node color - indicates source study")
                        ),
                        
                        h4("4. Filtering Options"),
                        tags$ul(
                          tags$li("PP.H4 threshold - filter by colocalization strength"),
                          tags$li("Study selection - show/hide specific studies"),
                          tags$li("Region selection - focus on specific genomic regions")
                        ),
                        
                        h4("5. Export"),
                        p("Click 'Download filtered data' to export the current view as CSV."),
                        
                        br(),
                        h3("Technical Details"),
                        p("For questions or issues, please contact the development team."),
                        br()
                 )
               )
      )
    )
  )
  
  server <- function(input, output, session) {
    
    # === Reactive Values ===
    
    # Reactive value to store the current loaded study
    current_study <- reactiveVal(NULL)
    
    # Add this new reactive value to track if we should use uploaded data
    use_uploaded_data <- reactiveVal(FALSE)
    
    # === Event Handlers ===
    
    # Handle study selection from home page
    observeEvent(input$selected_study, {
      req(input$selected_study)
      
      # Get the folder name from the study name
      folder_name <- available_studies[[input$selected_study]]
      
      # Construct file path
      file_path <- file.path(STUDY_BASE_PATH, folder_name, "annot", "annot_filt.RDS")
      
      # Check if file exists
      if (file.exists(file_path)) {
        # Store current study (now storing the name, not the folder)
        current_study(input$selected_study)
        
        # Clear the uploaded data flag
        use_uploaded_data(FALSE)
        
        # Switch to Analysis tab
        updateTabsetPanel(session, "main_tabs", selected = "Analysis")
        
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
    
    # Handle manual navigation to analysis
    observeEvent(input$go_to_analysis_manual, {
      current_study(NULL)
      updateTabsetPanel(session, "main_tabs", selected = "Analysis")
    })
    
    # Handle back to home button
    observeEvent(input$back_to_home, {
      updateTabsetPanel(session, "main_tabs", selected = "Home")
    })
    
    # Clear current study when user uploads a file
    observeEvent(input$upload_file, {
      if (!is.null(input$upload_file)) {
        current_study(NULL)  # Clear the auto-loaded study
        use_uploaded_data(TRUE)  # Set flag to use uploaded data
        showNotification(
          paste("Loaded uploaded file:", input$upload_file$name),
          type = "message",
          duration = 3
        )
      }
    })
    
    # Select all studies
    observeEvent(input$select_all, {
      updateCheckboxGroupInput(session, "filter_studies", 
                               selected = names(study_colors))
    })
    
    # Select no studies
    observeEvent(input$select_none, {
      updateCheckboxGroupInput(session, "filter_studies", 
                               selected = character(0))
    })
    
    # Update region selectors
    observe({
      req(regions_data())
      regions <- regions_data()
      
      # Update gene-based selector
      gene_choices <- setNames(
        paste0(regions$CHR_var, ":", regions$BP_START_var, "-", regions$BP_STOP_var),
        paste0(regions$nearest_gene_1, " (", regions$N, " coloc)")
      )
      updateSelectInput(session, "selected_region", choices = gene_choices)
      
      # Update coordinate-based selector  
      coord_choices <- setNames(
        paste0(regions$CHR_var, ":", regions$BP_START_var, "-", regions$BP_STOP_var),
        regions$region_id
      )
      updateSelectInput(session, "selected_region_coord", choices = coord_choices)
    })
    
    # Synchronize selectors
    observeEvent(input$selected_region, {
      if (!is.null(input$selected_region) && input$selected_region != "") {
        updateSelectInput(session, "selected_region_coord", selected = input$selected_region)
      }
    })
    
    observeEvent(input$selected_region_coord, {
      if (!is.null(input$selected_region_coord) && input$selected_region_coord != "") {
        updateSelectInput(session, "selected_region", selected = input$selected_region_coord)
      }
    })
    
    # === Reactive Data Processing ===
    
    # Get current selected region
    current_region <- reactive({
      if (!is.null(input$selected_region) && input$selected_region != "") {
        input$selected_region
      } else if (!is.null(input$selected_region_coord) && input$selected_region_coord != "") {
        input$selected_region_coord
      } else {
        NULL
      }
    })
    
    # Load and process data
    coloc_data <- reactive({
      # First check if we should use uploaded data
      if (use_uploaded_data() && !is.null(input$upload_file)) {
        ext <- tools::file_ext(input$upload_file$datapath)
        if (ext == "gz") {
          return(fread(input$upload_file$datapath))
        } else if (ext == "RDS") {
          return(readRDS(input$upload_file$datapath))
        }
      }
      
      # Second priority: auto-loaded study from home page
      if (!is.null(current_study())) {
        folder_name <- available_studies[[current_study()]]
        file_path <- file.path(STUDY_BASE_PATH, folder_name, "annot", "annot_filt.RDS")
        if (file.exists(file_path)) {
          data <- readRDS(file_path)
          showNotification(
            paste("Successfully loaded:", current_study()),
            type = "message",
            duration = 3
          )
          return(data)
        }
      }
      
      # Third priority: demo data
      if (input$use_demo > 0) {
        demo_file <- ""
        if (file.exists(demo_file)) {
          return(fread(demo_file))
        } else {
          showNotification("Demo file not found at specified path", type = "error")
          return(NULL)
        }
      }
      
      return(NULL)
    })
    
    # Get unique regions
    regions_data <- reactive({
      req(coloc_data())
      dt <- coloc_data()
      
      # Get unique regions with their info
      regions <- unique(dt[, .(
        CHR_var, BP_START_var, BP_STOP_var,
        region_center_pos, nearest_gene_1,
        region_id = paste0("chr", CHR_var, ":", BP_START_var, "-", BP_STOP_var)
      )])
      
      # Count colocalizations per region
      region_counts <- dt[, .N, by = .(CHR_var, BP_START_var, BP_STOP_var)]
      regions <- merge(regions, region_counts, by = c("CHR_var", "BP_START_var", "BP_STOP_var"))
      
      setorder(regions, CHR_var, BP_START_var)
      regions
    })
    
    # Filter data for selected region
    filtered_data <- reactive({
      req(coloc_data(), current_region())
      
      # Parse region
      region_parts <- strsplit(current_region(), ":")[[1]]
      chr <- region_parts[1]
      pos_parts <- strsplit(region_parts[2], "-")[[1]]
      start_pos <- as.numeric(pos_parts[1])
      end_pos <- as.numeric(pos_parts[2])
      
      # Filter data
      dt <- coloc_data()[CHR_var == chr & BP_START_var == start_pos & BP_STOP_var == end_pos]
      
      # Apply PP.H4 filter
      dt <- dt[PP.H4.abf >= input$min_pph4]
      
      # Apply nlog10 filter
      dt <- dt[sumstats_2_max_nlog10P >= input$min_nlog10P]
      
      # Apply study filter
      dt <- dt[source_study %in% input$filter_studies]
      
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
      region_info <- unique(dt[, .(nearest_gene_1, CHR_var, BP_START_var, BP_STOP_var)])
      central_node <- data.frame(
        id = "region",
        label = region_info$nearest_gene_1[1],
        title = paste0("Region: chr", region_info$CHR_var[1], ":", 
                       region_info$BP_START_var[1], "-", region_info$BP_STOP_var[1]),
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
                             "CKDGen_r5" = CKDGen_r5_ckdgen_r5_name,
                             "FinnGen_r9" = FinnGen_r9_phenotype,
                             "GCKD_mGWAS_plasma" = GCKD_mGWAS_plasma_BIOCHEMICAL,
                             "GCKD_mGWAS_urine" = GCKD_mGWAS_urine_BIOCHEMICAL,
                             "Icelanders_pGWAS" = Icelanders_pGWAS_Protein..short.name.,
                             "MVP_R4" = MVP_R4_Analyzed.variable,
                             "UKB_PPP_EUR" = UKB_PPP_EUR_olink_target_fullname,
                             "UKB_TOPMed" = UKB_TOPMed_phenostring,
                             "pho_ca" = pho_ca_pho_ca_name,
                             "eQTLGen" =  eQTLGen_gene_name,
                             "Kidney_eQTL" = Kidney_eQTL_gene_name,
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
          size = 10 + 20 * PP.H4.abf,  # Size by PP.H4
          color = study_colors[source_study]
        )
      }, by = .I]
      
      # Convert to data.frame with same columns as central_node
      trait_nodes <- data.frame(
        id = trait_nodes_list$id,
        label = trait_nodes_list$label,
        title = trait_nodes_list$title,
        group = trait_nodes_list$group,
        shape = trait_nodes_list$shape,
        size = trait_nodes_list$size,
        color = trait_nodes_list$color,
        stringsAsFactors = FALSE
      )
      
      nodes <- rbind(central_node, trait_nodes)
      
      # Create edges
      edges <- data.frame(
        from = rep("region", nrow(trait_nodes)),
        to = trait_nodes$id,
        width = 2,  # Fixed width since PP.H4 is pre-filtered
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
          h5("Currently Loaded:", style = "margin: 0;"),
          p(current_study(), style = "margin: 5px 0; font-weight: bold;"),
          p(paste("Folder:", available_studies[[current_study()]]), style = "margin: 0; font-size: 11px; color: #666;")
        )
      } else {
        div(
          style = "background-color: #f5f5f5; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
          p("No study auto-loaded. Upload your own data below.", style = "margin: 0; color: #666;")
        )
      }
    })
    
    # Create color legend
    output$color_legend <- renderUI({
      legend_html <- paste0(
        lapply(names(study_colors), function(study) {
          paste0(
            '<div style="margin: 2px 0;">',
            '<span style="display: inline-block; width: 12px; height: 12px; ',
            'background-color: ', study_colors[study], '; margin-right: 5px;"></span>',
            '<span style="font-size: 11px;">', study, '</span>',
            '</div>'
          )
        }),
        collapse = ""
      )
      HTML(legend_html)
    })
    
    # Render network
    output$network <- renderVisNetwork({
      net_data <- network_data()
      
      if (nrow(net_data$nodes) == 0) {
        visNetwork(data.frame(id = 1, label = "No data"), data.frame()) %>%
          visNodes(shape = "text", font = list(size = 20))
      } else {
        visNetwork(net_data$nodes, net_data$edges) %>%
          visOptions(
            highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
            selectedBy = "group",
            nodesIdSelection = TRUE
          ) %>%
          visPhysics(enabled = input$physics, stabilization = FALSE) %>%
          visLayout(randomSeed = 123) %>%
          visInteraction(
            dragNodes = TRUE,
            dragView = TRUE,
            zoomView = TRUE
          ) %>%
          visEvents(select = "function(nodes) {
          Shiny.setInputValue('selected_node', nodes.nodes[0]);
        }")
      }
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
    
    # Region info
    output$region_info <- renderPrint({
      req(filtered_data())
      dt <- filtered_data()
      
      cat("Chromosome:", unique(dt$CHR_var), "\n")
      cat("Start:", format(unique(dt$BP_START_var), big.mark = ","), "\n")
      cat("End:", format(unique(dt$BP_STOP_var), big.mark = ","), "\n")
      cat("Nearest gene:", unique(dt$nearest_gene_1), "\n")
      cat("Total colocalizations:", nrow(dt), "\n")
      cat("Studies represented:", length(unique(dt$source_study)), "\n")
    })
    
    # Data table
    output$data_table <- renderDT({
      req(filtered_data())
      
      # Select key columns for display
      display_cols <- c("source_study", "nearest_gene_1", "PP.H4.abf", 
                        "sumstats_2_max_nlog10P", "directionality")
      
      # Check which columns exist
      available_cols <- intersect(display_cols, names(filtered_data()))
      
      dt_display <- filtered_data()[, ..available_cols]
      
      datatable(dt_display,
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE) %>%
        formatRound(c("PP.H4.abf", "sumstats_2_max_nlog10P"), 3)
    })
    
    # Summary plot
    output$summary_plot <- renderPlot({
      req(filtered_data())
      dt <- filtered_data()
      
      # Count by study
      study_counts <- dt[, .N, by = source_study]
      
      par(mar = c(8, 4, 2, 2))
      barplot(
        study_counts$N,
        names.arg = study_counts$source_study,
        col = study_colors[study_counts$source_study],
        las = 2,
        main = "Colocalizations by Study",
        ylab = "Count"
      )
    })
    
    # Summary statistics
    output$summary_stats <- renderPrint({
      req(filtered_data())
      dt <- filtered_data()
      
      cat("Summary Statistics:\n")
      cat("==================\n\n")
      cat("Mean PP.H4:", round(mean(dt$PP.H4.abf), 3), "\n")
      cat("Median PP.H4:", round(median(dt$PP.H4.abf), 3), "\n")
      cat("Range PP.H4:", round(min(dt$PP.H4.abf), 3), "-", 
          round(max(dt$PP.H4.abf), 3), "\n\n")
      
      cat("Colocalizations by study:\n")
      print(dt[, .N, by = source_study][order(-N)])
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
    
    # LLM interpret button handler
    observeEvent(input$llm_interpret, {
      req(filtered_data())
      
      # Show loading message
      output$llm_response <- renderUI({
        tags$div(
          class = "text-muted",
          tags$i(class = "fa fa-spinner fa-spin"),
          " Analyzing with RAG system..."
        )
      })
      
      # Prepare context
      dt <- filtered_data()
      region_info <- unique(dt[, .(nearest_gene_1, CHR_var, BP_START_var, BP_STOP_var)])
      
      # Get top traits
      top_traits <- dt[order(-PP.H4.abf)][1:min(5, nrow(dt)), ]
      trait_names <- sapply(1:nrow(top_traits), function(i) {
        trait_name <- switch(top_traits$source_study[i],
                             "CKDGen_r5" = top_traits$CKDGen_r5_ckdgen_r5_name[i],
                             "FinnGen_r9" = top_traits$FinnGen_r9_phenotype[i],
                             "GCKD_mGWAS_plasma" = top_traits$GCKD_mGWAS_plasma_BIOCHEMICAL[i],
                             "GCKD_mGWAS_urine" = top_traits$GCKD_mGWAS_urine_BIOCHEMICAL[i],
                             "Icelanders_pGWAS" = top_traits$Icelanders_pGWAS_Protein..short.name.[i],
                             "MVP_R4" = top_traits$MVP_R4_Analyzed.variable[i],
                             "UKB_PPP_EUR" = top_traits$UKB_PPP_EUR_olink_target_fullname[i],
                             "UKB_TOPMed" = top_traits$UKB_TOPMed_phenostring[i],
                             "pho_ca" = top_traits$pho_ca_pho_ca_name[i],
                             "eQTLGen" = top_traits$eQTLGen_gene_name[i],
                             "Kidney_eQTL" = top_traits$Kidney_eQTL_gene_name[i],
                             "Unknown trait")
        if (is.na(trait_name) || trait_name == "") trait_name <- "Unknown trait"
        return(trait_name)
      })
      
      trait_names <- unique(trait_names[!is.na(trait_names)])
      
      additional_context <- paste0(
        "GWAS colocalization analysis for chr", region_info$CHR_var[1], ":", 
        format(region_info$BP_START_var[1], big.mark = ","), "-", 
        format(region_info$BP_STOP_var[1], big.mark = ","), 
        " with ", nrow(dt), " colocalizations (mean PP.H4=", 
        round(mean(dt$PP.H4.abf), 3), ")."
      )
      
      # Call RAG system
      response <- interpret_with_ragnar(
        gene = region_info$nearest_gene_1[1],
        traits = trait_names,
        additional_context = additional_context,
        store_path = STORE_PATH,
        ollama_model = ollama_model
      )
      
      # Display response
      output$llm_response <- renderUI({
        tags$div(
          style = "max-height: 400px; overflow-y: auto; padding: 10px; background-color: #f5f5f5; border-radius: 5px;",
          HTML(gsub("\n", "<br>", response))
        )
      })
    })
    
    # Custom query button
    observeEvent(input$llm_submit, {
      req(input$llm_query, filtered_data())
      
      dt <- filtered_data()
      region_info <- unique(dt[, .(nearest_gene_1)])
      
      print(input$llm_query)
      # Use RAG for custom queries too
      response <- interpret_with_ragnar(
        gene = region_info$nearest_gene_1[1],
        traits = paste0("Question: ", input$llm_query),
        additional_context = input$llm_query,
        store_path = STORE_PATH,
        ollama_model = ollama_model
      )    
      
      output$llm_response <- renderUI({
        tags$div(
          style = "max-height: 400px; overflow-y: auto; padding: 10px; background-color: #f5f5f5; border-radius: 5px;",
          HTML(gsub("\n", "<br>", response))
        )
      })
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
  
  # Run app
  app <- shinyApp(ui = ui, server = server)
  cat("Starting Shiny app on port", port, "\n")
  runApp(app, port = port, host = "127.0.0.1")
  
}
