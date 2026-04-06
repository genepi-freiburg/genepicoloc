# LLM-powered region interpretation (Claude API via ellmer + shinychat)

# Build context string from colocalization results for a region
build_coloc_context <- function(region_data, gene_annotation = NULL) {
  if (is.null(region_data) || nrow(region_data) == 0) return(NULL)

  # Region info
  region_str <- paste0("chr", region_data$CHR_var[1], ":",
                        region_data$BP_START_var[1], "-", region_data$BP_STOP_var[1])
  nearest_gene <- region_data$nearest_gene_1[1]

  sections <- paste0("Genomic region: ", region_str, " (nearest gene: ", nearest_gene, ")")

  # Group colocalizations by study type
  study_groups <- list(
    eQTL = c("eQTLGen", "Kidney_eQTL", "GTEXv8_eQTL"),
    pQTL = c("UKB_PPP_EUR", "Icelanders_pGWAS"),
    mQTL = c("GCKD_mGWAS_plasma", "GCKD_mGWAS_urine"),
    PheWAS = c("MVP_R4_EUR", "MVP_R4_AFR", "MVP_R4_AMR", "MVP_R4_EAS", "MVP_R4_META",
               "FinnGen_r9", "UKB_TOPMed", "CKDGen_r4"),
    Imaging = c("UKB_kidney_vol")
  )

  for (group_name in names(study_groups)) {
    studies <- study_groups[[group_name]]
    group_data <- region_data[region_data$source_study %in% studies, ]
    if (nrow(group_data) == 0) next

    # Extract trait names based on study
    traits <- sapply(seq_len(nrow(group_data)), function(i) {
      row <- group_data[i, ]
      study <- row$source_study
      trait <- tryCatch({
        if (study == "eQTLGen" && "eQTLGen_gene_name" %in% names(row))
          row$eQTLGen_gene_name
        else if (study == "Kidney_eQTL" && "Kidney_eQTL_gene_name" %in% names(row))
          row$Kidney_eQTL_gene_name
        else if (study == "UKB_PPP_EUR" && "UKB_PPP_EUR_olink_target_fullname" %in% names(row))
          row$UKB_PPP_EUR_olink_target_fullname
        else if (study == "Icelanders_pGWAS" && "Icelanders_pGWAS_Protein..full.name." %in% names(row))
          row$`Icelanders_pGWAS_Protein..full.name.`
        else if (study == "FinnGen_r9" && "FinnGen_r9_phenotype" %in% names(row))
          row$FinnGen_r9_phenotype
        else if (grepl("^MVP_R4", study)) {
          col <- paste0(study, "_Analyzed.variable")
          if (col %in% names(row)) row[[col]] else basename(row$sumstats_2_file)
        }
        else if (study == "UKB_TOPMed" && "UKB_TOPMed_phenostring" %in% names(row))
          row$UKB_TOPMed_phenostring
        else if (study == "CKDGen_r4" && "CKDGen_r4_Name" %in% names(row))
          row$CKDGen_r4_Name
        else if (study == "UKB_kidney_vol")
          gsub("model1_qnorm_(.*)_chr.*", "\\1", basename(row$sumstats_2_file))
        else if (study == "GCKD_mGWAS_plasma" && "GCKD_mGWAS_plasma_BIOCHEMICAL" %in% names(row))
          row$GCKD_mGWAS_plasma_BIOCHEMICAL
        else if (study == "GCKD_mGWAS_urine" && "GCKD_mGWAS_urine_BIOCHEMICAL" %in% names(row))
          row$GCKD_mGWAS_urine_BIOCHEMICAL
        else
          basename(row$sumstats_2_file)
      }, error = function(e) "Unknown")
      paste0(trait, " (PP.H4=", round(row$PP.H4.abf, 2), ")")
    })

    # Limit to top 20 per group
    if (length(traits) > 20) {
      traits <- c(head(traits, 20), paste0("... and ", length(traits) - 20, " more"))
    }

    studies_used <- paste(unique(group_data$source_study), collapse = ", ")
    sections <- c(sections, paste0(
      group_name, " colocalizations (", studies_used, "): ",
      paste(traits, collapse = "; ")
    ))
  }

  # Add gene annotation if available
  if (!is.null(gene_annotation)) {
    gene_name_clean <- gsub("\\(.*\\)", "", nearest_gene)
    gene_name_clean <- trimws(gsub("INTERGENIC:.*", "", gene_name_clean))
    gene_row <- gene_annotation[gene_annotation$gene_name == gene_name_clean, ]
    if (nrow(gene_row) > 0) {
      gene_row <- gene_row[1, ]
      if ("description" %in% names(gene_row) && !is.na(gene_row$description))
        sections <- c(sections, paste0("Gene description: ", gene_row$description))
      if ("reactome_pathways" %in% names(gene_row) && !is.na(gene_row$reactome_pathways))
        sections <- c(sections, paste0("Pathways: ", gene_row$reactome_pathways))
      if ("hpo_terms" %in% names(gene_row) && !is.na(gene_row$hpo_terms))
        sections <- c(sections, paste0("HPO phenotypes: ", gene_row$hpo_terms))
    }
  }

  paste(sections, collapse = "\n")
}

# Base system prompt for chat conversations
COLOC_CHAT_SYSTEM_PROMPT <- paste0(
  "You are a helpful genomics research assistant with expertise in genetic ",
  "colocalization analysis, GWAS interpretation, and kidney biology. ",
  "You are embedded in the Kidney Genomics Colocalization Atlas - a resource ",
  "that tests whether kidney trait GWAS signals share causal variants with ",
  "eQTL, pQTL, metabolite QTL, disease/phenotype GWAS, and kidney MRI traits. ",
  "PP.H4 > 0.8 indicates strong evidence of a shared causal variant.\n\n",
  "The user is looking at a specific genomic region. The context for this ",
  "region is provided below. Engage conversationally: answer questions about ",
  "the likely causal gene and mechanism, the eQTL/pQTL/metabolite evidence, ",
  "disease connections, and clinical relevance. Be evidence-based; distinguish ",
  "established biology from plausible hypotheses. Keep answers focused and ",
  "readable. Use hyphens instead of em dashes."
)

# Build a system prompt with embedded region context
build_chat_system_prompt <- function(context) {
  if (is.null(context) || nchar(context) == 0) return(COLOC_CHAT_SYSTEM_PROMPT)
  paste0(COLOC_CHAT_SYSTEM_PROMPT,
         "\n\n--- CURRENT REGION CONTEXT ---\n",
         context)
}

# Build the opening greeting prompt (sent as first "user" turn to get an
# initial assistant message that acknowledges the region)
build_greeting_prompt <- function(region_str, nearest_gene) {
  paste0(
    "Please greet the user briefly (1-2 sentences). Mention that you see ",
    "they are looking at region ", region_str,
    if (!is.null(nearest_gene) && !is.na(nearest_gene) && nchar(nearest_gene) > 0) {
      paste0(" near ", gsub("\\(.*\\)", "", nearest_gene))
    } else "",
    ", summarize in one sentence the groups of colocalizations you see, and ",
    "invite them to ask questions. Do not repeat the full context; keep it short."
  )
}
