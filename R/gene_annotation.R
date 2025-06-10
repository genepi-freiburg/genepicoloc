#' Setup Gene Annotation from Processed Gene Data
#'
#' Load processed gene annotations and create functions for annotating genomic positions.
#' This function uses either the packaged gene data or a user-provided file.
#'
#' @description
#' This function loads gene annotation data and returns functions for annotating
#' genomic positions. By default, it uses the pre-processed gene data included
#' with the package, but users can provide their own updated gene files.
#'
#' @param gene_file Character. Path to processed gene file. If NULL, uses packaged data.
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{gene_map}{Data frame with processed gene annotations}
#'   \item{annotate_position}{Function to annotate genomic positions}
#'   \item{metadata}{Information about the gene data source}
#'
#' @examples
#' \dontrun{
#' # Method 1: Use packaged gene data (recommended for most users)
#' result <- setup_gene_annotation()
#' 
#' # Annotate single positions
#' gene1 <- result$annotate_position(chr = "2", pos = 169000000)  # ABCB11
#' gene2 <- result$annotate_position(chr = "2", pos = 169127111)  # LRP2
#' gene3 <- result$annotate_position(chr = "1", pos = 1000000)    # Intergenic
#' 
#' print(gene1)  # "ABCB11"
#' print(gene2)  # "LRP2"
#' print(gene3)  # "INTERGENIC: upstream=..."
#' 
#' # Method 2: Use updated gene data
#' # First download and process with download_gene_data()
#' result_updated <- setup_gene_annotation(gene_file = "my_updated_genes.txt.gz")
#' 
#' # Check gene data info
#' print(result$metadata)
#' cat("Total genes:", nrow(result$gene_map), "\n")
#' head(result$gene_map)
#' }
#'
#' @importFrom utils read.table
#' @author Generated from workflow by Matthias Wuttke
#' @export

setup_gene_annotation <- function(gene_file = NULL, verbose = TRUE) {
  
  if (verbose) cat("Setting up gene annotation...\n")
  
  # Determine gene file to use
  if (is.null(gene_file)) {
    # Try to use packaged data
    packaged_file <- system.file("extdata", "genes_chr.txt.gz", package = "genepicoloc")
    if (file.exists(packaged_file)) {
      gene_file <- packaged_file
      data_source <- "packaged"
      if (verbose) cat("Using packaged gene data from genepicoloc package\n")
    } else {
      stop("No gene file specified and packaged data not found.",
           "\nPlease either:",
           "\n1. Provide gene_file parameter with path to processed gene data, or",
           "\n2. Install package with bundled gene data, or", 
           "\n3. Use download_gene_data() to create a gene file first")
    }
  } else {
    data_source <- "user_provided"
    if (verbose) cat("Using user-provided gene file:", gene_file, "\n")
  }
  
  # Check if file exists
  if (!file.exists(gene_file)) {
    stop("Gene file not found: ", gene_file,
         "\nIf you need to create a gene file, use download_gene_data() first.")
  }
  
  # Read gene data (handle both .gz and regular files)
  if (verbose) cat("Reading gene data from:", gene_file, "\n")
  
  if (grepl("\\.gz$", gene_file)) {
    # Read compressed file
    e <- read.table(gzfile(gene_file), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else {
    # Read regular file  
    e <- read.table(gene_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }
  
  # Validate data structure
  expected_cols <- c("HGNC", "GeneStart", "GeneEnd", "Chromosome")
  if (!all(expected_cols %in% colnames(e))) {
    stop("Invalid gene file format. Expected columns: ", paste(expected_cols, collapse = ", "),
         "\nFound columns: ", paste(colnames(e), collapse = ", "),
         "\nPlease use download_gene_data() to create a properly formatted file.")
  }
  
  # Convert chromosome to character for consistency
  e$Chromosome <- as.character(e$Chromosome)
  
  if (verbose) {
    cat("Loaded", nrow(e), "protein-coding genes\n")
    cat("Chromosomes:", paste(sort(unique(e$Chromosome)), collapse = ", "), "\n")
  }
  
  # Create metadata
  metadata <- list(
    source = data_source,
    file_path = gene_file,
    total_genes = nrow(e),
    chromosomes = sort(unique(e$Chromosome)),
    assembly = "GRCh38",
    created = Sys.time()
  )
  
  # Create the annotation function
  annotate_position <- function(chr, pos) {
    # Annotate a genomic position with gene information
    # 
    # chr: Character or numeric. Chromosome (1-22, X)
    # pos: Numeric. Genomic position in base pairs
    # Returns: Character string with gene annotation
    
    chr <- as.character(chr)
    pos <- as.numeric(pos)
    
    # Validate inputs
    if (is.na(chr) || is.na(pos)) {
      return("ERROR: Invalid chromosome or position")
    }
    
    # Filter for the specific chromosome
    map_chr <- e[e$Chromosome == chr, ]
    
    if (nrow(map_chr) == 0) {
      return(paste0("ERROR: Unknown chromosome: ", chr))
    }
    
    # Check if position is within any gene (intronic/exonic)
    intronic <- map_chr[map_chr$GeneStart <= pos & map_chr$GeneEnd >= pos, "HGNC"]
    
    if (length(intronic) > 0) {
      if (length(intronic) > 1) {
        # Multiple overlapping genes
        return(paste(intronic, collapse = ", "))
      } else {
        return(as.character(intronic[1]))
      }
    } else {
      # Position is intergenic - find closest upstream and downstream genes
      upstream_diff <- map_chr$GeneStart - pos
      upstream_diff[upstream_diff < 0] <- Inf  # Genes that start before position
      
      downstream_diff <- pos - map_chr$GeneEnd  
      downstream_diff[downstream_diff < 0] <- Inf  # Genes that end after position
      
      # Find closest genes
      closest_upstream_idx <- which.min(upstream_diff)
      closest_downstream_idx <- which.min(downstream_diff)
      
      upstream_gene <- map_chr[closest_upstream_idx, "HGNC"]
      downstream_gene <- map_chr[closest_downstream_idx, "HGNC"]
      upstream_dist <- upstream_diff[closest_upstream_idx]
      downstream_dist <- downstream_diff[closest_downstream_idx]
      
      # Handle edge cases
      if (is.infinite(upstream_dist)) upstream_dist <- "NA"
      if (is.infinite(downstream_dist)) downstream_dist <- "NA"
      
      return(paste0("INTERGENIC: upstream=", upstream_gene, "(+", upstream_dist, "bp), ",
                    "downstream=", downstream_gene, "(-", downstream_dist, "bp)"))
    }
  }
  
  if (verbose) {
    cat("Gene annotation setup complete!\n")
    cat("Testing with example positions...\n")
    cat("Chr 2, pos 169000000:", annotate_position("2", 169000000), "\n")
    cat("Chr 2, pos 169127111:", annotate_position("2", 169127111), "\n")
    cat("Chr 1, pos 1000000:", annotate_position("1", 1000000), "\n")
  }
  
  # Return annotation tools and data
  return(list(
    gene_map = e,
    annotate_position = annotate_position,
    metadata = metadata
  ))
}

#' Batch Annotate Multiple Genomic Positions
#'
#' Efficiently annotate multiple genomic positions using gene annotation data.
#' This function adds gene annotations to a data frame containing genomic coordinates.
#'
#' @param data Data frame with columns for chromosome and position
#' @param chr_col Character. Name of chromosome column (default: "Chr")  
#' @param pos_col Character. Name of position column (default: "Pos")
#' @param annotation_function Function returned by setup_gene_annotation()$annotate_position
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return Original data frame with added "Gene_Annotation" column
#'
#' @examples
#' \dontrun{
#' # First set up gene annotation
#' result <- setup_gene_annotation()
#' 
#' # Example 1: Simple batch annotation
#' my_variants <- data.frame(
#'   Chr = c("1", "2", "7", "17"),
#'   Pos = c(230000000, 169000000, 140453136, 41276113),
#'   Variant_ID = c("var1", "var2", "var3", "var4")
#' )
#' 
#' annotated <- batch_annotate_positions(
#'   data = my_variants,
#'   annotation_function = result$annotate_position
#' )
#' print(annotated)
#' 
#' # Example 2: Custom column names
#' gwas_data <- data.frame(
#'   chromosome = c("2", "4", "X"),
#'   bp_position = c(169127111, 88000000, 60000000),
#'   rsid = c("rs123", "rs456", "rs789"),
#'   pvalue = c(1e-8, 5e-7, 2e-6)
#' )
#' 
#' gwas_annotated <- batch_annotate_positions(
#'   data = gwas_data,
#'   chr_col = "chromosome",
#'   pos_col = "bp_position", 
#'   annotation_function = result$annotate_position
#' )
#' print(gwas_annotated$Gene_Annotation)
#' 
#' # Example 3: Large dataset processing
#' # For datasets with >10,000 variants, consider data.table for speed
#' if (require(data.table)) {
#'   large_data <- data.table(
#'     Chr = sample(c(1:22, "X"), 1000, replace = TRUE),
#'     Pos = sample(1:250000000, 1000),
#'     ID = paste0("variant_", 1:1000)
#'   )
#'   
#'   # This will be faster for large datasets
#'   large_data[, Gene := result$annotate_position(Chr, Pos), by = 1:nrow(large_data)]
#' }
#' }
#'
#' @export

batch_annotate_positions <- function(data, chr_col = "Chr", pos_col = "Pos", 
                                     annotation_function, verbose = TRUE) {
  
  # Validate inputs
  if (missing(annotation_function)) {
    stop("Please provide annotation_function from setup_gene_annotation()$annotate_position")
  }
  
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  if (!chr_col %in% colnames(data)) {
    stop("Column '", chr_col, "' not found in data. Available columns: ", 
         paste(colnames(data), collapse = ", "))
  }
  
  if (!pos_col %in% colnames(data)) {
    stop("Column '", pos_col, "' not found in data. Available columns: ", 
         paste(colnames(data), collapse = ", "))
  }
  
  if (verbose) cat("Annotating", nrow(data), "positions...\n")
  
  # Add gene annotations
  data$Gene_Annotation <- apply(data, 1, function(row) {
    annotation_function(row[chr_col], as.numeric(row[pos_col]))
  })
  
  if (verbose) {
    cat("Annotation complete!\n")
    
    # Summary statistics
    annotations <- data$Gene_Annotation
    n_genic <- sum(!grepl("^INTERGENIC:", annotations) & !grepl("^ERROR:", annotations))
    n_intergenic <- sum(grepl("^INTERGENIC:", annotations))
    n_errors <- sum(grepl("^ERROR:", annotations))
    
    cat("Summary:\n")
    cat("- Genic positions:", n_genic, sprintf("(%.1f%%)", 100 * n_genic / nrow(data)), "\n")
    cat("- Intergenic positions:", n_intergenic, sprintf("(%.1f%%)", 100 * n_intergenic / nrow(data)), "\n")
    if (n_errors > 0) {
      cat("- Errors:", n_errors, sprintf("(%.1f%%)", 100 * n_errors / nrow(data)), "\n")
    }
  }
  
  return(data)
}#' Download and Process Protein-Coding Genes from Ensembl BioMart
#'
#' This function provides a complete workflow to download and process protein-coding genes
#' from Ensembl BioMart. Use this function to update the gene annotations with the latest
#' data from Ensembl.
#'
#' @description
#' The workflow consists of three main steps:
#' 1. Download gene data from Ensembl BioMart (manual step)
#' 2. Process and clean the data 
#' 3. Save processed gene annotations
#'
#' @details
#' Manual BioMart Download Steps:
#' 1. Go to: http://www.ensembl.org/biomart/martview/
#' 2. Select Database: Ensembl Genes
#' 3. Select Dataset: Human genes (GRCh38.p14)
#' 4. Go to Filters panel:
#'    - Gene type: Select "protein_coding"
#' 5. Go to Attributes panel:
#'    - Gene stable ID: Uncheck
#'    - Gene name (HGNC symbol): Check
#'    - Gene start (bp): Check  
#'    - Gene end (bp): Check
#'    - Chromosome/scaffold name: Check
#' 6. Click "Results" and download as TSV file
#'
#' Alternative direct BioMart URL (may need updating):
#' http://www.ensembl.org/biomart/martview/62ec57e3bef4d3400759e631a7778f3e?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.hgnc_symbol|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.chromosome_name&FILTERS=hsapiens_gene_ensembl.default.filters.biotype."protein_coding"&VISIBLEPANEL=filterpanel
#'
#' @param biomart_file Character. Path to the downloaded BioMart TSV file
#' @param output_file Character. Path for the processed gene file (default: "genes_chr.txt")
#' @param compress Logical. Whether to compress the output file with gzip (default: TRUE)
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return Character. Path to the created gene annotation file
#'
#' @examples
#' \dontrun{
#' # Step 1: Download data from BioMart manually
#' # Follow the steps in function documentation
#' # Save as "mart_export.txt"
#' 
#' # Step 2: Process the downloaded data
#' gene_file <- download_gene_data(
#'   biomart_file = "mart_export.txt",
#'   output_file = "my_genes_chr.txt",
#'   compress = TRUE
#' )
#' 
#' # The processed file is now ready for use with setup_gene_annotation()
#' print(paste("Gene data saved to:", gene_file))
#' 
#' # Example: Update package data
#' # This creates a file that could replace the packaged data
#' updated_genes <- download_gene_data(
#'   biomart_file = "fresh_biomart_download.txt",
#'   output_file = "genes_chr_updated.txt.gz"
#' )
#' }
#'
#' @importFrom utils read.table write.table
#' @author Generated from workflow by Matthias Wuttke
#' @export

download_gene_data <- function(biomart_file, 
                               output_file = "genes_chr.txt", 
                               compress = TRUE,
                               verbose = TRUE) {
  
  if (verbose) cat("Processing gene data from Ensembl BioMart...\n")
  
  # Check if input file exists
  if (!file.exists(biomart_file)) {
    stop(paste("BioMart file not found:", biomart_file, 
               "\nPlease download from Ensembl BioMart first (see function documentation)",
               "\nURL: http://www.ensembl.org/biomart/martview/"))
  }
  
  # Read BioMart data
  if (verbose) cat("Reading BioMart data from:", biomart_file, "\n")
  d <- read.table(biomart_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  if (verbose) cat("Original dataset contains", nrow(d), "entries\n")
  
  # Validate expected columns
  expected_cols <- c("HGNC.symbol", "Gene.start..bp.", "Gene.end..bp.", "Chromosome.scaffold.name")
  if (!all(expected_cols %in% colnames(d))) {
    stop(paste("Expected columns not found in BioMart file:",
               "\nExpected:", paste(expected_cols, collapse = ", "),
               "\nFound:", paste(colnames(d), collapse = ", "),
               "\nPlease check your BioMart download settings."))
  }
  
  # Filter for main chromosomes and non-empty gene symbols
  if (verbose) cat("Filtering for chromosomes 1-22, X and non-empty HGNC symbols...\n")
  e <- subset(d, d$Chromosome.scaffold.name %in% c(1:22, "X"))
  e <- subset(e, e$HGNC.symbol != "")
  
  if (verbose) cat("Filtered dataset contains", nrow(e), "protein-coding genes\n")
  
  # Standardize column names
  colnames(e) <- c("HGNC", "GeneStart", "GeneEnd", "Chromosome")
  
  # Convert chromosome to character for consistency
  e$Chromosome <- as.character(e$Chromosome)
  
  # Sort by chromosome and position for efficient access
  e <- e[order(as.numeric(ifelse(e$Chromosome == "X", 23, e$Chromosome)), e$GeneStart), ]
  
  # Add compression suffix if requested
  if (compress && !grepl("\\.gz$", output_file)) {
    output_file <- paste0(output_file, ".gz")
  }
  
  # Write processed file
  if (verbose) cat("Writing processed data to:", output_file, "\n")
  
  if (compress && grepl("\\.gz$", output_file)) {
    # Write compressed file
    gz_con <- gzfile(output_file, "w")
    write.table(e, gz_con, row.names = FALSE, col.names = TRUE, sep = "\t")
    close(gz_con)
  } else {
    # Write regular file
    write.table(e, output_file, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
  
  if (verbose) {
    cat("Gene data processing complete!\n")
    cat("Summary:\n")
    cat("- Total genes:", nrow(e), "\n")
    cat("- Chromosomes:", paste(sort(unique(e$Chromosome)), collapse = ", "), "\n")
    cat("- Output file:", output_file, "\n")
    cat("- File size:", round(file.info(output_file)$size / 1024^2, 2), "MB\n")
  }
  
  return(output_file)
}

#' Get Current Ensembl Release Information
#'
#' Fetch information about the current Ensembl release to help users
#' determine if they need to update their gene annotations.
#'
#' @param verbose Logical. Print release information (default: TRUE)
#' @return List with release information, or NULL if unable to fetch
#'
#' @examples
#' \dontrun{
#' # Check current Ensembl release
#' release_info <- get_ensembl_release()
#' print(release_info$release)
#' }
#'
#' @export

get_ensembl_release <- function(verbose = TRUE) {
  
  if (verbose) cat("Checking current Ensembl release...\n")
  
  tryCatch({
    # Try to get release info from Ensembl REST API
    if (requireNamespace("httr", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE)) {
      
      response <- httr::GET("https://rest.ensembl.org/info/data/?content-type=application/json")
      
      if (httr::status_code(response) == 200) {
        content <- httr::content(response, as = "text", encoding = "UTF-8")
        release_data <- jsonlite::fromJSON(content)
        
        if (verbose) {
          cat("Current Ensembl release:", release_data$releases[[1]], "\n")
          cat("Assembly:", release_data$species$homo_sapiens$assembly, "\n")
        }
        
        return(list(
          release = release_data$releases[[1]],
          assembly = release_data$species$homo_sapiens$assembly,
          url = "http://www.ensembl.org/biomart/martview/"
        ))
      }
    }
    
    # Fallback message
    if (verbose) {
      cat("Unable to fetch current release info automatically.\n")
      cat("Please check manually at: http://www.ensembl.org/biomart/martview/\n")
    }
    
    return(list(
      release = "Unknown",
      assembly = "GRCh38",
      url = "http://www.ensembl.org/biomart/martview/"
    ))
    
  }, error = function(e) {
    if (verbose) {
      cat("Error checking Ensembl release:", e$message, "\n")
      cat("Please check manually at: http://www.ensembl.org/biomart/martview/\n")
    }
    return(NULL)
  })
}