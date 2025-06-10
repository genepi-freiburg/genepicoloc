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
#'   \item{metadata}{Information about the gene data source}
#'
#' @examples
#' \dontrun{
#' # Method 1: Use packaged gene data (recommended for most users)
#' result <- setup_gene_annotation()
#' 
#' # Annotate single positions
#' gene1 <- annotate_position(chr = "2", pos = 169000000, gene_map = result$gene_map)  # ABCB11
#' gene2 <- annotate_position(chr = "2", pos = 169127111, gene_map = result$gene_map)  # LRP2
#' gene3 <- annotate_position(chr = "1", pos = 1000000, gene_map = result$gene_map)    # Intergenic
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
  
  # Return annotation tools and data
  return(list(
    gene_map = e,
    metadata = metadata
  ))
}

#' Annotate a genomic position with gene information
#'
#' @param chr Character or numeric. Chromosome (1-22, X)
#' @param pos Numeric. Genomic position in base pairs
#' @param gene_map Data frame. Gene map from setup_gene_annotation()
#' @param n_nearest Integer. Number of nearest genes to return for intergenic positions (default: 2)
#' @param output_format Character. Output format: "default", "simple", or "detailed" (default: "default")
#' @return Character string with gene annotation
#' 
#' @details
#' For genic positions (position within a gene), returns the gene name(s).
#' For intergenic positions, behavior depends on n_nearest:
#' - n_nearest = 1: Returns closest gene only
#' - n_nearest = 2: Returns upstream and downstream genes (default behavior)  
#' - n_nearest > 2: Returns N closest genes with distances
#' 
#' Output formats:
#' - "default": Standard format with distances (e.g., "INTERGENIC: upstream=GENE1(+1000bp), downstream=GENE2(-500bp)")
#' - "simple": Gene names only, comma-separated (e.g., "GENE1, GENE2, GENE3")
#' - "detailed": Genes with distances, comma-separated (e.g., "GENE1(+1000bp), GENE2(-500bp), GENE3(+2000bp)")
#'
#' @examples
#' \dontrun{
#' # Setup gene annotation first
#' result <- setup_gene_annotation()
#' 
#' # Standard usage
#' annotate_position("1", 1000000, gene_map = result$gene_map)
#' 
#' # Or in one line
#' annotate_position("1", 1000000, gene_map = setup_gene_annotation()$gene_map)
#' 
#' # Get 5 nearest genes in simple format
#' annotate_position("1", 1000000, gene_map = result$gene_map, 
#'                   n_nearest = 5, output_format = "simple")
#' 
#' # Get 3 nearest genes with detailed distances
#' annotate_position("1", 1000000, gene_map = result$gene_map, 
#'                   n_nearest = 3, output_format = "detailed")
#' 
#' # Get 10 nearest genes (default format)
#' annotate_position("1", 1000000, gene_map = result$gene_map, n_nearest = 10)
#' }
#' 
#' @export
annotate_position <- function(chr, pos, gene_map, n_nearest = 2, output_format = "default") {
  # Validate gene_map parameter
  if (missing(gene_map) || !is.data.frame(gene_map)) {
    stop("gene_map must be a data frame from setup_gene_annotation()$gene_map")
  }
  
  # Validate parameters
  if (!is.numeric(n_nearest) || n_nearest < 1) {
    stop("n_nearest must be a positive integer")
  }
  
  if (!output_format %in% c("default", "simple", "detailed")) {
    stop("output_format must be 'default', 'simple', or 'detailed'")
  }
  
  chr <- as.character(chr)
  pos <- as.numeric(pos)
  n_nearest <- as.integer(n_nearest)
  
  # Validate inputs
  if (is.na(chr) || is.na(pos)) {
    return("ERROR: Invalid chromosome or position")
  }
  
  # Filter for the specific chromosome
  map_chr <- gene_map[gene_map$Chromosome == chr, ]
  
  if (nrow(map_chr) == 0) {
    return(paste0("ERROR: Unknown chromosome: ", chr))
  }
  
  # Check if position is within any gene (intronic/exonic)
  intronic <- map_chr[map_chr$GeneStart <= pos & map_chr$GeneEnd >= pos, "HGNC"]
  
  if (length(intronic) > 0) {
    # Position is within gene(s)
    if (output_format == "default" && n_nearest == 2) {
      # Default behavior for backward compatibility
      if (length(intronic) > 1) {
        return(paste(intronic, collapse = ", "))
      } else {
        return(as.character(intronic[1]))
      }
    } else {
      # For other cases, still find nearest genes but mark the genic ones
      return(find_nearest_genes(map_chr, pos, n_nearest, output_format, genic_genes = intronic))
    }
  } else {
    # Position is intergenic - find nearest genes
    return(find_nearest_genes(map_chr, pos, n_nearest, output_format))
  }
}

#' Helper function to find nearest genes for intergenic positions
#'
#' @param map_chr Data frame. Gene map filtered for specific chromosome
#' @param pos Numeric. Genomic position
#' @param n_nearest Integer. Number of nearest genes to find
#' @param output_format Character. Output format
#' @param genic_genes Character vector. Genes that contain the position (optional)
#' @return Character string with nearest genes
find_nearest_genes <- function(map_chr, pos, n_nearest, output_format, genic_genes = NULL) {
  
  # Calculate distances to all genes
  gene_distances <- data.frame(
    gene = map_chr$HGNC,
    start = map_chr$GeneStart,
    end = map_chr$GeneEnd,
    stringsAsFactors = FALSE
  )
  
  # Calculate minimum distance to each gene
  gene_distances$distance <- apply(gene_distances, 1, function(row) {
    start <- as.numeric(row["start"])
    end <- as.numeric(row["end"])
    
    if (pos < start) {
      # Position is upstream of gene
      return(start - pos)
    } else if (pos > end) {
      # Position is downstream of gene  
      return(pos - end)
    } else {
      # Position is within gene - distance is 0
      return(0)
    }
  })
  
  # Add direction indicator
  gene_distances$direction <- ifelse(
    pos < gene_distances$start, "upstream",
    ifelse(pos > gene_distances$end, "downstream", "within")
  )
  
  # Sort by distance and get top N
  gene_distances <- gene_distances[order(gene_distances$distance), ]
  nearest_genes <- head(gene_distances, n_nearest)
  
  # Handle special case for default format with n_nearest = 2 and intergenic position
  if (output_format == "default" && n_nearest == 2 && is.null(genic_genes)) {
    return(format_default_output(nearest_genes, pos))
  }
  
  # Format output based on requested format
  if (output_format == "simple") {
    return(paste(nearest_genes$gene, collapse = ", "))
  } else if (output_format == "detailed") {
    gene_strings <- mapply(function(gene, dist, direction) {
      if (direction == "within") {
        paste0(gene, "(within)")
      } else {
        sign <- if (direction == "upstream") "+" else "-"
        paste0(gene, "(", sign, dist, "bp)")
      }
    }, nearest_genes$gene, nearest_genes$distance, nearest_genes$direction)
    return(paste(gene_strings, collapse = ", "))
  } else {
    # Default format for n_nearest != 2 or genic positions
    gene_strings <- mapply(function(gene, dist, direction) {
      if (direction == "within") {
        paste0(gene, "(within)")
      } else {
        sign <- if (direction == "upstream") "+" else "-"
        paste0(gene, "(", sign, dist, "bp)")
      }
    }, nearest_genes$gene, nearest_genes$distance, nearest_genes$direction)
    
    if (is.null(genic_genes)) {
      return(paste0("INTERGENIC: ", paste(gene_strings, collapse = ", ")))
    } else {
      return(paste(gene_strings, collapse = ", "))
    }
  }
}

#' Helper function to format default output (backward compatibility)
#'
#' @param nearest_genes Data frame with nearest genes
#' @param pos Numeric. Position
#' @return Character string in default format
format_default_output <- function(nearest_genes, pos) {
  if (nrow(nearest_genes) == 0) {
    return("INTERGENIC: no nearby genes")
  }
  
  # Try to find one upstream and one downstream gene for classic format
  upstream_genes <- nearest_genes[nearest_genes$direction == "upstream", ]
  downstream_genes <- nearest_genes[nearest_genes$direction == "downstream", ]
  
  if (nrow(upstream_genes) > 0 && nrow(downstream_genes) > 0) {
    # Classic format: one upstream, one downstream
    upstream_gene <- upstream_genes[1, ]
    downstream_gene <- downstream_genes[1, ]
    
    upstream_dist <- upstream_gene$distance
    downstream_dist <- downstream_gene$distance
    
    return(paste0("INTERGENIC: upstream=", upstream_gene$gene, "(+", upstream_dist, "bp), ",
                  "downstream=", downstream_gene$gene, "(-", downstream_dist, "bp)"))
  } else {
    # Fallback: use two closest genes regardless of direction
    if (nrow(nearest_genes) >= 2) {
      gene1 <- nearest_genes[1, ]
      gene2 <- nearest_genes[2, ]
      
      sign1 <- if (gene1$direction == "upstream") "+" else "-"
      sign2 <- if (gene2$direction == "upstream") "+" else "-"
      
      return(paste0("INTERGENIC: ", gene1$gene, "(", sign1, gene1$distance, "bp), ",
                    gene2$gene, "(", sign2, gene2$distance, "bp)"))
    } else {
      # Only one gene found
      gene1 <- nearest_genes[1, ]
      sign1 <- if (gene1$direction == "upstream") "+" else "-"
      return(paste0("INTERGENIC: ", gene1$gene, "(", sign1, gene1$distance, "bp)"))
    }
  }
}
#' Batch Annotate Multiple Genomic Positions
#'
#' Efficiently annotate multiple genomic positions using gene annotation data.
#' This function adds gene annotations to a data frame containing genomic coordinates.
#'
#' @param data Data frame with columns for chromosome and position
#' @param chr_col Character. Name of chromosome column (default: "Chr")  
#' @param pos_col Character. Name of position column (default: "Pos")
#' @param annotation_function Function to annotate positions. Create with:
#'   annotation_function <- function(chr, pos) { annotate_position(chr, pos, gene_map) }
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return Original data frame with added "Gene_Annotation" column
#'
#' @examples
#' \dontrun{
#' # First set up gene annotation
#' result <- setup_gene_annotation()
#' 
#' # Create annotation function
#' annotation_func <- function(chr, pos) {
#'   annotate_position(chr, pos, result$gene_map)
#' }
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
#'   annotation_function = annotation_func
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
#'   annotation_function = annotation_func
#' )
#' print(gwas_annotated$Gene_Annotation)
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
}

#' Download and Process Protein-Coding Genes from Ensembl BioMart
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
  if (!all(expected_cols == colnames(d))) {
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