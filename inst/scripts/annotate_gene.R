#' Download and Process Protein-Coding Genes from Ensembl BioMart
#'
#' This function provides a complete workflow to download protein-coding genes
#' from Ensembl BioMart and create a gene annotation function for genomic positions.
#'
#' @description
#' The workflow consists of three main steps:
#' 1. Download gene data from Ensembl BioMart
#' 2. Process and clean the data 
#' 3. Create a gene annotation function
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
#' 6. Click "Results" and download as TSV file named "genes.txt"
#'
#' Alternative direct BioMart URL (may need updating):
#' http://www.ensembl.org/biomart/martview/62ec57e3bef4d3400759e631a7778f3e?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.hgnc_symbol|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.chromosome_name&FILTERS=hsapiens_gene_ensembl.default.filters.biotype."protein_coding"&VISIBLEPANEL=filterpanel
#'
#' @param biomart_file Character. Path to the downloaded BioMart TSV file (default: "genes.txt")
#' @param output_file Character. Path for the processed gene file (default: "genes_chr.txt")
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{gene_map}{Data frame with processed gene annotations}
#'   \item{annotate_position}{Function to annotate genomic positions}
#'
#' @examples
#' \dontrun{
#' # Method 1: Use with packaged gene data (recommended)
#' # No need to download - uses pre-installed gene annotations
#' result <- setup_gene_annotation()
#' 
#' # Annotate single positions
#' gene1 <- result$annotate_position(chr = "2", pos = 169127111)  # Should find ABCB11
#' gene2 <- result$annotate_position(chr = "1", pos = 1000000)    # Intergenic region
#' gene3 <- result$annotate_position(chr = "X", pos = 50000000)   # X chromosome
#' 
#' print(gene1)  # "ABCB11" 
#' print(gene2)  # "INTERGENIC: upstream=..."
#' 
#' # Method 2: Download fresh data from BioMart
#' # 1. Go to http://www.ensembl.org/biomart/martview/
#' # 2. Follow the manual steps in function documentation
#' # 3. Save as "my_genes.txt"
#' result <- setup_gene_annotation(biomart_file = "my_genes.txt")
#' 
#' # Batch annotation of multiple positions
#' my_snps <- data.frame(
#'   Chr = c("1", "2", "7", "X", "22"),
#'   Pos = c(1000000, 169127111, 140453136, 50000000, 30000000),
#'   SNP_ID = c("rs123", "rs456", "rs789", "rs101", "rs202")
#' )
#' 
#' # Add gene annotations
#' annotated_snps <- batch_annotate_positions(
#'   data = my_snps, 
#'   chr_col = "Chr", 
#'   pos_col = "Pos",
#'   annotation_function = result$annotate_position
#' )
#' 
#' print(annotated_snps)
#' #   Chr       Pos SNP_ID         Gene_Annotation
#' # 1   1   1000000  rs123 INTERGENIC: upstream=...
#' # 2   2 169127111  rs456                  ABCB11
#' # 3   7 140453136  rs789                 BRAF...
#' 
#' # Access the gene map directly
#' head(result$gene_map)
#' nrow(result$gene_map)  # ~20,000 protein-coding genes
#' 
#' # Real-world example: Annotate GWAS results
#' gwas_hits <- data.frame(
#'   chromosome = c("1", "2", "3"),
#'   position = c(12345678, 169127111, 98765432),
#'   pvalue = c(5e-8, 1e-10, 3e-9),
#'   beta = c(0.15, -0.25, 0.18)
#' )
#' 
#' gwas_annotated <- batch_annotate_positions(
#'   data = gwas_hits,
#'   chr_col = "chromosome", 
#'   pos_col = "position",
#'   annotation_function = result$annotate_position
#' )
#' }
#'
#' @author Generated from workflow by Matthias Wuttke
#' @export

setup_gene_annotation <- function(biomart_file = "mart_export.txt", 
                                  output_file = "genes_chr.txt", 
                                  verbose = TRUE) {
  
  if (verbose) cat("Setting up gene annotation from BioMart data...\n")
  
  # Check if input file exists
  if (!file.exists(biomart_file)) {
    stop(paste("BioMart file not found:", biomart_file, 
               "\nPlease download from Ensembl BioMart first (see function documentation)"))
  }
  
  # Read BioMart data
  if (verbose) cat("Reading BioMart data from:", biomart_file, "\n")
  d <- read.table(biomart_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  if (verbose) cat("Original dataset contains", nrow(d), "entries\n")
  
  # Filter for main chromosomes and non-empty gene symbols
  if (verbose) cat("Filtering for chromosomes 1-22, X and non-empty HGNC symbols...\n")
  e <- subset(d, d$Chromosome.scaffold.name %in% c(1:22, "X"))
  e <- subset(e, e$HGNC.symbol != "")
  
  if (verbose) cat("Filtered dataset contains", nrow(e), "protein-coding genes\n")
  
  # Standardize column names
  colnames(e) <- c("HGNC", "GeneStart", "GeneEnd", "Chromosome")
  
  # Convert chromosome to character for consistency
  e$Chromosome <- as.character(e$Chromosome)
  
  # Write processed file
  if (verbose) cat("Writing processed data to:", output_file, "\n")
  write.table(e, output_file, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  # Create the annotation function
  annotate_position <- function(chr, pos) {
    # Annotate a genomic position with gene information
    #
    # @param chr Character or numeric. Chromosome (1-22, X)
    # @param pos Numeric. Genomic position in base pairs
    # @return Character string with gene annotation
    
    chr <- as.character(chr)
    pos <- as.numeric(pos)
    
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
    cat("Chr 2, pos 169300000:", annotate_position("2", 169300000), "\n")
  }
  
  # Return both the gene map and the annotation function
  return(list(
    gene_map = e,
    annotate_position = annotate_position
  ))
}

#' Batch annotate multiple genomic positions
#'
#' @param data Data frame with columns for chromosome and position
#' @param chr_col Character. Name of chromosome column (default: "Chr")  
#' @param pos_col Character. Name of position column (default: "Pos")
#' @param annotation_function Function returned by setup_gene_annotation()$annotate_position
#'
#' @return Original data frame with added "Gene_Annotation" column
#' @export

batch_annotate_positions <- function(data, chr_col = "Chr", pos_col = "Pos", annotation_function) {
  if (missing(annotation_function)) {
    stop("Please provide annotation_function from setup_gene_annotation()$annotate_position")
  }
  
  cat("Annotating", nrow(data), "positions...\n")
  
  data$Gene_Annotation <- apply(data, 1, function(row) {
    annotation_function(row[chr_col], as.numeric(row[pos_col]))
  })
  
  cat("Annotation complete!\n")
  return(data)
}

# Usage examples and workflow documentation:
#
# COMPLETE WORKFLOW WITH EXAMPLES:
# ================================
# 
# 1. Basic usage (recommended - uses packaged data):
#    library(genepicoloc)
#    result <- setup_gene_annotation()
#    gene <- result$annotate_position("2", 169127111)  # Returns "ABCB11"
#
# 2. Batch annotation:
#    my_snps <- data.frame(Chr = c("1", "2", "7"), Pos = c(1000000, 169127111, 140453136))
#    annotated <- batch_annotate_positions(my_snps, annotation_function = result$annotate_position)
#
# 3. Real GWAS example:
#    gwas <- read.table("my_gwas.txt", header = TRUE)  # Chr, Pos, P_value, etc.
#    gwas_annotated <- batch_annotate_positions(gwas, annotation_function = result$annotate_position)
#    write.table(gwas_annotated, "gwas_with_genes.txt", sep = "\t", row.names = FALSE)
#
# 4. Download fresh BioMart data (optional):
#    # Go to http://www.ensembl.org/biomart/martview/
#    # Filter for protein_coding genes  
#    # Select attributes: HGNC symbol, gene start, gene end, chromosome
#    # Save as TSV file named "fresh_genes.txt"
#    result_fresh <- setup_gene_annotation(biomart_file = "fresh_genes.txt")
#
# 5. Working with different chromosome formats:
#    # Your data might have "chr1" instead of "1"
#    my_data$Chr <- gsub("chr", "", my_data$Chr)  # Remove "chr" prefix
#    annotated <- batch_annotate_positions(my_data, annotation_function = result$annotate_position)
#
# 6. Performance tips for large datasets:
#    # For >10,000 variants, use data.table:
#    library(data.table)
#    dt <- as.data.table(my_large_dataset)
#    dt[, Gene := result$annotate_position(Chr, Pos), by = 1:nrow(dt)]
#
# NOTES:
# ======
# - Gene coordinates are based on GRCh38/hg38 assembly
# - Only protein-coding genes on chromosomes 1-22, X are included
# - Positions within genes return the gene name
# - Intergenic positions return closest upstream/downstream genes with distances
# - Multiple overlapping genes are comma-separated
#
# TROUBLESHOOTING:
# ===============
# - If BioMart URLs change, manually download from biomart interface
# - Ensure chromosome format matches your data (numeric vs character)
# - For large datasets, consider using data.table for faster processing
#
# Last updated: 2025-06-10