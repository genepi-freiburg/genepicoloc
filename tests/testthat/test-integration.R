# Integration tests using packaged data
# tests/testthat/test-integration.R

library(testthat)

test_that("packaged gene data integration works", {
  # Skip if packaged data not available (e.g., during development)
  packaged_file <- system.file("extdata", "genes_chr.txt.gz", package = "genepicoloc")
  skip_if(packaged_file == "", "Packaged gene data not available")
  
  # Test setup with packaged data
  result <- setup_gene_annotation(verbose = FALSE)
  
  # Basic structure checks
  expect_type(result, "list")
  expect_named(result, c("gene_map", "metadata"))
  expect_equal(result$metadata$source, "packaged")
  
  # Check reasonable number of genes (should be ~19,000-20,000)
  expect_gt(nrow(result$gene_map), 15000)
  expect_lt(nrow(result$gene_map), 25000)
  
  # Test known gene positions
  # These are real coordinates that should be stable
  
  # Test ABCB11 (chr2:168,915,498-169,031,324)
  abcb11_result <- annotate_position("2", 169000000, result$gene_map)
  expect_equal(abcb11_result, "ABCB11")
  
  # Test LRP2 (chr2:169,127,109-169,362,534) 
  lrp2_result <- annotate_position("2", 169127111, result$gene_map)
  expect_equal(lrp2_result, "LRP2")
  
  # Test intergenic region
  intergenic_result <- annotate_position("1", 500000, result$gene_map)
  expect_true(grepl("INTERGENIC|ERROR", intergenic_result))
})

test_that("real-world genomic positions work correctly", {
  packaged_file <- system.file("extdata", "genes_chr.txt.gz", package = "genepicoloc")
  skip_if(packaged_file == "", "Packaged gene data not available")
  
  result <- setup_gene_annotation(verbose = FALSE)
  gene_map <- result$gene_map
  
  # Test various chromosome formats
  expect_type(annotate_position("1", 1000000, gene_map), "character")
  expect_type(annotate_position("X", 50000000, gene_map), "character")
  expect_type(annotate_position("22", 30000000, gene_map), "character")
  
  # Test edge cases
  expect_true(grepl("ERROR", annotate_position("25", 1000000, gene_map)))  # Invalid chr
  expect_true(grepl("ERROR", annotate_position("1", NA, gene_map)))         # Invalid pos
})

test_that("batch processing with real data works", {
  packaged_file <- system.file("extdata", "genes_chr.txt.gz", package = "genepicoloc")
  skip_if(packaged_file == "", "Packaged gene data not available")
  
  result <- setup_gene_annotation(verbose = FALSE)
  
  # Create annotation function
  annotation_func <- function(chr, pos) {
    annotate_position(chr, pos, result$gene_map)
  }
  
  # Create realistic test dataset
  test_gwas <- data.frame(
    Chr = c("1", "2", "2", "7", "X"),
    Pos = c(230000000, 169000000, 169127111, 140453136, 50000000),
    SNP_ID = paste0("rs", 1:5),
    P_value = c(1e-8, 5e-9, 2e-7, 3e-6, 8e-5),
    stringsAsFactors = FALSE
  )
  
  # Test batch annotation
  annotated <- batch_annotate_positions(
    data = test_gwas,
    annotation_function = annotation_func,
    verbose = FALSE
  )
  
  # Check structure
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), 5)
  expect_true("Gene_Annotation" %in% colnames(annotated))
  
  # Check that we got some gene annotations
  annotations <- annotated$Gene_Annotation
  n_genic <- sum(!grepl("^INTERGENIC:", annotations) & !grepl("^ERROR:", annotations))
  expect_gt(n_genic, 0)  # Should find at least some genes
  
  # Test specific known positions
  expect_equal(annotated$Gene_Annotation[2], "ABCB11")  # chr2:169000000
  expect_equal(annotated$Gene_Annotation[3], "LRP2")    # chr2:169127111
})

test_that("performance with larger datasets", {
  packaged_file <- system.file("extdata", "genes_chr.txt.gz", package = "genepicoloc")
  skip_if(packaged_file == "", "Packaged gene data not available")
  
  result <- setup_gene_annotation(verbose = FALSE)
  
  # Create annotation function
  annotation_func <- function(chr, pos) {
    annotate_position(chr, pos, result$gene_map)
  }
  
  # Create moderately large test dataset (100 positions)
  set.seed(123)  # For reproducible tests
  large_test <- data.frame(
    Chr = sample(c(1:22, "X"), 100, replace = TRUE),
    Pos = sample(1:250000000, 100),
    ID = paste0("variant_", 1:100),
    stringsAsFactors = FALSE
  )
  
  # Test that it completes in reasonable time
  start_time <- Sys.time()
  annotated_large <- batch_annotate_positions(
    data = large_test,
    annotation_function = annotation_func,
    verbose = FALSE
  )
  end_time <- Sys.time()
  
  # Should complete in under 10 seconds for 100 positions
  expect_lt(as.numeric(end_time - start_time), 10)
  
  # Check all positions were annotated
  expect_equal(nrow(annotated_large), 100)
  expect_true(all(!is.na(annotated_large$Gene_Annotation)))
})

test_that("chromosome format handling", {
  packaged_file <- system.file("extdata", "genes_chr.txt.gz", package = "genepicoloc")
  skip_if(packaged_file == "", "Packaged gene data not available")
  
  result <- setup_gene_annotation(verbose = FALSE)
  gene_map <- result$gene_map
  
  # Test different chromosome input formats
  expect_type(annotate_position(1, 1000000, gene_map), "character")      # Numeric
  expect_type(annotate_position("1", 1000000, gene_map), "character")    # Character
  expect_type(annotate_position("01", 1000000, gene_map), "character")   # Zero-padded
  
  # Test X chromosome
  expect_type(annotate_position("X", 50000000, gene_map), "character")
  expect_type(annotate_position("x", 50000000, gene_map), "character")   # Lowercase
  
  # Test invalid chromosomes
  expect_true(grepl("ERROR", annotate_position("chr1", 1000000, gene_map)))  # With chr prefix
  expect_true(grepl("ERROR", annotate_position("Y", 1000000, gene_map)))     # Y chromosome
  expect_true(grepl("ERROR", annotate_position("MT", 1000000, gene_map)))    # Mitochondrial
})

test_that("gene data metadata is accurate", {
  packaged_file <- system.file("extdata", "genes_chr.txt.gz", package = "genepicoloc")
  skip_if(packaged_file == "", "Packaged gene data not available")
  
  result <- setup_gene_annotation(verbose = FALSE)
  meta <- result$metadata
  
  # Check metadata structure
  expect_type(meta, "list")
  expect_true(all(c("source", "file_path", "total_genes", "chromosomes", "assembly") %in% names(meta)))
  
  # Check metadata values
  expect_equal(meta$source, "packaged")
  expect_equal(meta$assembly, "GRCh38")
  expect_equal(meta$total_genes, nrow(result$gene_map))
  
  # Check chromosomes are reasonable
  expected_chrs <- c(as.character(1:22), "X")
  expect_true(all(expected_chrs %in% meta$chromosomes))
  expect_equal(length(meta$chromosomes), 23)  # 1-22 plus X
})
