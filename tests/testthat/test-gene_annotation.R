# Test file for gene annotation functions
# tests/testthat/test-gene_annotation.R

library(testthat)

# Helper function to create test gene data
create_test_gene_data <- function() {
  data.frame(
    HGNC = c("GENE1", "GENE2", "GENE3", "GENE4"),
    GeneStart = c(1000, 5000, 10000, 50000),
    GeneEnd = c(2000, 6000, 12000, 52000),
    Chromosome = c("1", "1", "2", "X"),
    stringsAsFactors = FALSE
  )
}

# Create temporary test file
create_test_gene_file <- function() {
  test_data <- create_test_gene_data()
  temp_file <- tempfile(fileext = ".txt")
  write.table(test_data, temp_file, sep = "\t", row.names = FALSE, col.names = TRUE)
  return(temp_file)
}

# Create compressed test file
create_test_gene_file_gz <- function() {
  test_data <- create_test_gene_data()
  temp_file <- tempfile(fileext = ".txt.gz")
  gz_con <- gzfile(temp_file, "w")
  write.table(test_data, gz_con, sep = "\t", row.names = FALSE, col.names = TRUE)
  close(gz_con)
  return(temp_file)
}

# Test setup_gene_annotation function
test_that("setup_gene_annotation works with user-provided file", {
  # Create test file
  test_file <- create_test_gene_file()
  
  # Test setup
  result <- setup_gene_annotation(gene_file = test_file, verbose = FALSE)
  
  # Check structure
  expect_type(result, "list")
  expect_named(result, c("gene_map", "metadata"))
  
  # Check gene_map
  expect_s3_class(result$gene_map, "data.frame")
  expect_equal(nrow(result$gene_map), 4)
  expect_named(result$gene_map, c("HGNC", "GeneStart", "GeneEnd", "Chromosome"))
  
  # Check metadata
  expect_type(result$metadata, "list")
  expect_equal(result$metadata$source, "user_provided")
  expect_equal(result$metadata$total_genes, 4)
  
  # Clean up
  file.remove(test_file)
})

test_that("setup_gene_annotation works with compressed files", {
  # Create compressed test file
  test_file_gz <- create_test_gene_file_gz()
  
  # Test setup
  result <- setup_gene_annotation(gene_file = test_file_gz, verbose = FALSE)
  
  # Check that it loads correctly
  expect_equal(nrow(result$gene_map), 4)
  expect_equal(result$metadata$source, "user_provided")
  
  # Clean up
  file.remove(test_file_gz)
})

test_that("setup_gene_annotation handles missing files", {
  expect_error(
    setup_gene_annotation(gene_file = "nonexistent_file.txt", verbose = FALSE),
    "Gene file not found"
  )
})

test_that("setup_gene_annotation validates file format", {
  # Create file with wrong columns
  bad_data <- data.frame(
    WrongCol1 = c("A", "B"),
    WrongCol2 = c(1, 2)
  )
  temp_file <- tempfile(fileext = ".txt")
  write.table(bad_data, temp_file, sep = "\t", row.names = FALSE)
  
  expect_error(
    setup_gene_annotation(gene_file = temp_file, verbose = FALSE),
    "Invalid gene file format"
  )
  
  file.remove(temp_file)
})

# Test annotate_position function (standalone)
test_that("annotate_position works correctly", {
  test_file <- create_test_gene_file()
  result <- setup_gene_annotation(gene_file = test_file, verbose = FALSE)
  gene_map <- result$gene_map
  
  # Test genic positions
  expect_equal(annotate_position("1", 1500, gene_map), "GENE1")  # Within GENE1
  expect_equal(annotate_position("2", 11000, gene_map), "GENE3") # Within GENE3
  expect_equal(annotate_position("X", 51000, gene_map), "GENE4") # Within GENE4
  
  # Test intergenic position
  intergenic_result <- annotate_position("1", 3000, gene_map)  # Between GENE1 and GENE2
  expect_true(grepl("INTERGENIC", intergenic_result))
  expect_true(grepl("GENE1", intergenic_result))
  expect_true(grepl("GENE2", intergenic_result))
  
  # Test unknown chromosome
  expect_true(grepl("ERROR: Unknown chromosome", annotate_position("99", 1000, gene_map)))
  
  # Test invalid inputs
  expect_true(grepl("ERROR: Invalid", annotate_position(NA, 1000, gene_map)))
  expect_true(grepl("ERROR: Invalid", annotate_position("1", NA, gene_map)))
  
  # Test missing gene_map parameter
  expect_error(annotate_position("1", 1000), "gene_map must be a data frame")
  
  file.remove(test_file)
})

test_that("annotate_position handles overlapping genes", {
  # Create data with overlapping genes
  overlap_data <- data.frame(
    HGNC = c("GENE_A", "GENE_B"),
    GeneStart = c(1000, 1500),
    GeneEnd = c(2000, 2500),
    Chromosome = c("1", "1"),
    stringsAsFactors = FALSE
  )
  temp_file <- tempfile(fileext = ".txt")
  write.table(overlap_data, temp_file, sep = "\t", row.names = FALSE)
  
  result <- setup_gene_annotation(gene_file = temp_file, verbose = FALSE)
  
  # Test position in overlapping region
  overlap_result <- annotate_position("1", 1750, result$gene_map)
  expect_true(grepl("GENE_A.*GENE_B|GENE_B.*GENE_A", overlap_result))
  
  file.remove(temp_file)
})

# Test batch_annotate_positions function
test_that("batch_annotate_positions works correctly", {
  test_file <- create_test_gene_file()
  result <- setup_gene_annotation(gene_file = test_file, verbose = FALSE)
  
  # Create annotation function
  annotation_func <- function(chr, pos) {
    annotate_position(chr, pos, result$gene_map)
  }
  
  # Create test data
  test_positions <- data.frame(
    Chr = c("1", "2", "X", "1"),
    Pos = c(1500, 11000, 51000, 3000),
    ID = c("var1", "var2", "var3", "var4")
  )
  
  # Test batch annotation
  annotated <- batch_annotate_positions(
    data = test_positions,
    annotation_function = annotation_func,
    verbose = FALSE
  )
  
  # Check structure
  expect_s3_class(annotated, "data.frame")
  expect_true("Gene_Annotation" %in% colnames(annotated))
  expect_equal(nrow(annotated), 4)
  
  # Check specific annotations
  expect_equal(annotated$Gene_Annotation[1], "GENE1")
  expect_equal(annotated$Gene_Annotation[2], "GENE3")
  expect_equal(annotated$Gene_Annotation[3], "GENE4")
  expect_true(grepl("INTERGENIC", annotated$Gene_Annotation[4]))
  
  file.remove(test_file)
})

test_that("batch_annotate_positions handles custom column names", {
  test_file <- create_test_gene_file()
  result <- setup_gene_annotation(gene_file = test_file, verbose = FALSE)
  
  # Create annotation function
  annotation_func <- function(chr, pos) {
    annotate_position(chr, pos, result$gene_map)
  }
  
  # Create test data with custom column names
  test_positions <- data.frame(
    chromosome = c("1", "2"),
    bp_position = c(1500, 11000),
    variant_id = c("rs1", "rs2")
  )
  
  # Test with custom column names
  annotated <- batch_annotate_positions(
    data = test_positions,
    chr_col = "chromosome",
    pos_col = "bp_position",
    annotation_function = annotation_func,
    verbose = FALSE
  )
  
  expect_equal(annotated$Gene_Annotation[1], "GENE1")
  expect_equal(annotated$Gene_Annotation[2], "GENE3")
  
  file.remove(test_file)
})

test_that("batch_annotate_positions validates inputs", {
  test_file <- create_test_gene_file()
  result <- setup_gene_annotation(gene_file = test_file, verbose = FALSE)
  
  test_data <- data.frame(Chr = "1", Pos = 1000)
  
  # Test missing annotation_function parameter
  expect_error(
    batch_annotate_positions(data = test_data, verbose = FALSE),
    "Please provide annotation_function from setup_gene_annotation"
  )
  
  # Create annotation function for other tests
  annotation_func <- function(chr, pos) {
    annotate_position(chr, pos, result$gene_map)
  }
  
  # Test wrong data type
  expect_error(
    batch_annotate_positions(
      data = "not_a_dataframe",
      annotation_function = annotation_func,
      verbose = FALSE
    ),
    "data must be a data frame"
  )
  
  # Test missing columns
  expect_error(
    batch_annotate_positions(
      data = test_data,
      chr_col = "NonexistentCol",
      annotation_function = annotation_func,
      verbose = FALSE
    ),
    "Column 'NonexistentCol' not found"
  )
  
  file.remove(test_file)
})

# Test download_gene_data function
test_that("download_gene_data validates input file", {
  expect_error(
    download_gene_data("nonexistent_biomart_file.txt", verbose = FALSE),
    "BioMart file not found"
  )
})

test_that("download_gene_data validates column structure", {
  # Create file with wrong BioMart columns
  wrong_biomart <- data.frame(
    WrongCol1 = c("Gene1", "Gene2"),
    WrongCol2 = c(1000, 2000)
  )
  temp_file <- tempfile(fileext = ".txt")
  write.table(wrong_biomart, temp_file, sep = "\t", row.names = FALSE)
  
  expect_error(
    download_gene_data(temp_file, verbose = FALSE),
    "Expected columns not found"
  )
  
  file.remove(temp_file)
})

test_that("download_gene_data processes BioMart data correctly", {
  # Create mock BioMart data
  biomart_data <- data.frame(
    HGNC.symbol = c("GENE1", "GENE2", "", "GENE3", "GENE4"),
    Gene.start..bp. = c(1000, 5000, 7000, 10000, 50000),
    Gene.end..bp. = c(2000, 6000, 8000, 12000, 52000),
    Chromosome.scaffold.name = c("1", "1", "1", "2", "MT"), # Include one to filter out
    stringsAsFactors = FALSE
  )
  
  biomart_file <- tempfile(fileext = ".txt")
  write.table(biomart_data, biomart_file, sep = "\t", row.names = FALSE)
  
  output_file <- tempfile(fileext = ".txt")
  
  # Test processing
  result_file <- download_gene_data(
    biomart_file = biomart_file,
    output_file = output_file,
    compress = FALSE,
    verbose = FALSE
  )
  
  # Check output
  expect_true(file.exists(result_file))
  
  # Read and validate processed data
  processed_data <- read.table(result_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Should have 3 genes (filtered out empty HGNC and MT chromosome)
  expect_equal(nrow(processed_data), 3)
  expect_named(processed_data, c("HGNC", "GeneStart", "GeneEnd", "Chromosome"))
  expect_true(all(processed_data$HGNC != ""))
  expect_true(all(processed_data$Chromosome %in% c("1", "2")))
  
  # Clean up
  file.remove(biomart_file)
  file.remove(result_file)
})

test_that("download_gene_data handles compression", {
  # Create minimal BioMart data
  biomart_data <- data.frame(
    HGNC.symbol = c("GENE1"),
    Gene.start..bp. = c(1000),
    Gene.end..bp. = c(2000),
    Chromosome.scaffold.name = c("1"),
    stringsAsFactors = FALSE
  )
  
  biomart_file <- tempfile(fileext = ".txt")
  write.table(biomart_data, biomart_file, sep = "\t", row.names = FALSE)
  
  output_file <- tempfile(fileext = ".txt")
  
  # Test with compression
  result_file <- download_gene_data(
    biomart_file = biomart_file,
    output_file = output_file,
    compress = TRUE,
    verbose = FALSE
  )
  
  # Should have .gz extension
  expect_true(grepl("\\.gz$", result_file))
  expect_true(file.exists(result_file))
  
  # Should be able to read compressed file
  processed_data <- read.table(gzfile(result_file), header = TRUE, sep = "\t")
  expect_equal(nrow(processed_data), 1)
  
  # Clean up
  file.remove(biomart_file)
  file.remove(result_file)
})

# Helper function to check if we're online
skip_if_offline <- function() {
  if (!capabilities("http/ftp")) {
    skip("HTTP capabilities not available")
  }
  
  # Try a simple connection test
  tryCatch({
    con <- url("https://www.google.com", open = "rb", timeout = 2)
    close(con)
  }, error = function(e) {
    skip("Internet connection not available")
  })
}
