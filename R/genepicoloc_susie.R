#' Create LD matrix from genotype data using PLINK2
#' 
#' @description
#' Calculates a linkage disequilibrium (LD) correlation matrix for a set of
#' variants using PLINK2. Handles reference allele alignment and missing data.
#' 
#' @param sumstats Data.frame containing variant information with columns
#'   specified by Name_var and optionally A1_var.
#' @param bfile Character string. Path to PLINK binary fileset (without extension).
#'   Files .bed, .bim, and .fam must exist.
#' @param fixref Logical. Whether to fix reference alleles to match sumstats
#'   (default: FALSE). If TRUE, ensures LD calculations use the same effect allele.
#' @param keep Character string. Path to file with sample IDs to include
#'   (PLINK --keep format). Default: NULL (use all samples).
#' @param plink2 Character string. Path to PLINK2 executable (default: "plink2").
#' @param Name_var Character string. Column name for variant identifiers
#'   (default: "Name").
#' @param A1_var Character string. Column name for effect allele when fixref=TRUE
#'   (default: "A1").
#' 
#' @return A numeric matrix with variants as both row and column names.
#'   Values are correlation coefficients between -1 and 1.
#'   Variants with missing data are removed.
#' 
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts specified variants from the PLINK dataset
#'   \item Optionally aligns reference alleles with summary statistics
#'   \item Calculates unphased correlation matrix
#'   \item Removes variants with any missing correlations
#'   \item Cleans up temporary files
#' }
#' 
#' Missing data handling:
#' \itemize{
#'   \item Warns if >20% of variants have missing correlations
#'   \item Stops if >80% of variants have missing correlations
#'   \item Suggests using --geno and --mind filters if missingness is high
#' }
#' 
#' @examples
#' \dontrun{
#' # Basic LD matrix calculation
#' ld_matrix <- create_LD_matrix(
#'   sumstats = my_sumstats,
#'   bfile = "data/genotypes/chr1"
#' )
#' 
#' # With reference allele correction
#' ld_matrix <- create_LD_matrix(
#'   sumstats = my_sumstats,
#'   bfile = "data/genotypes/chr1",
#'   fixref = TRUE
#' )
#' 
#' # Using specific samples
#' ld_matrix <- create_LD_matrix(
#'   sumstats = my_sumstats,
#'   bfile = "data/genotypes/chr1",
#'   keep = "data/samples_european.txt"
#' )
#' }
#' 
#' @export
create_LD_matrix <- function(sumstats, bfile, fixref = FALSE, keep = NULL,
                             plink2 = "plink2", Name_var = "Name", A1_var = "A1") {
  
  # Input validation
  if (!is.data.frame(sumstats)) {
    stop("sumstats must be a data.frame")
  }
  
  if (!Name_var %in% colnames(sumstats)) {
    stop("Column '", Name_var, "' not found in sumstats")
  }
  
  if (fixref && !A1_var %in% colnames(sumstats)) {
    stop("Column '", A1_var, "' not found in sumstats (required when fixref=TRUE)")
  }
  
  # Check PLINK files exist
  plink_files <- paste0(bfile, c(".bed", ".bim", ".fam"))
  missing_files <- plink_files[!file.exists(plink_files)]
  if (length(missing_files) > 0) {
    stop("PLINK files not found: ", paste(missing_files, collapse = ", "))
  }
  
  # Check PLINK2 executable
  if (system(paste(plink2, "--version"), ignore.stdout = TRUE) != 0) {
    stop("PLINK2 not found at: ", plink2)
  }
  
  # Create temporary files
  tmp_base <- tempfile()
  variant_file <- paste0(tmp_base, "_variants.txt")
  
  # Ensure cleanup on exit
  on.exit({
    unlink(paste0(tmp_base, "*"))
  }, add = TRUE)
  
  # Build base PLINK command
  plink2_cmd <- paste(plink2, "--bfile", shQuote(bfile))
  
  # Add sample filter if provided
  if (!is.null(keep)) {
    if (!file.exists(keep)) {
      stop("Keep file not found: ", keep)
    }
    plink2_cmd <- paste(plink2_cmd, "--keep", shQuote(keep))
  }
  
  # Handle reference allele alignment
  if (fixref) {
    message("Aligning reference alleles with summary statistics")
    
    # Write variant and allele information
    write.table(
      sumstats[, c(Name_var, A1_var)], 
      variant_file, 
      row.names = FALSE, 
      col.names = FALSE, 
      quote = FALSE
    )
    
    # Create new dataset with aligned alleles
    align_cmd <- paste(
      plink2_cmd,
      "--extract", shQuote(variant_file),
      "--ref-allele", shQuote(variant_file),
      "--make-bed",
      "--out", shQuote(paste0(tmp_base, "_ref"))
    )
    
    if (system(align_cmd, ignore.stdout = TRUE) != 0) {
      stop("Failed to align reference alleles")
    }
    
    # Update file paths
    bfile <- paste0(tmp_base, "_ref")
    plink2_cmd <- paste(plink2, "--bfile", shQuote(bfile))
    
  } else {
    # Just write variant list
    write.table(
      sumstats[, Name_var], 
      variant_file, 
      row.names = FALSE, 
      col.names = FALSE, 
      quote = FALSE
    )
  }
  
  # Calculate LD matrix
  message("Calculating LD matrix for ", nrow(sumstats), " variants")
  
  ld_cmd <- paste(
    plink2_cmd,
    "--extract", shQuote(variant_file),
    "--r-unphased square ref-based",
    "--out", shQuote(paste0(tmp_base, "_r_mat"))
  )
  
  if (system(ld_cmd, ignore.stdout = TRUE) != 0) {
    stop("Failed to calculate LD matrix")
  }
  
  # Read results
  ld_file <- paste0(tmp_base, "_r_mat.unphased.vcor1")
  var_file <- paste0(tmp_base, "_r_mat.unphased.vcor1.vars")
  
  if (!file.exists(ld_file) || !file.exists(var_file)) {
    stop("LD calculation did not produce expected output files")
  }
  
  # Load LD matrix
  r_mat <- as.matrix(read.table(ld_file, header = FALSE))
  r_mat_vars <- read.table(var_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  
  # Set dimension names
  dimnames(r_mat) <- list(r_mat_vars, r_mat_vars)
  
  # Check for and remove variants with missing data
  na_per_row <- rowSums(is.na(r_mat))
  has_na <- na_per_row > 0
  n_removed <- sum(has_na)
  
  if (n_removed > 0) {
    prop_removed <- n_removed / nrow(r_mat)
    
    if (prop_removed > 0.8) {
      stop("More than 80% of variants have missing LD values. ",
           "Check genotype quality and consider filtering with --geno and --mind")
    }
    
    if (prop_removed > 0.2) {
      warning("Removed ", n_removed, " variants (", 
              round(prop_removed * 100, 1), "%) due to missing LD values. ",
              "Consider filtering genotypes with --geno and --mind")
    }
    
    message("Removed ", n_removed, " variants with missing LD values")
    r_mat <- r_mat[!has_na, !has_na]
  }
  
  # Final validation
  if (nrow(r_mat) == 0) {
    stop("No variants remaining after filtering")
  }
  
  message("LD matrix calculated successfully: ", 
          nrow(r_mat), " x ", ncol(r_mat))
  
  return(r_mat)
}


#' Run SuSiE with individual-level data
#' 
#' @description
#' Performs fine-mapping using SuSiE (Sum of Single Effects) with individual-level
#' genotype and phenotype data. Handles data alignment and missing values.
#' 
#' @param raw_gt Data.frame. PLINK raw format genotypes with columns:
#'   FID, IID, PAT, MAT, SEX, PHENOTYPE, then genotype columns (0/1/2 coding).
#' @param pheno Data.frame. Phenotype data with at least IID column and phenotype column.
#' @param pheno_name Character string. Name of the phenotype column to analyze.
#' @param coverage Numeric. Coverage for credible sets (default: 0.95).
#' @param scale_genotypes Logical. Whether to standardize genotypes (default: FALSE).
#' 
#' @return A susie object containing:
#'   \itemize{
#'     \item alpha: Posterior inclusion probabilities
#'     \item mu: Posterior mean effects
#'     \item sets: Credible sets at specified coverage
#'     \item pip: Posterior inclusion probability for each variant
#'   }
#' 
#' @details
#' The function:
#' \enumerate{
#'   \item Validates presence of IID columns
#'   \item Removes samples with missing phenotypes
#'   \item Removes variants with any missing genotypes
#'   \item Aligns samples between genotype and phenotype data
#'   \item Optionally scales genotypes to mean 0, variance 1
#'   \item Runs SuSiE fine-mapping
#' }
#' 
#' @importFrom susieR susie
#' 
#' @examples
#' \dontrun{
#' # Load genotype and phenotype data
#' raw_gt <- read.table("genotypes.raw", header = TRUE)
#' pheno <- read.table("phenotypes.txt", header = TRUE)
#' 
#' # Run SuSiE fine-mapping
#' susie_results <- genepicoloc_susie(
#'   raw_gt = raw_gt,
#'   pheno = pheno,
#'   pheno_name = "LDL_cholesterol",
#'   coverage = 0.95
#' )
#' 
#' # Extract credible sets
#' print(susie_results$sets)
#' }
#' 
#' @export
genepicoloc_susie <- function(raw_gt, pheno, pheno_name,
                              coverage = 0.95,
                              scale_genotypes = FALSE) {
  
  # Input validation
  if (!is.data.frame(raw_gt)) {
    stop("raw_gt must be a data.frame")
  }
  
  if (!is.data.frame(pheno)) {
    stop("pheno must be a data.frame")
  }
  
  if (!"IID" %in% colnames(raw_gt)) {
    stop("No IID column found in genotype data")
  }
  
  if (!"IID" %in% colnames(pheno)) {
    stop("No IID column found in phenotype data")
  }
  
  if (!pheno_name %in% colnames(pheno)) {
    stop("Phenotype '", pheno_name, "' not found in pheno data")
  }
  
  if (coverage <= 0 || coverage >= 1) {
    stop("coverage must be between 0 and 1")
  }
  
  # Prepare phenotype data
  message("Preparing phenotype data")
  pheno <- pheno[, c("IID", pheno_name)]
  n_total <- nrow(pheno)
  
  # Remove missing phenotypes
  pheno <- subset(pheno, !is.na(pheno[[pheno_name]]))
  n_missing <- n_total - nrow(pheno)
  
  if (n_missing > 0) {
    message("Removed ", n_missing, " samples with missing phenotypes")
  }
  
  if (nrow(pheno) == 0) {
    stop("No samples with non-missing phenotypes")
  }
  
  # Prepare genotype data
  message("Preparing genotype data")
  raw_gt <- as.data.frame(raw_gt)
  
  # Identify genotype columns (after the 6 PLINK columns)
  if (ncol(raw_gt) <= 6) {
    stop("No genotype columns found in raw_gt")
  }
  
  geno_cols <- 7:ncol(raw_gt)
  
  # Remove variants with any missing genotypes
  n_variants_total <- length(geno_cols)
  has_missing <- sapply(geno_cols, function(col) any(is.na(raw_gt[[col]])))
  
  if (any(has_missing)) {
    n_removed <- sum(has_missing)
    message("Removed ", n_removed, " variants with missing genotypes")
    keep_cols <- c(1:6, geno_cols[!has_missing])
    raw_gt <- raw_gt[, keep_cols]
    geno_cols <- 7:ncol(raw_gt)
  }
  
  # Align samples
  message("Aligning samples between genotype and phenotype data")
  
  # Subset to common samples
  common_iids <- intersect(raw_gt$IID, pheno$IID)
  
  if (length(common_iids) == 0) {
    stop("No common samples between genotype and phenotype data")
  }
  
  if (length(common_iids) < length(pheno$IID)) {
    message("Using ", length(common_iids), " samples present in both datasets")
  }
  
  # Subset and order both datasets
  raw_gt <- raw_gt[raw_gt$IID %in% common_iids, ]
  pheno <- pheno[pheno$IID %in% common_iids, ]
  
  # Match order
  pheno <- pheno[match(raw_gt$IID, pheno$IID), ]
  
  # Verify alignment
  if (!all(pheno$IID == raw_gt$IID)) {
    stop("Failed to align samples between datasets")
  }
  
  # Extract matrices
  pheno_vec <- pheno[[pheno_name]]
  geno_mat <- as.matrix(raw_gt[, geno_cols])
  
  # Scale genotypes if requested
  if (scale_genotypes) {
    message("Scaling genotypes to mean 0, variance 1")
    geno_mat <- scale(geno_mat)
    
    # Check for constant genotypes
    const_vars <- which(is.na(geno_mat[1, ]))
    if (length(const_vars) > 0) {
      warning("Removed ", length(const_vars), " constant variants after scaling")
      geno_mat <- geno_mat[, -const_vars]
    }
  }
  
  # Run SuSiE
  message("Running SuSiE with ", nrow(geno_mat), " samples and ", 
          ncol(geno_mat), " variants")
  
  susie_out <- susieR::susie(
    X = geno_mat, 
    y = pheno_vec,
    coverage = coverage,
    verbose = TRUE
  )
  
  # Add variant names if available
  if (!is.null(colnames(geno_mat))) {
    names(susie_out$pip) <- colnames(geno_mat)
  }
  
  message("SuSiE analysis complete")
  
  return(susie_out)
}


#' Scale genotype matrix
#' 
#' @description
#' Standardizes genotype matrix to have mean 0 and variance 1 for each variant.
#' This is often required for statistical methods that assume standardized inputs.
#' 
#' @param raw_gt_mat Numeric matrix or data.frame. Genotype matrix with
#'   samples in rows and variants in columns.
#' 
#' @return Numeric matrix with scaled genotypes. Constant variants (zero variance)
#'   will have NA values.
#' 
#' @details
#' The function uses R's scale() function which:
#' \itemize{
#'   \item Centers by subtracting the column mean
#'   \item Scales by dividing by the column standard deviation
#'   \item Returns NA for constant columns (zero variance)
#' }
#' 
#' @export
scale_genotypes <- function(raw_gt_mat) {
  
  # Convert to matrix if needed
  if (!is.matrix(raw_gt_mat)) {
    if (is.data.frame(raw_gt_mat)) {
      raw_gt_mat <- as.matrix(raw_gt_mat)
    } else {
      stop("Input must be a matrix or data.frame")
    }
  }
  
  # Check numeric
  if (!is.numeric(raw_gt_mat)) {
    stop("Genotype matrix must be numeric")
  }
  
  # Scale
  scaled_mat <- scale(raw_gt_mat)
  
  # Check for constant variants
  const_vars <- colSums(is.na(scaled_mat[1, , drop = FALSE])) > 0
  
  if (any(const_vars)) {
    warning(sum(const_vars), " variants have zero variance and cannot be scaled")
  }
  
  return(scaled_mat)
}


#' Run SuSiE RSS with summary statistics
#' 
#' @description
#' Performs fine-mapping using SuSiE RSS (Regression with Summary Statistics)
#' which requires only summary statistics and an LD matrix.
#' 
#' @param sumstats Data.frame. Summary statistics with columns specified by
#'   the name parameters.
#' @param LD_matrix Numeric matrix. LD correlation matrix with variants in the
#'   same order as sumstats. Row and column names must match sumstats Name column.
#' @param Name_name Character string. Column name for variant identifiers
#'   (default: "Name").
#' @param BETA_name Character string. Column name for effect sizes
#'   (default: "BETA").
#' @param SE_name Character string. Column name for standard errors
#'   (default: "SE").
#' @param N_name Character string. Column name for sample sizes
#'   (default: "N").
#' @param N Numeric. Sample size override. If not provided, uses median from
#'   the N column.
#' @param coverage Numeric. Coverage for credible sets (default: 0.95).
#' @param ... Additional arguments passed to susie_rss().
#' 
#' @return A susie object containing fine-mapping results.
#' 
#' @details
#' Key parameters that can be passed via ...:
#' \itemize{
#'   \item L: Maximum number of causal variants (default: 10)
#'   \item max_iter: Maximum iterations (default: 100)
#'   \item estimate_residual_variance: Whether to estimate residual variance
#'   \item refine: Whether to refine the model
#' }
#' 
#' The function requires that variant order in sumstats matches the order
#' in the LD matrix exactly.
#' 
#' @importFrom susieR susie_rss
#' 
#' @examples
#' \dontrun{
#' # Fine-mapping with summary statistics
#' susie_results <- genepicoloc_susie_rss(
#'   sumstats = my_sumstats,
#'   LD_matrix = my_ld_matrix,
#'   N = 10000,
#'   coverage = 0.95,
#'   L = 10
#' )
#' }
#' 
#' @export
genepicoloc_susie_rss <- function(sumstats, LD_matrix,
                                  Name_name = "Name", 
                                  BETA_name = "BETA",
                                  SE_name = "SE", 
                                  N_name = "N",
                                  N, 
                                  coverage = 0.95, 
                                  ...) {
  
  # Input validation
  if (!is.data.frame(sumstats)) {
    stop("sumstats must be a data.frame")
  }
  
  if (!is.matrix(LD_matrix)) {
    stop("LD_matrix must be a matrix")
  }
  
  # Check required columns exist
  required_cols <- c(Name_name, BETA_name, SE_name)
  missing_cols <- setdiff(required_cols, colnames(sumstats))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check dimensions
  if (nrow(sumstats) != nrow(LD_matrix) || nrow(LD_matrix) != ncol(LD_matrix)) {
    stop("Dimension mismatch: sumstats has ", nrow(sumstats), " variants, ",
         "LD matrix is ", nrow(LD_matrix), " x ", ncol(LD_matrix))
  }
  
  # Check variant order matches
  if (!identical(sumstats[[Name_name]], colnames(LD_matrix)) ||
      !identical(sumstats[[Name_name]], rownames(LD_matrix))) {
    stop("Variant order mismatch between sumstats and LD matrix. ",
         "Ensure sumstats[[", Name_name, "]] matches LD matrix row/column names")
  }
  
  # Handle sample size
  if (missing(N)) {
    if (!N_name %in% colnames(sumstats)) {
      stop("N not provided and column '", N_name, "' not found in sumstats")
    }
    
    N <- median(sumstats[[N_name]], na.rm = TRUE)
    
    if (is.na(N)) {
      stop("Could not determine sample size from sumstats")
    }
    
    message("Using median sample size: N = ", round(N))
  }
  
  # Validate sample size
  if (!is.numeric(N) || length(N) != 1 || N <= 0) {
    stop("N must be a single positive number")
  }
  
  # Check for missing values
  if (any(is.na(sumstats[[BETA_name]]))) {
    stop("Missing values in effect sizes (", BETA_name, ")")
  }
  
  if (any(is.na(sumstats[[SE_name]]))) {
    stop("Missing values in standard errors (", SE_name, ")")
  }
  
  if (any(is.na(LD_matrix))) {
    stop("Missing values in LD matrix")
  }
  
  # Check LD matrix properties
  if (!isSymmetric(LD_matrix)) {
    warning("LD matrix is not symmetric, using (M + t(M))/2")
    LD_matrix <- (LD_matrix + t(LD_matrix)) / 2
  }
  
  # Run SuSiE RSS
  message("Running SuSiE RSS with ", nrow(sumstats), " variants and N = ", round(N))
  
  fitted_rss <- susieR::susie_rss(
    bhat = sumstats[[BETA_name]],
    shat = sumstats[[SE_name]],
    n = N,
    R = LD_matrix,
    coverage = coverage,
    verbose = TRUE,
    ...
  )
  
  # Add variant names if available
  if (!is.null(sumstats[[Name_name]])) {
    names(fitted_rss$pip) <- sumstats[[Name_name]]
  }
  
  message("SuSiE RSS analysis complete")
  
  return(fitted_rss)
}


#' Wrapper for SuSiE RSS analysis with automatic LD calculation
#' 
#' @description
#' Convenience function that combines LD matrix calculation and SuSiE RSS
#' fine-mapping in a single step. Handles reference allele alignment and
#' variant filtering.
#' 
#' @param sumstats Data.frame. Summary statistics in standardized format.
#' @param bfile Character string. Path to PLINK binary fileset.
#' @param N Numeric. Sample size (default: median from sumstats).
#' @param coverage Numeric. Coverage for credible sets (default: 0.95).
#' @param rm_indels Logical. Whether to remove insertions/deletions
#'   (default: FALSE).
#' 
#' @return A susie object containing fine-mapping results.
#' 
#' @details
#' This wrapper:
#' \enumerate{
#'   \item Optionally removes indels (non-SNP variants)
#'   \item Calculates LD matrix with reference allele correction
#'   \item Ensures variants match between sumstats and LD matrix
#'   \item Runs SuSiE RSS fine-mapping
#' }
#' 
#' Indels are identified as variants where A1 or A2 is not in {A,T,G,C}.
#' 
#' @examples
#' \dontrun{
#' # Basic usage
#' susie_results <- susie_rss_wrapper(
#'   sumstats = my_sumstats,
#'   bfile = "data/genotypes/chr1",
#'   coverage = 0.95
#' )
#' 
#' # Remove indels and specify sample size
#' susie_results <- susie_rss_wrapper(
#'   sumstats = my_sumstats,
#'   bfile = "data/genotypes/chr1",
#'   N = 50000,
#'   rm_indels = TRUE
#' )
#' }
#' 
#' @export
susie_rss_wrapper <- function(sumstats, bfile, N = NULL, coverage = 0.95,
                              rm_indels = FALSE) {
  
  # Input validation
  if (!is.data.frame(sumstats)) {
    stop("sumstats must be a data.frame")
  }
  
  required_cols <- c("Name", "A1", "A2", "BETA", "SE", "N")
  missing_cols <- setdiff(required_cols, colnames(sumstats))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Remove indels if requested
  if (rm_indels) {
    n_before <- nrow(sumstats)
    
    # Identify SNPs (single nucleotide variants)
    is_snp <- sumstats$A1 %in% c("A", "T", "G", "C") & 
      sumstats$A2 %in% c("A", "T", "G", "C")
    
    sumstats <- sumstats[is_snp, ]
    
    n_removed <- n_before - nrow(sumstats)
    if (n_removed > 0) {
      message("Removed ", n_removed, " indels, ", 
              nrow(sumstats), " SNPs remaining")
    }
    
    if (nrow(sumstats) == 0) {
      stop("No SNPs remaining after removing indels")
    }
  }
  
  # Determine sample size
  if (is.null(N)) {
    N <- median(sumstats$N, na.rm = TRUE)
    
    if (is.na(N)) {
      stop("Could not determine sample size from sumstats")
    }
    
    message("Using median sample size: N = ", round(N))
  }
  
  # Calculate LD matrix with reference allele alignment
  message("Calculating LD matrix")
  LD_matrix <- create_LD_matrix(
    sumstats = sumstats, 
    bfile = bfile, 
    fixref = TRUE
  )
  
  # Ensure sumstats matches LD matrix
  common_variants <- intersect(sumstats$Name, colnames(LD_matrix))
  
  if (length(common_variants) == 0) {
    stop("No common variants between sumstats and LD matrix")
  }
  
  if (length(common_variants) < nrow(sumstats)) {
    n_missing <- nrow(sumstats) - length(common_variants)
    message("Note: ", n_missing, " variants in sumstats not found in LD matrix")
  }
  
  # Subset and reorder sumstats to match LD matrix
  sumstats <- sumstats[match(colnames(LD_matrix), sumstats$Name), ]
  
  # Verify alignment
  if (!identical(sumstats$Name, colnames(LD_matrix))) {
    stop("Failed to align sumstats with LD matrix")
  }
  
  # Run SuSiE RSS
  susie_obj <- genepicoloc_susie_rss(
    sumstats = sumstats,
    LD_matrix = LD_matrix,
    N = N,
    coverage = coverage
  )
  
  return(susie_obj)
}