# Harmonize and merge summary statistics from different sources
#
# Handles allele matching, strand flipping, and effect direction alignment
# when combining GWAS data from different studies.


#' Harmonize two sets of summary statistics by genomic position
#'
#' Merges two summary statistics data.tables by chromosome and position,
#' aligning alleles and flipping effect directions where necessary.
#' Handles cases where alleles are in different order or on different strands.
#'
#' @param sumstats_1 data.table with columns: CHR, POS, A1, A2, BETA, SE, and optionally AF, Name, rsID
#' @param sumstats_2 data.table with columns: CHR, POS, A1, A2, BETA, SE, and optionally AF, Name, rsID
#' @param match_by Character. Matching strategy: "position" (default), "rsid", or "name".
#' @param action Character. What to return: "merge" returns merged data.table,
#'   "datasets" returns list of two coloc-ready dataset lists.
#' @param type_1 Character. Type of sumstats_1: "quant" or "cc" (for action="datasets")
#' @param type_2 Character. Type of sumstats_2: "quant" or "cc" (for action="datasets")
#' @param sdY_1 Numeric. Trait SD for sumstats_1 (for action="datasets")
#' @param sdY_2 Numeric. Trait SD for sumstats_2 (for action="datasets")
#'
#' @return Depending on action:
#'   "merge": merged data.table with suffixed columns (.1, .2)
#'   "datasets": list with d1 and d2 ready for coloc_abf()
#'
#' @details
#' Allele matching logic:
#' \enumerate{
#'   \item Direct match: A1_1==A1_2 and A2_1==A2_2 -> keep as-is
#'   \item Flipped: A1_1==A2_2 and A2_1==A1_2 -> flip BETA_2 sign, swap AF_2
#'   \item Strand complement: try complement (A<->T, C<->G) then repeat 1-2
#'   \item No match: variant excluded
#' }
#'
#' @export
harmonize_sumstats <- function(sumstats_1, sumstats_2,
                                match_by = "position",
                                action = "datasets",
                                type_1 = "quant", type_2 = "quant",
                                sdY_1 = 1, sdY_2 = 1) {

  if (!data.table::is.data.table(sumstats_1)) sumstats_1 <- data.table::as.data.table(sumstats_1)
  if (!data.table::is.data.table(sumstats_2)) sumstats_2 <- data.table::as.data.table(sumstats_2)

  # Validate required columns
  req <- c("BETA", "SE")
  for (col in req) {
    if (!col %in% names(sumstats_1)) stop("sumstats_1 missing column: ", col)
    if (!col %in% names(sumstats_2)) stop("sumstats_2 missing column: ", col)
  }

  # --- Build merge key ---
  if (match_by == "position") {
    for (col in c("CHR", "POS")) {
      if (!col %in% names(sumstats_1)) stop("sumstats_1 missing column: ", col)
      if (!col %in% names(sumstats_2)) stop("sumstats_2 missing column: ", col)
    }
    sumstats_1[, .merge_key := paste0(CHR, ":", POS)]
    sumstats_2[, .merge_key := paste0(CHR, ":", POS)]
  } else if (match_by == "rsid") {
    rsid_col_1 <- intersect(c("rsID", "rsid", "RSID", "rs"), names(sumstats_1))[1]
    rsid_col_2 <- intersect(c("rsID", "rsid", "RSID", "rs"), names(sumstats_2))[1]
    if (is.na(rsid_col_1) || is.na(rsid_col_2)) stop("rsID column not found")
    sumstats_1[, .merge_key := get(rsid_col_1)]
    sumstats_2[, .merge_key := get(rsid_col_2)]
  } else if (match_by == "name") {
    name_col_1 <- intersect(c("Name", "name", "SNP"), names(sumstats_1))[1]
    name_col_2 <- intersect(c("Name", "name", "SNP"), names(sumstats_2))[1]
    if (is.na(name_col_1) || is.na(name_col_2)) stop("Name column not found")
    sumstats_1[, .merge_key := get(name_col_1)]
    sumstats_2[, .merge_key := get(name_col_2)]
  }

  # Remove duplicates on merge key (keep first)
  s1 <- sumstats_1[!duplicated(.merge_key)]
  s2 <- sumstats_2[!duplicated(.merge_key)]

  # Merge
  shared_keys <- intersect(s1$.merge_key, s2$.merge_key)
  if (length(shared_keys) == 0) {
    warning("No shared variants found")
    return(if (action == "datasets") list(d1 = NULL, d2 = NULL) else data.table::data.table())
  }

  idx1 <- match(shared_keys, s1$.merge_key)
  idx2 <- match(shared_keys, s2$.merge_key)

  # --- Allele harmonization ---
  has_alleles <- all(c("A1", "A2") %in% names(s1)) && all(c("A1", "A2") %in% names(s2))

  flip_beta <- rep(1, length(shared_keys))
  allele_match <- rep("none", length(shared_keys))

  if (has_alleles) {
    a1_1 <- toupper(s1$A1[idx1])
    a2_1 <- toupper(s1$A2[idx1])
    a1_2 <- toupper(s2$A1[idx2])
    a2_2 <- toupper(s2$A2[idx2])

    # Complement lookup
    comp <- c(A = "T", T = "A", C = "G", G = "C")
    complement <- function(x) {
      sapply(strsplit(x, ""), function(bases) paste(comp[bases], collapse = ""))
    }

    # Direct match
    direct <- a1_1 == a1_2 & a2_1 == a2_2
    allele_match[direct] <- "direct"

    # Flipped match (A1<->A2)
    flipped <- !direct & a1_1 == a2_2 & a2_1 == a1_2
    allele_match[flipped] <- "flipped"
    flip_beta[flipped] <- -1

    # Strand complement
    remaining <- allele_match == "none"
    if (any(remaining)) {
      ca1_2 <- complement(a1_2[remaining])
      ca2_2 <- complement(a2_2[remaining])

      comp_direct <- a1_1[remaining] == ca1_2 & a2_1[remaining] == ca2_2
      comp_flipped <- a1_1[remaining] == ca2_2 & a2_1[remaining] == ca1_2

      idx_rem <- which(remaining)
      allele_match[idx_rem[comp_direct]] <- "complement"
      allele_match[idx_rem[comp_flipped]] <- "complement_flipped"
      flip_beta[idx_rem[comp_flipped]] <- -1
    }

    # Filter to matched variants
    matched <- allele_match != "none"
    n_removed <- sum(!matched)
    if (n_removed > 0) {
      message("  Removed ", n_removed, " variants with incompatible alleles (",
              round(n_removed / length(shared_keys) * 100, 1), "%)")
    }

    shared_keys <- shared_keys[matched]
    idx1 <- idx1[matched]
    idx2 <- idx2[matched]
    flip_beta <- flip_beta[matched]
    allele_match <- allele_match[matched]
  }

  n_flipped <- sum(flip_beta == -1)
  message("  Harmonized: ", length(shared_keys), " variants (",
          sum(allele_match == "direct"), " direct, ",
          n_flipped, " flipped",
          if (any(grepl("complement", allele_match))) paste0(", ", sum(grepl("complement", allele_match)), " strand-complemented") else "",
          ")")

  # --- Build output ---
  # Create variant identifier for coloc
  if ("Name" %in% names(s1)) {
    snp_ids <- s1$Name[idx1]
  } else if (match_by == "position") {
    snp_ids <- shared_keys
  } else {
    snp_ids <- shared_keys
  }

  if (action == "datasets") {
    d1 <- list(
      beta = s1$BETA[idx1],
      varbeta = s1$SE[idx1]^2,
      snp = snp_ids,
      type = type_1
    )
    if (type_1 == "quant") d1$sdY <- sdY_1

    d2 <- list(
      beta = s2$BETA[idx2] * flip_beta,
      varbeta = s2$SE[idx2]^2,
      snp = snp_ids,
      type = type_2
    )
    if (type_2 == "quant") d2$sdY <- sdY_2

    # Clean up temp columns
    sumstats_1[, .merge_key := NULL]
    sumstats_2[, .merge_key := NULL]

    return(list(d1 = d1, d2 = d2, n_harmonized = length(shared_keys),
                n_flipped = n_flipped))
  }

  # action == "merge"
  merged <- data.table::data.table(
    snp = snp_ids,
    CHR = s1$CHR[idx1],
    POS = s1$POS[idx1],
    A1 = s1$A1[idx1],
    A2 = s1$A2[idx1],
    BETA.1 = s1$BETA[idx1],
    SE.1 = s1$SE[idx1],
    BETA.2 = s2$BETA[idx2] * flip_beta,
    SE.2 = s2$SE[idx2],
    allele_match = allele_match
  )

  # Add optional columns
  if ("AF" %in% names(s1)) merged[, AF.1 := s1$AF[idx1]]
  if ("AF" %in% names(s2)) merged[, AF.2 := ifelse(flip_beta == -1, 1 - s2$AF[idx2], s2$AF[idx2])]
  if ("nlog10P" %in% names(s1)) merged[, nlog10P.1 := s1$nlog10P[idx1]]
  if ("nlog10P" %in% names(s2)) merged[, nlog10P.2 := s2$nlog10P[idx2]]

  sumstats_1[, .merge_key := NULL]
  sumstats_2[, .merge_key := NULL]

  merged
}
