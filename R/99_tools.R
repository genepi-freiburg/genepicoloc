# Tools (supplementary functions) to facilitate other genepicoloc processes

handle_underflow <- function(P_vector) {
  if (is.character(P_vector)) {
    exp_indeces <- grepl("e-", P_vector, ignore.case = T)
    if (any(exp_indeces)) {
      significand_val <- as.numeric(gsub("(.*)[eE]-[0-9]+$", "\\1", P_vector[exp_indeces]))
      exponent_val <- as.numeric(gsub(".*[eE]-([0-9]+)$", "\\1", P_vector[exp_indeces]))
      P_vector[exp_indeces] <- as.character(-log10(significand_val) + exponent_val)
      P_vector[!exp_indeces] <- -log10(as.numeric(P_vector[!exp_indeces]))
      P_vector <- as.numeric(P_vector)
      return(P_vector)
    } else {
      stop("P column is character but scientific notation patterns 'e-' or 'E-' were not found")
    }
  } else {
    -log10(P_vector, na.rm=T)
  }
}


flip_alleles <- function(vec) {
  vec_out <- vec
  vec_out <- toupper(vec_out)
  vec_out[vec == "A"] <- "T"
  vec_out[vec == "T"] <- "A"
  vec_out[vec == "C"] <- "G"
  vec_out[vec == "G"] <- "C"
  return(vec_out)
}

cis_trans_annotation <- function(region_CHR_vec, region_BP_START_vec, region_BP_STOP_vec,
                                 gene_chr_vec, gene_start_vec,
                                 suggestive_window = 1e6) {
  cis_condition <- (gene_chr_vec == region_CHR_vec) &
    (gene_start_vec >= region_BP_START_vec & gene_start_vec <= region_BP_STOP_vec)
  suggestive_cis_condition <- (gene_chr_vec == region_CHR_vec) &
    (gene_start_vec >= region_BP_START_vec-suggestive_window &
       gene_start_vec <= region_BP_STOP_vec+suggestive_window)
  trans_condition <- ((gene_chr_vec != region_CHR_vec) | 
                        (gene_start_vec < region_BP_START_vec-suggestive_window) |
                        (gene_start_vec > region_BP_STOP_vec+suggestive_window))
  cis_trans <- ifelse(cis_condition, "cis",
                      ifelse(suggestive_cis_condition, "suggestive_cis",
                             ifelse(trans_condition, "trans", NA)))
  return(cis_trans)
}

# Temporary functions, not actively used
# params_df_to_list <- function(params_df) {
#   split(params_df, seq(nrow(params_df)))
# }
# params_df_to_chunks <- function(params_df,
#                                 chunk_size = 10000) {
#   nrow_var <- nrow(params_df)
#   N_chunks <- ceiling(nrow_var / chunk_size)
#   params_chunks <- lapply(1:N_chunks, function(i) {
#     i_end <- i*chunk_size
#     if (i == N_chunks) { i_end <- nrow_var }
#     params_df[((i-1)*chunk_size+1):i_end,]
#   })
#   return(params_chunks)
# }



