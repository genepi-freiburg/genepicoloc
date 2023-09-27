#' @param olink_protein_map_3k_v1_file path to olink_protein_map_3k_v1.
#' @return data.frame with processed olink_protein_map_3k_v1.
#' @examples
#' olink_protein_map_3k_v1(olink_protein_map_3k_v1_path)
#' @export
olink_annotation <- function(olink_protein_map_3k_v1_file, coloc_out) {
  olink_protein_map_3k_v1 <- read.delim(olink_protein_map_3k_v1_file, sep="\t")
  olink_protein_map_3k_v1$UKBPPP_ProteinID <- gsub(":", "_", olink_protein_map_3k_v1$UKBPPP_ProteinID)
  olink_protein_map_3k_v1$multiple_genes_per_OID <- 0
  olink_protein_map_3k_v1[which(duplicated(olink_protein_map_3k_v1$OlinkID)),]$multiple_genes_per_OID <- 1
  selected_cols <- c("OlinkID", "olink_target_fullname", "UniProt", "Assay",
                     "HGNC.symbol", "ensembl_id", "chr", "gene_start", "gene_end",
                     "multiple_genes_per_OID")
  olink_protein_map_3k_v1 <- olink_protein_map_3k_v1[,selected_cols]
  colnames(olink_protein_map_3k_v1)[-1] <- paste0("Olink_", colnames(olink_protein_map_3k_v1)[-1])
  coloc_out$OlinkID <- gsub(".*(OID[0-9]+).*", "\\1",coloc_out$sumstats_2_file)
  coloc_out <- merge(coloc_out, by.x="OlinkID",
                     olink_protein_map_3k_v1, by.y="OlinkID",
                     sort = FALSE)[, union(names(coloc_out), names(olink_protein_map_3k_v1))]
  # stopifnot(unique(coloc_out$CHR_var) %in% unique(coloc_out$Olink_chr))
  coloc_out$cis_trans <- "trans"
  cis_condition <- (coloc_out$CHR_var == coloc_out$Olink_chr) & ((coloc_out$Olink_gene_start >= coloc_out$BP_START_var & coloc_out$Olink_gene_start <= coloc_out$BP_STOP_var) | (coloc_out$Olink_gene_end >= coloc_out$BP_START_var & coloc_out$Olink_gene_end <= coloc_out$BP_STOP_var))
  if (any(cis_condition)) {
    coloc_out[cis_condition, ]$cis_trans <- "cis"
  }
  suggestive_cis_condition <- (coloc_out$CHR_var == coloc_out$Olink_chr) & ((coloc_out$Olink_gene_start >= coloc_out$BP_START_var-1e6 & coloc_out$Olink_gene_start <= coloc_out$BP_STOP_var+1e6) | (coloc_out$Olink_gene_end >= coloc_out$BP_START_var-1e6 & coloc_out$Olink_gene_end <= coloc_out$BP_STOP_var+1e6))
  if (any(suggestive_cis_condition)) {
    coloc_out[(!cis_condition) & suggestive_cis_condition, ]$cis_trans <- "suggestive_cis"
  }
  return(coloc_out)
}

