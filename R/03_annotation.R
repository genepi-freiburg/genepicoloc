#' @param annotation_file path to annotation file
#' @return data.frame with processed annotation file.
#' @examples
#' Under development
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
  coloc_out$OlinkID <- gsub(".*(OID[0-9]+).*", "\\1", coloc_out$sumstats_2_file)
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

#' @param annotation_file path to annotation file
#' @return data.frame with processed annotation file.
#' @examples
#' Under development
#' @export
somascan_annotation <- function(annotation_file, coloc_out) {
  annotation_df <- read.csv(annotation_file)
  annotation_df$multiple_genes_per_protein <- 0
  # annotation_df[which(duplicated(annotation_df$SeqId)),]$multiple_genes_per_protein <- 1
  # selected_cols <- c("SeqId", "Protein..short.name.", "Protein..full.name.", "Gene",
  #                    "UniProt", "Type", "Ensembl.Gene.ID", "chr", "gene_start", "gene_end",
  #                    "multiple_genes_per_protein")
  # annotation_df <- annotation_df[,selected_cols]
  prefix <- "Soma_"
  colnames(annotation_df)[-1] <- paste0(prefix, colnames(annotation_df)[-1])
  # coloc
  coloc_out$SeqId <- gsub(".*NG2021/([0-9]+_[0-9]+)_.*", "\\1", coloc_out$sumstats_2_file)
  coloc_out <- merge(coloc_out, by.x="SeqId",
                     annotation_df, by.y="SeqId",
                     sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  # stopifnot(unique(coloc_out$CHR_var) %in% unique(coloc_out$Olink_chr))
  coloc_out$cis_trans <- "trans"
  cis_condition <- (coloc_out[[paste0(prefix, "chr")]] == coloc_out$CHR_var) & ((coloc_out[[paste0(prefix, "gene_start")]] >= coloc_out$BP_START_var & coloc_out[[paste0(prefix, "gene_start")]] <= coloc_out$BP_STOP_var) | (coloc_out[[paste0(prefix, "gene_end")]] >= coloc_out$BP_START_var & coloc_out[[paste0(prefix, "gene_end")]] <= coloc_out$BP_STOP_var))
  if (any(cis_condition)) {
    coloc_out[cis_condition, ]$cis_trans <- "cis"
  }
  suggestive_cis_condition <- (coloc_out$CHR_var == coloc_out[[paste0(prefix, "chr")]]) & ((coloc_out[[paste0(prefix, "gene_start")]] >= coloc_out$BP_START_var-1e6 & coloc_out[[paste0(prefix, "gene_start")]] <= coloc_out$BP_STOP_var+1e6) | (coloc_out[[paste0(prefix, "gene_end")]] >= coloc_out$BP_START_var-1e6 & coloc_out[[paste0(prefix, "gene_end")]] <= coloc_out$BP_STOP_var+1e6))
  if (any(suggestive_cis_condition)) {
    coloc_out[(!cis_condition) & suggestive_cis_condition, ]$cis_trans <- "suggestive_cis"
  }
  return(coloc_out)
}



