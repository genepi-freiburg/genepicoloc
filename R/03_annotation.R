#' Olink_pGWAS annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
Olink_pGWAS <- function(annotation_file, coloc_out) {
  annotation_df <- read.delim(annotation_file, sep="\t")
  annotation_df$multiple_genes_per_OID <- 0
  annotation_df[which(duplicated(annotation_df$OlinkID)),]$multiple_genes_per_OID <- 1
  selected_cols <- c("OlinkID", "olink_target_fullname", "UniProt", "Assay",
                     "HGNC.symbol", "ensembl_id", "chr", "gene_start", "gene_end",
                     "multiple_genes_per_OID")
  annotation_df <- annotation_df[,selected_cols]
  colnames(annotation_df)[-1] <- paste0("Olink_", colnames(annotation_df)[-1])
  coloc_out$OlinkID <- gsub(".*(OID[0-9]+).*", "\\1", coloc_out$sumstats_2_file)
  coloc_out <- merge(coloc_out, by.x="OlinkID",
                     annotation_df, by.y="OlinkID",
                     sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
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

#' somascan annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
somascan_annotation <- function(annotation_file, coloc_out) {
  annotation_df <- read.csv(annotation_file)
  annotation_df$chr <- gsub("chr", "", annotation_df$chr)
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

#' ARIC pGWAS annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
ARIC_pGWAS_annotation <- function(annotation_file, coloc_out) {
  annotation_df <- read.delim(annotation_file)
  colnames(annotation_df)[colnames(annotation_df) == "chromosome_name"] <- "chr"
  colnames(annotation_df)[colnames(annotation_df) == "transcription_start_site"] <- "gene_start"
  prefix <- "ARIC_"
  colnames(annotation_df)[-1] <- paste0(prefix, colnames(annotation_df)[-1])
  # coloc
  coloc_out$seqid_in_sample <- gsub(".PHENO1.glm.linear.gz", "", basename(coloc_out$sumstats_2_file))
  coloc_out <- merge(coloc_out, by="seqid_in_sample",
                     annotation_df,
                     sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  coloc_out$cis_trans <- "trans"
  cis_condition <- (coloc_out[[paste0(prefix, "chr")]] == coloc_out$CHR_var) & ((coloc_out[[paste0(prefix, "gene_start")]] >= coloc_out$BP_START_var & coloc_out[[paste0(prefix, "gene_start")]] <= coloc_out$BP_STOP_var))
  if (any(cis_condition)) {
    coloc_out[cis_condition, ]$cis_trans <- "cis"
  }
  suggestive_cis_condition <- (coloc_out$CHR_var == coloc_out[[paste0(prefix, "chr")]]) & ((coloc_out[[paste0(prefix, "gene_start")]] >= coloc_out$BP_START_var-1e6 & coloc_out[[paste0(prefix, "gene_start")]] <= coloc_out$BP_STOP_var+1e6))
  if (any(suggestive_cis_condition)) {
    coloc_out[(!cis_condition) & suggestive_cis_condition, ]$cis_trans <- "suggestive_cis"
  }
  coloc_out[is.na(coloc_out$PP.H4.abf), ]$cis_trans <- NA
  return(coloc_out)
}

#' GTEXv8 annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
GTEXv8_annotation <- function(annotation_file, coloc_out, ...) {
  annotation_df <- read.delim(annotation_file, sep="\t", header=F)
  annotation_df$gene_id <- gsub(".*gene_id (ENSG[0-9]+.[0-9]+); .*", "\\1", annotation_df$V9)
  annotation_df$gene_type <- gsub(".*gene_type ([^;]+); .*", "\\1", annotation_df$V9)
  annotation_df$gene_name <- gsub(".*gene_name ([^;]+); .*", "\\1", annotation_df$V9)
  colnames(annotation_df)[colnames(annotation_df) == "V1"] <- "chr"
  colnames(annotation_df)[colnames(annotation_df) == "V4"] <- "gene_start"
  colnames(annotation_df)[colnames(annotation_df) == "V5"] <- "gene_end"
  annotation_df$chr <- gsub("chr", "", annotation_df$chr)
  selected_cols <- c("chr", "gene_start", "gene_end",
                     "gene_id", "gene_type", "gene_name")
  annotation_df <- annotation_df[,selected_cols]
  # coloc
  if (any(grepl("Susztak", coloc_out$sumstats_2_file))) {
    prefix <- "Kidney_eQTL_"
    colnames(annotation_df) <- paste0(prefix, colnames(annotation_df))
    annotation_df[[paste0(prefix, "gene_id_no_dot")]] <- gsub("(ENSG[0-9]+).[0-9]+", "\\1", annotation_df[[paste0(prefix, "gene_id")]])
    colnames(annotation_df)[colnames(annotation_df) == paste0(prefix, "gene_id_no_dot")] <- "ensembl_gene_id"
    coloc_out$ensembl_gene_id <- gsub(".*(ENSG[0-9]+).?.*", "\\1", basename(coloc_out$sumstats_2_file))
    coloc_out$Tissue <- gsub("_ENSG.*", "", gsub("Formated", "" , gsub("_hg38.txt.gz", "", basename(coloc_out$sumstats_2_file))))
  } else {
    prefix <- "GTExV8_"
    colnames(annotation_df) <- paste0(prefix, colnames(annotation_df))
    colnames(annotation_df)[colnames(annotation_df) == "GTExV8_gene_id"] <- "ensembl_gene_id"
    coloc_out$ensembl_gene_id <- gsub(".*(ENSG[0-9]+)", "\\1", coloc_out$sumstats_2_file)
    coloc_out$Tissue <- gsub(".*all_associations-(.*).v8.*", "\\1", coloc_out$sumstats_2_file)
  }
  coloc_out <- merge(coloc_out, by.x="ensembl_gene_id", all.x=T,
                     annotation_df, by.y="ensembl_gene_id",
                     sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  # stopifnot(unique(coloc_out$CHR_var) %in% unique(coloc_out$Olink_chr))
  coloc_out$cis_trans <- NA
  coloc_out[is.na(coloc_out[[paste0(prefix, "gene_name")]]),]$cis_trans <- NA
  coloc_out[!is.na(coloc_out[[paste0(prefix, "gene_name")]]),]$cis_trans <- "trans"
  cis_condition <- (coloc_out[[paste0(prefix, "chr")]] == coloc_out$CHR_var) & ((coloc_out[[paste0(prefix, "gene_start")]] >= coloc_out$BP_START_var & coloc_out[[paste0(prefix, "gene_start")]] <= coloc_out$BP_STOP_var) | (coloc_out[[paste0(prefix, "gene_end")]] >= coloc_out$BP_START_var & coloc_out[[paste0(prefix, "gene_end")]] <= coloc_out$BP_STOP_var))
  if (any(cis_condition)) {
    coloc_out$cis_trans[cis_condition] <- "cis"
  }
  suggestive_cis_condition <- (coloc_out$CHR_var == coloc_out[[paste0(prefix, "chr")]]) & ((coloc_out[[paste0(prefix, "gene_start")]] >= coloc_out$BP_START_var-1e6 & coloc_out[[paste0(prefix, "gene_start")]] <= coloc_out$BP_STOP_var+1e6) | (coloc_out[[paste0(prefix, "gene_end")]] >= coloc_out$BP_START_var-1e6 & coloc_out[[paste0(prefix, "gene_end")]] <= coloc_out$BP_STOP_var+1e6))
  if (any(suggestive_cis_condition)) {
    coloc_out$cis_trans[(!cis_condition) & suggestive_cis_condition] <- "suggestive_cis"
  }
  return(coloc_out)
}


#' mGWAS annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
mGWAS_annotation <- function(annotation_file, coloc_out, ...) {
  annotation_df <- read.csv(annotation_file)
  colnames(annotation_df)[colnames(annotation_df) == "ID"] <- "metabolite"
  coloc_out$metabolite <- gsub(".*[0-9]+-[0-9]+-[0-9]+_(.*_.*_.*)_.*_.*_.*_.*", "\\1", coloc_out$sumstats_2_file)
  nrow_before <- nrow(coloc_out)
  coloc_out <- merge(coloc_out, by = "metabolite",
                     annotation_df[,c("metabolite", "BIOCHEMICAL")],
                     sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  if(nrow_before != nrow(coloc_out)) { stop("Merge produced different number of rows, check duplicates or missing annotations") }
  return(coloc_out)
}

#' UKB TOPMed annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
UKB_TOPMed_annotation <- function(annotation_file, coloc_out, ...) {
  annotation_df <- read.delim(annotation_file, colClasses = "character")
  coloc_out$phenocode <- gsub(".gz", "", basename(coloc_out$sumstats_2_file))
  nrow_before <- nrow(coloc_out)
  coloc_out <- merge(coloc_out, by = "phenocode",
                     annotation_df, sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  if(nrow_before != nrow(coloc_out)) { stop("Merge produced different number of rows, check duplicates or missing annotations") }
  return(coloc_out)
}

#' FinnGen r9 annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
FinnGen_r9_annotation <- function(annotation_file, coloc_out, ...) {
  annotation_df <- read.delim(annotation_file, colClasses = "character")
  coloc_out$NAME <- gsub("finngen_R9_(.*).gz", "\\1", basename(coloc_out$sumstats_2_file))
  nrow_before <- nrow(coloc_out)
  coloc_out <- merge(coloc_out, by = "NAME", all.x=T,
                     annotation_df, sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  if(nrow_before != nrow(coloc_out)) { stop("Merge produced different number of rows, check duplicates or missing annotations") }
  return(coloc_out)
}
# 22 phenotypes are not present in the annotation file
# BMI_IRN, HEIGHT_IRN, K11_ANOMALI_DENTAL_ARCH_RELATIONS1_INCLAVO, K11_CROWDI_TEETH_INCLAVO, K11_DEEP_BITE_INCLAVO, K11_DISTAL_BITE_INCLAVO, K11_EMBED_TEETH_INCLAVO, K11_EMBIMPACT_TEETH_INCLAVO, K11_ERUPTION_INCLAVO, K11_HYPO_ONLY_INCLAVO, K11_HYPOLASIA_ENAMEL, K11_IMPACTED_TEETH_INCLAVO, K11_MAJOR_ANOMALI_JAW_SIZE_INCLAVO, K11_MIH_INCLAVO, K11_MN_PROGN_INCLAVO, K11_OPEN_BITE_INCLAVO, K11_ORAL_LICHEN_PLANUS_WIDE, K11_RESORBTION, K11_RETROG_MAXILLAE_INCLAVO, K11_SCISS_BITE_INCLAVO, Q17_CLEFT_AND_CARIES_INCLAVO, WEIGHT_IRN

