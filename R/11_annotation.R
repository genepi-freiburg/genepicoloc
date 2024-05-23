#' transcriptomics annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
transcriptomics_annotation <- function(study, annotation_file,
                                       sumstats_file = "sumstats_2_file",
                                       coloc_out,
                                       CHR_var = "CHR_var", BP_START_var = "BP_START_var",
                                       BP_STOP_var = "BP_STOP_var") {
  annotation_df <- read.delim(annotation_file)
  colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
  if (study == "GTEXv8") {
    merge_column <- paste0(study, "_gene_id")
    coloc_out[[merge_column]] <- gsub(".*(ENSG[0-9]+)", "\\1", basename(coloc_out[[sumstats_file]]))
    coloc_out[[paste0(study, "_Tissue")]] <- gsub(".*all_associations-(.*).v8.*", "\\1", basename(coloc_out[[sumstats_file]]))
  }
  if (study %in% c("Kidney_eQTL", "eQTLGen")) {
    merge_column <- paste0(study, "_gene_id_no_dot")
    annotation_df[[merge_column]] <- gsub("(ENSG[0-9]+).?.*", "\\1", annotation_df[[paste0(study, "_gene_id")]])
    coloc_out[[merge_column]] <- gsub(".*(ENSG[0-9]+).?.*", "\\1", basename(coloc_out[[sumstats_file]]))
    coloc_out[[paste0(study, "_Tissue")]] <- gsub("_ENSG.*", "", gsub("Formated", "" , gsub("_hg38.txt.gz", "", basename(coloc_out[[sumstats_file]]))))
  }
  coloc_out <- merge(coloc_out, annotation_df, by = merge_column, all.x=T,
                     sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  coloc_out[[paste0(study, "_cis_trans")]] <- cis_trans_annotation(region_CHR_vec = coloc_out[[CHR_var]],
                                                                   region_BP_START_vec = coloc_out[[BP_START_var]],
                                                                   region_BP_STOP_vec = coloc_out[[BP_STOP_var]],
                                                                   gene_chr_vec = coloc_out[[paste0(study, "_chr")]],
                                                                   gene_start_vec = coloc_out[[paste0(study, "_gene_start")]],
                                                                   suggestive_window = 1e6)
  return(coloc_out)
}

#' proteomics annotation
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
proteomics_annotation <- function(study, annotation_file,
                                  sumstats_file = "sumstats_2_file",
                                  coloc_out, CHR_var = "CHR_var", BP_START_var = "BP_START_var",
                                  BP_STOP_var = "BP_STOP_var") {
  if (study %in% c("UKB_PPP_EUR", "GCKD_pGWAS")) {
    merge_column <- paste0(study, "_OlinkID")
    annotation_df <- read.delim(annotation_file, sep="\t")
    colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
    annotation_df[[paste0(study, "_multiple_genes_per_OID")]] <- 0
    selected_cols <- paste0(study, "_", c("OlinkID", "olink_target_fullname", "UniProt", "Assay",
                                          "HGNC.symbol", "ensembl_id", "chr", "gene_start", "gene_end",
                                          "multiple_genes_per_OID"))
    annotation_df <- annotation_df[,selected_cols]
    if (sumstats_file == "sumstats_2_file") {
      coloc_out[[merge_column]] <- gsub(".*(OID[0-9]+).*", "\\1", coloc_out[[sumstats_file]])
    }
  }
  if (study == "Icelanders_pGWAS") {
    merge_column <- paste0(study, "_SeqId")
    annotation_df <- read.csv(annotation_file)
    colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
    annotation_df[[paste0(study, "_chr")]] <- gsub("chr", "", annotation_df[[paste0(study, "_chr")]])
    coloc_out[[merge_column]] <- gsub(".*NG2021/([0-9]+_[0-9]+)_.*", "\\1", coloc_out[[sumstats_file]])
  }
  if (study == "ARIC_pGWAS") {
    merge_column <- paste0(study, "_seqid_in_sample")
    annotation_df <- read.delim(annotation_file)
    colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
    colnames(annotation_df)[colnames(annotation_df) == paste0(study, "_chromosome_name")] <- paste0(study, "_chr")
    colnames(annotation_df)[colnames(annotation_df) == paste0(study, "_transcription_start_site")] <- paste0(study, "_gene_start")
    coloc_out[[merge_column]] <- gsub(".PHENO1.glm.linear.gz", "", basename(coloc_out[[sumstats_file]]))
  }
  cis_trans_column <- paste0(study, "_cis_trans")
  gene_chr_vec_name <- paste0(study, "_chr")
  gene_start_vec_name <- paste0(study, "_gene_start")
  if (sumstats_file == "sumstats_1_file") {
    if (study != "GCKD_pGWAS") {stop("Only GCKD_pGWAS as sumstats_1 is supported so far")}
    # rename sumstats_2 columns
    colnames(coloc_out)[grep(paste0(study, "_"), colnames(coloc_out))] <- 
      paste0("sumstats_2_", colnames(coloc_out)[grep(paste0(study, "_"), colnames(coloc_out))])
    merge_column <- paste0("sumstats_1_", merge_column)
    colnames(annotation_df) <- paste0("sumstats_1_", colnames(annotation_df))
    coloc_out[[merge_column]] <- gsub(".*(OID[0-9]+).*", "\\1", coloc_out[[sumstats_file]])
    # update cis_trans columns
    cis_trans_column <- paste0("sumstats_1_", study, "_cis_trans")
    gene_chr_vec_name <- paste0("sumstats_1_", study, "_chr")
    gene_start_vec_name <- paste0("sumstats_1_", study, "_gene_start")
  }
  coloc_out <- merge(coloc_out, annotation_df, by=merge_column, all.x=T,
                     sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  coloc_out[[cis_trans_column]] <- cis_trans_annotation(region_CHR_vec = coloc_out[[CHR_var]],
                                                        region_BP_START_vec = coloc_out[[BP_START_var]],
                                                        region_BP_STOP_vec = coloc_out[[BP_STOP_var]],
                                                        gene_chr_vec = coloc_out[[gene_chr_vec_name]],
                                                        gene_start_vec = coloc_out[[gene_start_vec_name]],
                                                        suggestive_window = 1e6)
  return(coloc_out)
}


#' standard annotation
#' mGWAS, UKB TOPMed, FinnGen r9
#' @param annotation_file path to annotation file
#' @return data.frame with annotated information.
#' @examples
#' Under development
#' @export
standard_annotation <- function(study, annotation_file,
                                sumstats_file = "sumstats_2_file",
                                coloc_out) {
  if (study %in% c("mGWAS_plasma", "mGWAS_urine")) {
    merge_column <- paste0(study, "_metabolite")
    annotation_df <- read.csv(annotation_file)
    colnames(annotation_df)[colnames(annotation_df) == "ID"] <- "metabolite"
    coloc_out[[merge_column]] <- gsub(".*[0-9]+-[0-9]+-[0-9]+_(.*_.*_.*)_.*_.*_.*_.*", "\\1", coloc_out[[sumstats_file]])
    colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
    nrow_before <- nrow(coloc_out)
    coloc_out <- merge(coloc_out, by = merge_column, all.x=T,
                       annotation_df[,c(paste0(study, "_metabolite"), paste0(study, "_BIOCHEMICAL"))],
                       sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  }
  if (study == "UKB_TOPMed") {
    merge_column <- paste0(study, "_phenocode")
    annotation_df <- read.delim(annotation_file, colClasses = "character")
    coloc_out[[merge_column]] <- gsub(".gz", "", basename(coloc_out[[sumstats_file]]))
    colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
    selected_cols <- paste0(study, "_", c("phenocode", "num_cases", "num_controls",
                                          "num_samples", "phenostring", "num_peaks",
                                          "gc_lambda_hundred", "category"))
    annotation_df <- annotation_df[,selected_cols]
    nrow_before <- nrow(coloc_out)
    coloc_out <- merge(coloc_out, by = merge_column, all.x=T,
                       annotation_df, sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  }
  if (study == "FinnGen_r9") {
    merge_column <- paste0(study, "_phenocode")
    annotation_df <- read.delim(annotation_file, colClasses = "character")
    coloc_out[[merge_column]] <- gsub("finngen_R9_(.*).gz", "\\1", basename(coloc_out[[sumstats_file]]))
    colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
    nrow_before <- nrow(coloc_out)
    coloc_out <- merge(coloc_out, by = merge_column, all.x=T,
                       annotation_df, sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  }
  if (study == "CKDGen" | study == "Infections23" | study == "RCC") {
    merge_column <- paste0(study, "_File")
    annotation_df <- read.delim(annotation_file, colClasses = "character", header = T)
    coloc_out[[merge_column]] <- basename(coloc_out[[sumstats_file]])
    colnames(annotation_df) <- paste0(study, "_", colnames(annotation_df))
    nrow_before <- nrow(coloc_out)
    coloc_out <- merge(coloc_out, by = merge_column, all.x=T,
                       annotation_df, sort = FALSE)[, union(names(coloc_out), names(annotation_df))]
  }
  if(nrow_before != nrow(coloc_out)) { stop("Merge produced different number of rows, check duplicates or missing annotations") }
  return(coloc_out)
}



