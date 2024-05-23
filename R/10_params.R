list_to_create_args_list <- list(
  Kidney_eQTL =
    list(EXPERIMENT = "Kidney_eQTL",
         sumstats_path = NULL,
         files = c("genepicoloc/data/Kidney_eQTL.txt.gz"),
         sumstats_pattern = "gz$",
         sumstats_function = "query_kidney_eQTL",
         sumstats_type = "quant",
         sumstats_sdY = 1,
         do_annotate = T,
         annotation_function = "transcriptomics_annotation",
         annotation_function_args = data.frame(annotation_file = "genepicoloc/data/gencode.v26.GRCh38.genes.gtf_genes_format.txt")
    ),
  UKB_PPP_EUR =
    list(EXPERIMENT = "UKB_PPP_EUR",
         sumstats_path = NULL,
         files = c("genepicoloc/data/UMOD_P07911_OID20237_v1_Cardiometabolic.txt.gz"),
         sumstats_pattern = "gz$",
         sumstats_function = "query_UKB_PPP_EUR",
         sumstats_type = "quant",
         sumstats_sdY = NA,
         do_annotate = T,
         annotation_function = "proteomics_annotation",
         annotation_function_args = data.frame(annotation_file = "genepicoloc/data/olink_protein_map_3k_v1.tsv")
         ),
  FinnGen_r9 =
    list(EXPERIMENT = "FinnGen_r9",
         sumstats_path = NULL,
         files = c("genepicoloc/data/finngen_R9_N14_CHRONKIDNEYDIS.gz"),
         sumstats_pattern = "gz$",
         sumstats_function = "query_finngen_GWAS",
         sumstats_type = "cc",
         sumstats_sdY = NA,
         do_annotate = T,
         annotation_function = "standard_annotation",
         annotation_function_args = list(annotation_file = "genepicoloc/data/endpoints.tsv")
    )
)

