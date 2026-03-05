#!/usr/bin/env Rscript
# Test end-to-end colocalization: eGFRcr (sumstats_1) vs UKB kidney volumes (sumstats_2)
# Run on cluster: Rscript test_kidney_vol_coloc.R

library(data.table)
devtools::load_all("/g/epi/data/programs/pipelines/genepicoloc/custom_scripts/genepicoloc_dev")

# --- sumstats_1: CKDGen r4 eGFR ---
sumstats_1_file <- "/g/epi/data/public_resources/CKDGen/preprocessing/Wuttke2019/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.txt_hg38.gz"
sumstats_1 <- fread(sumstats_1_file)

# Rename to standard columns
setnames(sumstats_1, c("CHR_hg38","POS_hg38","A1_hg38","A2_hg38","Name_hg38","RSID","Freq1","Effect","StdErr","P-value","n_total_sum"),
                      c("CHR","POS","A1","A2","Name","rsID","AF","BETA","SE","P","N"))

# Calculate nlog10P
sumstats_1[, nlog10P := -log10(P)]
sumstats_1 <- sumstats_1[is.finite(nlog10P) & !is.na(BETA) & !is.na(SE) & SE > 0]
message("Loaded ", nrow(sumstats_1), " variants")

# --- Get regions ---
coloc_regions_list <- get_coloc_regions(
  sumstats = sumstats_1,
  CHR_name = "CHR", POS_name = "POS", nlog10p_value_name = "nlog10P",
  nlogP_threshold = 7.30103, halfwindow = 500000
)

n_regions <- nrow(coloc_regions_list$coloc_regions_PASS)
message("Found ", n_regions, " regions")

# Limit to first 3 regions for quick test
if (n_regions > 3) {
  coloc_regions_list$coloc_regions_PASS <- coloc_regions_list$coloc_regions_PASS[1:3, ]
  message("Using first 3 regions for testing")
}
print(coloc_regions_list$coloc_regions_PASS)

# Write regions
output_dir <- "/tmp/test_kidney_vol_coloc"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sumstats_1_subset_file <- write_regions(
  coloc_regions_list,
  sumstats_name = file.path(output_dir, "eGFR_r4")
)

# --- sumstats_2: UKB kidney volumes via YAML config ---
yaml_content <- '
UKB_kidney_vol:
  description: "UKB kidney volumes (model1, BSA-adjusted)"
  query_function: tabix_UKB_kidney_vol
  type: quant
  sdY: 1
  files:
    - /g/epi/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/10Feb2024_final_volumes/maf001/model1_qnorm_tkv.bsa_chr1-22_maf001_liftOver_hg38_dedup_sorted.txt.gz
    - /g/epi/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/10Feb2024_final_volumes/maf001/model1_qnorm_cortex.bsa_chr1-22_maf001_liftOver_hg38_dedup_sorted.txt.gz
    - /g/epi/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/10Feb2024_final_volumes/maf001/model1_qnorm_medulla.bsa_chr1-22_maf001_liftOver_hg38_dedup_sorted.txt.gz
    - /g/epi/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/10Feb2024_final_volumes/maf001/model1_qnorm_hilus.bsa_chr1-22_maf001_liftOver_hg38_dedup_sorted.txt.gz
  expected_count: 4
'
config_file <- file.path(output_dir, "kidney_vol_config.yaml")
writeLines(yaml_content, config_file)

args_df <- create_args_df_from_config(config_file)
print(args_df)

# --- Run colocalization ---
sumstats_1_args <- list(
  coloc_regions_PASS = coloc_regions_list$coloc_regions_PASS,
  sumstats_1_function = "retrieve_sumstats_tabix",
  sumstats_1_file = sumstats_1_subset_file,
  sumstats_1_type = "quant",
  sumstats_1_sdY = NA
)

message("\nRunning colocalization: 3 regions x 4 kidney volumes...")
genepicoloc_wrapper(
  dir_out = output_dir,
  sumstats_1_args = sumstats_1_args,
  args_df = args_df,
  mc_cores = 4
)

# --- Check results ---
message("\nResults:")
result_files <- list.files(output_dir, pattern = "\\.RDS$", recursive = TRUE)
message("  Output files: ", length(result_files))
for (f in result_files) message("  ", f)

message("\nDone! Output: ", output_dir)
