# Test script for cluster - run interactively in R
# Step 1: Load package functions
devtools::load_all("/g/epi/data/programs/pipelines/genepicoloc/custom_scripts/genepicoloc_dev")

# Step 2: Verify new functions exist
exists("tabix_UKB_kidney_vol", mode = "function")
exists("format_UKB_kidney_vol", mode = "function")
exists("create_args_df_from_config", mode = "function")
exists("list_available_studies", mode = "function")
exists("load_studies_config", mode = "function")

# Step 3: Create a test YAML config with real kidney volume paths
# Mount info: real = /dsk/gedata, symlink = /g/epi/data, old /data is gone
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
writeLines(yaml_content, "/tmp/test_config.yaml")

# Step 4: Test config loading - should find 4 files, no warnings
list_available_studies("/tmp/test_config.yaml")

# Step 5: Test args_df creation - should show no "files not found" warning
args_df <- create_args_df_from_config("/tmp/test_config.yaml")
str(args_df)
print(args_df)
