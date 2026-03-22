# genepicoloc

**Scalable genetic colocalization analysis across molecular and clinical trait datasets**

**Authors:** Oleg Borisov, Zulema Rodriguez Hernandez, Sara Monteiro Martins,
Inga Steinbrenner, Pascal Schlosser, et al.
**Affiliation:** Institute of Epidemiology and Prevention, Faculty of Medicine and Medical Center, University of Freiburg, Germany
**Contact:** oleg.borisov [at] uniklinik-freiburg [dot] de

## Features

- **Standalone coloc implementation** - built-in approximate Bayes factor (ABF) colocalization with zero external dependencies
- **Preprocessing tools** - liftOver (hg19->hg38), variant naming (rsID + chr:pos:ref:alt), p-value underflow handling
- **Scalable pipeline** - parallel processing with tar-based storage for millions of colocalization tests
- **Public data integration** - eQTL Catalogue API, GWAS Catalog harmonised files
- **Flexible input** - YAML-based study configuration or programmatic setup

## Installation

```r
remotes::install_github("genepi-freiburg/genepicoloc")
```

**System requirements:** `bgzip`, `tabix` (from [htslib](http://www.htslib.org/)), and `tar`.

Optional for preprocessing: UCSC [liftOver](https://hgdownload.soe.ucsc.edu/admin/exe/) binary and chain file.

## Quick Start

```r
library(genepicoloc)
library(data.table)

# 1. Load and format your GWAS (must be GRCh38)
sumstats_1 <- fread("my_gwas_hg38.tsv.gz")

# 2. Identify significant regions
coloc_regions_list <- get_coloc_regions(sumstats = sumstats_1)

# 3. Write indexed regions
sumstats_file <- write_regions(coloc_regions_list, sumstats_name = "output/my_gwas")

# 4. Define target datasets via YAML config
args_df <- create_args_df_from_config("studies.yaml")

# 5. Run colocalization
genepicoloc_wrapper(
  dir_out = "output/coloc_results",
  sumstats_1_args = list(
    coloc_regions_PASS = coloc_regions_list$coloc_regions_PASS,
    sumstats_1_function = "retrieve_sumstats_tabix",
    sumstats_1_file = sumstats_file,
    sumstats_1_type = "quant",
    sumstats_1_sdY = NA
  ),
  args_df = args_df,
  mc_cores = 10
)
```

### Standalone Colocalization

For simple pairwise tests:

```r
result <- coloc_abf(
  dataset1 = list(beta = d1$BETA, varbeta = d1$SE^2, snp = d1$Name, type = "quant", sdY = 1),
  dataset2 = list(beta = d2$BETA, varbeta = d2$SE^2, snp = d2$Name, type = "quant", sdY = 1)
)
result$summary["PP.H4.abf"]  # posterior probability of shared causal variant
```

## Example: CKDGen eGFR Kidney Atlas

Using CKDGen round 4 eGFR GWAS (Wuttke et al. 2019) as sumstats_1 against 11 molecular and clinical datasets:

- **259 genomic regions** tested against **15,930 datasets**
- **4.1 million colocalization tests** completed in ~17 hours
- **6,439 high-confidence colocalizations** (PP.H4 > 0.8)

Top results by study:

| Study | Colocalizations | Type |
|-------|-----------------|------|
| MVP EUR | 2,504 | PheWAS |
| UKB PPP | 1,639 | Plasma proteomics |
| Icelanders pGWAS | 787 | Plasma proteomics |
| Kidney eQTL | 498 | Kidney gene expression |
| CKDGen r4 | 333 | Kidney function traits |
| FinnGen r9 | 281 | Clinical endpoints |
| GCKD metabolomics | 173 | Plasma + urine metabolites |
| eQTLGen | 89 | Blood eQTL |
| UKB kidney MRI | 66 | Kidney volumes |

See the [vignette](vignettes/genepicoloc-workflow.qmd) for the full workflow.

## Preprocessing

If your GWAS is in hg19:

```r
# LiftOver hg19 -> hg38
sumstats_hg38 <- genepi_liftover(
  sumstats = my_gwas, CHR_name = "Chr", POS_name = "Pos_b37",
  A1_name = "A1", A2_name = "A2",
  liftOver_bin = "liftOver", liftOver_chain = "hg19ToHg38.over.chain.gz"
)

# Add rsID and standardized variant names
sumstats_named <- name_by_position(
  sumstats = sumstats_hg38,
  CHR_name = "CHR_hg38", POS_name = "POS_hg38",
  A1_name = "A1_hg38", A2_name = "A2_hg38",
  dbSNP_dir = "/path/to/ensembl_vcf/"
)
```

## Required Columns

| Column | Description | Required |
|--------|-------------|----------|
| Name | chr:pos:ref:alt | yes |
| rsID | dbSNP identifier | optional |
| CHR, POS | Chromosome and position (GRCh38) | yes |
| A1, A2 | Effect and other allele | yes |
| BETA, SE | Effect size and standard error | yes |
| nlog10P | -log10(p-value) | yes |
| AF, N | Allele frequency, sample size | quantitative traits only |

## Citation

If you use genepicoloc, please cite:

> Borisov et al. (2026) genepicoloc: Scalable genetic colocalization analysis
> for kidney genomics. *In preparation.*

## License

MIT
