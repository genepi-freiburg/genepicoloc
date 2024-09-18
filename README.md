genepicoloc
==============

***R package to facilitate genetic colocalization analysis.***

**Authors:** Oleg Borisov, Zulema Rodriguez Hernandez, Sara Monteiro Martins,
Inga Steinbrenner, Pascal Schlosser, et al.  
**Affiliation:** *Institute of Genetic Epidemiology, Medical Center - University of Freiburg, Germany*  
**Contact:** *oleg.borisov [at] uniklinik-freiburg [dot] de*  

Source code: https://github.com/genepi-freiburg/genepicoloc

# ! Maintenance in progress
Maintenance is in progress as of September 2024, some functions can be limited.

# ----

# Overview

tl;dr - please go directly to a [typical use case](#typical-use-case) below (ideally after checking [system requirements](#system-requirements)).

Genetic colocalization is a statistical approach used to assess whether two complex traits (e.g., chronic kidney disease and type II diabetes, or chronic kidney disease and a molecular trait such as gene expression) share a causal genetic variant. As input, it uses the results of genome-wide association studies (GWASs) - summary statistics. Establishing shared genetic signals can reveal common biological mechanisms underlying such complex traits and provide insights into their etiology.

A frequently used tool for genetic colocalization is [coloc](https://chr1swallace.github.io/coloc/) which *genepicoloc* is based on. In particular, it uses the enumeration method implemented in the "coloc.abf" function. This is a Bayesian approach which calculates posterior probabilities (PPs) for colocalization between each pair of traits and tests whether they show evidence of colocalization driven by the same genetic variant (corresponds to the high PP of Hypothesis 4, H4).

In essence, *genepicoloc* is aimed at colocalizing input summary statistics (sumstats_1) against a set of available traits and phenotypes (sumstats_2) which include molecular traits (i.e., transcriptomics, proteomics, metabolomics, and other omics data) and diseases (clinical outcomes from phenome-wide association studies).

For more details on the method, please refer to the following resources:  
- Genome-wide association studies: [methods primer](https://www.nature.com/articles/s43586-021-00056-9)  
- Colocalization: [coloc website and references](https://chr1swallace.github.io/coloc/)  


# Installation

## System requirements
- R: tested with version 4, should be compatible with earlier versions as well (see https://www.r-project.org/)  
- *coloc* R package (see https://chr1swallace.github.io/coloc/)  
- *remotes* R package (see https://cran.r-project.org/web/packages/remotes/index.html)  
- *tabix* (see http://www.htslib.org/doc/tabix.html)  
- *git* (see https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)  

The code was tested on Linux systems (Ubuntu, Debian) but should work on other systems as well.  

To install *genepicoloc*, please run:
```
remotes::install_github("https://github.com/genepi-freiburg/genepicoloc.git")
library(genepicoloc)
```

Next, we need to load coloc package and ensure that tabix works:
```
library(coloc)
# Test tabix installation
tabix_test_fail <- suppressWarnings(system("tabix --version", intern=F, ignore.stdout=T, ignore.stderr=T))
if (tabix_test_fail) stop("Please check that tabix is installed and accessible")
```

# Typical use case

Please note our recommendations for colocalization analysis: [good genepicoloc practice](#good-genepicoloc-practice).

To illustrate the functionality of *genepicoloc*, we will perform colocalization analysis for a complex trait. As an example, we will use a summary statistics of kidney function (eGFR, estimated glomerular filtration rate). We will refer to this summary statistics as "sumstats_1". The GWAS of eGFR was performed by the CKDGen consortium as described in the [article](https://www.nature.com/articles/s41588-019-0407-x) and the source summary statistics file is available from the [CKDGen website](https://ckdgen.imbi.uni-freiburg.de/datasets/Wuttke_2019).

For illustration purpose, we use a subset of this file with chromosome 16 and a 4-megabase window around one of the significant regions (`Chr == 16 & Pos_b37 > 18e6 & Pos_b37 < 22e6`). In addition, the liftOver was performed to map coordinated to build GRCh38 (using [liftOver tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver)).

The workflow will consist of two parts:  
I. Format and identify significant regions.  
II. Run colocalization analysis.  


## I. Format and identify significant regions

We will start with three following steps:  
1. Read and format input summary statistics (`read_sumstats()`).
2. Identify significant regions (`get_coloc_regions()`).
3. Save and index the data (`write_regions()`).

```
# declare mandatory variables
sumstats_1_name <- "eGFR_sumstats"
sumstats_1_type <- "quant"
sumstats_1_sdY <- NA
output_folder <- "output"
# read input summary statistics into memory
test_file <- system.file("data/test_sumstats.RDS", package="genepicoloc")
sumstats <- readRDS(test_file)
# Function 1: format input sumstats
sumstats_1 <- read_sumstats(sumstats = sumstats,
                            Name = "Name_hg38",
                            rsID = "rs",
                            CHR = "CHR_hg38",
                            POS = "POS_hg38",
                            A1 = "A1_hg38",
                            A2 = "A2_hg38",
                            BETA = "Effect",
                            SE = "StdErr",
                            nlog10p_value = "nlog10P",
                            AF = "Freq1",
                            N = "n_total_sum")
head(sumstats_1, 2)
```

Using column names provided by the user, `read_sumstats()` automatically adapts the sumstats table for the required input format.  
- Name has the format of CHR:POS:REF:ALT (e.g., chr16:17906244:T:C).  
- rsID column is optional and can be filled with NA.  
- Note that p-value is represented as negative log10 (this is required for handling underflow issues in case of very low p-values, p < 1e-324, which is not uncommon in, e.g., large size pQTL studies).  
- For more information on formatting, please check [good genepicoloc practice](#good-genepicoloc-practice).  


Next, we will identify significant regions using `get_coloc_regions()` with the following parameters:  
- Standard GWAS p-value threshold (5e-8),  
- and 1-megabase window around the index variant with lowest p-value (i.e., halfwindow = 500000)
```
coloc_regions_list <- get_coloc_regions(sumstats = sumstats_1)
```

While running, the `get_coloc_regions()` function outputs logs with the identified window. There is only 1-megabase window containing significant variants.  

Finally, we save the obtained results using `write_regions()` and the `sumstats_name` variable that we declared in the beginning (e.g., "eGFR_sumstats"). At the same time, we will create a data.frame with arguments for colocalization.  

```
sumstats_1_file <- write_regions(coloc_regions_list, sumstats_1_name)
sumstats_1_args <- create_sumstats_1_args(sumstats_1_file=sumstats_1_file,
                                          sumstats_1_type=sumstats_1_type,
                                          sumstats_1_sdY=sumstats_1_sdY)
```

After running this command, you will find several new files in the working directory (assuming sumstats_name="eGFR_sumstats"):  
- `eGFR_sumstats_subset.tsv.gz` (formatted sumstats file) and its index `eGFR_sumstats_subset.tsv.gz.tbi` (created with tabix),  
- `eGFR_sumstats_coloc_regions_PASS.tsv` - file listing significant regions that will be used for colocalization (and its extended version with leading variants for each region `eGFR_sumstats_coloc_regions.tsv`)  
- `eGFR_sumstats_log.txt` - log file.  

After these preprocessed files ("sumstats_1") have been created, we are ready to start actual colocalization analysis.

## II. Run colocalization analysis

We are going to colocalize eGFR summary statistics against GTEx eQTLs from eQTL Catalog.  
```
eQTL_Catalogue <- make_eQTL_Catalogue_args()
args_df <- do.call(create_args_df, c(coloc_regions_list$coloc_regions_PASS,
                                     sumstats_1_args,
                                     eQTL_Catalogue))
```

To run colocalization analysis we use a wrapper function `map_over_args()` that takes care of all data wrangling. We select only two GTEx tissues for illustration purpose `c(30,49)`. `annotate_eQTL_Catalog()` will annotate the results.  

```
coloc_out <- map_over_args(args_df[c(30,49),], mc_cores=2)
annotation_file <- system.file("data/ENSG_HGNC.RDS", package="genepicoloc")
annotation_df <- readRDS(annotation_file)

datasets_eQTL_Catalogue <- get_datasets_eQTL_Catalogue()
coloc_out_annot <- annotate_eQTL_Catalog(coloc_out=coloc_out,
                                         annotation_df=annotation_df,
                                         datasets_eQTL_Catalogue=datasets_eQTL_Catalogue)
```

An example on how to create a data.frame with arguments manually.  

```
args_df_example <- create_args_df(CHR_var="16",
                                  BP_START_var=19850119,
                                  BP_STOP_var=20850119,
                                  sumstats_1_file="eGFR_sumstats_subset.tsv.gz",
                                  sumstats_1_function="query_sumstats_1",
                                  sumstats_1_type="quant",
                                  sumstats_1_sdY=NA,
                                  sumstats_2_files="QTD000356",
                                  sumstats_2_function="query_eQTL_Catalogue",
                                  sumstats_2_type="quant",
                                  sumstats_2_sdY=1)
```


# Additional information
## Good genepicoloc practice
- Genomic build. Input sumstats should have build GRCh38. If the build is different, please perform liftOver first, *genepicoloc::genepi_liftOver()* may be helpful.
- Variant naming. Input sumstats should contain variant name in the following format: CHR:POS:REF:ALT. Some summary statistics are coming from [metal software](https://genome.sph.umich.edu/wiki/METAL) and allele order may be arbitrary, i.e., allele 1 and 2 correspond to ALT and REF only in approximately 50% of cases. If that is the case, *genepicoloc::Name_by_position()* function may be helpful to find proper variant names. If liftOver needs to be performed, then genepi_liftOver() will run Name_by_position() automatically.
- The p-values. Input sumstats should contain -log10(P) column. If there is a normal "P" column, ensure there is not underflow issues (P<1e-324). If it is the case, please use genepicoloc::handle_underflow() to properly format input sumstats (use return_nlog10P=T if needed).
- Duplicated variant names. Input sumstats should contain only unique variant names, otherwise coloc will throw an error. Please remove duplicates from the input sumstats when performing QC.
- BETA and SE. Input sumstats should not contain any missing or NA values in BETA and SE column, otherwise coloc will throw an error. Please remove any missing values from the input sumstats when performing QC.
- MAF should be a numeric, strictly >0 & <1.



## Large-scale colocalization framework

A number of summary statistics are suppored by the *genepicoloc*. Due to the very large size of the input files, we are not able to upload them on github. However, we have implemented a colocalization framework at our local computing cluster in Freiburg (Germany). We are open for collaborations, feel free to reach out in case you are interested to run colocalization analysis for your GWAS summary statistics.

