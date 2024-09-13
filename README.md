genepicoloc
==============

***R package to facilitate genetic colocalization analysis.***

**Authors:** Oleg Borisov, Zulema Rodriguez Hernandez, Sara Monteiro Martins,
Inga Steinbrenner, Pascal Schlosser, et al.  
**Affiliation:** *Institute of Genetic Epidemiology, Medical Center - University of Freiburg, Germany*  
**Contact:** *oleg.borisov [at] uniklinik-freiburg [dot] de*  

Source code: https://github.com/genepi-freiburg/genepicoloc

# Overview

## Maintenance in progress
Maintenance is in progress as of September 2024, some functionalities can be limited or not available. Please check for updates in October 2024.

tl;dr - please go directly to a [typical use case](#typical-use-case) below (ideally after checking [system requirements](#system-requirements)).

Genetic colocalization is a statistical approach used to assess whether two complex traits (e.g., chronic kidney disease and type II diabetes, or chronic kidney disease and a molecular trait such as gene expression) share a causal genetic variant. As input, it uses summary statistics - the results of genome-wide association studies (GWASs). Establishing shared genetic signals can reveal common biological mechanisms underlying such complex traits and provide insights into their etiology.

A frequently used tool for genetic colocalization is [coloc](https://chr1swallace.github.io/coloc/) which *genepicoloc* is based on. In particular, it uses the enumeration method implemented in the "coloc.abf" function. This is a Bayesian approach which calculates posterior probabilities (PPs) for colocalization between each pair of traits and tests whether they show evidence of colocalization driven by the same genetic variant (corresponds to the high PP of Hypothesis 4, H4).

In essence, *genepicoloc* colocalizes input summary statistics (sumstats_1) against a set of available traits and phenotypes (sumstats_2) which include molecular traits (i.e., transcriptomics, proteomics, metabolomics, and other omics data) and diseases (clinical outcomes from phenome-wide association studies).

For more details on the method, please refer to the following resources:  
- Genome-wide association studies: [methods primer](https://www.nature.com/articles/s43586-021-00056-9)  
- Colocalization: [coloc website and references](https://chr1swallace.github.io/coloc/)  


# Installation

## System requirements
- R: tested with version 4, should be compatible with earlier versions as well (see https://www.r-project.org/)  
- *coloc* R package (see https://chr1swallace.github.io/coloc/)  
- *devtools* R package (see https://cran.r-project.org/web/packages/devtools/index.html)  
  - Usually can be installed with `install.packages("devtools")`  
- *tabix* (see http://www.htslib.org/doc/tabix.html)  
- *git* (see https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)  

The code was tested on Linux systems (Ubuntu, Debian) but should work on other systems as well.

To install *genepicoloc*, please clone the github repo first (all the following code can be run using R).

```
system("git clone https://github.com/genepi-freiburg/genepicoloc.git")
```

Next, load *genepicoloc* R using *devtools* (if the latter is available or can be installed).
```
devtools::load_all("genepicoloc")
```
Alternatively (if *devtools* is not available), simply source the scripts in "R" subfolder.
```
sapply(list.files("genepicoloc/R", full.names = T), source)
```

Next, we need to load coloc package and ensure that tabix works:
```
library(coloc)
# Test tabix installation
tabix_test <- suppressWarnings(system("tabix --version", intern=F, ignore.stdout=T, ignore.stderr=T))
if (tabix_test != 0) {
  stop("tabix command produced an error, please check that tabix is installed and accessible")
} else {"tabix is found and can be accessed successfully."}
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
3. Save and index the data (`save_coloc_regions()`).

```
# declare output name variable
sumstats_name <- "eGFR_sumstats"
# read input summary statistics into memory
sumstats <- read.csv(gzfile("genepicoloc/data/eGFR_sumstats_subset.csv.gz"))
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
coloc_regions_list <- get_coloc_regions(sumstats = sumstats_1,
                                        nlogP_threshold = -log10(5e-8),
                                        halfwindow = 500000)
```

While running, the `get_coloc_regions()` function outputs logs with the identified window. There is only 1-megabase window containing significant variants.   

To test the functionality, we can modify the parameters as following:
```
coloc_regions_list_test <- suppressMessages(get_coloc_regions(sumstats = sumstats_1,
                                        nlogP_threshold = -log10(1e-5),
                                        halfwindow = 20000))
print(coloc_regions_list_test$regions_log)
```
As we can see, four "significant" regions were identified, and two of them were merged to the first one as they were located closer than a half window to the border of the latter.  

Finally, we save the obtained results using `save_coloc_regions()` and the `sumstats_name` variable that we declared in the beginning (e.g., "eGFR_sumstats").

```
save_coloc_regions(coloc_regions_list, sumstats_name)
```

After running this command, you will find several new files in the working directory (assuming sumstats_name="eGFR_sumstats"):  
- `eGFR_sumstats_subset.tsv.gz` (formatted sumstats file) and its index `eGFR_sumstats_subset.tsv.gz.tbi` (created with tabix),  
- `eGFR_sumstats_coloc_regions_PASS.tsv` - file listing significant regions that will be used for colocalization (and its extended version with leading variants for each region `eGFR_sumstats_coloc_regions.tsv`)  
- `eGFR_sumstats_log.txt` - log file.  

After these preprocessed files ("sumstats_1") have been created, we are ready to start actual colocalization analysis.

## II. Run colocalization analysis

We are going to colocalize eGFR summary statistics against three other traits ("sumstats_2"):  
- Transcriptomics data from [the Human Kidney eQTL Atlas](https://susztaklab.com/Kidney_eQTL/)  
- Proteomics data (UMOD protein) from the [UK Biobank Pharma Proteomics Project](https://metabolomips.org/ukbbpgwas/)  
- GWAS of chronic kidney disease from the [FinnGen PheWAS, release 9](https://r9.finngen.fi/pheno/N14_CHRONKIDNEYDIS)  

We start by defining a data.frame with parameters (each line correspond to a single colocalization analysis).
```
sumstats_1_args <- data.frame(sumstats_1_file = paste0(sumstats_name, "_subset.tsv.gz"),
                              sumstats_1_function = "query_sumstats_1",
                              sumstats_1_type = "quant",
                              sumstats_1_sdY = NA)
sumstats_1_args <- data.frame(sumstats_1_args, coloc_regions_list[["coloc_regions_PASS"]])
list_of_args <- lapply(list_to_create_args_list, function(x) {
  do.call(create_coloc_params_df, 
          c(x, list(sumstats_1_args = sumstats_1_args)))
})
```

To run colocalization analysis we use a wrapper function that takes care of all data wrangling.

```
coloc_out <- suppressWarnings(Map(parallel_wrapper, list_of_args, debug_mode=T))
```

Finally, we summarize the results in the "output" folder.

```
summarize_coloc(selected_studies=selected_studies,
                output_folder = "output",
                remove_dirname = F,
                do_summary=F, do_xlsx=F)
```

As you see in the output files, our significant locus (around *UMOD* gene) in the eGFR summary statistics colocalized with:  
- *UMOD* gene expression in kidneys,  
- UMOD protein levels in plasma,  
- Chronic kidney disease (FinnGen).


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

### List of supported studies

#### Molecular traits  

**Transcriptome**
- GTEx V8 <https://www.gtexportal.org/home/>
  - 49 tissues, 39867 unique transcripts (phenotypes)
- Kidney eQTL <https://susztaklab.com/Kidney_eQTL/>
  - 5 summary statistics, 14694 unique transcripts (phenotypes)

**Proteome**

- UKB PPP <https://www.nature.com/articles/s41586-023-06592-6>
  - 2940 unique protein analytes
- Icelanders pGWAS <https://pubmed.ncbi.nlm.nih.gov/34857953/>
  - 4909 unique protein analytes
- Atherosclerosis Risk in Communities (ARIC) study <http://nilanjanchatterjeelab.org/pwas/>
  - 4657 unique protein analytes

#### Biobanks
- UKB-TOPMed <https://pheweb.org/UKB-TOPMed/about>
  - 1419 clinical outcomes (diseases)
- FinnGen release 9 <https://r9.finngen.fi/>
  - 2272 clinical outcomes (diseases)
