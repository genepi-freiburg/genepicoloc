# R package to facilitate genetic colocalization analysis.

# Installation
```
system("git clone https://github.com/genepi-freiburg/genepicoloc.git")
devtools::load_all("genepicoloc")
```

# Showcases

## 1. Identify significant regions
In this example, we will use `get_coloc_regions()` function to define significant regions for a given summary statistics.  We will use an example eGFR summary statistics file from the CKDGen consortium.  
The code was tested Linux (Ubuntu and Debian GNU) but it should also work with other platforms. For other platforms, please adjust the download command ("wget") or download the file manually.

```
library(data.table) # not mandatory, baseR functions can be used instead of `fread` below
# We use eGFR sumstats from CKDGen website
eGFR_link <- "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.txt.gz"
sumstats_name <- "eGFR_sumstats.gz"
# download sumstats
system(paste0("wget ", eGFR_link, " -O ", sumstats_name))
# read sumstats into R
sumstats <- fread(sumstats_name, fill = T)
# convert p-value to -log10
sumstats[["nlog10P"]] <- -log10(sumstats[["P-value"]])
# read_sumstats: format sumstats
sumstats_1 <- read_sumstats(sumstats = sumstats,
                            Name = "RSID",
                            CHR = "Chr",
                            POS = "Pos_b37",
                            A1 = "Allele1",
                            A2 = "Allele2",
                            BETA = "Effect",
                            SE = "StdErr",
                            nlog10p_value = "nlog10P",
                            AF = "Freq1",
                            N = "n_total_sum")
# find significant regions
## nlogP_threshold is set to -log10(5e-100) for test purposes, usually it should be -log10(5e-8)
## halfwindow is set to 500000 meaning that the window size will be 1 megabase
coloc_regions_list <- get_coloc_regions(sumstats = sumstats_1,
                                        nlogP_threshold = -log10(5e-100),
                                        halfwindow = 500000)
# inspect output
head(coloc_regions_list[["coloc_regions"]]) # identified regions
head(coloc_regions_list[["coloc_regions_PASS"]]) # identified regions excluding the merged ones
head(coloc_regions_list[["regions_log"]]) # function log
head(coloc_regions_list[["sumstats_filt"]]) # input sumstats subset to identified regions
dim(coloc_regions_list[["sumstats_filt"]]) # input sumstats subset to identified regions
```


## Supported studies

### Molecular traits

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

### Biobanks
- UKB-TOPMed <https://pheweb.org/UKB-TOPMed/about>
  - 1419 clinical outcomes (diseases)
- FinnGen release 9 <https://r9.finngen.fi/>
  - 2272 clinical outcomes (diseases)
