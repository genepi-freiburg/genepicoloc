# Atlas Colocalization Strategy

## Goal
Build a kidney genomics colocalization atlas with multi-ancestry support.
Three trait_1 categories: kidney diseases (CKDGen r4), kidney MRI (UKB), urine metabolomics (GCKD VAE).

## sumstats_2 Selection

### Tier 1: Core (all trait_1s)
| Study | Type | Traits | Ancestry | Status |
|-------|------|--------|----------|--------|
| MVP_R4_EUR | PheWAS | 2,068 | EUR | On cluster, ready |
| UKB_PPP_EUR | pQTL | 2,940 proteins | EUR | On cluster, ready |
| eQTLGen | eQTL | 19,250 transcripts | EUR | On cluster, ready |
| CKDGen_r4 | Kidney traits | eGFR, BUN, UACR, urate, gout, MA | EUR | On cluster, ready |
| UKB_kidney_vol | Imaging | 4 MRI volumes | EUR | On cluster, ready |

### Tier 2: Multi-ancestry
| Study | Type | Traits | Ancestry | Status |
|-------|------|--------|----------|--------|
| MVP_R4_AFR | PheWAS | ~2,000 | AFR | On cluster, ready |
| MVP_R4_EAS | PheWAS | ~2,000 | EAS | On cluster, ready |
| MVP_R4_AMR | PheWAS | ~2,000 | AMR | On cluster, ready |
| UKB_PPP non-EUR | pQTL | 2,940 | AFR/CSA/EAS/MID | Small N but publishable |
| Pan-UKB | PheWAS | 7,228 x 6 ancestries | Multi | On cluster, 7.2TB processed, tabix-ready |
| BioBank Japan | PheWAS | 220 | EAS | To download |

### Tier 3: Deferred
| Study | Reason |
|-------|--------|
| GTEx eQTL | Too many datasets, low kidney relevance |
| Kidney eQTL | p<1e-5 limitation, unreliable for coloc |
| Icelanders pGWAS | EUR only, overlaps UKB_PPP |
| FinnGen | EUR only, overlaps MVP/Pan-UKB |
| UKB_TOPMed | EUR only, replaced by Pan-UKB |
| GCKD mGWAS | Internal, low N (4,850) |

### Watch list
| Study | Note |
|-------|------|
| CKDGen R5 | Multi-ancestry kidney, await public data release |
| KidneyGenAfrica | AFR kidney eGFR, await data deposition |
| eQTLGen Phase 2 | Massive power upgrade, EUR only |
| ARIC pQTL | AA + EA proteomics on Synapse |

## Execution Plan

### Phase 1: Quick wins (Tier 1 only)
Run MRI and uMet against: eQTLGen, UKB_PPP_EUR, MVP_R4_EUR, CKDGen_r4, UKB_kidney_vol
- Reuse 02_run_coloc.R pattern from CKDGen_r4
- ~5 studies x ~10k datasets = manageable

### Phase 2: Multi-ancestry
Add MVP_R4_AFR/EAS/AMR to all trait_1s (MRI, uMet, and existing CKDGen)
- Same pipeline, just different selected_studies parameter

### Phase 3: Pan-UKB
Add Pan-UKB across all trait_1s and ancestries
- Needs tabix_PanUKB() and format_PanUKB() in genepicoloc or genepicoloc_internal
- 7.2 TB on cluster at /g/epi/data/public_resources/UKBB_public/09_PanUKB/sumstats_hg38/
- phenotype_manifest.csv has metadata for all traits
- Files: {category}-{phenocode}-{sex}-{modifier}_hg38_sorted_hg38.tsv.gz + .tbi

### Phase 4: Additional non-EUR resources
- BioBank Japan (EAS)
- ARIC pQTL (AFR)
- UKB_PPP non-EUR

## Pan-UKB Details
- Path: /g/epi/data/public_resources/UKBB_public/09_PanUKB/sumstats_hg38/
- Format: tabix-indexed TSV, hg38, all ancestries in one file
- P-values: stored as ln(P), need conversion to nlog10P = -ln(P) / log(10)
- Phenotype manifest: phenotype_manifest.csv (trait metadata, sample sizes, QC flags)
- QC: phenotype_qc_{ancestry} column - filter to "PASS" for each ancestry
- 7,228 phenotypes across 6 ancestries
- Processed by rodriguz, ready for genepicoloc integration

## Shiny App Categories
1. Kidney Disease Traits: eGFR, BUN, UACR, urate, gout, MA (CKDGen r4)
2. Kidney MRI Volumes: TKV, cortex, medulla, hilus (UKB)
3. Urine Metabolomics: ME19, ME30 (GCKD VAE)
