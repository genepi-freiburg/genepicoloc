# genepicoloc atlas - user guide

## What is this?

A region-centric multi-omics atlas of kidney genomics. For each
colocalization region in a source GWAS (e.g. CKDGen eGFR, MVP BUN),
the atlas shows which traits colocalize with the signal and renders
their regional association plots for spatial comparison.

**Data is curated; interpretation is yours.** The app shows what
colocalizes, not what genes to prioritize or what hypotheses to form.

## Datasets on the landing page

### Featured atlases

- **CKDGen Round 4** - Kidney function and disease phenotypes
  (eGFR creatinine, BUN, UACR, urate, gout, microalbuminuria) from
  the CKDGen Round 4 meta-analysis (Wuttke et al., Nat Genet 2019 and
  later updates). European ancestry primarily.
- **MVP Kidney (multi-ancestry)** - Million Veteran Program kidney
  phenotypes colocalized across AFR, AMR, EAS, EUR and a trans-ancestry
  META meta-analysis. Each trait loads all available ancestries at once
  and overlays them in the Manhattan and the regional plots. See the
  "Multi-ancestry view" section below.

### Additional phenotypes

- **Kidney MRI Volumes** - UK Biobank kidney MRI (total kidney volume,
  cortex, medulla, hilus), BSA-adjusted, EUR.
- **Urine Metabolomics** - GCKD urine metabolome GWAS (Schlosser et al.,
  Nat Genet 2023), 1,409 metabolites.

## How to read a region

Click a trait on the landing page to enter the **Region View**. The
top half shows a mini-Manhattan across the whole atlas; click a dot
to zoom into that region.

For the selected region you get, top to bottom:

1. **Trait of interest** - the source GWAS regional plot. One scatter
   per ancestry if the study is multi-ancestry.
2. **Selected trait** - when you click a trait tile below, its
   regional association plot appears here, sharing the x-axis with
   the trait of interest so you can visually check whether the peaks
   colocalize in space as well as statistically.
3. **Genes in region** - protein-coding genes in the window. Click a
   gene to get NCBI description, HGNC family, Reactome pathways and
   HPO phenotype annotations.
4. **Colocalized traits** - category cards on the left (kidney
   diseases, PheWAS, imaging, proteins, transcripts, metabolites)
   and a tile grid of the top traits in the selected category on the
   right. Click a tile to populate the Selected trait plot above.

## Multi-ancestry view (MVP Kidney)

When the source study has multiple ancestries, the Manhattan dot
colors reflect **how many ancestries hit each locus** (1-4 coverage
gradient) and tooltips list which ones. Per-ancestry regions are
automatically merged into consensus clusters based on interval
overlap.

When you pick a region:

- Both the **trait of interest** and the **selected trait** plots
  render one colored scatter per ancestry. Colors follow the
  [Okabe-Ito palette](https://jfly.uni-koeln.de/color/): AFR orange,
  AMR green, EAS yellow, EUR blue, META pink.
- Ancestries whose regional sumstats are available but whose
  colocalization did **not** pass the PP.H4 threshold are drawn with
  reduced opacity and labelled `<ANC> (no coloc.)` in the legend.
  This is deliberate: it shows the raw signal as context without
  claiming it colocalized.
- The trait tile badge (e.g. `[AMR,EUR]`) always matches the solid
  traces. A faded trace means "same signal, different coloc status".

## Filters (left sidebar)

- **Min PP.H4** (default 0.8) - minimum posterior probability for
  shared causal variant from `coloc.abf`.
- **Min -log10(P)** - minimum sumstats_2 lead-variant strength.
- **Include studies** - toggle individual sumstats_2 studies on/off.

## MHC exclusion

The extended MHC locus (chr6:28,510,120-33,480,577 in hg38) has
extreme LD complexity that breaks `coloc.abf`'s single-causal-variant
assumption. MHC colocalizations are dropped at extraction time and
will not appear in the atlas.

## Data sources and versions

| Category              | Source                          | Reference |
|-----------------------|---------------------------------|-----------|
| Source GWAS           | CKDGen Round 4                  | PMID 31152163 |
| Source GWAS           | MVP Release 4                   | MVP consortium |
| Imaging               | UK Biobank kidney MRI           | in-house |
| Metabolites           | GCKD urine metabolomics         | PMID 37277652 |
| eQTL                  | eQTL Catalogue (GTEx, eQTLGen)  | https://www.ebi.ac.uk/eqtl/ |
| eQTL                  | Susztak kidney eQTL             | PMID 30578417 |
| pQTL                  | UKB PPP (Sun et al.)            | PMID 37794183 |
| pQTL                  | deCODE Icelanders               | PMID 34857953 |
| PheWAS                | FinnGen r9                      | https://finngen.fi/ |
| PheWAS                | UK Biobank TOPMed               | Pan-UKB |
| PheWAS                | MVP_R4 (per ancestry)           | MVP consortium |
| Kidney traits         | CKDGen_r4                       | PMID 31152163 |

## Technical notes

- Regional sumstats are bundled in a slim format (`p < 1e-3`) to keep
  the atlas under a few hundred MB. Full sumstats live on the cluster
  and can be rerun through the `genepicoloc` R package if needed.
- The Manhattan gene label is the `nearest_gene_1` from the
  genepicoloc annotator (nearest protein-coding gene to the region
  center). The hover tooltip shows up to 10 neighbouring genes so you
  can sanity-check the label against known kidney genes.
- The MVP Kidney consensus-region logic is documented at length in
  `docs/multi-ancestry.md` in the genepicoloc repository.

## Source code and citation

- Package: https://github.com/genepi-freiburg/genepicoloc
- Internal study definitions:
  https://github.com/genepi-freiburg/genepicoloc_internal
- Please cite the respective source GWAS and sumstats_2 studies as
  well as the genepicoloc package when using results from this atlas.

## Feedback

Open an issue at
https://github.com/genepi-freiburg/genepicoloc/issues.
