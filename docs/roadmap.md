# genepicoloc atlas roadmap

Living roadmap for the public genepicoloc package + Shiny atlas app.
Companion design docs: [`convergence-panel.md`](convergence-panel.md),
[`multi-ancestry.md`](multi-ancestry.md).

## Status at a glance

| Area                               | State |
|------------------------------------|-------|
| CKDGen r4 kidney traits            | shipped |
| UKB MRI kidney volumes             | shipped |
| GCKD urine metabolomics            | shipped |
| MVP kidney multi-ancestry          | in production run (test set in app) |
| Convergence tab (Region View)      | shipped (replaces old Network + RAP tabs) |
| Multi-ancestry overlay RAPs        | shipped |
| Podman local deploy                | shipped (`inst/shiny/deploy.sh`) |
| Remote VPS deploy (`--remote`)     | shipped (bundle to Nextcloud -> scp) |
| In-app LLM chat                    | removed (Datenschutz; preserved on branch) |

## sumstats_2 source tiers

### Tier 1: Core (shipped on cluster, in atlas)
| Study              | Type        | Traits          | Ancestry | Notes |
|--------------------|-------------|-----------------|----------|-------|
| MVP_R4_EUR         | PheWAS      | ~2,000          | EUR      | ancestry-split |
| UKB_PPP_EUR        | pQTL        | 2,940 proteins  | EUR      | |
| eQTLGen            | eQTL        | 19,250          | EUR      | |
| GTEx V8            | eQTL        | 49 tissues      | EUR      | 1,127 datasets |
| Kidney_eQTL        | eQTL        | 5               | EUR      | Susztak |
| CKDGen r4          | Kidney GWAS | 6               | EUR      | |
| UKB_kidney_vol     | Imaging     | 4 MRI volumes   | EUR      | |
| GCKD mGWAS (u/p)   | metabolites | 2,706           | EUR      | |
| Icelanders pGWAS   | pQTL        | 4,907           | EUR      | deCODE |
| FinnGen r9         | PheWAS      | 2,272           | FIN      | |
| UKB_TOPMed         | PheWAS      | 1,419           | EUR      | |

### Tier 2: Multi-ancestry (running now)
| Study              | Type    | Ancestry          | State |
|--------------------|---------|-------------------|-------|
| MVP_R4_{AFR,AMR,EAS,EUR,META} | PheWAS  | 5 ancestries | **Running on cluster** (full genome, two-stage orchestrator per combo; chr16 test bundle already live in the app) |
| UKB_PPP non-EUR    | pQTL    | AFR, CSA, EAS, MID | on cluster, not yet in atlas |
| Pan-UKB            | PheWAS  | 6 ancestries      | 7.2 TB preprocessed on cluster, tabix-ready; not yet wrapped |
| BioBank Japan      | PheWAS  | EAS               | to download |

### Tier 3: Deferred / watch list
- CKDGen R5 (multi-ancestry kidney) - awaiting public release
- KidneyGenAfrica (AFR eGFR) - awaiting data deposition
- eQTLGen Phase 2 (large EUR) - watching for release
- ARIC pQTL (AFR + EUR) - Synapse
- BBJ non-kidney PheWAS

## Release milestones

### Shipped
- [x] Landing page with category cards (HTML/CSS/JS inside Shiny, no TS/React)
- [x] Interactive mini-Manhattan with gene labels + click-to-select
- [x] Multi-category atlas layout (`<category>/coloc/`, `<category>/regional/`)
- [x] Bundled slim regional RDS format (one file per region, p < 1e-3 only)
- [x] Convergence tab: category cards + trait tiles + stacked RAPs + gene track
- [x] Gene info panel (NCBI, Reactome, HPO, HGNC) via bundled annotation TSV
- [x] CKDGen r4 atlas (6 kidney traits, ~259 regions)
- [x] UKB kidney MRI volumes atlas (4 traits)
- [x] GCKD urine metabolomics atlas (auto-populated, empty-metabolite filter)
- [x] Bundled slim annotation DB (`inst/shiny/extdata/shiny_annotation_db.RDS`)
- [x] Podman dev / local / remote deploy (`deploy.sh`)
- [x] MVP_kidney virtual multi-ancestry studies: consensus Manhattan,
      overlay RAPs, consensus trait dropdown, "no coloc." faded traces
      (see [`multi-ancestry.md`](multi-ancestry.md))

### In progress
- [ ] Full MVP_kidney run across all 16 traits x 5 ancestries on cluster
      (chr16 test bundle already validated in the app)
- [ ] Extract + sync MVP_kidney atlas (genome-wide) to Nextcloud mirror
- [ ] Wire Pan-UKB into `create_args_df()` and atlas (`tabix_PanUKB` +
      `format_PanUKB` helpers; metadata from `phenotype_manifest.csv`;
      per-ancestry `phenotype_qc_*` PASS filter; ln(P) -> nlog10P conversion)

### Next
- [ ] Trait View tab multi-ancestry support (currently single-ancestry)
- [ ] Vignette `vignettes/genepicoloc-workflow.qmd` filled with demo data
      and end-to-end walkthrough
- [ ] Demo dataset bundled in `inst/extdata/` (3 regions: GATM, SLC34A1, UMOD)
- [ ] Benchmark + comparison section for publication
- [ ] Tile grid for traits on Convergence (cleaner than old network)
- [ ] Export current view as PNG/PDF (base RAP + selected trait + gene track)

### Later
- [ ] BioBank Japan (EAS)
- [ ] ARIC pQTL (AFR + EUR)
- [ ] UKB_PPP non-EUR
- [ ] Cross-region comparison (two regions side by side)
- [ ] SLURM array-job mode for very large wrapper runs (see HPC section in
      `CLAUDE.md`)

## Deployment

- **Local dev:** `cd inst/shiny && ./deploy.sh --dev` (bind-mount code,
  fast restart)
- **Local full build:** `./deploy.sh --local`
- **Remote (VPS):** `./deploy.sh --remote` bundles the app (+ atlas slim)
  to `~/Work/bioinfo/nextcloud/Downloads/epi-vps/apps/genepicoloc/` which
  syncs to the VPS; then `~/epi-vps/apps/genepicoloc/deploy-remote.sh` on
  the VPS rebuilds the container.

## Tech stack

- R Shiny server
- HTML/CSS/JS (no TypeScript/React) for the landing page
- plotly for Manhattans and RAPs
- Posit PPM precompiled binaries in Podman (Docker-compatible)
- Okabe-Ito palette for ancestry overlays (see `config.R::ANCESTRY_COLORS`)
