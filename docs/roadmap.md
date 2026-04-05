# Shiny App Roadmap

## Phase 1: Landing page + RAPs for CKDGen traits (current)

- Landing page with category cards (pure HTML/CSS/JS in Shiny)
  - Kidney Diseases: eGFR, BUN, UACR, urate, gout, MA
  - Kidney MRI: TKV, cortex, medulla, hilus (UKB volumes)
  - Urine Metabolomics: GCKD VAE modules (ME19, ME30)
- Click category -> shows available traits -> click trait -> network/region view
- Regional association plots (RAPs) for all 6 CKDGen traits
- Pull sumstats from cluster for RAP data

## Phase 2: MRI volumes as trait_1

- Run genepicoloc with UKB kidney MRI volumes as trait_1 on cluster
- Add MRI results to atlas data
- Enable MRI category on landing page

## Phase 3: Urine metabolomics

- GCKD VAE analysis (ME19/ME30) results from cluster
- Add metabolomics results to atlas data
- Enable metabolomics category on landing page

## Phase 4: Cross-module links

- Trait_1 vs trait_1 colocalization (if needed)
- Cross-category network views
- Shared region highlighting across categories

## Deployment

- Podman container with mounted data volume
- deploy.sh for local dev and full builds
- Remote deploy to VPS (to be added)

## Tech stack

- R Shiny (server)
- HTML/CSS/JS for landing page (no TypeScript/React needed)
- visNetwork for network plots
- plotly for regional association plots
- Posit PPM precompiled binaries in Docker/Podman
