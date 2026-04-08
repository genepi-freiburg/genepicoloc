# Convergence tab (Region View)

## Goal

A **region-centric multi-omics view** that answers: "What traits
colocalize with the GWAS signal in this region, grouped by evidence
type, and what does the spatial overlap look like?"

This **is** the primary region view in the current app - the old
separate Network + Regional Plot tabs are gone.

**Philosophy:** present evidence agnostically, no gene prioritization.
The atlas shows what colocalizes; the user decides what it means.

For multi-ancestry virtual studies (e.g. MVP_kidney) the same layout
renders overlaid per-ancestry traces in both RAPs - see
[`multi-ancestry.md`](multi-ancestry.md) for the overlay semantics,
consensus regions, and "no coloc." trace rule.

## Layout (top to bottom)

1. **Trait of interest** - regional association plot of the study-1 GWAS
   (e.g. eGFR / BUN) for the current region. Uses
   `plot_regional_association_interactive()`.

2. **Selected trait RAP** (conditional) - when the user clicks a trait
   tile in the drilldown, its regional plot appears directly below for
   spatial comparison. Shares x-range with the base plot. Disappears
   when no trait is selected.

3. **Genes in region** - gene track from `plot_gene_track()` in
   `R/gene_track.R`. Protein-coding genes in the window; clickable for
   the gene info panel (NCBI / Reactome / HPO / HGNC annotations).

4. **Colocalized traits** - category cards on the left, trait drilldown
   on the right.
   - 6 category cards stacked vertically, clickable `actionButton`s
     with colored stripes and counts.
   - Clicking a card selects that category (click again to close).
   - Right side shows a tile grid of up to 50 top traits in the
     category, sorted by sumstats_2 -log10P, with per-trait header,
     study badge and -log10P.
   - Clicking a tile triggers the selected-trait RAP above.

## Categories

Defined in `R/convergence.R::TRAIT_CATEGORIES`. Six top-level buckets:

| id | label | studies |
|----|-------|---------|
| `kidney_disease` | Kidney diseases & traits | CKDGen_r4, CKDGen_r5 |
| `phewas` | PheWAS (diseases & traits) | MVP_R4*, FinnGen_r9, UKB_TOPMed, PanUKB |
| `imaging` | Kidney imaging | UKB_kidney_vol |
| `proteins` | Proteins (pQTL) | UKB_PPP_EUR, Icelanders_pGWAS, GCKD_pGWAS |
| `transcripts` | Transcripts (eQTL) | eQTLGen, Kidney_eQTL, GTEXv8_eQTL |
| `metabolites` | Metabolites | GCKD_mGWAS_plasma/urine, GCKD_uMet |

Kidney diseases and kidney imaging are intentionally split from general
PheWAS because this is a kidney-focused atlas and those signals should
stand out.

## Key reactives

- `conv_region_coords()` - parses `current_region()` into chr/start/end.
- `conv_cat_counts()` - named int vector of trait counts per category
  for the current filtered region data (uses `filtered_region_data()`,
  NOT the `max_traits`-limited `filtered_data()`).
- `conv_selected_cat` - reactiveVal, currently selected category id or
  NULL.
- `conv_selected_trait` - reactiveVal, currently clicked trait's bundle
  key (or consensus key in virtual mode).
- `conv_drilldown_data()` - filtered to the selected category. For
  virtual studies it also collapses rows per consensus trait and
  attaches `ancestries` / `n_ancestries` / `consensus_key`.
- `conv_trait_sumstats()` - loads the regional RDS(s) for the selected
  trait. Returns a single data.table for normal studies, a named list
  of per-ancestry data.tables for virtual studies (see
  `multi-ancestry.md`).

## Click flow

```
User clicks category card
  -> observer sets conv_selected_cat(cat_id)
  -> conv_drilldown_header + conv_trait_tiles re-render
  -> top trait auto-selected from the new drilldown

User clicks trait tile
  -> Shiny.setInputValue('conv_trait_click', tile_id)
  -> observer resolves tile_id -> bundle_key via conv_node_keys()
  -> sets conv_selected_trait(bundle_key)
  -> conv_trait_sumstats loads from regional RDS bundle(s)
  -> output$conv_trait_plot renders (single or overlaid)
```

## Design decisions

- **No gene prioritization.** Earlier iterations had a gene-attributed
  heatmap showing which genes cis eQTL/pQTL pointed to. User feedback:
  that overstates the evidence and feels like prioritization. Removed.
- **`filtered_region_data()` vs `filtered_data()`**: `filtered_data()`
  caps at `input$max_traits` (default 50) for performance of the legacy
  network. That cap can drop cis molecular colocs when a region has 50+
  top-PP.H4 trans signals. The Convergence tab uses the uncapped
  `filtered_region_data()` so all evidence is counted; the drilldown
  still caps per-category tiles at 50 by signal strength for render
  speed.
- **Tile ids use index + separate `conv_node_keys` map.** Using bundle
  keys as DOM ids is unreliable because of dots / slashes / special
  characters. Indexed ids (`conv_t_1`, `conv_t_2`, ...) sidestep the
  issue; the server-side `conv_node_keys` reactiveVal holds the id ->
  bundle_key mapping.
- **Drilldown right, categories left.** Lets the user flip categories
  without losing visual context.
- **De-emphasize PP.H4 in tooltips.** Every trait on screen already
  satisfies PP.H4 >= 0.8 by atlas definition, so the number is not
  discriminating.

## Known quirks

- The classic "Missing value where TRUE/FALSE needed" warnings at
  startup come from an unrelated legacy code path; not functional.

## Files

- `inst/shiny/R/convergence.R` - `TRAIT_CATEGORIES`,
  `category_for_study()`, `count_by_category()`.
- `inst/shiny/app.R` - Convergence tab UI (search for
  `conv_category_cards`, `conv_trait_tiles`, `conv_trait_plot`) and
  server logic.
- `inst/shiny/R/plots.R` - `plot_regional_association_interactive()`
  and (for virtual studies) `plot_regional_overlay()`.
- `inst/shiny/www/landing.css` - landing page CSS (category grid,
  trait buttons).
