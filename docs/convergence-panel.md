# Convergence tab

## Goal

Give the user a **region-centric multi-omics view** that answers:
"What traits colocalize with the GWAS signal in this region, grouped by
evidence type, and what does the spatial overlap look like?"

The tab is intended to eventually replace the Network + Regional Plot tabs
once users are comfortable with it.

**Philosophy:** present evidence agnostically, no gene prioritization.
The atlas shows what colocalizes; the user decides what it means.

## Layout (top to bottom, single column)

1. **Trait of interest** - regional association plot of the study-1 GWAS
   (e.g. eGFR) for the current region. Reuses
   `plot_regional_association_interactive()`.

2. **Selected trait RAP** (conditional) - when the user clicks a trait node
   in the category drill-down, its regional plot appears directly below the
   trait of interest for spatial comparison. Shares x-range with the base
   plot. Disappears when no trait is selected.

3. **Genes in region** - gene track from `plot_gene_track()` in
   `R/gene_track.R`. Shows all protein-coding genes in the 1Mb coloc window.

4. **Colocalized traits** - category cards (left, `column(5)`) + drill-down
   network (right, `column(7)`).

   - 6 category cards stacked vertically on the left. Each is a clickable
     Shiny `actionButton` with colored stripe, icon, label, and count.
   - Clicking a card selects that category (toggle to close).
   - Right side shows a `visNetwork` of up to 50 top traits in the selected
     category, colored by study, sized by sumstats_2 -log10P.
   - Clicking a trait node triggers the RAP to render above (panel 2).

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

Kidney diseases and kidney imaging are kept separate from general PheWAS
because this is a kidney-focused atlas and we want those signals to stand out.

## Key reactives

- `conv_region_coords()` - parses `current_region()` into chr/start/end
- `conv_cat_counts()` - named int vector of trait counts per category for the
  current filtered region data (uses `filtered_region_data()`, NOT the
  max_traits-limited `filtered_data()`)
- `conv_selected_cat` - reactiveVal, currently selected category id or NULL
- `conv_selected_trait` - reactiveVal, currently clicked trait's bundle_key
- `conv_node_keys` - reactiveVal, map of network node id -> bundle_key
- `conv_drilldown_data()` - filtered to the selected category, top 50 by
  sumstats_2 -log10P
- `conv_trait_sumstats()` - loads the regional RDS for the selected trait

## Click flow

```
User clicks category card
  -> observer sets conv_selected_cat(cat_id)
  -> conv_drilldown_header + conv_drilldown_network re-render
  -> visNetwork shows traits in that category

User clicks trait node
  -> visEvents(select) fires JS that sets input$conv_trait_click = node_id
  -> observer resolves node_id -> bundle_key via conv_node_keys()
  -> sets conv_selected_trait(bundle_key)
  -> conv_trait_sumstats loads from regional RDS bundle
  -> output$conv_trait_plot renders
  -> conditionalPanel shows the RAP below the base plot
```

## Completed

- [x] Base plot of trait of interest
- [x] Gene track in region
- [x] 6 category overview cards (actionButton with inline style)
- [x] Drill-down network per category (top 50 by sumstats_2 -log10P)
- [x] Click trait node -> RAP stacked under trait of interest
- [x] Uses filtered_region_data (no max_traits cap) so cis signals aren't dropped
- [x] State resets on region change and category change
- [x] MHC excluded at extraction time

## Roadmap

### Short term
- [ ] **Tile grid for traits** instead of physics network: cleaner layout,
      easier to skim when there are 30-50 traits. Sorted by -log10P.
      Each tile shows trait name + study badge + -log10P bar. Click to select.
- [ ] **Add Pan-UKB** to `study_colors` in `config.R` once the atlas includes it
- [ ] **Hide the Network + Regional Plot tabs** once the Convergence tab
      covers all use cases (feature flag or delete)
- [ ] Trace the "missing value where TRUE/FALSE needed" warnings to their
      root cause (not blocking, but noisy in logs)
- [ ] **De-emphasize PP.H4** in tooltips (all traits are already PP.H4 >= 0.8
      by atlas definition, so the value isn't discriminating)

### Medium term
- [ ] **Category filter inside drill-down** - when a category has hundreds of
      traits (e.g. UKB_PPP trans in a hub locus), add a search box to filter
      by trait name
- [ ] **Sort controls** for the drill-down: by -log10P, by study, alphabetical
- [ ] **Highlight cis molecular colocs** in the gene track (red outline on
      genes that have cis eQTL/pQTL in the selected category)
- [ ] **Export current view** as PNG/PDF (full-page: base + selected trait +
      gene track + category summary)

### Long term
- [ ] **Cross-region comparison**: select two regions and compare their
      category profiles side by side
- [ ] **LLM "interpret this convergence"** button that pre-fills the chat
      with the full convergence context (base + category counts + top traits)

## Design decisions (for future reference)

- **No gene prioritization**. Earlier iterations had a "gene-attributed
  heatmap" showing which genes cis eQTL/pQTL pointed to. User feedback: this
  overstates the evidence and feels like prioritization. Removed.
- **filtered_region_data() vs filtered_data()**: `filtered_data()` caps at
  `max_traits` (default 50) for network rendering performance. That cap drops
  cis molecular colocs when a region has 50+ top-PP.H4 trans signals. The
  Convergence tab uses the uncapped `filtered_region_data()` so all evidence
  is counted; the drill-down network still caps at 50 per category for
  rendering speed, but top-by-signal-strength within each category.
- **Node IDs use index + separate bundle_key map**. Early attempt used the
  bundle_key directly as node ID, but visNetwork/plotly JS events don't
  reliably handle IDs with dots, slashes, and special characters. Indexed
  IDs (`conv_t_1`, `conv_t_2`, ...) avoid the issue; the server-side
  `conv_node_keys` reactiveVal stores the id -> bundle_key mapping.
- **Drill-down stays on the right while categories stay on the left** so the
  user can click between categories without losing their place.
- **"Missing value where TRUE/FALSE needed" warnings** at startup are from an
  unrelated legacy code path and don't affect functionality. Not fixed yet.

## Files

- `inst/shiny/R/convergence.R` - `TRAIT_CATEGORIES`, `category_for_study()`,
  `count_by_category()`
- `inst/shiny/app.R` - Convergence tab UI (~lines 260-285) and server logic
  (~lines 2070-2320)
- `inst/shiny/www/landing.css` - CSS for the old `.conv-cat-card` class
  (currently unused since cards are inline-styled `actionButton`s, but kept
  in case we switch back)
