# Shiny app architecture

Overview of how `inst/shiny/` is structured after the Phase 1-4
refactor (April 2026).

## File map

```
inst/shiny/
├── app.R                (~1430 lines)
│   ├── UI: fluidPage with 3 tabs (Atlas, Region View, Documentation)
│   └── server: reactive values, core data reactives, observers,
│       convergence outputs, gene info panel, download handler
├── R/                   (auto-sourced before app.R by Shiny)
│   ├── annotation_db.R  (66)   trait_label(), umet_label()
│   ├── bundle_keys.R    (79)   make_bundle_key(), make_consensus_trait_key(),
│   │                            resolve_consensus_in_bundle(), ANC_CODES
│   ├── config.R         (194)  DATA_PATH, study_colors, discover_studies(),
│   │                            discover_virtual_studies(),
│   │                            ANCESTRY_COLORS, ANCESTRY_COVERAGE_COLORS
│   ├── convergence.R    (74)   TRAIT_CATEGORIES, category_for_study(),
│   │                            count_by_category()
│   ├── gene_track.R     (200)  plot_gene_track(), resolve_annot_path(),
│   │                            is_virtual_study(), gene_annotation loading
│   ├── landing.R        (180)  atlas_categories data, build_card(),
│   │                            expand_card_traits(), is_known_study()
│   ├── manhattan.R      (256)  manhattan_consensus_table(),
│   │                            manhattan_plot()
│   ├── plots.R          (292)  plot_regional_association_interactive(),
│   │                            plot_regional_overlay()
│   ├── region_bundles.R (135)  load_region_bundle(), load_base_sumstats(),
│   │                            load_trait_sumstats(),
│   │                            load_multi_region_bundles()
│   └── regions.R        (63)   parse_region_key(), parse_region_id(),
│                                parse_regional_filename(),
│                                reconstruct_name()
├── www/
│   ├── docs.md                 User guide rendered in Documentation tab
│   └── landing.css             Landing page styles
├── extdata/
│   └── shiny_annotation_db.RDS Bundled slim annotation DB
├── Dockerfile                  Podman/Docker container
└── deploy.sh                   Local/remote deployment script
```

## Design principles

1. **R/*.R files are pure** - they never read Shiny inputs or reactives
   directly. Everything they need is passed as an argument by the
   server function in app.R. This makes them testable, reviewable, and
   safe to modify without Shiny-specific knowledge.

2. **Shiny state lives in app.R** - all reactiveVals are declared at
   the top of the server function. Core reactives (coloc_data,
   regions_data, filtered_region_data, current_virtual_info) follow
   immediately. Every observer and output references these by name.

3. **One concern per file** - bundle key logic (string manipulation)
   is separate from region bundle I/O (file reading + caching), which
   is separate from the Manhattan (plotly assembly). This means a bug
   in consensus-key resolution doesn't require reading 1400 lines of
   app.R to find.

4. **Process-wide caches** in `region_bundles.R` (plain environments,
   not reactiveVals). Region RDS files are immutable between deploys,
   so a per-process cache is correct and simpler than per-session
   reactiveVal caches.

## app.R server structure (top to bottom)

```
server <- function(input, output, session) {
  Landing page render               output$landing_categories
  Reactive values (all at top)      current_study, selected_gene, ...
  Event handlers
    - study card click              observeEvent(selected_study)
    - gene / region selectize sync  observe({ regions_data() ... })
    - Manhattan chr selector        observe({ coloc_data() ... })
    - Manhattan click               observeEvent(plotly_click)
  Core data reactives
    - coloc_data                    loads + merges coloc RDS files
    - regions_data                  unique regions from coloc_data
    - current_virtual_info          resolve virtual study metadata
    - filtered_region_data          coloc rows for current region
    - filtered_data                 alias for filtered_region_data
  Mini Manhattan                    output$mini_manhattan (delegates to
                                    manhattan_consensus_table + manhattan_plot)
  Convergence view
    - conv_region_coords            current region as list(chr, start, end)
    - conv_base_plot                upper RAP (base study)
    - conv_gene_track + click       gene track + gene info panel click
    - conv_cat_counts               category counts for card badges
    - conv_category_cards           6-button horizontal grid
    - conv_drilldown_header         category label above DT
    - conv_drilldown_data           per-category trait list (consensus-
                                    collapsed for virtual studies)
    - conv_trait_table              data.frame + keys for DT
    - conv_trait_list (DT output)   searchable, sortable, filterable
    - region change observer        resets state, auto-selects first cat
    - conv_trait_sumstats           loads trait RAP data (overlay for multi)
    - conv_selected_trait_label     display name for the lower RAP title
    - conv_trait_plot               lower RAP output
  Regional bundle data
    - regional_data_path            resolves study -> filesystem dir
    - regional_base_sumstats        base RAP data (overlay list for multi)
  Trait name helpers
    - get_trait_name                filename-safe trait identifier
    - get_trait_display_name        human-readable trait label
  Gene info panel                   output$gene_info_panel + close button
  Download handler                  output$download_data (CSV export)
}
```

## Data flow

```
Landing card click
  -> current_study()
  -> coloc_data()              reads per-ancestry RDS, rbinds, unifies
  -> regions_data()            unique (CHR, BP_START, BP_STOP) + gene meta
  -> mini Manhattan            manhattan_consensus_table() + manhattan_plot()

Manhattan dot click
  -> current_region()
  -> filtered_region_data()    rows for the region (consensus-widened if multi)
  -> conv_cat_counts()         category badge numbers
  -> conv_drilldown_data()     per-category trait list

Category card click
  -> conv_selected_cat()
  -> conv_drilldown_data()     updates for the new category
  -> conv_trait_table()        DT refresh

DT row click
  -> conv_selected_trait()
  -> conv_trait_sumstats()     load trait from per-ancestry bundles
  -> conv_trait_plot            overlay RAP

Gene track click
  -> selected_gene()
  -> output$gene_info_panel    NCBI / Reactome / HPO / HGNC annotations
```

## Multi-ancestry virtual studies

Full design in `docs/multi-ancestry.md`. High-level:

- `discover_virtual_studies()` in `config.R` collapses
  `<stem>_{AFR,AMR,EUR,META}` siblings into one virtual trait id.
- `coloc_data()` virtual branch rbinds all ancestries, adds `ancestry`
  column, rewrites `MVP_R4_<ANC>_*` metadata to unified `MVP_R4_*`.
- Manhattan: `manhattan_consensus_table()` clusters overlapping
  per-ancestry regions via single-pass interval merge, colors dots by
  ancestry coverage (1-4 gradient). Click routing uses the
  representative per-ancestry region key.
- RAP overlay: `load_multi_region_bundles()` re-derives consensus from
  the filesystem. `plot_regional_overlay()` draws one Okabe-Ito trace
  per ancestry; non-colocalized ancestries rendered faded with
  "(no coloc.)" legend tag.
- Consensus trait keys: `make_consensus_trait_key()` strips `_<ANC>`
  prefix and `.<ANC>.` filename token. `resolve_consensus_in_bundle()`
  inverts the mapping per ancestry.

## Removed features (preserved on branches)

| Feature | Branch | Reason |
|---------|--------|--------|
| Network (visNetwork) subtab | `feature/network-view` | Replaced by Convergence tile grid |
| Trait View tab | `feature/trait-view` | Replaced by Convergence drilldown |
| In-app LLM chat | `feature/llm-chat` | Datenschutz |

## What's NOT modularized yet (and why)

The Convergence server reactives + outputs (~400 lines) are tightly
coupled to each other and to the shared reactiveVals. Extracting them
into a proper `shiny::moduleServer` would be the "textbook" next step,
but it changes every consumer's call signature (module namespace
prefixes) and requires testing parity. Deferred to a future pass
when/if the convergence panel grows further.
