# Multi-ancestry atlas view (virtual studies)

## Goal

Present colocalization results that come from the **same trait analyzed in
several ancestries** as a single, unified study on the landing page, the
mini-Manhattan and the regional association plots (RAPs), while preserving
full visibility of which ancestries contributed each signal.

First producer: the MVP_kidney analysis
(`genepicoloc_internal/analyses/MVP_kidney/`) which runs each curated kidney
trait (e.g. BUN_BSP_Mean_INT) against ancestry-matched MVP_R4 sumstats_2 in
AFR, AMR, EAS, EUR and META separately. That emits five per-ancestry "real"
studies sharing a stem name:

```
MVPkid_BUN_BSP_Mean_INT_AFR
MVPkid_BUN_BSP_Mean_INT_AMR
MVPkid_BUN_BSP_Mean_INT_EAS
MVPkid_BUN_BSP_Mean_INT_EUR
MVPkid_BUN_BSP_Mean_INT_META
```

Without any UI work a user would see five separate study tiles, five
Manhattans, five independent RAPs. The multi-ancestry layer collapses all
of that into one **virtual study** `MVPkid_BUN_BSP_Mean_INT` and renders it
with ancestry-aware overlays.

## Virtual study discovery

`inst/shiny/R/config.R :: discover_virtual_studies()` scans
`DEFAULT_AVAILABLE_STUDIES` for ids matching `<stem>_<ANC>` where `ANC` is
one of `AFR`, `AMR`, `EAS`, `EUR`, `META`, and collapses siblings (>=2) that
share the same atlas category into a single entry in
`DEFAULT_VIRTUAL_STUDIES`. Each entry records:

```r
list(
  category      = "MVP_kidney",
  ancestries    = c("AFR","AMR","EUR","META"),
  coloc_files   = list(AFR = "...AFR.RDS", ...),
  regional_dirs = list(AFR = "<atlas>/MVP_kidney/regional/MVPkid_..._AFR",
                       ...),
  real_ids      = c("MVPkid_..._AFR", "MVPkid_..._AMR", ...)
)
```

Only *virtual* ids are added to the 4th landing card ("MVP Kidney
(multi-ancestry)"). The per-ancestry real ids are not shown on the landing
page, but they remain loadable if an external link points at one of them.

`is_virtual_study(id)` in `R/gene_track.R` is the dispatch flag used by the
rest of the app.

## coloc_data() unification

When `current_study()` is virtual, `coloc_data()` in `app.R`:

1. Reads each ancestry's `annot_filt.RDS`, tags rows with an `ancestry`
   column.
2. Rewrites the ancestry-specific metadata prefix
   (`MVP_R4_EUR_*`, `MVP_R4_AMR_*`, ...) to a unified `MVP_R4_*`, so every
   downstream display helper sees a single schema.
3. Sets `source_study = "MVP_R4"` so existing filter/category logic works
   without branching.
4. `rbindlist(fill = TRUE)` into one data.table.

The `ancestry` column is the only surviving per-ancestry tag and drives the
overlay/collapse logic described below.

## Region overlap: forward and backward

Each ancestry's `get_coloc_regions()` picks its own window around
genome-wide-significant lead variants, so the "same" locus arrives with
different boundaries. Example PDILT/UMOD:

```
AMR   19,881,010 - 20,881,010   (1.0 Mb)
META  19,372,042 - 20,881,010   (1.5 Mb)
EUR   19,330,554 - 21,586,583   (2.3 Mb)
```

The app computes the consensus window **twice** - from two independent
sources - so the Manhattan and the RAPs stay in sync even if the atlas is
rebuilt.

### Forward: coloc-table -> Manhattan consensus dot

Inside `output$mini_manhattan` (the `is_multi` branch, `app.R`):

1. Aggregate per-ancestry rows by `(CHR_var, BP_START_var, BP_STOP_var)`.
2. Sort by `BP_START_var` within chromosome.
3. **Single-pass interval merge**: walk the sorted list keeping a running
   `end = max(stop_so_far)`. If the next region's `start <= end`, it
   belongs to the same cluster; otherwise start a new cluster. Every
   region gets a `cluster_id`.
4. Collapse each cluster to one consensus row with
   `BP_START = min(start)`, `BP_STOP = max(stop)`,
   `ancs = sort(unique(ancestry))`, `n_ancestries = length(ancs)`.

That row becomes one Manhattan dot, colored by `ANCESTRY_COVERAGE_COLORS`
(`"1"` = light orange ... `"4"` = dark red), with `anc_label` in the
tooltip listing which ancestries hit.

Clicks need a `region_key` that downstream code understands, so each
consensus dot also carries a **representative per-ancestry region key** -
the in-cluster row with the most colocs (widest as tiebreak). When the
user clicks the dot, `current_region()` becomes e.g. `16:19330554-21586583`
(EUR's window), and every region-exact match reactive (`filtered_data`,
Convergence drilldown, regional RDS lookup) still works.

### Backward: clicked region -> all overlapping regional bundles

`load_multi_region_bundles(region_str)` in `app.R` re-derives the
consensus window **from the filesystem**:

1. Take the clicked `region_str` as a seed.
2. List every `.RDS` in every `regional_dirs[[anc]]` and parse filenames
   `chr16_<start>_<stop>.RDS` back into numeric coordinates.
3. Grow the window iteratively: any file overlapping the running
   `[start, stop]` extends it; repeat until stable. This yields the exact
   same consensus window as the forward pass (when the atlas is
   consistent), but it works directly from the bundles so the RAPs never
   desync with the Manhattan.
4. Within each ancestry keep the widest overlapping bundle (usually
   exactly one per ancestry per locus).
5. Read each RDS, attach `$.consensus = list(chr, start, stop)`, cache.
6. Return `list(AMR = bundle, EUR = bundle, META = bundle)`.

The two overlap passes use the same "grow window until stable" rule on
purpose - they're derived from different inputs (coloc table vs.
filesystem) but always produce identical coordinates in practice.

## Consensus trait keys

Region bundles use per-ancestry keys like

```
MVP_R4_EUR__MVP_R4.1000G_AGR.A1C_Max_INT.EUR.GIA.dbGaP.txt.gz
MVP_R4_AMR__MVP_R4.1000G_AGR.A1C_Max_INT.AMR.GIA.dbGaP.txt.gz
MVP_R4_META__MVP_R4.1000G_AGR.A1C_Max_INT.META.GIA.dbGaP.txt.gz
```

`make_consensus_trait_key()` strips `_<ANC>` from the prefix and `.<ANC>.`
from the filename body to produce one canonical key per trait:

```
MVP_R4__MVP_R4.1000G_AGR.A1C_Max_INT.GIA.dbGaP.txt.gz
```

`resolve_consensus_in_bundle()` inverts the mapping for a given ancestry:
it synthesizes the expected `<prefix>_<anc>__...<anc>.GIA...` key and, if
not found, falls back to fuzzy matching (scan keys starting with
`<prefix>_<anc>__` and reduce each through `make_consensus_trait_key`).

The Convergence tiles and the header dropdown use these consensus keys as
their selection values, so clicking one trait immediately drives an
overlay across every ancestry that stored data for it.

## RAP overlay (plot_regional_overlay)

`inst/shiny/R/plots.R :: plot_regional_overlay()` is triggered when
`plot_regional_association_interactive()` receives a **named list** of
data.tables (one per ancestry). It draws one scatter trace per ancestry
using `ANCESTRY_COLORS` (Okabe-Ito). Markers >= -log10P 5 are solid with
hover tooltips; background points are faded.

Both the **base study** (upper panel) and the **selected trait** (lower
panel) get this treatment whenever `current_study()` is virtual:

- `regional_base_sumstats()` returns a named list of per-ancestry base
  sumstats for virtual studies.
- `conv_trait_sumstats()` resolves the consensus trait key per ancestry
  bundle and returns a named list.

The shared x-axis is pinned to the consensus window
(`bundles[[1]]$.consensus$start..stop`), so the upper and lower panels line
up perfectly across ancestries.

## "No coloc." traces

A subtle but important semantic: the region **bundles** were extracted
with `save_sumstats = TRUE`, which writes raw per-region sumstats for every
trait with any signal in the window - not only the traits that passed the
PP.H4 threshold. So when the user selects e.g. `Creat_BSP_Max_INT`, the
loader may find raw sumstats in **all 3** ancestries even though only
`[AMR, EUR]` actually colocalized.

The overlay handles this honestly instead of hiding it:

1. `conv_trait_sumstats()` attaches
   `attr(out, "coloc_ancestries") = c("AMR","EUR")` based on the consensus
   row in `conv_drilldown_data()`.
2. `plot_regional_overlay()` reads that attribute. Ancestries present in
   the list but missing from `coloc_ancestries` are rendered with
   **reduced opacity** and relabelled in the legend as
   `META (no coloc.)`.
3. The dropdown/header label (`[AMR,EUR]`) shows only the truly
   colocalized ancestries, matching the solid traces.

This lets the user see "the signal is there in META too, it just didn't
pass our PP.H4 threshold here" without misrepresenting the coloc result.
Non-colocalized ancestries in the lower panel are context, not evidence.

## Drilldown (filtered_region_data / conv_drilldown_data)

`filtered_region_data()` defaults to region-exact matching
(`CHR == chr & BP_START == start & BP_STOP == end`), which would restrict
the Convergence tiles to the **representative** ancestry's colocs only.
For virtual studies it is widened to an **overlap match** against the
consensus window:

```r
dt <- coloc_data()[CHR_var == cons$chr &
                   BP_START_var <= cons$stop &
                   BP_STOP_var  >= cons$start]
```

`conv_drilldown_data()` then collapses rows per consensus trait:

- `consensus_key` via `make_consensus_trait_key(make_bundle_key(...))`
- `ancestries` / `n_ancestries` from the set of rows that share the
  consensus key
- the strongest-signal row's metadata (after sorting by
  `-sumstats_2_max_nlog10P`) is kept as the representative for display
- .SD scoping **gotcha**: inside the `dt[, { ... }, by = .ck]` block, any
  nested data.table op (e.g. `first[, foo := ...]`) re-binds `.SD` to the
  nested frame. Capture group-level vectors (`group_ancs`, `group_ck`) in
  plain R variables **before** touching `first`, otherwise `ancestries`
  silently collapses to a single element.

The tile UI reads `ancestries` as the badge (`[AMR,EUR,META]`) and uses
`consensus_key` as the click value, so clicking a tile triggers the
multi-ancestry overlay instead of a single-ancestry one.

## Manhattan gene label

The consensus dot's label is picked by:

1. For each cluster, find the per-ancestry `regions_data()` row whose
   **midpoint is closest** to the cluster midpoint.
2. Use that row's `nearest_gene_1` (or `Prioritized_Gene` when available).
3. Expose `nearest_genes_10` in the hover tooltip as `"Nearby: ..."` so
   the user sees the full 10-gene neighborhood, not just the single
   nearest gene.

This is deliberately simple and honest: the annotator call
(`coloc_out_annotate`) is already based on the region center, so we do
not second-guess it with ancestry voting or strongest-signal overrides.
The tooltip lets the user sanity-check the label against neighbors - e.g.
for the classic PDILT/UMOD locus the label may render as `ACSM2A(within)`
because the EUR lead variant sits inside ACSM2A, with `Nearby: ACSM2A,
ACSM5, PDILT, ACSM2B, UMOD, GP2, ...` showing PDILT/UMOD right next door.

## Known limitations

- **Trait View tab** (trait-centered drilldown from the Network tab) is
  not yet multi-ancestry aware. It operates on the per-ancestry key for
  the representative region only.
- **Convergence category card counts** for the left sidebar are computed
  from `filtered_region_data()` and therefore reflect the widened
  consensus window, but they do not yet deduplicate by consensus trait -
  multi-ancestry traits count once per ancestry. Cosmetic.
- **EAS** is a first-class ancestry in the discovery list but currently
  absent from the chr16 test bundle because no EAS trait colocalized in
  the `--test` run.
- Virtual discovery only recognizes the `AFR|AMR|EAS|EUR|META` suffix
  set. Adding SAS or sub-cohorts requires editing `ANC_CODES` and
  `ANCESTRY_COLORS` in `config.R`.

## Why two separate merge passes?

It's tempting to compute consensus clusters once and carry them forward.
The two-pass design is intentional:

- **Manhattan** should work from `coloc_data()` only, so category filters
  (PP.H4, nlog10P, selected studies) apply consistently and the dot count
  matches the tile count.
- **RAP overlay** should work from the filesystem only, because the
  regional bundles are the ground truth for what can actually be plotted.
  If the coloc table has a row whose regional RDS was not extracted, the
  forward pass would lie about what the overlay can show.

Keeping them independent means a broken or partial extraction shows up
as "fewer overlaid traces" rather than "empty plots with confusing
headers". They converge because both passes implement the same interval
merge rule.
