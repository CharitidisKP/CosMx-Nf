# shiny

Interactive explorer for a run's results. It runs where the run lives, on the server, off a
lightweight bundle; it never loads Giotto or the merged object. Label edits are saved into the
run's labels folder, which the pipeline reads, so nothing has to be copied between machines.

In Positron on the server, open `shiny/app.R` and press Run Shiny App, the ▶ at the top left of
the editor. Positron runs the app in its R console and shows it in the Viewer; the console's stop
button ends it. That R needs the app's packages (shiny, bslib, scattermore), for example an R
running inside the pipeline's image. From a console in the pipeline folder the same is

```r
shiny::runApp("shiny")
```

The app reads the run the pipeline ran last, among the runs beside the pipeline that have an app
bundle; the sidebar names it. To read another, set `COSMX_NF` to its results folder first:

```r
Sys.setenv(COSMX_NF = normalizePath("../Runs/<run_id>/results"))
```

If the Viewer stays blank, the app's port has not reached your machine: forward it in Positron's
Ports panel, or run `ssh -N -L <port>:127.0.0.1:<port> <server>` and open
`http://localhost:<port>`.

The page is drawn at 85% (`font-size-root` in `theme.R`), as a browser zoom of 85% would draw
it. The plot fills the card above the Tables bar and never runs past it; the legend beside it
scrolls on its own. **Figure size** (− 100% + next to Reset zoom, or Cmd − and Cmd +, Cmd 0 to
reset) scales the plot's text, points and lines and leaves the page and the legend as they are.
It works in Positron's Viewer, which cannot zoom the page; click the app first so the keys reach
it. For more room, hide the left sidebar with the ‹ at its top, or press the expand icon at the
plot card's bottom-right corner for full screen.

The Spatial tab shows one sample at a time, picked next to Jump to FOV, which lists that
sample's FOVs. All samples stacks every sample ticked under Samples; it draws every cell's
outline, so it is slow on a whole run. The other tabs always show every ticked sample.

## Data

The app reads a run's results folder, the run launched last by default (set `COSMX_NF` for
another), and the run's labels from `labels/` beside it (set `COSMX_LABELS` for another):

| file | written by | carries |
|---|---|---|
| `local_export/cells_umap_clusters.rds` | `export_app_bundle.sh` | one row per cell: `x`/`y`, UMAP, and every categorical annotation |
| `local_export/marker_expr_normalized.rds` | `export_app_bundle.sh` | sparse normalized expression, every gene (or markers + top DE genes, export option) |
| `local_export/cell_polygons.rds` | `export_app_bundle.sh` | cell boundaries as `i`/`px`/`py` |
| `tables/<step>/markers/cluster_labels.<col>.template.csv` | DIAGNOSTICS | per-cluster triage tables (steps 03, 04, 05) |
| `tables/<step>/markers/<col>_de_top.csv` | DIAGNOSTICS | top DE genes per cluster |
| `tables/05_subclustering/sub_<res>_cells.csv` | SUBCLUSTER | subcluster per cell |
| `tables/04_reference_typing/{hieratype,supervised}/*` | HIERATYPE, SUPERVISED | per-cell scores |
| `tables/06_annotation/labels/cluster_labels.<basis>.csv` | the annotation notebooks | draft labels, joined on `cluster_id` |
| `../labels/*` | the app and you | the run's curated labels, see Labels below |

The bundle changes only when the object does. After a stage that writes a new object
(discover, subcluster, annotate, spatial), run `export_app_bundle.sh <run_id> all xz <object>`
on the server, e.g. `objects/spatial` for the niche columns, and restart the app. Label edits
need no new bundle: the app builds its label columns from the labels folder.

Without any bundle the app falls back to `demo_*.rds`, which are generated fixtures
(`make_demo.R`) and are gitignored — 17 MB of synthetic data is not worth versioning.

## Embeddings

Tabs come from whichever coordinate pairs the bundle carries: `x`/`y` (Spatial),
`UMAP1`/`UMAP2`, and `UMAPpre1`/`UMAPpre2` when the run used Harmony. The last is the same
cells embedded on PCA alone, so switching between the two UMAP tabs is a before/after of
batch correction with colour, gene and sample filters held fixed. It disappears on runs
that clustered on PCA, because then there is no pair.

## Selecting and editing labels

The mode switch above the plot sets what a drag does:

- **Inspect** (default): click a cell to see its details, and drag a box to zoom to it.
- **Select**: drag a box. The bar above the plot filters every view to those cells or assigns
  them a label.
- **Draw**: press and drag to outline any shape. The outline stays in place while you zoom,
  and the same bar filters or labels the cells inside it.

In every mode, pinch or Cmd/Ctrl + scroll zooms and double-click resets.

### Labels

The app edits the run's labels where the pipeline reads them: `labels/` beside the run's
results (`Runs/<run>/labels/`), or the folder in `COSMX_LABELS`. Every edit is saved as it is
made, and the next `annotate` applies it with no flags. A later stage notices labels that
changed since its object was built and runs annotate (and spatial) again first.

| File | Holds | Written by |
|---|---|---|
| `cluster_labels.<clustering>.csv` | `cluster_id,label` for every cluster | Tables > Cluster labels |
| `label_overrides.csv` | `cell_ID,label,confidence,note` per edited cell | Annotate |
| `cluster_basis.txt` | the clustering `cell_type` comes from | Use for cell_type |

**Tables > Cluster labels** labels whole clusters. Pick the clustering and a cluster, pick or
type a label, then press Set. A clustering starts from its saved labels, else the pipeline's
draft in `tables/06_annotation/labels/`, else the majority label of its cells. **Use for
cell_type** makes annotate take `cell_type` from the clustering shown. **Discard edits** puts
the clustering's labels back as they were when the app started.

**Annotate** labels the selected cells. Pick a label or type a new one, pick a confidence, then
press Assign there or in the bar. `(cluster label)` hands cells back to their cluster's label,
and Undo steps back one assignment.

- **High** confidence sets `cell_type`, the label every analysis uses.
- **Low** confidence leaves `cell_type` alone. The label goes into `cell_type_display`, which
  figures draw with white stripes, and the cell is marked `low` in `label_confidence`.

Of a cell's edits the later one wins, except that a low-confidence label never replaces a
newer high-confidence one. `bin/label_state.R` holds these rules for both the app and annotate.

Three colourings show what annotate will write: `custom_annotation` (`cell_type`),
`custom_display` (`cell_type_display`) and `custom_confidence` (`label_confidence`). They head
the Colour by list, under Labels - live; `cell_type`, `cell_type_display` and
`cell_type_confidence` further down are what the last annotate wrote into the bundle. The app
opens on `custom_annotation`. Coloured by `custom_display`, each low-confidence cell carries a
white slash, as the figures stripe it. Coloured by a label, Cluster triage and Top DE show the
clustering picked in Cluster labels; highlight labels to list the top DE genes of the clusters
labelled with them. The pipeline computes DE per clustering, not per label.

Tables opens under the plot: click its bar to open or close it, and drag the bar up for more
table or down for more plot. Each table scrolls inside the panel, with its header row and first
column kept in view.

Each label keeps one colour, from `assets/label_colours.csv`, which also gives the name figures
show (for example `Plasma_cell` as B cell lineage). Kidney cell types are pastel and immune
types vivid, drawn on top; untyped labels are grey and drawn underneath. A label not in the
file takes a spare colour.

The Download buttons save copies of the two CSVs, e.g. to pass to another run with `--labels`
and `--label_overrides`. `annotate` writes `label_overrides_applied.csv`, which counts the
from-to changes, and `label_source.csv`, which names the files it used.

## Contract

Every **factor** column in the cells table automatically becomes a "Colour by" option. New
annotations need no app changes; they only need to reach the export. The export carries
Leiden resolutions, `ht_class`/`ht_call` (HieraType), `insitutype_unsup`,
`insitutype_sup_<profile>` and the sample metadata; numeric confidence columns are
deliberately left behind.

`poly$i` is a **row index** into the cells table, not a join key, so the two files are
valid only as the pair they were written in.

Measured on a 148,116-cell run with every gene: 3,107,462 polygon vertices (21 per
cell), `cell_polygons.rds` 5.9 MB at `compress = "xz"`. Each expression matrix is ~125 MB and
the whole bundle ~295 MB.

Columns that never vary within a sample — treatment, timepoint, patient, batch — are
withdrawn from "Colour by" on the Spatial view, because colouring tissue by them just
recolours the block of space each sample occupies. They stay available in "Split by", which
is the control that actually makes a multi-sample tissue view readable.

## Files

| file | |
|---|---|
| `app.R` | the application |
| `load_cosmx_nf.R` | reads a run's results tree into the contract |
| `export_app_bundle.sh` | server-side exporter, Giotto object to bundle |
| `from_giotto.R` | minimal statement of the contract `app.R` depends on; superseded for pipeline use by `export_app_bundle.sh` |
| `make_demo.R` | generates the synthetic fallback fixtures |
| `theme.R` | palettes and the bslib theme |

## Gene combination panel

Pick 2+ genes in the sidebar to switch from single-gene expression to co-detection. Cells
are coloured by how many of the picked genes are expressed (or All / Any), and the
**Gene combination** tab reports observed vs expected cells at each threshold.

The expected column is the point. With markers detected in a few percent of cells, chance
co-detection is common and it *rises with sequencing depth*, so "N cells express all three"
means nothing on its own. The baseline is computed within total-count deciles, which holds
depth fixed and shrinks the apparent enrichment a naive baseline reports. Genes that
truly co-express, such as the MHC-II set in antigen-presenting cells, should climb steeply
with the threshold; chance co-detection stays near 1x.

Needs `total_counts` in the bundle for the depth-matched baseline; without it the table
falls back to a naive baseline and says so. Re-run `export_app_bundle.sh` to get it.
