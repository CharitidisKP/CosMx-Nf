# shiny

Interactive explorer for the discovery-stage output. Runs on the laptop off a lightweight
bundle; it never loads Giotto or the merged object.

```r
setwd("shiny"); shiny::runApp(".")
```

## Data

The app reads a run's results folder, `../results/` by default (set `COSMX_NF` for another):

| file | written by | carries |
|---|---|---|
| `local_export/cells_umap_clusters.rds` | `export_app_bundle.sh` | one row per cell: `x`/`y`, UMAP, and every categorical annotation |
| `local_export/marker_expr_normalized.rds` | `export_app_bundle.sh` | sparse normalized expression, markers + top DE genes |
| `local_export/cell_polygons.rds` | `export_app_bundle.sh` | cell boundaries as `i`/`px`/`py` |
| `tables/<step>/markers/cluster_labels.<col>.template.csv` | DIAGNOSTICS | per-cluster triage tables (steps 03, 04, 05) |
| `tables/<step>/markers/<col>_de_top.csv` | DIAGNOSTICS | top DE genes per cluster |
| `tables/05_subclustering/sub_<res>_cells.csv` | SUBCLUSTER | subcluster per cell |
| `tables/04_reference_typing/{hieratype,supervised}/*` | HIERATYPE, SUPERVISED | per-cell scores |
| `tables/06_annotation/labels/cluster_labels.<basis>.csv` | the annotation notebooks | manual labels, joined on `cluster_id` |

Refresh the bundle by running `export_app_bundle.sh <run_id>` on the server (it reads
`Runs/<run_id>/results/objects/discover/giotto` by default) and rsyncing `local_export/` and
`tables/` down.

Without any bundle the app falls back to `demo_*.rds`, which are generated fixtures
(`make_demo.R`) and are gitignored — 17 MB of synthetic data is not worth versioning.

## Embeddings

Tabs come from whichever coordinate pairs the bundle carries: `x`/`y` (Spatial),
`UMAP1`/`UMAP2`, and `UMAPpre1`/`UMAPpre2` when the run used Harmony. The last is the same
cells embedded on PCA alone, so switching between the two UMAP tabs is a before/after of
batch correction with colour, gene and sample filters held fixed. It disappears on runs
that clustered on PCA, because then there is no pair.

## Contract

Every **factor** column in the cells table automatically becomes a "Colour by" option. New
annotations need no app changes; they only need to reach the export. The export carries
Leiden resolutions, `ht_class`/`ht_call` (HieraType), `insitutype_unsup`,
`insitutype_sup_<profile>` and the sample metadata; numeric confidence columns are
deliberately left behind.

`poly$i` is a **row index** into the cells table, not a join key, so the two files are
valid only as the pair they were written in.

Measured on the real object: 126,098 cells, 2,675,746 polygon vertices (21.2 per cell),
`cell_polygons.rds` 5.4 MB at `compress = "xz"`. The whole bundle is ~24 MB.

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
