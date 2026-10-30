# CosMx-Nf

A Nextflow (DSL2) pipeline for single-cell spatial transcriptomics from the Bruker (NanoString)
CosMx Spatial Molecular Imager. It takes each slide's flat-file export through QC, integration,
clustering and reference-based cell typing. A human-in-the-loop annotation step follows, then
spatial statistics, cell-cell communication and neighbourhood-aware differential expression.
Every step is an R script in `bin/` built on [Giotto](https://giottosuite.com).

It was written for multi-sample FFPE kidney biopsy studies on the CosMx 6K Discovery panel. The
marker panel (`assets/canonical_markers.csv`) and the default reference profiles are
kidney-specific; the rest of the pipeline is tissue-agnostic.

## How a run proceeds

The pipeline runs in stages. Each stage is its own launch (`--stage`) and picks up the Giotto
object the previous stage published, so the manual annotation step can sit between them.

```mermaid
flowchart LR
    SS[samplesheet.csv] --> D[discover]
    D --> SC["subcluster<br/>(optional)"]
    D --> MA["manual annotation<br/>notebooks/render.sh"]
    SC --> MA
    MA -- "draft labels" --> LB["labels folder<br/>(by hand or in the app)"]
    LB --> AN[annotate]
    AN --> SP[spatial]
    SP --> CCI[cci]
    SP --> DE[de]
```

| Stage | Processes | Methods |
|---|---|---|
| `discover` (default) | LOAD, QC, MERGE, INTEGRATE, CLUSTER, HIERATYPE, INSITUTYPE, SUPERVISED, DIAGNOSTICS | Giotto `importCosMx`, with negative probes as their own feature type. Per-sample QC: the count floor is `min(count_cap, count_quantile)`, so it never removes more than a set share of a sample; plus an area cap, an FOV flag and border-cell flags. Normalisation, then the top HVGs unioned with the marker panel, then PCA, then Harmony over the batch columns that vary. A shared-nearest-neighbour graph is clustered with Leiden over a resolution ladder, with UMAPs for the chosen resolutions. HieraType marker typing. Unsupervised InSituType with cohorts from the immunofluorescence channels and neighbourhood expression. Supervised InSituType against reference profiles. Marker-panel diagnostics per chosen column write a label template with a best guess per cluster. |
| `subcluster` | SUBCLUSTER, DIAGNOSTICS | Re-embeds each cluster of a basis column on its own (own HVGs, PCA, Harmony, graph, Leiden) and adds one `sub_<res>` column. Clusters below `subcluster_min_cells` are left whole. |
| manual annotation | `notebooks/render.sh` | Knits `manual_annotation.Rmd` in the container. It covers resolution overview, marker heatmaps and volcano plots, marker-score-guided draft labels, lineage and immunofluorescence checks, and a validation against the annotate contract. It writes a `cluster_id,label` CSV for you to curate. |
| `annotate` | APPLY_LABELS, COMPOSITION | Joins the run's curated cluster labels onto their clustering (`cell_type`, `cell_type_l1`) and draws a marker dot plot. The per-cell edits in the same folder go on top. A low-confidence edit leaves `cell_type` alone: it changes only the label figures draw (`cell_type_display`, striped) and marks the cell in `label_confidence`. Composition per sample, with the per-subject change against a baseline inside each treatment group. |
| `spatial` | SPATIAL_NET, BANKSY_PREP, BANKSY_SAMPLE, BANKSY, BCELL, IMMUNE_DENSITY | A Delaunay network with a maximum edge length, and cell-type proximity enrichment per sample against a permutation null. BANKSY niches per sample, one task each, on the variable genes that rise well above the sample's background, then matched across samples into niche types by their cell-type mix (`niche_type_lam<lambda>`). The settings are the `banksy_settings` row of `assets/banksy_variants.csv`. B-lineage subclustering. Immune cells per mm² of occupied tissue. |
| `cci` | CCI_LIANA, CCI_NICHENET, CCI_MISTY, CCI_MISTY_COLLECT, CCI_SPARKX, CCI_NNSVG, CCI_INSITUCOR | LIANA ligand-receptor consensus weighted by Delaunay contact. NicheNet ligand activity per comparison and receiver type. MISTy intra, juxta and para views per sample. SPARK-X spatially variable genes beyond cell-type composition, refined with nnSVG. InSituCor spatial co-expression modules. |
| `de` | SMIDE, SMIDE_OVERLAP, SMIDE_META, PSEUDOBULK_DE | smiDE per sample x cell type (a GLMM on niche type, with neighbour-expression and spatial terms), pooled across samples by fixed-effect meta-analysis; hits whose gene is far higher in neighbouring cells of that niche are flagged as possible transcript bleed. Pseudobulk edgeR per cell type, each later timepoint against baseline within subject, with the same bleed flag. |
| `banksy_sweep` | BANKSY_SWEEP_PREP, BANKSY_SWEEP, BANKSY_SWEEP_COLLECT | Optional. Runs each biopsy through every row of `assets/banksy_variants.csv`, with the same code as the `spatial` stage, and scores the niches: agreement with FOVs (low is better), with cell types, and spatial coherence. For choosing `banksy_settings`. |
| `diagnose` | DIAGNOSTICS | Re-runs the diagnostics on existing columns without recomputing anything upstream. |
| `supervised` | SUPERVISED, DIAGNOSTICS | Re-runs the supervised InSituType passes against the published discover object. |
| `report` | REPORT | Rebuilds `run_report.md` from the published results. |

Every stage ends by rewriting `run_report.md`, a Markdown summary of the run so far. It covers
samples and QC, integration settings, clusters per resolution, reference typing, subclusters,
cell-type composition, spatial statistics, CCI and DE hits, and links to the main figures.

## Requirements

- Nextflow 24.04 or later (the configuration uses `resourceLimits`). Developed on 26.04, where
  `nextflow lint` reports no errors.
- Apptainer. Every task runs inside one image with R 4.5 and Giotto 4.2. The other R packages
  can live in the image or in an R library bound into it (`r_libs` in `conf/site.yaml`):
  data.table, ggplot2, ggrepel, Matrix, matrixStats, irlba, harmony, scran,
  SingleCellExperiment, SpatialExperiment, InSituType, HieraType, Banksy, liana, nichenetr,
  mistyR, future, SPARK, nnSVG, InSituCor, smiDE, jsonlite, sf, ragg and patchwork. The app
  also needs shiny, bslib and scattermore. Giotto's Leiden goes through
  reticulate to a Python with leidenalg and python-igraph (`Environment/environment.yml`
  lists the Python side).
- Reference profile matrices for the supervised passes (genes x cell types, `<name>.csv` or
  `.rds`) are not shipped. Stage them in `assets/profiles/`, or set `supervised_profiles` to
  `""` to skip the supervised passes.
- NicheNet ligand-target model files, only for the NicheNet step (`nichenet_model_dir`).

## Running

```bash
git clone https://github.com/CharitidisKP/CosMx-Nf.git
cd CosMx-Nf
cp conf/site.example.yaml conf/site.yaml   # image, R libraries, bind mounts, samplesheet
./run.sh my_run -profile server            # discover
```

`run.sh <run_id> [nextflow options]` is the entry point. It keeps each run in one place,
`../Runs/<run_id>/{results,logs}`, with `work/` on local scratch (`$NXF_SCRATCH_ROOT`, default
`/tmp/$USER-nf`). It passes `conf/site.yaml` as `-params-file` and rejects any `--param` no
config declares, because Nextflow otherwise accepts a mistyped parameter silently. Without it:

```bash
nextflow run main.nf -profile server -params-file conf/site.yaml --outdir results --run_id my_run
```

A full pass, then:

```bash
# Optional: split every cluster of a resolution into subclusters
./run.sh my_run -profile server --stage subcluster --subcluster_basis leiden_clus_res0.3

# Draft labels for a clustering, written to results/tables/06_annotation/labels/
./notebooks/render.sh my_run sub_res0.3

# Curate them in the run's labels folder, by hand or in the Shiny app
RUN=../Runs/my_run
mkdir -p $RUN/labels
cp $RUN/results/tables/06_annotation/labels/cluster_labels.sub_res0.3.csv $RUN/labels/

./run.sh my_run -profile server --stage annotate
./run.sh my_run -profile server --stage spatial
./run.sh my_run -profile server --stage cci
./run.sh my_run -profile server --stage de
```

- **Labels** live in `Runs/<run_id>/labels/`, the folder the Shiny app edits:
  - `cluster_labels.<clustering>.csv` (`cluster_id,label`) for each clustering you label;
  - `label_overrides.csv` (`cell_ID,label,confidence,note`) for per-cell edits;
  - `cluster_basis.txt`, naming the clustering `cell_type` comes from. It is needed only when
    more than one clustering has labels.

  `annotate` needs no flags. `--labels`, `--label_overrides` and `--cluster_basis` each replace
  one file, and `--labels_dir` points at another folder.
- **Changed labels** are picked up on their own. `annotate` and `spatial` record the labels
  their objects were built from (`label_state.txt`). When a later stage finds the labels
  changed, it reruns `annotate`, and `spatial` if needed, first in the same launch. `spatial`
  also reruns when `banksy_settings` changes.
- **`de`** tests every typed label in `cell_type`, one smiDE job per sample and cell type, unless
  `smide_celltypes` names a subset.
- **`-resume`** works on any stage. Tasks use `cache = 'lenient'` (path and size, not timestamps), so a
  shared file system with unreliable timestamps does not invalidate the cache.
- **`--solo <id,id|all>`** also takes the named samples through discovery on their own, without
  MERGE or Harmony (`--solo_only` skips the merged run).
- **`--merged_giotto <dir>`** points any stage after `discover` at a different Giotto object
  than the one the previous stage published.
- **`assets/samplesheet.smoke.csv`** loads three FOVs per sample, for a quick end-to-end check
  of `discover`.

### Profiles

| Profile | What it sets |
|---|---|
| `server` | Local executor on one large machine, every task in the Apptainer image, per-task and total CPU/memory caps. |
| `slurm` | The same container settings with the SLURM executor, a time limit per resource label and `--partition`. |

Resources are set per label in `conf/base.config`. A task that dies with a kill or
out-of-memory exit code is retried once with twice the memory.

## Samplesheet

One row per sample; see `assets/samplesheet.csv`. The real sheets, which carry data paths,
belong outside git: point `input` in `conf/site.yaml` at them (`assets/private/` is gitignored).

| Column | Meaning |
|---|---|
| `sample_id` | Unique sample name, used in every output and by `--solo`. |
| `include` | `TRUE`/`FALSE`. `FALSE` drops the sample from every stage without deleting its row. A missing column counts as `TRUE`. |
| `file_prefix` | Prefix of the flat files in `data_dir` (`<prefix>_exprMat_file.csv.gz`, `<prefix>_metadata_file.csv.gz`, ...). |
| `data_dir` | The sample's flat-file export folder: expression matrix, cell metadata, FOV positions and polygons. Absolute, or relative to the launch directory. |
| `patient_id` | Specimen or patient code, stored on every cell. Use a study code, never a hospital or accession number. |
| `treatment` | Study group. Comparisons stay inside each group (`compare_within`). |
| `timepoint` | Compared against `compare_baseline` (`compare_by`). |
| `slide_id`, `batch` | Physical slide and instrument run; with `batch_column` these are the Harmony covariates. |
| `subset_fovs` | Optional. Semicolon-separated FOV ids to load, e.g. when several samples share one slide's export. |
| `subject_id` | Optional. The person; paired changes in `annotate` and NicheNet are computed within a subject (`pair_by`). |

## Configuration

Every parameter is declared, with a comment, in `nextflow.config`, grouped by step. The ones
most runs touch:

| Parameter | Default | |
|---|---|---|
| `input` | `assets/samplesheet.csv` | Samplesheet. |
| `chosen_resolutions` | `0.3,0.5` | Leiden resolutions that get a UMAP and diagnostics; every resolution in `leiden_res_ladder` is still clustered. |
| `markers` | `assets/canonical_markers.csv` | Marker panel (`cell_type,compartment,gene`) behind the diagnostics and HVG union. |
| `supervised_profiles` | `HCA_Kidney,KPMP` | Reference profiles for supervised InSituType, resolved in `profile_dir`. |
| `compare_by`, `compare_baseline`, `compare_within`, `pair_by` | `timepoint`, `T0`, `treatment`, `subject_id` | The study design, read from the samplesheet by composition and NicheNet. |
| `celltype_column` | `cell_type` | Label used by `spatial`, `cci` and `de`; can be any typing column instead. |
| `untyped_celltypes` | `Low_count` | Labels without a cell identity: kept in the object, left out of LIANA, NicheNet and the typed proportions. |
| `label_colours` | `assets/label_colours.csv` | Colour, group (immune, non_immune or untyped) and display name of each label, for every figure and the app. |
| `fig_format` | `png` | Formats of the label and niche figures: `png`, `pdf`, `jpeg`, a comma list, or `all`. Tissue maps are always PNG. |
| `banksy_settings` | `docs_above_background` | The row of `assets/banksy_variants.csv` the `spatial` stage runs. |

Machine-specific settings (image, R libraries, Python, bind mounts, profile and NicheNet folders,
the real samplesheet) go in `conf/site.yaml`, copied from `conf/site.example.yaml`. It is a
params file rather than a `-c` config because params set in a `-c` file do not reach the
profile settings built from them.

## Outputs

```
Runs/<run_id>/
├── labels/                            curated labels: cluster_labels.<clustering>.csv, label_overrides.csv, cluster_basis.txt
├── results/
│   ├── run_report.md                  rewritten at the end of every stage
│   ├── figures/<step>/                01_qc ... 09_de
│   ├── tables/<step>/                 CSVs by the same steps; draft labels in 06_annotation/labels/
│   ├── objects/<stage>/giotto         discover, subcluster, annotate, spatial, solo/<id>
│   ├── pipeline_info/                 <stage>_{report.html,trace.txt,timeline.html,params.json}, stage_log.tsv
│   └── local_export/                  Shiny bundle (shiny/export_app_bundle.sh)
├── logs/nextflow.log
└── work -> <scratch>/<run_id>
```

Steps: `01_qc`, `02_integration`, `03_clustering`, `04_reference_typing`, `05_subclustering`,
`06_annotation`, `07_spatial`, `08_cci`, `09_de`; the BANKSY sweep writes under
`07_spatial/banksy_sweep`. Each launch also records its parameters (`<stage>_params.json`) and a
line in `stage_log.tsv` with the time, the pipeline commit and the command line.

## Repository layout

```
main.nf            entry workflow: stage checks, one named workflow per --stage, run records
nextflow.config    parameters, profiles, reports, manifest
conf/              base.config (resources), server.config, slurm.config, site.example.yaml
subworkflows/      discovery, solo, annotation (ANNOTATION, SPATIAL, DE), cci; samplesheet helpers
modules/           one process per file (insitutype_supervised.nf holds SUPERVISED and ATTACH_SUPERVISED)
bin/               the R scripts the processes call
assets/            marker panel and example samplesheets
notebooks/         manual_annotation.Rmd and render.sh, the manual annotation step
shiny/             interactive explorer and label editor for a run (see shiny/README.md)
tests/             testthat tests for the label rules and the app's label handling
Environment/       conda environment for interactive work on the same R version
run.sh             launcher: run folder, scratch work dir, log, site file, parameter check
```

## Shiny explorer

`shiny/` is an app for browsing a run and curating its labels. It has UMAP and tissue views
coloured by any annotation, gene expression, depth-matched gene co-detection, and the
diagnostics and label tables side by side. Cluster labels and per-cell edits (high or low
confidence) are saved to the run's labels folder as you make them, and the next stage applies
them. The app runs where the results are, for example from Positron or RStudio on the server,
and opens the run launched last. `shiny/export_app_bundle.sh <run_id>` writes a bundle of the
cells, expression and cell outlines next to the results, so the app never loads Giotto. Without
a bundle it falls back to synthetic data from `make_demo.R`.

## Limitations

- CosMx exports are large and patient-derived, so the repository ships no test data. The
  smoke samplesheet shows the shape of a minimal run.
- There is no CI yet. `tests/testthat` covers the label rules and the app's label handling
  (`Rscript -e 'testthat::test_dir("tests/testthat")'`). The rest is checked with
  `nextflow lint`, a parse check of every `bin/` script, and by running the affected stages.
- The marker panel and default profiles assume kidney tissue.

## Acknowledgements

Built on Giotto, InSituType, HieraType, Harmony, BANKSY, LIANA, NicheNet, MISTy, SPARK-X,
nnSVG, InSituCor and smiDE. Please cite the methods you use.

## License

MIT, see [LICENSE](LICENSE).
