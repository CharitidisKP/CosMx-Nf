#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  sample_id = arg("sample_id"),
  file_prefix = arg("file_prefix"),
  data_dir = arg("data_dir"),
  patient_id = arg("patient_id"),
  treatment = arg("treatment"),
  timepoint = arg("timepoint"),
  slide_id = arg("slide_id"),
  batch = arg("batch"),
  subset_fovs = arg("subset_fovs"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

## Slide number within the cell_id prefix. Samples split out of one slide's export share it ##
slide <- as.integer(sub(
  "^c_(\\d+)_.*",
  "\\1",
  fread(
    file.path(opt$data_dir, paste0(opt$file_prefix, "_metadata_file.csv.gz")),
    nrows = 1,
    select = "cell_id"
  )$cell_id[1]
))

## FOV ids are per directory and do not start at 1. See the samplesheet ##
fovs <- if (nzchar(opt$subset_fovs)) {
  as.integer(strsplit(opt$subset_fovs, ";")[[1]])
} else {
  NULL
}

## Negative probes get their own feat type so qc.R reads them directly ##
g <- importCosMx(
  cosmx_dir = opt$data_dir,
  slide = slide,
  fovs = fovs
)$create_gobject(
  load_expression = TRUE,
  load_transcripts = FALSE,
  feat_type = c("rna", "negprobes"),
  split_keyword = list("Negative|SystemControl|FalseCode")
)

## Positions from each cell's own outline, matched by cell ID. Giotto's own centroids are
## paired with IDs by row, and terra::centroids() drops the centre of a one-point outline,
## so every cell after one would carry the next cell's centre ##
centres_by_id <- function(sv) {
  g <- as.data.table(terra::geom(sv))[hole == 0]
  g[,
    `:=`(
      x2 = data.table::shift(x, -1L, fill = x[1L]),
      y2 = data.table::shift(y, -1L, fill = y[1L])
    ),
    by = .(geom, part)
  ]
  g[, cr := x * y2 - x2 * y]
  ct <- g[,
    .(
      a = sum(cr) / 2,
      cx = sum((x + x2) * cr),
      cy = sum((y + y2) * cr),
      mx = mean(x),
      my = mean(y)
    ),
    by = .(geom, part)
  ]
  ## Largest part; an outline without area falls back to the mean of its vertices ##
  ct <- ct[order(geom, -abs(a))][, .SD[1L], by = geom]
  ct[,
    `:=`(
      sdimx = fifelse(abs(a) > 1e-9, cx / (6 * a), mx),
      sdimy = fifelse(abs(a) > 1e-9, cy / (6 * a), my)
    )
  ]
  out <- data.table(
    cell_ID = as.character(sv$poly_ID),
    sdimx = NA_real_,
    sdimy = NA_real_
  )
  out[ct$geom, `:=`(sdimx = ct$sdimx, sdimy = ct$sdimy)]
  out[]
}
old <- getSpatialLocations(g, output = "spatLocsObj")
xy <- centres_by_id(
  getPolygonInfo(g, polygon_name = "cell", return_giottoPolygon = TRUE)[]
)[!is.na(sdimx)]
g <- setSpatialLocations(
  g,
  createSpatLocsObj(
    coordinates = xy[, .(sdimx, sdimy, cell_ID)],
    name = objName(old),
    spat_unit = spatUnit(old),
    provenance = prov(old),
    verbose = FALSE
  ),
  verbose = FALSE,
  initialize = FALSE
)

## A cell without an outline has no position and cannot be kept ##
sl <- as.character(getSpatialLocations(g, output = "data.table")$cell_ID)
cid <- as.character(pDataDT(g)$cell_ID)
lost <- setdiff(cid, sl)
if (length(lost)) {
  message(
    "load.R: ",
    length(lost),
    " cell(s) with no polygon derived position, dropped: ",
    paste(utils::head(lost, 10), collapse = ", ")
  )
  g <- subsetGiotto(g, cell_ids = intersect(cid, sl))
}

## Dont open a plot device through the pipeline ##
instructions(g) <- createGiottoInstructions(
  save_plot = FALSE,
  show_plot = FALSE,
  return_plot = FALSE
)

## Add new sample metadata ##
pd <- pDataDT(g)
g <- addCellMetadata(
  g,
  by_column = TRUE,
  column_cell_ID = "cell_ID",
  new_metadata = data.frame(
    cell_ID = pd$cell_ID,
    sample_id = opt$sample_id,
    patient_id = opt$patient_id,
    treatment = opt$treatment,
    timepoint = opt$timepoint,
    slide_id = opt$slide_id,
    batch = opt$batch
  )
)

## Write ##
saveGiotto(g, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
write.csv(
  data.frame(
    sample_id = opt$sample_id,
    n_cells = nrow(pDataDT(g)),
    n_rna = nrow(fDataDT(g, feat_type = "rna")),
    n_negprobe = nrow(fDataDT(g, feat_type = "negprobes")),
    fovs = if (is.null(fovs)) "all" else paste(fovs, collapse = ";")
  ),
  file.path(opt$outdir, paste0(opt$sample_id, "_load_summary.csv")),
  row.names = FALSE
)
