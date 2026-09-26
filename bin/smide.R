#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON")),
  sample = arg("sample"),
  cell_type = arg("cell_type"),
  sample_col = arg("sample_col", "sample_id"),
  celltype_column = arg("celltype_column", "cell_type"),
  groupvar = arg("groupvar", "niche_id_lam0.8"),
  radius = as.numeric(arg("radius", "415")),
  overlap_threshold = as.numeric(arg("overlap_threshold", "0.8")),
  min_detection = as.numeric(arg("min_detection", "0.01")),
  min_cells_per_niche = as.integer(arg("min_cells_per_niche", "30")),
  family = arg("family", "nbinom2"),
  spatial_model = arg("spatial_model", "none,kmeans"),
  k_prop_n = as.numeric(arg("k_prop_n", "0.05")),
  formula = arg("formula", ""),
  comparisons = arg("comparisons", "pairwise,emmeans,one.vs.rest,one.vs.all"),
  cores = as.integer(arg("cores", "1")),
  seed = as.integer(arg("seed", "42"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(Matrix)
  library(smiDE)
})
setDTthreads(1)
set.seed(opt$seed)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}
stub <- paste0(
  gsub("[^A-Za-z0-9]+", "_", opt$sample),
  "__",
  gsub("[^A-Za-z0-9]+", "_", opt$cell_type)
)
out_file <- function(x) file.path(opt$outdir, paste0(stub, "_smide_", x))

## A gate that fails ends the task with a row saying why, not with an error ##
skips <- data.table(
  sample = character(),
  cell_type = character(),
  stage = character(),
  reason = character()
)
skip <- function(stage, reason) {
  skips <<- rbind(
    skips,
    data.table(
      sample = opt$sample,
      cell_type = opt$cell_type,
      stage = stage,
      reason = reason
    )
  )
}
finish <- function() {
  fwrite(skips, out_file("skipped.csv"))
  quit(save = "no", status = 0)
}
try_smide <- function(what, expr) {
  tryCatch(expr, error = function(e) {
    message("smide.R: ", what, " failed: ", conditionMessage(e))
    NULL
  })
}

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
for (v in c(opt$celltype_column, opt$groupvar)) {
  if (!v %in% names(pd)) {
    stop("smide.R: no column '", v, "' on the object")
  }
}
ids <- pd$cell_ID[pd[[opt$sample_col]] == opt$sample]
if (!length(ids)) {
  skip("sample", "no cells for sample")
  finish()
}

## Every cell of the sample. The focal cells' neighbours are read from here ##
m <- pd[match(ids, pd$cell_ID)]
locs <- getSpatialLocations(merged, output = "data.table")
setkey(locs, cell_ID)
raw <- getExpression(merged, values = "raw", output = "matrix")[, ids]
meta <- data.table(
  cell_ID = ids,
  sdimx = locs[ids]$sdimx,
  sdimy = locs[ids]$sdimy,
  smp = opt$sample,
  ct = m[[opt$celltype_column]],
  grp = m[[opt$groupvar]],
  nCount_RNA = pmax(Matrix::colSums(raw), 1)
)
setnames(
  meta,
  c("smp", "ct", "grp"),
  c(opt$sample_col, opt$celltype_column, opt$groupvar)
)
meta <- meta[!is.na(get(opt$celltype_column)) & !is.na(get(opt$groupvar))]
raw <- raw[, meta$cell_ID]

rng <- max(diff(range(meta$sdimx)), diff(range(meta$sdimy)))
if (opt$radius > rng) {
  skip(
    "radius",
    sprintf("radius %.0f exceeds the coordinate range %.0f", opt$radius, rng)
  )
  finish()
}

## Focal cells, in the niches that hold enough of them ##
focal <- meta[get(opt$celltype_column) == opt$cell_type]
if (nrow(focal) < opt$min_cells_per_niche) {
  skip("celltype", sprintf("%d focal cells", nrow(focal)))
  finish()
}
nsz <- table(focal[[opt$groupvar]])
keep <- sort(names(nsz)[nsz >= opt$min_cells_per_niche])
if (length(keep) < 2L) {
  skip(
    "niche",
    sprintf(
      "%d niche(s) with >= %d cells, need 2",
      length(keep),
      opt$min_cells_per_niche
    )
  )
  finish()
}
focal <- focal[get(opt$groupvar) %in% keep]
focal[, (opt$groupvar) := factor(get(opt$groupvar), levels = keep)]

## Genes: little signal from the neighbours, and detected in the focal cells ##
orm <- try_smide(
  "overlap_ratio_metric",
  overlap_ratio_metric(
    assay_matrix = raw,
    metadata = as.data.frame(meta),
    cluster_col = opt$celltype_column,
    cellid_col = "cell_ID",
    sdimx_col = "sdimx",
    sdimy_col = "sdimy",
    radius = opt$radius
  )
)
if (is.null(orm)) {
  skip("overlap", "overlap_ratio_metric errored")
  finish()
}
orm <- as.data.table(orm)[get(opt$celltype_column) == opt$cell_type]
fwrite(
  orm[, `:=`(sample = opt$sample, cell_type = opt$cell_type)],
  out_file("overlap_ratio.csv")
)
det <- Matrix::rowMeans(raw[, focal$cell_ID] > 0)
targets <- intersect(
  orm[ratio <= opt$overlap_threshold, target],
  names(det)[det >= opt$min_detection]
)
message(
  sprintf(
    "smide.R: %d of %d genes pass overlap <= %.2f and detection >= %.3f",
    length(targets),
    nrow(raw),
    opt$overlap_threshold,
    opt$min_detection
  )
)
if (length(targets) < 2L) {
  skip("genes", sprintf("%d genes pass", length(targets)))
  finish()
}

pre <- try_smide(
  "pre_de",
  pre_de(
    metadata = as.data.frame(meta),
    ref_celltype = opt$cell_type,
    cell_type_metadata_colname = opt$celltype_column,
    cellid_colname = "cell_ID",
    sdimx_colname = "sdimx",
    sdimy_colname = "sdimy",
    split_neighbors_by_colname = opt$sample_col,
    mm_radius = opt$radius,
    verbose = FALSE
  )
)
if (is.null(pre)) {
  skip("pre_de", "pre_de errored")
  finish()
}

## smi_de builds otherct_expr from the cells in assay_matrix, so it holds the whole
## sample. Neighbour counts are scaled to the mean library size first ##
assay <- as.matrix(raw[targets, ])
scalefactor <- setNames(mean(meta$nCount_RNA) / meta$nCount_RNA, meta$cell_ID)
fmla <- if (nzchar(opt$formula)) {
  as.formula(opt$formula)
} else {
  as.formula(sprintf(
    "~ RankNorm(otherct_expr) + %s + offset(log(nCount_RNA))",
    opt$groupvar
  ))
}

res <- list()
for (sm in split_arg(opt$spatial_model)) {
  dat <- copy(focal)
  f <- fmla
  if (sm == "kmeans") {
    kc <- as.data.table(xy_kmeans_clusters(
      metadata = as.data.frame(dat),
      x_coord_col = "sdimx",
      y_coord_col = "sdimy",
      cluster_name = "k_cluster",
      k_prop_n = opt$k_prop_n,
      seed = opt$seed
    ))
    dat[, k_cluster := factor(kc$k_cluster[match(cell_ID, kc$cell_ID)])]
    f <- update(fmla, ~ . + (1 | k_cluster))
  } else if (sm != "none") {
    stop("smide.R: unknown spatial_model '", sm, "', use none or kmeans")
  }

  t0 <- Sys.time()
  fit <- try_smide(
    paste0("smi_de [", sm, "]"),
    smi_de(
      assay_matrix = assay,
      metadata = as.data.frame(dat),
      formula = f,
      pre_de_obj = pre,
      groupVar = opt$groupvar,
      groupVar_levels = keep,
      nCores = opt$cores,
      family = opt$family,
      targets = targets,
      cellid_colname = "cell_ID",
      neighbor_expr_cell_type_metadata_colname = opt$celltype_column,
      neighbor_expr_overlap_agg = "sum",
      neighbor_expr_totalcount_normalize = TRUE,
      neighbor_expr_totalcount_scalefactor = scalefactor
    )
  )
  if (is.null(fit)) {
    skip("fit", paste0("smi_de failed, spatial_model = ", sm))
    next
  }
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  r <- try_smide(
    paste0("results [", sm, "]"),
    results(
      fit,
      comparisons = split_arg(opt$comparisons),
      variable = opt$groupvar
    )
  )
  if (is.null(r)) {
    skip("results", paste0("results() failed, spatial_model = ", sm))
    next
  }
  d <- rbindlist(
    lapply(names(r), function(nm) {
      if (nrow(r[[nm]])) r[[nm]][, result_component := nm]
    }),
    fill = TRUE,
    use.names = TRUE
  )
  if (!nrow(d)) {
    skip("results", paste0("every comparison empty, spatial_model = ", sm))
    next
  }
  d[, `:=`(
    sample = opt$sample,
    cell_type = opt$cell_type,
    spatial_model = sm,
    fit_secs = secs
  )]
  if ("p.value" %in% names(d)) {
    d[, fdr := p.adjust(p.value, "BH"), by = result_component]
  }
  res[[sm]] <- d
}
if (!length(res)) {
  skip("fit", "no model produced results")
  finish()
}

## Write ##
de <- rbindlist(res, fill = TRUE, use.names = TRUE)
fwrite(de, out_file("de_results.csv"))
writeLines(setdiff(targets, de$target), out_file("failed_genes.txt"))
message(sprintf("smide.R: wrote %d rows", nrow(de)))
finish()
