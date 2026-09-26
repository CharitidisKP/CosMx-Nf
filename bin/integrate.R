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
  markers = arg("markers"),
  scalefactor = arg("scalefactor", "auto"),
  n_hvgs = as.integer(arg("n_hvgs", "4000")),
  n_pcs = as.integer(arg("n_pcs", "50")),
  batch_column = arg("batch_column"),
  batch_correct = arg("batch_correct", "auto"),
  dims_use = arg("dims_use", "1:30"),
  seed = as.integer(arg("seed", "42"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(scran)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

parse_dims <- function(s) {
  p <- strsplit(s, ":")[[1]]
  if (length(p) == 2) {
    as.integer(p[1]):as.integer(p[2])
  } else {
    as.integer(strsplit(s, ",")[[1]])
  }
}
dims <- parse_dims(opt$dims_use)

merged <- loadGiotto(opt$input)

## Make sure features are used downstream, not negative probes ##
raw <- getExpression(
  merged,
  feat_type = "rna",
  values = "raw",
  output = "matrix"
)
lib_med <- round(median(Matrix::colSums(raw)))
sf <- if (opt$scalefactor == "auto") lib_med else as.numeric(opt$scalefactor)

merged <- normalizeGiotto(merged, feat_type = "rna", scalefactor = sf)
merged <- addStatistics(merged, feat_type = "rna")

## Top HVGs unioned with the canonical panel, so markers cannot be ranked out ##
nm <- as.matrix(
  getExpression(
    merged,
    feat_type = "rna",
    values = "normalized",
    output = "matrix"
  )
)
dec <- suppressWarnings(scran::modelGeneVar(nm))
hvf <- suppressWarnings(scran::getTopHVGs(dec, n = opt$n_hvgs))
canon <- unique(read.csv(opt$markers)$gene)
pca_feats <- union(hvf, intersect(canon, rownames(nm)))
rm(nm)

merged <- runPCA(
  merged,
  feat_type = "rna",
  feats_to_use = pca_feats,
  ncp = opt$n_pcs,
  name = "pca"
)

## batch_column is a comma separated list. Drop any entry that is absent or constant ##
pd <- pDataDT(merged)
bvars <- trimws(strsplit(opt$batch_column, ",")[[1]])
bvars <- bvars[nzchar(bvars)]

miss <- setdiff(bvars, names(pd))
if (length(miss)) {
  message(
    "integrate.R: batch column(s) absent, dropped: ",
    paste(miss, collapse = ", ")
  )
}
have <- intersect(bvars, names(pd))
levs <- vapply(have, function(v) length(unique(pd[[v]])), integer(1L))
use <- have[levs >= 2L]
if (length(use) < length(have)) {
  message(
    "integrate.R: batch column(s) with < 2 levels, dropped: ",
    paste(setdiff(have, use), collapse = ", ")
  )
}

## A solo object has one level of everything, so it falls back to pca.
## cluster.R reads reduction.txt rather than assuming harmony ##
do_bc <- switch(
  tolower(opt$batch_correct),
  "true" = TRUE,
  "false" = FALSE,
  length(use) > 0L
)
if (do_bc && !length(use)) {
  message("integrate.R: no usable batch covariate, reduction = pca")
  do_bc <- FALSE
}
red <- if (do_bc) "harmony" else "pca"
if (do_bc) {
  merged <- runGiottoHarmony(
    merged,
    vars_use = use,
    dim_reduction_to_use = "pca",
    dim_reduction_name = "pca",
    dimensions_to_use = dims,
    name = "harmony",
    seed_number = opt$seed
  )
}
writeLines(red, file.path(opt$outdir, "reduction.txt"))

## Write ##
p <- screePlot(
  merged,
  dim_reduction_name = "pca",
  ncp = opt$n_pcs,
  save_plot = FALSE,
  return_plot = TRUE
)
ggsave(
  file.path(opt$outdir, "scree.png"),
  p,
  width = 7,
  height = 4,
  dpi = 150,
  bg = "white"
)
writeLines(pca_feats, file.path(opt$outdir, "pca_feats.txt"))
write.csv(
  data.frame(
    scalefactor = sf,
    lib_median = lib_med,
    n_hvg = length(hvf),
    n_pca_feats = length(pca_feats),
    batch_requested = opt$batch_column,
    batch_used = if (length(use)) paste(use, collapse = ",") else NA_character_,
    batch_levels = if (length(use)) {
      paste(levs[use], collapse = ",")
    } else {
      NA_character_
    },
    batch_corrected = do_bc,
    reduction = red,
    dims = paste0(min(dims), ":", max(dims))
  ),
  file.path(opt$outdir, "integrate_summary.csv"),
  row.names = FALSE
)
saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
