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
  reduction = arg("reduction", "harmony"),
  reduction_file = arg("reduction_file", ""),
  dims_use = arg("dims_use", "1:30"),
  knn_k = as.integer(arg("knn_k", "20")),
  batch_column = arg("batch_column", "slide_id"),
  call_top = as.numeric(arg("call_top", "0.9")),
  call_second = as.numeric(arg("call_second", "0.5")),
  seed = as.integer(arg("seed", "42"))
)
if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(HieraType)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

parse_dims <- function(s) {
  p <- strsplit(s, ":")[[1]]
  if (length(p) == 2) {
    as.integer(p[1]):as.integer(p[2])
  } else {
    as.integer(strsplit(s, ",")[[1]])
  }
}
dims <- parse_dims(opt$dims_use)

## Default to embedding from integrate.R, else use the first ##
red <- if (nzchar(opt$reduction_file) && file.exists(opt$reduction_file)) {
  trimws(readLines(opt$reduction_file)[1])
} else {
  opt$reduction
}

## Index markers define a class, predictors refine it. Both are intersected with the panel ##
ht_classes <- list(
  bcell = list(
    i = c("MS4A1", "CD79A"),
    p = c(
      "CD79B",
      "CD19",
      "BANK1",
      "PAX5",
      "FCRLA",
      "VPREB3",
      "TCL1A",
      "CR2",
      "FCER2",
      "CD22",
      "CD24",
      "BLK",
      "TNFRSF13B",
      "TNFRSF13C",
      "IGHM",
      "IGHD"
    )
  ),
  plasma = list(
    i = c("MZB1", "DERL3", "TNFRSF17"),
    p = c(
      "JCHAIN",
      "XBP1",
      "PRDM1",
      "CD38",
      "SDC1",
      "IGKC",
      "IGHA1",
      "POU2AF1",
      "SEC11C"
    )
  ),
  tcell = list(
    i = c("CD3D", "CD3E", "CD3G"),
    p = c(
      "CD2",
      "TRAC",
      "IL7R",
      "CD8A",
      "CD8B",
      "CD4",
      "CD40LG",
      "GZMK",
      "CCL5",
      "GZMA",
      "LTB"
    )
  ),
  nk = list(
    i = c("NKG7", "GNLY", "KLRD1"),
    p = c(
      "NCAM1",
      "GZMB",
      "FCGR3A"
    )
  ),
  myeloid = list(
    i = c("CD68", "C1QA"),
    p = c(
      "C1QB",
      "C1QC",
      "CD163",
      "LYZ",
      "CD14",
      "ITGAX",
      "FCN1",
      "S100A8",
      "S100A9",
      "VCAN",
      "LYVE1",
      "CD1C",
      "FCER1A"
    )
  ),
  epithelial = list(
    i = c("PAX8", "LRP2"),
    p = c(
      "CUBN",
      "SLC34A1",
      "ALDOB",
      "GPX3",
      "UMOD",
      "SLC12A1",
      "AQP2",
      "AQP3",
      "SLC4A1",
      "CA2",
      "SLC12A3",
      "CALB1",
      "HNF4A"
    )
  ),
  endothelial = list(
    i = c("PECAM1", "EMCN"),
    p = c(
      "CD34",
      "KDR",
      "FLT1",
      "VWF",
      "CLDN5",
      "ENG",
      "PLVAP"
    )
  ),
  stromal = list(
    i = c("PDGFRB", "COL1A1"),
    p = c(
      "PDGFRA",
      "COL3A1",
      "DCN",
      "LUM",
      "ACTA2",
      "TAGLN",
      "RGS5",
      "NOTCH3",
      "POSTN"
    )
  )
)

## Exclusions ##
ddp <- list(
  bcell = "plasma",
  plasma = "bcell",
  epithelial = c("tcell", "bcell", "plasma", "myeloid"),
  endothelial = c("tcell", "bcell", "plasma", "myeloid"),
  stromal = c("tcell", "bcell", "plasma", "myeloid")
)

merged <- loadGiotto(opt$input)
panel <- fDataDT(merged, feat_type = "rna")$feat_ID
counts <- t(getExpression(
  merged,
  feat_type = "rna",
  values = "raw",
  output = "matrix"
))
md <- pDataDT(merged, feat_type = "rna")

if (!opt$batch_column %in% names(md)) {
  stop("hieratype.R: batch column '", opt$batch_column, "' not in metadata")
}
if (length(unique(md[[opt$batch_column]])) < 2L) {
  message(
    "hieratype.R: batch column '",
    opt$batch_column,
    "' is constant, fitting as one batch"
  )
}

idx <- lapply(ht_classes, function(x) intersect(x$i, panel))
prd <- lapply(ht_classes, function(x) intersect(x$p, panel))

## A class with no index marker on the panel can never be called, so fail loudly ##
write.csv(
  data.frame(
    class = names(idx),
    n_index = lengths(idx),
    n_predictor = lengths(prd),
    dropped = vapply(
      ht_classes,
      function(x) paste(setdiff(c(x$i, x$p), panel), collapse = ";"),
      character(1L)
    )
  ),
  file.path(opt$outdir, "hieratype_markers.csv"),
  row.names = FALSE
)
if (any(lengths(idx) == 0L)) {
  stop(
    "hieratype.R: no index marker on the panel for class(es): ",
    paste(names(idx)[lengths(idx) == 0L], collapse = ", ")
  )
}

mkl <- HieraType::make_markerslist(index_marker = idx, predictors = prd)

## HieraType builds its own graph and does not reuse the clustering one ##
emb <- getDimReduction(
  merged,
  feat_type = "rna",
  reduction_method = red,
  name = red,
  output = "matrix"
)
adj <- HieraType::jaccard_adjacency_matrix(
  cell_embeddings = emb,
  k_nearest_neighbors = opt$knn_k,
  npcs = max(dims)
)

mg <- HieraType::fit_metagene_scores(
  markerslist = mkl,
  counts_matrix = counts,
  adjacency_matrix = adj$snn,
  obs = as.data.frame(md),
  batch_variable = opt$batch_column,
  cellid_colname = "cell_ID"
)

set.seed(opt$seed)
cm <- HieraType::cluster_metagenes(
  metagenes = mg,
  to_model = "yhat",
  discourage_double_positive = ddp,
  seed = opt$seed
)

pp <- cm$post_probs
CLS <- setdiff(names(pp), c("cell_ID", "best_class", "best_score"))
P <- as.matrix(pp[, ..CLS])
top2 <- apply(P, 1, function(x) sort(x, decreasing = TRUE)[1:2])

## A confident call needs a clear winner and a distant runner up. Everything else is unknown ##
pp[,
  ht_call := ifelse(
    top2[1, ] >= opt$call_top & top2[2, ] <= opt$call_second,
    CLS[max.col(P)],
    "unknown"
  )
]

calls <- data.frame(
  cell_ID = pp$cell_ID,
  ht_class = pp$best_class,
  ht_score = pp$best_score,
  ht_call = pp$ht_call
)
merged <- addCellMetadata(
  merged,
  by_column = TRUE,
  column_cell_ID = "cell_ID",
  new_metadata = calls
)

## Write ##
write.csv(
  calls,
  file.path(opt$outdir, "hieratype_calls.csv"),
  row.names = FALSE
)

n_called <- sum(pp$ht_call != "unknown")
message(
  "hieratype.R: ",
  n_called,
  " of ",
  nrow(pp),
  " cells called (",
  round(100 * n_called / nrow(pp), 1),
  "%)"
)

tab <- table(pp$ht_call)
write.csv(
  data.frame(
    ht_call = names(tab),
    n_cells = as.integer(tab),
    pct = round(100 * as.integer(tab) / nrow(pp), 2)
  ),
  file.path(opt$outdir, "hieratype_summary.csv"),
  row.names = FALSE
)

saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
