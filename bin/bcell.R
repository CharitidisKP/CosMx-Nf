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
  celltype_column = arg("celltype_column", "cell_type"),
  bcell_regex = arg("bcell_regex"),
  min_cells = as.integer(arg("min_cells", "30")),
  markers = arg("markers"),
  n_hvgs = as.integer(arg("n_hvgs", "4000")),
  n_pcs = as.integer(arg("n_pcs", "50")),
  dims_use = arg("dims_use", "1:20"),
  batch_column = arg("batch_column"),
  batch_correct = arg("batch_correct", "auto"),
  knn_k = as.integer(arg("knn_k", "10")),
  resolution = as.numeric(arg("resolution", "0.4")),
  n_iterations = as.integer(arg("n_iterations", "1000")),
  top_n = as.integer(arg("top_n", "25")),
  sample_col = arg("sample_col", "sample_id"),
  seed = as.integer(arg("seed", "42"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(ggplot2)
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

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
b_ids <- pd$cell_ID[grepl(opt$bcell_regex, pd[[opt$celltype_column]])]
summ <- data.frame(
  n_cells = length(b_ids),
  min_cells = opt$min_cells,
  reduction = NA_character_,
  n_subclusters = NA_integer_
)
if (length(b_ids) < opt$min_cells) {
  message(
    "bcell.R: ",
    length(b_ids),
    " B-lineage cells < min_cells ",
    opt$min_cells,
    ", not subclustered"
  )
  write.csv(
    summ,
    file.path(opt$outdir, "bcell_summary.csv"),
    row.names = FALSE
  )
  quit(save = "no", status = 0)
}

## Own HVGs and PCA, as subcluster.R does per cluster. The global embedding is
## shaped by the abundant types, not by what varies inside this lineage ##
gb <- subsetGiotto(merged, cell_ids = b_ids)
nm <- getExpression(
  gb,
  feat_type = "rna",
  values = "normalized",
  output = "matrix"
)
dec <- suppressWarnings(scran::modelGeneVar(nm))
hvf <- suppressWarnings(scran::getTopHVGs(dec, n = opt$n_hvgs))
feats <- union(hvf, intersect(unique(read.csv(opt$markers)$gene), rownames(nm)))
rm(nm)
ncp <- min(opt$n_pcs, length(feats) - 1L, length(b_ids) - 1L)
dims <- parse_dims(opt$dims_use)
dims <- dims[dims <= ncp]
gb <- runPCA(
  gb,
  feat_type = "rna",
  feats_to_use = feats,
  ncp = ncp,
  name = "pca_b"
)

## Same batch switch as integrate.R, re-tested on these cells ##
bd <- pDataDT(gb)
bvars <- trimws(strsplit(opt$batch_column, ",")[[1]])
use <- Filter(
  function(v) v %in% names(bd) && length(unique(bd[[v]])) >= 2L,
  bvars[nzchar(bvars)]
)
do_bc <- length(use) > 0L && tolower(opt$batch_correct) != "false"
red <- c(type = "pca", name = "pca_b")
if (do_bc) {
  h <- tryCatch(
    runGiottoHarmony(
      gb,
      vars_use = use,
      dim_reduction_to_use = "pca",
      dim_reduction_name = "pca_b",
      dimensions_to_use = dims,
      name = "harmony_b",
      seed_number = opt$seed
    ),
    error = function(e) {
      message("bcell.R: harmony failed (", conditionMessage(e), "), using PCA")
      NULL
    }
  )
  hit <- if (is.null(h)) {
    character()
  } else {
    Filter(
      function(n) {
        !is.null(tryCatch(
          getDimReduction(
            h,
            reduction = "cells",
            reduction_method = "harmony",
            name = n,
            output = "matrix"
          ),
          error = function(e) NULL
        ))
      },
      c("harmony_b", "harmony")
    )
  }
  if (length(hit)) {
    gb <- h
    red <- c(type = "harmony", name = hit[[1]])
  }
}

gb <- createNearestNetwork(
  gb,
  type = "sNN",
  dim_reduction_to_use = red[["type"]],
  dim_reduction_name = red[["name"]],
  dimensions_to_use = dims,
  k = min(opt$knn_k, length(b_ids) - 1L),
  name = "NN.bcell"
)
gb <- doLeidenCluster(
  gb,
  nn_network_to_use = "sNN",
  network_name = "NN.bcell",
  resolution = opt$resolution,
  n_iterations = opt$n_iterations,
  name = "bcell_sub",
  seed_number = opt$seed
)
gb <- runUMAP(
  gb,
  dim_reduction_to_use = red[["type"]],
  dim_reduction_name = red[["name"]],
  dimensions_to_use = dims,
  name = "umap_bcell",
  seed_number = opt$seed
)
de <- suppressWarnings(findMarkers_one_vs_all(
  gb,
  method = "scran",
  expression_values = "normalized",
  cluster_column = "bcell_sub",
  min_feats = 5
))

## Write ##
gpd <- pDataDT(gb)
summ$reduction <- red[["name"]]
summ$n_subclusters <- uniqueN(gpd$bcell_sub)
fwrite(
  gpd[,
    c("cell_ID", opt$sample_col, opt$celltype_column, "bcell_sub"),
    with = FALSE
  ],
  file.path(opt$outdir, "bcell_subcluster_cells.csv")
)
fwrite(
  gpd[, .(n_cells = .N), by = .(bcell_sub, sample = get(opt$sample_col))],
  file.path(opt$outdir, "bcell_subcluster_by_sample.csv")
)
fwrite(
  de[, head(.SD[order(ranking)], opt$top_n), by = cluster],
  file.path(opt$outdir, "bcell_subcluster_markers.csv")
)
p <- plotUMAP(
  gb,
  dim_reduction_name = "umap_bcell",
  cell_color = "bcell_sub",
  point_size = 1,
  show_plot = FALSE,
  return_plot = TRUE
)
ggsave(
  file.path(opt$outdir, "umap_bcell_subclusters.png"),
  p,
  width = 7,
  height = 6,
  dpi = 150,
  bg = "white"
)
write.csv(summ, file.path(opt$outdir, "bcell_summary.csv"), row.names = FALSE)
