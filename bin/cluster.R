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
  knn_k = as.integer(arg("knn_k", "15")),
  n_iterations = as.integer(arg("n_iterations", "1000")),
  leiden_res_ladder = arg("leiden_res_ladder"),
  umap_resolutions = arg("umap_resolutions"),
  dims_use = arg("dims_use", "1:30"),
  reduction = arg("reduction", "harmony"),
  reduction_file = arg("reduction_file", ""),
  seed = as.integer(arg("seed", "42"))
)
if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
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
dims <- parse_dims(opt$dims_use)
res_tok <- trimws(strsplit(opt$leiden_res_ladder, ",")[[1]])
## Every resolution is clustered, but only the chosen ones get a UMAP. Empty: all ##
umap_tok <- trimws(strsplit(opt$umap_resolutions, ",")[[1]])
umap_tok <- if (length(umap_tok)) intersect(res_tok, umap_tok) else res_tok

## Which embedding to cluster on from integrate.R. --reduction is the fallback. ##
red <- if (nzchar(opt$reduction_file) && file.exists(opt$reduction_file)) {
  trimws(readLines(opt$reduction_file)[1])
} else {
  opt$reduction
}
nn <- paste0("NN.", red)

merged <- loadGiotto(opt$input)

## One shared nearest neighbour graph, reused by every resolution ##
merged <- createNearestNetwork(
  merged,
  feat_type = "rna",
  type = "sNN",
  dim_reduction_to_use = red,
  dim_reduction_name = red,
  dimensions_to_use = dims,
  k = opt$knn_k,
  name = nn
)

## Every resolution is kept. The choice between them is made downstream ##
for (rt in res_tok) {
  merged <- doLeidenCluster(
    merged,
    feat_type = "rna",
    nn_network_to_use = "sNN",
    network_name = nn,
    resolution = as.numeric(rt),
    n_iterations = opt$n_iterations,
    name = paste0("leiden_clus_res", rt),
    seed_number = opt$seed
  )
}

merged <- runUMAP(
  merged,
  feat_type = "rna",
  dim_reduction_to_use = red,
  dim_reduction_name = red,
  dimensions_to_use = dims,
  name = "umap",
  seed_number = opt$seed
)
if (red != "pca") {
  merged <- runUMAP(
    merged,
    feat_type = "rna",
    dim_reduction_to_use = "pca",
    dim_reduction_name = "pca",
    dimensions_to_use = dims,
    name = "umap_pca",
    seed_number = opt$seed
  )
}

save_umap <- function(colr, nm, umap_name = "umap") {
  p <- plotUMAP(
    merged,
    feat_type = "rna",
    dim_reduction_name = umap_name,
    cell_color = colr,
    point_size = 0.6,
    show_plot = FALSE,
    return_plot = TRUE
  )
  ggsave(
    file.path(opt$outdir, nm),
    p,
    width = 7,
    height = 6,
    dpi = 150,
    bg = "white"
  )
}

## Same covariates on both embeddings. This pair is the only evidence Harmony did anything ##
for (v in c("sample_id", "batch", "slide_id", "treatment")) {
  if (red != "pca") {
    save_umap(v, paste0("umap_pca_by_", v, ".png"), "umap_pca")
  }
  save_umap(v, paste0("umap_by_", v, ".png"))
}
for (rt in umap_tok) {
  save_umap(
    paste0("leiden_clus_res", rt),
    paste0("umap_leiden_res", rt, ".png")
  )
}

## Write ##
md <- pDataDT(merged, feat_type = "rna")
for (rt in res_tok) {
  col <- paste0("leiden_clus_res", rt)
  mix <- round(
    prop.table(
      table(cluster = md[[col]], sample_id = md$sample_id),
      margin = 1
    ),
    3
  )
  write.csv(
    as.data.frame.matrix(mix),
    file.path(opt$outdir, paste0("sample_mix_res", rt, ".csv"))
  )
}
write.csv(
  data.frame(
    resolution = res_tok,
    reduction = red,
    n_clusters = sapply(
      res_tok,
      function(rt) length(unique(md[[paste0("leiden_clus_res", rt)]]))
    )
  ),
  file.path(opt$outdir, "cluster_summary.csv"),
  row.names = FALSE
)
saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
