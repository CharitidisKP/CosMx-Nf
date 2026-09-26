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
  n_hvgs = as.integer(arg("n_hvgs", "4000")),
  k_geom = arg("k_geom", "15,30"),
  lambda_ladder = arg("lambda_ladder", "0.2,0.8"),
  use_agf = as.logical(arg("use_agf", "true")),
  n_pcs = as.integer(arg("n_pcs", "50")),
  k_nn = as.integer(arg("k_nn", "30")),
  banksy_res = as.numeric(arg("banksy_res", "0.4")),
  batch_column = arg("batch_column"),
  sample_col = arg("sample_col", "sample_id"),
  celltype_column = arg("celltype_column", "cell_type"),
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
  library(SpatialExperiment)
  library(Banksy)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}

## k_geom is one k per harmonic: the mean, then the AGF ##
k_geom <- as.integer(split_arg(opt$k_geom))
lambdas <- as.numeric(split_arg(opt$lambda_ladder))
if (!length(k_geom) || anyNA(k_geom) || !length(lambdas) || anyNA(lambdas)) {
  stop("banksy.R: --k_geom and --lambda_ladder must be comma separated numbers")
}
if (!opt$use_agf) {
  k_geom <- k_geom[1]
}

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)

## Same features as the discovery PCA: top HVGs unioned with the canonical panel ##
norm <- getExpression(merged, values = "normalized", output = "matrix")
dec <- suppressWarnings(scran::modelGeneVar(norm))
hvf <- suppressWarnings(scran::getTopHVGs(dec, n = opt$n_hvgs))
canon <- intersect(unique(read.csv(opt$markers)$gene), rownames(norm))
norm <- norm[union(hvf, canon), ]

## Neighbourhoods are built per sample, so no cell borrows from another tissue ##
locs <- getSpatialLocations(merged, output = "data.table")
setkey(locs, cell_ID)
samples <- sort(unique(pd[[opt$sample_col]]))
spes <- lapply(samples, function(s) {
  ids <- pd$cell_ID[pd[[opt$sample_col]] == s]
  xy <- as.matrix(locs[ids, .(sdimx, sdimy)])
  rownames(xy) <- ids
  spe <- SpatialExperiment(
    assays = list(normalized = norm[, ids, drop = FALSE]),
    spatialCoords = xy,
    sample_id = s
  )
  computeBanksy(
    spe,
    assay_name = "normalized",
    compute_agf = opt$use_agf,
    k_geom = k_geom,
    seed = opt$seed,
    verbose = FALSE
  )
})
aug <- do.call(cbind, spes)
rm(spes, norm)
aug <- runBanksyPCA(
  aug,
  use_agf = opt$use_agf,
  lambda = lambdas,
  npcs = opt$n_pcs,
  seed = opt$seed
)

## Harmony on the same batch columns as discovery, minus any absent or constant ##
bvars <- intersect(split_arg(opt$batch_column), names(pd))
bvars <- bvars[vapply(bvars, function(v) uniqueN(pd[[v]]) > 1L, logical(1L))]
message(
  "banksy.R: Harmony on ",
  if (length(bvars)) paste(bvars, collapse = ", ") else "nothing"
)
bmeta <- as.data.frame(pd[match(colnames(aug), cell_ID), ..bvars])

## One niche column per lambda. 0.2 follows cell type, 0.8 follows tissue domain ##
niches <- data.table(cell_ID = colnames(aug))
for (lam in lambdas) {
  emb <- reducedDim(aug, sprintf("PCA_M%d_lam%s", as.integer(opt$use_agf), lam))
  if (length(bvars)) {
    set.seed(opt$seed)
    emb <- harmony::RunHarmony(
      emb,
      meta_data = bmeta,
      vars_use = bvars,
      verbose = FALSE
    )
  }
  red <- paste0("banksy_lam", lam)
  reducedDim(aug, red) <- emb
  before <- names(colData(aug))
  aug <- clusterBanksy(
    aug,
    dimred = red,
    algo = "leiden",
    k_neighbors = opt$k_nn,
    resolution = opt$banksy_res,
    seed = opt$seed
  )
  cl <- setdiff(names(colData(aug)), before)
  niches[[paste0("niche_id_lam", lam)]] <- paste0("niche_", colData(aug)[[cl]])
}
merged <- addCellMetadata(
  merged,
  by_column = TRUE,
  column_cell_ID = "cell_ID",
  new_metadata = as.data.frame(niches)
)

## Write ##
pdc <- pDataDT(merged)
cols <- setdiff(names(niches), "cell_ID")
long <- rbindlist(lapply(cols, function(cl) {
  pdc[, .(
    lambda = sub("niche_id_lam", "", cl),
    sample = get(opt$sample_col),
    cell_type = get(opt$celltype_column),
    niche = get(cl)
  )]
}))
fwrite(
  long[, .(n_cells = .N), by = .(lambda, sample, niche)],
  file.path(opt$outdir, "niche_by_sample.csv")
)
fwrite(
  long[, .(n_cells = .N), by = .(lambda, niche, cell_type)],
  file.path(opt$outdir, "niche_by_celltype.csv")
)
comp <- long[, .(n = .N), by = .(lambda, niche, cell_type)]
comp[, share := n / sum(n), by = .(lambda, niche)]
for (l in unique(comp$lambda)) {
  d <- comp[lambda == l]
  d[,
    niche := factor(
      niche,
      levels = unique(niche[order(as.numeric(sub("^niche_", "", niche)))])
    )
  ]
  h <- ggplot(d, aes(cell_type, niche, fill = share)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1)) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Share",
      title = paste0("Cell types per niche, lambda ", l)
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(
    file.path(opt$outdir, paste0("niche_celltype_heatmap_lam", l, ".png")),
    h,
    width = max(6, 2 + 0.35 * uniqueN(d$cell_type)),
    height = max(4, 1.5 + 0.3 * uniqueN(d$niche)),
    dpi = 150,
    bg = "white"
  )
}
write.csv(
  data.frame(
    lambda = sub("niche_id_lam", "", cols),
    n_niches = vapply(cols, function(cl) uniqueN(pdc[[cl]]), integer(1L)),
    n_feats = nrow(aug),
    k_geom = paste(k_geom, collapse = ","),
    use_agf = opt$use_agf,
    batch_used = paste(bvars, collapse = ",")
  ),
  file.path(opt$outdir, "banksy_summary.csv"),
  row.names = FALSE
)
for (s in samples) {
  gs <- subsetGiotto(merged, cell_ids = pdc$cell_ID[pdc[[opt$sample_col]] == s])
  for (cl in cols) {
    p <- spatInSituPlotPoints(
      gs,
      show_polygon = TRUE,
      polygon_feat_type = "cell",
      polygon_fill = cl,
      polygon_fill_as_factor = TRUE,
      polygon_line_size = 0.05,
      feats = NULL,
      show_plot = FALSE,
      return_plot = TRUE
    )
    ggsave(
      file.path(opt$outdir, paste0("niche_spatial_", s, "_", cl, ".png")),
      p + labs(title = paste(s, cl)),
      width = 7,
      height = 7,
      dpi = 150,
      bg = "white"
    )
  }
}
saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
