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
  cell_types = arg("cell_types"),
  sample_col = arg("sample_col", "sample_id"),
  celltype_column = arg("celltype_column", "cell_type"),
  groupvar = arg("groupvar", "niche_type_lam0.8"),
  radius = as.numeric(arg("radius", "415"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(Matrix)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(
  opt$outdir,
  paste0(gsub("[^A-Za-z0-9]+", "_", opt$sample), "_smide_niche_overlap.csv")
)
empty <- data.table(
  sample = character(),
  cell_type = character(),
  niche = character(),
  target = character(),
  avg_neighbor_othercluster = numeric(),
  avg_cluster = numeric(),
  ratio = numeric()
)

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
for (v in c(opt$celltype_column, opt$groupvar)) {
  if (!v %in% names(pd)) {
    stop("smide_overlap.R: no column '", v, "' on the object")
  }
}
ids <- pd$cell_ID[pd[[opt$sample_col]] == opt$sample]
locs <- getSpatialLocations(merged, output = "data.table")
setkey(locs, cell_ID)
raw <- getExpression(merged, values = "raw", output = "matrix")[, ids]

## Only these tables are used below ##
rm(merged)
invisible(gc())

m <- pd[match(ids, pd$cell_ID)]
meta <- data.table(
  cell_ID = ids,
  sdimx = locs[ids]$sdimx,
  sdimy = locs[ids]$sdimy,
  ct = m[[opt$celltype_column]],
  grp = m[[opt$groupvar]]
)
meta <- meta[!is.na(ct) & !is.na(grp)]
if (nrow(meta) < 2L) {
  message("smide_overlap.R: ", opt$sample, " has no labelled cells")
  fwrite(empty, out_file)
  quit(save = "no", status = 0)
}
raw <- raw[, meta$cell_ID]

## smiDE's overlap ratio, restricted to focal cells in one niche: the mean in their
## neighbouring cells of other types (from any niche, so bleed across a niche border
## counts) over the mean in the focal cells. Same neighbours and averaging as smiDE's
## overlap_ratio_metric: radius search, each focal cell's other-type neighbours averaged,
## then averaged over the focal cells (a cell with none contributes zero) ##
xy <- as.matrix(meta[, .(sdimx, sdimy)])
nn <- RANN::nn2(
  xy,
  xy,
  k = min(200L, nrow(xy)),
  searchtype = "radius",
  radius = opt$radius
)$nn.idx
i <- rep(seq_len(nrow(nn)), ncol(nn))
j <- as.vector(nn)
keep <- j > 0 & j != i
adj <- sparseMatrix(i = i[keep], j = j[keep], x = 1, dims = c(nrow(xy), nrow(xy)))

focal <- trimws(strsplit(opt$cell_types, ",")[[1]])
focal <- intersect(focal[nzchar(focal)], unique(meta$ct))
orm <- rbindlist(lapply(focal, function(fct) {
  other <- which(meta$ct != fct)
  rbindlist(lapply(sort(unique(meta$grp[meta$ct == fct])), function(g) {
    f <- which(meta$ct == fct & meta$grp == g)
    w <- adj[f, other, drop = FALSE]
    rs <- Matrix::rowSums(w)
    rs[rs == 0] <- Inf
    w <- Diagonal(x = 1 / rs) %*% w
    nb <- Matrix::rowMeans(raw[, other, drop = FALSE] %*% Matrix::t(w))
    cl <- Matrix::rowMeans(raw[, f, drop = FALSE])
    data.table(
      ct = fct,
      grp = g,
      target = rownames(raw),
      avg_neighbor_othercluster = as.numeric(nb),
      avg_cluster = as.numeric(cl),
      ratio = as.numeric(nb) / as.numeric(cl)
    )
  }))
}))
if (!nrow(orm)) {
  message("smide_overlap.R: ", opt$sample, " has none of the focal cell types")
  fwrite(empty, out_file)
  quit(save = "no", status = 0)
}

## Write ##
fwrite(
  orm[, .(
    sample = opt$sample,
    cell_type = ct,
    niche = grp,
    target,
    avg_neighbor_othercluster,
    avg_cluster,
    ratio
  )],
  out_file
)
message(sprintf("smide_overlap.R: %s, %d rows", opt$sample, nrow(orm)))
