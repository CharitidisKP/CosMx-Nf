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
  sample_col = arg("sample_col", "sample_id"),
  k = as.integer(arg("k", "30"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(Matrix)
  library(InSituCor)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
raw <- getExpression(merged, values = "raw", output = "matrix")
ids <- colnames(raw)
m <- pd[match(ids, pd$cell_ID)]
locs <- getSpatialLocations(merged, output = "data.table")
setkey(locs, cell_ID)
xy <- as.matrix(locs[ids, .(sdimx, sdimy)])
rownames(xy) <- ids

## InSituCor wants cells x genes, scaled to the same total per cell ##
total <- Matrix::colSums(raw)
counts <- Matrix::t(raw) * (mean(total) / pmax(total, 1))

## Conditioned on cell type, tissue, signal and background, per the package docs ##
conditionon <- data.frame(
  celltype = m[[opt$celltype_column]],
  tissue = m[[opt$sample_col]],
  totalcounts = total,
  negmean = m$neg_mean,
  row.names = ids
)
res <- insitucor(
  counts = counts,
  conditionon = conditionon,
  celltype = m[[opt$celltype_column]],
  xy = xy,
  k = opt$k,
  tissue = m[[opt$sample_col]]
)

## Write ##
fwrite(
  as.data.table(res$modules),
  file.path(opt$outdir, "insitucor_modules.csv")
)
write.csv(
  as.data.frame(res$celltypeinvolvement),
  file.path(opt$outdir, "insitucor_celltype_involvement.csv")
)
fwrite(
  data.table(cell_ID = ids, as.matrix(res$scores_env)),
  file.path(opt$outdir, "insitucor_scores_env.csv")
)
