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
  sample = arg("sample"),
  sample_col = arg("sample_col", "sample_id"),
  juxta_thr = as.numeric(arg("juxta_thr", "150")),
  para_ls = arg("para_ls", "400,800,1600"),
  min_cells = as.integer(arg("min_cells", "50")),
  cores = as.integer(arg("cores", "1")),
  seed = as.integer(arg("seed", "42")),
  reuse_dir = arg("reuse_dir")
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(mistyR)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

## Views and models run in parallel only under a future plan ##
future::plan(future::multisession, workers = opt$cores)

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
norm <- getExpression(merged, values = "normalized", output = "matrix")
locs <- getSpatialLocations(merged, output = "data.table")
setkey(locs, cell_ID)

## Only these tables are used below. The object would hold memory for the whole run ##
rm(merged)
invisible(gc())
para_ls <- as.numeric(strsplit(opt$para_ls, ",")[[1]])

## Targets become formula terms, so IGHG1/2 and the like need syntactic names ##
feats <- intersect(unique(read.csv(opt$markers)$gene), rownames(norm))
gene_of <- setNames(feats, make.names(feats))

## One joint model for this sample: intra, juxta and every paraview scale ##
ids <- pd$cell_ID[pd[[opt$sample_col]] == opt$sample]
if (length(ids) < opt$min_cells) {
  message(
    "cci_misty.R: ",
    opt$sample,
    " has ",
    length(ids),
    " cells < min_cells, skipped"
  )
  quit(save = "no", status = 0)
}
ex <- t(as.matrix(norm[feats, ids, drop = FALSE]))
colnames(ex) <- names(gene_of)
rm(norm)
invisible(gc())

## run_misty stops on a target with no variance ##
ex <- as.data.frame(ex[, matrixStats::colVars(ex) > 0, drop = FALSE])
f <- file.path(
  opt$outdir,
  paste0("misty_", gsub("[^A-Za-z0-9]+", "_", opt$sample))
)

## A folder from an earlier run is copied instead of refitted when it is finished with these
## settings: one performance row per target and an importance file per target for each view
## this run asks for, named as mistyR names them (intra, juxta.<thr>, para.<l>) ##
prev <- if (nzchar(opt$reuse_dir)) file.path(opt$reuse_dir, basename(f)) else ""
views_ab <- c("intra", paste0("juxta.", opt$juxta_thr), paste0("para.", para_ls))
imp <- if (nzchar(prev)) list.files(prev, "^importances_") else character()
finished <- nzchar(prev) &&
  file.exists(file.path(prev, "performance.txt")) &&
  nrow(fread(file.path(prev, "performance.txt"))) == ncol(ex) &&
  length(imp) == ncol(ex) * length(views_ab) &&
  all(vapply(views_ab, function(v) sum(endsWith(imp, paste0("_", v, ".txt"))), 0L) == ncol(ex))
if (finished) {
  message("cci_misty.R: ", opt$sample, " reused from ", prev)
  invisible(file.copy(prev, opt$outdir, recursive = TRUE))
} else {
  xy <- as.data.frame(locs[ids, .(x = sdimx, y = sdimy)])
  views <- create_initial_view(ex) |>
    add_juxtaview(xy, neighbor.thr = opt$juxta_thr, verbose = FALSE)
  for (l in para_ls) {
    views <- add_paraview(
      views,
      xy,
      l = l,
      zoi = opt$juxta_thr,
      verbose = FALSE
    )
  }
  invisible(run_misty(views, results.folder = f, seed = opt$seed))
}

## Write ##
## The key takes the collect step from the folder back to the sample, and from each
## syntactic target name back to the gene ##
fwrite(
  data.table(
    folder = basename(f),
    sample = opt$sample,
    name = names(gene_of),
    gene = unname(gene_of)
  ),
  paste0(f, "_key.csv")
)
