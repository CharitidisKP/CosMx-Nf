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
  sample_col = arg("sample_col", "sample_id"),
  svg_genes = arg("svg_genes", ""),
  max_genes = as.integer(arg("max_genes", "0")),
  order = arg("order", "maxmin"),
  cores = as.integer(arg("cores", "1"))
)
stopifnot(opt$order %in% c("maxmin", "AMMD", "Sum_coords"))

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(SpatialExperiment)
  library(nnSVG)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

## Exact maximum-minimum-distance order: start at the cell nearest the centre, then always take
## the cell farthest from every cell already ordered. With "maxmin" it replaces BRISC's grid
## approximation ("AMMD"), which drops cells from overfull grid boxes and then never finishes ##
order_maxmin_exact <- function(locs) {
  t0 <- Sys.time()
  x <- locs[, 1]
  y <- locs[, 2]
  ord <- integer(length(x))
  ord[1] <- which.min((x - mean(x))^2 + (y - mean(y))^2)
  d <- (x - x[ord[1]])^2 + (y - y[ord[1]])^2
  d[ord[1]] <- -1
  for (k in seq_along(ord)[-1]) {
    i <- which.max(d)
    ord[k] <- i
    d <- pmin(d, (x - x[i])^2 + (y - y[i])^2)
    d[i] <- -1
  }
  message(sprintf(
    "cci_nnsvg.R: maxmin order of %d cells in %.0f s",
    length(x),
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ))
  ord
}
if (opt$order == "maxmin") {
  utils::assignInNamespace("orderMaxMinLocal", order_maxmin_exact, ns = "BRISC")
}

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
norm <- getExpression(merged, values = "normalized", output = "matrix")
locs <- getSpatialLocations(merged, output = "data.table")
setkey(locs, cell_ID)

## Freed before nnSVG forks its workers, so no fork copies it ##
rm(merged)
invisible(gc())

## Genes: the SPARK-X hits of each sample, else the top max_genes by variance, else all ##
svg <- if (nzchar(opt$svg_genes)) {
  fread(opt$svg_genes)[significant == TRUE]
} else {
  NULL
}
genes <- rownames(norm)
if (is.null(svg) && opt$max_genes > 0 && opt$max_genes < length(genes)) {
  v <- Matrix::rowMeans(norm^2) - Matrix::rowMeans(norm)^2
  genes <- names(sort(v, decreasing = TRUE))[seq_len(opt$max_genes)]
}

res <- rbindlist(lapply(sort(unique(pd[[opt$sample_col]])), function(s) {
  ids <- pd$cell_ID[pd[[opt$sample_col]] == s]
  g <- if (is.null(svg)) genes else intersect(svg[sample == s]$gene, genes)
  x <- as.matrix(norm[g, ids, drop = FALSE])

  ## BRISC cannot fit a gene with no variance in the sample ##
  x <- x[matrixStats::rowVars(x) > 0, , drop = FALSE]
  if (nrow(x) < 2) {
    message("cci_nnsvg.R: ", s, " has < 2 testable genes, skipped")
    return(NULL)
  }
  spe <- SpatialExperiment(
    assays = list(logcounts = x),
    spatialCoords = as.matrix(locs[ids, .(sdimx, sdimy)])
  )
  message(sprintf(
    "cci_nnsvg.R: %s, %d genes x %d cells, %s order, started %s",
    s,
    nrow(x),
    ncol(x),
    opt$order,
    format(Sys.time(), "%H:%M")
  ))
  t0 <- Sys.time()
  spe <- tryCatch(
    nnSVG(
      spe,
      assay_name = "logcounts",
      order = if (opt$order == "Sum_coords") "Sum_coords" else "AMMD",
      n_threads = opt$cores
    ),
    error = function(e) {
      message("cci_nnsvg.R: nnSVG failed for ", s, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(spe)) {
    return(NULL)
  }
  message(sprintf(
    "cci_nnsvg.R: %s done in %.1f min",
    s,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
  rd <- as.data.table(as.data.frame(rowData(spe)), keep.rownames = "gene")
  rd[, sample := s]
}))

## Write ##
fwrite(res, file.path(opt$outdir, "nnsvg_spatially_variable_genes.csv"))
