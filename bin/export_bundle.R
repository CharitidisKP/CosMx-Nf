#!/usr/bin/env Rscript
## Writes the app bundle that shiny/load_cosmx_nf.R reads: cells, expression and outlines ##
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input"),
  outdir = arg("outdir", "."),
  log = arg("log"),
  genes = arg("genes", "all"),
  poly_compress = arg("poly_compress", "xz"),
  markers = arg("markers"),
  tables = arg("tables"),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON"))
)
if (opt$genes != "all" && !grepl("^[0-9]+$", opt$genes)) {
  stop("export_bundle.R: genes must be 'all' or a number, not '", opt$genes, "'")
}
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(Matrix)
})
if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
  options(giotto.py_path = opt$python)
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
if (nzchar(opt$log)) {
  con <- file(opt$log, "w")
  sink(con, split = TRUE)
  sink(con, type = "message")
}

g <- loadGiotto(opt$input)
cat("object:", opt$input, "\n")
pd <- as.data.table(pDataDT(g))
loc <- getSpatialLocations(g, output = "data.table")
## A UMAP by name, NULL when the object has none of that name ##
get_umap <- function(nm) {
  tryCatch(
    getDimReduction(g, reduction = "cells", reduction_method = "umap", name = nm, output = "matrix"),
    error = function(e) NULL
  )
}
umap <- get_umap("umap")
upre <- get_umap("umap_pca")
if (is.null(umap)) stop("export_bundle.R: no dim reduction named 'umap' on this object")
have <- tryCatch(as.character(GiottoClass::list_dim_reductions(g)$name), error = function(e) character())
cat("UMAP: exporting 'umap' and 'umap_pca'. Available:", if (length(have)) paste(have, collapse = ", ") else "(not listed)", "\n")
if (!is.null(upre)) {
  cat("  umap_pca is there, so 'umap' is after batch correction and 'umap_pca' before. Both are exported.\n")
} else {
  cat("  no umap_pca, so this run clustered on PCA and there is no before and after pair\n")
}

## Cells: positions, both UMAPs and every metadata column but the fit internals ##
i <- match(pd$cell_ID, loc$cell_ID)
cells <- data.frame(cell_ID = pd$cell_ID, x = loc$sdimx[i], y = loc$sdimy[i])
j <- match(pd$cell_ID, rownames(umap))
cells$UMAP1 <- umap[j, 1]
cells$UMAP2 <- umap[j, 2]
if (!is.null(upre)) {
  jp <- match(pd$cell_ID, rownames(upre))
  cells$UMAPpre1 <- upre[jp, 1]
  cells$UMAPpre2 <- upre[jp, 2]
}
drop <- grep("loglik|_delta|_entropy|_margin|^cell_ID$", names(pd), value = TRUE)
keep <- setdiff(names(pd), c(drop, names(cells)))
## Column names written as text, and those written as numbers ##
catrx <- "^leiden_clus_res|^sub_res|^insitutype|^ht_|^cell_type|^monaco_|^fov$|^list_ID$"
numrx <- "_UMAP[0-9]$|_score$"
for (nm in keep) {
  cells[[nm]] <- if (grepl(numrx, nm)) {
    as.numeric(pd[[nm]])
  } else if (!is.numeric(pd[[nm]]) || grepl(catrx, nm)) {
    as.character(pd[[nm]])
  } else {
    as.numeric(pd[[nm]])
  }
}
cat("dropped as fit internals:", if (length(drop)) paste(drop, collapse = ", ") else "(none)", "\n")
rw_mat <- getExpression(g, values = "raw", output = "matrix")
cells$total_counts <- as.numeric(Matrix::colSums(rw_mat))[match(cells$cell_ID, colnames(rw_mat))]
cat("cells:", nrow(cells), "rows,", length(keep), "metadata columns\n")
cat("  ", paste(keep, collapse = ", "), "\n")
saveRDS(cells, file.path(opt$outdir, "cells_umap_clusters.rds"))

## Expression: every gene, or the markers and the top N genes per cluster of each diagnosed column ##
nm_mat <- getExpression(g, values = "normalized", output = "matrix")
if (opt$genes == "all") {
  genes <- rownames(nm_mat)
  cat("expression: all", length(genes), "genes\n")
} else {
  nde <- as.integer(opt$genes)
  mk <- tryCatch(read.csv(opt$markers)$gene, error = function(e) character())
  de_genes <- character()
  for (f in list.files(opt$tables, "_de_top\\.csv$", recursive = TRUE, full.names = TRUE)) {
    d <- fread(f)
    de_genes <- c(de_genes, d[, .(feats = head(feats[order(ranking)], nde)), by = cluster]$feats)
  }
  genes <- intersect(unique(c(mk, de_genes)), rownames(nm_mat))
  if (!length(de_genes)) cat("WARNING: no DE genes under", opt$tables, "\n")
  cat("expression:", length(genes), "genes,", length(intersect(mk, rownames(nm_mat))), "markers and",
      length(setdiff(genes, mk)), "DE genes\n")
}
sub <- nm_mat[genes, match(cells$cell_ID, colnames(nm_mat)), drop = FALSE]
saveRDS(Matrix(sub, sparse = TRUE), file.path(opt$outdir, "marker_expr_normalized.rds"))
graw <- intersect(genes, rownames(rw_mat))
if (length(graw) != length(genes)) {
  cat("WARNING: the raw matrix lacks", length(genes) - length(graw), "of the exported genes\n")
}
subr <- rw_mat[graw, match(cells$cell_ID, colnames(rw_mat)), drop = FALSE]
saveRDS(Matrix(subr, sparse = TRUE), file.path(opt$outdir, "marker_expr_raw.rds"))
cat("raw counts:", nrow(subr), "genes beside the normalized matrix\n")

## Outlines: one row per vertex, i the row of its cell in cells ##
ok <- tryCatch(
  {
    gp <- getPolygonInfo(g, polygon_name = "cell", return_giottoPolygon = TRUE)
    crd <- terra::geom(gp@spatVector)
    cat("polygons:", nrow(crd), "vertices,", round(nrow(crd) / length(unique(gp@spatVector$poly_ID)), 1), "per cell\n")
    cat("  about", round(nrow(crd) * 2.0 / 1e6, 1), "MB on disk with xz,", round(nrow(crd) * 2.8 / 1e6, 1), "MB with gzip\n")
    if (identical(opt$poly_compress, "none")) {
      cat("  poly_compress none: vertices counted, nothing written\n")
      FALSE
    } else {
      ids <- gp@spatVector$poly_ID[crd[, "geom"]]
      k <- !is.na(match(ids, cells$cell_ID))
      poly <- data.frame(i = match(ids[k], cells$cell_ID), px = crd[k, "x"], py = crd[k, "y"])
      whole <- all(poly$px == round(poly$px)) && all(poly$py == round(poly$py)) &&
        max(abs(c(poly$px, poly$py))) < .Machine$integer.max
      if (whole) {
        storage.mode(poly$px) <- "integer"
        storage.mode(poly$py) <- "integer"
        cat("  whole-pixel coordinates, stored as integers\n")
      } else {
        cat("  sub-pixel coordinates, kept as doubles\n")
      }
      cat("  kept", nrow(poly), "vertices for", length(unique(poly$i)), "cells,", sum(!k), "dropped with their cells\n")
      cat("  writing with compress", opt$poly_compress, "\n")
      saveRDS(poly, file.path(opt$outdir, "cell_polygons.rds"), compress = opt$poly_compress)
      TRUE
    }
  },
  error = function(e) {
    cat("polygons FAILED:", conditionMessage(e), "\n")
    FALSE
  }
)

cat("\nsizes\n")
for (f in list.files(opt$outdir, "\\.rds$", full.names = TRUE)) {
  cat(sprintf("%8.1f MB  %s\n", file.size(f) / 1e6, basename(f)))
}
if (nzchar(opt$log)) {
  sink(type = "message")
  sink()
  close(con)
}
