# Export a Giotto object to the two files app.R reads. Run this beside your real
# object (server-side); the app itself never loads Giotto.
#
# CONTRACT — app.R needs exactly this and nothing else:
#   demo_cells.rds  one row per cell: cell_id, x, y, <embedding cols>, and any
#                   number of FACTOR columns. Every factor becomes a "Colour by"
#                   option automatically, so resolutions, annotations, sample and
#                   group all arrive just by being factors.
#   demo_poly.rds   vertex table: i (row index into cells), px, py.
#                   Index, not a join key — a 1.5M-row merge per render is not free.

library(Giotto)

export_giotto <- function(gobj,
                          out        = ".",
                          spat_unit  = "cell",
                          feat_type  = "rna",
                          embeddings = c(UMAP = "umap", PCA = "pca"),
                          meta_cols  = NULL,   # NULL = every non-numeric metadata column
                          poly_name  = "cell") {

  md  <- data.table::as.data.table(pDataDT(gobj, spat_unit = spat_unit, feat_type = feat_type))
  loc <- getSpatialLocations(gobj, spat_unit = spat_unit, output = "data.table")

  cells <- data.frame(
    cell_id = md$cell_ID,
    x = loc$sdimx[match(md$cell_ID, loc$cell_ID)],
    y = loc$sdimy[match(md$cell_ID, loc$cell_ID)]
  )

  for (nm in names(embeddings)) {
    dr <- getDimReduction(gobj, spat_unit = spat_unit, feat_type = feat_type,
                          reduction = "cells", reduction_method = embeddings[[nm]],
                          output = "matrix")
    dr <- dr[match(md$cell_ID, rownames(dr)), 1:2, drop = FALSE]
    cells[[paste0(nm, "_1")]] <- dr[, 1]
    cells[[paste0(nm, "_2")]] <- dr[, 2]
  }

  if (is.null(meta_cols)) {
    meta_cols <- setdiff(names(md)[!vapply(md, is.numeric, logical(1))], "cell_ID")
  }
  for (nm in meta_cols) cells[[nm]] <- factor(md[[nm]])

  # Polygons: terra SpatVector -> flat vertex table, i indexing into cells
  gp  <- getPolygonInfo(gobj, polygon_name = poly_name, return_giottoPolygon = TRUE)
  crd <- terra::geom(gp@spatVector)                       # geom, part, x, y, hole
  ids <- gp@spatVector$poly_ID[crd[, "geom"]]
  keep <- !is.na(match(ids, cells$cell_id))
  poly <- data.frame(i  = match(ids[keep], cells$cell_id),
                     px = crd[keep, "x"], py = crd[keep, "y"])

  saveRDS(cells, file.path(out, "demo_cells.rds"))
  saveRDS(poly,  file.path(out, "demo_poly.rds"), compress = "xz")
  message(nrow(cells), " cells, ", nrow(poly), " vertices, ",
          length(meta_cols), " colour-by columns")
  invisible(list(cells = cells, poly = poly))
}

# export_giotto(my_gobj, out = "path/to/CosMx-Nf/results/local_export")

# VERIFIED against a real Giotto object (CosMx 6K, 8 samples). sdimx/sdimy, poly_ID and
# the "umap" reduction name all held; terra::geom + the i-index contract produced
# 2,675,746 vertices for 126,098 cells (21.2 per cell) with 0 unmatched.
# Sizes measured: cell_polygons.rds 5.4 MB at compress = "xz" -- 2.0 bytes/vertex, because
# segmentation vertices land on a pixel grid and xz crushes low-entropy mantissas. gzip is
# ~40% larger. terra::geom holds an nverts x 5 double matrix (~107 MB here) before slicing.
#
# For the pipeline this is superseded by export_app_bundle.sh, which runs server-side in the
# container and also writes x/y, the expression matrix and every annotation column. Keep this
# file as the minimal statement of the contract app.R depends on.
# If a reduction is missing: showGiottoDimRed(gobj) lists what exists.
