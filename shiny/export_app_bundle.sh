#!/usr/bin/env bash
# Export a Giotto object to the Rshiny/spatial explorer contract.
# Writes <results>/local_export/, which load_cosmx_nf.R reads.
#
#   ./export_app_bundle.sh [run_id] [n_de_genes] [poly_compress] [obj_subdir]
# defaults: run_id=full  n_de_genes=25  poly_compress=xz  obj_subdir=objects/discover
# poly_compress=none skips the polygon write but still reports the vertex count.
# obj_subdir=objects/subcluster exports the --stage subcluster object instead, which carries the
# sub_<res> columns; objects/annotate also carries cell_type. The bundle lands in local_export either way.
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN=${1:-full}
NDE=${2:-25}
PCZ=${3:-xz}
OBJ=${4:-objects/discover}
RES="$PIPE/../Runs/$RUN/results"
[ -d "$RES/$OBJ/giotto" ] || { echo "no object at $RES/$OBJ/giotto" >&2; exit 1; }
# Same image, libraries and binds as the pipeline, read from conf/site.yaml
site() { sed -n "s/^$1:[[:space:]]*//p" "$PIPE/conf/site.yaml" 2>/dev/null | head -n 1 | tr -d "\"'"; }
SIF=$(site sif)
[ -n "$SIF" ] || { echo "no sif in $PIPE/conf/site.yaml (see conf/site.example.yaml)" >&2; exit 1; }
OPTS=(--cleanenv)
[ -n "$(site r_libs_user)" ]      && OPTS+=(--env "R_LIBS_USER=$(site r_libs_user)")
[ -n "$(site r_libs)" ]           && OPTS+=(--env "R_LIBS=$(site r_libs)")
[ -n "$(site container_python)" ] && OPTS+=(--env "RETICULATE_PYTHON=$(site container_python)")
IFS=, read -r -a BINDS <<<"$(site bind)"
for b in ${BINDS[@]+"${BINDS[@]}"}; do OPTS+=(-B "$b"); done
LOG="$PIPE/../Runs/$RUN/logs/export_app_bundle.txt"
mkdir -p "$(dirname "$LOG")" "$RES/local_export"

apptainer exec "${OPTS[@]}" \
  "$SIF" Rscript --vanilla - "$RES" "$NDE" "$LOG" "$PCZ" "$PIPE" "$OBJ" <<'RS'
suppressPackageStartupMessages({ library(Giotto); library(data.table); library(Matrix) })
a   <- commandArgs(trailingOnly = TRUE)
res <- a[1]; nde <- as.integer(a[2]); log <- a[3]; pcz <- a[4]; pipe <- a[5]
obj <- if (length(a) >= 6 && nzchar(a[6])) a[6] else "objects/discover"
ex  <- file.path(res, "local_export")   # where load_cosmx_nf() reads from
con <- file(log, "w"); sink(con, split = TRUE); sink(con, type = "message")

g   <- loadGiotto(file.path(res, obj, "giotto"))
cat("object:", file.path(res, obj, "giotto"), "\n")
pd  <- as.data.table(pDataDT(g))
loc <- getSpatialLocations(g, output = "data.table")
# cluster.R writes TWO umaps: "umap" (built on harmony when the object is merged) and
# "umap_pca" (raw PCA, only when harmony ran). Ask for "umap" BY NAME -- leaving name=NULL
# takes whichever Giotto treats as default, which is not something to leave to chance.
get_umap <- function(nm) tryCatch(
  getDimReduction(g, reduction = "cells", reduction_method = "umap",
                  name = nm, output = "matrix"),
  error = function(e) NULL)
umap <- get_umap("umap")          # post batch correction when harmony ran
upre <- get_umap("umap_pca")      # same cells, PCA only -- the uncorrected view
if (is.null(umap)) stop("no dim reduction named 'umap' on this object")
# which reductions exist tells us whether harmony ran: cluster.R only writes "umap_pca"
# when red != "pca". Wrapped -- a listing helper is not worth failing the export over.
have <- tryCatch(as.character(GiottoClass::list_dim_reductions(g)$name),
                 error = function(e) character())
cat("UMAP: exporting 'umap' and 'umap_pca'.  available:",
    if (length(have)) paste(have, collapse = ", ") else "(could not list)", "\n")
if (!is.null(upre)) {
  cat("  umap_pca present -> the merged run used HARMONY; 'umap' is the corrected",
      "embedding and 'umap_pca' the uncorrected one. Both are exported.\n")
} else {
  cat("  no umap_pca -> this run clustered on PCA, not harmony;",
      "there is no before/after pair to export\n")
}

# --- cells: coords + every non-numeric metadata column ------------------------------------
i <- match(pd$cell_ID, loc$cell_ID)
cells <- data.frame(cell_ID = pd$cell_ID, x = loc$sdimx[i], y = loc$sdimy[i])
j <- match(pd$cell_ID, rownames(umap))
cells$UMAP1 <- umap[j, 1]; cells$UMAP2 <- umap[j, 2]
# Before batch correction, as a SECOND embedding rather than a second bundle -- one cells
# table means the same colouring, gene and sample filters apply to both without reloading.
if (!is.null(upre)) {
  jp <- match(pd$cell_ID, rownames(upre))
  cells$UMAPpre1 <- upre[jp, 1]; cells$UMAPpre2 <- upre[jp, 2]
}

# Everything on the object except per-cell fit internals. A whitelist here silently cost us
# the IF channels; metadata is a few hundred KB against a 28 MB expression matrix, and the
# loader decides what is a category -- numerics stay numeric (IF intensity, Area).
drop <- grep("loglik|_delta|_entropy|_margin|^cell_ID$", names(pd), value = TRUE)
keep <- setdiff(names(pd), c(drop, names(cells)))
# Cluster ids arrive numeric but are categories; everything else keeps its natural type
catrx <- "^leiden_clus_res|^sub_res|^insitutype|^ht_|^cell_type|^monaco_|^fov$|^list_ID$"
numrx <- "_UMAP[0-9]$|_score$"
for (nm in keep) cells[[nm]] <- if (grepl(numrx, nm)) as.numeric(pd[[nm]])
  else if (!is.numeric(pd[[nm]]) || grepl(catrx, nm)) as.character(pd[[nm]])
  else as.numeric(pd[[nm]])
cat("dropped as fit internals:", if (length(drop)) paste(drop, collapse = ", ") else "(none)", "\n")
# Total counts per cell. The explorer needs it to build a DEPTH-MATCHED independence
# baseline for the gene-combination panel: detection tracks depth, so a naive baseline
# reports co-expression that is only sequencing depth.
rw_mat <- getExpression(g, values = "raw", output = "matrix")
cells$total_counts <- as.numeric(Matrix::colSums(rw_mat))[match(cells$cell_ID, colnames(rw_mat))]
cat("cells:", nrow(cells), "rows,", length(keep), "metadata columns\n")
cat("  ", paste(keep, collapse = ", "), "\n")
saveRDS(cells, file.path(ex, "cells_umap_clusters.rds"))

# --- expression: canonical markers + top-N DE per cluster of every diagnosed column --------
nm_mat <- getExpression(g, values = "normalized", output = "matrix")
mk <- tryCatch(read.csv(file.path(pipe, "assets/canonical_markers.csv"))$gene,
               error = function(e) character())
de_genes <- character()
for (f in list.files(file.path(res, "tables"), "_de_top\\.csv$", recursive = TRUE, full.names = TRUE)) {
  d <- fread(f)
  # name the j expression: an unnamed one comes back as V1 and $feats is silently NULL
  de_genes <- c(de_genes, d[, .(feats = head(feats[order(ranking)], nde)), by = cluster]$feats)
}
genes <- intersect(unique(c(mk, de_genes)), rownames(nm_mat))
if (!length(de_genes))
  cat("WARNING: no DE genes found under", file.path(res, "tables"), "\n")
cat("expression:", length(genes), "genes (", length(intersect(mk, rownames(nm_mat))),
    "markers +", length(setdiff(genes, mk)), "DE )\n")
sub <- nm_mat[genes, match(cells$cell_ID, colnames(nm_mat)), drop = FALSE]
saveRDS(Matrix(sub, sparse = TRUE), file.path(ex, "marker_expr_normalized.rds"))

# Raw counts as well. Detection questions ("is this cell CD19+?") must be asked on counts:
# normalization is library-size scaled, so ONE transcript in a low-depth cell outranks one
# in a high-depth cell -- the opposite of what you want when judging positivity.
graw   <- intersect(genes, rownames(rw_mat))   # rw_mat read above, with total_counts
if (length(graw) != length(genes))
  cat("WARNING: raw matrix is missing", length(genes) - length(graw), "of the exported genes\n")
subr <- rw_mat[graw, match(cells$cell_ID, colnames(rw_mat)), drop = FALSE]
saveRDS(Matrix(subr, sparse = TRUE), file.path(ex, "marker_expr_raw.rds"))
cat("raw counts:", nrow(subr), "genes exported alongside the normalized matrix\n")

# --- polygons: terra SpatVector -> flat vertex table, i indexing into `cells` --------------
# i is a ROW INDEX into the cells file written above, and `cells` is not touched after this
# point. Computing i before any reordering of cells is what keeps the two files aligned;
# a merge that re-sorts cells would silently point every polygon at the wrong cell.
ok <- tryCatch({
  gp  <- getPolygonInfo(g, polygon_name = "cell", return_giottoPolygon = TRUE)
  # terra::geom materialises an nverts x 5 double matrix before we slice it: ~40 bytes
  # per vertex transient, on top of the data.frame built from it.
  crd <- terra::geom(gp@spatVector)
  cat("polygons:", nrow(crd), "vertices,",
      round(nrow(crd) / length(unique(gp@spatVector$poly_ID)), 1), "per cell;",
      "terra matrix held", round(nrow(crd) * 40 / 1e6), "MB\n")
  # 2.0 B/vertex measured on real CosMx segmentation (2.68M verts -> 5.4 MB, xz). Pixel-grid
  # coordinates compress far better than the synthetic-ring figure this used to quote.
  cat("  estimated on disk:", round(nrow(crd) * 2.0 / 1e6, 1), "MB (xz),",
      round(nrow(crd) * 2.8 / 1e6, 1), "MB (gzip)\n")
  if (identical(pcz, "none")) {
    cat("  poly_compress=none -- vertex count reported, nothing written\n"); FALSE
  } else {
    ids <- gp@spatVector$poly_ID[crd[, "geom"]]
    # QC-dropped cells keep their polygons; without this they become NA indices
    k   <- !is.na(match(ids, cells$cell_ID))
    poly <- data.frame(i = match(ids[k], cells$cell_ID), px = crd[k, "x"], py = crd[k, "y"])
    # Segmentation vertices come off a raster mask, so they are usually integral. Storing them
    # as integer halves the in-memory table. Guarded: a sub-pixel double would be silently
    # truncated, which is a far worse trade than the RAM it saves.
    if (all(poly$px == round(poly$px)) && all(poly$py == round(poly$py)) &&
        max(abs(c(poly$px, poly$py))) < .Machine$integer.max) {
      storage.mode(poly$px) <- "integer"; storage.mode(poly$py) <- "integer"
      cat("  coordinates are integral -> stored as integer\n")
    } else cat("  coordinates are sub-pixel -> kept as double\n")
    cat("  kept", nrow(poly), "vertices for", length(unique(poly$i)), "cells;",
        sum(!k), "dropped (not in cells)\n")
    cat("  writing with compress =", pcz, "-- xz can take minutes at this size\n")
    saveRDS(poly, file.path(ex, "cell_polygons.rds"), compress = pcz)
    TRUE
  }
}, error = function(e) { cat("polygons FAILED:", conditionMessage(e), "\n"); FALSE })

cat("\n--- sizes ---\n")
for (f in list.files(ex, full.names = TRUE))
  cat(sprintf("%8.1f MB  %s\n", file.size(f) / 1e6, basename(f)))
sink(type = "message"); sink(); close(con)
RS
echo "log: $LOG"
