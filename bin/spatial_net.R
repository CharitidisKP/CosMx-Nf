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
  max_delaunay_dist = as.numeric(arg("max_delaunay_dist", "400")),
  prox_sim = as.integer(arg("prox_sim", "250")),
  sample_col = arg("sample_col", "sample_id"),
  seed = as.integer(arg("seed", "42")),
  um_per_px = as.numeric(arg("um_per_px", "0.12028")),
  immune_types = arg("immune_types"),
  untyped = arg("untyped", "Low_count,Unassigned"),
  label_colours = arg("label_colours")
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
## The figure helpers sit beside this script ##
bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(bin, "cosmx_figures.R"))
set_fig_format(arg("fig_format", "png"))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
if (!opt$celltype_column %in% names(pd)) {
  stop("spatial_net.R: no column '", opt$celltype_column, "' on the object")
}

## Positions from each cell's own outline, matched by cell ID. terra::centroids() drops the
## centre of a one-point outline but keeps its ID, and Giotto pairs centres with IDs by row,
## so every cell after such an outline had the next cell's centre (seen in one sample) ##
centres_by_id <- function(sv) {
  g <- as.data.table(terra::geom(sv))[hole == 0]
  g[,
    `:=`(
      x2 = data.table::shift(x, -1L, fill = x[1L]),
      y2 = data.table::shift(y, -1L, fill = y[1L])
    ),
    by = .(geom, part)
  ]
  g[, cr := x * y2 - x2 * y]
  ct <- g[,
    .(
      a = sum(cr) / 2,
      cx = sum((x + x2) * cr),
      cy = sum((y + y2) * cr),
      mx = mean(x),
      my = mean(y)
    ),
    by = .(geom, part)
  ]
  ## Largest part; an outline without area falls back to the mean of its vertices ##
  ct <- ct[order(geom, -abs(a))][, .SD[1L], by = geom]
  ct[,
    `:=`(
      sdimx = fifelse(abs(a) > 1e-9, cx / (6 * a), mx),
      sdimy = fifelse(abs(a) > 1e-9, cy / (6 * a), my)
    )
  ]
  out <- data.table(
    cell_ID = as.character(sv$poly_ID),
    sdimx = NA_real_,
    sdimy = NA_real_
  )
  out[ct$geom, `:=`(sdimx = ct$sdimx, sdimy = ct$sdimy)]
  out[]
}
old <- getSpatialLocations(merged, output = "spatLocsObj")
cur <- getSpatialLocations(merged, output = "data.table")
xy <- centres_by_id(
  getPolygonInfo(merged, polygon_name = "cell", return_giottoPolygon = TRUE)[]
)
xy <- xy[cur[, .(cell_ID)], on = "cell_ID"]
no_outline <- is.na(xy$sdimx)
xy[no_outline, `:=`(sdimx = cur$sdimx[no_outline], sdimy = cur$sdimy[no_outline])]
moved <- sqrt((xy$sdimx - cur$sdimx)^2 + (xy$sdimy - cur$sdimy)^2)
pos <- data.table(
  sample = pd[[opt$sample_col]][match(cur$cell_ID, pd$cell_ID)],
  moved,
  no_outline
)[,
  .(
    cells = .N,
    ## 20 px is about 2.4 um, well inside a cell ##
    moved_over_20px = sum(moved > 20),
    max_moved_px = round(max(moved)),
    no_outline = sum(no_outline)
  ),
  by = sample
]
merged <- setSpatialLocations(
  merged,
  createSpatLocsObj(
    coordinates = xy[, .(sdimx, sdimy, cell_ID)],
    name = objName(old),
    spat_unit = spatUnit(old),
    provenance = prov(old),
    verbose = FALSE
  ),
  verbose = FALSE
)
message(
  "spatial_net.R: positions from the cell outlines; moved > 20 px: ",
  paste0(pos$sample, " ", pos$moved_over_20px, collapse = ", ")
)
fwrite(pos, file.path(opt$outdir, "positions_check.csv"))

merged <- createSpatialNetwork(
  merged,
  method = "Delaunay",
  maximum_distance_delaunay = opt$max_delaunay_dist,
  name = "Delaunay_network",
  verbose = FALSE
)

## Samples are shift-joined, so an edge between two samples means the distance is too long ##
net <- getSpatialNetwork(
  merged,
  name = "Delaunay_network",
  output = "networkDT"
)
smp <- setNames(pd[[opt$sample_col]], pd$cell_ID)
cross <- sum(smp[net$from] != smp[net$to], na.rm = TRUE)
if (cross > 0) {
  warning(
    "spatial_net.R: ",
    cross,
    " cross-sample Delaunay edges, lower max_delaunay_dist"
  )
}

## Maps draw the display label, striped where its confidence is low, when the analysis label
## is cell_type; any other typing column is drawn as it is ##
style <- read_label_style(opt$label_colours)
use_display <- opt$celltype_column == "cell_type" &&
  all(c("cell_type_display", "label_confidence") %in% names(pd))
map_cells <- pd[, .(
  cell_ID,
  sample = get(opt$sample_col),
  group = as.character(get(if (use_display) "cell_type_display" else opt$celltype_column)),
  low = if (use_display) label_confidence == "low" else FALSE,
  fov = if ("fov" %in% names(pd)) as.character(fov) else "all"
)]
lv <- label_levels(map_cells$group, style)
outlines <- cell_outlines(merged)

## The null shuffles labels over the whole network. Pooled, differences in
## composition between samples would read as enrichment, so it runs per sample ##
samples <- sort(unique(pd[[opt$sample_col]]))
enr <- list()
for (s in samples) {
  gs <- subsetGiotto(
    merged,
    cell_ids = pd$cell_ID[pd[[opt$sample_col]] == s]
  )
  cpe <- cellProximityEnrichment(
    gs,
    cluster_column = opt$celltype_column,
    spatial_network_name = "Delaunay_network",
    number_of_simulations = opt$prox_sim,
    adjust_method = "fdr",
    seed_number = opt$seed
  )
  enr[[s]] <- cpe$enrichm_res[, sample := s]

  ## Every cell as its outline, as the presentation's maps; immune cells on top ##
  d <- map_cells[sample == s]
  n <- d[, .N, by = group]
  slv <- intersect(lv, n$group)
  d[, group := factor(group, levels = slv)]
  m <- biopsy_map(
    d,
    outlines[cell_ID %in% d$cell_ID],
    label_colours(slv, style),
    stats::setNames(
      sprintf(
        "%s (%s)",
        display_name(slv, style),
        prettyNum(n$N[match(slv, n$group)], big.mark = ",")
      ),
      slv
    ),
    opt$um_per_px,
    "Cell type",
    title = s,
    on_top = slv[label_layer(slv, style, split_arg(opt$untyped), split_arg(opt$immune_types)) == 3L]
  )
  save_fig(
    m$plot,
    file.path(opt$outdir, paste0("spatial_celltype_", s, ".png")),
    m$w,
    m$h,
    map = TRUE
  )
}
enr <- rbindlist(enr)
enr[, unified_int := as.character(unified_int)]
setcolorder(enr, "sample")

## Heatmap: log2 enrichment per cell-type pair and sample, both orders of each pair ##
pair <- tstrsplit(enr$unified_int, "--", fixed = TRUE)
hm <- data.table(
  sample = enr$sample,
  a = pair[[1]],
  b = pair[[2]],
  enrichm = enr$enrichm,
  sig = (enr$p.adj_higher < 0.05 | enr$p.adj_lower < 0.05) %in% TRUE
)
hm <- unique(rbind(hm, hm[, .(sample, a = b, b = a, enrichm, sig)]))
n_col <- min(3L, length(samples))
q <- ggplot(hm, aes(a, b, fill = enrichm)) +
  geom_tile() +
  geom_point(data = hm[sig == TRUE], size = 0.6) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
  scale_x_discrete(labels = function(x) display_name(x, style)) +
  scale_y_discrete(labels = function(x) display_name(x, style)) +
  facet_wrap(~sample, ncol = n_col) +
  labs(
    x = NULL,
    y = NULL,
    fill = "log2 enrichment",
    caption = "Dots: FDR < 0.05"
  ) +
  theme_cosmx(9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Write ##
fwrite(enr, file.path(opt$outdir, "proximity_enrichment.csv"))
save_fig(
  q,
  file.path(opt$outdir, "proximity_heatmap.png"),
  1 + 4.5 * n_col,
  1 + 4.5 * ceiling(length(samples) / n_col)
)
write.csv(
  data.frame(
    n_edges = nrow(net),
    cross_sample_edges = cross,
    max_delaunay_dist = opt$max_delaunay_dist,
    prox_sim = opt$prox_sim
  ),
  file.path(opt$outdir, "spatial_net_summary.csv"),
  row.names = FALSE
)
saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
