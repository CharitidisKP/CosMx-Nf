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
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
if (!opt$celltype_column %in% names(pd)) {
  stop("spatial_net.R: no column '", opt$celltype_column, "' on the object")
}

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

  p <- spatInSituPlotPoints(
    gs,
    show_polygon = TRUE,
    polygon_feat_type = "cell",
    polygon_fill = opt$celltype_column,
    polygon_fill_as_factor = TRUE,
    polygon_line_size = 0.05,
    feats = NULL,
    show_plot = FALSE,
    return_plot = TRUE
  )
  ggsave(
    file.path(opt$outdir, paste0("spatial_celltype_", s, ".png")),
    p + labs(title = s),
    width = 7,
    height = 7,
    dpi = 150,
    bg = "white"
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
  facet_wrap(~sample, ncol = n_col) +
  labs(
    x = NULL,
    y = NULL,
    fill = "log2 enrichment",
    caption = "Dots: FDR < 0.05"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Write ##
fwrite(enr, file.path(opt$outdir, "proximity_enrichment.csv"))
ggsave(
  file.path(opt$outdir, "proximity_heatmap.png"),
  q,
  width = 1 + 4.5 * n_col,
  height = 1 + 4.5 * ceiling(length(samples) / n_col),
  dpi = 150,
  bg = "white"
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
