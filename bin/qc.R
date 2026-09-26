#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  sample_id = arg("sample_id"),
  input = arg("input"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON")),
  gene_min_cells = as.integer(arg("gene_min_cells", "5")),
  cell_min_genes = as.integer(arg("cell_min_genes", "5")),
  count_cap = as.integer(arg("count_cap", "200")),
  count_quantile = as.numeric(arg("count_quantile", "0.10")),
  area_max = as.numeric(arg("area_max", "30000")),
  fov_min_count = as.numeric(arg("fov_min_count", "100")),
  split_ratio_threshold = as.numeric(arg("split_ratio_threshold", "0.5"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(Matrix)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

g <- loadGiotto(opt$input)

## Mean negprobe counts per cell. estimateBackground() expects exactly this ##
neg <- getExpression(
  g,
  feat_type = "negprobes",
  values = "raw",
  output = "matrix"
)
neg_mean <- Matrix::colMeans(neg[
  grep("^Negative", rownames(neg)),
  ,
  drop = FALSE
])
g <- addCellMetadata(
  g,
  by_column = TRUE,
  column_cell_ID = "cell_ID",
  new_metadata = data.frame(
    cell_ID = names(neg_mean),
    neg_mean = as.numeric(neg_mean)
  )
)

## Set counts as feat_type. Negprobes need to be excluded from total expr ##
g <- addStatistics(g, feat_type = "rna", expression_values = "raw")
m <- as.data.frame(pDataDT(g))
loc <- getSpatialLocations(g, output = "data.table")
setkey(loc, cell_ID)

## Count floor is self limiting: can never drop more than count_quantile of a sample ##
min_count <- min(opt$count_cap, quantile(m$total_expr, opt$count_quantile))

flag_min_count <- m$total_expr < min_count
flag_min_genes <- m$nr_feats < opt$cell_min_genes
flag_area <- if (is.null(m$Area)) rep(FALSE, nrow(m)) else m$Area > opt$area_max

## Bruker's FOV call when the export carries it, otherwise mean counts per cell ##
flag_fov <- if (!is.null(m$qcFlagsFOV)) {
  m$qcFlagsFOV != "Pass"
} else {
  ave(m$total_expr, m$fov, FUN = mean) < opt$fov_min_count
}

## Border cells are flagged, never dropped. 0 means the cell is not on a border ##
flag_border <- if (is.null(m$SplitRatioToLocal)) {
  rep(FALSE, nrow(m))
} else {
  m$SplitRatioToLocal > 0 & m$SplitRatioToLocal < opt$split_ratio_threshold
}

keep <- !flag_min_count & !flag_min_genes & !flag_area & !flag_fov

## Flags are written for every cell, with coordinates, so exclusions can be plotted in space ##
write.csv(
  data.frame(
    cell_ID = m$cell_ID,
    flag_min_count,
    flag_min_genes,
    flag_area,
    flag_fov,
    flag_border,
    keep,
    total_expr = m$total_expr,
    nr_feats = m$nr_feats,
    neg_mean = m$neg_mean,
    Area = if (is.null(m$Area)) NA_real_ else m$Area,
    fov = m$fov,
    sdimx = loc[m$cell_ID]$sdimx,
    sdimy = loc[m$cell_ID]$sdimy,
    min_count = min_count
  ),
  file.path(opt$outdir, paste0(opt$sample_id, "_qc_flags.csv")),
  row.names = FALSE
)

## Filter: cells first, so the gene threshold counts detections only in kept cells ##
g <- subsetGiotto(g, cell_ids = m$cell_ID[keep])
g <- filterGiotto(
  g,
  feat_type = "rna",
  expression_threshold = 1,
  feat_det_in_min_cells = opt$gene_min_cells,
  min_det_feats_per_cell = 0
)

## The two thresholds on their distributions. A cell can carry several flags ##
dist <- rbind(
  data.table(
    metric = "Total counts (log10)",
    value = log10(m$total_expr + 1),
    kept = keep
  ),
  if (!is.null(m$Area)) data.table(metric = "Area", value = m$Area, kept = keep)
)
cut <- rbind(
  data.table(metric = "Total counts (log10)", value = log10(min_count + 1)),
  if (!is.null(m$Area)) data.table(metric = "Area", value = opt$area_max)
)
p <- ggplot(dist, aes(value, fill = ifelse(kept, "kept", "removed"))) +
  geom_histogram(bins = 60) +
  geom_vline(data = cut, aes(xintercept = value), linetype = "dashed") +
  facet_wrap(~metric, scales = "free") +
  labs(
    title = opt$sample_id,
    subtitle = sprintf(
      "%d cells, %d kept. Flagged: counts %d, genes %d, area %d, FOV %d, border %d",
      nrow(m),
      sum(keep),
      sum(flag_min_count),
      sum(flag_min_genes),
      sum(flag_area),
      sum(flag_fov),
      sum(flag_border)
    ),
    x = NULL,
    y = "Cells",
    fill = NULL
  ) +
  theme_minimal(base_size = 11)

## Write ##
ggsave(
  file.path(opt$outdir, paste0("qc_", opt$sample_id, ".png")),
  p,
  width = 10,
  height = 4,
  dpi = 150,
  bg = "white"
)
write.csv(
  data.frame(
    sample_id = opt$sample_id,
    cells_before = nrow(m),
    cells_after = sum(keep),
    genes_after = nrow(fDataDT(g, feat_type = "rna")),
    min_count = min_count,
    flagged_min_count = sum(flag_min_count),
    flagged_min_genes = sum(flag_min_genes),
    flagged_area = sum(flag_area),
    flagged_fov = sum(flag_fov),
    flagged_border = sum(flag_border),
    median_counts_kept = stats::median(m$total_expr[keep]),
    median_genes_kept = stats::median(m$nr_feats[keep])
  ),
  file.path(opt$outdir, paste0(opt$sample_id, "_qc_summary.csv")),
  row.names = FALSE
)
saveGiotto(g, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
