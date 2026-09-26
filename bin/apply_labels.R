#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input"),
  labels = arg("labels"),
  cluster_basis = arg("cluster_basis"),
  immune_types = arg("immune_types"),
  markers = arg("markers"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON"))
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

merged <- loadGiotto(opt$input)

## colClasses: read.csv types "3.1" and "3.10" as numeric and both collapse to 3.1 ##
lab <- read.csv(opt$labels, colClasses = "character")
if (!all(c("cluster_id", "label") %in% names(lab))) {
  stop("apply_labels.R: ", opt$labels, " needs columns cluster_id and label")
}
dup <- unique(lab$cluster_id[duplicated(lab$cluster_id)])
if (length(dup)) {
  stop(
    "apply_labels.R: duplicate cluster_id in ",
    opt$labels,
    ": ",
    paste(dup, collapse = ", ")
  )
}
map <- setNames(lab$label, lab$cluster_id)

clusters <- as.character(unique(pDataDT(merged)[[opt$cluster_basis]]))
bad <- setdiff(clusters, names(map)[!is.na(map) & nzchar(map)])
if (length(bad)) {
  stop(
    "apply_labels.R: unmapped or blank clusters in ",
    opt$cluster_basis,
    ": ",
    paste(bad, collapse = ", ")
  )
}

merged <- annotateGiotto(
  merged,
  annotation_vector = map,
  cluster_column = opt$cluster_basis,
  name = "cell_type"
)

## Level 1 is read off immune_types, so it follows whatever the labels are called ##
immune <- trimws(strsplit(opt$immune_types, ",")[[1]])
immune <- immune[nzchar(immune)]
absent <- setdiff(immune, map)
if (length(absent)) {
  message(
    "apply_labels.R: immune_types not among the labels: ",
    paste(absent, collapse = ", ")
  )
}
l1 <- setNames(ifelse(map %in% immune, "immune", "non_immune"), names(map))
merged <- annotateGiotto(
  merged,
  annotation_vector = l1,
  cluster_column = opt$cluster_basis,
  name = "cell_type_l1"
)

## Write ##
p <- plotUMAP(
  merged,
  dim_reduction_name = "umap",
  cell_color = "cell_type",
  point_size = 0.6,
  show_plot = FALSE,
  return_plot = TRUE
)
ggsave(
  file.path(opt$outdir, "umap_cell_type.png"),
  p,
  width = 8,
  height = 6,
  dpi = 150,
  bg = "white"
)

## Dotplot: the three panel markers most specific to each label, i.e. highest detection in
## the label minus the highest detection in any other label ##
if (nzchar(opt$markers) && file.exists(opt$markers)) {
  pd <- pDataDT(merged)
  expr <- getExpression(merged, values = "normalized", output = "matrix")
  genes <- intersect(unique(read.csv(opt$markers)$gene), rownames(expr))
  expr <- expr[genes, pd$cell_ID, drop = FALSE]
  lab_f <- factor(pd$cell_type)
  ind <- sparseMatrix(
    i = seq_along(lab_f),
    j = as.integer(lab_f),
    x = 1,
    dims = c(length(lab_f), nlevels(lab_f))
  )
  n_lab <- Matrix::colSums(ind)
  det <- as.matrix((expr > 0) %*% ind) / rep(n_lab, each = length(genes))
  avg <- as.matrix(expr %*% ind) / rep(n_lab, each = length(genes))
  dimnames(det) <- dimnames(avg) <- list(genes, levels(lab_f))
  top <- unique(unlist(lapply(levels(lab_f), function(l) {
    other <- apply(det[, setdiff(levels(lab_f), l), drop = FALSE], 1, max)
    spec <- (det[, l] - other)[det[, l] >= 0.1]
    head(names(sort(spec, decreasing = TRUE)), 3)
  })))
  z <- t(scale(t(avg[top, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  dd <- data.table(
    gene = factor(rep(top, times = nlevels(lab_f)), levels = top),
    label = factor(
      rep(levels(lab_f), each = length(top)),
      levels = rev(levels(lab_f))
    ),
    detection = as.vector(det[top, ]),
    z = as.vector(z)
  )
  p <- ggplot(dd, aes(gene, label, size = detection, colour = z)) +
    geom_point() +
    scale_colour_gradient2(
      low = "steelblue",
      mid = "grey90",
      high = "firebrick"
    ) +
    scale_size(range = c(0, 5), limits = c(0, 1)) +
    labs(x = NULL, y = NULL, size = "Detected", colour = "Scaled mean") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(
    file.path(opt$outdir, "marker_dotplot.png"),
    p,
    width = max(8, 2 + 0.25 * length(top)),
    height = max(4, 1.5 + 0.35 * nlevels(lab_f)),
    dpi = 150,
    bg = "white"
  )
}
n <- table(pDataDT(merged)$cell_type)
write.csv(
  data.frame(
    cell_type = names(n),
    cell_type_l1 = ifelse(names(n) %in% immune, "immune", "non_immune"),
    n_cells = as.integer(n)
  ),
  file.path(opt$outdir, "cell_type_counts.csv"),
  row.names = FALSE
)
saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
