#!/usr/bin/env Rscript
## The curated labels onto the object: each cell's cluster label in --cluster_basis, then the
## per-cell edits in --overrides, by the rules in label_state.R. Writes cell_type (every analysis
## reads it), cell_type_l1, cell_type_source, cell_type_display (figures draw it) and
## label_confidence ##
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input"),
  labels = arg("labels"),
  overrides = arg("overrides"),
  cluster_basis = arg("cluster_basis"),
  label_state = arg("label_state"),
  labels_from = arg("labels_from"),
  overrides_from = arg("overrides_from"),
  immune_types = arg("immune_types"),
  untyped = arg("untyped", "Low_count,Unassigned"),
  markers = arg("markers"),
  label_colours = arg("label_colours"),
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
## The helpers sit beside this script ##
bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(bin, "label_state.R"))
source(file.path(bin, "cosmx_figures.R"))
set_fig_format(arg("fig_format", "png"))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
if (!opt$cluster_basis %in% names(pd)) {
  stop("apply_labels.R: no column '", opt$cluster_basis, "' on the object")
}

map <- read_cluster_labels(opt$labels)
clusters <- as.character(unique(pd[[opt$cluster_basis]]))
bad <- setdiff(clusters, names(map)[!is.na(map) & nzchar(map)])
if (length(bad)) {
  stop(
    "apply_labels.R: unmapped or blank clusters in ",
    opt$cluster_basis,
    ": ",
    paste(bad, collapse = ", ")
  )
}

ov <- read_label_overrides(opt$overrides)
miss <- setdiff(ov$cell_ID, pd$cell_ID)
if (length(miss)) {
  message(
    "apply_labels.R: ",
    length(unique(miss)),
    " override cells are not in the object and were ignored"
  )
}
st <- as.data.table(apply_label_state(pd$cell_ID, pd[[opt$cluster_basis]], map, ov))

## Level 1 is read off immune_types, so it follows whatever the labels are called ##
immune <- split_arg(opt$immune_types)
absent <- setdiff(immune, st$cell_type)
if (length(absent)) {
  message(
    "apply_labels.R: immune_types not among the labels: ",
    paste(absent, collapse = ", ")
  )
}
st[, cell_type_l1 := fifelse(cell_type %in% immune, "immune", "non_immune")]
merged <- addCellMetadata(
  merged,
  new_metadata = st[, .(
    cell_ID,
    cell_type,
    cell_type_l1,
    cell_type_source,
    cell_type_display,
    label_confidence
  )],
  by_column = TRUE,
  column_cell_ID = "cell_ID"
)

## What the edits changed: high edits from the cluster label to cell_type, low edits from
## cell_type to the display label ##
st[, cluster_label := unname(map[as.character(pd[[opt$cluster_basis]])])]
fwrite(
  rbind(
    st[cell_type_source == "manual", .(confidence = "high", from = cluster_label, to = cell_type)],
    st[label_confidence == "low", .(confidence = "low", from = cell_type, to = cell_type_display)]
  )[, .(n_cells = .N), by = .(confidence, from, to)],
  file.path(opt$outdir, "label_overrides_applied.csv")
)
n_hi <- sum(st$cell_type_source == "manual")
n_lo <- sum(st$label_confidence == "low")
if (n_hi + n_lo) {
  message(
    "apply_labels.R: ",
    n_hi,
    " cells relabelled from overrides, ",
    n_lo,
    " given a low-confidence display label"
  )
}

## Write ##
n <- table(st$cell_type)
write.csv(
  data.frame(
    cell_type = names(n),
    cell_type_l1 = ifelse(names(n) %in% immune, "immune", "non_immune"),
    n_cells = as.integer(n)
  ),
  file.path(opt$outdir, "cell_type_counts.csv"),
  row.names = FALSE
)
fwrite(
  st[, .(n_cells = .N), by = .(cell_type_display, label_confidence)][
    order(cell_type_display, label_confidence)
  ],
  file.path(opt$outdir, "cell_type_display_counts.csv")
)
## The label state this object was built from, which the later stages compare against, and
## the files it came from ##
if (nzchar(opt$label_state)) {
  writeLines(opt$label_state, file.path(opt$outdir, "label_state.txt"))
}
fwrite(
  data.table(
    cluster_basis = opt$cluster_basis,
    labels = opt$labels_from,
    overrides = opt$overrides_from,
    override_rows = nrow(ov),
    relabelled = n_hi,
    low_confidence = n_lo,
    label_state = opt$label_state
  ),
  file.path(opt$outdir, "label_source.csv")
)

style <- read_label_style(opt$label_colours)
untyped <- split_arg(opt$untyped)

## Cell types on the UMAP: the display label, low-confidence labels
## slashed, immune cells drawn last ##
um <- tryCatch(
  getDimReduction(
    merged,
    reduction = "cells",
    reduction_method = "umap",
    name = "umap",
    output = "matrix"
  ),
  error = function(e) NULL
)
if (is.null(um)) {
  message("apply_labels.R: no UMAP named 'umap' on the object, so no UMAP figure")
} else {
  j <- match(st$cell_ID, rownames(um))
  u <- st[!is.na(j), .(
    UMAP1 = um[j[!is.na(j)], 1],
    UMAP2 = um[j[!is.na(j)], 2],
    group = cell_type_display,
    low = label_confidence == "low"
  )]
  lv <- label_levels(u$group, style)
  u[, `:=`(
    group = factor(group, levels = lv),
    layer = label_layer(group, style, untyped, immune)
  )]
  p <- umap_plot(u, label_colours(lv, style), stats::setNames(display_name(lv, style), lv))
  save_fig(p, file.path(opt$outdir, "umap_cell_type.png"), 10, 7.5)
}

## Dotplot: the three panel markers most specific to each label, i.e. highest detection in
## the label minus the highest detection in any other label ##
if (nzchar(opt$markers) && file.exists(opt$markers)) {
  expr <- getExpression(merged, values = "normalized", output = "matrix")
  genes <- intersect(unique(read.csv(opt$markers)$gene), rownames(expr))
  expr <- expr[genes, st$cell_ID, drop = FALSE]
  lab_f <- factor(st$cell_type, levels = label_levels(st$cell_type, style))
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
    scale_y_discrete(labels = function(x) display_name(x, style)) +
    labs(x = NULL, y = NULL, size = "Detected", colour = "Scaled mean") +
    theme_cosmx() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  save_fig(
    p,
    file.path(opt$outdir, "marker_dotplot.png"),
    max(8, 2 + 0.25 * length(top)),
    max(4, 1.5 + 0.35 * nlevels(lab_f))
  )
}
saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
