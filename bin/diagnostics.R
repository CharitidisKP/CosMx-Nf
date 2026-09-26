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
  cluster_column = arg("cluster_column"),
  markers = arg("markers"),
  top_n = as.integer(arg("top_n", "25")),
  min_detection = as.numeric(arg("min_detection", "0.05")),
  min_panel_genes = as.integer(arg("min_panel_genes", "3")),
  min_score = as.numeric(arg("min_score", "0.5"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(ggplot2)
  library(Matrix)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
col <- opt$cluster_column

merged <- loadGiotto(opt$input)
mk <- read.csv(opt$markers)
markers <- split(mk$gene, mk$cell_type)
norm_mat <- getExpression(merged, values = "normalized", output = "matrix")
pd <- pDataDT(merged)
if (!col %in% names(pd)) {
  clus <- grep(
    "^leiden_clus_res|^sub_|insitutype|niche_id|cell_type",
    names(pd),
    value = TRUE
  )
  stop(
    "diagnostics.R: no column '",
    col,
    "'. Cluster columns present: ",
    if (length(clus)) paste(clus, collapse = ", ") else "none"
  )
}
cl <- as.character(setNames(pd[[col]], pd$cell_ID)[colnames(norm_mat)])

## Cluster ids sort numerically, and subcluster ids on both parts: 3.2 before 3.10 ##
ul <- local({
  u <- unique(cl)
  k <- vapply(
    strsplit(u, ".", fixed = TRUE),
    function(z) suppressWarnings(as.numeric(c(z, NA_character_))[1:2]),
    numeric(2)
  )
  if (any(is.na(k[1, ]))) sort(u) else u[order(k[1, ], k[2, ], na.last = FALSE)]
})

de <- suppressWarnings(findMarkers_one_vs_all(
  merged,
  method = "scran",
  expression_values = "normalized",
  cluster_column = col,
  min_feats = 10
))
top <- de[, head(.SD[order(ranking)], opt$top_n), by = cluster]

## Marker scores: z per gene across clusters, then averaged over each panel.
## Genes below min_detection are left out, their small sd turns noise into specificity ##
genes_all <- intersect(unique(mk$gene), rownames(norm_mat))
detected <- Matrix::rowMeans(norm_mat[genes_all, , drop = FALSE] > 0)
genes_use <- genes_all[detected >= opt$min_detection]
gm <- t(sapply(genes_all, function(g) tapply(norm_mat[g, ], cl, mean)))
gz <- t(scale(t(gm[genes_use, , drop = FALSE])))
gz[is.na(gz)] <- 0

panel_score <- function(M, gs) {
  gs <- intersect(gs, rownames(M))
  if (!length(gs)) {
    return(setNames(rep(NA_real_, ncol(M)), colnames(M)))
  }
  round(colMeans(M[gs, , drop = FALSE]), 3)
}
scores_mean <- t(sapply(markers, function(gs) panel_score(gm, gs)))
scores <- t(sapply(markers, function(gs) panel_score(gz, gs)))
n_used <- sapply(markers, function(gs) length(intersect(gs, genes_use)))

## Level 1: a compartment's genes score as one panel, not as the mean of its cell types ##
comps <- lapply(split(mk$gene, mk$compartment), unique)
scores_comp <- t(sapply(comps, function(gs) panel_score(gz, gs)))
if (!"immune_infiltrate" %in% rownames(scores_comp)) {
  stop(
    "diagnostics.R: no 'immune_infiltrate' compartment in ",
    opt$markers,
    ", present: ",
    paste(rownames(scores_comp), collapse = ", ")
  )
}
comp_used <- sapply(comps, function(gs) length(intersect(gs, genes_use)))
other <- scores_comp[
  setdiff(names(comps)[comp_used >= opt$min_panel_genes], "immune_infiltrate"),
  ul,
  drop = FALSE
]
if (!nrow(other)) {
  stop(
    "diagnostics.R: no non-immune compartment has >= ",
    opt$min_panel_genes,
    " usable markers"
  )
}
immune_score <- scores_comp["immune_infiltrate", ul]
immune_margin <- round(immune_score - apply(other, 2, max, na.rm = TRUE), 3)
top_nonimmune <- rownames(other)[apply(other, 2, which.max)]

## PTPRC detection: a high immune score without it is shared-gene bleed ##
frac_ptprc <- if ("PTPRC" %in% rownames(norm_mat)) {
  round(tapply(norm_mat["PTPRC", ] > 0, cl, mean)[ul], 3)
} else {
  NA_real_
}

## Only types with >= min_panel_genes usable markers can win ##
eligible <- names(markers)[n_used[names(markers)] >= opt$min_panel_genes]
if (!length(eligible)) {
  stop(
    "diagnostics.R: no cell type has >= ",
    opt$min_panel_genes,
    " markers detected above ",
    opt$min_detection
  )
}
message(
  "diagnostics.R [",
  col,
  "]: ",
  length(eligible),
  "/",
  length(markers),
  " cell types eligible for best_guess"
)
sc_e <- scores[eligible, , drop = FALSE]
rank_k <- function(x, k) {
  if (all(is.na(x))) NA_character_ else rownames(sc_e)[order(-x)][k]
}
best <- apply(sc_e, 2, rank_k, 1)
second <- apply(sc_e, 2, rank_k, 2)
best_score <- apply(sc_e, 2, function(x) {
  if (all(is.na(x))) NA_real_ else round(max(x, na.rm = TRUE), 3)
})
margin <- apply(sc_e, 2, function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  s <- sort(x, decreasing = TRUE, na.last = NA)
  round(s[1] - s[2], 3)
})

## No call when the winner is only the least negative score ##
no_call <- !is.na(best_score) & best_score < opt$min_score
if (any(no_call)) {
  message(
    "diagnostics.R [",
    col,
    "]: ",
    sum(no_call),
    " cluster(s) below --min_score, best_guess blank: ",
    paste(colnames(sc_e)[no_call], collapse = ", ")
  )
}
best[no_call] <- ""
topmk <- top[,
  .(top_markers = paste(head(feats[order(ranking)], 5), collapse = ";")),
  by = cluster
]
n_cell <- table(cl)

## Write ##
fwrite(de, file.path(opt$outdir, paste0(col, "_de_all.csv")))
fwrite(top, file.path(opt$outdir, paste0(col, "_de_top.csv")))
write.csv(
  scores_mean,
  file.path(opt$outdir, paste0(col, "_marker_scores_mean.csv"))
)
write.csv(scores, file.path(opt$outdir, paste0(col, "_marker_scores.csv")))
write.csv(
  data.frame(
    cell_type = names(markers),
    n_markers = lengths(markers)[names(markers)],
    n_present = sapply(markers, function(gs) length(intersect(gs, genes_all))),
    n_used = n_used[names(markers)]
  ),
  file.path(opt$outdir, paste0(col, "_marker_panel_coverage.csv")),
  row.names = FALSE
)
write.csv(
  scores_comp,
  file.path(opt$outdir, paste0(col, "_compartment_scores.csv"))
)
## cluster.R and subcluster.R already draw the Leiden and subcluster UMAPs ##
if (!grepl("^leiden_clus_res|^sub_", col)) {
  p <- plotUMAP(
    merged,
    dim_reduction_name = "umap",
    cell_color = col,
    point_size = 0.6,
    show_plot = FALSE,
    return_plot = TRUE
  )
  ggsave(
    file.path(opt$outdir, paste0("umap_", col, ".png")),
    p,
    width = 7,
    height = 6,
    dpi = 150,
    bg = "white"
  )
}
write.csv(
  data.frame(
    cluster_id = ul,
    n_cells = as.integer(n_cell[ul]),
    best_guess = best[ul],
    second_guess = second[ul],
    best_score = best_score[ul],
    margin = margin[ul],
    n_genes_used = ifelse(nzchar(best[ul]), n_used[best[ul]], NA_integer_),
    top_markers = topmk$top_markers[match(ul, as.character(topmk$cluster))],
    top_nonimmune = top_nonimmune,
    immune_score = round(immune_score, 3),
    immune_margin = immune_margin,
    immune_call = ifelse(immune_margin > 0, "immune", "non_immune"),
    frac_ptprc = frac_ptprc,
    label = ""
  ),
  file.path(opt$outdir, paste0("cluster_labels.", col, ".template.csv")),
  row.names = FALSE
)
