#!/usr/bin/env Rscript
## Checks on discover's choices: the Leiden resolution, how many PCs and which genes ##
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
on <- function(x) tolower(x) %in% c("true", "1", "yes")
opt <- list(
  input = arg("input"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON")),
  reduction = arg("reduction"),
  leiden_res_ladder = arg("leiden_res_ladder"),
  chosen_resolutions = arg("chosen_resolutions"),
  dims_use = arg("dims_use", "1:30"),
  knn_k = as.integer(arg("knn_k", "15")),
  n_iterations = as.integer(arg("n_iterations", "1000")),
  markers = arg("markers"),
  n_hvgs = as.integer(arg("n_hvgs", "4000")),
  hvg_method = arg("hvg_method", "variance"),
  hvg_file = arg("hvg_file"),
  batch_column = arg("batch_column"),
  cells = as.integer(arg("cells", "20000")),
  stability = on(arg("stability", "false")),
  dims = on(arg("dims", "false")),
  genes = on(arg("genes", "false")),
  seed = as.integer(arg("seed", "42")),
  cores = as.integer(arg("cores", "1"))
)
if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
  options(giotto.py_path = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(ggplot2)
  library(ggraph)
})
## The helpers sit beside this script ##
bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(bin, "gene_selection.R"))
source(file.path(bin, "check_metrics.R"))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(opt$outdir, f)

parse_dims <- function(s) {
  p <- strsplit(s, ":")[[1]]
  if (length(p) == 2) {
    as.integer(p[1]):as.integer(p[2])
  } else {
    as.integer(strsplit(s, ",")[[1]])
  }
}
dims <- parse_dims(opt$dims_use)
res_tok <- trimws(strsplit(opt$leiden_res_ladder, ",")[[1]])
res_tok <- res_tok[nzchar(res_tok)]
chosen <- trimws(strsplit(opt$chosen_resolutions, ",")[[1]])
chosen <- chosen[nzchar(chosen)]

g <- loadGiotto(opt$input)

## Batch columns that vary, as one blocking factor over the given cells ##
batch_block <- function(meta) {
  v <- intersect(trimws(strsplit(opt$batch_column, ",")[[1]]), names(meta))
  v <- v[vapply(v, function(x) uniqueN(meta[[x]]) > 1L, logical(1))]
  if (length(v)) interaction(as.list(meta[, ..v]), drop = TRUE) else NULL
}

## The embedding CLUSTER used: --reduction, else harmony when the object has it, else pca ##
have_red <- tryCatch(
  as.character(GiottoClass::list_dim_reductions(g)$name),
  error = function(e) character()
)
red <- if (nzchar(opt$reduction)) {
  opt$reduction
} else if ("harmony" %in% have_red) {
  "harmony"
} else {
  "pca"
}
emb <- getDimReduction(
  g,
  reduction = "cells",
  reduction_method = red,
  name = red,
  output = "matrix"
)[, dims, drop = FALSE]
pd <- pDataDT(g)
pd <- pd[match(rownames(emb), cell_ID)]
message("cluster_checks.R: ", nrow(emb), " cells on ", red, ", dimensions ", opt$dims_use)

## The ladder's Leiden columns that the object carries ##
cols <- paste0("leiden_clus_res", res_tok)
res_tok <- res_tok[cols %in% names(pd)]
cols <- cols[cols %in% names(pd)]
if (length(cols) < 2L) {
  stop("cluster_checks.R: needs two or more leiden_clus_res columns, found ", length(cols))
}

## A seeded sample of cells for the silhouette, PC and gene checks ##
set.seed(opt$seed)
idx <- if (nrow(emb) > opt$cells) sort(sample(nrow(emb), opt$cells)) else seq_len(nrow(emb))

## Clustree across the ladder, nodes coloured by their SC3 stability ##
tree <- clustree::clustree(
  as.data.frame(pd[, ..cols]),
  prefix = "leiden_clus_res",
  node_colour = "sc3_stability"
)
ggsave(out("clustree.png"), tree, width = 10, height = 12, dpi = 150, bg = "white")
nodes <- as.data.table(as.data.frame(tryCatch(tree$data, error = function(e) tree@data)))
nodes[, resolution := res_tok[match(as.numeric(as.character(leiden_clus_res)), as.numeric(res_tok))]]
nodes <- nodes[, .(resolution, cluster, size, sc3_stability = round(sc3_stability, 4))]
fwrite(nodes, out("clustree_nodes.csv"))

## Mean silhouette width per resolution, in the clustering embedding, on the sampled cells ##
sil <- vapply(
  cols,
  function(cc) {
    mean(bluster::approxSilhouette(emb[idx, , drop = FALSE], factor(pd[[cc]][idx]))$width)
  },
  numeric(1)
)
res_tab <- data.table(
  resolution = res_tok,
  n_clusters = vapply(cols, function(cc) uniqueN(pd[[cc]]), integer(1)),
  clusters_under_100_cells = vapply(cols, function(cc) sum(table(pd[[cc]]) < 100L), integer(1)),
  median_sc3_stability = vapply(
    res_tok,
    function(r) round(stats::median(nodes[resolution == r, sc3_stability]), 3),
    numeric(1)
  ),
  silhouette = round(sil, 3)
)

## Gene sets by each method on the sampled cells, and how far each one moves the PCA ##
if (opt$genes) {
  sub <- subsetGiotto(g, cell_ids = rownames(emb)[idx])
  nm <- as.matrix(getExpression(sub, feat_type = "rna", values = "normalized", output = "matrix"))
  raw <- getExpression(sub, feat_type = "rna", values = "raw", output = "matrix")[, colnames(nm)]
  spd <- pDataDT(sub)[match(colnames(nm), cell_ID)]
  canon <- intersect(unique(read.csv(opt$markers)$gene), rownames(nm))
  ## NanoString's method runs inside each cluster of the first chosen resolution ##
  grp <- intersect(paste0("leiden_clus_res", c(chosen, res_tok)), cols)[1]
  sets <- list(
    variance = select_genes(nm, "variance", opt$n_hvgs),
    variance_batch = select_genes(nm, "variance_batch", opt$n_hvgs, block = batch_block(spd)),
    all = select_genes(nm, "all"),
    per_group = per_group_genes(raw, spd[[grp]], min_cells = 200L)
  )
  if (nzchar(opt$hvg_file)) sets$file <- select_genes(nm, "file", gene_file = opt$hvg_file)
  run_set <- sets[[opt$hvg_method]]
  for (m in names(sets)) writeLines(sets[[m]], out(paste0("genes_", m, ".txt")))
  ## PCA scores on the sampled cells from a gene set with the markers added, as the pipeline runs it ##
  pcs_for <- function(f) {
    f <- union(f, canon)
    k <- min(max(dims), length(f) - 1L)
    s2 <- runPCA(sub, feat_type = "rna", feats_to_use = f, ncp = k, name = "pca_genes")
    getDimReduction(s2, reduction = "cells", reduction_method = "pca", name = "pca_genes", output = "matrix")
  }
  ref_pc <- pcs_for(run_set)
  emb_tab <- rbindlist(lapply(names(sets), function(m) {
    if (length(union(sets[[m]], canon)) < 50L) return(NULL)
    cs <- subspace_cosines(ref_pc, pcs_for(sets[[m]])[rownames(ref_pc), , drop = FALSE])
    data.table(method = m, mean_cos = round(mean(cs), 3), min_cos = round(min(cs), 3), shared_dims = sum(cs >= 0.9))
  }))
  tab <- as.data.table(gene_set_table(sets, run_set, canon))
  if (nrow(emb_tab)) tab <- merge(tab, emb_tab, by = "method", all.x = TRUE, sort = FALSE)
  tab[, run_method := method == opt$hvg_method]
  fwrite(tab, out("genes_methods.csv"))
  message("cluster_checks.R: gene sets ", paste(sprintf("%s %d", tab$method, tab$n_genes), collapse = ", "))
  rm(sub, nm, raw, ref_pc)
  invisible(gc())
}

## How many PCs the data supports: a 100-PC PCA on the run's genes, with noise from its tail ##
if (opt$dims) {
  pca <- getDimReduction(g, reduction = "cells", reduction_method = "pca", name = "pca", output = "dimObj")
  feats <- tryCatch(rownames(pca@misc$loadings), error = function(e) NULL)
  ## Without stored loadings, the run's genes again, picked as integrate.R picked them ##
  if (!length(feats)) {
    nm <- as.matrix(getExpression(g, feat_type = "rna", values = "normalized", output = "matrix"))
    blk <- batch_block(pDataDT(g)[match(colnames(nm), cell_ID)])
    feats <- union(
      select_genes(nm, opt$hvg_method, opt$n_hvgs, block = blk, gene_file = opt$hvg_file),
      intersect(unique(read.csv(opt$markers)$gene), rownames(nm))
    )
    rm(nm)
  }
  ncp <- min(100L, length(feats) - 1L)
  g <- runPCA(g, feat_type = "rna", feats_to_use = feats, ncp = ncp, name = "pca_check")
  pcs <- getDimReduction(g, reduction = "cells", reduction_method = "pca", name = "pca_check", output = "matrix")
  pc_var <- apply(pcs, 2, stats::var)
  sigma <- sqrt(mean(pc_var[seq(max(1L, ncp - 19L), ncp)]))
  set.seed(opt$seed)
  at <- sort(sample(nrow(pcs), min(opt$cells, nrow(pcs))))
  ## Noise sd for distances with noise at both ends and at one end ##
  sigmas <- c(both_ends = sqrt(2) * sigma, one_end = sigma)
  dims_tab <- as.data.table(dims_estimate(pcs, at, c(10L, 20L, 50L), sigmas, opt$cores))
  dims_tab[, `:=`(estimate = round(estimate, 1), noise_sd = round(sigma, 4), n_pcs = ncp, dims_use = max(dims))]
  fwrite(dims_tab, out("dims_estimate.csv"))
  pv <- data.table(pc = seq_along(pc_var), variance = pc_var)
  fwrite(pv, out("dims_pc_variance.csv"))
  p <- ggplot(pv, aes(pc, variance)) +
    annotate("rect", xmin = ncp - 19.5, xmax = ncp + 0.5, ymin = min(pc_var), ymax = Inf, alpha = 0.12) +
    geom_vline(xintercept = max(dims), linetype = "dashed", colour = "grey40") +
    geom_point(colour = "#4C7A8C", size = 1.4) +
    scale_y_log10() +
    labs(x = "PC", y = "Variance") +
    theme_minimal(base_size = 12)
  ggsave(out("dims_pc_variance.png"), p, width = 7, height = 4, dpi = 150, bg = "white")
  message("cluster_checks.R: intrinsic dimension ", paste(dims_tab[noise == "both_ends", estimate], collapse = ", "))
  rm(pcs)
  invisible(gc())
}

## Stability: three 80% samples, each with its own graph, clustered at every resolution ##
if (opt$stability) {
  n_full <- stats::setNames(lapply(cols, function(cc) table(as.character(pd[[cc]]))), res_tok)
  jac <- stats::setNames(vector("list", length(res_tok)), res_tok)
  for (b in 1:3) {
    set.seed(opt$seed + b)
    s <- subsetGiotto(g, cell_ids = sample(pd$cell_ID, round(0.8 * nrow(pd))))
    s <- createNearestNetwork(
      s,
      feat_type = "rna",
      type = "sNN",
      dim_reduction_to_use = red,
      dim_reduction_name = red,
      dimensions_to_use = dims,
      k = opt$knn_k,
      name = "NN.check"
    )
    for (i in seq_along(res_tok)) {
      nm_col <- paste0("check_res", res_tok[i])
      s <- doLeidenCluster(
        s,
        feat_type = "rna",
        nn_network_to_use = "sNN",
        network_name = "NN.check",
        resolution = as.numeric(res_tok[i]),
        n_iterations = opt$n_iterations,
        name = nm_col,
        seed_number = opt$seed + b
      )
      sp <- pDataDT(s)
      j <- best_jaccard(sp[[cols[i]]], sp[[nm_col]])
      jac[[i]] <- cbind(jac[[i]], j[names(n_full[[i]])])
    }
    message("cluster_checks.R: stability sample ", b, " of 3 done")
    rm(s)
    invisible(gc())
  }
  mean_j <- lapply(jac, function(m) rowMeans(m, na.rm = TRUE))
  fwrite(
    rbindlist(lapply(seq_along(res_tok), function(i) {
      data.table(
        resolution = res_tok[i],
        cluster = names(mean_j[[i]]),
        n_cells = as.integer(n_full[[i]][names(mean_j[[i]])]),
        mean_jaccard = round(mean_j[[i]], 3)
      )
    })),
    out("stability_clusters.csv")
  )
  res_tab <- cbind(
    res_tab,
    rbindlist(lapply(seq_along(res_tok), function(i) {
      j <- mean_j[[i]][is.finite(mean_j[[i]])]
      stability_row(j, n_full[[i]])
    }))
  )
  ## The resolution later steps take when no clustering is named ##
  pick <- pick_resolution(res_tab$resolution, res_tab$median_jaccard)
  res_tab[, picked := resolution %in% pick]
  if (!is.na(pick)) {
    writeLines(cols[match(pick, res_tok)], out("picked_resolution.txt"))
    message("cluster_checks.R: picked resolution ", pick)
  }
}
fwrite(res_tab, out("resolution_checks.csv"))

## Each measure against the resolution, the chosen resolutions dashed and the pick solid ##
long <- melt(
  res_tab[, intersect(
    names(res_tab),
    c("resolution", "silhouette", "median_sc3_stability", "median_jaccard")
  ), with = FALSE],
  id.vars = "resolution",
  variable.name = "measure"
)
long[, resolution := as.numeric(resolution)]
vl <- as.numeric(intersect(chosen, res_tok))
pk <- numeric(0)
if ("picked" %in% names(res_tab)) pk <- as.numeric(res_tab[picked == TRUE, resolution])
p <- ggplot(long, aes(resolution, value)) +
  (if (length(vl)) geom_vline(xintercept = vl, linetype = "dashed", colour = "grey60")) +
  (if (length(pk)) geom_vline(xintercept = pk, colour = "grey30")) +
  geom_line(colour = "#4C7A8C") +
  geom_point(colour = "#4C7A8C", size = 2) +
  facet_wrap(~measure, ncol = 1, scales = "free_y") +
  labs(x = "Leiden resolution", y = NULL) +
  theme_minimal(base_size = 12)
ggsave(out("resolution_checks.png"), p, width = 7, height = 2.4 * uniqueN(long$measure) + 0.6, dpi = 150, bg = "white")

message("cluster_checks.R: done")
