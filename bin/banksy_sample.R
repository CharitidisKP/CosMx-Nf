#!/usr/bin/env Rscript
## One biopsy through BANKSY with one settings set, a row of assets/banksy_variants.csv, at every
## lambda and resolution of that set. The spatial stage runs it with the pipeline's row, the sweep
## with every row. Writes one niche column per (lambda, resolution), and the genes it used ##
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  raw = arg("raw"),
  meta = arg("meta"),
  hvg = arg("hvg"),
  sample = arg("sample"),
  variant = arg("variant"),
  genes = arg("genes", "hvg"),
  values = arg("values", "counts"),
  k_geom = arg("k_geom", "15;30"),
  agf = as.logical(arg("agf", "true")),
  npcs = as.integer(arg("npcs", "50")),
  group = arg("group", "none"),
  lambdas = arg("lambdas", "0.8"),
  resolutions = arg("resolutions", "0.4"),
  k_nn = as.integer(arg("k_nn", "50")),
  cluster = arg("cluster", "banksy"),
  bg_mult = as.numeric(arg("bg_mult", "10")),
  untyped = arg("untyped", "Low_count,Unassigned"),
  min_cells = as.integer(arg("min_cells", "30")),
  seed = as.integer(arg("seed", "42")),
  outdir = arg("outdir", ".")
)
## Lists in the sweep table are separated by ; (commas belong to the CSV) ##
split_arg <- function(s) {
  v <- trimws(strsplit(s, "[;,]")[[1]])
  v[nzchar(v)]
}
k_geom <- as.integer(split_arg(opt$k_geom))
lambdas <- as.numeric(split_arg(opt$lambdas))
res <- as.numeric(split_arg(opt$resolutions))
stopifnot(
  length(k_geom) > 0,
  !anyNA(k_geom),
  length(lambdas) > 0,
  !anyNA(lambdas),
  length(res) > 0,
  !anyNA(res),
  !is.na(opt$agf)
)
if (!opt$agf) {
  k_geom <- k_geom[1]
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SpatialExperiment)
  library(Banksy)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[["elapsed"]]

raw <- readRDS(opt$raw)
meta <- fread(opt$meta, colClasses = list(character = c("cell_ID", "fov")))
stopifnot(identical(colnames(raw), meta$cell_ID))
hvg <- intersect(readLines(opt$hvg), rownames(raw))

## Genes: all, the pooled variable genes, or those of them whose best cell type reaches bg_mult
## times the biopsy's background (mean count per negative probe per cell) ##
genes <- switch(
  opt$genes,
  all = rownames(raw),
  hvg = hvg,
  hvg_bg = {
    untyped <- split_arg(opt$untyped)
    ct <- meta$cell_type
    n_ct <- table(ct[!is.na(ct) & !ct %in% untyped])
    types <- names(n_ct)[n_ct >= opt$min_cells]
    if (!length(types)) {
      stop("banksy_sample.R: ", opt$sample, " has no cell type with ", opt$min_cells, " cells")
    }
    means <- vapply(
      types,
      function(t) Matrix::rowMeans(raw[, which(ct == t), drop = FALSE]),
      numeric(nrow(raw))
    )
    fold <- apply(means, 1L, max) / mean(meta$neg_mean, na.rm = TRUE)
    kept <- intersect(hvg, rownames(raw)[fold >= opt$bg_mult])
    best_type <- types[apply(means[kept, , drop = FALSE], 1L, which.max)]
    x_background <- round(fold[kept], 1)
    kept
  },
  stop(
    "banksy_sample.R: unknown genes '",
    opt$genes,
    "', use all, hvg or hvg_bg"
  )
)
if (length(genes) < 50) {
  stop("banksy_sample.R: only ", length(genes), " genes for ", opt$sample, " with ", opt$variant)
}
message("banksy_sample.R: ", opt$sample, " ", length(genes), " genes")
## The genes used, with their best cell type and its multiple of background when hvg_bg chose them ##
gene_table <- data.table(
  sample = opt$sample,
  gene = genes,
  best_cell_type = if (opt$genes == "hvg_bg") best_type else NA_character_,
  x_background = if (opt$genes == "hvg_bg") unname(x_background) else NA_real_
)

## Values: counts over library size scaled to the mean library size without a log
## (normalizeCounts(log = FALSE)), or log2 with the median library size or 6000 as scale factor.
## The library size sums every gene ##
ls <- pmax(Matrix::colSums(raw), 1)
m <- raw[genes, , drop = FALSE]
sf <- switch(
  opt$values,
  counts = mean(ls) / ls,
  log = median(ls) / ls,
  log6000 = 6000 / ls,
  stop(
    "banksy_sample.R: unknown values '",
    opt$values,
    "', use counts, log or log6000"
  )
)
m@x <- m@x * rep(sf, diff(m@p))
if (opt$values != "counts") {
  m@x <- log2(m@x + 1)
}
rm(raw)

## BANKSY features and one PCA per lambda. group = fov z-scores each gene within each FOV ##
xy <- as.matrix(meta[, .(x, y)])
rownames(xy) <- meta$cell_ID
spe <- SpatialExperiment(
  assays = list(values = m),
  spatialCoords = xy,
  sample_id = opt$sample
)
spe$fov <- meta$fov
spe <- computeBanksy(
  spe,
  assay_name = "values",
  compute_agf = opt$agf,
  k_geom = k_geom,
  seed = opt$seed,
  verbose = FALSE
)
npcs <- min(opt$npcs, ncol(spe) - 1L)
spe <- runBanksyPCA(
  spe,
  use_agf = opt$agf,
  lambda = lambdas,
  npcs = npcs,
  group = if (opt$group == "none") NULL else opt$group,
  seed = opt$seed
)

## Leiden at every resolution: BANKSY's own shared-neighbour graph (clusterBanksy), or the old
## pipeline's unweighted k-nearest-neighbour graph with igraph (modularity) ##
labels <- data.table(cell_ID = meta$cell_ID)
col_of <- function(l, r) paste0(opt$variant, "_lam", l, "_res", r)
for (l in lambdas) {
  if (opt$cluster == "banksy") {
    before <- names(colData(spe))
    spe <- clusterBanksy(
      spe,
      use_agf = opt$agf,
      lambda = l,
      npcs = npcs,
      algo = "leiden",
      k_neighbors = opt$k_nn,
      resolution = res,
      seed = opt$seed
    )
    cl <- setdiff(names(colData(spe)), before)
    if (length(cl) != length(res)) {
      stop(
        "banksy_sample.R: clusterBanksy added ",
        length(cl),
        " columns, expected ",
        length(res)
      )
    }
    for (i in seq_along(res)) {
      labels[[col_of(l, res[i])]] <- as.character(colData(spe)[[cl[i]]])
    }
  } else if (opt$cluster == "knn_igraph") {
    pcs <- reducedDim(spe, sprintf("PCA_M%d_lam%s", as.integer(opt$agf), l))
    knn <- BiocNeighbors::findKNN(pcs, k = opt$k_nn)$index
    g <- igraph::simplify(igraph::graph_from_edgelist(
      cbind(rep(seq_len(nrow(knn)), each = ncol(knn)), as.vector(t(knn))),
      directed = FALSE
    ))
    for (r in res) {
      set.seed(opt$seed)
      labels[[col_of(l, r)]] <- as.character(igraph::membership(
        igraph::cluster_leiden(
          g,
          objective_function = "modularity",
          resolution = r
        )
      ))
    }
  } else {
    stop(
      "banksy_sample.R: unknown cluster '",
      opt$cluster,
      "', use banksy or knn_igraph"
    )
  }
  message(
    "banksy_sample.R: ",
    opt$sample,
    " ",
    opt$variant,
    " lambda ",
    l,
    " done"
  )
}

stem <- paste0(opt$sample, "__", opt$variant)
fwrite(labels, file.path(opt$outdir, paste0(stem, "_labels.csv.gz")))
fwrite(gene_table, file.path(opt$outdir, paste0(stem, "_genes.csv")))
fwrite(
  data.table(
    sample = opt$sample,
    variant = opt$variant,
    n_cells = nrow(meta),
    n_genes = length(genes),
    seconds = round(proc.time()[["elapsed"]] - t0)
  ),
  file.path(opt$outdir, paste0(stem, "_time.csv"))
)
