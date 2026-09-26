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
  exclude = arg("exclude"),
  min_cells = as.integer(arg("min_cells", "10")),
  cores = as.integer(arg("cores", "1")),
  seed = as.integer(arg("seed", "42"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(SingleCellExperiment)
  library(liana)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
norm <- getExpression(merged, values = "normalized", output = "matrix")
raw <- getExpression(merged, values = "raw", output = "matrix")
raw <- raw[, colnames(norm)]

## liana reads counts and logcounts and assumes log2, which is Giotto's default base ##
sce <- SingleCellExperiment(
  assays = list(counts = raw, logcounts = norm),
  colData = S4Vectors::DataFrame(
    cell_type = pd[[opt$celltype_column]][match(colnames(norm), pd$cell_ID)],
    row.names = colnames(norm)
  )
)
## Excluded labels carry no identity, so they neither send nor receive ##
exclude <- trimws(strsplit(opt$exclude, ",")[[1]])
sce <- sce[, !is.na(sce$cell_type) & !sce$cell_type %in% exclude]

## Contact weight: the share of target cells with a source cell among their Delaunay
## neighbours. Pairs that never touch get 0, which liana drops ##
ct_of <- setNames(pd[[opt$celltype_column]], pd$cell_ID)
net <- getSpatialNetwork(
  merged,
  name = "Delaunay_network",
  output = "networkDT"
)
edges <- rbind(
  net[, .(cell = from, nb = to)],
  net[, .(cell = to, nb = from)]
)[, .(cell, target = ct_of[cell], source = ct_of[nb])]
edges <- unique(edges[!is.na(target) & !is.na(source)])
types <- sort(unique(sce$cell_type))
adj <- merge(
  CJ(source = types, target = types),
  merge(
    edges[, .(touching = .N), by = .(target, source)],
    as.data.table(table(target = sce$cell_type)),
    by = "target"
  )[, .(source, target, adjacency = touching / N)],
  by = c("source", "target"),
  all.x = TRUE
)
adj[is.na(adjacency), adjacency := 0]

## Workers only spread the CellPhoneDB permutations ##
res <- liana_wrap(
  sce,
  idents_col = "cell_type",
  min_cells = opt$min_cells,
  cell.adj = adj,
  parallelize = opt$cores > 1,
  workers = opt$cores,
  seed = opt$seed,
  verbose = FALSE
)

agg <- as.data.table(liana_aggregate(res, verbose = FALSE))

## Write ##
fwrite(adj, file.path(opt$outdir, "liana_cell_adjacency.csv"))
fwrite(agg, file.path(opt$outdir, "liana_ligand_receptor.csv"))
need <- c(
  "source",
  "target",
  "ligand.complex",
  "receptor.complex",
  "aggregate_rank"
)
if (nrow(agg) && all(need %in% names(agg))) {
  top <- head(agg[order(aggregate_rank)], 30)
  top[, `:=`(
    pair = paste(source, "->", target),
    lr = paste(ligand.complex, "-", receptor.complex)
  )]
  p <- ggplot(top, aes(pair, lr, size = -log10(aggregate_rank))) +
    geom_point(colour = "firebrick") +
    labs(
      x = NULL,
      y = NULL,
      size = "-log10 rank",
      title = "Top 30 interactions by LIANA aggregate rank"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(
    file.path(opt$outdir, "liana_top_interactions.png"),
    p,
    width = max(7, 3 + 0.4 * uniqueN(top$pair)),
    height = max(5, 1.5 + 0.25 * uniqueN(top$lr)),
    dpi = 150,
    bg = "white"
  )
}
