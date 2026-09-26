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
  markers = arg("markers"),
  sample_col = arg("sample_col", "sample_id"),
  juxta_thr = as.numeric(arg("juxta_thr", "150")),
  para_ls = arg("para_ls", "400,800,1600"),
  min_cells = as.integer(arg("min_cells", "50")),
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
  library(mistyR)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

## Views and models run in parallel only under a future plan ##
future::plan(future::multisession, workers = opt$cores)

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
norm <- getExpression(merged, values = "normalized", output = "matrix")
locs <- getSpatialLocations(merged, output = "data.table")
setkey(locs, cell_ID)
para_ls <- as.numeric(strsplit(opt$para_ls, ",")[[1]])

## Targets become formula terms, so IGHG1/2 and the like need syntactic names ##
feats <- intersect(unique(read.csv(opt$markers)$gene), rownames(norm))
gene_of <- setNames(feats, make.names(feats))

## One joint model per sample: intra, juxta and every paraview scale ##
folders <- character()
for (s in sort(unique(pd[[opt$sample_col]]))) {
  ids <- pd$cell_ID[pd[[opt$sample_col]] == s]
  if (length(ids) < opt$min_cells) {
    message(
      "cci_misty.R: ",
      s,
      " has ",
      length(ids),
      " cells < min_cells, skipped"
    )
    next
  }
  ex <- t(as.matrix(norm[feats, ids, drop = FALSE]))
  colnames(ex) <- names(gene_of)

  ## run_misty stops on a target with no variance ##
  ex <- as.data.frame(ex[, matrixStats::colVars(ex) > 0, drop = FALSE])
  xy <- as.data.frame(locs[ids, .(x = sdimx, y = sdimy)])
  views <- create_initial_view(ex) |>
    add_juxtaview(xy, neighbor.thr = opt$juxta_thr, verbose = FALSE)
  for (l in para_ls) {
    views <- add_paraview(
      views,
      xy,
      l = l,
      zoi = opt$juxta_thr,
      verbose = FALSE
    )
  }
  f <- file.path(opt$outdir, paste0("misty_", gsub("[^A-Za-z0-9]+", "_", s)))
  run_misty(views, results.folder = f, seed = opt$seed)
  folders[basename(f)] <- s
}
if (!length(folders)) {
  stop("cci_misty.R: no sample had >= ", opt$min_cells, " cells")
}

## Back to sample ids and real gene names ##
r <- collect_results(file.path(opt$outdir, names(folders)))
tidy <- function(x, gene_cols) {
  d <- as.data.table(x)
  d[, sample := folders[basename(sample)]]
  for (gc in intersect(gene_cols, names(d))) {
    d[, (gc) := unname(gene_of[get(gc)])]
  }
  d[]
}

imp <- tidy(r$improvements, "target")

## Write ##
fwrite(imp, file.path(opt$outdir, "misty_improvements.csv"))
gain <- imp[measure == "gain.R2"]
if (nrow(gain)) {
  top <- gain[, .(m = mean(value, na.rm = TRUE)), by = target][order(-m)]
  top <- head(top$target, 25)
  p <- ggplot(
    gain[target %in% top],
    aes(sample, factor(target, levels = rev(top)), fill = value)
  ) +
    geom_tile() +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
    labs(
      x = NULL,
      y = NULL,
      fill = "Gain in R2",
      title = "Expression explained by the tissue context, top 25 targets"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(
    file.path(opt$outdir, "misty_gain_r2.png"),
    p,
    width = max(6, 3 + 0.6 * uniqueN(gain$sample)),
    height = 7,
    dpi = 150,
    bg = "white"
  )
}
fwrite(
  tidy(r$contributions, "target"),
  file.path(opt$outdir, "misty_contributions.csv")
)
fwrite(
  tidy(r$importances, c("Predictor", "Target")),
  file.path(opt$outdir, "misty_importances.csv")
)
