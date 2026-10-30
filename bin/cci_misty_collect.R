#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input", "."),
  outdir = arg("outdir", ".")
)

suppressPackageStartupMessages({
  library(data.table)
  library(mistyR)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

## One key per sample, written beside its results folder by cci_misty.R ##
keys <- list.files(opt$input, "^misty_.*_key\\.csv$", full.names = TRUE)
if (!length(keys)) {
  stop("cci_misty_collect.R: no misty_*_key.csv in ", opt$input)
}
key <- rbindlist(lapply(keys, fread, colClasses = "character"))
fs <- unique(key[, .(folder, sample)])
folders <- setNames(fs$sample, fs$folder)
gs <- unique(key[, .(name, gene)])
gene_of <- setNames(gs$gene, gs$name)

## Back to sample ids and real gene names ##
r <- collect_results(file.path(opt$input, names(folders)))
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
