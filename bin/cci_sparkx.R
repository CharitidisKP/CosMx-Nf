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
  sample_col = arg("sample_col", "sample_id"),
  fdr = as.numeric(arg("fdr", "0.05")),
  cores = as.integer(arg("cores", "1"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(SPARK)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
raw <- getExpression(merged, values = "raw", output = "matrix")
locs <- getSpatialLocations(merged, output = "data.table")
setkey(locs, cell_ID)

## One test per sample, the kernels assume one continuous tissue. Cell type and
## log depth are regressed out, so a hit varies in space beyond composition ##
res <- rbindlist(lapply(sort(unique(pd[[opt$sample_col]])), function(s) {
  ids <- pd[
    get(opt$sample_col) == s & !is.na(get(opt$celltype_column)),
    cell_ID
  ]
  x <- raw[, ids, drop = FALSE]
  ct <- factor(pd[match(ids, cell_ID), get(opt$celltype_column)])
  covars <- cbind(
    stats::model.matrix(~ 0 + ct),
    log_depth = log(Matrix::colSums(x))
  )
  sk <- sparkx(
    x,
    as.matrix(locs[ids, .(sdimx, sdimy)]),
    X_in = covars,
    numCores = opt$cores,
    option = "mixture",
    verbose = FALSE
  )
  data.table(
    sample = s,
    gene = rownames(sk$res_mtest),
    combinedPval = sk$res_mtest$combinedPval,
    adjustedPval = sk$res_mtest$adjustedPval
  )
}))
res[, significant := adjustedPval < opt$fdr]
setorder(res, sample, adjustedPval)

## Write ##
fwrite(res, file.path(opt$outdir, "sparkx_svgs.csv"))
p <- ggplot(res[, .(n = sum(significant)), by = sample], aes(sample, n)) +
  geom_col(fill = "grey40") +
  labs(
    x = NULL,
    y = paste0("Genes at FDR < ", opt$fdr),
    title = "Spatially variable genes beyond cell-type composition"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(opt$outdir, "sparkx_svg_counts.png"),
  p,
  width = max(5, 2 + 0.6 * uniqueN(res$sample)),
  height = 4,
  dpi = 150,
  bg = "white"
)
