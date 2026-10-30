#!/usr/bin/env Rscript
## Inputs for banksy_sample.R, in the spatial stage and the sweep alike: the top variable genes
## of each biopsy, pooled over all biopsies (Banksy multi-sample vignette), and per biopsy its raw
## counts and cell metadata, cells in the object's order ##
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON")),
  samples = arg("samples"),
  hvgs = as.integer(arg("hvgs", "2000")),
  celltype_column = arg("celltype_column", "cell_type"),
  sample_col = arg("sample_col", "sample_id")
)
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
samples <- sort(unique(pd[[opt$sample_col]]))
if (nzchar(opt$samples)) {
  miss <- setdiff(split_arg(opt$samples), samples)
  if (length(miss)) {
    stop(
      "banksy_prep.R: not in the object: ",
      paste(miss, collapse = ", ")
    )
  }
  samples <- split_arg(opt$samples)
}

## Top variable genes of each biopsy on its own, pooled over every biopsy in the object ##
norm <- getExpression(merged, values = "normalized", output = "matrix")
hvg <- Reduce(
  union,
  lapply(split(pd$cell_ID, pd[[opt$sample_col]]), function(ids) {
    suppressWarnings(scran::getTopHVGs(
      scran::modelGeneVar(norm[, ids]),
      n = opt$hvgs
    ))
  })
)
writeLines(hvg, file.path(opt$outdir, "hvg_pooled.txt"))
message("banksy_prep.R: ", length(hvg), " pooled variable genes")
rm(norm)
invisible(gc())

## Per biopsy: raw counts of every gene, and position, FOV, cell type, background and the
## stored per-biopsy niches ##
raw <- getExpression(merged, values = "raw", output = "matrix")
locs <- getSpatialLocations(merged, output = "data.table")
stored <- intersect(c("niche_sample_lam0.2", "niche_sample_lam0.8"), names(pd))
for (s in samples) {
  m <- pd[get(opt$sample_col) == s]
  meta <- data.table(
    cell_ID = m$cell_ID,
    fov = as.character(m$fov),
    cell_type = m[[opt$celltype_column]],
    neg_mean = m$neg_mean
  )
  meta[locs, on = "cell_ID", `:=`(x = i.sdimx, y = i.sdimy)]
  for (cl in stored) {
    meta[[sub("^niche_sample_", "stored_", cl)]] <- m[[cl]]
  }
  saveRDS(
    as(raw[, meta$cell_ID, drop = FALSE], "CsparseMatrix"),
    file.path(opt$outdir, paste0(s, "_raw.rds"))
  )
  fwrite(meta, file.path(opt$outdir, paste0(s, "_meta.csv.gz")))
  message("banksy_prep.R: ", s, " ", nrow(meta), " cells")
}
