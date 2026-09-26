#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

merged <- loadGiotto(opt$input)
ids <- pDataDT(merged)$cell_ID
files <- sort(Sys.glob("insitutype_sup_*_cells.csv"))
if (!length(files)) {
  stop("attach_supervised.R: no insitutype_sup_*_cells.csv staged")
}

## Profile <name> adds insitutype_sup_<name> for the call and
## insitutype_sup_<name>_<column> for every other per-cell column ##
for (f in files) {
  nm <- sub("^insitutype_sup_(.*)_cells\\.csv$", "\\1", basename(f))
  d <- read.csv(f, stringsAsFactors = FALSE)
  i <- match(ids, d$cell_ID)
  if (anyNA(i)) {
    message(
      "attach_supervised.R [",
      nm,
      "]: ",
      sum(is.na(i)),
      " of ",
      length(ids),
      " cells absent from the CSV, set NA"
    )
  }
  new <- data.frame(cell_ID = ids)
  for (cl in setdiff(names(d), "cell_ID")) {
    nn <- switch(
      cl,
      celltype = paste0("insitutype_sup_", nm),
      celltype_refined = paste0("insitutype_sup_", nm, "_refined"),
      paste0("insitutype_sup_", nm, "_", cl)
    )
    new[[nn]] <- d[[cl]][i]
  }
  merged <- addCellMetadata(
    merged,
    by_column = TRUE,
    column_cell_ID = "cell_ID",
    new_metadata = new
  )
}

## Write ##
saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
