#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  sample_ids = arg("sample_ids"),
  x_padding = as.numeric(arg("x_padding", "1000")),
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
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

sids <- strsplit(opt$sample_ids, ",")[[1]]
objs <- lapply(paste0("input_giotto_", seq_along(sids)), loadGiotto)

merged <- joinGiottoObjects(
  gobject_list = objs,
  gobject_names = sids,
  join_method = "shift",
  x_padding = opt$x_padding
)

## Covariates are already loaded on each object ##
pd <- pDataDT(merged)

## If list_ID and sample_id disagree, the staging order did not match sample_ids ##
if (!all(pd$sample_id == pd$list_ID)) {
  stop("merge.R: sample_id does not match list_ID, staging order mismatch")
}
if (!"neg_mean" %in% names(pd)) {
  stop("merge.R: neg_mean did not survive the join")
}

## Write ##
saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
write.csv(
  pd[,
    .(n_cells = .N),
    by = .(sample_id, patient_id, treatment, timepoint, slide_id, batch)
  ],
  file.path(opt$outdir, "merge_cells_per_sample.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    n_samples = length(sids),
    n_cells = nrow(pd),
    n_rna = nrow(fDataDT(merged, feat_type = "rna")),
    n_negprobe = nrow(fDataDT(merged, feat_type = "negprobes")),
    sample_ids = paste(sids, collapse = ";")
  ),
  file.path(opt$outdir, "merge_summary.csv"),
  row.names = FALSE
)
