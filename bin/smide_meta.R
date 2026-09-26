#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input", "."),
  outdir = arg("outdir", "."),
  fdr = as.numeric(arg("fdr", "0.05"))
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
setDTthreads(1)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(opt$input, "_smide_de_results\\.csv$", full.names = TRUE)
if (!length(files)) {
  stop("smide_meta.R: no *_smide_de_results.csv in ", opt$input)
}
de <- rbindlist(lapply(files, fread), fill = TRUE, use.names = TRUE)

## emmeans rows are level estimates, not contrasts. Only contrasts pool ##
de <- de[
  result_component %in%
    c("pairwise", "one.vs.rest", "one.vs.all") &
    is.finite(ratio) &
    ratio > 0 &
    is.finite(SE) &
    SE > 0
]
if (!nrow(de)) {
  stop("smide_meta.R: no contrast rows with a finite ratio and SE")
}

## Fixed-effect pool on the log scale. Delta method: se(log ratio) = SE / ratio ##
if (!"ncells_1" %in% names(de)) {
  de[, ncells_1 := NA_real_]
}
de[, `:=`(beta = log(ratio), se = SE / ratio)]
de[, w := 1 / se^2]
m <- de[,
  {
    bp <- sum(w * beta) / sum(w)
    sp <- sqrt(1 / sum(w))
    q <- sum(w * (beta - bp)^2)
    dfq <- .N - 1L
    .(
      n_samples = .N,
      samples = paste(sort(unique(sample)), collapse = ";"),
      n_cells = sum(ncells_1, na.rm = TRUE),
      log_fc = bp,
      se = sp,
      z = bp / sp,
      p_value = 2 * stats::pnorm(-abs(bp / sp)),
      fold_change = exp(bp),
      Q = q,
      df_Q = dfq,
      p_het = if (dfq > 0L) {
        stats::pchisq(q, dfq, lower.tail = FALSE)
      } else {
        NA_real_
      },
      I2 = if (dfq > 0L && q > 0) max(0, (q - dfq) / q) * 100 else NA_real_
    )
  },
  by = .(cell_type, spatial_model, result_component, contrast, target)
]
m[,
  fdr := p.adjust(p_value, "BH"),
  by = .(cell_type, spatial_model, result_component)
]
setorder(m, cell_type, spatial_model, result_component, fdr)
hits <- m[
  fdr < opt$fdr,
  .(hits = .N, heterogeneous = sum(p_het < 0.05, na.rm = TRUE)),
  by = .(cell_type, spatial_model, result_component)
]

## Volcano on the niche-versus-rest contrasts, one panel per cell type ##
v <- if (any(m$result_component == "one.vs.rest")) {
  m[result_component == "one.vs.rest"]
} else {
  m
}
p <- ggplot(
  v,
  aes(log_fc / log(2), -log10(fdr), colour = fdr < opt$fdr)
) +
  geom_point(size = 0.4) +
  scale_colour_manual(
    values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
    guide = "none"
  ) +
  facet_wrap(~cell_type) +
  labs(
    x = "Pooled log2 fold change",
    y = "-log10 FDR",
    title = "smiDE meta-analysis, niche against the rest"
  ) +
  theme_minimal(base_size = 10)

## Write ##
fwrite(m, file.path(opt$outdir, "smide_meta.csv"))
fwrite(hits, file.path(opt$outdir, "smide_meta_summary.csv"))
ggsave(
  file.path(opt$outdir, "smide_meta_volcano.png"),
  p,
  width = 10,
  height = 3 + 2.5 * ceiling(uniqueN(v$cell_type) / 3),
  dpi = 150,
  bg = "white"
)
message(sprintf(
  "smide_meta.R: %d pooled rows, %d at FDR < %.2f",
  nrow(m),
  nrow(m[fdr < opt$fdr]),
  opt$fdr
))
