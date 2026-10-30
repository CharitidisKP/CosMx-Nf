#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input", "."),
  outdir = arg("outdir", "."),
  fdr = as.numeric(arg("fdr", "0.05")),
  niche_overlap_flag = as.numeric(arg("niche_overlap_flag", "2"))
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
setDTthreads(1)
## The figure helpers sit beside this script ##
bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(bin, "cosmx_figures.R"))
set_fig_format(arg("fig_format", "png"))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(opt$input, "_smide_de_results\\.csv$", full.names = TRUE)
if (!length(files)) {
  stop("smide_meta.R: no *_smide_de_results.csv in ", opt$input)
}
## fread can type a numeric column as text in one file (a long run of blanks, say), and
## rbindlist then makes it text for every file. The pooled columns are converted
## explicitly, and anything that is not a number is reported rather than silently dropped ##
num_cols <- c("ratio", "SE", "p.value", "ncells_1", "ncells_2")
read_one <- function(f) {
  d <- fread(f)
  txt <- intersect(num_cols, names(d)[vapply(d, is.character, logical(1))])
  for (cl in txt) {
    v <- d[[cl]]
    x <- suppressWarnings(as.numeric(v))
    bad <- sum(is.na(x) & !is.na(v) & nzchar(v))
    message(sprintf(
      "smide_meta.R: %s: %s read as text, converted%s",
      basename(f),
      cl,
      if (bad) sprintf(" (%d values not numeric, set to NA)", bad) else ""
    ))
    set(d, j = cl, value = x)
  }
  d
}
de <- rbindlist(lapply(files, read_one), fill = TRUE, use.names = TRUE)

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

## Bleed inside the tested niche. A contrast names one niche (against the rest or all) or
## two (pairwise); a row takes the highest overlap ratio among them, in its own sample ##
de[, niche_overlap := NA_real_]
ovf <- list.files(opt$input, "_smide_niche_overlap\\.csv$", full.names = TRUE)
if (length(ovf)) {
  ov <- rbindlist(lapply(ovf, fread, colClasses = "character"), fill = TRUE)
  ov[, ratio := suppressWarnings(as.numeric(ratio))]
  ctr <- unique(de[, .(sample, cell_type, contrast)])
  ## Whole-name tokens, so niche_1 never matches inside niche_10 ##
  tok <- regmatches(ctr$contrast, gregexpr("[A-Za-z0-9_.]+", ctr$contrast))
  pairs <- rbindlist(lapply(unique(ov$niche), function(n) {
    hit <- vapply(tok, function(t) n %in% t, logical(1))
    ctr[hit][, niche := n]
  }))
  if (nrow(pairs)) {
    x <- pairs[
      de[, .(sample, cell_type, contrast, target)],
      on = .(sample, cell_type, contrast),
      allow.cartesian = TRUE,
      nomatch = 0
    ]
    x[ov, r := i.ratio, on = .(sample, cell_type, niche, target)]
    x <- x[,
      .(niche_overlap = if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE)),
      by = .(sample, cell_type, contrast, target)
    ]
    de[x, niche_overlap := i.niche_overlap, on = .(sample, cell_type, contrast, target)]
  }
  message(sprintf(
    "smide_meta.R: niche overlap found for %d of %d rows",
    sum(!is.na(de$niche_overlap)),
    nrow(de)
  ))
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
      I2 = if (dfq > 0L && q > 0) max(0, (q - dfq) / q) * 100 else NA_real_,
      niche_overlap = if (all(is.na(niche_overlap))) {
        NA_real_
      } else {
        max(niche_overlap, na.rm = TRUE)
      }
    )
  },
  by = .(cell_type, spatial_model, result_component, contrast, target)
]
m[, niche_bleed := niche_overlap > opt$niche_overlap_flag]
m[,
  fdr := p.adjust(p_value, "BH"),
  by = .(cell_type, spatial_model, result_component)
]
setorder(m, cell_type, spatial_model, result_component, fdr)
hits <- m[
  fdr < opt$fdr,
  .(
    hits = .N,
    heterogeneous = sum(p_het < 0.05, na.rm = TRUE),
    flagged = sum(niche_bleed %in% TRUE)
  ),
  by = .(cell_type, spatial_model, result_component)
]

## Volcano on the niche-versus-rest contrasts, one panel per cell type ##
v <- if (any(m$result_component == "one.vs.rest")) {
  m[result_component == "one.vs.rest"]
} else {
  copy(m)
}
v[,
  call := fcase(
    fdr < opt$fdr & niche_bleed %in% TRUE,
    "hit, possible bleed",
    fdr < opt$fdr,
    "hit",
    default = "not significant"
  )
]
p <- ggplot(v, aes(log_fc / log(2), -log10(fdr), colour = call)) +
  geom_point(size = 0.4) +
  scale_colour_manual(
    values = c(
      hit = "firebrick",
      `hit, possible bleed` = "steelblue",
      `not significant` = "grey70"
    ),
    name = NULL
  ) +
  facet_wrap(~cell_type) +
  labs(x = "Pooled log2 fold change", y = "-log10 FDR") +
  theme_cosmx(10)

## Write ##
fwrite(m, file.path(opt$outdir, "smide_meta.csv"))
fwrite(hits, file.path(opt$outdir, "smide_meta_summary.csv"))
save_fig(
  p,
  file.path(opt$outdir, "smide_meta_volcano.png"),
  10,
  3 + 2.5 * ceiling(uniqueN(v$cell_type) / 3)
)
message(sprintf(
  "smide_meta.R: %d pooled rows, %d at FDR < %.2f",
  nrow(m),
  nrow(m[fdr < opt$fdr]),
  opt$fdr
))
