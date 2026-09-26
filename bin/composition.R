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
  samplesheet = arg("samplesheet"),
  compare_by = arg("compare_by", "timepoint"),
  compare_baseline = arg("compare_baseline", "T0"),
  compare_within = arg("compare_within", "treatment"),
  pair_by = arg("pair_by", "subject_id"),
  exclude = arg("exclude"),
  sample_col = arg("sample_col", "sample_id")
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)

## Sample facts come from the samplesheet, so correcting one there needs no rerun upstream ##
ss <- fread(opt$samplesheet, colClasses = "character")
ss <- ss[sample_id %in% pd[[opt$sample_col]]]
within <- split_arg(opt$compare_within)
miss <- setdiff(c(opt$compare_by, within), names(ss))
if (length(miss)) {
  stop(
    "composition.R: samplesheet has no column ",
    paste(miss, collapse = ", ")
  )
}
ss[, grp := "all"]
if (length(within)) {
  ss[, grp := do.call(paste, c(.SD, sep = "/")), .SDcols = within]
}
ss[,
  subject := if (opt$pair_by %in% names(ss)) get(opt$pair_by) else NA_character_
]
ss[!nzchar(subject), subject := NA_character_]

## One row per sample and cell type, zeros included. Untyped labels count towards a sample's
## cells but not towards the typed shares, so a change in their share cannot move the rest ##
tab <- as.data.table(table(
  sample_id = pd[[opt$sample_col]],
  cell_type = pd$cell_type
))
setnames(tab, "N", "n_cells")
tab[, typed := !cell_type %in% split_arg(opt$exclude)]
tab[, share_of_all := round(n_cells / sum(n_cells), 5), by = sample_id]
tab[, share_of_typed := round(n_cells / sum(n_cells[typed]), 5), by = sample_id]
tab[(!typed), share_of_typed := NA_real_]
tab[,
  log2_prop := log2((n_cells + 0.5) / sum((n_cells + 0.5)[typed])),
  by = sample_id
]
tab <- merge(
  ss[, .(sample_id, group = grp, level = get(opt$compare_by), subject)],
  tab,
  by = "sample_id"
)

## Each level against the baseline, inside each group. A handful of subjects per group is too
## few for a test, so the change is per subject. Half a cell keeps zeros finite ##
chg <- rbindlist(lapply(unique(tab$group), function(g) {
  s <- tab[group == g & !is.na(subject) & typed]
  base <- s[level == opt$compare_baseline]
  rbindlist(lapply(setdiff(unique(s$level), opt$compare_baseline), function(l) {
    d <- merge(
      s[level == l],
      base[, .(
        subject,
        cell_type,
        baseline_sample = sample_id,
        base = log2_prop
      )],
      by = c("subject", "cell_type")
    )
    d[, .(
      comparison = paste0(g, ": ", l, " vs ", opt$compare_baseline),
      subject,
      cell_type,
      baseline_sample,
      sample_id,
      log2_change = round(log2_prop - base, 3)
    )]
  }))
}))

## Write ##
fwrite(
  tab[, .(
    sample_id,
    group,
    level,
    subject,
    cell_type,
    n_cells,
    share_of_all,
    share_of_typed
  )],
  file.path(opt$outdir, "composition_by_sample.csv")
)
fwrite(chg, file.path(opt$outdir, "composition_paired_change.csv"))
p <- ggplot(tab, aes(sample_id, share_of_all, fill = cell_type)) +
  geom_col() +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Share of all cells", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(opt$outdir, "composition_barplot.png"),
  p,
  width = max(6, 3 + 0.6 * uniqueN(tab$sample_id)),
  height = 6,
  dpi = 150,
  bg = "white"
)
if (nrow(chg)) {
  q <- ggplot(chg, aes(log2_change, cell_type, colour = subject)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_point(size = 2) +
    facet_wrap(~comparison) +
    labs(x = "log2 change in share of typed cells", y = NULL, colour = NULL) +
    theme_minimal(base_size = 11)
  ggsave(
    file.path(opt$outdir, "composition_paired_change.png"),
    q,
    width = max(6, 3 + 3 * uniqueN(chg$comparison)),
    height = max(4, 1 + 0.3 * uniqueN(chg$cell_type)),
    dpi = 150,
    bg = "white"
  )
}
