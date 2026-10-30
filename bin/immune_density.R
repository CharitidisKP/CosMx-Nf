#!/usr/bin/env Rscript
## Immune cells per mm2 of tissue and as a share of all cells, for each sample. Tissue area
## counts the bin_um squares that hold at least one cell centre: gaps narrower than a square
## count as tissue, the empty margins of a FOV do not ##
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
  celltype_column = arg("celltype_column", "cell_type"),
  sample_col = arg("sample_col", "sample_id"),
  immune_types = arg("immune_types"),
  um_per_px = as.numeric(arg("um_per_px", "0.12028")),
  bin_um = as.numeric(arg("bin_um", "50")),
  label_colours = arg("label_colours")
)
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}
immune <- split_arg(opt$immune_types)
stopifnot(length(immune) > 0, opt$um_per_px > 0, opt$bin_um > 0)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(ggplot2)
})
## The figure helpers sit beside this script ##
bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(bin, "cosmx_figures.R"))
set_fig_format(arg("fig_format", "png"))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
xy <- getSpatialLocations(merged, output = "data.table")
cells <- data.table(
  sample = pd[[opt$sample_col]],
  cell_type = as.character(pd[[opt$celltype_column]]),
  x = xy$sdimx[match(pd$cell_ID, xy$cell_ID)] * opt$um_per_px,
  y = xy$sdimy[match(pd$cell_ID, xy$cell_ID)] * opt$um_per_px
)[!is.na(x) & !is.na(y)]

## Tissue area per sample ##
area <- cells[,
  .(
    all_cells = .N,
    tissue_mm2 = uniqueN(paste(floor(x / opt$bin_um), floor(y / opt$bin_um))) *
      opt$bin_um^2 /
      1e6
  ),
  by = sample
]

## Count, density and share of each immune type, and of all of them together ##
types <- c(immune, "all_immune")
counts <- rbind(
  cells[cell_type %in% immune, .(n_cells = .N), by = .(sample, cell_type)],
  cells[cell_type %in% immune, .(cell_type = "all_immune", n_cells = .N), by = sample]
)
res <- CJ(sample = area$sample, cell_type = types, sorted = FALSE)
res <- merge(res, counts, by = c("sample", "cell_type"), all.x = TRUE)
res[is.na(n_cells), n_cells := 0L]
res <- merge(res, area, by = "sample")
res[, `:=`(
  per_mm2 = n_cells / tissue_mm2,
  pct_of_all = 100 * n_cells / all_cells,
  bin_um = opt$bin_um
)]
res <- res[order(sample, match(cell_type, types))]
setnames(res, "sample", opt$sample_col)

message(
  "immune_density.R: tissue mm2 ",
  paste0(area$sample, " ", round(area$tissue_mm2, 2), collapse = ", ")
)
fwrite(res, file.path(opt$outdir, "immune_density.csv"))

## Sample facts come from the samplesheet, as in composition.R. A person's samples are joined
## by a line; a sample without a partner (e.g. a control) is a grey point ##
ss <- fread(opt$samplesheet, colClasses = "character")
within <- split_arg(opt$compare_within)
miss <- setdiff(c("sample_id", opt$compare_by, within), names(ss))
if (length(miss)) {
  stop("immune_density.R: samplesheet has no column ", paste(miss, collapse = ", "))
}
ss[, grp := "all"]
if (length(within)) {
  ss[, grp := do.call(paste, c(.SD, sep = "/")), .SDcols = within]
}
ss[, subject := if (opt$pair_by %in% names(ss)) get(opt$pair_by) else sample_id]
ss[is.na(subject) | !nzchar(subject), subject := sample_id]
pl <- merge(
  res,
  ss[, .(sample_id, group = grp, level = get(opt$compare_by), subject)],
  by.x = opt$sample_col,
  by.y = "sample_id"
)
pl[, level := factor(level, levels = unique(c(opt$compare_baseline, sort(unique(level)))))]
style <- read_label_style(opt$label_colours)
names_ <- c(stats::setNames(display_name(immune, style), immune), all_immune = "All immune")
pl[, panel := factor(names_[cell_type], levels = names_[types])]

## One line per person from baseline to repeat biopsy ##
draw <- function(value, ylab, file, log) {
  p <- subject_lines(pl[, .(panel, subject, level, value = get(value))], ylab, log = log)
  save_fig(p, file.path(opt$outdir, file), 1.6 + 2.4 * length(types), 4)
}
draw("per_mm2", expression("Cells per mm"^2 * " of tissue"), "immune_density_per_mm2.png", TRUE)
draw("pct_of_all", "% of all cells", "immune_density_pct.png", FALSE)
