#!/usr/bin/env Rscript
## BANKSY niches into the object. banksy_sample.R runs BANKSY on each sample with the pipeline's
## settings, the row of assets/banksy_variants.csv that banksy_settings names, and writes its
## niche labels. Here they are named <sample>:<niche>, matched across samples into niche types,
## written to the object, counted and drawn ##
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input"),
  labels = arg("labels", "labels"),
  genes = arg("genes", "genes"),
  variant = arg("variant"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON")),
  bg_mult = as.numeric(arg("bg_mult", "10")),
  untyped = arg("untyped", "Low_count,Unassigned"),
  k_geom = arg("k_geom", "15;30"),
  lambda_ladder = arg("lambda_ladder", "0.8"),
  use_agf = as.logical(arg("use_agf", "true")),
  k_nn = as.integer(arg("k_nn", "50")),
  sample_res = as.numeric(arg("sample_res", "0.4")),
  match_cut = as.numeric(arg("match_cut", "0.3")),
  sample_col = arg("sample_col", "sample_id"),
  celltype_column = arg("celltype_column", "cell_type"),
  um_per_px = as.numeric(arg("um_per_px", "0.12028")),
  focus = arg("focus"),
  label_colours = arg("label_colours")
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
## The figure helpers sit beside this script ##
bin <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(bin, "cosmx_figures.R"))
set_fig_format(arg("fig_format", "png"))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
## Lists from the settings table are separated by ; and lists from params by , ##
split_arg <- function(s) {
  v <- trimws(strsplit(s, "[;,]")[[1]])
  v[nzchar(v)]
}

## k_geom is one k per harmonic: the mean, then the AGF ##
k_geom <- as.integer(split_arg(opt$k_geom))
lambdas <- as.numeric(split_arg(opt$lambda_ladder))
if (!length(k_geom) || anyNA(k_geom) || !length(lambdas) || anyNA(lambdas)) {
  stop("banksy.R: --k_geom and --lambda_ladder must be lists of numbers")
}
if (!opt$use_agf) {
  k_geom <- k_geom[1]
}
untyped <- split_arg(opt$untyped)

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
samples <- sort(unique(pd[[opt$sample_col]]))
ct_of <- setNames(pd[[opt$celltype_column]], pd$cell_ID)

## Each sample's niches at the pipeline's resolution, one column per lambda. Rows in the
## object's order: samples sorted, cells as banksy_prep.R wrote them ##
per_sample <- rbindlist(lapply(samples, function(s) {
  f <- file.path(opt$labels, paste0(s, "__", opt$variant, "_labels.csv.gz"))
  if (!file.exists(f)) {
    stop("banksy.R: no niche labels for ", s, " (", basename(f), ")")
  }
  l <- fread(f, colClasses = "character")
  out <- data.table(cell_ID = l$cell_ID)
  for (lam in lambdas) {
    col <- paste0(opt$variant, "_lam", lam, "_res", opt$sample_res)
    if (!col %in% names(l)) {
      stop("banksy.R: no column ", col, " in ", basename(f))
    }
    out[[paste0("niche_sample_lam", lam)]] <- paste0(s, ":", l[[col]])
  }
  out
}))
genes_used <- rbindlist(lapply(samples, function(s) {
  fread(
    file.path(opt$genes, paste0(s, "__", opt$variant, "_genes.csv")),
    colClasses = c(sample = "character", gene = "character", best_cell_type = "character")
  )
}))

## Per-sample niches matched across samples by their typed cells: average linkage on
## 1 - correlation of their cell-type shares, cut at match_cut. Untyped labels are left out,
## because their share is large everywhere and makes unrelated niches correlate. A niche
## without typed cells stays a type of its own ##
types <- rbindlist(lapply(lambdas, function(lam) {
  col <- paste0("niche_sample_lam", lam)
  d <- per_sample[, .(niche = get(col), ct = ct_of[cell_ID])]
  n_all <- d[, .N, by = niche]
  d <- d[!is.na(ct) & !ct %in% untyped]
  m <- unclass(table(factor(d$niche, levels = n_all$niche), d$ct))
  sh <- m / pmax(rowSums(m), 1)
  typ <- if (nrow(sh) < 2L) {
    setNames(1L, rownames(sh))
  } else {
    dd <- 1 - suppressWarnings(cor(t(sh)))
    dd[is.na(dd)] <- 1
    cutree(hclust(as.dist(dd), method = "average"), h = opt$match_cut)
  }
  ## Types numbered by size, largest first ##
  size <- sort(tapply(n_all$N, typ[n_all$niche], sum), decreasing = TRUE)
  rank <- setNames(seq_along(size), names(size))
  data.table(
    lambda = lam,
    niche_sample = rownames(sh),
    sample = sub(":[^:]*$", "", rownames(sh)),
    type = paste0("type_", rank[as.character(typ)]),
    n_cells = n_all$N[match(rownames(sh), n_all$niche)],
    top_cell_types = apply(sh, 1L, function(v) {
      o <- order(-v)[seq_len(min(3L, sum(v > 0)))]
      paste0(colnames(sh)[o], " ", round(100 * v[o]), "%", collapse = ", ")
    })
  )
}))
for (lam in lambdas) {
  tl <- types[lambda == lam]
  per_sample[[paste0("niche_type_lam", lam)]] <- tl$type[match(
    per_sample[[paste0("niche_sample_lam", lam)]],
    tl$niche_sample
  )]
}

merged <- addCellMetadata(
  merged,
  by_column = TRUE,
  column_cell_ID = "cell_ID",
  new_metadata = as.data.frame(per_sample)
)

## Write ##
pdc <- pDataDT(merged)
cols <- setdiff(names(per_sample), "cell_ID")
mode_of <- function(cl) {
  fifelse(startsWith(cl, "niche_sample_"), "sample", "type")
}
lambda_of <- function(cl) sub("^niche_[a-z]+_lam([0-9.]+).*$", "\\1", cl)
long <- rbindlist(lapply(cols, function(cl) {
  pdc[, .(
    column = cl,
    lambda = lambda_of(cl),
    sample = get(opt$sample_col),
    cell_type = get(opt$celltype_column),
    niche = get(cl)
  )]
}))
fwrite(
  long[, .(n_cells = .N), by = .(column, lambda, sample, niche)],
  file.path(opt$outdir, "niche_by_sample.csv")
)
fwrite(
  long[, .(n_cells = .N), by = .(column, lambda, niche, cell_type)],
  file.path(opt$outdir, "niche_by_celltype.csv")
)
fwrite(types, file.path(opt$outdir, "niche_types.csv"))
fwrite(genes_used, file.path(opt$outdir, "banksy_genes.csv"))

## Typed cell types per niche, rows in the given order from the top ##
style <- read_label_style(opt$label_colours)
heat <- function(d, rows, file) {
  d <- d[!is.na(cell_type) & !cell_type %in% untyped]
  d[, share := n / sum(n), by = niche]
  d[, niche := factor(niche, levels = rev(rows))]
  d[, cell_type := factor(cell_type, levels = label_levels(cell_type, style))]
  h <- ggplot(d, aes(cell_type, niche, fill = share)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1)) +
    scale_x_discrete(labels = function(x) display_name(x, style)) +
    labs(x = NULL, y = NULL, fill = "Share") +
    theme_cosmx(10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  save_fig(
    h,
    file.path(opt$outdir, file),
    max(6, 2 + 0.35 * uniqueN(d$cell_type)),
    max(4, 1.5 + 0.3 * length(rows))
  )
}
for (lam in lambdas) {
  cs <- paste0("niche_sample_lam", lam)
  ns <- unique(per_sample[[cs]])
  heat(
    long[column == cs, .(n = .N), by = .(niche, cell_type)],
    ns[order(sub(":[^:]*$", "", ns), as.integer(sub("^.*:", "", ns)))],
    paste0("niche_celltype_heatmap_lam", lam, ".png")
  )
  heat(
    long[
      column == paste0("niche_type_lam", lam),
      .(n = .N),
      by = .(niche, cell_type)
    ],
    paste0("type_", seq_len(uniqueN(types[lambda == lam, type]))),
    paste0("niche_type_heatmap_lam", lam, ".png")
  )
}
write.csv(
  data.frame(
    column = cols,
    settings = opt$variant,
    mode = mode_of(cols),
    lambda = lambda_of(cols),
    resolution = ifelse(mode_of(cols) == "sample", opt$sample_res, NA_real_),
    n_niches = vapply(cols, function(cl) uniqueN(pdc[[cl]]), integer(1L)),
    k_geom = paste(k_geom, collapse = ","),
    use_agf = opt$use_agf,
    k_nn = opt$k_nn,
    bg_mult = opt$bg_mult,
    match_cut = ifelse(mode_of(cols) == "type", opt$match_cut, NA_real_)
  ),
  file.path(opt$outdir, "banksy_summary.csv"),
  row.names = FALSE
)
## Each sample's niches as the presentation draws them: every cell as its outline, a niche in
## its type's colour, named by its number and commonest cell types, with the focus cell types
## drawn over the niches in their own colour and a dark edge ##
focus <- split_arg(opt$focus)
outlines <- cell_outlines(merged)
for (lam in lambdas) {
  cs <- paste0("niche_sample_lam", lam)
  tl <- types[lambda == lam][, k := as.integer(sub("^.*:", "", niche_sample))]
  for (s in samples) {
    d <- pdc[get(opt$sample_col) == s, .(
      cell_ID,
      niche = get(cs),
      ct = as.character(get(opt$celltype_column)),
      fov = if ("fov" %in% names(pdc)) as.character(fov) else "all",
      low = FALSE
    )]
    key <- tl[sample == s][order(k)]
    shown <- intersect(focus, d$ct)
    d[, group := niche]
    d[ct %in% shown, group := ct]
    n_on <- d[ct %in% shown, .N, by = ct]
    d[, group := factor(group, levels = c(shown, key$niche_sample))]
    m <- biopsy_map(
      d,
      outlines[cell_ID %in% d$cell_ID],
      c(
        label_colours(shown, style),
        stats::setNames(niche_type_colour(key$type), key$niche_sample)
      ),
      c(
        stats::setNames(
          sprintf(
            "%s (%s cells)",
            display_name(shown, style),
            prettyNum(n_on$N[match(shown, n_on$ct)], big.mark = ",")
          ),
          shown
        ),
        stats::setNames(
          sprintf(
            "%d \u00b7 %s (%s cells)",
            key$k,
            display_in_text(key$top_cell_types, style),
            prettyNum(key$n_cells, big.mark = ",")
          ),
          key$niche_sample
        )
      ),
      opt$um_per_px,
      "BANKSY niche",
      title = s,
      on_top = shown,
      top_line = "grey10",
      top_lwd = 0.1,
      legend_w = 5.5
    )
    save_fig(
      m$plot,
      file.path(opt$outdir, paste0("niche_spatial_", s, "_", cs, ".png")),
      m$w,
      m$h,
      map = TRUE
    )
  }
}
saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
