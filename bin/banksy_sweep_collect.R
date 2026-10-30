#!/usr/bin/env Rscript
## Scores every niche labelling of the BANKSY sweep, per biopsy, and draws them on the tissue.
## ari_vs_fov: agreement of niches with FOVs (1 = the niches are the FOVs, 0 = unrelated).
## fovs_one_niche: FOVs with at least 90% of their cells in one niche.
## celltype_fov_overlap: per typed cell type, how alike its niche mix is in each FOV and in the
## whole biopsy (1 = the same everywhere, 0 = FOV-specific niches).
## ari_vs_celltype: agreement of niches with cell types among typed cells (high = cell typing).
## coherence: share of each cell's nearest spatial neighbours in its niche (high = patches).
## focal_*: the niche holding the most cells of the focal cell type ##
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  labels = arg("labels", "labels"),
  times = arg("times", "times"),
  meta = arg("meta", "meta"),
  map_res = arg("map_res", "0.4,0.8"),
  untyped = arg("untyped", "Low_count,Unassigned"),
  focal = arg("focal", "Plasma_cell"),
  min_cells = as.integer(arg("min_cells", "30")),
  k_coherence = as.integer(arg("k_coherence", "15")),
  outdir = arg("outdir", ".")
)
split_arg <- function(s) {
  v <- trimws(strsplit(s, "[;,]")[[1]])
  v[nzchar(v)]
}
untyped <- split_arg(opt$untyped)
map_res <- as.numeric(split_arg(opt$map_res))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

## Adjusted Rand index in doubles ##
ari <- function(x, y) {
  tab <- table(x, y)
  s <- sum(choose(tab, 2))
  r <- sum(choose(rowSums(tab), 2))
  c <- sum(choose(colSums(tab), 2))
  e <- r * c / choose(sum(tab), 2)
  round((s - e) / ((r + c) / 2 - e), 3)
}

## Metrics of one labelling. typed: cells with a cell type outside the untyped labels ##
score <- function(niche, meta, typed, nn) {
  f <- data.table(fov = meta$fov, niche = niche)[, .N, by = .(fov, niche)][,
    .(top = max(N) / sum(N)),
    by = fov
  ]
  d <- data.table(ct = meta$cell_type, fov = meta$fov, niche = niche)[typed]
  big <- d[, .N, by = ct][N >= opt$min_cells, ct]
  pa <- d[ct %in% big, .N, by = .(ct, niche)][, p_all := N / sum(N), by = ct]
  pf <- d[ct %in% big, .N, by = .(ct, fov, niche)][,
    n_f := sum(N),
    by = .(ct, fov)
  ][
    n_f >= 10
  ][, p_f := N / n_f]
  ov <- pf[pa, on = .(ct, niche), nomatch = NULL][,
    .(o = sum(pmin(p_f, p_all)), n_f = n_f[1]),
    by = .(ct, fov)
  ][, .(o = weighted.mean(o, n_f), w = sum(n_f)), by = ct]
  fi <- which(meta$cell_type == opt$focal)
  ft <- if (length(fi)) sort(table(niche[fi]), decreasing = TRUE) else NULL
  data.table(
    niches = uniqueN(niche),
    ari_vs_fov = ari(niche, meta$fov),
    n_fovs = nrow(f),
    fovs_one_niche = sum(f$top >= 0.9),
    celltype_fov_overlap = if (nrow(ov)) {
      round(weighted.mean(ov$o, ov$w), 3)
    } else {
      NA_real_
    },
    ari_vs_celltype = ari(niche[typed], meta$cell_type[typed]),
    coherence = round(
      mean(rowMeans(matrix(niche[nn], nrow = nrow(nn)) == niche)),
      3
    ),
    focal_cells = length(fi),
    focal_top_pct = if (length(fi)) {
      round(100 * ft[[1]] / length(fi))
    } else {
      NA_real_
    },
    focal_niche_cells = if (length(fi)) {
      sum(niche == names(ft)[1])
    } else {
      NA_integer_
    }
  )
}

meta_files <- list.files(
  opt$meta,
  pattern = "_meta\\.csv\\.gz$",
  full.names = TRUE
)
label_files <- list.files(
  opt$labels,
  pattern = "_labels\\.csv\\.gz$",
  full.names = TRUE
)
samples <- sub("_meta\\.csv\\.gz$", "", basename(meta_files))
parse_col <- function(cl) {
  data.table(
    labelling = cl,
    variant = sub("_lam.*$", "", cl),
    lambda = as.numeric(sub("^.*_lam([0-9.]+).*$", "\\1", cl)),
    res = suppressWarnings(as.numeric(sub("^.*_res", "", cl)))
  )
}

metrics <- list()
for (i in seq_along(samples)) {
  s <- samples[i]
  meta <- fread(
    meta_files[i],
    colClasses = list(character = c("cell_ID", "fov"))
  )
  labs <- meta[,
    c("cell_ID", grep("^stored_lam", names(meta), value = TRUE)),
    with = FALSE
  ]
  for (f in label_files[startsWith(basename(label_files), paste0(s, "__"))]) {
    d <- fread(f, colClasses = "character")
    labs <- cbind(labs, d[match(labs$cell_ID, d$cell_ID), !"cell_ID"])
  }
  cols <- setdiff(names(labs), "cell_ID")
  typed <- !is.na(meta$cell_type) & !meta$cell_type %in% untyped
  nn <- dbscan::kNN(as.matrix(meta[, .(x, y)]), k = opt$k_coherence)$id
  metrics[[s]] <- rbindlist(lapply(cols, function(cl) {
    cbind(
      data.table(sample = s),
      parse_col(cl),
      score(labs[[cl]], meta, typed, nn)
    )
  }))
  message(
    "banksy_sweep_collect.R: ",
    s,
    " ",
    length(cols),
    " labellings scored"
  )

  ## Every variant and lambda at each map resolution, the stored niches on top, FOVs outlined ##
  box <- meta[,
    .(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y)),
    by = fov
  ]
  aspect <- diff(range(meta$y)) / diff(range(meta$x))
  for (r in map_res) {
    pc <- parse_col(cols)[variant == "stored" | res == r]
    if (!nrow(pc)) {
      next
    }
    long <- rbindlist(lapply(seq_len(nrow(pc)), function(j) {
      data.table(
        x = meta$x,
        y = meta$y,
        variant = pc$variant[j],
        lambda = pc$lambda[j],
        niche_no = as.integer(factor(labs[[pc$labelling[j]]]))
      )
    }))
    long[,
      variant := factor(
        variant,
        levels = c("stored", sort(setdiff(unique(variant), "stored")))
      )
    ]
    p <- ggplot(long, aes(x, y, colour = factor(niche_no))) +
      geom_point(size = 0.05) +
      geom_rect(
        data = box,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        colour = "black",
        linewidth = 0.2
      ) +
      facet_grid(variant ~ lambda, labeller = label_both) +
      coord_equal() +
      guides(colour = "none") +
      labs(
        x = NULL,
        y = NULL,
        title = sprintf("%s, niches at resolution %s, FOVs outlined", s, r)
      ) +
      theme_minimal(base_size = 11) +
      theme(axis.text = element_blank(), panel.grid = element_blank())
    ggsave(
      file.path(opt$outdir, sprintf("banksy_sweep_map_%s_res%s.png", s, r)),
      p,
      width = 4 * uniqueN(long$lambda) + 2,
      height = uniqueN(long$variant) * max(1, 4 * aspect) + 1.5,
      dpi = 110,
      bg = "white",
      limitsize = FALSE
    )
  }
}
met <- rbindlist(metrics)

times <- rbindlist(lapply(
  list.files(opt$times, pattern = "_time\\.csv$", full.names = TRUE),
  fread
))
if (nrow(times)) {
  met <- times[, .(sample, variant, n_genes, seconds)][
    met,
    on = .(sample, variant)
  ]
}
setcolorder(met, c("sample", "variant", "lambda", "res"))
setorder(met, sample, variant, lambda, res)
fwrite(met, file.path(opt$outdir, "banksy_sweep_metrics.csv"))
fwrite(times, file.path(opt$outdir, "banksy_sweep_times.csv"))

## Across biopsies, per variant, lambda and resolution ##
summ <- met[,
  .(
    biopsies = .N,
    ari_vs_fov = round(mean(ari_vs_fov), 3),
    fovs_one_niche_pct = round(100 * sum(fovs_one_niche) / sum(n_fovs)),
    celltype_fov_overlap = round(mean(celltype_fov_overlap, na.rm = TRUE), 3),
    ari_vs_celltype = round(mean(ari_vs_celltype), 3),
    coherence = round(mean(coherence), 3),
    niches_median = as.numeric(median(niches)),
    focal_top_pct = round(
      sum(focal_top_pct * focal_cells, na.rm = TRUE) /
        sum(focal_cells[!is.na(focal_top_pct)])
    )
  ),
  by = .(variant, lambda, res)
]
fwrite(summ, file.path(opt$outdir, "banksy_sweep_summary.csv"))

## One heatmap per metric: lambda against resolution, one panel per variant ##
for (mt in c(
  "ari_vs_fov",
  "celltype_fov_overlap",
  "ari_vs_celltype",
  "coherence",
  "niches_median"
)) {
  d <- summ[variant != "stored"]
  d[, value := get(mt)]
  p <- ggplot(d, aes(factor(res), factor(lambda), fill = value)) +
    geom_tile() +
    geom_text(aes(label = value), size = 2.6) +
    scale_fill_viridis_c(name = mt) +
    facet_wrap(~variant, ncol = 4) +
    labs(
      x = "Leiden resolution",
      y = "lambda",
      title = paste(mt, "(mean over biopsies)")
    ) +
    theme_minimal(base_size = 11)
  ggsave(
    file.path(opt$outdir, paste0("banksy_sweep_", mt, ".png")),
    p,
    width = 16,
    height = 2.4 * ceiling(uniqueN(d$variant) / 4) + 1.5,
    dpi = 130,
    bg = "white"
  )
}
message(
  "banksy_sweep_collect.R: ",
  nrow(met),
  " labellings over ",
  length(samples),
  " biopsies"
)
