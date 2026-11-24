## The report's own figures: ggplot2 in the pipeline's deck style, drawn again only when their inputs change ##
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
if (nzchar(RUN$fig_style) && file.exists(RUN$fig_style)) {
  source(RUN$fig_style)
} else {
  ## Without the pipeline's style file, a plain theme and default colours ##
  theme_cosmx <- function(base = 11) theme_classic(base_size = base)
  read_label_style <- function(file) {
    data.table(label = character(), colour = character(), group = character(), display = character())
  }
  display_name <- function(x, style) gsub("_", " ", as.character(x))
  label_levels <- function(x, style) sort(unique(as.character(x[!is.na(x)])))
  label_colours <- function(lv, style) stats::setNames(scales::hue_pal(l = 70, c = 60)(length(lv)), lv)
  label_layer <- function(x, style, untyped = character(), immune = character()) rep(2L, length(x))
  SUBJECT_COLOURS <- c("#0072B2", "#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#F0E442")
  NICHE_COLOURS <- scales::hue_pal(l = 85, c = 40)(16)
}
STYLE <- read_label_style(RUN$label_colours)

## Colours the figures share ##
TEAL <- "#4C7A8C"
WARM <- "#E07B7B"
COOL <- "#6FA8DC"
PALE <- "grey82"
KEPT <- "#9CC3DA"
SITE_DPI <- 200

theme_report <- function(base = 11) theme_cosmx(base)
tilt_x <- theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
no_axis_lines <- theme(axis.line = element_blank(), axis.ticks = element_blank())

## The figure helpers, the shared helpers and the pipeline's style, so a change to any draws every figure again ##
FIG_CODE <- md5_text(c(readLines("_figures.R"), readLines("_functions.R"), if (file.exists(RUN$fig_style)) readLines(RUN$fig_style)))

## A figure in results/report_figures, drawn again only when its inputs, code or size change ##
rfig <- function(name, inputs, draw, w = 7, h = 4.5) {
  inputs <- unique(inputs[!is.na(inputs) & nzchar(inputs) & file.exists(inputs)])
  if (!length(inputs)) {
    return(character())
  }
  w <- min(max(w, 3), 16)
  h <- min(max(h, 2), 30)
  stem <- gsub("[^A-Za-z0-9_.-]+", "_", name)
  key <- substr(fig_key(name, inputs, c(FIG_CODE, deparse(draw)), w, h), 1, 12)
  out <- res_path("report_figures", paste0(stem, "_", key, ".png"))
  if (file.exists(out)) {
    return(out)
  }
  p <- draw()
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  unlink(list.files(dirname(out), paste0("^", rx(stem), "_[0-9a-f]{12}\\.png$"), full.names = TRUE))
  if (is.null(p)) {
    return(character())
  }
  ## Drawn into a part file that takes the figure's name only once it is whole ##
  part <- paste0(out, ".part")
  on.exit(unlink(part), add = TRUE)
  ggsave(part, p, width = w, height = h, dpi = SITE_DPI, bg = "white", device = ragg::agg_png, limitsize = FALSE)
  file.rename(part, out)
  out
}

## Figures as cards, captioned by the names of a list whose empty entries are dropped ##
rgallery <- function(items, ..., empty = NULL) {
  items <- items[lengths(items) > 0]
  gallery(unlist(items, use.names = FALSE), names(items), trim = FALSE, empty = empty, scale = 100 / SITE_DPI)
}

## One figure at a time from a named list, its empty entries dropped ##
rpicker <- function(items, empty = "No figures yet.") {
  items <- items[lengths(items) > 0]
  figure_picker(unlist(items, use.names = FALSE), names(items), trim = FALSE, empty = empty, scale = 100 / SITE_DPI)
}

## The figures the pipeline wrote for a page, folded ##
pipeline_figures <- function(files) {
  files <- unique(files[file.exists(files)])
  if (!length(files)) {
    return(NULL)
  }
  tags$details(
    class = "made", tags$summary(sprintf("Figures the pipeline wrote (%d)", length(files))),
    gallery(files, nice(files), wide = TRUE)
  )
}

## ---- Colours ----

## A colour per cluster ##
cluster_colours <- function(lv) stats::setNames(scales::hue_pal(l = 68, c = 55)(length(lv)), lv)

## A colour per sample, Okabe-Ito first ##
sample_colours <- function(lv) {
  more <- max(0L, length(lv) - length(SUBJECT_COLOURS))
  stats::setNames(c(SUBJECT_COLOURS, if (more) scales::hue_pal(l = 65, c = 60)(more))[seq_along(lv)], lv)
}

## Samples in samplesheet order ##
sample_levels <- function(x) {
  ss <- samples()
  s <- if (!is.null(ss)) ss$sample_id else character()
  u <- unique(as.character(x[!is.na(x)]))
  c(intersect(s, u), sort(setdiff(u, s)))
}

## The colour file's colours for labels, grey for those without a cell identity ##
cell_colours <- function(lv) {
  v <- label_colours(lv, STYLE)
  v[lv %in% c(RUN$untyped, "unknown", "none", "Unassigned")] <- "grey78"
  v
}

## A cell-type fill or colour scale with display names ##
scale_cells <- function(x, aesthetics = "fill", ...) {
  lv <- label_levels(x, STYLE)
  scale_discrete_manual(aesthetics, values = cell_colours(lv), breaks = lv, labels = display_name(lv, STYLE), ...)
}

## Niche types in their numeric order, type_1 first ##
niche_levels <- function(x) {
  u <- unique(as.character(x[!is.na(x)]))
  u[order(suppressWarnings(as.numeric(sub("^type_", "", u))), u)]
}

## Niche types in order, coloured as the pipeline's niche maps ##
niche_colours <- function(lv) {
  i <- suppressWarnings(as.integer(sub("^type_", "", lv)))
  v <- rep(NICHE_COLOURS, length.out = max(c(i, length(lv)), na.rm = TRUE))
  stats::setNames(ifelse(is.na(i), "grey80", v[pmax(i, 1L)]), lv)
}

## ---- Data from the app bundle ----

## The bundle's cells, this stage's own copy first, then results/local_export ##
bundle_file <- function() {
  cur <- if (nzchar(RUN$current)) Sys.glob(paste0(RUN$current, "*/cells_umap_clusters.rds")) else character()
  f <- c(cur, res_path("local_export", "cells_umap_clusters.rds"))
  f <- f[file.exists(f)]
  if (length(f)) f[1] else NA_character_
}
bundle_cells <- function() {
  if (is.null(.cache$bundle)) .cache$bundle <- as.data.table(readRDS(bundle_file()))
  .cache$bundle
}
bundle_has <- function(cols) {
  f <- bundle_file()
  if (is.na(f)) {
    return(FALSE)
  }
  if (is.null(.cache$bundle_cols)) .cache$bundle_cols <- names(bundle_cells())
  all(cols %in% .cache$bundle_cols)
}

## ---- Plots ----

## Cells on a UMAP coloured by one column, with a legend of counts or the labels on the plot ##
umap_plot <- function(d, col, values = NULL, labels = NULL, x = "UMAP1", y = "UMAP2",
                      centroids = FALSE, layer = NULL, size = NULL, legend_cols = 1L) {
  v <- as.character(d[[col]])
  v[is.na(v) | !nzchar(v)] <- "none"
  ok <- is.finite(d[[x]]) & is.finite(d[[y]])
  present <- unique(v[ok])
  lv <- if (is.null(values)) cluster_levels(present) else c(intersect(names(values), present), sort(setdiff(present, names(values))))
  if (is.null(values)) values <- cluster_colours(lv)
  miss <- setdiff(lv, names(values))
  values <- c(values, stats::setNames(rep("grey78", length(miss)), miss))
  lab <- if (is.null(labels)) lv else ifelse(is.na(labels[lv]), lv, labels[lv])
  n <- as.integer(table(factor(v[ok], lv)))
  if (is.null(size)) size <- min(1.5, max(0.08, 40 / sqrt(sum(ok))))
  pd <- data.frame(x = d[[x]][ok], y = d[[y]][ok], g = factor(v[ok], lv))
  o <- sample.int(nrow(pd))
  if (!is.null(layer)) o <- o[order(layer[ok][o])]
  pd <- pd[o, ]
  p <- ggplot(pd, aes(x, y, colour = g)) +
    geom_point(size = size, stroke = 0, shape = 16) +
    scale_colour_manual(
      values = values[lv], breaks = lv, labels = paste0(lab, " (", fmt_n(n), ")"),
      guide = if (legend_cols > 0) guide_legend(ncol = legend_cols, override.aes = list(size = 2.6)) else "none"
    ) +
    coord_equal() +
    theme_void(base_size = 11) +
    theme(legend.title = element_blank(), legend.key.spacing.y = unit(1, "pt"), plot.margin = margin(6, 6, 6, 6))
  if (centroids || legend_cols == 0) {
    cen <- stats::aggregate(cbind(x, y) ~ g, data = pd, FUN = stats::median)
    p <- p + if (nrow(cen) > 40) {
      ggrepel::geom_text_repel(
        data = cen, aes(x, y, label = g), inherit.aes = FALSE, size = 2.2, colour = "grey15", max.overlaps = Inf,
        min.segment.length = 0.3, segment.size = 0.15, segment.colour = "grey50", seed = 1
      )
    } else {
      geom_text(data = cen, aes(x, y, label = g), inherit.aes = FALSE, size = 3, colour = "grey15")
    }
  }
  p
}

## A UMAP figure of one bundle column, cached on the bundle ##
umap_fig <- function(name, col, ..., w = 7.5, h = 6, x = "UMAP1", y = "UMAP2") {
  if (!bundle_has(c(col, x, y))) {
    return(character())
  }
  args <- list(...)
  k <- legend_columns(length(unique(bundle_cells()[[col]])))
  rfig(name, bundle_file(), function() {
    do.call(umap_plot, c(list(bundle_cells(), col, x = x, y = y, legend_cols = k), args))
  }, w + 1.7 * max(0, k - 1), h)
}

## A UMAP of cell types in the colour file's colours, untyped cells underneath ##
umap_cells_fig <- function(name, col) {
  if (!bundle_has(c(col, "UMAP1", "UMAP2"))) {
    return(character())
  }
  k <- legend_columns(length(unique(bundle_cells()[[col]])))
  rfig(name, c(bundle_file(), RUN$label_colours), function() {
    d <- bundle_cells()
    lv <- label_levels(d[[col]], STYLE)
    umap_plot(d, col,
      values = cell_colours(lv), labels = stats::setNames(display_name(lv, STYLE), lv),
      layer = label_layer(d[[col]], STYLE, RUN$untyped), legend_cols = k
    )
  }, 8 + 1.7 * max(0, k - 1), 6)
}

## Genes by effect and significance, the strongest of each panel labelled, up to three panels to a row, or NULL without a gene to show ##
volcano_plot <- function(d, lfc, p, cut = 0.05, gene = NULL, facet = NULL, flag = NULL,
                         n_label = 6L, xlab = "log2 fold change", ylab = "-log10 FDR", ncol = 3L) {
  d <- as.data.frame(d)
  d$x <- d[[lfc]]
  d$y <- -log10(pmax(d[[p]], 1e-300))
  d$sig <- factor(
    ifelse(d[[p]] < cut & d$x > 0, "higher", ifelse(d[[p]] < cut & d$x < 0, "lower", "not significant")),
    c("higher", "lower", "not significant")
  )
  d$flagged <- if (!is.null(flag)) d[[flag]] %in% TRUE else FALSE
  d <- d[is.finite(d$x) & is.finite(d$y), ]
  if (!nrow(d)) {
    return(NULL)
  }
  ## Each panel stops a quarter past its largest hit, or at three times its upper quartile, with genes past it at the edge ##
  panel <- if (!is.null(facet)) as.character(d[[facet]]) else rep("all", nrow(d))
  edge <- vapply(split(seq_len(nrow(d)), panel), function(i) {
    max(1.25 * abs(d$x[i][d$sig[i] != "not significant"]), 3 * stats::quantile(abs(d$x[i]), 0.75, names = FALSE), 1)
  }, numeric(1))
  d$x <- pmax(pmin(d$x, edge[panel]), -edge[panel])
  d <- d[order(d$sig != "not significant"), ]
  g <- ggplot(d, aes(x, y)) +
    geom_hline(yintercept = -log10(cut), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(aes(colour = sig, shape = flagged), size = 1, stroke = 0.4) +
    scale_colour_manual(values = c(higher = WARM, lower = COOL, `not significant` = PALE), name = NULL, drop = FALSE) +
    scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 2), name = NULL,
      labels = c(`FALSE` = "clean", `TRUE` = "may be bleed from neighbours"), guide = if (is.null(flag)) "none" else "legend"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
    labs(x = xlab, y = ylab) +
    theme_report()
  if (!is.null(gene) && n_label > 0L) {
    top <- d[d$sig != "not significant", ]
    if (nrow(top)) {
      top <- top[order(top[[p]]), ]
      if (!is.null(facet)) top <- do.call(rbind, lapply(split(top, top[[facet]]), head, n_label)) else top <- head(top, n_label)
      g <- g + ggrepel::geom_text_repel(
        data = top, aes(label = .data[[gene]]), size = 2.6, max.overlaps = 30,
        min.segment.length = 0, segment.colour = "grey60", segment.size = 0.2
      )
    }
  }
  if (!is.null(facet)) g <- g + facet_wrap(facet, scales = "free", ncol = min(ncol, length(unique(d[[facet]]))))
  g
}

## A tiled heatmap, diverging around 0 or sequential from white, with dots marking a column ##
heat_plot <- function(d, x, y, fill, diverging = FALSE, name = fill, dot = NULL, limits = NULL) {
  d <- as.data.frame(d)
  g <- ggplot(d, aes(.data[[x]], .data[[y]], fill = .data[[fill]])) +
    geom_tile(colour = "white", linewidth = 0.3)
  g <- g + if (diverging) {
    lim <- limits %||% (c(-1, 1) * max(abs(d[[fill]]), na.rm = TRUE))
    scale_fill_gradient2(low = COOL, mid = "white", high = WARM, midpoint = 0, name = name, limits = lim, oob = scales::squish)
  } else {
    scale_fill_gradient(low = "white", high = TEAL, name = name, limits = limits, na.value = "grey95")
  }
  if (!is.null(dot)) {
    g <- g + geom_point(data = d[d[[dot]] %in% TRUE, , drop = FALSE], shape = 16, size = 0.9, colour = "grey15", show.legend = FALSE)
  }
  g + labs(x = NULL, y = NULL) + theme_report() + tilt_x + no_axis_lines
}

## Width and height in inches for a heatmap of nx columns and ny rows ##
heat_size <- function(nx, ny, extra_w = 2.6, extra_h = 1.6) {
  c(w = extra_w + 0.24 * nx, h = extra_h + 0.2 * ny)
}

## Size of a cluster by panel heatmap, the clusters across while they fit 16 inches and down the side past that ##
cluster_heat_size <- function(n_clusters, n_panels) {
  across <- heat_size(n_clusters, n_panels)
  if (across[["w"]] <= 16) c(across, down = 0) else c(heat_size(n_panels, n_clusters), down = 1)
}
