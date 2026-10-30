## The pipeline's figure style: one pinned colour and a display name per label
## (assets/label_colours.csv), immune cells drawn last, low-confidence labels striped, tissue
## maps drawn from the cell outlines with the legend on the right and its counts, 600 dpi PNGs
## with a vector PDF beside every figure that is not a map, no plot titles except the sample on
## a map. Sourced by the bin/ scripts that draw cell types, niches or per-subject changes ##
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## Resolution of every PNG (dots per inch), and the base text size (pt) ##
FIG_DPI <- 600
FIG_PT <- 11
## Low-confidence stripes on maps: spacing (um) and line width ##
STRIPE_UM <- 3
STRIPE_LWD <- 0.12
## Niche types in the presentation's colours, type_1 first; later types take the rest ##
NICHE_COLOURS <- c(
  "#C3E7FF", "#F7DD9E", "#95C7C4", "#E0D2F2", "#A5FDDB", "#C4B3E1", "#D5C6B9", "#B5DDD8",
  "#D499ED", "#80FDCC", "#FFD1DC", "#D8E2A8", "#F5C6A5", "#B8C9F2", "#E8E3B0", "#C9E4CA"
)
## Subjects seen at more than one level, in Okabe-Ito order. A subject seen once is grey ##
SUBJECT_COLOURS <- c("#0072B2", "#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#F0E442")
SINGLE_COLOUR <- "grey55"
## Labels the colour file does not list take these, in order ##
EXTRA_COLOURS <- c(
  "#85660D", "#F99379", "#FE00FA", "#008856", "#90AD1C", "#3283FE", "#882D17", "#D55E00",
  "#E68FAC", "#3B00FB", "#E41A1C", "#325A9B", "#B10DA1", "#1CBE4F", "#AA0DFE", "#F28E2B"
)

## The label colour file: label, colour, group (immune, non_immune or untyped) and display, the
## name figures show. A label without a display name shows with spaces for underscores ##
read_label_style <- function(file) {
  s <- if (nzchar(file) && file.exists(file)) {
    fread(file, colClasses = "character")
  } else {
    data.table(label = character(), colour = character())
  }
  if (!"group" %in% names(s)) s[, group := "non_immune"]
  if (!"display" %in% names(s)) s[, display := NA_character_]
  s[is.na(display) | !nzchar(display), display := gsub("_", " ", label)]
  s[]
}

display_name <- function(x, style) {
  x <- as.character(x)
  d <- style$display[match(x, style$label)]
  ifelse(is.na(d), gsub("_", " ", x), d)
}

## Display names inside a text such as "PT 71%, T_NK 32%" ##
display_in_text <- function(x, style) {
  for (k in style$label[order(-nchar(style$label))]) {
    x <- gsub(paste0("\\b", k, "\\b"), style$display[style$label == k][1], x)
  }
  gsub("_", " ", x)
}

## A niche type's colour by its number: NICHE_COLOURS, then pale hues past their end ##
niche_type_colour <- function(type) {
  i <- as.integer(sub("^type_", "", type))
  more <- max(0L, max(i, na.rm = TRUE) - length(NICHE_COLOURS))
  c(NICHE_COLOURS, if (more) scales::hue_pal(l = 85, c = 40)(more))[i]
}

## A colour per label: the file's, then EXTRA_COLOURS for labels it does not list ##
label_colours <- function(lv, style) {
  col <- style$colour[match(lv, style$label)]
  new <- is.na(col)
  col[new] <- rep_len(setdiff(EXTRA_COLOURS, style$colour), sum(new))
  stats::setNames(col, lv)
}

## Drawing order: untyped labels underneath, immune labels on top. untyped and immune add to the
## groups in the colour file ##
label_layer <- function(x, style, untyped = character(), immune = character()) {
  x <- as.character(x)
  g <- style$group[match(x, style$label)]
  fifelse(
    x %in% untyped | g %in% "untyped",
    1L,
    fifelse(x %in% immune | g %in% "immune", 3L, 2L)
  )
}

## Labels in legend order: the colour file's order, then any others alphabetically ##
label_levels <- function(x, style) {
  x <- unique(as.character(x[!is.na(x)]))
  c(intersect(style$label, x), sort(setdiff(x, style$label)))
}

theme_cosmx <- function(base = FIG_PT) {
  theme_classic(base_size = base) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      legend.key.spacing.y = unit(1, "pt")
    )
}

## Formats every figure is saved in, set from the scripts' --fig_format: png, pdf and jpeg,
## comma separated, or all ##
FIG_FORMATS <- "png"
set_fig_format <- function(x) {
  f <- tolower(trimws(strsplit(x, ",")[[1]]))
  f <- f[nzchar(f)]
  if (!length(f)) f <- "png"
  if ("all" %in% f) f <- c("png", "pdf", "jpeg")
  bad <- setdiff(f, c("png", "pdf", "jpeg"))
  if (length(bad)) {
    stop("fig_format takes png, pdf, jpeg or all, not ", paste(bad, collapse = ", "))
  }
  FIG_FORMATS <<- unique(f)
}

## A figure in every format of FIG_FORMATS, beside file (named .png). Sizes in inches. A map is
## never a PDF, one path per cell makes it huge; asked only for PDFs, it comes as a PNG ##
save_fig <- function(p, file, w, h, map = FALSE) {
  formats <- FIG_FORMATS
  if (map) formats <- setdiff(formats, "pdf")
  if (!length(formats)) formats <- "png"
  stem <- sub("\\.png$", "", file)
  for (f in formats) {
    switch(
      f,
      png = ggsave(
        paste0(stem, ".png"), p, width = w, height = h, dpi = FIG_DPI, bg = "white",
        device = ragg::agg_png, limitsize = FALSE
      ),
      jpeg = ggsave(
        paste0(stem, ".jpg"), p, width = w, height = h, dpi = FIG_DPI, bg = "white",
        device = ragg::agg_jpeg, limitsize = FALSE
      ),
      pdf = ggsave(
        paste0(stem, ".pdf"), p, width = w, height = h, bg = "white",
        device = grDevices::pdf, limitsize = FALSE
      )
    )
  }
  invisible(file)
}

## Radius (mm) of a ggplot point of this size and stroke, as grid draws shape 16 ##
dot_r <- function(size, stroke = 0.5) 0.375 * (size * .pt + stroke * .stroke / 2) / .pt

## Legend entry for low-confidence labels: a grey dot with white stripes ##
draw_key_low <- function(data, params, size) {
  tile <- grid::grobTree(
    grid::rectGrob(gp = grid::gpar(fill = "grey55", col = NA)),
    grid::segmentsGrob(
      x0 = c(0, -0.5, 0.5), y0 = c(0, 0.5, -0.5), x1 = c(1, 0.5, 1.5), y1 = c(1, 1.5, 0.5),
      gp = grid::gpar(col = "white", lwd = 1.5, lineend = "square")
    )
  )
  stripes <- grid::pattern(
    tile,
    width = unit(1.2, "mm"),
    height = unit(1.2, "mm"),
    extend = "repeat"
  )
  grid::circleGrob(r = unit(dot_r(4, 0), "mm"), gp = grid::gpar(fill = stripes, col = NA))
}

## The legend's low-confidence entry, as an invisible point carrying draw_key_low. pos: one
## row with x and y anywhere on the plot ##
low_key_layers <- function(pos, n_low) {
  key <- data.table(x = pos$x[1], y = pos$y[1], conf = sprintf(
    "Low confidence (%s)",
    prettyNum(n_low, big.mark = ",")
  ))
  list(
    geom_point(
      data = key, aes(x, y, shape = conf), inherit.aes = FALSE, size = 0,
      key_glyph = draw_key_low
    ),
    scale_shape_manual(values = 16, name = NULL),
    guides(shape = guide_legend(order = 2))
  )
}

## Cell outlines of a Giotto object: cell_ID, part (one ring) and x, y in pixels. Holes are left
## out; a cell drawn in several parts keeps them all ##
cell_outlines <- function(g, cell_ids = NULL) {
  sv <- Giotto::getPolygonInfo(g, polygon_name = "cell", return_giottoPolygon = TRUE)[]
  crd <- terra::geom(sv)
  d <- data.table(
    cell_ID = as.character(sv$poly_ID)[crd[, "geom"]],
    part = paste(crd[, "geom"], crd[, "part"]),
    x = crd[, "x"],
    y = crd[, "y"]
  )[crd[, "hole"] == 0]
  if (!is.null(cell_ids)) d <- d[cell_ID %in% cell_ids]
  d[]
}

## White diagonal lines every spacing, clipped to the union of some outlines (part, x, y) ##
stripe_paths <- function(v, spacing = STRIPE_UM) {
  empty <- data.table(X = numeric(), Y = numeric(), g = character())
  if (!nrow(v) || !requireNamespace("sf", quietly = TRUE)) return(empty)
  rings <- lapply(split(v, by = "part"), function(p) {
    xy <- as.matrix(p[, .(x, y)])
    if (nrow(unique(xy)) < 3) return(NULL)
    sf::st_polygon(list(rbind(xy, xy[1, ])))
  })
  rings <- Filter(Negate(is.null), rings)
  if (!length(rings)) return(empty)
  shape <- sf::st_union(sf::st_make_valid(sf::st_sfc(rings)))
  b <- sf::st_bbox(shape)
  cs <- seq(b[["ymin"]] - b[["xmax"]], b[["ymax"]] - b[["xmin"]], by = spacing * sqrt(2))
  lines <- sf::st_sfc(lapply(cs, function(k) {
    sf::st_linestring(rbind(c(b[["xmin"]], b[["xmin"]] + k), c(b[["xmax"]], b[["xmax"]] + k)))
  }))
  cut <- sf::st_intersection(lines, shape)
  cut <- cut[sf::st_geometry_type(cut) %in% c("LINESTRING", "MULTILINESTRING")]
  if (!length(cut)) return(empty)
  xy <- as.data.table(sf::st_coordinates(sf::st_cast(sf::st_sfc(cut), "MULTILINESTRING")))
  xy[, g := paste(L1, L2)][, .(X, Y, g)]
}

## A scale bar about a fifth of the width, in a round number of um ##
scale_bar_um <- function(width_um) {
  nice <- c(50, 100, 250, 500, 1000, 2000, 5000)
  nice[max(1L, findInterval(width_um / 5, nice))]
}

## One tissue map, as the presentation draws them: every cell as its outline, filled by group;
## groups in on_top drawn last, with an edge when top_line is set; low-confidence cells striped;
## a scale bar; the legend on the right in one column. d: cell_ID, group (a factor), low and
## optionally fov, whose cells are striped together. Returns the plot and its size (inches) ##
biopsy_map <- function(d, outlines, values, labels, um_per_px, legend_title, title = NULL,
                       on_top = character(), top_line = NA, top_lwd = 0, legend_w = 3.5) {
  d <- copy(d)
  if (!"fov" %in% names(d)) d[, fov := "all"]
  v <- outlines[d, on = "cell_ID", nomatch = NULL]
  v[, `:=`(x = x * um_per_px, y = y * um_per_px)]
  st <- rbindlist(lapply(split(v[low == TRUE], by = "fov"), function(w) {
    stripe_paths(w)[, g := paste(w$fov[1], g)]
  }))
  if (!nrow(st)) st <- data.table(X = numeric(), Y = numeric(), g = character())
  wx <- diff(range(v$x))
  aspect <- diff(range(v$y)) / max(wx, 1)
  map_w <- if (aspect < 0.5) 16 else 9
  bar <- scale_bar_um(wx)
  y0 <- min(v$y) - 0.04 * diff(range(v$y))
  p <- ggplot(v, aes(x, y, group = part, fill = group)) +
    geom_polygon(data = v[!group %in% on_top], colour = NA, key_glyph = "point") +
    geom_path(
      data = st, aes(X, Y, group = g), inherit.aes = FALSE,
      colour = "white", linewidth = STRIPE_LWD
    ) +
    geom_polygon(
      data = v[group %in% on_top], colour = top_line, linewidth = top_lwd,
      key_glyph = "point"
    ) +
    annotate("segment", x = max(v$x) - bar, xend = max(v$x), y = y0, yend = y0, linewidth = 0.8) +
    annotate(
      "text", x = max(v$x) - bar / 2, y = y0, label = paste(bar, "µm"),
      vjust = -0.6, size = 0.9 * FIG_PT / .pt
    ) +
    scale_fill_manual(
      values = values, labels = labels, breaks = names(values), limits = names(values),
      name = legend_title
    ) +
    guides(fill = guide_legend(
      override.aes = list(shape = 21, size = 4, colour = NA), ncol = 1, order = 1
    )) +
    coord_fixed(clip = "off") +
    labs(title = title) +
    theme_void(base_size = FIG_PT + 2) +
    theme(
      plot.title = element_text(face = "bold", margin = margin(b = 8)),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      legend.key.spacing.y = unit(2, "pt"),
      plot.margin = margin(10, 12, 10, 12)
    )
  n_low <- sum(d$low)
  if (n_low > 0) p <- p + low_key_layers(v[1, .(x, y)], n_low)
  n_keys <- length(values) + (n_low > 0)
  list(plot = p, w = map_w + legend_w, h = max(map_w * aspect + 1.2, 0.3 * n_keys + 1))
}

## A white line across a dot at 45 degrees, sized in mm so it stays inside the dot ##
GeomSlash <- ggproto("GeomSlash", Geom,
  required_aes = c("x", "y"),
  default_aes = aes(colour = "white", linewidth = 0.13),
  draw_key = draw_key_blank,
  draw_panel = function(data, panel_params, coord, half = 0.2) {
    d <- coord$transform(data, panel_params)
    h <- unit(half, "mm")
    grid::segmentsGrob(
      x0 = unit(d$x, "npc") - h, y0 = unit(d$y, "npc") - h,
      x1 = unit(d$x, "npc") + h, y1 = unit(d$y, "npc") + h,
      gp = grid::gpar(col = d$colour, lwd = d$linewidth * .pt, lineend = "butt")
    )
  }
)

## Cells on the UMAP, coloured by group, low-confidence labels drawn first
## with a white slash across each dot, then the rest by layer (untyped first, immune last); axes
## as two arrows. d: UMAP1, UMAP2, group (a factor), low, layer ##
umap_plot <- function(d, values, labels, size = 0.4, seed = 1) {
  set.seed(seed)
  d <- d[order(!low, layer, sample(.N))]
  half <- 0.85 * dot_r(size) / sqrt(2)
  slashes <- layer(
    geom = GeomSlash, stat = "identity", position = "identity", data = d[low == TRUE],
    mapping = aes(UMAP1, UMAP2), inherit.aes = FALSE, show.legend = FALSE,
    params = list(half = half, linewidth = 0.13, na.rm = FALSE)
  )
  ax <- d[, .(x = min(UMAP1) - 1, y = min(UMAP2) - 1, len = 0.1 * diff(range(UMAP1)))]
  tip <- arrow(length = unit(2.5, "mm"), type = "closed")
  p <- ggplot(d, aes(UMAP1, UMAP2, colour = group)) +
    geom_point(data = d[low == TRUE], size = size, shape = 16) +
    slashes +
    geom_point(data = d[low == FALSE], size = size, shape = 16) +
    annotate(
      "segment", x = ax$x, xend = ax$x + ax$len, y = ax$y, yend = ax$y,
      arrow = tip, linewidth = 0.5
    ) +
    annotate(
      "segment", x = ax$x, xend = ax$x, y = ax$y, yend = ax$y + ax$len,
      arrow = tip, linewidth = 0.5
    ) +
    annotate(
      "text", x = ax$x + ax$len / 2, y = ax$y, label = "UMAP 1",
      vjust = 1.6, size = 0.9 * FIG_PT / .pt
    ) +
    annotate(
      "text", x = ax$x, y = ax$y + ax$len / 2, label = "UMAP 2",
      angle = 90, vjust = -0.8, size = 0.9 * FIG_PT / .pt
    ) +
    scale_colour_manual(
      values = values, labels = labels, breaks = names(values), limits = names(values),
      name = NULL
    ) +
    guides(colour = guide_legend(
      ncol = 1, override.aes = list(size = 4, stroke = 0), order = 1
    )) +
    coord_fixed(clip = "off") +
    theme_void(base_size = FIG_PT + 2) +
    theme(legend.key.spacing.y = unit(1, "pt"), plot.margin = margin(10, 12, 10, 12))
  n_low <- sum(d$low)
  if (n_low > 0) p <- p + low_key_layers(d[1, .(x = UMAP1, y = UMAP2)], n_low)
  p
}

## Fold change from the first level to the last for each subject seen at both, as "32x lower" ##
fold_labels <- function(d) {
  f <- d[order(level), if (.N >= 2 && value[1] > 0 && value[.N] > 0) {
    .(level = level[.N], value = value[.N], r = value[.N] / value[1])
  }, by = .(panel, subject)]
  if (!nrow(f)) return(data.table(panel = d$panel[0], subject = d$subject[0], level = d$level[0],
                                  value = numeric(), lab = character()))
  f[, lab := sprintf(
    "%s× %s",
    formatC(pmax(r, 1 / r), digits = 2, format = "fg"),
    fifelse(r < 1, "lower", "higher")
  )]
  f[]
}

## One line per subject across the levels (e.g. timepoints), one panel per measure.
## Subjects seen at one level only are grey, with a dotted line at their value
## through every level. A log axis when every value is above zero; fold labels at the last
## level. d: panel, subject, level (a factor, baseline first), value ##
subject_lines <- function(d, ylab, log = TRUE, fold = TRUE, subject_names = NULL) {
  d <- copy(d)
  seen <- d[, .(n = uniqueN(level)), by = subject]
  multi <- seen[n > 1, subject]
  subjects <- c(sort(as.character(multi)), sort(setdiff(as.character(seen$subject), multi)))
  cols <- stats::setNames(rep(SINGLE_COLOUR, length(subjects)), subjects)
  cols[subjects %in% multi] <- rep_len(SUBJECT_COLOURS, sum(subjects %in% multi))
  d[, subject := factor(subject, levels = subjects)]
  log <- log && all(d$value > 0, na.rm = TRUE)
  single <- d[!subject %in% multi]
  leg_lines <- ifelse(subjects %in% multi, "solid", "dotted")
  p <- ggplot(d, aes(level, value, colour = subject, group = subject)) +
    geom_hline(
      data = single, aes(yintercept = value, colour = subject),
      linetype = "dotted", linewidth = 0.8, show.legend = FALSE
    ) +
    # Every subject, so the legend draws a line for each; one point draws no line
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.6) +
    scale_colour_manual(
      values = cols, breaks = subjects, name = NULL,
      labels = if (is.null(subject_names)) subjects else subject_names[subjects]
    ) +
    scale_x_discrete(expand = expansion(add = c(0.3, 0.6))) +
    labs(x = NULL, y = ylab) +
    guides(colour = guide_legend(ncol = 1, override.aes = list(linetype = leg_lines))) +
    theme_cosmx() +
    theme(legend.key.width = unit(2.4, "lines"))
  if (uniqueN(d$panel) > 1) p <- p + facet_wrap(vars(panel), scales = "free_y", nrow = 1)
  p <- if (log) {
    p + scale_y_log10(
      breaks = scales::breaks_log(6),
      labels = scales::label_number(drop0trailing = TRUE)
    )
  } else {
    p + expand_limits(y = 0)
  }
  if (fold) {
    fd <- fold_labels(d[subject %in% multi])
    if (nrow(fd)) {
      p <- p + geom_text(
        data = fd, aes(label = lab), hjust = 0, nudge_x = 0.12,
        size = 0.8 * FIG_PT / .pt, show.legend = FALSE
      ) + coord_cartesian(clip = "off")
    }
  }
  p
}
