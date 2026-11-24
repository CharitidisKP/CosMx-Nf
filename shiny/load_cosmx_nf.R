# Reads a run's local_export bundle and published tables for the explorer

load_cosmx_nf <- function(root) {
  ex    <- file.path(root, "local_export")
  cells <- as.data.frame(readRDS(file.path(ex, "cells_umap_clusters.rds")))

  # Makes a factor with levels in number order when all parse as numbers
  as_ordered_factor <- function(v) {
    lv <- unique(v[!is.na(v)])
    n  <- suppressWarnings(as.numeric(lv))
    factor(v, levels = if (!anyNA(n)) lv[order(n)] else sort(lv))
  }
  # Text columns and cluster ids, whatever their stored type, become factors
  CATRX <- "^leiden_clus_res|^sub_res|^insitutype|^ht_call$|^cell_type|^monaco_|^fov$|^list_ID$"
  cat_nm <- union(names(cells)[vapply(cells, is.character, logical(1))],
                  grep(CATRX, names(cells), value = TRUE))
  cat_nm <- setdiff(cat_nm, c("cell_ID", grep("_score$|_UMAP[0-9]$", cat_nm, value = TRUE)))
  # Embedding coordinates are numbers whatever the export stored them as
  for (nm in grep("_UMAP[0-9]$", names(cells), value = TRUE))
    cells[[nm]] <- suppressWarnings(as.numeric(as.character(cells[[nm]])))
  for (nm in cat_nm) cells[[nm]] <- as_ordered_factor(as.character(cells[[nm]]))

  # Drops typing diagnostics: _anchor, _second_type, _second_score and ht_class
  drop <- grep("_anchor$|_second_type$|_second_score$|^ht_class$", names(cells), value = TRUE)
  cells <- cells[, setdiff(names(cells), drop), drop = FALSE]

  const <- names(Filter(function(z) is.factor(z) && nlevels(droplevels(z)) < 2, cells))
  cells <- cells[, setdiff(names(cells), const), drop = FALSE]

  # Orders samples by treatment, then timepoint, then id
  samp <- intersect(c("sample_id", "sample"), names(cells))[1]
  if (!is.na(samp) && "treatment" %in% names(cells)) {
    k <- unique(cells[, c(samp, "treatment",
                          intersect("timepoint", names(cells))), drop = FALSE])
    tr <- tolower(as.character(k$treatment))
    ord <- order(tr,
                 if (!is.null(k$timepoint)) as.character(k$timepoint) else character(nrow(k)),
                 as.character(k[[samp]]))
    cells[[samp]] <- factor(cells[[samp]], levels = as.character(k[[samp]][ord]))
  }

  # Moves each sample to the origin, and divides x_rel and y_rel by its width
  if (!is.na(samp) && all(c("x", "y") %in% names(cells))) {
    g  <- cells[[samp]]
    x0 <- ave(cells$x, g, FUN = min); y0 <- ave(cells$y, g, FUN = min)
    sx <- ave(cells$x, g, FUN = function(v) max(diff(range(v)), 1))
    cells$x_rel <- (cells$x - x0) / sx     # fit-width: every tissue drawn to the same width
    cells$y_rel <- (cells$y - y0) / sx
    cells$x_abs <- cells$x - x0            # true scale: translated only, one um-per-pixel
    cells$y_abs <- cells$y - y0
    attr(cells, "rel") <- list(x0 = x0, y0 = y0, sx = sx)   # polygons need the same transform
  }

  # Adds each cluster_labels.<basis>.csv as a label_ column, joined on cluster_id
  for (f in list.files(file.path(root, "tables", "06_annotation", "labels"),
                       "^cluster_labels\\..+\\.csv$", full.names = TRUE)) {
    basis <- sub("^cluster_labels\\.(.+)\\.csv$", "\\1", basename(f))
    if (!basis %in% names(cells)) next
    # colClasses: "3.1" and "3.10" both parse to the number 3.1 and the join collides
    m  <- utils::read.csv(f, colClasses = "character")
    nm <- paste0("label_", sub("^leiden_clus_", "", basis))
    cells[[nm]] <- factor(m$label[match(as.character(cells[[basis]]),
                                        as.character(m$cluster_id))])
  }

  # Adds the subcluster columns from the sub_*_cells.csv tables, not the bundle
  for (f in list.files(file.path(root, "tables", "05_subclustering"), "^sub_.+_cells\\.csv$",
                       full.names = TRUE)) {
    # As text: read as numbers, "3.10" becomes 3.1 and two subclusters merge
    d  <- utils::read.csv(f, colClasses = "character")
    for (uc in grep("_UMAP[12]$", names(d), value = TRUE)) d[[uc]] <- as.numeric(d[[uc]])
    nm <- sub("^(sub_.+)_cells\\.csv$", "\\1", basename(f))
    if (!all(c("cell_ID", nm) %in% names(d))) next
    v <- as.character(d[[nm]])[match(cells$cell_ID, d$cell_ID)]
    # "3.10" is not the number 3.1: sort on parent and child separately or they collide
    u <- unique(v[!is.na(v)])
    ps <- strsplit(u, ".", fixed = TRUE); dp <- max(1L, lengths(ps))
    k <- vapply(ps, function(z) suppressWarnings(as.numeric(c(z, rep(NA_character_, dp))[seq_len(dp)])),
                numeric(dp))
    k <- matrix(k, nrow = dp)
    cells[[nm]] <- factor(v, levels = u[do.call(order, c(lapply(seq_len(dp), function(j) k[j, ]),
                                                          na.last = FALSE))])
    # Adds the per-parent UMAPs from the same CSV, for the Subclusters tab
    for (uc in grep(paste0("^", nm, "_c.+_UMAP[12]$"), names(d), value = TRUE))
      if (!uc %in% names(cells)) cells[[uc]] <- d[[uc]][match(cells$cell_ID, d$cell_ID)]
  }

  # Adds per-cell scores from the HieraType and supervised typing tables
  score_src <- c(
    list(list(f = file.path(root, "tables/04_reference_typing/hieratype/hieratype_calls.csv"),
              col = "ht_score", nm = "ht")),
    lapply(list.files(file.path(root, "tables/04_reference_typing/supervised"), "_cells\\.csv$",
                      full.names = TRUE),
           function(f) list(f = f, col = "score",
                            nm = sub("_cells\\.csv$", "", basename(f)))))
  for (z in score_src) {
    if (!file.exists(z$f)) next
    d <- utils::read.csv(z$f)
    if (!all(c("cell_ID", z$col) %in% names(d))) next
    v <- d[[z$col]][match(cells$cell_ID, d$cell_ID)]
    cells[[paste0(z$nm, "_score")]] <- v
  }

  # Which methods call each cell B-lineage, one logical column per method
  LINEAGE <- list(
    B = list(
      HT   = list(col = "ht_call",                     val = c("bcell", "plasma")),
      KPMP = list(col = "insitutype_sup_KPMP",         val = c("B", "PL")),
      HCA  = list(col = "insitutype_sup_HCA_Kidney",   val = "B_cell"),
      io   = list(col = "insitutype_sup_HCA_io",       val = c("B-cell", "plasmablast")),
      RCC  = list(col = "insitutype_sup_KidneyRCC_6k",
                  val = c("B_cell", "Plasma", "Plasmablast"))
    )
  )
  lineage <- list()
  for (ln in names(LINEAGE)) {
    d <- Filter(function(z) z$col %in% names(cells), LINEAGE[[ln]])
    if (!length(d)) next
    hit <- vapply(d, function(z) as.character(cells[[z$col]]) %in% z$val,
                  logical(nrow(cells)))
    dim(hit) <- c(nrow(cells), length(d)); colnames(hit) <- names(d)
    lineage[[ln]] <- hit
  }

  # Normalised and raw expression, each NULL when the bundle lacks it
  read_expr <- function(fn) {
    f <- file.path(ex, fn)
    if (!file.exists(f)) return(NULL)
    # Attaches Matrix, for `[` on a dgCMatrix to dispatch under Rscript
    library(Matrix)
    m <- readRDS(f)
    m[, match(cells$cell_ID, colnames(m)), drop = FALSE]
  }
  expr     <- read_expr("marker_expr_normalized.rds")
  expr_raw <- read_expr("marker_expr_raw.rds")

  # Adds numeric metadata to both matrices as rows prefixed conf:, IF: or meta:
  NUM <- setdiff(names(cells)[vapply(cells, is.numeric, logical(1))],
                 c("x", "y", "UMAP1", "UMAP2", "UMAPpre1", "UMAPpre2", "PCA_1", "PCA_2",
                   "x_rel", "y_rel", "x_abs", "y_abs", "total_counts"))
  if (length(NUM)) {
    rows <- t(as.matrix(cells[, NUM, drop = FALSE]))
    rownames(rows) <- ifelse(grepl("_score$", NUM), paste0("conf:", sub("_score$", "", NUM)),
                      ifelse(grepl("^(Mean|Max)\\.", NUM), paste0("IF:", NUM),
                             paste0("meta:", NUM)))
    colnames(rows) <- cells$cell_ID
    rows[!is.finite(rows)] <- 0
    add_if <- function(m) if (is.null(m)) NULL else rbind(m, rows[, colnames(m), drop = FALSE])
    expr     <- add_if(expr)
    expr_raw <- add_if(expr_raw)
  }

  # Cell polygons, where i is a row of the cells file from the same bundle
  f <- file.path(ex, "cell_polygons.rds")
  poly <- if (file.exists(f)) as.data.frame(readRDS(f)) else NULL

  # Each step's per-cluster diagnostics and top DE genes, keyed by cluster column
  dd <- file.path(root, "tables", c("03_clustering", "04_reference_typing", "05_subclustering"),
                  "markers")
  diag <- de <- list()
  for (f in list.files(dd, "^cluster_labels\\..+\\.template\\.csv$", full.names = TRUE)) {
    col <- sub("^cluster_labels\\.(.+)\\.template\\.csv$", "\\1", basename(f))
    t <- utils::read.csv(f, colClasses = c(cluster_id = "character"))
    diag[[col]] <- t[order(-t$n_cells), ]
    g <- file.path(dirname(f), paste0(col, "_de_top.csv"))
    if (file.exists(g)) {
      x <- utils::read.csv(g, colClasses = c(cluster = "character"))
      de[[col]] <- x
    }
  }

  list(cells = cells, poly = poly, expr = expr, expr_raw = expr_raw,
       diag = diag, de = de, lineage = lineage)
}
