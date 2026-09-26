# Read a CosMx-Nf results tree into the explorer's contract: the local_export bundle that
# shiny/export_app_bundle.sh writes, plus the published tables. Rsync both down; the
# Giotto objects stay on the server.

load_cosmx_nf <- function(root) {
  ex    <- file.path(root, "local_export")
  cells <- as.data.frame(readRDS(file.path(ex, "cells_umap_clusters.rds")))

  # Leiden resolutions come through numeric. They are categories, not quantities —
  # left numeric they would get a continuous scale and sort as 1,10,11,2.
  # Cluster ids arrive as character, so plain sort() gives 1, 10, 11, 2. Order numerically
  # whenever every level parses as a number; fall back to lexicographic for real labels.
  as_ordered_factor <- function(v) {
    lv <- unique(v[!is.na(v)])
    n  <- suppressWarnings(as.numeric(lv))
    factor(v, levels = if (!anyNA(n)) lv[order(n)] else sort(lv))
  }
  # IF / morphology channels are quantities, not categories. Held out of the factor sweep and
  # served through the gene machinery below, so they get the continuous scale and the
  # threshold filter for free.
  # Characters are categories; so are cluster ids, whatever type they were stored as. The
  # loader decides this rather than trusting the export to have typed them.
  CATRX <- "^leiden_clus_res|^sub_res|^insitutype|^ht_call$|^cell_type|^monaco_|^fov$|^list_ID$"
  cat_nm <- union(names(cells)[vapply(cells, is.character, logical(1))],
                  grep(CATRX, names(cells), value = TRUE))
  cat_nm <- setdiff(cat_nm, c("cell_ID", grep("_score$|_UMAP[0-9]$", cat_nm, value = TRUE)))
  # Embedding coordinates are numbers whatever the export stored them as
  for (nm in grep("_UMAP[0-9]$", names(cells), value = TRUE))
    cells[[nm]] <- suppressWarnings(as.numeric(as.character(cells[[nm]])))
  for (nm in cat_nm) cells[[nm]] <- as_ordered_factor(as.character(cells[[nm]]))

  # One colour-by column per typing method: the call. `_anchor` flags the exemplars used to
  # rescale the reference, `_second_type` is the runner-up, `ht_class` the coarse tier --
  # all fitting diagnostics that read as peer annotations in a dropdown.
  drop <- grep("_anchor$|_second_type$|_second_score$|^ht_class$", names(cells), value = TRUE)
  cells <- cells[, setdiff(names(cells), drop), drop = FALSE]

  const <- names(Filter(function(z) is.factor(z) && nlevels(droplevels(z)) < 2, cells))
  cells <- cells[, setdiff(names(cells), const), drop = FALSE]

  # Sample order: treatment, then timepoint, then id, so each subject's samples sit together
  samp <- intersect(c("sample_id", "sample"), names(cells))[1]
  if (!is.na(samp) && "treatment" %in% names(cells)) {
    k <- unique(cells[, c(samp, "treatment",
                          intersect("timepoint", names(cells))), drop = FALSE])
    tr <- tolower(as.character(k$treatment))
    ord <- order(tr,
                 if (!is.null(k$timepoint)) as.character(k$timepoint) else "",
                 as.character(k[[samp]]))
    cells[[samp]] <- factor(cells[[samp]], levels = as.character(k[[samp]][ord]))
  }

  # merge.R shifts each sample along x, so faceted tissue panels would each show one tissue
  # in a sliver of empty space. Re-centring per sample lets the panels share a scale, which
  # is what coord_fixed() needs -- free scales and a fixed aspect ratio are incompatible.
  # Translating alone is not enough: the tissues differ ~4x in extent, so with the shared
  # scale that coord_fixed() requires, the small ones become slivers in an empty panel.
  # Divide both axes by the SAME per-sample factor -- each tissue then fills the panel width
  # while keeping its true shape. Units become arbitrary, so the app hides these axes.
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

  # Manual annotation: one cluster_labels.<basis>.csv per cluster basis, joined on
  # cluster_id. Each becomes its own colour-by column, so auto vs manual sit
  # side by side and you can switch basis without losing the other.
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

  # Read from the results tree, not the bundle: a new subcluster run then needs only its
  # own small CSV rsynced.
  for (f in list.files(file.path(root, "tables", "05_subclustering"), "^sub_.+_cells\\.csv$",
                       full.names = TRUE)) {
    d  <- utils::read.csv(f)
    nm <- sub("^(sub_.+)_cells\\.csv$", "\\1", basename(f))
    if (!all(c("cell_ID", nm) %in% names(d))) next
    v <- as.character(d[[nm]])[match(cells$cell_ID, d$cell_ID)]
    # "3.10" is not the number 3.1: sort on parent and child separately or they collide
    u <- unique(v[!is.na(v)])
    k <- vapply(strsplit(u, ".", fixed = TRUE),
                function(z) suppressWarnings(as.numeric(c(z, NA_character_))[1:2]), numeric(2))
    cells[[nm]] <- factor(v, levels = u[order(k[1, ], k[2, ], na.last = FALSE)])
  }

  # Per-cell confidence. The bundle deliberately carries no numeric score columns, but each
  # stage publishes them per cell and those CSVs come down with the rest of the tree. Read
  # them here rather than widening the export: no server round trip, and the join is on
  # cell_ID, which cannot silently misalign the way a row index would.
  # Banded, not left continuous: the app colours by factor, and >=0.9 is where HieraType's
  # own gate sits, so the top band IS ht_call and the rest is what the gate threw away.
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

  # Which methods call a cell B-lineage. Five vocabularies for one compartment, so the
  # disagreement is invisible until they are put side by side. Handed back as a logical
  # matrix rather than as finished columns: WHICH methods count as "the others" is a
  # question the user answers in the app, and precomputing every subset is 2^n columns.
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

  # Both scales when the bundle carries them. Normalized is the default for looking at
  # relative expression; raw is what detection questions must use, because normalization
  # is library-size scaled and inflates a single transcript in a low-depth cell.
  read_expr <- function(fn) {
    f <- file.path(ex, fn)
    if (!file.exists(f)) return(NULL)
    # attach, don't just load: `[` on a dgCMatrix is an S4 method, and under
    # Rscript the namespace alone leaves it dispatching to base `[`
    library(Matrix)
    m <- readRDS(f)
    m[, match(cells$cell_ID, colnames(m)), drop = FALSE]
  }
  expr     <- read_expr("marker_expr_normalized.rds")
  expr_raw <- read_expr("marker_expr_raw.rds")

  # IF channels ride in both matrices under an "IF:" prefix, which groups them in the gene
  # picker and cannot collide with a gene symbol. Identical on both scales -- an intensity is
  # not a count, so library-size normalization does not apply to it.
  # Every numeric metadata column is a quantity, so it is served through the gene machinery
  # and gets the gradient, the threshold filter and the scale toggle for free.
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

  # poly$i indexes cells BY ROW, so it is only valid against the cells file it
  # shipped with. Both are written in the same chunk; never mix bundles.
  f <- file.path(ex, "cell_polygons.rds")
  poly <- if (file.exists(f)) as.data.frame(readRDS(f)) else NULL

  # Per-cluster triage tables, keyed by the cluster column they describe. Each step keeps
  # its own under markers/: clusters, reference calls and subclusters.
  dd <- file.path(root, "tables", c("03_clustering", "04_reference_typing", "05_subclustering"),
                  "markers")
  diag <- de <- list()
  for (f in list.files(dd, "^cluster_labels\\..+\\.template\\.csv$", full.names = TRUE)) {
    col <- sub("^cluster_labels\\.(.+)\\.template\\.csv$", "\\1", basename(f))
    t <- utils::read.csv(f)
    t$cluster_id <- as.character(t$cluster_id)
    diag[[col]] <- t[order(-t$n_cells), ]
    g <- file.path(dirname(f), paste0(col, "_de_top.csv"))
    if (file.exists(g)) {
      x <- utils::read.csv(g)
      x$cluster <- as.character(x$cluster)
      de[[col]] <- x
    }
  }

  list(cells = cells, poly = poly, expr = expr, expr_raw = expr_raw,
       diag = diag, de = de, lineage = lineage)
}
