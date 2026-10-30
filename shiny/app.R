library(shiny)
library(ggplot2)
library(scattermore)
source("theme.R")

# The label rules annotate applies, so the live columns here are what it will write
source("../bin/label_state.R")
# The figures' white slash across a low-confidence label (GeomSlash, dot_r)
source("../bin/cosmx_figures.R")

# Which results tree to read: COSMX_NF when it is set, otherwise the run the pipeline ran last
# among those beside the pipeline (../../Runs) with an app bundle. Its stage log gains a line at
# every launch; a run without one counts from its bundle. Set COSMX_NF to read another run
latest_results <- function(runs = file.path("..", "..", "Runs")) {
  b <- Sys.glob(file.path(runs, "*", "results", "local_export", "cells_umap_clusters.rds"))
  if (!length(b)) return(NULL)
  res <- dirname(dirname(b))
  last <- vapply(res, function(r) {
    f <- file.path(r, c("pipeline_info/stage_log.tsv", "local_export/cells_umap_clusters.rds"))
    as.numeric(max(file.mtime(f[file.exists(f)])))
  }, numeric(1))
  normalizePath(res[which.max(last)])
}
NF <- Sys.getenv("COSMX_NF")
if (!nzchar(NF)) NF <- latest_results() %||% "../results"
# The run's name, shown at the top of the sidebar
RUN_NAME <- basename(dirname(normalizePath(NF, mustWork = FALSE)))
message("CosMx explorer: ", normalizePath(NF, mustWork = FALSE))

src <- if (dir.exists(file.path(NF, "local_export"))) {
  source("load_cosmx_nf.R"); load_cosmx_nf(NF)
} else if (file.exists("demo_cells.rds")) {
  list(cells = readRDS("demo_cells.rds"), poly = readRDS("demo_poly.rds"),
       expr = NULL, diag = list(), de = list())
} else {
  stop("No data. Either rsync a bundle into ", normalizePath(NF, mustWork = FALSE),
       "/local_export (build it with shiny/export_app_bundle.sh <run_id> on the server), ",
       "or run make_demo.R for fixtures.")
}
cells <- src$cells; poly <- src$poly; expr <- src$expr; expr_raw <- src$expr_raw
DIAG  <- src$diag;  DE   <- src$de

# Lineage agreement: a logical cells x methods matrix per lineage, built by the loader.
# The colour column itself is assembled in the server from the primary method and whichever
# others are ticked, so it is a reactive vector that is not, and cannot be, in `cells`.
LIN     <- src$lineage %||% list()
LIN_NM  <- if (length(LIN)) names(LIN)[1] else NULL
LIN_COL <- if (is.null(LIN_NM)) NULL else paste0(LIN_NM, "_agreement")
LIN_M   <- if (is.null(LIN_NM)) character() else colnames(LIN[[LIN_NM]])

pick <- function(...) { for (n in c(...)) if (n %in% names(cells)) return(n); NULL }
emb  <- function(a, b) if (all(c(a, b) %in% names(cells))) c(a, b)

# One tab per subclustered parent. Its cells were re-embedded on their own PCA, so they do not
# belong on the harmony UMAP -- the structure there is the parent's, not theirs.
# Per-parent subcluster embeddings, named "<basis> · cluster <n>" and ordered by basis, then
# parent. They share one Subclusters tab, with the parent picked from a list.
SUB_EMB <- {
  u1 <- grep("_c[^_]+_UMAP1$", names(cells), value = TRUE)
  basis <- sub("_c[^_]+_UMAP1$", "", u1)
  par   <- sub("^.*_c([^_]+)_UMAP1$", "\\1", u1)
  ord   <- order(basis, suppressWarnings(as.numeric(par)), par)
  stats::setNames(lapply(u1[ord], function(a) c(a, sub("1$", "2", a))),
                  paste0(basis[ord], " · cluster ", par[ord], recycle0 = TRUE))
}
SUB_EMB <- SUB_EMB[vapply(SUB_EMB, function(z) all(z %in% names(cells)), logical(1))]
SUB_TAB <- "Subclusters"
SUB_COLS <- unique(sub(" · cluster .*$", "", names(SUB_EMB)))

EMB <- Filter(Negate(is.null), c(list(
  Spatial = emb("x", "y"),
  UMAP    = emb("UMAP1", "UMAP2") %||% emb("UMAP_1", "UMAP_2"),
  # The same cells embedded on PCA instead of harmony. Present only when the run corrected
  # for batch, and carried in the same cells table so colour, gene and sample filters hold
  # across the pair -- switching tabs is then a clean before/after of the correction alone.
  "UMAP pre-batch" = emb("UMAPpre1", "UMAPpre2"),
  PCA     = emb("PCA_1", "PCA_2")
), if (length(SUB_EMB)) stats::setNames(list(SUB_EMB[[1]]), SUB_TAB)))
# Output ids, brush ids and the per-page facet memory are keyed by a slug: a display name
# with a space or a hyphen in it cannot be part of an input id.
VID <- c(Spatial = "spatial", UMAP = "umap",
         "UMAP pre-batch" = "umap_pre", PCA = "pca", Subclusters = "subclusters")[names(EMB)]
names(VID) <- names(EMB)
# Manual edits from the Annotate menu are laid over this column and shown as <column>_edited.
# "(cluster label)" as a label hands a cell back to its cluster's label.
BASE_LAB <- pick("cell_type", grep("^label_", names(cells), value = TRUE))
EDIT_COL <- if (is.null(BASE_LAB)) NULL else paste0(BASE_LAB, "_edited")
KEEP_LAB <- KEEP_LABEL

# The run's curated labels, beside its results: what annotate applies (main.nf's labels_dir).
# Cluster labels, per-cell edits and the chosen clustering are saved here as they are made
LABEL_DIR  <- normalizePath(Sys.getenv("COSMX_LABELS", file.path(NF, "..", "labels")),
                            mustWork = FALSE)
OV_FILE    <- file.path(LABEL_DIR, "label_overrides.csv")
BASIS_FILE <- file.path(LABEL_DIR, "cluster_basis.txt")
cl_file    <- function(b) file.path(LABEL_DIR, paste0("cluster_labels.", b, ".csv"))
read_labs  <- function(dir) {
  f <- list.files(dir, "^cluster_labels\\..+\\.csv$", full.names = TRUE)
  f <- f[!grepl("\\.(draft|template)\\.csv$", f)]
  stats::setNames(lapply(f, read_cluster_labels),
                  sub("^cluster_labels\\.(.+)\\.csv$", "\\1", basename(f)))
}

# Cluster-level labels, edited in Tables > Cluster labels. Each clustering starts from its saved
# labels, else from the pipeline's draft in tables/06_annotation/labels, else from the majority
# current label of its cells.
CL_BASES <- grep("^sub_|^leiden_clus_res", names(Filter(is.factor, cells)), value = TRUE)
LABS_SAVED <- read_labs(LABEL_DIR)
LABS <- utils::modifyList(read_labs(file.path(NF, "tables", "06_annotation", "labels")), LABS_SAVED)
# The clustering cell_type comes from: cluster_basis.txt, else the only saved labels
CHOSEN0 <- if (file.exists(BASIS_FILE)) trimws(readLines(BASIS_FILE, n = 1L, warn = FALSE)) else
  if (length(LABS_SAVED) == 1L) names(LABS_SAVED) else NULL
CL_DEFAULT <- if (isTRUE(CHOSEN0 %in% CL_BASES)) CHOSEN0 else
  if ("sub_res0.3" %in% CL_BASES) "sub_res0.3" else CL_BASES[1]
# The labels in Tables > Cluster labels laid over the cells of the clustering picked there, with
# the Annotate edits on top: the cell_type annotate will write. The colouring the app opens on.
# Beside it the label figures will draw (low-confidence edits included) and its confidence
CUSTOM_COL  <- if (length(CL_BASES)) "custom_annotation" else NULL
DISPLAY_COL <- if (length(CL_BASES)) "custom_display" else NULL
CONF_COL    <- if (length(CL_BASES)) "custom_confidence" else NULL
LIVE_COLS   <- c(CUSTOM_COL, DISPLAY_COL, CONF_COL)
# Columns of cell-type labels, coloured from assets/label_colours.csv
LAB_COLS <- unique(c(CUSTOM_COL, DISPLAY_COL, BASE_LAB, EDIT_COL,
                     intersect("cell_type_display", names(cells)),
                     grep("^label_", names(cells), value = TRUE)))
# The marker-panel types each label covers. A cluster whose marker best guess falls outside its
# label's set is flagged for review.
LABEL_OK <- list(
  PT = "PT", PT_injured = c("PT_injured", "PT"), TAL_DCT = c("TAL", "DCT", "Macula_densa"),
  CNT_CD = c("CNT", "Principal_cell", "IC"), Podocyte = "Podocyte",
  Endothelial = c("Endothelial", "Glomerular_endothelium", "Peritubular_capillary_EC",
                  "Vasa_recta_EC", "Arteriolar_EC", "Lymphatic_EC"),
  Fibroblast = "Fibroblast", Myofibroblast = "Myofibroblast", Pericyte = "Pericyte",
  Mural = c("vSMC", "Pericyte", "Mesangial", "Renin_JG"),
  Macrophage = c("Macrophage", "Monocyte", "cDC", "pDC"),
  Lymphocyte = c("T_cell", "CD4_T", "CD8_T", "Treg", "NK", "B_cell"),
  T_NK = c("T_cell", "CD4_T", "CD8_T", "Treg", "NK"),
  Plasma_cell = "Plasma")
CAT   <- c(CUSTOM_COL, DISPLAY_COL, CONF_COL, names(Filter(is.factor, cells)), LIN_COL, EDIT_COL)
SAMP  <- pick("sample_id", "sample")
# A column that never varies within a sample is sample-level metadata. Colouring tissue by
# it just recolours the block of space each sample occupies, so offer those on embeddings
# only. Derived rather than listed by name, so condition columns travel whatever they are
# called. Split by is deliberately NOT filtered -- see the emb observer.
# nlevels > 1 guards against a degenerate clustering column (res0.02 collapses to one
# cluster, so it is constant within every sample and would be misread as metadata).
GRP <- if (is.null(SAMP)) character() else {
  const_in_sample <- function(z)
    nlevels(droplevels(z)) > 1 &&
      all(tapply(z, cells[[SAMP]], function(v) length(unique(v[!is.na(v)])) <= 1))
  # intersect with names(cells): CAT also carries the custom, lineage and edited columns,
  # which are reactive and have no column to test
  c(SAMP, names(Filter(const_in_sample, cells[intersect(setdiff(CAT, SAMP), names(cells))])))
}
# GRP still drives Split by; on Spatial the samples are stacked, so colouring by treatment or
# timepoint blocks them visibly and is worth having.
cat_for <- function(view) CAT

# selectInput renders a named list as <optgroup>s. Group by what the column IS, so the
# annotation columns stop being an undifferentiated wall of names.
GROUPS <- list(
  # The labels as the next annotate will write them, edits included
  "Labels - live (labels folder)" = function(x) x %in% LIVE_COLS,
  "Metadata"                 = function(x) x %in% GRP,
  "Clustering"               = function(x) grepl("^leiden_clus_res", x),
  "Clustering - subclusters" = function(x) grepl("^sub_res", x),
  "Annotation - manual"      = function(x) grepl("^label_|^cell_type", x),
  "Annotation - HieraType"   = function(x) x == "ht_call",
  "Annotation - InSituType unsupervised" = function(x) grepl("^insitutype_unsup", x),
  "Annotation - InSituType supervised"   = function(x) grepl("^insitutype_sup", x),
  "Annotation - SingleR"     = function(x) grepl("^monaco_|^singler", x),
  "Derived - lineage agreement" = function(x) grepl("_agreement$", x))

grouped_choices <- function(cols) {
  out <- list(); left <- cols
  for (g in names(GROUPS)) {
    hit <- left[GROUPS[[g]](left)]
    # as.list(): Shiny flattens a length-1 group into a plain option, so a group holding a
    # single column would silently lose its heading
    if (length(hit)) { out[[g]] <- as.list(hit); left <- setdiff(left, hit) }
  }
  if (length(left)) out[["Other"]] <- as.list(left)
  out
}
SPLIT <- c("Default" = "none", intersect(c("sample_id", "treatment", "timepoint", "batch"), names(cells)))

# FOV is read from the cell_ID, <sample>-c_<slide>_<fov>_<cell> in the merged object
if (!"fov" %in% names(cells) && all(grepl("c_\\d+_\\d+_\\d+$", cells$cell_ID)))
  cells$fov <- as.integer(sub("^.*c_\\d+_(\\d+)_\\d+$", "\\1", cells$cell_ID))
FOVS <- if (is.null(SAMP) || !"fov" %in% names(cells)) NULL else {
  u <- unique(data.frame(s = as.character(cells[[SAMP]]), f = cells$fov))
  u <- u[order(match(u$s, levels(cells[[SAMP]])), u$f), ]
  lapply(split(u, factor(u$s, levels = unique(u$s))),
         function(z) stats::setNames(paste(z$s, z$f, sep = "|"), paste("FOV", z$f)))
}
# Spatial shows one sample at a time; this choice stacks them all instead
ALL_SAMPLES <- ".all"
GENES <- if (is.null(expr)) character() else rownames(expr)
GREY  <- "#DCE3E0"
# --- marker evidence -------------------------------------------------------------------
# A label is a claim; the markers are the evidence for it. Panels come from the repo's own
# canonical_markers.csv so the app scores against the same vocabulary DIAGNOSTICS does.
MK <- tryCatch(utils::read.csv("../assets/canonical_markers.csv"), error = function(e) NULL)
mk_panels <- function(by) {
  if (is.null(MK) || is.null(expr_raw)) return(list())
  g <- lapply(split(MK$gene, MK[[by]]),
              function(z) intersect(unique(z), rownames(expr_raw)))
  g[lengths(g) >= 2]                      # one gene is an anecdote, not a panel
}
MK_SET <- if (is.null(MK)) list() else
  list(compartment = mk_panels("compartment"), cell_type = mk_panels("cell_type"))
MK_GENES <- unique(unlist(MK_SET, use.names = FALSE))
# Detection tracks sequencing depth, so a naive fold reports depth as biology. Deciles of
# total counts hold depth fixed; the expected count for a group is then the background rate
# in each decile weighted by how that group is distributed across them.
DEC <- if (!"total_counts" %in% names(cells) || anyNA(cells$total_counts)) NULL else {
  br <- unique(stats::quantile(cells$total_counts, seq(0, 1, 0.1), na.rm = TRUE))
  if (length(br) > 2) cut(cells$total_counts, br, include.lowest = TRUE, labels = FALSE)
  else NULL
}
NDEC <- if (is.null(DEC)) 1L else max(DEC)
# Second recessive tier. "others only" is not the subject but is not background either, and
# at 59,296 cells it takes over the panel if it is given a palette colour.
GREY2 <- "#A8B6C0"
# Per-cell posteriors, joined from the published stage CSVs. Offered as a threshold rather
# than as bands on the legend: where to cut is the question, so it should be a control, not
# a decision baked into a column. The legend is then free to carry the cell types.
SCORE <- grep("_score$", names(cells), value = TRUE)
# ht_call's score is published as ht_score; everything else is <column>_score
score_of <- function(col) {
  if (is.null(col)) return(NULL)
  s <- if (identical(col, "ht_call")) "ht_score" else paste0(col, "_score")
  if (s %in% SCORE) s else NULL
}

# a label_<res> column is a renaming of leiden_clus_<res>; show that basis's triage
basis_of <- function(c) if (grepl("^label_", c)) paste0("leiden_clus_", sub("^label_", "", c)) else c

pal_n <- function(k) if (k <= length(pal_cat)) pal_cat[seq_len(k)] else
  c(pal_cat, scales::hue_pal(h = c(200, 560), l = 62, c = 72)(k - length(pal_cat)))

# keyed by LEVEL NAME, rebuilt per column: cluster 3 at res0.3 is not cluster 3 at
# res0.8, so the two must not be forced to share a colour.
# A level that means "not this" is background and gets GREY rather than a palette slot:
# 94,179 cells painted spring green swamp the 514 that are the point, and burn a colour
# the real classes need. Named by convention -- "none", "no", or anything "not ...".
pal_for <- function(f) {
  lv <- levels(f)
  bg  <- lv %in% c("none", "no", "unknown") | grepl("^not ", lv)
  bg2 <- lv %in% "others only"
  p   <- setNames(character(length(lv)), lv)
  p[!bg & !bg2] <- pal_n(sum(!bg & !bg2))
  p[bg]  <- GREY
  p[bg2] <- GREY2
  p
}

# One colour per cell-type label, shared with the figures: kidney types pastel, immune types
# vivid, untyped grey. Keyed by label, so a label keeps its colour whatever else is present.
# A label missing from the file takes the next pal_wide colour the file does not use.
LAB_PAL <- tryCatch(utils::read.csv("../assets/label_colours.csv", colClasses = "character"),
                    error = function(e) NULL)
pal_labels <- function(lv) {
  p <- if (is.null(LAB_PAL)) rep(NA_character_, length(lv)) else LAB_PAL$colour[match(lv, LAB_PAL$label)]
  new <- is.na(p)
  p[new] <- rep_len(setdiff(pal_wide, LAB_PAL$colour), sum(new))
  stats::setNames(p, lv)
}
# Draw order for labels: untyped cells underneath, immune cells on top
lab_layer <- function(v) {
  if (is.null(LAB_PAL)) return(rep(1L, length(v)))
  g <- LAB_PAL$group[match(as.character(v), LAB_PAL$label)]
  ifelse(g %in% "untyped", 0L, ifelse(g %in% "immune", 2L, 1L))
}

# Pinch or Cmd/Ctrl + scroll zooms; plain scrolling is left to the page. Deltas are summed and
# sent once the gesture pauses, so a flick costs one re-render instead of thirty. Chrome
# reports pinch as a ctrl-wheel, Safari as gesture events.
ZOOM_JS <- "
(function () {
  function flush(el) {
    clearTimeout(el._t);
    el._t = setTimeout(function () {
      if (el._dy) Shiny.setInputValue(el.id + '_wheel', {dy: el._dy, n: Date.now()}, {priority: 'event'});
      el._dy = 0;
    }, 120);
  }
  function flushPan(el) {
    clearTimeout(el._pt);
    el._pt = setTimeout(function () {
      if (el._px || el._py) Shiny.setInputValue(el.id + '_pan',
        {dx: el._px || 0, dy: el._py || 0, w: el.clientWidth, h: el.clientHeight, n: Date.now()},
        {priority: 'event'});
      el._px = 0; el._py = 0;
    }, 60);
  }
  var zoomed = false;
  // The height the plot may take is its card body's, measured: the window, the Tables panel and
  // full screen all change it. A hidden tab measures 0 and is skipped.
  function sendAvail(area) {
    clearTimeout(area._t);
    area._t = setTimeout(function () {
      if (area.offsetParent === null) return;
      var cs = getComputedStyle(area);
      var h = area.clientHeight - parseFloat(cs.paddingTop) - parseFloat(cs.paddingBottom);
      if (h > 50) Shiny.setInputValue('cosmx_avail', Math.floor(h) - 4);
    }, 150);
  }
  var ro = new ResizeObserver(function (es) { es.forEach(function (e) { sendAvail(e.target); }); });
  // Figure size: Cmd/Ctrl + and - scale the plot's text, points and lines, Cmd/Ctrl 0 resets.
  // The page, the legend and the controls keep their size.
  var SC = [0.5, 0.67, 0.75, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 1.75, 2], si = 5;
  try { var v = SC.indexOf(parseFloat(localStorage.getItem('cosmx_fig'))); if (v >= 0) si = v; } catch (e) {}
  function setFig(i) {
    si = Math.max(0, Math.min(SC.length - 1, i));
    Shiny.setInputValue('fig_scale', SC[si]);
    document.querySelectorAll('.cosmx-figpct').forEach(function (el) {
      el.textContent = Math.round(SC[si] * 100) + '%';
    });
    try { localStorage.setItem('cosmx_fig', SC[si]); } catch (e) {}
  }
  window.cosmxFig = function (d) { setFig(d === 0 ? 5 : si + d); };
  document.addEventListener('keydown', function (e) {
    if (!(e.metaKey || e.ctrlKey) || e.altKey) return;
    var k = e.key, c = e.code;
    var d = (k === '-' || k === '_' || c === 'Minus' || c === 'NumpadSubtract') ? -1 :
            (k === '=' || k === '+' || c === 'Equal' || c === 'NumpadAdd') ? 1 :
            (k === '0' || c === 'Digit0' || c === 'Numpad0') ? 0 : null;
    if (d === null) return;
    e.preventDefault(); e.stopPropagation();
    window.cosmxFig(d);
  }, true);
  $(document).on('shiny:connected', function () {
    Shiny.addCustomMessageHandler('cosmx-zoomed', function (z) { zoomed = z; });
    document.querySelectorAll('.cosmx-plot-area').forEach(function (a) { ro.observe(a); });
    setFig(si);
  });
  function target(e) { return e.target.closest ? e.target.closest('.cosmx-zoom') : null; }
  document.addEventListener('wheel', function (e) {
    var el = target(e); if (!el) return;
    var k = e.deltaMode === 1 ? 16 : 1;
    if (e.ctrlKey || e.metaKey) {
      e.preventDefault();
      if (Date.now() - (el._g || 0) < 200) return;
      el._dy = (el._dy || 0) + e.deltaY * k;
      flush(el);
      return;
    }
    // Zoomed: sideways swipes pan; up/down pans only when the plot fits its box,
    // otherwise it stays with the page
    if (!zoomed) return;
    var box = el.closest('.card-body');
    var fits = !box || el.offsetHeight <= box.clientHeight + 4;
    var horiz = Math.abs(e.deltaX) > Math.abs(e.deltaY);
    if (!horiz && !fits) return;
    e.preventDefault();
    if (horiz) el._px = (el._px || 0) + e.deltaX * k;
    else el._py = (el._py || 0) + e.deltaY * k;
    flushPan(el);
  }, {passive: false});
  ['gesturestart', 'gesturechange'].forEach(function (type) {
    document.addEventListener(type, function (e) {
      var el = target(e); if (!el) return;
      e.preventDefault();
      el._g = Date.now();
      if (type === 'gesturestart') { el._s = 1; return; }
      el._dy = (el._dy || 0) - Math.log(e.scale / (el._s || 1)) * 400;
      el._s = e.scale;
      flush(el);
    }, {passive: false});
  });
})();
"

# Draw mode: press and drag on the plot to outline cells. The path is kept in CSS pixels of
# the plot image and sent on release; the server maps it to data coordinates with the
# plot's own coordmap. The outline is cleared when the plot redraws, which is when the
# server draws it again in data coordinates.
LASSO_JS <- "
(function () {
  var st = null;
  function add(e) {
    var x = e.clientX - st.r.left, y = e.clientY - st.r.top, p = st.pts;
    if (p.length && Math.abs(p[p.length - 1][0] - x) + Math.abs(p[p.length - 1][1] - y) < 3) return;
    p.push([x, y]);
    st.line.setAttribute('points', p.map(function (q) { return q.join(','); }).join(' '));
  }
  document.addEventListener('mousedown', function (e) {
    var box = e.target.closest ? e.target.closest('.cosmx-draw') : null;
    if (!box || e.button !== 0) return;
    var im = box.querySelector('img'); if (!im) return;
    e.preventDefault(); e.stopPropagation();
    var svg = box.querySelector('svg.cosmx-lasso');
    if (!svg) {
      svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
      svg.setAttribute('class', 'cosmx-lasso');
      box.appendChild(svg);
    }
    svg.style.cssText = 'position:absolute;pointer-events:none;left:' + im.offsetLeft + 'px;top:' +
      im.offsetTop + 'px;width:' + im.clientWidth + 'px;height:' + im.clientHeight + 'px;';
    svg.innerHTML = '<polyline fill=\\\"rgba(0,0,0,0.06)\\\" stroke=\\\"#000\\\" ' +
      'stroke-width=\\\"1.5\\\" stroke-dasharray=\\\"4 3\\\" points=\\\"\\\"/>';
    st = {box: box, im: im, r: im.getBoundingClientRect(), pts: [], line: svg.firstChild};
    add(e);
  }, true);
  document.addEventListener('mousemove', function (e) { if (st) add(e); }, true);
  document.addEventListener('mouseup', function () {
    if (!st) return;
    var s = st; st = null;
    if (s.pts.length > 2)
      Shiny.setInputValue(s.box.id + '_lasso', {pts: s.pts, n: Date.now()}, {priority: 'event'});
    else s.line.setAttribute('points', '');
  }, true);
  $(document).on('shiny:value', function (e) {
    var el = document.getElementById(e.name);
    var svg = el ? el.querySelector('svg.cosmx-lasso') : null;
    if (svg) svg.innerHTML = '';
  });
})();
"

# The Tables panel under the plot. A click on its bar opens or closes it; dragging the bar
# resizes it, up for more table and down for more plot. Shiny is told when its outputs show.
TABLES_JS <- "
(function () {
  var drag = null;
  function open(p, on) {
    p.classList.toggle('open', on);
    Shiny.setInputValue('tables_open', on);
    $(p).trigger(on ? 'shown' : 'hidden');
  }
  document.addEventListener('mousedown', function (e) {
    var bar = e.target.closest ? e.target.closest('.cosmx-tables-bar') : null;
    if (!bar || e.button !== 0) return;
    var p = bar.parentElement, body = p.querySelector('.cosmx-tables-body');
    var h0 = p.classList.contains('open') ? body.offsetHeight : 0;
    // The tables may grow only as far as leaves the plot 280 px
    var area = Array.prototype.find.call(document.querySelectorAll('.cosmx-plot-area'),
                                         function (a) { return a.offsetParent !== null; });
    drag = {p: p, y: e.clientY, h: h0, moved: false,
            max: h0 + Math.max(0, (area ? area.clientHeight : 400) - 280)};
    e.preventDefault();
  });
  document.addEventListener('mousemove', function (e) {
    if (!drag) return;
    var dy = drag.y - e.clientY;
    if (!drag.moved && Math.abs(dy) < 4) return;
    drag.moved = true;
    var h = Math.min(drag.max, drag.h + dy);
    if (h < 60) { if (drag.p.classList.contains('open')) open(drag.p, false); return; }
    drag.p.style.setProperty('--cosmx-tables-h', Math.round(h) + 'px');
    if (!drag.p.classList.contains('open')) open(drag.p, true);
  });
  document.addEventListener('mouseup', function () {
    if (!drag) return;
    var d = drag; drag = null;
    if (!d.moved) open(d.p, !d.p.classList.contains('open'));
    else $(d.p).trigger('shown');
  });
  document.addEventListener('keydown', function (e) {
    var bar = e.target.closest ? e.target.closest('.cosmx-tables-bar') : null;
    if (!bar || (e.key !== 'Enter' && e.key !== ' ')) return;
    e.preventDefault();
    open(bar.parentElement, !bar.parentElement.classList.contains('open'));
  });
})();
"
CHEVRON <- paste0('<svg width="14" height="14" viewBox="0 0 16 16" aria-hidden="true">',
                  '<path d="M3 6l5 5 5-5" fill="none" stroke="currentColor" stroke-width="1.6"/></svg>')

# Point in polygon by ray casting, vectorised over points and looped over edges
in_poly <- function(px, py, vx, vy) {
  inside <- logical(length(px))
  j <- length(vx)
  for (i in seq_along(vx)) {
    cr <- ((vy[i] > py) != (vy[j] > py)) &
      (px < (vx[j] - vx[i]) * (py - vy[i]) / (vy[j] - vy[i]) + vx[i])
    inside <- xor(inside, cr %in% TRUE)
    j <- i
  }
  inside
}

# An overrides CSV: cell_ID and label, confidence and note optional. Cells not in this bundle
# are dropped.
read_overrides <- function(f) {
  d <- read_label_overrides(f)
  d[d$cell_ID %in% cells$cell_ID, , drop = FALSE]
}

APP_CSS <- "
.cosmx-draw { position: relative; cursor: crosshair; }
html, body { overscroll-behavior-x: none; }
.selbar { position: sticky; top: 0; z-index: 20; display: flex; gap: 8px; align-items: center;
          background: #fff; padding: 4px 6px; border-bottom: 1px solid #E4E9E7; }
.plot-toolbar .dropdown-menu { font-size: .9rem; }
.legend-row { white-space: nowrap; font-size: .88rem; margin-bottom: 2px; }
.legend-sw { display: inline-block; width: 11px; height: 11px; border-radius: 2px;
             margin-right: 7px; vertical-align: -1px; }
.legend-low { background: repeating-linear-gradient(135deg, #9A9A9A 0 2px, #FFFFFF 2px 4px); }
/* The plot's card body scrolls on its own (All samples, split views); the legend stays put */
.cosmx-plot-area { overflow: auto; }
/* The Tables panel: a bar to click or drag, over a body whose height the drag sets */
.cosmx-tables { flex: 0 0 auto; }
.cosmx-tables-bar { display: flex; align-items: center; gap: .75rem; padding: .5rem .9rem;
  border: 1px solid var(--bs-border-color); border-radius: .5rem; background: var(--bs-body-bg);
  cursor: row-resize; user-select: none; }
.cosmx-tables-bar:hover { border-color: #00A870; }
.cosmx-grip { width: 2.5rem; height: 4px; border-radius: 2px; background: #C9D3CF; margin: 0 auto; }
.cosmx-chev { display: inline-flex; transition: transform .15s; }
.cosmx-tables.open .cosmx-chev { transform: rotate(180deg); }
.cosmx-tables-body { display: none; height: var(--cosmx-tables-h, 40vh); margin-top: .5rem; }
.cosmx-tables.open .cosmx-tables-body { display: flex; flex-direction: column; }
.cosmx-tables-body > .card { flex: 1 1 auto; min-height: 0; margin-bottom: 0; }
/* Each table fills its tab and scrolls inside it */
.cosmx-scroll { overflow: auto; font-size: .9em; }
.cosmx-scroll table { margin-bottom: 0; }
.cosmx-figscale .btn { padding: 0 .5rem; }
.cosmx-figpct { min-width: 3.6em; }
.cosmx-scroll table thead th { position: sticky; top: 0; z-index: 2; background: #FFFFFF;
                               box-shadow: inset 0 -1px 0 #DEE2E6; }
.cosmx-scroll table td:first-child, .cosmx-scroll table th:first-child {
  position: sticky; left: 0; z-index: 1; background: #FFFFFF; }
.cosmx-scroll table thead th:first-child { z-index: 3; }
"

scroll_box <- function(...) htmltools::bindFillRole(div(class = "cosmx-scroll", ...), item = TRUE)

# Bootstrap dropdown that stays open while its inputs are used. Fixed positioning lets the
# menu escape the card's overflow clipping.
# Widths are in rem, so they shrink with the page's root font size (theme.R)
dd <- function(label, ..., width = "20rem") div(
  class = "dropdown",
  tags$button(class = "btn btn-sm btn-outline-primary dropdown-toggle", type = "button",
              `data-bs-toggle` = "dropdown", `data-bs-auto-close` = "outside",
              `data-bs-popper-config` = '{"strategy":"fixed"}', label),
  div(class = "dropdown-menu p-3", style = paste0("width:", width, ";"), ...))

toolbar <- div(
  class = "plot-toolbar d-flex flex-wrap gap-2 align-items-center",
  # What a drag on the plot does. Kept as an explicit switch so a click meant to inspect a cell
  # never starts a box, and a box meant to select never zooms.
  tagAppendAttributes(
    div(style = "margin-bottom: -1rem;",
        radioButtons("mode", NULL, c("Inspect" = "inspect", "Select" = "select", "Draw" = "draw"),
                     "inspect", inline = TRUE)),
    title = paste("Inspect: click a cell; drag a box to zoom.",
                  "Select: drag a box to filter or label the cells in it.",
                  "Draw: draw around cells to filter or label them.",
                  "In every mode pinch or Cmd/Ctrl + scroll zooms and double-click resets.")),
  if (!is.null(SAMP)) dd(
    "Samples",
    checkboxGroupInput("samples", NULL, levels(cells[[SAMP]]), selected = levels(cells[[SAMP]])),
    actionLink("samples_all", "All"), " · ", actionLink("samples_none", "None"),
    width = "16.25rem"),
  if (length(SUB_EMB)) conditionalPanel(
    sprintf("input.view == '%s'", SUB_TAB),
    div(style = "width: 14.4rem; margin-bottom: -1rem;",
        selectInput("sub_parent", NULL, names(SUB_EMB), selectize = FALSE))),
  if (!is.null(SAMP)) conditionalPanel(
    "input.view == 'Spatial'",
    tagAppendAttributes(
      div(style = "width: 12.5rem; margin-bottom: -1rem;",
          selectInput("one_sample", NULL,
                      c(stats::setNames(levels(cells[[SAMP]]), levels(cells[[SAMP]])),
                        "All samples" = ALL_SAMPLES),
                      levels(cells[[SAMP]])[1], selectize = FALSE)),
      title = "The sample Spatial shows. All samples stacks them, which draws far more cells")),
  if (!is.null(FOVS)) conditionalPanel(
    "input.view == 'Spatial'",
    div(style = "width: 13.1rem; margin-bottom: -1rem;",
        selectizeInput("fov", NULL, choices = c("Jump to FOV" = "", FOVS),
                       # in body, so the list opens over the plot, not under it
                       options = list(placeholder = "Jump to FOV", dropdownParent = "body")))),
  dd(
    "Filter",
    selectizeInput("keep_lv", "Keep levels of Colour by", NULL, multiple = TRUE,
                   options = list(placeholder = "all", plugins = list("remove_button"))),
    if (length(GENES)) tagList(
      selectizeInput("fgene", "Keep cells expressing", NULL,
                     options = list(placeholder = "gene", allowEmptyOption = TRUE)),
      numericInput("fmin", "above", 0, step = 0.5)),
    uiOutput("conf_ui"),
    uiOutput("sel_status"),
    actionButton("filter_reset", "Clear all filters", class = "btn-sm btn-outline-secondary w-100")),
  dd(
    "Highlight",
    selectizeInput("hl_lv", "Levels of Colour by", NULL, multiple = TRUE,
                   options = list(placeholder = "none", plugins = list("remove_button"))),
    if (length(GENES)) tagList(
      selectizeInput("gene", "Gene", NULL,
                     options = list(placeholder = "none", allowEmptyOption = TRUE)),
      conditionalPanel("input.gene", numericInput("emin", "Colour only cells above", NA, step = 0.5)),
      # One stray transcript looks like real expression; co-detection against an
      # independence baseline is the question a single gene cannot answer
      selectizeInput("genes", "Gene combination (2+)", NULL, multiple = TRUE,
                     options = list(placeholder = "none", plugins = list("remove_button"))),
      conditionalPanel(
        "input.genes && input.genes.length > 1",
        numericInput("gthr", "Expressed if value >", 0, min = 0, step = 0.5),
        radioButtons("gmode", "Show",
                     c("How many are expressed" = "count",
                       "Each +/- combination"   = "pattern",
                       "All of them"            = "all",
                       "Any of them"            = "any"), "count"),
        # Pattern mode enumerates 2^n groups, so it is capped at 3 genes (8 groups)
        conditionalPanel("input.gmode == 'pattern'",
                         selectInput("hlpat", "Highlight one combination",
                                     choices = "(show all)", selected = "(show all)")))),
    width = "21.25rem"),
  if (!is.null(EDIT_COL)) dd(
    "Annotate",
    selectizeInput("ann_label", paste("Label for", BASE_LAB),
                   c("", levels(cells[[BASE_LAB]]), KEEP_LAB),
                   options = list(create = TRUE, placeholder = "pick or type a new label")),
    radioButtons("ann_conf", "Confidence",
                 c("High: sets cell_type" = "high",
                   "Low: figures only, striped" = "low"), "high"),
    textInput("ann_note", "Note (optional)"),
    actionButton("ann_assign", "Assign to selected cells", class = "btn-sm btn-success w-100 mb-2"),
    uiOutput("ann_summary"),
    div(class = "mb-2", actionLink("ann_undo", "Undo last"), " · ", actionLink("ann_clear", "Remove all")),
    downloadButton("ann_dl", "Download overrides CSV", class = "btn-sm btn-outline-primary w-100 mb-2"),
    fileInput("ann_up", "Load an overrides CSV", accept = ".csv"),
    div(class = "text-muted small",
        "Edits save as you go to ", code(OV_FILE), ", which annotate applies on its next run. ",
        "A low-confidence label leaves cell_type alone: figures draw it, striped."),
    width = "21.25rem"),
  tagAppendAttributes(
    actionButton("zoom_reset", "Reset zoom", class = "btn-sm btn-outline-secondary"),
    title = "Zoom: drag a box in Inspect mode, pinch, or Cmd/Ctrl + scroll. Double-click to reset."),
  # Figure size: the plot's text, points and lines. The page and the legend keep theirs
  div(class = "btn-group btn-group-sm cosmx-figscale", role = "group",
      title = paste("Figure size: the plot's text, points and lines (Cmd - and Cmd +;",
                    "Cmd 0 resets). The legend and the page keep their size."),
      tags$button(type = "button", class = "btn btn-outline-secondary", onclick = "cosmxFig(-1)",
                  HTML("&minus;")),
      tags$button(type = "button", class = "btn btn-outline-secondary cosmx-figpct",
                  onclick = "cosmxFig(0)", "100%"),
      tags$button(type = "button", class = "btn btn-outline-secondary", onclick = "cosmxFig(1)", "+")),
  span(class = "text-muted small", textOutput("filter_summary", inline = TRUE)))

# The card fills the height the Tables panel leaves. The plot's body and the legend's sidebar
# scroll separately, so neither carries the other along
plot_card <- function(id) bslib::card(
  full_screen = TRUE, fill = TRUE,
  bslib::layout_sidebar(
    fillable = TRUE, fill = TRUE,
    sidebar = bslib::sidebar(
      position = "right", width = "18.75rem", open = "always",
      div(class = "cosmx-legend",
          uiOutput(paste0("cellinfo_", id)),
          uiOutput(paste0("legend_", id)))),
    bslib::card_body(fillable = FALSE, fill = TRUE, class = "cosmx-plot-area",
                     style = "padding: 6px;",
                     uiOutput(paste0("selbar_", id)),
                     uiOutput(paste0("ui_", id)))))

ui <- bslib::page_sidebar(
  title = NULL,
  window_title = "CosMx explorer",
  theme = app_theme,
  tags$head(tags$style(HTML(APP_CSS)), tags$script(HTML(ZOOM_JS)), tags$script(HTML(LASSO_JS)),
            tags$script(HTML(TABLES_JS))),
  sidebar = bslib::sidebar(
    width = "19.4rem",
    div(class = "text-muted small", "Run: ", strong(RUN_NAME)),
    selectInput("cby", "Colour by", grouped_choices(cat_for("Spatial")),
                selected = if (!is.null(CUSTOM_COL)) CUSTOM_COL
                           else pick("label_res0.8", "insitutype_unsup", "anno_auto")),
    uiOutput("conf_colour_ui"),
    if (is.null(LIN_COL)) NULL else conditionalPanel(
      sprintf("input.cby == '%s'", LIN_COL),
      checkboxGroupInput("lin_methods", "Methods", LIN_M,
                         selected = intersect(c("HT", "KPMP"), LIN_M), inline = TRUE)),
    # Raw is the correct scale for "is this cell positive?": normalization is library-size
    # scaled, so one transcript in a shallow cell outranks one in a deep cell
    if (length(GENES) && !is.null(expr_raw))
      selectInput("escale", "Expression scale", c("Normalized" = "norm", "Raw counts" = "raw"),
                  "norm", selectize = FALSE),
    conditionalPanel("input.view == 'Spatial'",
                     selectInput("layout", "Sample layout",
                                 c("True scale" = "abs", "Fit width" = "fit"), "abs",
                                 selectize = FALSE)),
    conditionalPanel("input.view != 'Spatial'",
                     selectInput("facet", "Split by", SPLIT, "none", selectize = FALSE))),

  do.call(bslib::navset_card_tab, c(
    list(id = "view", header = toolbar),
    lapply(names(EMB), function(v) bslib::nav_panel(v, plot_card(VID[[v]]))))),

  # Closed by default: the plot is the point, these are a reference. A click on the bar opens or
  # closes them; dragging it gives the tables more room (up) or the plot (down)
  div(
    id = "cosmx_tables", class = "cosmx-tables",
    div(class = "cosmx-tables-bar", role = "button", tabindex = "0",
        title = "Click to open or close. Drag up or down to resize.",
        span(class = "fw-semibold", "Tables"), span(class = "cosmx-grip"),
        span(class = "cosmx-chev", HTML(CHEVRON))),
    div(
      class = "cosmx-tables-body",
      bslib::navset_card_tab(
        id = "tables_tab",
        if (length(CL_BASES)) bslib::nav_panel(
          "Cluster labels",
          div(class = "p-2",
              div(class = "d-flex flex-wrap gap-2 align-items-end",
                  selectInput("cl_basis", "Clustering", CL_BASES, CL_DEFAULT,
                              selectize = FALSE, width = "10.6rem"),
                  selectInput("cl_pick", "Cluster", NULL, selectize = FALSE, width = "18.75rem"),
                  selectizeInput("cl_label", "Label", NULL, width = "12.5rem",
                                 options = list(create = TRUE, placeholder = "pick or type",
                                                dropdownParent = "body")),
                  div(class = "mb-3",
                      actionButton("cl_set", "Set", class = "btn-sm btn-success"),
                      actionButton("cl_next", "Next to review", class = "btn-sm btn-outline-primary"),
                      actionLink("cl_discard", "Discard edits", class = "ms-2 small")),
                  div(class = "mb-1", checkboxInput("cl_show", "Show on plot", TRUE)),
                  div(class = "mb-3",
                      actionButton("cl_use", "Use for cell_type", class = "btn-sm btn-outline-success"),
                      downloadButton("cl_dl", "Download labels CSV",
                                     class = "btn-sm btn-outline-primary"))),
              div(class = "text-muted small", textOutput("cl_hdr"))),
          scroll_box(tableOutput("cl_tbl"))),
        bslib::nav_panel("Composition", scroll_box(textOutput("comp_hdr"), tableOutput("comp"))),
        bslib::nav_panel("Cluster triage", scroll_box(textOutput("tri_hdr"), tableOutput("triage"))),
        bslib::nav_panel("Top DE", scroll_box(textOutput("de_hdr"), tableOutput("de_tbl"))),
        bslib::nav_panel("Gene combination",
                         scroll_box(textOutput("combo_hdr"), tableOutput("combo_tbl"))),
        bslib::nav_panel(
          "Marker evidence",
          div(class = "p-2",
              # Cell type by default: immune_infiltrate spans T, NK, myeloid, B and plasma
              # over 74 genes, so a B cell scores 1.00 on it
              radioButtons("mkby", NULL, c("Cell type" = "cell_type",
                                           "Compartment" = "compartment"),
                           "cell_type", inline = TRUE)),
          scroll_box(textOutput("mk_hdr"), tableOutput("mk_tbl")))))))

server <- function(input, output, session) {
  # The axes of a view. The Subclusters tab takes them from the parent picked in its list.
  emb_of <- function(view) {
    view <- view %||% names(EMB)[1]
    if (identical(view, SUB_TAB) && isTruthy(input$sub_parent) && input$sub_parent %in% names(SUB_EMB))
      SUB_EMB[[input$sub_parent]]
    else EMB[[view]]
  }
  # Picking a parent colours by its subclusters and starts unzoomed
  observeEvent(input$sub_parent, {
    b <- sub(" · cluster .*$", "", input$sub_parent)
    if (b %in% names(cells) && !identical(input$cby, b)) updateSelectInput(session, "cby", selected = b)
  }, ignoreInit = TRUE)
  if (length(GENES)) {
    updateSelectizeInput(session, "gene",  choices = c("", GENES), server = TRUE)
    updateSelectizeInput(session, "fgene", choices = c("", GENES), server = TRUE)
    updateSelectizeInput(session, "genes", choices = GENES, server = TRUE)
  }

  # Tick methods; every intersection among them that actually occurs becomes a colour.
  # Two methods gives three: each alone and the pair. Cells no ticked method calls are NA,
  # so they are drawn grey and never reach the legend -- a "not B" entry of 153,137 cells
  # is not an annotation, it is the rest of the object.
  # k methods gives at most 2^k - 1 levels, which is why the control is a tick list: the
  # readable range is yours to choose rather than something the column decides for you.
  lin_vec <- reactive({
    h  <- LIN[[LIN_NM]]
    ms <- LIN_M[LIN_M %in% (input$lin_methods %||% character())]
    if (!length(ms)) return(factor(rep(NA_character_, nrow(h))))
    bit  <- 2^(seq_along(ms) - 1)
    code <- as.vector((h[, ms, drop = FALSE] * 1) %*% bit)   # one integer per cell
    seen <- setdiff(sort(unique(code)), 0)
    if (!length(seen)) return(factor(rep(NA_character_, nrow(h))))
    lbl  <- vapply(seen, function(k) paste(ms[bitwAnd(k, bit) > 0], collapse = "+"), "")
    # singletons first, then pairs, then triples: the agreements sit after the methods
    # they are agreements between, which is the order the legend is read in
    ord  <- order(vapply(seen, function(k) sum(bitwAnd(k, bit) > 0), 1L), seen)
    i <- match(code, seen)
    factor(ifelse(is.na(i), NA_character_, lbl[i]), levels = lbl[ord])
  })

  # Every read of the colour column goes through here: it is a real column of `cells` for
  # every annotation except the lineage one, which exists only for as long as the ticks do.
  # Confidence belongs to the annotation being shown, so it appears with it rather than as a
  # second column picker. Colouring by it is the "conf:" row in the gene list.
  output$conf_ui <- renderUI({
    sc <- score_of(input$cby)
    if (is.null(sc)) return(NULL)
    nm <- sub("^ht$", "HieraType", sub("_score$", "", sub("^insitutype_sup_", "", sc)))
    tagList(
      tags$label(class = "control-label", paste0("Minimum ", nm, " confidence")),
      sliderInput("conf_min", NULL, 0, 1, isolate(input$conf_min %||% 0),
                  step = 0.05, ticks = FALSE))
  })

  output$conf_colour_ui <- renderUI({
    sc <- score_of(input$cby)
    if (is.null(sc)) return(NULL)
    nm <- sub("^ht$", "HieraType", sub("_score$", "", sub("^insitutype_sup_", "", sc)))
    checkboxInput("conf_show", paste0("Colour by ", nm, " confidence"), FALSE)
  })

  cby_all <- reactive({
    req(input$cby)
    if (!is.null(LIN_COL) && identical(input$cby, LIN_COL)) lin_vec()
    else if (!is.null(EDIT_COL) && identical(input$cby, EDIT_COL)) edit_vec()
    else if (input$cby %in% LIVE_COLS) live_vec(input$cby)
    else cells[[input$cby]]
  })

  spatial <- reactive(identical(input$view, "Spatial"))
  vid     <- reactive({ v <- unname(VID[input$view %||% ""]); if (is.na(v)) "spatial" else v })
  zoom    <- reactiveValues(x = NULL, y = NULL)
  sel_ids <- reactiveVal(NULL)   # cell_IDs kept by "Filter cells to selected"
  mode    <- reactive(input$mode %||% "inspect")
  # The last outline drawn: its view, the outline in that view's data coordinates (so it stays
  # put through zooming) and the cells inside it
  lasso   <- reactiveVal(NULL)
  # Each view's latest pixel-to-data mapping, read off its hover events
  hv_map  <- reactiveValues()

  # Manual label edits, one batch per Assign so Undo steps back one action. Every row is kept in
  # order: of a cell's edits the later one wins, by label_state.R's rules
  batches <- reactiveVal(list())
  edits <- reactive({
    b <- batches()
    if (!length(b)) return(empty_overrides())
    do.call(rbind, b)
  })
  # The edits that still count, one or two rows per cell
  edits_now <- reactive(compact_overrides(edits()))
  # Levels start from the cell types and new labels follow, so a label that arrives leaves the
  # others in place
  lab_factor <- function(v) {
    base <- if (is.null(BASE_LAB)) character() else levels(cells[[BASE_LAB]])
    factor(v, levels = c(base, setdiff(sort(unique(v[!is.na(v)])), base)))
  }
  # The edits over the bundle's own cell_type, high-confidence ones only
  edit_vec <- reactive({
    base <- as.character(cells[[BASE_LAB]])
    lab_factor(apply_label_state(cells$cell_ID, cells$cell_ID,
                                 stats::setNames(base, cells$cell_ID), edits())$cell_type)
  })
  # What annotate will write: each cell's cluster label from Tables > Cluster labels, then its
  # Annotate edits
  live <- reactive({
    d <- cl_tab()
    lab <- stats::setNames(d$label, d$cluster)
    lab[!nzchar(lab)] <- NA
    apply_label_state(cells$cell_ID, cells[[input$cl_basis]], lab, edits())
  })
  live_vec <- function(col) {
    l <- live()
    if (identical(col, CONF_COL)) factor(l$label_confidence, levels = c("high", "low"))
    else lab_factor(if (identical(col, DISPLAY_COL)) l$cell_type_display else l$cell_type)
  }
  custom_vec <- reactive(live_vec(CUSTOM_COL))

  # Plots are always width 100%, so the width Shiny reports back is the container width
  base_w <- reactiveVal(900)
  observe({
    w <- session$clientData[[paste0("output_plot_", vid(), "_width")]]
    if (!is.null(w) && w > 200) base_w(round(w / 50) * 50)
  })

  keep_s <- reactive({
    if (is.null(SAMP)) return(rep(TRUE, nrow(cells)))
    req(input$samples)
    k <- cells[[SAMP]] %in% input$samples
    # Spatial shows the one sample picked above the plot, unless All samples is picked
    one <- input$one_sample %||% ALL_SAMPLES
    if (spatial() && one != ALL_SAMPLES) k <- k & cells[[SAMP]] == one
    k
  })
  # Greying 170,000 cells still draws them, and the grey is what you end up looking at.
  # Dropping them is the only way to see the 200 that are left.
  keep <- reactive({
    k <- keep_s()
    cc <- score_of(input$cby)
    if (!is.null(cc) && (input$conf_min %||% 0) > 0)
      k <- k & !is.na(cells[[cc]]) & cells[[cc]] >= input$conf_min
    # A subcluster embedding only holds its own parent's cells; the rest have no coordinates
    ax <- emb_of(input$view)
    if (!is.null(ax) && all(ax %in% names(cells))) k <- k & !is.na(cells[[ax[1]]])
    if (length(input$keep_lv)) k <- k & cby_all() %in% input$keep_lv
    fg <- input$fgene
    if (isTruthy(fg) && !is.null(input$fmin) && !is.na(input$fmin)) {
      m <- if (identical(input$escale, "raw") && !is.null(expr_raw) && fg %in% rownames(expr_raw))
        expr_raw else expr
      if (fg %in% rownames(m)) k <- k & as.numeric(m[fg, ]) > input$fmin
    }
    if (!is.null(sel_ids())) k <- k & cells$cell_ID %in% sel_ids()
    k
  })
  # Colour by confidence: reuses the gene gradient, driven by a switch next to Colour by
  conf_row <- reactive({
    sc <- score_of(input$cby)
    if (is.null(sc)) return(NULL)
    r <- paste0("conf:", sub("_score$", "", sc))
    if (r %in% GENES) r else NULL
  })
  gene  <- reactive({
    if (isTRUE(input$conf_show) && !is.null(conf_row())) return(conf_row())
    if (length(GENES) && isTruthy(input$gene)) input$gene else NULL
  })
  emat  <- reactive({
    m <- if (identical(input$escale, "raw") && !is.null(expr_raw)) expr_raw else expr
    # a gene missing from the raw export would error on `[`; fall back rather than crash
    if (!is.null(gene()) && !gene() %in% rownames(m)) expr else m
  })
  # 2+ genes selected takes precedence over the single-gene view
  gcombo <- reactive({
    g <- intersect(input$genes, rownames(emat()))
    if (length(g) > 1) g else NULL
  })
  gthr <- reactive(if (is.null(input$gthr) || is.na(input$gthr)) 0 else input$gthr)

  # every +/- combination of the picked genes, ordered fewest "+" first
  pat_levels <- reactive({
    g <- gcombo(); req(!is.null(g), length(g) <= 3)
    grid <- expand.grid(rep(list(c("+", "-")), length(g)), stringsAsFactors = FALSE)
    lv <- apply(grid, 1, function(v) paste(paste0(g, v), collapse = " "))
    lv[order(sapply(lv, function(x) lengths(regmatches(x, gregexpr("\\+", x)))))]
  })
  observeEvent(list(gcombo(), input$gmode), {
    if (identical(input$gmode, "pattern") && !is.null(gcombo()) && length(gcombo()) <= 3)
      updateSelectInput(session, "hlpat", choices = c("(show all)", pat_levels()),
                        selected = "(show all)")
  })

  # attach expression as a column rather than re-slicing per consumer: brushedPoints and
  # the composition table then carry it for free, with no realignment step
  shown <- reactive({
    d <- cells[keep(), ]
    if (!is.null(LIN_COL)) d[[LIN_COL]] <- lin_vec()[keep()]
    if (!is.null(EDIT_COL)) d[[EDIT_COL]] <- edit_vec()[keep()]
    # Only while shown, so a Set does not redraw a plot coloured by something else
    if (input$cby %in% LIVE_COLS) d[[input$cby]] <- live_vec(input$cby)[keep()]
    if (!is.null(gcombo())) {
      mm <- as.matrix(emat()[gcombo(), keep(), drop = FALSE] > gthr())
      d$.nhit <- as.integer(colSums(mm))
      # "CD79A+ SDC1-" per cell. Built with vectorised paste, not apply() over 126k columns.
      if (identical(input$gmode, "pattern") && length(gcombo()) <= 3)
        d$.pat <- do.call(paste, lapply(seq_along(gcombo()),
                          function(j) paste0(gcombo()[j], ifelse(mm[j, ], "+", "-"))))
    } else if (!is.null(gene())) {
      d$.expr <- as.numeric(emat()[gene(), keep()])
      # NA rather than a second layer: scale_colour_gradientn paints NA with na.value,
      # so below-threshold cells go grey and the gradient still spans only the kept ones
      if (isTruthy(input$emin) && !is.na(input$emin)) d$.expr[d$.expr <= input$emin] <- NA_real_
    }
    if (stacked()) d <- add_stack(d)
    d
  })
  hl <- reactive(if (length(input$hl_lv)) input$hl_lv else NULL)

  # Tissue is always drawn as ONE panel with the samples stacked vertically, not as facets.
  # facet_wrap + coord_fixed forces every panel to the tallest sample's y range, so a flat
  # tissue sits in 14% of its box; stacking packs them and keeps brush and zoom single-valued.
  stacked <- reactive(spatial() && !is.null(SAMP) && all(c("x_rel", "y_abs") %in% names(cells)))
  axes <- reactive({
    if (stacked()) c(".sx", ".sy") else emb_of(input$view)
  })

  # Band geometry, first sample at the top. Built from the sample selection only, so filters
  # never reflow the layout. Each band has a label row above it: names inside the panel
  # instead of on a y axis, which narrowed the panel and left coord_fixed dead space.
  bands <- reactive({
    req(stacked())
    d  <- cells[keep_s(), ]
    fit <- identical(input$layout, "fit")
    yc <- if (fit) "y_rel" else "y_abs"
    xc <- if (fit) "x_rel" else "x_abs"
    lv <- levels(droplevels(d[[SAMP]]))
    h  <- vapply(lv, function(s) max(d[[yc]][d[[SAMP]] == s], na.rm = TRUE), numeric(1))
    xr   <- range(d[[xc]], na.rm = TRUE)
    pad  <- 0.01 * diff(xr)                                # polygons reach past their centroids
    xlim <- xr + c(-pad, pad)
    lab  <- 30 * fig_s() * diff(xlim) / max(200, base_w() - 24)   # ~30 px label row, in data units
    n   <- length(lv); top <- numeric(n); names(top) <- lv
    if (n > 1) for (i in seq(n - 1L, 1L)) top[i] <- top[i + 1L] + h[i + 1L] + lab
    total <- top[1] + h[1] + lab
    list(lv = lv, h = h, off = top, lab_y = top + h + lab / 2,
         xlim = xlim, ylim = c(-pad, total), yc = yc, xc = xc)
  })

  # Unzoomed extent of the current view, in data units
  full_range <- reactive({
    if (stacked()) return(list(x = bands()$xlim, y = bands()$ylim))
    d <- shown(); xy <- axes()
    list(x = range(d[[xy[1]]], na.rm = TRUE), y = range(d[[xy[2]]], na.rm = TRUE))
  })

  add_stack <- function(d) {
    b <- bands()
    d$.sx <- d[[b$xc]]
    d$.sy <- d[[b$yc]] + b$off[as.character(d[[SAMP]])]
    d
  }
  n_panels <- reactive({
    if (spatial() || !isTruthy(input$facet) || input$facet == "none") 1L
    else max(1L, nlevels(droplevels(shown()[[input$facet]])))
  })
  brush <- reactive(input[[paste0("br_", vid())]])

  observeEvent(input$view, {
    zoom$x <- NULL; zoom$y <- NULL
    ch <- cat_for(input$view)
    updateSelectInput(session, "cby", choices = grouped_choices(ch),
                      selected = if (input$cby %in% ch) input$cby else ch[1])
  }, ignoreInit = TRUE)

  # Level pickers follow Colour by; a level name from another column would filter to nothing
  # A highlight asked for before its Colour by is in place waits here (Cluster labels panel)
  pend_hl <- reactiveVal(NULL)
  observeEvent(cby_all(), {
    lv <- levels(droplevels(cby_all()))
    for (id in c("keep_lv", "hl_lv")) {
      sel <- if (id == "hl_lv" && !is.null(pend_hl())) intersect(pend_hl(), lv)
             else intersect(input[[id]], lv)
      updateSelectizeInput(session, id, choices = lv, selected = sel)
    }
    pend_hl(NULL)
  })

  reset_zoom <- function() {
    zoom$x <- NULL; zoom$y <- NULL
    if (isTruthy(isolate(input$fov))) updateSelectizeInput(session, "fov", selected = "")
  }
  # Changes only when zoom switches on or off, so the plot container is not rebuilt per pan
  zoomed <- reactiveVal(FALSE)
  observe(zoomed(!is.null(zoom$x)))
  # The page only hands scroll gestures to the plot while zoomed
  observe(session$sendCustomMessage("cosmx-zoomed", zoomed()))

  # A zoomed Spatial view is a screen-sized viewport, not the full-height stack: it pans in
  # every direction and each re-render draws ~10x fewer pixels
  # The height of the plot's card body, measured in the page (ZOOM_JS), so a view fits the room
  # the window and the Tables panel leave it
  avail_h  <- reactive(max(220, input$cosmx_avail %||% 600))
  view_h   <- reactive(avail_h())
  legend_px <- reactive(if (is.null(gene()) && is.null(gcombo())) 0 else 90)
  view_asp <- reactive((view_h() - 24 - legend_px()) / max(100, base_w() - 24))

  # Window of the given width centred on (cx, cy), shaped to the viewport
  viewport <- function(cx, cy, w) list(x = cx + c(-1, 1) * w / 2, y = cy + c(-1, 1) * w * view_asp() / 2)

  # Pan by a swipe of (dx, dy) screen px, kept inside the data's extent
  pan_by <- function(p) {
    req(!is.null(zoom$x))
    fr <- full_range()
    keep_in <- function(r, f) {
      if (diff(r) >= diff(f)) return(r)
      r - (min(0, r[1] - f[1]) + max(0, r[2] - f[2]))
    }
    zoom$x <- keep_in(zoom$x + p$dx / max(1, p$w) * diff(zoom$x), fr$x)
    zoom$y <- keep_in(zoom$y - p$dy / max(1, p$h) * diff(zoom$y), fr$y)
  }

  # Zoom to a box. On Spatial it is widened to the viewport's shape
  zoom_to <- function(b) {
    x <- c(b$xmin, b$xmax); y <- c(b$ymin, b$ymax)
    if (spatial()) {
      v <- viewport(mean(x), mean(y), max(diff(x), diff(y) / view_asp()))
      x <- v$x; y <- v$y
    }
    zoom$x <- x; zoom$y <- y
  }
  observeEvent(input$zoom_reset, reset_zoom())
  observeEvent(list(input$samples, input$one_sample, input$layout, input$sub_parent), reset_zoom(),
               ignoreInit = TRUE)

  # Zoom about the cursor. dy > 0 (scroll down, pinch in) zooms out. Zooming out past the
  # full width returns to the overview.
  zoom_by <- function(dy, hv) {
    fr <- full_range()
    f  <- exp(max(-2, min(2, dy * 0.002)))
    if (is.null(zoom$x)) {
      if (f >= 1) return()
      px <- hv$x %||% mean(fr$x); py <- hv$y %||% mean(fr$y)
      start <- if (spatial()) viewport(mean(fr$x), py, diff(fr$x)) else fr
      cx <- start$x; cy <- start$y
    } else {
      cx <- zoom$x; cy <- zoom$y
      px <- hv$x %||% mean(cx); py <- hv$y %||% mean(cy)
    }
    nx <- px + (cx - px) * f; ny <- py + (cy - py) * f
    if (diff(nx) >= diff(fr$x) && (spatial() || diff(ny) >= diff(fr$y))) return(reset_zoom())
    if (diff(nx) < diff(fr$x) / 500) return()
    zoom$x <- nx; zoom$y <- ny
  }

  # Jump to an FOV: its cells' extent in stacked coordinates, with a small margin
  observeEvent(input$fov, {
    if (!nzchar(input$fov)) return(reset_zoom())
    req(stacked())
    s <- sub("\\|.*$", "", input$fov); f <- as.integer(sub("^.*\\|", "", input$fov))
    b <- bands(); req(s %in% b$lv)
    i <- which(cells[[SAMP]] == s & cells$fov == f)
    x <- range(cells[[b$xc]][i]); y <- range(cells[[b$yc]][i]) + b$off[[s]]
    m <- 0.03 * max(diff(x), diff(y))
    zoom_to(list(xmin = x[1] - m, xmax = x[2] + m, ymin = y[1] - m, ymax = y[2] + m))
  }, ignoreInit = TRUE)

  # Clicked cell: inspected in the side panel and outlined on the plot
  sel_cell <- reactiveVal(NULL)

  cell_info_ui <- function(id) {
    cid <- sel_cell(); req(cid)
    i <- match(cid, cells$cell_ID); req(!is.na(i))
    r <- cells[i, , drop = FALSE]
    val <- function(col) if (col %in% names(r)) as.character(r[[col]]) else
      if (col %in% LIVE_COLS) as.character(live_vec(col)[i]) else NA_character_
    meta <- c(Sample = if (!is.null(SAMP)) val(SAMP), FOV = val("fov"),
              Patient = val("patient_id"), Treatment = val("treatment"),
              Timepoint = val("timepoint"),
              `Total counts` = if ("total_counts" %in% names(r)) format(r$total_counts, big.mark = ","))
    ann <- unique(c(input$cby, grep("^label_", names(cells), value = TRUE),
                    intersect(c("insitutype_unsup", "ht_call"), names(cells)),
                    grep("^insitutype_sup_", names(cells), value = TRUE) |>
                      grep(pattern = "_(anchor|second_type|score|conf|refined)$", invert = TRUE,
                           value = TRUE)))
    ann <- setdiff(ann, LIN_COL)
    rows <- c(meta, stats::setNames(vapply(ann, val, ""), ann))
    if (!is.null(EDIT_COL)) {
      e <- edits_now()[edits_now()$cell_ID == cid, , drop = FALSE]
      if (nrow(e)) rows <- c(rows, `Edited label` = paste0(e$label, " (", e$confidence, ")",
                                                           collapse = ", "))
    }
    rows <- rows[!is.na(rows) & nzchar(rows)]
    top <- if (!is.null(expr_raw)) {
      # Genes only: the IF:/conf:/meta: rows are intensities and scores, not counts
      v <- expr_raw[!grepl(":", rownames(expr_raw)), i]; v <- sort(v[v > 0], decreasing = TRUE)
      if (length(v)) paste(sprintf("%s %g", names(v)[seq_len(min(10, length(v)))],
                                   v[seq_len(min(10, length(v)))]), collapse = ", ")
    }
    div(class = "border rounded p-2 mb-3", style = "font-size: .85rem;",
        div(class = "d-flex justify-content-between",
            strong(cid), actionLink(paste0("cell_close_", id), "×")),
        tags$table(class = "table table-sm mb-1",
                   lapply(names(rows), function(k) tags$tr(tags$td(class = "text-muted", k),
                                                           tags$td(rows[[k]])))),
        if (!is.null(top)) div(span(class = "text-muted", "Top genes (raw counts): "), top))
  }

  if (!is.null(SAMP)) {
    observeEvent(input$samples_all,
                 updateCheckboxGroupInput(session, "samples", selected = levels(cells[[SAMP]])))
    observeEvent(input$samples_none,
                 updateCheckboxGroupInput(session, "samples", selected = character()))
    # The sample Spatial shows is one of those ticked under Samples
    observeEvent(input$samples, {
      cur <- isolate(input$one_sample) %||% ALL_SAMPLES
      ch <- c(stats::setNames(input$samples, input$samples), "All samples" = ALL_SAMPLES)
      updateSelectInput(session, "one_sample", choices = ch,
                        selected = if (cur %in% ch) cur else ch[[1]])
    }, ignoreInit = TRUE, ignoreNULL = FALSE)
  }
  # Jump to FOV offers the FOVs of the sample on screen
  if (!is.null(FOVS)) observeEvent(input$one_sample, {
    one <- input$one_sample
    updateSelectizeInput(session, "fov", selected = "",
                         choices = c("Jump to FOV" = "", if (one == ALL_SAMPLES) FOVS else FOVS[one]))
  })

  observeEvent(input$filter_reset, {
    updateSelectizeInput(session, "keep_lv", selected = character())
    if (length(GENES)) updateSelectizeInput(session, "fgene", choices = c("", GENES),
                                            selected = "", server = TRUE)
    updateSliderInput(session, "conf_min", value = 0)
    sel_ids(NULL)
  })
  observeEvent(input$sel_clear, sel_ids(NULL))

  output$sel_status <- renderUI({
    req(!is.null(sel_ids()))
    div(class = "d-flex justify-content-between align-items-center mb-2",
        span(sprintf("Selection: %s cells", format(length(sel_ids()), big.mark = ","))),
        actionLink("sel_clear", "Clear"))
  })

  output$filter_summary <- renderText({
    on <- c(if (length(input$keep_lv)) "levels",
            if (isTruthy(input$fgene)) paste0(input$fgene, " > ", input$fmin),
            if (!is.null(score_of(input$cby)) && (input$conf_min %||% 0) > 0) "confidence",
            if (!is.null(sel_ids())) "selection")
    sprintf("%s of %s cells%s", format(sum(keep()), big.mark = ","),
            format(sum(keep_s()), big.mark = ","),
            if (length(on)) paste0("  ·  filtered by ", paste(on, collapse = ", ")) else "")
  })

  pal_hl <- reactive({
    f <- cby_all()
    p <- pal_for(f)
    # Subclusters take their colours from the ones in view: one parent's, or those kept by the
    # level filter. Spread over all 134 of sub_res0.3 the hues sit ~3 degrees apart and a
    # parent's dozen came out as one colour. No parent has more than 24.
    if (input$cby %in% SUB_COLS) {
      k  <- !is.na(f)
      ax <- emb_of(input$view)
      if (!is.null(ax) && all(ax %in% names(cells))) k <- k & !is.na(cells[[ax[1]]])
      if (length(input$keep_lv)) k <- k & f %in% input$keep_lv
      lv <- levels(droplevels(f[k]))
      if (length(lv) <= length(pal_wide)) p[lv] <- pal_wide[seq_along(lv)]
    }
    # Labels arrive one Set at a time. pal_n spreads its hues over however many levels there
    # are, so each new label would recolour the rest; a colour per label does not.
    if (input$cby %in% LAB_COLS) p <- pal_labels(levels(f))
    if (!is.null(hl())) p[!names(p) %in% hl()] <- GREY
    p
  })

  legend_ui <- function() {
    req(input$cby)
    f <- droplevels(factor(shown()[[input$cby]]))
    n <- table(f); p <- pal_hl()
    low <- low_cells(input$cby)
    n_low <- if (is.null(low)) 0L else sum(low[match(shown()$cell_ID, cells$cell_ID)], na.rm = TRUE)
    tagList(
      lapply(levels(f), function(l) div(
        class = "legend-row",
        span(class = "legend-sw", style = paste0("background:", p[[l]], ";")), l,
        span(style = "color:#8a8a8a;font-size:.85em", HTML("&nbsp;"),
             format(as.integer(n[[l]]), big.mark = ",")))),
      if (n_low > 0) div(
        class = "legend-row mt-2",
        span(class = "legend-sw legend-low"), "Low confidence (white slash)",
        span(style = "color:#8a8a8a;font-size:.85em", HTML("&nbsp;"),
             format(n_low, big.mark = ","))))
  }

  # Cells whose label in this colouring is a low-confidence guess, which the figures stripe and
  # the views here slash: only for the display labels, whose guesses are what is drawn
  low_cells <- function(col) {
    if (!is.null(DISPLAY_COL) && identical(col, DISPLAY_COL)) return(live()$label_confidence == "low")
    if (identical(col, "cell_type_display")) {
      cc <- intersect(c("label_confidence", "cell_type_confidence"), names(cells))[1]
      if (!is.na(cc)) return(as.character(cells[[cc]]) == "low")
    }
    NULL
  }

  size_ui <- function(id) {
    # Draw mode hands every drag to the outline script, so it has no box, click or double-click.
    # In Inspect the box zooms as soon as the drag pauses, so it waits longer than in Select.
    md <- mode()
    out <- function(h) tagAppendAttributes(
      plotOutput(paste0("plot_", id), height = h,
                 brush = if (md != "draw")
                   brushOpts(paste0("br_", id), resetOnNew = TRUE,
                             delay = if (md == "inspect") 700 else 300),
                 click = if (md != "draw") paste0("ck_", id),
                 dblclick = if (md != "draw") paste0("dbl_", id),
                 hover = hoverOpts(paste0("hv_", id), delay = 60, delayType = "throttle",
                                   nullOutside = TRUE)),
      class = paste("cosmx-zoom", if (md == "draw") "cosmx-draw"))
    if (!spatial()) {
      if (n_panels() == 1) return(out(paste0(avail_h(), "px")))
      return(out(paste0(min(20000, ceiling(n_panels() / 2) * round(base_w() / 2 * 0.8) + 40), "px")))
    }
    if (zoomed()) return(out(paste0(view_h(), "px")))
    # coord_fixed() pins the aspect, so the height is exactly what the panel needs at this
    # width: 24 px of plot margin, plus the colour bar when a gene view puts one on top
    fr  <- full_range()
    asp <- diff(fr$y) / max(1e-9, diff(fr$x))
    leg <- legend_px()
    # One sample fits the window; All samples is a tall stack the page scrolls through
    cap <- if (identical(input$one_sample %||% ALL_SAMPLES, ALL_SAMPLES)) 20000 else avail_h()
    out(paste0(max(220, min(cap, round((base_w() - 24) * asp + 24 + leg))), "px"))
  }

  # geom_point below RASTER_ABOVE cells, scattermore above it. Measured on this bundle,
  # geom_point beat scattermore at every canvas size tried, by 5-16x:
  #   126k cells, 1188x792 px   scatter 2.35 s   point 0.34 s
  #   126k cells, 1320x3300 px  scatter 5.82 s   point 0.48 s
  # scattermore is kept only as a guard for datasets far larger than this one.
  RASTER_ABOVE <- 400000L
  PT <- 1   # point size; points are only used for the embeddings
  # Figure size, from the page's Cmd -/+ and its buttons: scales the plot's text, points and lines
  fig_s <- reactive(input$fig_scale %||% 1)

  pts <- function(n, size) {
    if (n > RASTER_ABOVE) geom_scattermore(pointsize = max(1, size * 6), pixels = c(1000, 1000))
    else geom_point(size = size, shape = 16)
  }

  # Cells inside the zoom window. Drawing only these is what makes zoomed-in polygons fast
  vis <- reactive({
    d <- shown()
    if (is.null(zoom$x)) return(d)
    xy <- axes()
    # A small margin keeps cells whose outline crosses the edge of the view
    mx <- 0.03 * diff(zoom$x); my <- 0.03 * diff(zoom$y)
    d[d[[xy[1]]] >= zoom$x[1] - mx & d[[xy[1]]] <= zoom$x[2] + mx &
      d[[xy[2]]] >= zoom$y[1] - my & d[[xy[2]]] <= zoom$y[2] + my, , drop = FALSE]
  })

  # Tissue is always drawn as segmented cells; ~2.3 s for all 2.7M vertices at full zoom-out
  use_poly <- reactive(!is.null(poly) && spatial())

  # Hoisted out of renderPlot: a 2.7M-row index per frame is the most expensive step here
  poly_sub <- reactive({
    req(use_poly())
    ids <- vis()$cell_ID
    pd  <- poly[cells$cell_ID[poly$i] %in% ids, , drop = FALSE]
    if (stacked()) {                          # same transform + band offset as the points
      r <- attr(cells, "rel"); b <- bands()
      den <- if (b$xc == "x_rel") r$sx[pd$i] else 1
      pd$px <- (pd$px - r$x0[pd$i]) / den
      pd$py <- (pd$py - r$y0[pd$i]) / den + b$off[as.character(cells[[SAMP]][pd$i])]
    }
    pd
  })

  build_plot <- function(id) {
    xy <- axes()
    s <- fig_s()

    # Each mode below picks the cells in draw order (d), the column that colours them (col)
    # and a scale; the geom is chosen once, after. Scales serve colour AND fill, so polygons
    # follow the gene views too -- they used to fill by Colour by whatever else was set.
    both <- c("colour", "fill")
    low <- NULL
    if (!is.null(gcombo()) && identical(input$gmode, "pattern") &&
        ".pat" %in% names(vis())) {
      lv <- pat_levels()
      d  <- vis()
      d$.pat <- factor(d$.pat, levels = lv)
      nplus <- sapply(lv, function(x) lengths(regmatches(x, gregexpr("\\+", x))))
      hp <- input$hlpat %||% "(show all)"
      cols <- if (hp %in% lv) {
        # one combination in red, everything else grey -- the "show me CD79A+/SDC1-" view
        stats::setNames(ifelse(lv == hp, "#D7263D", "#E4E9E7"), lv)
      } else {
        # all-negative stays background grey; the rest get separable hues
        stats::setNames(c("#E4E9E7",
                          grDevices::colorRampPalette(
                            c("#4C9BE8", "#00A870", "#F2A93B", "#D7263D"))(length(lv) - 1)
                          )[rank(nplus, ties.method = "first")], lv)
      }
      # draw the rarer combinations last so they are not buried
      d <- d[order(-table(d$.pat)[as.character(d$.pat)]), ]
      cnt  <- table(factor(d$.pat, levels = lv))
      labs <- stats::setNames(sprintf("%s  (%s)", lv, format(as.integer(cnt[lv]), big.mark = ",")), lv)
      col <- ".pat"
      sc  <- scale_colour_manual(values = cols, labels = labs, drop = FALSE, aesthetics = both,
                                 limits = names(cols),
                                 name = paste0("> ", gthr(),
                                               if (identical(input$escale, "raw")) " (raw counts)" else ""))

    } else if (!is.null(gcombo())) {
      n  <- length(gcombo())
      d  <- vis()[order(vis()$.nhit), ]          # co-expressing cells drawn on top
      lv <- 0:n
      d$.combo <- factor(switch(input$gmode %||% "count",
                                all   = ifelse(d$.nhit == n, "all", "not all"),
                                any   = ifelse(d$.nhit >= 1, "any",  "none"),
                                as.character(d$.nhit)),
                         levels = switch(input$gmode %||% "count",
                                         all = c("not all", "all"),
                                         any = c("none", "any"),
                                         as.character(lv)))
      # 0 gets a near-background grey so it cannot be mistaken for expression; the ramp
      # starts at the first real colour. Without this, "1 of 3" reads as bright green and
      # looks like a hit.
      cols <- switch(input$gmode %||% "count",
                     all = c("not all" = "#E4E9E7", "all" = "#04251C"),
                     any = c("none"    = "#E4E9E7", "any" = "#00A870"),
                     stats::setNames(c("#E4E9E7", grDevices::colorRampPalette(pal_seq[-1])(n)),
                                     as.character(lv)))
      # counts in the legend: "how many cells is this colour" is the first question asked
      # of this plot, and reading it off the picture is exactly how you get it wrong
      cnt  <- table(factor(d$.combo, levels = names(cols)))
      labs <- stats::setNames(sprintf("%s  (%s)", names(cols),
                                      format(as.integer(cnt[names(cols)]), big.mark = ",")),
                              names(cols))
      col <- ".combo"
      sc  <- scale_colour_manual(values = cols, labels = labs, drop = FALSE, aesthetics = both,
                                 limits = names(cols),
                                 name = paste0(paste(gcombo(), collapse = " + "), "\n> ", gthr(),
                                               if (identical(input$escale, "raw")) " (raw counts)" else ""))

    } else if (!is.null(gene())) {
      d <- vis()[order(vis()$.expr, na.last = FALSE), ]   # expressing cells on top
      npass <- sum(!is.na(d$.expr))
      ttl <- paste0(gene(),
                    if (isTruthy(input$emin) && !is.na(input$emin))
                      paste0(" > ", input$emin) else "",
                    if (identical(input$escale, "raw")) "\n(raw counts)" else "",
                    if (isTruthy(input$emin) && !is.na(input$emin))
                      paste0("\n", format(npass, big.mark = ","), " cells") else "")
      col <- ".expr"
      sc  <- scale_colour_gradientn(colours = pal_seq, na.value = "#E4E9E7", name = ttl,
                                    aesthetics = both)

    } else {
      d <- vis()
      if (input$cby %in% LAB_COLS) d <- d[order(lab_layer(d[[input$cby]])), ]
      if (!is.null(hl())) d <- d[order(d[[input$cby]] %in% hl()), ]   # highlighted on top
      col <- input$cby
      sc  <- scale_colour_manual(values = pal_hl(), limits = names(pal_hl()), drop = FALSE, na.value = GREY,
                                 name = input$cby, aesthetics = both)
      # One flag per row of d: its label is a low-confidence guess
      low <- low_cells(col)
      if (!is.null(low)) low <- low[match(d$cell_ID, cells$cell_ID)] %in% TRUE
    }

    p <- if (use_poly()) {
      pd <- poly_sub()
      # group = the cell's row in d: GeomPolygon draws in group order, so expressing and
      # highlighted cells land on top exactly as they do as points
      pd$.g <- match(pd$i, match(d$cell_ID, cells$cell_ID))
      pd <- pd[!is.na(pd$.g), , drop = FALSE]
      pd$.v <- d[[col]][pd$.g]
      p <- ggplot(pd, aes(px, py, group = .g, fill = .v)) + geom_polygon(linewidth = 0) + sc
      # A white slash across each low-confidence cell, corner to corner inside its outline's box
      if (any(low)) {
        q <- pd[low[pd$.g], , drop = FALSE]
        x0 <- tapply(q$px, q$.g, min); x1 <- tapply(q$px, q$.g, max)
        y0 <- tapply(q$py, q$.g, min); y1 <- tapply(q$py, q$.g, max)
        sl <- data.frame(x = x0 + 0.2 * (x1 - x0), xend = x1 - 0.2 * (x1 - x0),
                         y = y0 + 0.2 * (y1 - y0), yend = y1 - 0.2 * (y1 - y0))
        p <- p + geom_segment(data = sl, aes(x = x, y = y, xend = xend, yend = yend),
                              inherit.aes = FALSE, colour = "#FFFFFF", linewidth = 0.35 * s)
      }
      hit <- pd[cells$cell_ID[pd$i] %in% sel_cell(), , drop = FALSE]
      if (nrow(hit)) p <- p + geom_polygon(data = hit, aes(px, py, group = .g), inherit.aes = FALSE,
                                           fill = NA, colour = "#000000", linewidth = 0.8 * s)
      p
    } else {
      p <- ggplot(d, aes(.data[[xy[1]]], .data[[xy[2]]], colour = .data[[col]])) +
        pts(nrow(d), PT * s) + sc
      if (col %in% c(".pat", ".combo"))
        p <- p + guides(colour = guide_legend(override.aes = list(size = 2.5 * s, alpha = 1)))
      # A white slash across each low-confidence dot, as on the pipeline's UMAP
      if (any(low))
        p <- p + layer(geom = GeomSlash, stat = "identity", position = "identity",
                       data = d[low, , drop = FALSE], mapping = aes(.data[[xy[1]]], .data[[xy[2]]]),
                       inherit.aes = FALSE, show.legend = FALSE,
                       params = list(half = 0.85 * dot_r(PT * s) / sqrt(2), linewidth = 0.2 * s,
                                     na.rm = TRUE))
      hit <- d[d$cell_ID %in% sel_cell(), , drop = FALSE]
      if (nrow(hit)) p <- p + geom_point(data = hit, aes(.data[[xy[1]]], .data[[xy[2]]]),
                                         inherit.aes = FALSE, shape = 21, size = 4 * s, stroke = 1.2 * s,
                                         colour = "#000000", fill = NA)
      p
    }

    lz <- lasso()
    if (!is.null(lz) && lz$id == id)
      p <- p + geom_path(data = lz$poly, aes(x, y), inherit.aes = FALSE,
                         colour = "#000000", linewidth = 0.6 * s, linetype = "22")

    # Figure size: text and grid lines scale with the points
    p <- p + theme(text = element_text(size = PLOT_PT * s),
                   line = element_line(linewidth = PLOT_PT / 22 * s))

    legend_on <- !is.null(gene()) || !is.null(gcombo())
    # Spatial is real geometry: aspect pinned, no axes, sample names inside the panel
    if (spatial()) {
      fr <- full_range()
      p <- p + coord_fixed(xlim = zoom$x %||% fr$x, ylim = zoom$y %||% fr$y, expand = FALSE) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              axis.title = element_blank(), panel.grid = element_blank(),
              plot.margin = margin(6, 6, 6, 6),
              legend.position = if (legend_on) "top" else "none")
      if (stacked()) {
        b <- bands()
        p <- p + annotate("text", x = b$xlim[1], y = b$lab_y, label = b$lv,
                          hjust = 0, vjust = 0.5, size = 3.6 * s, colour = "#4A4A4A")
      }
      return(p)
    }
    # Embedding axes have no metric meaning, so they stretch to fill the panel
    if (input$facet != "none")
      p <- p + facet_wrap(vars(.data[[input$facet]]), ncol = min(2L, n_panels()))
    # Zooming changes tick label widths ("-20" -> "-2.5"); padding pins the panel edge
    padlab <- function(v) ifelse(is.na(v), "", formatC(v, width = 6, format = "g"))
    p + scale_x_continuous(labels = padlab) + scale_y_continuous(labels = padlab) +
      coord_cartesian(xlim = zoom$x, ylim = zoom$y) +
      theme(legend.position = if (legend_on) "right" else "none") +
      labs(x = xy[1], y = xy[2])
  }

  # The cells an Assign or a Filter acts on: the outline drawn on this view, else its box
  pending_ids <- reactive({
    lz <- lasso()
    if (!is.null(lz) && lz$id == vid()) return(lz$ids)
    b <- brush()
    if (!is.null(b) && mode() == "select") brushedPoints(vis(), b, axes()[1], axes()[2])$cell_ID
  })
  clear_pending <- function() {
    lasso(NULL)
    session$resetBrush(paste0("br_", vid()))
  }
  observeEvent(input$mode, clear_pending(), ignoreInit = TRUE)
  # An outline is in the coordinates it was drawn on: a new parent, sample set or layout moves them
  observeEvent(list(input$sub_parent, input$samples, input$one_sample, input$layout), lasso(NULL),
               ignoreInit = TRUE)

  assign_now <- function() {
    ids <- pending_ids(); lab <- trimws(input$ann_label %||% "")
    if (!length(ids))
      return(showNotification("Select or draw around some cells first.", type = "warning"))
    if (!nzchar(lab))
      return(showNotification("Pick or type a label in Annotate first.", type = "warning"))
    conf <- input$ann_conf %||% "high"
    batches(c(batches(), list(data.frame(cell_ID = ids, label = lab, confidence = conf,
                                         note = trimws(input$ann_note %||% "")))))
    clear_pending()
    # A low-confidence label shows in the display column, not in cell_type
    show <- if (is.null(CUSTOM_COL)) EDIT_COL else if (conf == "low") DISPLAY_COL else CUSTOM_COL
    if (!input$cby %in% c(EDIT_COL, LIVE_COLS)) updateSelectInput(session, "cby", selected = show)
    showNotification(sprintf("%s cells set to %s (%s confidence)",
                             format(length(ids), big.mark = ","), lab, conf), duration = 4)
  }

  if (!is.null(EDIT_COL)) {
    observeEvent(input$ann_assign, assign_now())
    observeEvent(input$ann_undo, {
      b <- batches()
      if (length(b)) batches(b[-length(b)])
    })
    observeEvent(input$ann_clear, showModal(modalDialog(
      paste("Remove every edited label? This empties", OV_FILE, "too."),
      footer = tagList(modalButton("Cancel"),
                       actionButton("ann_clear_ok", "Remove all", class = "btn-danger")))))
    observeEvent(input$ann_clear_ok, { batches(list()); removeModal() })
    observeEvent(input$ann_up, {
      d <- tryCatch(read_overrides(input$ann_up$datapath), error = function(e) {
        showNotification(conditionMessage(e), type = "error"); NULL })
      req(d)
      batches(c(batches(), list(d)))
      showNotification(sprintf("Loaded %s edited cells", format(nrow(d), big.mark = ",")))
    })
    output$ann_dl <- downloadHandler(
      filename = function() paste0("label_overrides_", format(Sys.Date()), ".csv"),
      content = function(f) utils::write.csv(edits_now(), f, row.names = FALSE))
    output$ann_summary <- renderUI({
      e <- edits_now()
      if (!nrow(e)) return(div(class = "text-muted small mb-2", "No edits yet."))
      n <- table(factor(e$label), factor(e$confidence, levels = c("high", "low")))
      n <- n[order(-rowSums(n)), , drop = FALSE]
      div(class = "small mb-2",
          strong(format(length(unique(e$cell_ID)), big.mark = ","), " cells edited"),
          tags$table(class = "table table-sm mb-0",
                     tags$tr(tags$th(""), tags$th("high"), tags$th("low")),
                     lapply(rownames(n), function(l)
                       tags$tr(tags$td(l), lapply(n[l, ], function(k)
                         tags$td(format(as.integer(k), big.mark = ",")))))))
    })
    # New labels join the picker so they can be reused
    observeEvent(edits(), updateSelectizeInput(
      session, "ann_label",
      choices = c("", union(levels(edit_vec()), sort(setdiff(edits()$label, KEEP_LAB))), KEEP_LAB),
      selected = isolate(input$ann_label)), ignoreInit = TRUE)

    # The run's edits: read at start, saved as they are made. Cut to the rows that count, so a
    # cell relabelled often stays one or two rows
    if (file.exists(OV_FILE)) {
      d0 <- tryCatch(read_overrides(OV_FILE), error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
      if (!is.null(d0) && nrow(d0)) {
        batches(list(d0))
        showNotification(sprintf("Loaded %s label edits from %s",
                                 format(nrow(d0), big.mark = ","), OV_FILE), duration = 8)
      }
    }
    observeEvent(edits(), {
      ok <- tryCatch({
        dir.create(LABEL_DIR, recursive = TRUE, showWarnings = FALSE)
        write_csv_atomic(edits_now(), OV_FILE)
        TRUE
      }, error = function(e) FALSE)
      if (!ok) showNotification(paste("Could not save the edits to", OV_FILE), type = "error")
    }, ignoreInit = TRUE)
  }

  for (v in unname(VID)) local({
    id <- v
    br <- paste0("br_", id)
    output[[paste0("legend_", id)]] <- renderUI(legend_ui())
    output[[paste0("cellinfo_", id)]] <- renderUI(cell_info_ui(id))

    # Click a cell to inspect it; clicking empty space clears the panel
    observeEvent(input[[paste0("ck_", id)]], {
      hit <- nearPoints(vis(), input[[paste0("ck_", id)]], axes()[1], axes()[2],
                        threshold = 40, maxpoints = 1)
      sel_cell(if (nrow(hit)) hit$cell_ID[1] else NULL)
    })
    observeEvent(input[[paste0("cell_close_", id)]], sel_cell(NULL))
    output[[paste0("ui_", id)]]     <- renderUI(size_ui(id))
    output[[paste0("plot_", id)]]   <- renderPlot(build_plot(id), res = 112)

    # A box in Select mode or an outline in Draw mode pops up the selection bar. Filtering
    # narrows every view to those cells; Assign gives them the label chosen in Annotate.
    output[[paste0("selbar_", id)]] <- renderUI({
      lz <- lasso()
      drawn <- !is.null(lz) && lz$id == id
      n <- if (drawn) length(lz$ids) else {
        b <- input[[br]]; req(b, mode() == "select")
        nrow(brushedPoints(vis(), b, axes()[1], axes()[2]))
      }
      div(class = "selbar",
          strong(format(n, big.mark = ","), " cells selected"),
          actionButton(paste0("sel_apply_", id), "Filter to these", class = "btn-sm btn-primary"),
          if (!is.null(EDIT_COL))
            actionButton(paste0("sel_assign_", id), "Assign label", class = "btn-sm btn-success"),
          actionButton(paste0("sel_cancel_", id), "Clear", class = "btn-sm btn-outline-secondary"),
          span(class = "text-muted small",
               if (drawn) "draw again to replace" else "drag again to replace"))
    })
    observeEvent(input[[paste0("sel_apply_", id)]], {
      ids <- pending_ids(); req(length(ids))
      sel_ids(if (is.null(sel_ids())) ids else intersect(sel_ids(), ids))
      clear_pending()
    })
    observeEvent(input[[paste0("sel_assign_", id)]], assign_now())
    observeEvent(input[[paste0("sel_cancel_", id)]], clear_pending())

    # Inspect mode: a box zooms straight away
    observeEvent(input[[br]], {
      b <- input[[br]]
      if (!is.null(b) && mode() == "inspect") {
        zoom_to(b)
        session$resetBrush(br)
      }
    })
    # Double-click inside a Select box zooms to it; anywhere else resets
    observeEvent(input[[paste0("dbl_", id)]], {
      b <- input[[br]]
      if (is.null(b)) return(reset_zoom())
      zoom_to(b)
      session$resetBrush(br)
    })

    # Hover events carry the plot's coordmap. The latest one maps a drawn outline to data units.
    observeEvent(input[[paste0("hv_", id)]], {
      h <- input[[paste0("hv_", id)]]
      if (!is.null(h$domain) && !is.null(h$range) && !is.null(h$img_css_ratio))
        hv_map[[id]] <- h[c("domain", "range", "img_css_ratio")]
    })
    observeEvent(input[[paste0("plot_", id, "_lasso")]], {
      m <- hv_map[[id]]
      if (is.null(m))
        return(showNotification("Move the pointer over the plot once, then draw again.", type = "warning"))
      if (n_panels() > 1)
        return(showNotification("Drawing needs one panel: set Split by to none.", type = "warning"))
      xy <- matrix(unlist(input[[paste0("plot_", id, "_lasso")]]$pts), ncol = 2, byrow = TRUE)
      if (nrow(xy) > 400) xy <- xy[unique(round(seq(1, nrow(xy), length.out = 400))), , drop = FALSE]
      dm <- m$domain; rg <- m$range
      px <- xy[, 1] * m$img_css_ratio$x; py <- xy[, 2] * m$img_css_ratio$y
      x <- dm$left + (px - rg$left) / (rg$right - rg$left) * (dm$right - dm$left)
      y <- dm$bottom + (rg$bottom - py) / (rg$bottom - rg$top) * (dm$top - dm$bottom)
      v <- vis(); ax <- axes()
      inside <- in_poly(v[[ax[1]]], v[[ax[2]]], x, y)
      lasso(list(id = id, poly = data.frame(x = c(x, x[1]), y = c(y, y[1])),
                 ids = v$cell_ID[inside]))
    })
    observeEvent(input[[paste0("plot_", id, "_wheel")]], {
      zoom_by(input[[paste0("plot_", id, "_wheel")]]$dy, input[[paste0("hv_", id)]])
    })
    observeEvent(input[[paste0("plot_", id, "_pan")]], pan_by(input[[paste0("plot_", id, "_pan")]]))
  })

  # --- marker evidence ------------------------------------------------------------------
  # Fold enrichment of each marker panel, per level of Colour by, for the CURRENT selection.
  # The background is the sample filter only, never the legend filter: scored against itself
  # a filtered selection is enriched 1.0 and the table says nothing.
  # One sparse pass. obs = cells in level L detecting gene g, summed over the panel;
  # exp = the same count predicted from background rates within depth deciles. Aggregating
  # counts before dividing, rather than averaging per-gene folds, keeps a gene that is near
  # zero in the background from dominating the panel.
  mk_eval <- reactive({
    req(length(MK_SET), !is.null(expr_raw), input$cby)
    pan <- MK_SET[[input$mkby %||% "cell_type"]]
    req(length(pan))
    d <- shown(); req(nrow(d) > 0)
    f <- droplevels(factor(d[[input$cby]]))
    # the lineage column is NA for cells no ticked method calls; they are canvas, not a
    # group, and fac2sparse has no column to put them in
    ok <- !is.na(f); d <- d[ok, , drop = FALSE]; f <- droplevels(f[ok])
    req(nrow(d) > 0, nlevels(f) > 0)
    i  <- match(d$cell_ID, cells$cell_ID)
    bg <- keep_s()
    dec <- if (is.null(DEC)) rep(1L, nrow(cells)) else DEC

    H  <- (expr_raw[MK_GENES, i, drop = FALSE] > 0) * 1        # genes x selected
    Lm <- Matrix::fac2sparse(f)                                 # levels x selected
    OBS <- as.matrix(Matrix::tcrossprod(H, Lm))                 # genes x levels
    NQ  <- as.matrix(Matrix::tcrossprod(
      Matrix::fac2sparse(factor(dec[i], levels = seq_len(NDEC))), Lm))   # deciles x levels
    bgd <- vapply(seq_len(NDEC), function(q) {
      k <- bg & dec == q
      if (!any(k)) rep(0, length(MK_GENES))
      else Matrix::rowMeans(expr_raw[MK_GENES, k, drop = FALSE] > 0)
    }, numeric(length(MK_GENES)))
    EXP <- bgd %*% NQ                                           # genes x levels

    fold <- vapply(pan, function(gs) {
      k <- match(gs, MK_GENES)
      o <- colSums(OBS[k, , drop = FALSE]); e <- colSums(EXP[k, , drop = FALSE])
      ifelse(e > 0, o / e, NA_real_)
    }, numeric(nlevels(f)))
    dim(fold) <- c(nlevels(f), length(pan))
    out <- data.frame(level = levels(f), n_cells = as.integer(table(f)),
                      round(as.data.frame(fold), 1), check.names = FALSE)
    names(out)[-(1:2)] <- sprintf("%s (%d)", names(pan), lengths(pan))
    out[order(-out$n_cells), ]
  })
  output$mk_tbl <- renderTable(mk_eval(), digits = 1)
  output$mk_hdr <- renderText({
    if (!length(MK_SET)) return("assets/canonical_markers.csv not found beside the app.")
    pan <- MK_SET[[input$mkby %||% "cell_type"]]
    sprintf("Fold over %s background, depth-matched on %d count deciles. Each column is a marker panel, its number of genes in brackets.",
            if (is.null(DEC)) "NAIVE (no total_counts)" else "sample-filtered", NDEC)
  })

  picked <- reactive({
    xy <- axes()
    if (is.null(brush())) vis() else brushedPoints(vis(), brush(), xy[1], xy[2])
  })

  output$comp_hdr <- renderText({
    sprintf("%s of %s cells%s%s", format(nrow(picked()), big.mark = ","),
            format(nrow(vis()), big.mark = ","),
            if (is.null(brush())) "" else " (brushed)",
            if (is.null(zoom$x)) "  --  Inspect mode: drag a box to zoom (or pinch, Cmd/Ctrl + scroll)"
            else "  --  zoomed; double-click to reset")
  })

  output$comp <- renderTable({
    d <- picked(); req(nrow(d) > 0)
    t <- sort(table(droplevels(d[[input$cby]])), decreasing = TRUE)
    out <- data.frame(level = names(t), n = as.integer(t),
                      pct = sprintf("%.1f%%", 100 * as.integer(t) / nrow(d)))
    if (!is.null(gene()) && ".expr" %in% names(d))
      out$mean_expr <- round(tapply(d$.expr, droplevels(d[[input$cby]]), mean)[out$level], 3)
    out
  })

  # Poisson-binomial: exact P(exactly k of n independent events), by convolution.
  # n is the number of picked genes, so this is trivially small.
  pb_dist <- function(p) { d <- 1; for (q in p) d <- c(d, 0) * (1 - q) + c(0, d) * q; d }

  output$combo_hdr <- renderText({
    if (is.null(gcombo())) "Pick 2 or more genes under Highlight."
    else sprintf("%s  (>%s, %s)  --  %s cells%s", paste(gcombo(), collapse = " + "),
                 format(gthr()), if (identical(input$escale, "raw")) "raw counts" else "normalized",
                 format(nrow(picked()), big.mark = ","),
                 if (is.null(brush())) "" else " (brushed)")
  })

  # The whole point of this panel. "N cells express all three" is meaningless on its own:
  # with markers detected in a few percent of cells, chance co-detection is common, and it
  # rises with sequencing depth. So compare against a baseline computed WITHIN depth
  # deciles, which holds depth fixed and destroys only the cell-level association.
  output$combo_tbl <- renderTable({
    req(!is.null(gcombo()))
    d <- picked(); req(nrow(d) > 20)
    # pattern mode: plain cross-tab of every +/- combination, no statistics
    if (identical(input$gmode, "pattern") && ".pat" %in% names(d)) {
      lv <- pat_levels(); cnt <- table(factor(d$.pat, levels = lv))
      out <- data.frame(combination = lv, cells = as.integer(cnt[lv]),
                        pct = sprintf("%.2f%%", 100*as.integer(cnt[lv])/nrow(d)))
      return(out[order(-out$cells), ])
    }
    gs <- gcombo(); n <- length(gs)
    ids <- match(d$cell_ID, colnames(emat()))
    hit <- emat()[gs, ids, drop = FALSE] > gthr()
    obs <- as.integer(Matrix::colSums(hit))

    strat <- rep(1L, nrow(d))
    depth_matched <- "total_counts" %in% names(d) && !anyNA(d$total_counts)
    if (depth_matched) {
      br <- unique(stats::quantile(d$total_counts, seq(0, 1, 0.1), na.rm = TRUE))
      if (length(br) > 2) strat <- cut(d$total_counts, br, include.lowest = TRUE, labels = FALSE)
    }
    exp_exact <- numeric(n + 1)
    for (b in unique(strat)) {
      i <- strat == b
      exp_exact <- exp_exact + sum(i) * pb_dist(Matrix::rowMeans(hit[, i, drop = FALSE]))
    }
    o_at <- rev(cumsum(rev(tabulate(obs + 1, n + 1))))
    e_at <- rev(cumsum(rev(exp_exact)))
    data.frame(genes_expressed = paste0(">= ", 0:n),
               cells = as.integer(o_at),
               expected = round(e_at, 1),
               ratio = round(o_at / pmax(e_at, 0.5), 2),
               baseline = c(if (depth_matched) "depth-matched" else "NAIVE (no total_counts in bundle)",
                            rep("", n)))[-1, ]
  })

  # The clustering whose per-cluster tables go with the colouring: the colouring itself, or for a
  # label column the clustering picked in Cluster labels
  label_cols <- c(LIVE_COLS, BASE_LAB, EDIT_COL, "cell_type_display")
  # Coloured by a label column: the tables show the clustering picked in Cluster labels
  on_labels <- reactive((input$cby %||% "") %in% label_cols && length(CL_BASES) > 0)
  tables_basis <- reactive(
    if (on_labels()) input$cl_basis %||% CL_DEFAULT else basis_of(input$cby %||% ""))
  # Its clusters to show: the highlighted levels, or for a label column the clusters labelled
  # with them in Cluster labels
  tables_levels <- reactive({
    h <- hl(); if (is.null(h)) return(NULL)
    if (on_labels()) { d <- cl_tab(); d$cluster[d$label %in% h] } else h
  })
  tri <- reactive(DIAG[[tables_basis()]])

  # --- cluster labels -----------------------------------------------------------------
  # This session's edits per clustering, as named vectors cluster_id -> label. Each Set saves the
  # clustering's labels to the labels folder, where annotate reads them
  cl_edits <- reactiveVal(list())
  # The clustering cell_type comes from, as cluster_basis.txt names it
  chosen <- reactiveVal(CHOSEN0)
  save_basis <- function(b) {
    dir.create(LABEL_DIR, recursive = TRUE, showWarnings = FALSE)
    writeLines(b, BASIS_FILE)
    chosen(b)
  }

  cl_tab <- reactive({
    b <- input$cl_basis; req(b %in% CL_BASES)
    f  <- droplevels(cells[[b]])
    ids <- levels(f)
    n  <- as.integer(table(f)[ids])
    maj <- share <- rep(NA, length(ids))
    if (!is.null(BASE_LAB) && BASE_LAB != b) {
      tb <- table(f, cells[[BASE_LAB]])[ids, , drop = FALSE]
      maj   <- colnames(tb)[max.col(tb, "first")]
      share <- round(100 * apply(tb, 1, max) / pmax(1, rowSums(tb)))
    }
    orig <- if (!is.null(LABS[[b]])) unname(LABS[[b]][ids]) else rep(NA_character_, length(ids))
    from <- ifelse(is.na(orig), "majority",
                   if (!is.null(LABS_SAVED[[b]])) "labels folder" else "pipeline draft")
    orig[is.na(orig)] <- maj[is.na(orig)]
    orig[is.na(orig)] <- ""
    cur <- orig
    ed <- cl_edits()[[b]]
    if (length(ed)) {
      k <- intersect(names(ed), ids)
      cur[match(k, ids)] <- ed[k]
    }
    d <- data.frame(cluster = ids, cells = n, label = cur, original = orig, from = from,
                    majority = maj, share = share, check.names = FALSE)
    t <- DIAG[[b]]
    if (!is.null(t)) {
      i <- match(ids, t$cluster_id)
      d$best_guess  <- t$best_guess[i]
      d$score       <- round(t$best_score[i], 2)
      d$top_markers <- t$top_markers[i]
    }
    # The majority share only means something when the colour labels did not come from this
    # clustering; the marker check is the independent one
    guess <- if (is.null(d$best_guess)) rep("", nrow(d)) else ifelse(is.na(d$best_guess), "", d$best_guess)
    off <- d$label %in% names(LABEL_OK) & nzchar(guess) &
      !mapply(function(l, g) g %in% LABEL_OK[[l]], d$label, guess)
    d$review <- (!is.na(d$share) & d$share < 80) | !nzchar(d$label) | off
    d
  })

  cl_vocab <- reactive(sort(unique(c(
    if (!is.null(BASE_LAB)) levels(cells[[BASE_LAB]]), unlist(LABS, use.names = FALSE),
    unlist(cl_edits(), use.names = FALSE), cl_tab()$label))))

  observeEvent(cl_tab(), {
    d <- cl_tab()
    ch <- stats::setNames(d$cluster, sprintf("%s  ·  %s  (%s)%s", d$cluster,
                                             ifelse(nzchar(d$label), d$label, "unlabelled"),
                                             format(d$cells, big.mark = ","),
                                             ifelse(d$review, "  !", "")))
    sel <- if (isTruthy(input$cl_pick) && input$cl_pick %in% d$cluster) input$cl_pick else d$cluster[1]
    updateSelectInput(session, "cl_pick", choices = ch, selected = sel)
  })

  observeEvent(list(input$cl_pick, input$cl_basis), {
    d <- cl_tab(); req(input$cl_pick %in% d$cluster)
    updateSelectizeInput(session, "cl_label", choices = c("", cl_vocab()),
                         selected = d$label[d$cluster == input$cl_pick])
    # Only while the panel is in view: the picker is also filled in at start-up
    panel_open <- isTRUE(input$tables_open) && identical(input$tables_tab, "Cluster labels")
    if (isTRUE(input$cl_show) && panel_open) {
      if (identical(input$cby, input$cl_basis)) {
        updateSelectizeInput(session, "hl_lv", selected = input$cl_pick)
      } else {
        pend_hl(input$cl_pick)
        updateSelectInput(session, "cby", selected = input$cl_basis)
      }
    }
  })

  observeEvent(input$cl_set, {
    lab <- trimws(input$cl_label %||% ""); req(isTruthy(input$cl_pick))
    if (!nzchar(lab)) return(showNotification("Pick or type a label first.", type = "warning"))
    e <- cl_edits(); b <- input$cl_basis
    v <- e[[b]] %||% character()
    v[input$cl_pick] <- lab
    e[[b]] <- v
    cl_edits(e)
    # Saving a second clustering's labels must not leave annotate without a choice: name the one
    # it was using first, or this one when there was none
    if (is.null(chosen())) save_basis(if (length(LABS_SAVED) == 1L) names(LABS_SAVED) else b)
    d <- cl_tab()
    ok <- tryCatch({
      dir.create(LABEL_DIR, recursive = TRUE, showWarnings = FALSE)
      write_csv_atomic(data.frame(cluster_id = d$cluster, label = d$label), cl_file(b))
      TRUE
    }, error = function(e) FALSE)
    if (!ok) showNotification(paste("Could not save", cl_file(b)), type = "error")
  })

  # Drops this session's edits to this clustering: its saved labels go back to what they were at
  # start, or away when there were none
  observeEvent(input$cl_discard, {
    e <- cl_edits(); b <- input$cl_basis
    e[[b]] <- NULL; cl_edits(e)
    if (!is.null(LABS_SAVED[[b]])) {
      write_csv_atomic(data.frame(cluster_id = names(LABS_SAVED[[b]]),
                                  label = unname(LABS_SAVED[[b]])), cl_file(b))
    } else {
      unlink(cl_file(b))
    }
    showNotification(paste("Discarded the edits for", b), duration = 3)
  })

  # cell_type comes from this clustering from the next annotate on
  observeEvent(input$cl_use, {
    b <- input$cl_basis
    if (!file.exists(cl_file(b))) {
      d <- cl_tab()
      write_csv_atomic(data.frame(cluster_id = d$cluster, label = d$label), cl_file(b))
    }
    save_basis(b)
    showNotification(paste("Annotate takes cell_type from", b), duration = 4)
  })

  observeEvent(input$cl_next, {
    d <- cl_tab(); r <- which(d$review)
    if (!length(r)) return(showNotification("Nothing left to review.", duration = 3))
    i <- match(input$cl_pick, d$cluster)
    nx <- r[r > (i %||% 0)][1]
    updateSelectInput(session, "cl_pick", selected = d$cluster[if (is.na(nx)) r[1] else nx])
  })

  output$cl_hdr <- renderText({
    d <- cl_tab(); b <- input$cl_basis
    sprintf("%d clusters: %d labelled from %s, %d edited, %d to review (marker guess outside the label, majority under 80%%, or no label). Set saves to %s. cell_type comes from %s.",
            nrow(d), sum(d$from != "majority"),
            if (any(d$from != "majority")) unique(d$from[d$from != "majority"]) else "nothing saved",
            sum(d$label != d$original), sum(d$review), cl_file(b),
            if (is.null(chosen())) "no clustering yet (Use for cell_type picks one)" else chosen())
  })

  output$cl_tbl <- renderTable({
    d <- cl_tab()
    d$edited <- ifelse(d$label != d$original, paste("was", d$original), "")
    d$share  <- ifelse(is.na(d$share), "", paste0(d$share, "%"))
    d$review <- ifelse(d$review, "!", "")
    keep <- intersect(c("cluster", "cells", "label", "edited", "majority", "share",
                        "best_guess", "score", "top_markers", "review"), names(d))
    d[, keep]
  })

  output$cl_dl <- downloadHandler(
    filename = function() paste0("cluster_labels.", input$cl_basis, ".csv"),
    content = function(f) {
      d <- cl_tab()
      if (any(!nzchar(d$label)))
        showNotification(sprintf("%d clusters have no label; annotate will stop on them.",
                                 sum(!nzchar(d$label))), type = "warning")
      utils::write.csv(data.frame(cluster_id = d$cluster, label = d$label), f, row.names = FALSE)
    })

  output$tri_hdr <- renderText({
    lv <- tables_levels()
    if (is.null(tri())) sprintf("No diagnostics template for %s.", tables_basis())
    else sprintf("%s — %d clusters%s", tables_basis(), nrow(tri()),
                 if (is.null(lv)) "" else sprintf(", showing %s", paste(lv, collapse = ", ")))
  })

  output$triage <- renderTable({
    d <- tri(); req(!is.null(d))
    lv <- tables_levels()
    if (!is.null(lv) && any(d$cluster_id %in% lv)) d[d$cluster_id %in% lv, ] else d
  })

  output$de_hdr <- renderText({
    d <- DE[[tables_basis()]]
    lv <- tables_levels()
    if (is.null(d)) sprintf("No DE table for %s.", tables_basis())
    # The pipeline's DE is per cluster. A label's genes are those of the clusters labelled with it
    else if (is.null(lv) && on_labels())
      sprintf(paste("Highlight labels to see the top DE genes of the %s clusters labelled with them.",
                    "The pipeline computes DE per clustering, not per label."), tables_basis())
    else if (is.null(lv)) sprintf("%s: pick levels under Highlight to see their top DE genes.", tables_basis())
    else if (!length(lv)) sprintf("No %s cluster carries the highlighted labels.", tables_basis())
    else sprintf("Top DE genes, %s cluster %s", tables_basis(), paste(lv, collapse = ", "))
  })

  output$de_tbl <- renderTable({
    d <- DE[[tables_basis()]]; lv <- tables_levels(); req(!is.null(d), !is.null(lv))
    x <- d[d$cluster %in% lv, ]; req(nrow(x) > 0)
    x <- do.call(rbind, lapply(split(x, x$cluster), function(z) head(z[order(z$ranking), ], 15)))
    out <- data.frame(cluster = x$cluster, gene = x$feats, logFC = round(x$logFC, 3),
                      FDR = signif(x$FDR, 3), rank = as.integer(x$ranking))
    if (on_labels()) {
      cl <- cl_tab()
      out <- cbind(out[1], label = cl$label[match(out$cluster, cl$cluster)], out[-1])
    }
    out
  })
}

shinyApp(ui, server)
