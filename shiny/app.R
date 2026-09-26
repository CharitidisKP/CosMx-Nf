library(shiny)
library(ggplot2)
library(scattermore)
source("theme.R")

# Which results tree to read. Override to compare runs without overwriting a bundle:
#   COSMX_NF=../../Runs/my_run/results R -e 'shiny::runApp("shiny")'
NF <- Sys.getenv("COSMX_NF", "../results")

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
SUB_EMB <- {
  u1 <- grep("_UMAP1$", names(cells), value = TRUE)
  stats::setNames(lapply(u1, function(a) c(a, sub("1$", "2", a))),
                  paste("Subcluster", sub("^.*_c(.+)_UMAP1$", "\\1", u1)))
}
SUB_EMB <- SUB_EMB[vapply(SUB_EMB, function(z) all(z %in% names(cells)), logical(1))]

EMB <- Filter(Negate(is.null), c(list(
  Spatial = emb("x", "y"),
  UMAP    = emb("UMAP1", "UMAP2") %||% emb("UMAP_1", "UMAP_2"),
  # The same cells embedded on PCA instead of harmony. Present only when the run corrected
  # for batch, and carried in the same cells table so colour, gene and sample filters hold
  # across the pair -- switching tabs is then a clean before/after of the correction alone.
  "UMAP pre-batch" = emb("UMAPpre1", "UMAPpre2"),
  PCA     = emb("PCA_1", "PCA_2")
), SUB_EMB))
# Output ids, brush ids and the per-page facet memory are keyed by a slug: a display name
# with a space or a hyphen in it cannot be part of an input id.
VID <- c(Spatial = "spatial", UMAP = "umap",
         "UMAP pre-batch" = "umap_pre", PCA = "pca")[names(EMB)]
# Subcluster tabs are discovered, not listed, so they need slugs generated
VID[is.na(VID)] <- paste0("sub", gsub("[^A-Za-z0-9]", "", names(EMB)[is.na(VID)]))
names(VID) <- names(EMB)
CAT   <- c(names(Filter(is.factor, cells)), LIN_COL)
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
  # intersect with names(cells): CAT also carries the lineage column, which is reactive
  # and has no column to test
  c(SAMP, names(Filter(const_in_sample, cells[intersect(setdiff(CAT, SAMP), names(cells))])))
}
# GRP still drives Split by; on Spatial the samples are stacked, so colouring by treatment or
# timepoint blocks them visibly and is worth having.
cat_for <- function(view) CAT

# selectInput renders a named list as <optgroup>s. Group by what the column IS, so the
# annotation columns stop being an undifferentiated wall of names.
GROUPS <- list(
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
  function sendVh() { Shiny.setInputValue('cosmx_vh', window.innerHeight); }
  $(document).on('shiny:connected', function () {
    Shiny.addCustomMessageHandler('cosmx-zoomed', function (z) { zoomed = z; });
    sendVh();
  });
  window.addEventListener('resize', function () {
    clearTimeout(window._vht); window._vht = setTimeout(sendVh, 300);
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

APP_CSS <- "
html, body { overscroll-behavior-x: none; }
.selbar { position: sticky; top: 0; z-index: 20; display: flex; gap: 8px; align-items: center;
          background: #fff; padding: 4px 6px; border-bottom: 1px solid #E4E9E7; }
.plot-toolbar .dropdown-menu { font-size: .9rem; }
.legend-row { white-space: nowrap; font-size: .88rem; margin-bottom: 2px; }
.legend-sw { display: inline-block; width: 11px; height: 11px; border-radius: 2px;
             margin-right: 7px; vertical-align: -1px; }
"

scroll_box <- function(...) div(style = "max-height: 260px; overflow-y: auto; font-size: .9em;",
                                ...)

# Bootstrap dropdown that stays open while its inputs are used. Fixed positioning lets the
# menu escape the card's overflow clipping.
dd <- function(label, ..., width = "320px") div(
  class = "dropdown",
  tags$button(class = "btn btn-sm btn-outline-primary dropdown-toggle", type = "button",
              `data-bs-toggle` = "dropdown", `data-bs-auto-close` = "outside",
              `data-bs-popper-config` = '{"strategy":"fixed"}', label),
  div(class = "dropdown-menu p-3", style = paste0("width:", width, ";"), ...))

toolbar <- div(
  class = "plot-toolbar d-flex flex-wrap gap-2 align-items-center",
  if (!is.null(SAMP)) dd(
    "Samples",
    checkboxGroupInput("samples", NULL, levels(cells[[SAMP]]), selected = levels(cells[[SAMP]])),
    actionLink("samples_all", "All"), " · ", actionLink("samples_none", "None"),
    width = "260px"),
  if (!is.null(FOVS)) conditionalPanel(
    "input.view == 'Spatial'",
    div(style = "width: 210px; margin-bottom: -1rem;",
        selectizeInput("fov", NULL, choices = c("Jump to FOV" = "", FOVS),
                       options = list(placeholder = "Jump to FOV")))),
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
    width = "340px"),
  tagAppendAttributes(
    actionButton("zoom_reset", "Reset zoom", class = "btn-sm btn-outline-secondary"),
    title = "Zoom: drag a box and double-click it, pinch, or Cmd/Ctrl + scroll. Double-click empty space to reset."),
  span(class = "text-muted small", textOutput("filter_summary", inline = TRUE)))

plot_card <- function(id) bslib::card(
  full_screen = TRUE, fill = FALSE,
  bslib::layout_sidebar(
    fillable = FALSE, fill = FALSE,
    sidebar = bslib::sidebar(
      position = "right", width = 300, open = "always",
      div(style = "max-height: 80vh; overflow: auto;",
          uiOutput(paste0("cellinfo_", id)),
          uiOutput(paste0("legend_", id)))),
    bslib::card_body(fillable = FALSE, fill = FALSE,
                     style = "overflow: auto; max-height: 86vh; padding: 6px;",
                     uiOutput(paste0("selbar_", id)),
                     uiOutput(paste0("ui_", id)))))

ui <- bslib::page_sidebar(
  title = "CosMx explorer",
  theme = app_theme,
  tags$head(tags$style(HTML(APP_CSS)), tags$script(HTML(ZOOM_JS))),
  sidebar = bslib::sidebar(
    width = 310,
    selectInput("cby", "Colour by", grouped_choices(cat_for("Spatial")),
                selected = pick("label_res0.8", "insitutype_unsup", "anno_auto")),
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

  # Collapsed by default and height-capped: the plot is the point, these are a reference
  bslib::accordion(
    open = FALSE, class = "mt-2",
    bslib::accordion_panel(
      "Tables", icon = NULL,
      bslib::navset_card_tab(
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
    if (!is.null(LIN_COL) && identical(input$cby, LIN_COL)) lin_vec() else cells[[input$cby]]
  })

  spatial <- reactive(identical(input$view, "Spatial"))
  vid     <- reactive({ v <- unname(VID[input$view %||% ""]); if (is.na(v)) "spatial" else v })
  zoom    <- reactiveValues(x = NULL, y = NULL)
  sel_ids <- reactiveVal(NULL)   # cell_IDs kept by "Filter cells to selected"

  # Plots are always width 100%, so the width Shiny reports back is the container width
  base_w <- reactiveVal(900)
  observe({
    w <- session$clientData[[paste0("output_plot_", vid(), "_width")]]
    if (!is.null(w) && w > 200) base_w(round(w / 50) * 50)
  })

  keep_s <- reactive({
    if (is.null(SAMP)) return(rep(TRUE, nrow(cells)))
    req(input$samples)
    cells[[SAMP]] %in% input$samples
  })
  # Greying 170,000 cells still draws them, and the grey is what you end up looking at.
  # Dropping them is the only way to see the 200 that are left.
  keep <- reactive({
    k <- keep_s()
    cc <- score_of(input$cby)
    if (!is.null(cc) && (input$conf_min %||% 0) > 0)
      k <- k & !is.na(cells[[cc]]) & cells[[cc]] >= input$conf_min
    # A subcluster embedding only holds its own parent's cells; the rest have no coordinates
    ax <- EMB[[input$view %||% names(EMB)[1]]]
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
    if (stacked()) c(".sx", ".sy") else EMB[[input$view %||% names(EMB)[1]]]
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
    lab  <- 30 * diff(xlim) / max(200, base_w() - 24)     # ~30 px label row, in data units
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
  observeEvent(cby_all(), {
    lv <- levels(droplevels(cby_all()))
    for (id in c("keep_lv", "hl_lv"))
      updateSelectizeInput(session, id, choices = lv, selected = intersect(input[[id]], lv))
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
  view_h   <- reactive(round(0.78 * (input$cosmx_vh %||% 900)))
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
  observeEvent(list(input$samples, input$layout), reset_zoom(), ignoreInit = TRUE)

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
    val <- function(col) if (col %in% names(r)) as.character(r[[col]]) else NA_character_
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
  }

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
    p <- pal_for(cby_all())
    if (!is.null(hl())) p[!names(p) %in% hl()] <- GREY
    p
  })

  legend_ui <- function() {
    req(input$cby)
    f <- droplevels(factor(shown()[[input$cby]]))
    n <- table(f); p <- pal_hl()
    tagList(lapply(levels(f), function(l) div(
      class = "legend-row",
      span(class = "legend-sw", style = paste0("background:", p[[l]], ";")), l,
      span(style = "color:#8a8a8a;font-size:.85em", HTML("&nbsp;"),
           format(as.integer(n[[l]]), big.mark = ",")))))
  }

  size_ui <- function(id) {
    out <- function(h) tagAppendAttributes(
      plotOutput(paste0("plot_", id), height = h,
                 brush = brushOpts(paste0("br_", id), resetOnNew = TRUE),
                 click = paste0("ck_", id), dblclick = paste0("dbl_", id),
                 hover = hoverOpts(paste0("hv_", id), delay = 60, delayType = "throttle",
                                   nullOutside = TRUE)),
      class = "cosmx-zoom")
    if (!spatial()) {
      if (n_panels() == 1) return(out("76vh"))
      return(out(paste0(min(20000, ceiling(n_panels() / 2) * round(base_w() / 2 * 0.8) + 40), "px")))
    }
    if (zoomed()) return(out(paste0(view_h(), "px")))
    # coord_fixed() pins the aspect, so the height is exactly what the panel needs at this
    # width: 24 px of plot margin, plus the colour bar when a gene view puts one on top
    fr  <- full_range()
    asp <- diff(fr$y) / max(1e-9, diff(fr$x))
    leg <- legend_px()
    out(paste0(max(220, min(20000, round((base_w() - 24) * asp + 24 + leg))), "px"))
  }

  # geom_point below RASTER_ABOVE cells, scattermore above it. Measured on this bundle,
  # geom_point beat scattermore at every canvas size tried, by 5-16x:
  #   126k cells, 1188x792 px   scatter 2.35 s   point 0.34 s
  #   126k cells, 1320x3300 px  scatter 5.82 s   point 0.48 s
  # scattermore is kept only as a guard for datasets far larger than this one.
  RASTER_ABOVE <- 400000L
  PT <- 1   # point size; points are only used for the embeddings

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

    # Each mode below picks the cells in draw order (d), the column that colours them (col)
    # and a scale; the geom is chosen once, after. Scales serve colour AND fill, so polygons
    # follow the gene views too -- they used to fill by Colour by whatever else was set.
    both <- c("colour", "fill")
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
      if (!is.null(hl())) d <- d[order(d[[input$cby]] %in% hl()), ]   # highlighted on top
      col <- input$cby
      sc  <- scale_colour_manual(values = pal_hl(), limits = names(pal_hl()), drop = FALSE, na.value = GREY,
                                 name = input$cby, aesthetics = both)
    }

    p <- if (use_poly()) {
      pd <- poly_sub()
      # group = the cell's row in d: GeomPolygon draws in group order, so expressing and
      # highlighted cells land on top exactly as they do as points
      pd$.g <- match(pd$i, match(d$cell_ID, cells$cell_ID))
      pd <- pd[!is.na(pd$.g), , drop = FALSE]
      pd$.v <- d[[col]][pd$.g]
      p <- ggplot(pd, aes(px, py, group = .g, fill = .v)) + geom_polygon(linewidth = 0) + sc
      hit <- pd[cells$cell_ID[pd$i] %in% sel_cell(), , drop = FALSE]
      if (nrow(hit)) p <- p + geom_polygon(data = hit, aes(px, py, group = .g), inherit.aes = FALSE,
                                           fill = NA, colour = "#000000", linewidth = 0.8)
      p
    } else {
      p <- ggplot(d, aes(.data[[xy[1]]], .data[[xy[2]]], colour = .data[[col]])) +
        pts(nrow(d), PT) + sc
      if (col %in% c(".pat", ".combo"))
        p <- p + guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
      hit <- d[d$cell_ID %in% sel_cell(), , drop = FALSE]
      if (nrow(hit)) p <- p + geom_point(data = hit, aes(.data[[xy[1]]], .data[[xy[2]]]),
                                         inherit.aes = FALSE, shape = 21, size = 4, stroke = 1.2,
                                         colour = "#000000", fill = NA)
      p
    }

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
                          hjust = 0, vjust = 0.5, size = 3.6, colour = "#4A4A4A")
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
    output[[paste0("plot_", id)]]   <- renderPlot(build_plot(id), res = 132)

    # Drawing a box pops up the selection bar; applying it narrows every view to those cells
    output[[paste0("selbar_", id)]] <- renderUI({
      b <- input[[br]]; req(b)
      n <- nrow(brushedPoints(vis(), b, axes()[1], axes()[2]))
      div(class = "selbar",
          strong(format(n, big.mark = ","), " cells selected"),
          actionButton(paste0("sel_apply_", id), "Filter cells to selected",
                       class = "btn-sm btn-primary"),
          actionButton(paste0("sel_cancel_", id), "Cancel", class = "btn-sm btn-outline-secondary"),
          span(class = "text-muted small", "or double-click the box to zoom"))
    })
    observeEvent(input[[paste0("sel_apply_", id)]], {
      b <- input[[br]]; req(b)
      ids <- brushedPoints(vis(), b, axes()[1], axes()[2])$cell_ID
      sel_ids(if (is.null(sel_ids())) ids else intersect(sel_ids(), ids))
      session$resetBrush(br)
    })
    observeEvent(input[[paste0("sel_cancel_", id)]], session$resetBrush(br))

    # Double-click inside a box zooms to it; anywhere else resets
    observeEvent(input[[paste0("dbl_", id)]], {
      b <- input[[br]]
      if (is.null(b)) return(reset_zoom())
      zoom_to(b)
      session$resetBrush(br)
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
    names(out)[-(1:2)] <- names(pan)
    out[order(-out$n_cells), ]
  })
  output$mk_tbl <- renderTable(mk_eval(), digits = 1)
  output$mk_hdr <- renderText({
    if (!length(MK_SET)) return("assets/canonical_markers.csv not found beside the app.")
    pan <- MK_SET[[input$mkby %||% "cell_type"]]
    sprintf("Fold over %s background, depth-matched on %d count deciles. Panels: %s",
            if (is.null(DEC)) "NAIVE (no total_counts)" else "sample-filtered",
            NDEC, paste(sprintf("%s %d", names(pan), lengths(pan)), collapse = ", "))
  })

  picked <- reactive({
    xy <- axes()
    if (is.null(brush())) vis() else brushedPoints(vis(), brush(), xy[1], xy[2])
  })

  output$comp_hdr <- renderText({
    sprintf("%s of %s cells%s%s", format(nrow(picked()), big.mark = ","),
            format(nrow(vis()), big.mark = ","),
            if (is.null(brush())) "" else " (brushed)",
            if (is.null(zoom$x)) "  --  drag a box, double-click it to zoom (or pinch, Cmd/Ctrl + scroll)"
            else "  --  zoomed; double-click outside a box to reset")
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

  tri <- reactive(DIAG[[basis_of(input$cby)]])

  output$tri_hdr <- renderText({
    if (is.null(tri())) sprintf("No diagnostics template for %s.", basis_of(input$cby))
    else sprintf("%s — %d clusters%s", basis_of(input$cby), nrow(tri()),
                 if (is.null(hl())) "" else sprintf(", showing %s", paste(hl(), collapse = ", ")))
  })

  output$triage <- renderTable({
    d <- tri(); req(!is.null(d))
    if (!is.null(hl()) && any(d$cluster_id %in% hl())) d[d$cluster_id %in% hl(), ] else d
  })

  output$de_hdr <- renderText({
    d <- DE[[basis_of(input$cby)]]
    if (is.null(d)) sprintf("No DE table for %s.", basis_of(input$cby))
    else if (is.null(hl())) "Pick levels under Highlight to see their top DE genes."
    else sprintf("Top DE genes — cluster %s", paste(hl(), collapse = ", "))
  })

  output$de_tbl <- renderTable({
    d <- DE[[basis_of(input$cby)]]; req(!is.null(d), !is.null(hl()))
    x <- d[d$cluster %in% hl(), ]; req(nrow(x) > 0)
    x <- do.call(rbind, lapply(split(x, x$cluster), function(z) head(z[order(z$ranking), ], 15)))
    data.frame(cluster = x$cluster, gene = x$feats, logFC = round(x$logFC, 3),
               FDR = signif(x$FDR, 3), rank = as.integer(x$ranking))
  })
}

shinyApp(ui, server)
