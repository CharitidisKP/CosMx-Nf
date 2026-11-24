## Shared by every page of the run report. The pages only read what the run wrote in results/ ##
suppressPackageStartupMessages({
  library(data.table)
  library(htmltools)
})
source("_functions.R")

if (!exists("%||%", baseenv())) `%||%` <- function(a, b) if (is.null(a)) b else a

## Settings the REPORT step passes, with defaults for a build by hand ##
env <- function(k, d = "") {
  v <- Sys.getenv(k)
  if (nzchar(v)) v else d
}
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}
RES <- normalizePath(env("COSMX_RESULTS", "../results"), mustWork = FALSE)
RUN <- list(
  id = env("COSMX_RUN_ID", basename(dirname(RES))),
  stage = env("COSMX_STAGE", "report"),
  commit = env("COSMX_COMMIT", "unknown"),
  updated = env("COSMX_UPDATED", format(Sys.time(), "%Y-%m-%d %H:%M")),
  chosen = split_arg(env("COSMX_CHOSEN", "")),
  untyped = split_arg(env("COSMX_UNTYPED", "")),
  current = env("COSMX_CURRENT", ""),
  sheet = env("COSMX_SAMPLESHEET", ""),
  config = env("COSMX_CONFIG", ""),
  history = env("COSMX_HISTORY", ""),
  no_overrides = env("COSMX_NO_OVERRIDES", ""),
  fig_style = env("COSMX_FIG_STYLE", ""),
  label_colours = env("COSMX_LABEL_COLOURS", ""),
  labels = env("COSMX_LABELS", file.path(dirname(RES), "labels"))
)
.cache <- new.env()
res_path <- function(...) file.path(RES, ...)
href <- function(path) rel_href(path, RES)
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
rx <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
missing_note <- function(...) tags$p(class = "missing-note", ...)

## ---- Tables ----

## Every published table, listed once per page ##
published_tables <- function() {
  if (is.null(.cache$pub)) {
    .cache$pub <- list.files(res_path("tables"), recursive = TRUE, full.names = TRUE)
  }
  .cache$pub
}

## Tables by file name, this stage's own first, then the published tree ##
find_tab <- function(pattern) {
  cur <- if (nzchar(RUN$current)) Sys.glob(paste0(RUN$current, "*")) else character()
  mine <- if (length(cur)) list.files(cur, pattern, full.names = TRUE) else character()
  pub <- published_tables()
  hits <- c(mine, pub[grepl(pattern, basename(pub))])
  hits[!duplicated(basename(hits))]
}
read_files <- function(f, ..., header = TRUE) {
  if (!length(f)) {
    return(NULL)
  }
  cc <- list(...)[["colClasses"]]
  d <- rbindlist(lapply(f, function(x) as_numbers(fread(x, header = header, ...), cc)), fill = TRUE)
  if (nrow(d)) d else NULL
}
read_tab <- function(pattern, ...) read_files(find_tab(pattern), ...)

## The published copy of a table, for its link ##
published_copy <- function(f) {
  if (startsWith(f, RES)) {
    return(f)
  }
  pub <- published_tables()
  m <- pub[basename(pub) == basename(f)]
  if (length(m)) m[1] else NA_character_
}

## Links to the files a table comes from ##
csv_links <- function(files) {
  files <- unique(stats::na.omit(vapply(files, published_copy, "")))
  files <- files[file.exists(files)]
  if (!length(files)) {
    return(NULL)
  }
  a <- vapply(files, function(f) as.character(tags$a(href = href(f), target = "_blank", basename(f))), "")
  tags$p(class = "csv-links", HTML(paste0(
    if (length(files) > 1L) "Full tables: " else "Full table: ", paste(a, collapse = ", ")
  )))
}

.ids <- new.env()
next_id <- function(prefix) {
  n <- (.ids[[prefix]] %||% 0L) + 1L
  .ids[[prefix]] <- n
  paste0(prefix, n)
}

## Searchable, sortable table, its p columns with threshold filters and the rest rounded for the page ##
hits_table <- function(df, page = 15L, order = NULL, cap = 5000L, csv = NULL, note = NULL) {
  links <- csv_links(csv)
  if (is.null(df) || !nrow(df)) {
    return(tagList(links, missing_note(if (is.null(df)) "Not in this run." else "No rows.")))
  }
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  pc <- names(df)[is_p_col(names(df)) & vapply(df, is.numeric, logical(1))]
  order <- order %||% if (length(pc)) pc[1]
  r <- cap_rows(df, cap, order)
  df <- r$df
  capped <- if (r$total > nrow(df)) {
    tags$p(class = "table-note", sprintf(
      "Showing %s of %s rows%s. The full table has every row.",
      fmt_n(nrow(df)), fmt_n(r$total),
      if (!is.null(r$by)) {
        paste0(", those with the ", if (startsWith(order, "-")) "largest " else "smallest ", r$by)
      } else {
        ""
      }
    ))
  }
  for (n in names(df)) {
    x <- df[[n]]
    if (is.character(x) || is.logical(x)) {
      v <- unique(x[!is.na(x)])
      if (length(v) >= 1L && length(v) <= 300L && max(nchar(as.character(v)), 0L) <= 60L) {
        df[[n]] <- factor(x)
      }
    }
  }
  dec <- vapply(df, function(x) is.numeric(x) && any(x != round(x), na.rm = TRUE), logical(1))
  num <- setdiff(names(df)[dec], pc)
  p_idx <- match(pc, names(df)) - 1L
  opts <- list(
    pageLength = page, lengthMenu = c(15, 50, 200, 1000), autoWidth = FALSE,
    language = list(search = "Search all columns:")
  )
  if (nrow(df) <= min(page, 15L)) opts$dom <- "t"
  if (length(p_idx)) opts$columnDefs <- list(list(searchable = FALSE, targets = p_idx))
  dt <- DT::datatable(
    df,
    rownames = FALSE, filter = "top", elementId = next_id("tbl"), style = "bootstrap5",
    class = "table table-sm table-striped table-hover compact", options = opts,
    callback = DT::JS(if (length(p_idx)) {
      sprintf("cosmxPFilter(table, %s); return table;", jsonlite::toJSON(p_idx))
    } else {
      "return table;"
    })
  )
  if (length(pc)) dt <- DT::formatSignif(dt, pc, digits = 2)
  if (length(num)) dt <- DT::formatRound(dt, num, digits = 3)
  tagList(links, if (!is.null(note)) tags$p(class = "table-note", note), capped, tags$div(class = "dt-wrap", dt))
}

## A table read by file name, shown with its source files linked ##
tab <- function(pattern, ..., cols = NULL, drop = NULL, page = 15L, order = NULL, cap = 5000L, note = NULL) {
  f <- find_tab(pattern)
  d <- read_files(f, ...)
  if (is.null(d) && length(f)) d <- data.frame()
  if (!is.null(d) && nrow(d)) {
    if (!is.null(cols)) d <- d[, intersect(cols, names(d)), with = FALSE]
    if (!is.null(drop)) d <- d[, setdiff(names(d), drop), with = FALSE]
  }
  hits_table(d, page = page, order = order, cap = cap, csv = f, note = note)
}

## A short table as plain HTML, numbers rounded for the page ##
plain_table <- function(df, digits = 3L) {
  if (is.null(df) || !nrow(df)) {
    return(missing_note("No rows."))
  }
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  for (n in names(df)) {
    x <- df[[n]]
    if (is.numeric(x)) {
      df[[n]] <- vapply(x, function(v) {
        if (is.na(v)) "" else if (v == round(v)) fmt_n(v) else fmt_n(signif(v, digits))
      }, "")
    }
  }
  html_table(df)
}

## A data frame as an HTML table built from tags, so Quarto leaves it and what follows alone ##
html_table <- function(df, class = "table table-sm plain") {
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  tags$table(
    class = class,
    tags$thead(tags$tr(lapply(names(df), tags$th))),
    tags$tbody(lapply(seq_len(nrow(df)), function(i) {
      tags$tr(lapply(df[i, , drop = FALSE], function(v) tags$td(if (is.na(v)) "" else as.character(v))))
    }))
  )
}

## ---- Figures ----

## Figures under a folder of figures/, sorted ##
figs <- function(dir, pattern = "\\.(png|jpe?g|gif)$", recursive = TRUE) {
  d <- res_path("figures", dir)
  if (!dir.exists(d)) {
    return(character())
  }
  sort(list.files(d, pattern, full.names = TRUE, recursive = recursive, ignore.case = TRUE))
}

## Figures by name under figures/, those that exist ##
fig_files <- function(...) {
  f <- res_path("figures", ...)
  f[file.exists(f)]
}
nice <- function(f, prefix = "") gsub("_", " ", sub(paste0("^", prefix), "", tools::file_path_sans_ext(basename(f))))

## Figures by file name with their captions, those that exist ##
gallery_named <- function(dir, named, ..., empty = NULL) {
  f <- res_path("figures", dir, names(named))
  ok <- file.exists(f)
  gallery(f[ok], unname(named[ok]), ..., empty = empty)
}

## A PNG copy at most max_px wide, its white margins trimmed on request, cached outside the site folder ##
display_copy <- function(f, max_px = 2400L, trim = FALSE) {
  if (!grepl("\\.png$", f, ignore.case = TRUE) || !file.exists(f)) {
    return(f)
  }
  w <- png_size(f)[1]
  if (!trim && (is.na(w) || w <= max_px)) {
    return(f)
  }
  key <- paste0(gsub("[^A-Za-z0-9]+", "_", tools::file_path_sans_ext(substring(f, nchar(RES) + 2L))), if (trim) "_t")
  out <- res_path("report_display", paste0(key, "_", format(file.mtime(f), "%Y%m%d%H%M%S"), ".png"))
  if (file.exists(out)) {
    return(out)
  }
  ok <- tryCatch(
    {
      dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
      unlink(list.files(dirname(out), paste0("^", key, "_[0-9]{14}\\.png$"), full.names = TRUE))
      img <- magick::image_read(f)
      if (trim) img <- magick::image_border(magick::image_trim(img, fuzz = 2), "white", "24x24")
      if (magick::image_info(img)$width > max_px) img <- magick::image_scale(img, as.character(max_px))
      magick::image_write(img, out, format = "png")
      TRUE
    },
    error = function(e) FALSE
  )
  if (ok) out else f
}

## Figures as cards that open the original. wide: two a row at most. big: one a row, as large as the page ##
gallery <- function(files, labels = nice(files), wide = FALSE, big = FALSE, trim = big, empty = "No figures yet.",
                    scale = NULL) {
  if (!length(files)) {
    return(if (!is.null(empty)) missing_note(empty))
  }
  tags$div(
    class = paste("gal", if (!is.null(scale)) "rep" else if (big) "big" else if (wide) "wide"),
    Map(function(f, l) {
      w <- if (!is.null(scale)) png_size(f)[1] else NA
      tags$figure(
        tags$a(
          href = href(f), target = "_blank",
          tags$img(
            src = href(display_copy(f, trim = trim)), loading = "lazy", alt = l,
            style = if (is.finite(w)) sprintf("width: %dpx", round(w * scale))
          )
        ),
        tags$figcaption(l)
      )
    }, files, labels)
  )
}

## One large figure at a time, with a drop-down of the rest and buttons for the previous and next ##
.pick_js <- "window.pickFigure = window.pickFigure || function(id, i) {
  var d = window.pickData[id], n = d.src.length; i = ((i % n) + n) % n;
  document.getElementById(id + '-img').src = d.src[i];
  document.getElementById(id + '-a').href = d.orig[i];
  document.getElementById(id + '-sel').value = i;
  document.getElementById(id + '-cap').textContent = d.label[i];
  if (d.w) document.getElementById(id + '-img').style.width = d.w[i] ? d.w[i] + 'px' : '';
};"
figure_picker <- function(files, labels = nice(files), trim = TRUE, empty = "No figures yet.", scale = NULL) {
  if (!length(files)) {
    return(missing_note(empty))
  }
  id <- next_id("pick")
  src <- unname(vapply(files, function(f) href(display_copy(f, trim = trim)), ""))
  orig <- unname(vapply(files, href, ""))
  wpx <- if (!is.null(scale)) unname(vapply(files, function(f) round(png_size(f)[1] * scale), 0)) else NULL
  step <- function(by) sprintf("pickFigure('%s', Number(document.getElementById('%s-sel').value) + %d)", id, id, by)
  tags$div(
    class = "picker",
    tags$div(
      class = "picker-bar",
      tags$button(type = "button", class = "btn btn-sm btn-outline-secondary", onclick = step(-1), HTML("&#9664;")),
      tags$select(
        id = paste0(id, "-sel"), class = "form-select form-select-sm picker-select",
        onchange = sprintf("pickFigure('%s', Number(this.value))", id),
        Map(function(i, l) tags$option(value = i - 1L, l), seq_along(files), labels)
      ),
      tags$button(type = "button", class = "btn btn-sm btn-outline-secondary", onclick = step(1), HTML("&#9654;")),
      tags$span(class = "picker-count", paste(length(files), if (length(files) == 1L) "figure" else "figures"))
    ),
    tags$a(
      id = paste0(id, "-a"), href = orig[1], target = "_blank",
      tags$img(
        id = paste0(id, "-img"), src = src[1], class = "picker-img", alt = labels[1],
        style = if (length(wpx) && is.finite(wpx[1])) sprintf("width: %dpx", wpx[1])
      )
    ),
    tags$div(id = paste0(id, "-cap"), class = "picker-cap", labels[1]),
    tags$script(HTML(paste0(
      .pick_js, "\nwindow.pickData = window.pickData || {}; window.pickData['", id, "'] = ",
      jsonlite::toJSON(list(src = src, orig = orig, label = unname(labels), w = wpx)), ";"
    )))
  )
}

## ---- Layout ----

## Bootstrap tabs from a named list of contents, the empty ones dropped ##
tabs <- function(x) {
  x <- x[lengths(x) > 0]
  if (!length(x)) {
    return(missing_note("Nothing to show yet."))
  }
  ids <- vapply(seq_along(x), function(i) next_id("rtab"), "")
  tagList(
    tags$ul(
      class = "nav nav-tabs", role = "tablist",
      Map(function(n, id, i) {
        tags$li(
          class = "nav-item", role = "presentation",
          tags$a(
            class = paste("nav-link", if (i == 1L) "active"), `data-bs-toggle` = "tab",
            href = paste0("#", id), role = "tab", n
          )
        )
      }, names(x), ids, seq_along(x))
    ),
    tags$div(
      class = "tab-content",
      Map(function(cont, id, i) {
        tags$div(class = paste("tab-pane fade", if (i == 1L) "show active"), id = id, role = "tabpanel", cont)
      }, x, ids, seq_along(x))
    )
  )
}

## Sections with real headings, which reach the contents list. A function is a section's content, a list holds subsections ##
render_tree <- function(tree, level = 2L) {
  leaves <- list()
  lines <- character()
  walk <- function(node, level) {
    for (i in seq_along(node)) {
      nm <- names(node)[i]
      item <- node[[i]]
      if (is.null(item)) next
      if (!startsWith(nm, ".")) lines <<- c(lines, paste(strrep("#", min(level, 6L)), nm), "")
      if (is.function(item)) {
        leaves[[length(leaves) + 1L]] <<- local({
          f <- item
          function() {
            out <- tryCatch(f(), error = function(e) missing_note("Not shown: ", conditionMessage(e)))
            if (is.null(out)) invisible(NULL) else out
          }
        })
        col <- attr(item, "column")
        lines <<- c(
          lines, "```{r}", if (!is.null(col)) paste0("#| column: ", col),
          sprintf(".leaves[[%d]]()", length(leaves)), "```", ""
        )
      } else {
        walk(item, if (startsWith(nm, ".")) level else level + 1L)
      }
    }
  }
  walk(tree, level)
  e <- new.env(parent = globalenv())
  e$.leaves <- leaves
  knitr::knit_child(text = lines, envir = e, quiet = TRUE)
}

## A section whose figures may also use the page's right margin ##
wide <- function(f, column = "page-right") structure(f, column = column)

## Figures as cards when the figure helpers are loaded ##
rgallery_if <- function(items) if (exists("rgallery")) rgallery(items, wide = TRUE)

## Per clustering column, its score heatmap, label template and top genes per cluster, one tab each ##
marker_tabs <- function(cols, labels = cols) {
  tmpl <- function(cc) find_tab(paste0("^cluster_labels\\.", rx(cc), "\\.template\\.csv$"))
  top <- function(cc) find_tab(paste0("^", rx(cc), "_de_top\\.csv$"))
  have <- vapply(cols, function(cc) length(tmpl(cc)) + length(top(cc)) > 0L, logical(1))
  if (!any(have)) {
    return(missing_note("No marker diagnostics yet."))
  }
  out <- lapply(cols[have], function(cc) {
    more <- find_tab(paste0(
      "^", rx(cc), "_(de_all|marker_scores|marker_scores_mean|marker_panel_coverage|compartment_scores)\\.csv$"
    ))
    scores <- find_tab(paste0("^", rx(cc), "_marker_scores_mean\\.csv$"))
    heat <- if (exists("rfig") && length(scores)) {
      sz <- cluster_heat_size(ncol(fread(scores[1], nrows = 0L, header = TRUE)) - 1L, nrow(fread(scores[1], select = 1L, header = TRUE)))
      rfig(paste0("marker_scores_", cc), scores, function() {
        d <- fread(scores[1], header = TRUE)
        setnames(d, 1, "panel")
        long <- melt(d, id.vars = "panel", variable.name = "cluster", value.name = "score")
        lv <- cluster_levels(as.character(long$cluster))
        panels <- sort(unique(long$panel))
        if (sz[["down"]] == 1) {
          long[, `:=`(cluster = factor(as.character(cluster), rev(lv)), panel = factor(panel, panels))]
          heat_plot(long, "panel", "cluster", "score", name = "Mean score")
        } else {
          long[, `:=`(cluster = factor(as.character(cluster), lv), panel = factor(panel, rev(panels)))]
          heat_plot(long, "cluster", "panel", "score", name = "Mean score")
        }
      }, sz[["w"]], sz[["h"]])
    }
    tagList(
      rgallery_if(list(`Mean marker-panel score per cluster` = heat)),
      tags$p(tags$b("Label template: "), "the best marker-panel guess per cluster, the draft the labels start from."),
      hits_table(read_files(tmpl(cc), colClasses = list(character = 1)), page = 50, csv = tmpl(cc)),
      tags$p(tags$b("Top genes per cluster")),
      hits_table(read_files(top(cc), colClasses = list(character = 1)), csv = top(cc)),
      csv_links(more)
    )
  })
  tabs(setNames(out, labels[have]))
}

## ---- Settings ----

## The settings in nextflow.config, with their defaults and comments ##
config_table <- function() {
  if (is.null(.cache$config)) {
    .cache$config <- if (nzchar(RUN$config) && file.exists(RUN$config)) {
      config_params(readLines(RUN$config, warn = FALSE))
    } else {
      config_params(character())
    }
  }
  .cache$config
}

## The settings a stage's last launch ran with ##
stage_params <- function(stage) {
  f <- res_path("pipeline_info", paste0(stage, "_params.json"))
  if (!file.exists(f)) {
    return(NULL)
  }
  tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
}

## The stage of a group that launched last, from the stage log ##
latest_stage <- function(stages) {
  f <- res_path("pipeline_info", "stage_log.tsv")
  if (!file.exists(f)) {
    return(stages[1])
  }
  lg <- tryCatch(fread(f, sep = "\t", quote = ""), error = function(e) NULL)
  if (is.null(lg) || !nrow(lg)) {
    return(stages[1])
  }
  lg <- lg[stage %in% stages]
  if (nrow(lg)) lg$stage[nrow(lg)] else stages[1]
}

## Each setting with the value the stage ran and what it does, from the comment in nextflow.config ##
params_table <- function(names, stage) {
  cfg <- config_table()
  p <- stage_params(stage)
  val <- vapply(names, function(n) {
    v <- p[[n]]
    if (!is.null(v)) {
      return(paste(as.character(unlist(v)), collapse = ","))
    }
    d <- cfg$default[match(n, cfg$name)]
    if (is.na(d)) "" else d
  }, "")
  doc <- cfg$doc[match(names, cfg$name)]
  df <- data.frame(Setting = names, Value = val, `What it does` = ifelse(is.na(doc), "", doc), check.names = FALSE)
  tagList(
    html_table(df, "table table-sm plain params"),
    tags$p(class = "table-note", if (is.null(p)) {
      sprintf("The %s stage has not run, so these are the defaults in nextflow.config.", stage)
    } else {
      sprintf("As the last %s launch ran them (pipeline_info/%s_params.json). What each does is its comment in nextflow.config.", stage, stage)
    })
  )
}

## ---- Run status ----

STAGES <- c("discover", "checks", "subcluster", "diagnose", "supervised", "annotate", "spatial", "cci", "de", "banksy_sweep", "report")

## The run's launches in Nextflow's history ##
launches <- function() {
  if (is.null(.cache$launches)) {
    h <- if (nzchar(RUN$history) && file.exists(RUN$history)) readLines(RUN$history, warn = FALSE) else character()
    .cache$launches <- history_launches(h, RUN$id)
  }
  .cache$launches
}

## True when the stage's object or tables are in results/ ##
stage_made <- function(s) {
  switch(s,
    discover = ,
    subcluster = ,
    annotate = ,
    spatial = file.exists(res_path("objects", s, "giotto", "gobject.RDS")),
    cci = dir.exists(res_path("tables", "08_cci")),
    de = dir.exists(res_path("tables", "09_de")),
    FALSE
  )
}

## The clustering the labels name, from cluster_basis.txt or the only labels file ##
labels_basis <- function() {
  f <- file.path(RUN$labels, "cluster_basis.txt")
  if (file.exists(f)) {
    return(trimws(readLines(f, n = 1L, warn = FALSE)))
  }
  l <- list.files(RUN$labels, "^cluster_labels\\..+\\.csv$")
  if (length(l) == 1L) sub("^cluster_labels\\.(.+)\\.csv$", "\\1", l) else NULL
}

## True when everything the stage wrote is older than the object it reads ##
stage_outdated <- function(s) {
  io <- switch(s,
    annotate = c("objects/annotate/giotto", if (startsWith(labels_basis() %||% "", "sub_")) "subcluster" else "discover"),
    spatial = c("objects/spatial/giotto", "annotate"),
    cci = c("tables/08_cci", "spatial"),
    de = c("tables/09_de", "spatial"),
    NULL
  )
  if (is.null(io)) {
    return(FALSE)
  }
  obj <- res_path("objects", io[2], "giotto", "gobject.RDS")
  out <- list.files(res_path(io[1]), recursive = TRUE, full.names = TRUE)
  outdated_by(if (file.exists(obj)) as.numeric(file.mtime(obj)) else NA, as.numeric(file.mtime(out)))
}

## A stage's state, its last launch and how the badge shows it ##
stage_status <- function(s) {
  L <- launches()
  state <- launch_state(s, L, made = stage_made(s), current = RUN$stage)
  if (state == "done" && stage_outdated(s)) state <- "outdated"
  mine <- L[L$stage == s, , drop = FALSE]
  last <- if (nrow(mine)) mine[nrow(mine), , drop = FALSE] else NULL
  list(
    stage = s, state = state,
    class = c(done = "ok", outdated = "stale", running = "running", failed = "failed", stopped = "failed", none = "missing")[[state]],
    label = c(done = "done", outdated = "outdated", running = "running now", failed = "failed", stopped = "stopped", none = "not run")[[state]],
    time = if (!is.null(last)) substr(last$time, 1L, 16L) else "",
    duration = if (!is.null(last)) last$duration else "",
    launches = nrow(mine),
    options = if (!is.null(last)) stage_options(last$command) else ""
  )
}
badge <- function(st) tags$span(class = paste("status", st$class), st$label)
status_detail <- function(st) {
  if (nzchar(st$time)) {
    paste0("last launch ", st$time, if (nzchar(st$duration) && st$duration != "-") paste0(", took ", st$duration))
  } else if (st$state == "done") {
    "outputs in results/, launched from another folder"
  } else {
    ""
  }
}

## Status lines for the folded box at the top of a page ##
status_lines <- function(stages) {
  tags$div(lapply(stages, function(s) {
    st <- stage_status(s)
    tags$div(class = "status-line", tags$code(s), " ", badge(st), " ", tags$small(status_detail(st)))
  }))
}

## Under the title, only when a stage behind the page failed or is out of date ##
stale_warning <- function(stages) {
  bad <- unlist(lapply(stages, function(s) {
    st <- stage_status(s)
    if (st$class %in% c("stale", "failed")) as.character(tagList(tags$code(s), " ", badge(st)))
  }))
  if (!length(bad)) {
    return(NULL)
  }
  tags$div(class = "stale-note", HTML(paste0(
    paste(bad, collapse = " "), " See <a href=\"status.html\">Run status</a>."
  )))
}

## The folded box with how the page's outputs were made, the stages behind them first ##
made <- function(..., stages = NULL) {
  items <- lapply(list(...), HTML)
  if (length(stages)) items <- c(list(Run = status_lines(stages)), items)
  tags$details(
    class = "made", tags$summary("How it was made"),
    tags$dl(Map(function(n, v) tagList(tags$dt(n), tags$dd(v)), names(items), items))
  )
}
page_top <- function(..., stages = NULL) tagList(stale_warning(stages), made(..., stages = stages))

## Whether annotate was built from the labels in the run's labels folder ##
labels_state <- function() {
  n <- length(list.files(RUN$labels, "^cluster_labels\\..+\\.csv$"))
  basis <- labels_basis()
  if (!n) {
    return(list(class = "missing", label = "no labels", text = "No labels in the run's labels folder yet."))
  }
  if (is.null(basis)) {
    return(list(class = "stale", label = "no basis", text = sprintf("%d clusterings have labels, and cluster_basis.txt names none of them.", n)))
  }
  lab <- file.path(RUN$labels, paste0("cluster_labels.", basis, ".csv"))
  if (!file.exists(lab)) {
    return(list(class = "failed", label = "no labels file", text = sprintf("cluster_basis.txt names %s, which has no cluster_labels file.", basis)))
  }
  ovr <- file.path(RUN$labels, "label_overrides.csv")
  edits <- file.exists(ovr)
  if (!edits) ovr <- RUN$no_overrides
  what <- sprintf("cell_type comes from %s, %s.", basis, if (edits) "with per-cell edits" else "with no per-cell edits")
  built_f <- res_path("objects", "annotate", "label_state.txt")
  if (!file.exists(built_f)) {
    return(list(class = "missing", label = "not used yet", text = paste(what, "annotate has not used them yet.")))
  }
  if (!file.exists(ovr)) {
    return(list(class = "ok", label = "labels", text = what))
  }
  same <- identical(trimws(readLines(built_f, n = 1L, warn = FALSE)), label_digest(basis, lab, ovr))
  if (same) {
    list(class = "ok", label = "up to date", text = paste(what, "annotate is up to date with them."))
  } else {
    list(class = "stale", label = "changed", text = paste(what, "They changed after annotate ran, so cci and de rerun annotate and spatial first."))
  }
}

## ---- Pages and samples ----

PAGES <- data.table(
  page = c(
    "01_qc.html", "02_integration.html", "03_clustering.html",
    "04_hieratype.html", "04_insitutype.html", "04_supervised.html",
    "05_subclustering.html", "06_annotation.html",
    "07_neighbours.html", "07_density.html", "07_niches.html", "07_bcell.html", "07_sweep.html",
    "08_liana.html", "08_nichenet.html", "08_misty.html", "08_svg.html", "08_insitucor.html",
    "09_pseudobulk.html", "09_smide.html"
  ),
  title = c(
    "01 - Quality control", "02 - Integration", "03 - Clustering",
    "04 - HieraType", "04 - InSituType, unsupervised", "04 - InSituType, supervised",
    "05 - Subclustering", "06 - Annotation",
    "07 - Neighbours", "07 - Immune density", "07 - Niches", "07 - B-lineage cells", "07 - BANKSY sweep",
    "08 - LIANA", "08 - NicheNet", "08 - MISTy", "08 - SPARK-X and nnSVG", "08 - InSituCor",
    "09 - Pseudobulk", "09 - smiDE"
  ),
  holds = c(
    "Cells per sample, what each filter removed, the counts, genes and area per cell, and the FOVs.",
    "Batch correction on the UMAP, the PCs the data supports and the gene-selection methods.",
    "Leiden over the resolution ladder, how well each resolution holds, the UMAPs and the marker diagnostics.",
    "HieraType's marker-based calls.",
    "InSituType's clusters and cohorts without a reference.",
    "InSituType's calls against each reference profile.",
    "Each cluster re-clustered on its own embedding.",
    "The curated labels on the cells, the cell types and their composition.",
    "The neighbour network, which cell types sit together, and the cell-type maps.",
    "Immune cells per mm2 of tissue and as a share of all cells.",
    "BANKSY niches per sample, matched into niche types across samples.",
    "The B-lineage cells re-clustered on their own.",
    "Every BANKSY settings row on every biopsy.",
    "Ligand-receptor pairs between neighbouring cell types.",
    "The ligands behind each comparison's changes in each receiver.",
    "How much of each gene the surroundings explain.",
    "Genes whose expression follows space.",
    "Gene modules that vary together in space, and the cell types behind them.",
    "Each comparison per cell type, paired on the subject.",
    "Each cell type's genes between niches, per sample and pooled."
  ),
  stage = c(
    "discover", "discover", "discover", "discover", "discover", "discover",
    "subcluster", "annotate", "spatial", "spatial", "spatial", "spatial", "banksy_sweep",
    "cci", "cci", "cci", "cci", "cci", "de", "de"
  )
)

## The samplesheet without paths or patient codes ##
samples <- function() {
  if (!nzchar(RUN$sheet) || !file.exists(RUN$sheet)) {
    return(NULL)
  }
  ss <- fread(RUN$sheet, colClasses = "character")
  ss[, intersect(c("sample_id", "include", "treatment", "timepoint", "subject_id", "slide_id", "batch"), names(ss)), with = FALSE]
}
