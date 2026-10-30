## The app saves its label edits into the run's labels folder, where annotate reads them.
## Runs the app's server on a made-up run of 60 cells in 3 clusters ##
skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("scattermore")

make_run <- function() {
  root <- tempfile("run_")
  res <- file.path(root, "results")
  dir.create(file.path(res, "local_export"), recursive = TRUE)
  dir.create(file.path(res, "tables", "06_annotation", "labels"), recursive = TRUE)
  set.seed(1)
  cl <- rep(c("1", "2", "3"), each = 20)
  cells <- data.frame(
    cell_ID = sprintf("c%02d", 1:60),
    x = runif(60, 0, 1000), y = runif(60, 0, 1000),
    UMAP1 = rnorm(60), UMAP2 = rnorm(60),
    sample_id = rep(c("s1", "s2"), 30), fov = "1",
    leiden_clus_res0.3 = cl,
    cell_type = c(PT = "PT", TAL = "TAL", Unassigned = "Unassigned")[
      c("1" = "PT", "2" = "TAL", "3" = "Unassigned")[cl]
    ]
  )
  saveRDS(cells, file.path(res, "local_export", "cells_umap_clusters.rds"))
  # The pipeline's draft for the clustering, as render.sh leaves it
  write.csv(
    data.frame(cluster_id = c("1", "2", "3"), label = c("PT", "TAL", "Unassigned")),
    file.path(res, "tables", "06_annotation", "labels", "cluster_labels.leiden_clus_res0.3.csv"),
    row.names = FALSE
  )
  # Square outlines around each cell, and the clustering's triage table
  saveRDS(
    data.frame(
      i = rep(seq_len(nrow(cells)), each = 4),
      px = rep(cells$x, each = 4) + c(-8, 8, 8, -8),
      py = rep(cells$y, each = 4) + c(-8, -8, 8, 8)
    ),
    file.path(res, "local_export", "cell_polygons.rds")
  )
  dir.create(file.path(res, "tables", "03_clustering", "markers"), recursive = TRUE)
  write.csv(
    data.frame(cluster_id = c("1", "2", "3"), n_cells = 20L, best_guess = c("PT", "TAL", ""),
               best_score = c(0.9, 0.8, NA), top_markers = ""),
    file.path(res, "tables", "03_clustering", "markers",
              "cluster_labels.leiden_clus_res0.3.template.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(cluster = c("1", "1", "2"), Top = c(1, 2, 1), p.value = 0, FDR = 0,
               logFC = c(1.3, 1.0, 2.1), feats = c("GPX3", "ALDOB", "UMOD"), ranking = c(1, 2, 1)),
    file.path(res, "tables", "03_clustering", "markers", "leiden_clus_res0.3_de_top.csv"),
    row.names = FALSE
  )
  root
}

app_dir <- normalizePath(testthat::test_path("..", "..", "shiny"))
source(testthat::test_path("..", "..", "bin", "label_state.R"))

test_that("Set saves the clustering's labels and names it for cell_type", {
  root <- make_run()
  withr::local_envvar(COSMX_NF = file.path(root, "results"))
  app <- shiny::shinyAppDir(app_dir)
  shiny::testServer(app, {
    session$setInputs(cl_basis = "leiden_clus_res0.3", cl_pick = "1", cl_label = "PT_injured")
    session$setInputs(cl_set = 1)
    labs <- read_cluster_labels(file.path(root, "labels", "cluster_labels.leiden_clus_res0.3.csv"))
    expect_equal(unname(labs[c("1", "2", "3")]), c("PT_injured", "TAL", "Unassigned"))
    expect_equal(readLines(file.path(root, "labels", "cluster_basis.txt")), "leiden_clus_res0.3")
    expect_equal(as.character(live_vec("custom_annotation"))[1], "PT_injured")
    # Discard: nothing was saved at start, so the file goes
    session$setInputs(cl_discard = 1)
    expect_false(file.exists(file.path(root, "labels", "cluster_labels.leiden_clus_res0.3.csv")))
  })
})

test_that("uploaded low-confidence edits change the display label only, and are saved", {
  root <- make_run()
  withr::local_envvar(COSMX_NF = file.path(root, "results"))
  up <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(
      cell_ID = c("c41", "c42", "c01"),
      label = c("PT", "Endothelial", "TAL"),
      confidence = c("low", "low", "high")
    ),
    up,
    row.names = FALSE
  )
  app <- shiny::shinyAppDir(app_dir)
  shiny::testServer(app, {
    session$setInputs(cl_basis = "leiden_clus_res0.3")
    session$setInputs(ann_up = data.frame(name = "up.csv", size = file.size(up), type = "text/csv",
                                          datapath = up))
    saved <- read_label_overrides(file.path(root, "labels", "label_overrides.csv"))
    expect_equal(saved$confidence, c("low", "low", "high"))
    l <- live()
    i <- match(c("c41", "c42", "c01"), l$cell_ID)
    expect_equal(l$cell_type[i], c("Unassigned", "Unassigned", "TAL"))
    expect_equal(l$cell_type_display[i], c("PT", "Endothelial", "TAL"))
    expect_equal(l$label_confidence[i], c("low", "low", "high"))
    expect_equal(as.character(live_vec("custom_confidence"))[i], c("low", "low", "high"))
  })
})

test_that("Use for cell_type names the clustering and saves its labels", {
  root <- make_run()
  withr::local_envvar(COSMX_NF = file.path(root, "results"))
  app <- shiny::shinyAppDir(app_dir)
  shiny::testServer(app, {
    session$setInputs(cl_basis = "leiden_clus_res0.3")
    session$setInputs(cl_use = 1)
    expect_equal(readLines(file.path(root, "labels", "cluster_basis.txt")), "leiden_clus_res0.3")
    expect_true(file.exists(file.path(root, "labels", "cluster_labels.leiden_clus_res0.3.csv")))
  })
})

test_that("Spatial shows one sample at a time, other views every sample", {
  root <- make_run()
  withr::local_envvar(COSMX_NF = file.path(root, "results"))
  app <- shiny::shinyAppDir(app_dir)
  shiny::testServer(app, {
    session$setInputs(view = "Spatial", cby = "custom_annotation", samples = c("s1", "s2"),
                      one_sample = "s1")
    expect_equal(sum(keep_s()), 30L)
    session$setInputs(one_sample = ".all")
    expect_equal(sum(keep_s()), 60L)
    session$setInputs(one_sample = "s2", view = "UMAP")
    expect_equal(sum(keep_s()), 60L)
  })
})

test_that("display colourings slash low-confidence cells; label colourings use Cluster labels", {
  root <- make_run()
  dir.create(file.path(root, "labels"))
  write.csv(
    data.frame(cell_ID = c("c41", "c42"), label = c("PT", "Endothelial"), confidence = "low"),
    file.path(root, "labels", "label_overrides.csv"),
    row.names = FALSE
  )
  withr::local_envvar(COSMX_NF = file.path(root, "results"))
  app <- shiny::shinyAppDir(app_dir)
  shiny::testServer(app, {
    session$setInputs(view = "UMAP", cby = "custom_display", samples = c("s1", "s2"),
                      one_sample = "s1", cl_basis = "leiden_clus_res0.3", mode = "inspect",
                      facet = "none", layout = "abs")
    expect_equal(sum(low_cells("custom_display")), 2L)
    expect_null(low_cells("custom_annotation"))
    is_layer <- function(p, cls) any(vapply(p$layers, function(l) inherits(l$geom, cls), logical(1)))
    expect_true(is_layer(build_plot("umap"), "GeomSlash"))
    session$setInputs(view = "Spatial", one_sample = ".all")
    expect_true(is_layer(build_plot("spatial"), "GeomSegment"))
    # Without low-confidence cells in view, no slashes
    session$setInputs(cby = "custom_annotation")
    expect_false(is_layer(build_plot("spatial"), "GeomSegment"))
    expect_equal(tables_basis(), "leiden_clus_res0.3")
    expect_equal(nrow(tri()), 3L)
  })
})

test_that("without COSMX_NF the app reads the run launched last that has a bundle", {
  ex <- parse(file.path(app_dir, "app.R"))
  i <- which(vapply(ex, function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("latest_results"))
  }, logical(1)))
  eval(ex[[i]])
  runs <- tempfile("Runs_")
  mk <- function(run, bundle_age, log_age = NA) {
    r <- file.path(runs, run, "results")
    dir.create(file.path(r, "local_export"), recursive = TRUE)
    dir.create(file.path(r, "pipeline_info"))
    b <- file.path(r, "local_export", "cells_umap_clusters.rds")
    saveRDS(1, b)
    Sys.setFileTime(b, Sys.time() - bundle_age)
    if (!is.na(log_age)) {
      l <- file.path(r, "pipeline_info", "stage_log.tsv")
      writeLines("stage", l)
      Sys.setFileTime(l, Sys.time() - log_age)
    }
  }
  mk("old", bundle_age = 3600, log_age = 3600)
  mk("recent", bundle_age = 7200, log_age = 60)
  dir.create(file.path(runs, "no_bundle", "results"), recursive = TRUE)
  expect_equal(basename(dirname(latest_results(runs))), "recent")
  expect_null(latest_results(tempfile("empty_")))
})

test_that("figure size scales the plot's text and points, and the plot takes the measured height", {
  root <- make_run()
  withr::local_envvar(COSMX_NF = file.path(root, "results"))
  app <- shiny::shinyAppDir(app_dir)
  shiny::testServer(app, {
    session$setInputs(view = "UMAP", cby = "custom_annotation", samples = c("s1", "s2"),
                      one_sample = "s1", cl_basis = "leiden_clus_res0.3", mode = "inspect",
                      facet = "none", layout = "abs")
    p <- build_plot("umap")
    expect_equal(p$theme$text$size, 13)
    expect_equal(p$layers[[1]]$aes_params$size, 1)
    session$setInputs(fig_scale = 0.5)
    p <- build_plot("umap")
    expect_equal(p$theme$text$size, 6.5)
    expect_equal(p$layers[[1]]$aes_params$size, 0.5)
    # The page measures the plot's card body; the plot is drawn that tall, 220 px at least
    session$setInputs(cosmx_avail = 480)
    expect_equal(avail_h(), 480)
    expect_match(output$ui_umap$html, "height:480px", fixed = TRUE)
    session$setInputs(cosmx_avail = 100)
    expect_equal(avail_h(), 220)
  })
})

test_that("on a label colouring, Top DE lists the clusters carrying the highlighted labels", {
  root <- make_run()
  withr::local_envvar(COSMX_NF = file.path(root, "results"))
  app <- shiny::shinyAppDir(app_dir)
  shiny::testServer(app, {
    session$setInputs(cby = "custom_annotation", cl_basis = "leiden_clus_res0.3")
    expect_match(output$de_hdr, "Highlight labels", fixed = TRUE)
    session$setInputs(hl_lv = "PT")
    expect_equal(tables_levels(), "1")
    expect_match(output$de_hdr, "cluster 1$")
    expect_match(output$de_tbl, "label", fixed = TRUE)
    expect_match(output$de_tbl, "GPX3", fixed = TRUE)
    expect_no_match(output$de_tbl, "UMOD", fixed = TRUE)
  })
})

test_that("the live label columns head the Colour by list", {
  ex <- parse(file.path(app_dir, "app.R"))
  defs <- function(name) Filter(function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(name))
  }, as.list(ex))[[1]]
  env <- new.env()
  env$LIVE_COLS <- c("custom_annotation", "custom_display", "custom_confidence")
  env$GRP <- "sample_id"
  eval(defs("GROUPS"), env)
  eval(defs("grouped_choices"), env)
  g <- env$grouped_choices(c(env$LIVE_COLS, "sample_id", "leiden_clus_res0.3", "cell_type",
                             "cell_type_confidence"))
  expect_equal(names(g)[1], "Labels - live (labels folder)")
  expect_equal(unlist(g[[1]]), env$LIVE_COLS)
  expect_equal(unlist(g[["Annotation - manual"]]), c("cell_type", "cell_type_confidence"))
})
