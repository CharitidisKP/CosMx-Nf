## The run report's figure helpers, loaded on an empty results folder ##
figures_env <- function(res, dir = testthat::test_path("..", "..", "report")) {
  withr::local_envvar(COSMX_RESULTS = res, COSMX_CURRENT = "", COSMX_FIG_STYLE = "", COSMX_LABEL_COLOURS = "")
  withr::local_dir(dir)
  e <- new.env()
  sys.source("_helpers.R", envir = e)
  sys.source("_figures.R", envir = e)
  e
}

test_that("a figure that fails to draw leaves no file, and one with nothing left to draw takes its old file away", {
  skip_if_not_installed("ragg")
  res <- withr::local_tempdir()
  e <- figures_env(res)
  input <- file.path(res, "genes.csv")
  writeLines("gene", input)
  figs <- function() list.files(file.path(res, "report_figures"), all.files = TRUE, no.. = TRUE)
  no_rows <- function() ggplot2::ggplot(data.frame(x = numeric()), ggplot2::aes(x, x)) + ggplot2::facet_wrap(~x)
  expect_error(e$rfig("volcano", input, no_rows), "at least one value")
  expect_length(figs(), 0)
  drawn <- e$rfig("volcano", input, function() ggplot2::ggplot(data.frame(x = 1), ggplot2::aes(x, x)) + ggplot2::geom_point())
  expect_equal(figs(), basename(drawn))
  expect_length(e$rfig("volcano", input, function() NULL), 0)
  expect_length(figs(), 0)
})

test_that("a change to the shared helpers draws every figure again", {
  src <- withr::local_tempdir()
  file.copy(file.path(testthat::test_path("..", "..", "report"), c("_helpers.R", "_functions.R", "_figures.R")), src)
  before <- figures_env(withr::local_tempdir(), src)$FIG_CODE
  cat("\n## A new helper ##\n", file = file.path(src, "_functions.R"), append = TRUE)
  expect_false(identical(figures_env(withr::local_tempdir(), src)$FIG_CODE, before))
})

test_that("a cluster heatmap puts its clusters down the side once they no longer fit across", {
  e <- figures_env(withr::local_tempdir())
  expect_equal(unname(e$cluster_heat_size(50, 41)), c(2.6 + 0.24 * 50, 1.6 + 0.2 * 41, 0))
  expect_equal(unname(e$cluster_heat_size(134, 41)), c(2.6 + 0.24 * 41, 1.6 + 0.2 * 134, 1))
  expect_named(e$cluster_heat_size(134, 41), c("w", "h", "down"))
})

test_that("a volcano leaves out genes without an FDR and wraps its panels up to three to a row", {
  e <- figures_env(withr::local_tempdir())
  d <- data.frame(
    gene = paste0("g", 1:20), cell_type = rep(c("PT", "TAL", "DCT", "CNT", "PC"), each = 4),
    logFC = c(2, -1, 0.5, 0.1), FDR = NA_real_
  )
  expect_null(e$volcano_plot(d, "logFC", "FDR", gene = "gene", facet = "cell_type"))
  d$FDR <- c(0.001, 0.02, 0.4, 0.9)
  d$FDR[d$cell_type == "PC"] <- NA
  layout <- function(d) ggplot2::ggplot_build(e$volcano_plot(d, "logFC", "FDR", gene = "gene", facet = "cell_type"))$layout
  four <- layout(d)
  expect_equal(nrow(four$layout), 4)
  expect_equal(c(max(four$layout$ROW), four$facet_params$ncol), c(2, 3))
  d$FDR[d$cell_type %in% c("DCT", "CNT")] <- NA
  expect_equal(layout(d)$facet_params$ncol, 2)
})
