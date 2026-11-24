## Gene selection for the PCA and NanoString's per-group method ##
source(testthat::test_path("..", "..", "bin", "gene_selection.R"))

## 200 genes x 300 cells of Poisson counts, genes g1 to g5 far more variable than their mean ##
toy_counts <- function() {
  set.seed(1)
  mu <- rep(c(0.1, 1, 4, 10), each = 50)
  m <- matrix(
    rpois(200 * 300, mu),
    nrow = 200,
    dimnames = list(paste0("g", 1:200), paste0("c", 1:300))
  )
  m[1:5, 1:60] <- m[1:5, 1:60] + rpois(5 * 60, 40)
  m
}
toy_norm <- function() log1p(toy_counts())

test_that("all returns every gene on the panel", {
  nm <- toy_norm()
  expect_identical(select_genes(nm, "all"), rownames(nm))
})

test_that("file keeps the listed genes on the panel, in file order", {
  nm <- toy_norm()
  f <- tempfile(fileext = ".txt")
  writeLines(c("g7", "", "g3", "not_a_gene", "g3"), f)
  expect_message(out <- select_genes(nm, "file", gene_file = f), "1 genes")
  expect_identical(out, c("g7", "g3"))
})

test_that("file stops when no listed gene is on the panel", {
  nm <- toy_norm()
  f <- tempfile(fileext = ".txt")
  writeLines(c("x1", "x2"), f)
  expect_error(select_genes(nm, "file", gene_file = f), "none of the genes")
  expect_error(select_genes(nm, "file"), "hvg_file")
})

test_that("an unknown method stops with the allowed ones", {
  expect_error(select_genes(toy_norm(), "fancy"), "variance, variance_batch, all, file")
})

test_that("variance is scran's top genes, as integrate.R picked them before", {
  skip_if_not_installed("scran")
  nm <- toy_norm()
  ref <- suppressWarnings(scran::getTopHVGs(suppressWarnings(scran::modelGeneVar(nm)), n = 20))
  expect_identical(select_genes(nm, "variance", n_hvgs = 20), ref)
})

test_that("variance_batch blocks on the batch, and falls back without one", {
  skip_if_not_installed("scran")
  nm <- toy_norm()
  b <- rep(c("a", "b"), each = 150)
  ref <- suppressWarnings(scran::getTopHVGs(
    suppressWarnings(scran::modelGeneVar(nm, block = b)),
    n = 20
  ))
  expect_identical(select_genes(nm, "variance_batch", n_hvgs = 20, block = b), ref)
  expect_message(
    out <- select_genes(nm, "variance_batch", n_hvgs = 20, block = rep("a", 300)),
    "falls back"
  )
  expect_identical(out, select_genes(nm, "variance", n_hvgs = 20))
})

test_that("loess genes find the overdispersed genes and skip low-mean ones", {
  counts <- toy_counts()
  out <- loess_genes(counts)
  expect_true(all(paste0("g", 1:5) %in% out))
  low <- rownames(counts)[rowMeans(counts) <= 0.2]
  expect_length(intersect(out, low), 0)
})

test_that("loess genes work on a sparse matrix and return nothing for too few genes", {
  counts <- toy_counts()
  expect_identical(
    loess_genes(Matrix::Matrix(counts, sparse = TRUE)),
    loess_genes(counts)
  )
  expect_identical(loess_genes(counts[1:5, ]), character())
})

test_that("per-group genes combine each large group's loess genes", {
  counts <- toy_counts()
  grp <- c(rep("A", 150), rep("B", 140), rep("tiny", 10))
  out <- per_group_genes(counts, grp, min_cells = 100)
  ref <- sort(unique(c(
    loess_genes(counts[, grp == "A"]),
    loess_genes(counts[, grp == "B"])
  )))
  expect_identical(out, ref)
  expect_identical(per_group_genes(counts, grp, min_cells = 1000), character())
})
