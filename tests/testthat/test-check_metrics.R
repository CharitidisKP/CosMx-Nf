## Scores the discovery checks report ##
source(testthat::test_path("..", "..", "bin", "check_metrics.R"))

test_that("best Jaccard is 1 for every cluster when only the names differ", {
  ref <- rep(c("1", "2", "3"), times = c(10, 20, 30))
  sub <- rep(c("b", "c", "a"), times = c(10, 20, 30))
  expect_equal(best_jaccard(ref, sub), c("1" = 1, "2" = 1, "3" = 1))
})

test_that("a cluster split in two halves scores 0.5", {
  ref <- rep(c("1", "2"), each = 20)
  sub <- c(rep("x", 10), rep("y", 10), rep("z", 20))
  expect_equal(best_jaccard(ref, sub), c("1" = 0.5, "2" = 1))
})

test_that("two merged clusters score by their share of the merge", {
  ref <- rep(c("1", "2"), times = c(30, 10))
  sub <- rep("m", 40)
  expect_equal(best_jaccard(ref, sub), c("1" = 0.75, "2" = 0.25))
})

test_that("stability rows count stable clusters and their cells", {
  j <- c("1" = 0.9, "2" = 0.8, "3" = 0.3)
  n <- c("1" = 100, "2" = 50, "3" = 50)
  out <- stability_row(j, n)
  expect_equal(out$median_jaccard, 0.8)
  expect_equal(out$stable_clusters, 2L)
  expect_equal(out$cells_in_stable_pct, 75)
})

test_that("the pick is the resolution with the highest median Jaccard", {
  res <- c("0.1", "0.2", "0.3", "0.4")
  expect_identical(pick_resolution(res, c(0.47, 0.756, 0.839, 0.771)), "0.3")
})

test_that("a tie goes to the coarser resolution, in any order", {
  expect_identical(pick_resolution(c("0.5", "0.2", "0.3"), c(0.8, 0.8, 0.7)), "0.2")
})

test_that("resolutions without a score are passed over", {
  expect_identical(pick_resolution(c("0.1", "0.2"), c(NA, 0.5)), "0.2")
  expect_identical(pick_resolution(c("0.1", "0.2"), c(NA_real_, NaN)), NA_character_)
})

test_that("the pick keeps the resolution as written", {
  expect_identical(pick_resolution(c("0.9", "1.0"), c(0.5, 0.6)), "1.0")
})

test_that("the pick needs one score per resolution", {
  expect_error(pick_resolution(c("0.1", "0.2"), 0.5), "one score per resolution")
})

test_that("the same subspace in another basis has cosines of 1", {
  set.seed(1)
  a <- matrix(rnorm(200 * 5), 200, 5)
  rot <- qr.Q(qr(matrix(rnorm(25), 5, 5)))
  expect_equal(subspace_cosines(a, a %*% rot), rep(1, 5), tolerance = 1e-8)
})

test_that("orthogonal subspaces have cosines of 0", {
  a <- diag(6)[, 1:3]
  b <- diag(6)[, 4:6]
  expect_equal(subspace_cosines(a, b), rep(0, 3), tolerance = 1e-12)
})

test_that("the gene set table reports size, overlap and markers kept", {
  sets <- list(run = c("a", "b", "c", "d"), other = c("c", "d", "e"))
  out <- gene_set_table(sets, ref = sets$run, markers = c("a", "e", "z"))
  expect_identical(out$method, c("run", "other"))
  expect_identical(out$n_genes, c(4L, 3L))
  expect_equal(out$jaccard_vs_run, c(1, 0.4))
  expect_equal(out$markers_kept_pct, c(33.3, 33.3))
  expect_identical(out$markers_dropped, c("e z", "a z"))
})

test_that("the dimension estimate is the median per noise level and k over every chunk", {
  per_cell <- function(x, ix, ks, sigmas) outer(x[ix, 1], seq_len(length(ks) * length(sigmas)))
  x <- matrix(c(4, 1, 7, 3, 9, 2, 8), ncol = 1)
  at <- c(1L, 3L, 4L, 6L, 7L)
  sig <- c(both_ends = 2, one_end = 1)
  for (cores in c(1L, 2L, 8L)) {
    out <- dims_estimate(x, at, c(10L, 20L), sig, cores, per_cell = per_cell)
    expect_identical(out$noise, c("both_ends", "both_ends", "one_end", "one_end"))
    expect_identical(out$k, c(10L, 20L, 10L, 20L))
    expect_equal(out$estimate, median(x[at, 1]) * 1:4)
  }
})

test_that("cells whose estimate failed are left out of the median", {
  per_cell <- function(x, ix, ks, sigmas) {
    v <- x[ix, 1]
    v[v == 7] <- NA
    matrix(v, nrow = length(ix), ncol = length(ks) * length(sigmas))
  }
  x <- matrix(c(4, 1, 7, 3, 9, 2, 8), ncol = 1)
  out <- dims_estimate(x, c(1L, 3L, 4L, 6L, 7L), 10L, c(one_end = 1), 2L, per_cell = per_cell)
  expect_equal(out$estimate, median(c(4, 3, 2, 8)))
})

test_that("the dimension estimate stops when a chunk fails", {
  per_cell <- function(x, ix, ks, sigmas) {
    if (3L %in% ix) stop("no neighbours")
    matrix(1, length(ix), 1)
  }
  x <- matrix(1:7, ncol = 1)
  expect_error(dims_estimate(x, 1:6, 10L, c(one_end = 1), 1L, per_cell = per_cell), "no neighbours")
  expect_error(
    suppressWarnings(dims_estimate(x, 1:6, 10L, c(one_end = 1), 2L, per_cell = per_cell)),
    "1 of 2 chunks failed"
  )
})

test_that("cell estimates from one neighbour search match the package's pointwise estimates", {
  skip_if_not_installed("intrinsicDimension")
  set.seed(1)
  x <- cbind(matrix(rnorm(300 * 3), 300), matrix(rnorm(300 * 5, sd = 0.05), 300))
  at <- sort(sample(300, 40))
  ks <- c(5L, 8L)
  sig <- c(both_ends = sqrt(2) * 0.05, one_end = 0.05)
  ref <- do.call(cbind, lapply(sig, function(s) {
    sapply(ks, function(k) {
      intrinsicDimension::maxLikPointwiseDimEst(
        x, k = k, dnoise = intrinsicDimension::dnoiseGaussH, sigma = s,
        n = ncol(x), indices = at, unbiased = TRUE
      )$dim.est
    })
  }))
  expect_equal(cell_dims(x, at, ks, sig), unname(ref))
  expect_equal(
    dims_estimate(x, at, ks, sig, 2L)$estimate,
    unname(apply(ref, 2, stats::median, na.rm = TRUE))
  )
})
