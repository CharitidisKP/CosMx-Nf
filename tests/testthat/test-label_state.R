source(testthat::test_path("..", "..", "bin", "label_state.R"))

ov <- function(...) check_overrides(data.frame(..., stringsAsFactors = FALSE))

cells <- paste0("c", 1:6)
clusters <- c("1", "1", "2", "2", "3.10", "3.1")
labels <- c("1" = "PT", "2" = "Unassigned", "3.1" = "TAL", "3.10" = "DCT")

test_that("cluster labels alone set every label, with high confidence", {
  s <- apply_label_state(cells, clusters, labels, empty_overrides())
  expect_equal(s$cell_type, c("PT", "PT", "Unassigned", "Unassigned", "DCT", "TAL"))
  expect_equal(s$cell_type_display, s$cell_type)
  expect_true(all(s$label_confidence == "high"))
  expect_true(all(s$cell_type_source == "cluster"))
})

test_that("a high edit sets cell_type, a low edit only the display label", {
  e <- ov(
    cell_ID = c("c1", "c3"),
    label = c("PT_injured", "Endothelial"),
    confidence = c("high", "low")
  )
  s <- apply_label_state(cells, clusters, labels, e)
  expect_equal(s$cell_type[1:3], c("PT_injured", "PT", "Unassigned"))
  expect_equal(s$cell_type_display[1:3], c("PT_injured", "PT", "Endothelial"))
  expect_equal(s$label_confidence[1:3], c("high", "high", "low"))
  expect_equal(s$cell_type_source[1:3], c("manual", "cluster", "cluster"))
})

test_that("a forced guess after a high relabel keeps the relabel in cell_type", {
  # As in label_overrides_forced.csv: plasma cells set Unassigned, then given a guess
  e <- ov(
    cell_ID = c("c1", "c1"),
    label = c("Unassigned", "PT"),
    confidence = c("high", "low")
  )
  s <- apply_label_state(cells, clusters, labels, e)
  expect_equal(s$cell_type[1], "Unassigned")
  expect_equal(s$cell_type_display[1], "PT")
  expect_equal(s$label_confidence[1], "low")
})

test_that("a later high edit replaces an earlier low guess", {
  e <- ov(
    cell_ID = c("c3", "c3"),
    label = c("Endothelial", "Podocyte"),
    confidence = c("low", "high")
  )
  s <- apply_label_state(cells, clusters, labels, e)
  expect_equal(s$cell_type[3], "Podocyte")
  expect_equal(s$cell_type_display[3], "Podocyte")
  expect_equal(s$label_confidence[3], "high")
})

test_that("(cluster label) undoes a cell's edits", {
  e <- ov(
    cell_ID = c("c1", "c1", "c3", "c3"),
    label = c("TAL", KEEP_LABEL, "PT", KEEP_LABEL),
    confidence = c("high", "high", "low", "low")
  )
  s <- apply_label_state(cells, clusters, labels, e)
  expect_equal(s$cell_type[c(1, 3)], c("PT", "Unassigned"))
  expect_equal(s$cell_type_display[c(1, 3)], c("PT", "Unassigned"))
  expect_equal(s$label_confidence[c(1, 3)], c("high", "high"))
})

test_that("the last of several edits wins", {
  e <- ov(cell_ID = rep("c2", 3), label = c("TAL", "DCT", "PEC"))
  s <- apply_label_state(cells, clusters, labels, e)
  expect_equal(s$cell_type[2], "PEC")
})

test_that("cluster ids are matched as text, so 3.1 and 3.10 stay apart", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("cluster_id,label", "3.1,TAL", "3.10,DCT"), f)
  m <- read_cluster_labels(f)
  expect_equal(unname(m[c("3.1", "3.10")]), c("TAL", "DCT"))
})

test_that("duplicate cluster ids and bad confidence values stop", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("cluster_id,label", "1,PT", "1,TAL"), f)
  expect_error(read_cluster_labels(f), "duplicate cluster_id")
  expect_error(ov(cell_ID = "c1", label = "PT", confidence = "maybe"), "high or low")
})

test_that("blank confidence reads as high and blank labels are dropped", {
  e <- ov(cell_ID = c("c1", "c2"), label = c("PT", ""), confidence = c("", "low"))
  expect_equal(nrow(e), 1L)
  expect_equal(e$confidence, "high")
})

test_that("edits of cells not in the object are ignored", {
  e <- ov(cell_ID = "elsewhere", label = "PT")
  s <- apply_label_state(cells, clusters, labels, e)
  expect_equal(s$cell_type, unname(labels[clusters]))
})

test_that("compacting edits never changes the labels they give", {
  set.seed(1)
  for (k in 1:300) {
    n <- sample(1:25, 1)
    e <- ov(
      cell_ID = sample(cells, n, replace = TRUE),
      label = sample(c("PT", "TAL", "Endothelial", KEEP_LABEL), n, replace = TRUE),
      confidence = sample(c("high", "low"), n, replace = TRUE)
    )
    full <- apply_label_state(cells, clusters, labels, e)
    cut <- compact_overrides(e)
    expect_lte(nrow(cut), nrow(e))
    expect_false(any(cut$label == KEEP_LABEL))
    expect_equal(apply_label_state(cells, clusters, labels, cut), full)
  }
})

test_that("a file of edits round-trips through write_csv_atomic", {
  e <- ov(
    cell_ID = c("c1", "c3"),
    label = c("PT", "TAL"),
    confidence = c("high", "low"),
    note = c("", "a, note")
  )
  f <- tempfile(fileext = ".csv")
  write_csv_atomic(e, f)
  expect_equal(read_label_overrides(f), e)
  expect_equal(nrow(read_label_overrides(tempfile())), 0L)
})
