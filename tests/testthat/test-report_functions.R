## The run report's helpers that need no results folder ##
source(testthat::test_path("..", "..", "report", "_functions.R"))

test_that("each setting gets its own comment, without the group headers", {
  lines <- c(
    "nextflow.enable.dsl = 2",
    "params {",
    "    input                       = \"x\"",
    "",
    "    // Quality control parameters",
    "    // Giotto vignette",
    "    gene_min_cells              = 5",
    "    // Effective floor is min(cap, quantile of the sample)",
    "    count_cap                   = 200",
    "",
    "    // Dimensionality reduction",
    "    scalefactor                 = \"auto\"",
    "",
    "    // Cell-cell communication",
    "    // LIANA",
    "    // Cells a cell type needs for ligand-receptor scoring",
    "    liana_min_cells             = 10",
    "",
    "    // Samples run through discovery on their own, a comma list of sample_id or all",
    "    solo                        = null",
    "    // Per-cell label edits, replacing labels_dir's label_overrides.csv",
    "    // Columns cell_ID, label, confidence and note",
    "    label_overrides             = null",
    "",
    "    // Spatial domains with BANKSY, each sample on its own",
    "    // Settings sets, one per row",
    "    //   variant      the set's name",
    "    banksy_variants             = \"${projectDir}/assets/banksy_variants.csv\"",
    "}",
    "process {",
    "    cpus = 2",
    "}"
  )
  d <- config_params(lines)
  doc <- setNames(d$doc, d$name)
  expect_equal(doc[["gene_min_cells"]], "Giotto vignette")
  expect_equal(doc[["count_cap"]], "Effective floor is min(cap, quantile of the sample)")
  expect_equal(doc[["scalefactor"]], "")
  expect_equal(doc[["liana_min_cells"]], "Cells a cell type needs for ligand-receptor scoring")
  expect_equal(doc[["solo"]], "Samples run through discovery on their own, a comma list of sample_id or all")
  expect_equal(
    doc[["label_overrides"]],
    "Per-cell label edits, replacing labels_dir's label_overrides.csv. Columns cell_ID, label, confidence and note"
  )
  expect_equal(doc[["banksy_variants"]], "Settings sets, one per row")
  expect_equal(doc[["input"]], "")
  expect_false("cpus" %in% d$name)
})

test_that("each setting's default is its value without quotes", {
  d <- config_params(c(
    "params {",
    "    count_cap = 200",
    "    scalefactor = \"auto\"",
    "    hvg_file = null",
    "}"
  ))
  expect_equal(setNames(d$default, d$name), c(count_cap = "200", scalefactor = "auto", hvg_file = "null"))
})

test_that("tables over the cap keep the rows with the smallest p", {
  d <- data.frame(g = letters[1:6], p = c(0.5, 0.01, NA, 0.2, 0.001, 0.05))
  r <- cap_rows(d, 3, "p")
  expect_equal(r$df$g, c("e", "b", "f"))
  expect_equal(r$total, 6L)
  expect_equal(r$by, "p")
})

test_that("a minus sign caps a table from its largest values", {
  d <- data.frame(g = letters[1:4], n = c(5, 50, 1, 20))
  expect_equal(cap_rows(d, 2, "-n")$df$g, c("b", "d"))
})

test_that("tables within the cap stay whole and in their order", {
  d <- data.frame(g = letters[1:3], p = c(0.3, 0.1, 0.2))
  r <- cap_rows(d, 5, "p")
  expect_equal(r$df$g, c("a", "b", "c"))
  expect_null(r$by)
})

test_that("without an order column the cap keeps the first rows", {
  d <- data.frame(g = letters[1:4])
  r <- cap_rows(d, 2, NULL)
  expect_equal(r$df$g, c("a", "b"))
  expect_null(r$by)
})

test_that("p-value columns are told apart from other statistics", {
  yes <- c(
    "p.value", "FDR", "PValue", "fdr", "p_value", "p.adj_higher", "p_lower_orig",
    "adjustedPval", "combinedPval", "pval", "padj", "wilcox_p", "p_het"
  )
  no <- c("Q", "df_Q", "logFC", "PI_value", "top_sample_pct", "prop_sv", "phi", "aupr", "pearson", "n_cells")
  expect_true(all(is_p_col(yes)))
  expect_false(any(is_p_col(no)))
})

test_that("links from the site go up one folder and encode each part of the path", {
  expect_equal(
    rel_href("/r/results/figures/03 clustering/a#1.png", "/r/results"),
    "../figures/03%20clustering/a%231.png"
  )
})

test_that("a PNG's size is read from its header", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 123, height = 45)
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  grDevices::dev.off()
  expect_equal(png_size(f), c(123L, 45L))
  g <- tempfile(fileext = ".png")
  writeLines("not a png", g)
  expect_equal(png_size(g), c(NA_integer_, NA_integer_))
})

test_that("the history gives each of the run's launches with its stage", {
  h <- c(
    paste(
      "2026-11-03 13:31:12", "5h 54m", "compassionate_albattani", "OK", "b69ccc31fa", "sess-1",
      "nextflow run main.nf --outdir /x --run_id r1 --stage checks --checks_stability true",
      sep = "\t"
    ),
    paste("2026-11-02 10:00:00", "1h", "other", "ERR", "abc", "sess-2", "nextflow run main.nf --run_id r10 --stage spatial", sep = "\t"),
    paste("2026-11-01 09:00:00", "3h", "first", "OK", "abc", "sess-3", "nextflow run main.nf --run_id r1", sep = "\t"),
    "a line that is not a launch"
  )
  d <- history_launches(h, "r1")
  expect_equal(d$stage, c("checks", "discover"))
  expect_equal(d$status, c("OK", "OK"))
  expect_equal(d$session, c("sess-1", "sess-3"))
  expect_equal(nrow(history_launches(character(), "r1")), 0L)
})

test_that("a stage's state comes from its last launch", {
  d <- data.frame(
    stage = c("discover", "spatial", "spatial", "cci"),
    status = c("OK", "OK", "ERR", "-"),
    stringsAsFactors = FALSE
  )
  expect_equal(launch_state("discover", d), "done")
  expect_equal(launch_state("spatial", d), "failed")
  expect_equal(launch_state("cci", d, current = "cci"), "running")
  expect_equal(launch_state("cci", d, current = "de"), "stopped")
  expect_equal(launch_state("de", d, made = TRUE), "done")
  expect_equal(launch_state("de", d, made = FALSE), "none")
})

test_that("a stage is outdated when nothing it wrote is newer than the object it read", {
  expect_true(outdated_by(10, c(5, 9)))
  expect_false(outdated_by(10, c(5, 11)))
  expect_false(outdated_by(NA, 5))
  expect_false(outdated_by(10, numeric()))
})

test_that("the label digest matches main.nf and the launcher", {
  d <- tempfile()
  dir.create(d)
  writeLines(c("cluster_id,label", "1,PT", "2,TAL"), file.path(d, "l.csv"))
  writeLines(c("cell_ID,label,confidence,note", "c1,Podocyte,low,checked"), file.path(d, "o.csv"))
  expect_equal(
    label_digest("leiden_clus_res0.3", file.path(d, "l.csv"), file.path(d, "o.csv")),
    "7f2e6431922f3a99a226014ac53b08bb"
  )
})

test_that("a launch's options leave out what cosmx-nf adds", {
  cmd <- paste(
    "nextflow -log /x/logs/nextflow.log run main.nf -work-dir /tmp/w --outdir /x/results",
    "--run_id r1 --stage checks -params-file /p/site.yaml",
    "-resume 6de67345-cc42-4582-94c2-633ecbce730e --checks_stability true -profile server"
  )
  expect_equal(stage_options(cmd), "--checks_stability true -profile server")
  expect_equal(stage_options("nextflow run main.nf --run_id r1 -resume"), "")
})

test_that("a figure's key changes with its inputs, its code and its size, and only then", {
  f <- tempfile(fileext = ".csv")
  writeLines("a,b", f)
  k <- fig_key("qc", f, "ggplot()", 7, 4)
  expect_equal(fig_key("qc", f, "ggplot()", 7, 4), k)
  expect_false(fig_key("qc", f, "ggplot() + geom_point()", 7, 4) == k)
  expect_false(fig_key("qc", f, "ggplot()", 8, 4) == k)
  expect_false(fig_key("other", f, "ggplot()", 7, 4) == k)
  writeLines(c("a,b", "1,2"), f)
  expect_false(fig_key("qc", f, "ggplot()", 7, 4) == k)
})

test_that("text columns that hold only numbers come back as numbers, unless asked for as text", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("fold_change,gene,id,lambda", "0.368,A,1234567890123456789012,0.8",
               "2300052571982320000,B,2,1.0", ",C,3,0.8"), f)
  raw <- data.table::fread(f, colClasses = c(lambda = "character"))
  expect_type(raw$fold_change, "character")
  out <- as_numbers(data.table::copy(raw), c(lambda = "character"))
  expect_type(out$fold_change, "double")
  expect_equal(out$fold_change, c(0.368, 2300052571982320000, NA))
  expect_type(out$gene, "character")
  expect_type(out$id, "character")
  expect_type(out$lambda, "character")
  expect_type(as_numbers(data.table::copy(raw))$lambda, "double")
  expect_type(as_numbers(data.table::copy(raw), list(character = "lambda"))$lambda, "character")
  expect_type(as_numbers(data.table::copy(raw), "character")$fold_change, "character")
  expect_type(as_numbers(as.data.frame(raw))$fold_change, "double")
})

test_that("git keeps every file the report site is built from", {
  root <- testthat::test_path("..", "..")
  skip_if_not(nzchar(Sys.which("git")) && file.exists(file.path(root, ".git")))
  files <- list.files(file.path(root, "report"), recursive = TRUE)
  files <- files[!grepl("^(_site|\\.quarto|_freeze|site_libs)/|_files/", files)]
  ignored <- suppressWarnings(system2(
    "git", c("-C", root, "check-ignore", "--", file.path("report", files)), stdout = TRUE, stderr = FALSE
  ))
  expect_length(ignored, 0)
})

test_that("the report step prints its warnings in the console", {
  nf <- readLines(testthat::test_path("..", "..", "modules", "report.nf"))
  expect_true(any(grepl("^\\s*debug\\s+true\\s*$", nf)))
  warn <- grep("REPORT: ", nf, value = TRUE, fixed = TRUE)
  expect_gte(length(warn), 5)
  expect_false(any(grepl(">&2", warn, fixed = TRUE)))
})

test_that("the report step writes the HTML report and nothing else", {
  nf <- readLines(testthat::test_path("..", "..", "modules", "report.nf"))
  expect_false(any(grepl("run_report", nf, fixed = TRUE)))
  outputs <- grep("^\\s*path \"", nf, value = TRUE)
  expect_equal(trimws(outputs), "path \"report\", optional: true")
})

test_that("cluster IDs sort by number at each level, names last", {
  expect_equal(cluster_levels(c("6.10", "6.2", "6.1", "10.1", "2.3")), c("2.3", "6.1", "6.2", "6.10", "10.1"))
  expect_equal(cluster_levels(c("10", "2", "1", NA, "2")), c("1", "2", "10"))
  expect_equal(cluster_levels(c("b", "3", "a")), c("3", "a", "b"))
  expect_equal(cluster_levels(c("5.5.10", "5.5.2", "5.4", "5.5.1", "10.1", "5.1.8", "2.1")), c("2.1", "5.1.8", "5.4", "5.5.1", "5.5.2", "5.5.10", "10.1"))
})

test_that("legends get a column per 22 entries and none past 66", {
  expect_equal(vapply(c(5, 22, 23, 44, 45, 66, 67, 100), legend_columns, 1L), c(1L, 1L, 2L, 2L, 3L, 3L, 0L, 0L))
})

test_that("a header of numbers is still read as the header", {
  f <- tempfile(fileext = ".csv")
  writeLines(c('"","1.1","1.10","2.1"', '"PT",0.5,0.2,0.1'), f)
  expect_true(startsWith(names(data.table::fread(f))[1], "V"))
  expect_equal(names(data.table::fread(f, header = TRUE))[-1], c("1.1", "1.10", "2.1"))
})
