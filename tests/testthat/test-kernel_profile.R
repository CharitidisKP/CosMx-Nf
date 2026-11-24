## The Positron kernel's R profile ##

kernel_profile <- function() {
  normalizePath(testthat::test_path("..", "..", "notebooks", "kernel.Rprofile"))
}

test_that("the kernel profile sets Giotto's Python, then reads the folder's profile", {
  prof <- kernel_profile()
  dir <- withr::local_tempdir()
  writeLines("kernel_profile_read <- TRUE", file.path(dir, ".Rprofile"))
  withr::local_envvar(RETICULATE_PYTHON = "/opt/py/bin/python")
  withr::local_options(giotto.py_path = NULL)
  withr::defer(suppressWarnings(rm("kernel_profile_read", envir = globalenv())))
  withr::with_dir(dir, source(prof))
  expect_equal(getOption("giotto.py_path"), "/opt/py/bin/python")
  expect_true(exists("kernel_profile_read", envir = globalenv(), inherits = FALSE))
})

test_that("without RETICULATE_PYTHON the kernel profile leaves Giotto's option alone", {
  prof <- kernel_profile()
  dir <- withr::local_tempdir()
  # An empty profile in the folder, read in place of the home one
  writeLines("", file.path(dir, ".Rprofile"))
  withr::local_envvar(RETICULATE_PYTHON = NA)
  withr::local_options(giotto.py_path = "/kept")
  withr::with_dir(dir, source(prof))
  expect_equal(getOption("giotto.py_path"), "/kept")
})
