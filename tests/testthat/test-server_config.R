## The server profile's container settings ##

test_that("each task gets its own cache folder in its container", {
  cfg <- readLines(testthat::test_path("..", "..", "conf", "server.config"))
  at <- grep("XDG_CACHE_HOME", cfg, fixed = TRUE)
  expect_length(at, 1)
  # Single quotes leave the variable for the task's shell to fill in
  expect_equal(trimws(cfg[at]), "' --env XDG_CACHE_HOME=$NXF_TASK_WORKDIR/.cache' +")
  # The option sits among the container's run options
  expect_gt(at, grep("runOptions", cfg, fixed = TRUE))
  expect_lt(at, grep("envWhitelist", cfg, fixed = TRUE))
})
