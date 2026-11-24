# Checks the image it runs in against a lockfile, pandoc and Giotto's Python
#   Rscript check_image.R <renv.lock> [giotto_dir]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
lock <- jsonlite::fromJSON(args[[1]], simplifyVector = FALSE)
giotto_dir <- if (length(args) >= 2L) args[[2]] else ""
problems <- character()
problem <- function(...) problems <<- c(problems, paste0(...))

# R and the libraries
if (!identical(as.character(getRversion()), lock$R$Version))
  problem("R is ", getRversion(), "; the lock records ", lock$R$Version)
home <- .libPaths()[startsWith(normalizePath(.libPaths()), normalizePath("~"))]
if (length(home)) problem("a library under your home reaches R: ", paste(home, collapse = ", "))
cat("R       :", as.character(getRversion()), "with", paste(.libPaths(), collapse = ", "), "\n")

# Every package at its recorded version, in the first library that has it
want <- vapply(lock$Packages, function(r) r$Version, "")
names(want) <- vapply(lock$Packages, function(r) r$Package, "")
ip <- installed.packages()
ip <- ip[!duplicated(ip[, "Package"]), , drop = FALSE]
have <- setNames(ip[, "Version"], ip[, "Package"])
absent <- setdiff(names(want), names(have))
if (length(absent)) problem("not installed: ", paste(absent, collapse = ", "))
present <- intersect(names(want), names(have))
other <- present[have[present] != want[present]]
if (length(other))
  problem("other versions: ", paste0(other, " ", have[other], " (lock ", want[other], ")", collapse = ", "))

# Each package loads, in fresh sessions of ten, below R's limit on loaded DLLs
child <- tempfile(fileext = ".R")
writeLines(c(
  "for (p in commandArgs(TRUE)) {",
  "  err <- tryCatch({ suppressPackageStartupMessages(loadNamespace(p)); NULL },",
  "                  error = function(e) conditionMessage(e))",
  "  if (!is.null(err)) cat('FAIL', p, gsub('[[:space:]]+', ' ', err), '\\n')",
  "}"), child)
rscript <- file.path(R.home("bin"), "Rscript")
# Raises R's limit of 100 loaded DLLs, kept under 60% of the open-file limit
fds <- suppressWarnings(as.numeric(system("ulimit -n", intern = TRUE)))
dlls <- if (is.na(fds)) 500 else min(500, floor(0.6 * fds) - 1)
dll_env <- if (dlls > 100) paste0("R_MAX_NUM_DLLS=", dlls) else character()
failed <- character()
for (batch in split(present, ceiling(seq_along(present) / 10))) {
  out <- suppressWarnings(system2(rscript, c(child, batch), stdout = TRUE, stderr = FALSE,
                                  env = dll_env))
  failed <- c(failed, sub("^FAIL ", "", grep("^FAIL ", out, value = TRUE)))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L)
    failed <- c(failed, paste0("a session stopped (status ", status, ") loading one of ",
                               paste(batch, collapse = ", ")))
}
if (length(failed)) problem("do not load:\n    ", paste(failed, collapse = "\n    "))
cat("packages:", length(present) - length(failed), "of the lock's", length(want), "installed and loading\n")

# rmarkdown renders a page, as `cosmx-nf <run> labels` does
if (!rmarkdown::pandoc_available()) {
  problem("rmarkdown finds no pandoc")
} else {
  rmd <- tempfile(fileext = ".Rmd")
  writeLines(c("---", "title: check", "output: html_document", "---", "", "```{r}", "1 + 1", "```"), rmd)
  rendered <- tryCatch({
    rmarkdown::render(rmd, quiet = TRUE)
    TRUE
  }, error = function(e) {
    problem("rmarkdown cannot render: ", conditionMessage(e))
    FALSE
  })
  if (rendered) cat("pandoc  :", as.character(rmarkdown::pandoc_version()), "renders a page\n")
}

# Giotto's Python: its modules import, and Giotto takes it from giotto.py_path as the scripts set it
py <- Sys.getenv("RETICULATE_PYTHON")
mods <- "leidenalg, igraph, community, sklearn, smfishHmrf, umap, numpy, pandas, scipy, matplotlib"
if (!nzchar(py)) {
  problem("RETICULATE_PYTHON is not set in the image")
} else {
  out <- suppressWarnings(system2(py, c("-c", shQuote(paste("import", mods))), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status"))) {
    problem("Python modules do not import: ", paste(out, collapse = " "))
  } else {
    cat("python  :", py, "imports", mods, "\n")
  }
  options(giotto.py_path = py)
  accepted <- tryCatch({
    suppressMessages(GiottoClass::set_giotto_python_path())
    TRUE
  }, error = function(e) {
    problem("Giotto does not accept ", py, ": ", conditionMessage(e))
    FALSE
  })
  if (accepted) cat("giotto  : accepts", py, "\n")
}

if (nzchar(giotto_dir)) {
  g <- tryCatch(suppressMessages(GiottoClass::loadGiotto(giotto_dir)), error = function(e) {
    problem("loadGiotto(", giotto_dir, ") fails: ", conditionMessage(e))
    NULL
  })
  if (!is.null(g)) cat("object  :", giotto_dir, "loads,", nrow(GiottoClass::pDataDT(g)), "cells\n")
}

if (length(problems)) {
  cat("\nPROBLEMS\n", paste0("- ", problems, "\n"), sep = "")
  quit(status = 1)
}
cat("\nThe image is complete for the pipeline\n")
