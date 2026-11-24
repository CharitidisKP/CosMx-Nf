# Installs every package in renv.lock, at its locked version, into one library
#   Rscript install.R <renv.lock> <library>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
lockfile <- args[[1]]
library_dir <- args[[2]]

# No renv cache, and each package stays in the library as soon as it installs
Sys.setenv(RENV_CONFIG_CACHE_ENABLED = "FALSE", RENV_CONFIG_INSTALL_TRANSACTIONAL = "FALSE",
           RENV_PATHS_ROOT = tempfile("renv-root-"))
# At most 8 parallel installs, so a build on a shared machine leaves room for other jobs
options(Ncpus = max(1L, min(8L, parallel::detectCores() - 1L)))

# renv only does the installing here, so whichever version the base image's repository has
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# CRAN repositories, a dated P3M snapshot for binaries and CRAN for sources
binaries <- c(CRAN = "https://p3m.dev/cran/__linux__/noble/2026-10-25")
sources <- c(CRAN = "https://cloud.r-project.org")

# A copy of the lock with each repository package's Repository set to CRAN
lock <- renv::lockfile_read(lockfile)
lock$Packages <- lapply(lock$Packages, function(r) {
  if (identical(r$Source, "Repository")) r$Repository <- "CRAN"
  r
})
lockfile <- tempfile(fileext = ".lock")
invisible(capture.output(renv::lockfile_write(lock, file = lockfile)))
records <- lock$Packages
needs <- function(r, fields = c("Depends", "Imports", "LinkingTo"))
  intersect(trimws(sub("\\(.*", "", unlist(r[intersect(fields, names(r))]))), names(records))
# The packages plus everything in the lock that needs them, directly or not
with_dependents <- function(pkgs) {
  repeat {
    more <- names(Filter(function(r) any(needs(r) %in% pkgs), records))
    if (all(more %in% pkgs)) return(pkgs)
    pkgs <- union(pkgs, more)
  }
}

# Built from source: packages linking to RcppParallel and their dependents
from_source <- with_dependents(
  names(Filter(function(r) "RcppParallel" %in% needs(r, "LinkingTo"), records)))

project <- tempfile("renv-project-")
dir.create(project)
restore <- function(...)
  renv::restore(project = project, lockfile = lockfile, library = library_dir,
                prompt = FALSE, ...)

done <- FALSE
for (pass in 1:10) {
  failed <- tryCatch({
    restore(repos = binaries, exclude = from_source)
    NULL
  }, error = function(e) {
    hit <- regmatches(conditionMessage(e),
                      regexec("testing if '([^']+)' can be loaded", conditionMessage(e)))[[1]]
    if (length(hit) < 2L || hit[[2]] %in% from_source) stop(e)
    hit[[2]]
  })
  if (is.null(failed)) { done <- TRUE; break }
  message("The binary of ", failed, " does not load; it will be built from source")
  from_source <- with_dependents(union(from_source, failed))
}
if (!done) stop("binaries still failing to load after 10 passes")

if (length(from_source)) {
  restore(repos = sources, packages = from_source)
  message("Built from source against the locked versions: ",
          paste(sort(from_source), collapse = ", "))
}
