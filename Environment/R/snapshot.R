# Writes a lockfile of the packages on .libPaths(), the first copy of each, to stdout

# Every library but R's own
libs <- setdiff(normalizePath(.libPaths()), normalizePath(.Library))
lock <- tempfile(fileext = ".lock")
invisible(capture.output(
  renv::lockfile_write(
    renv::lockfile_create(type = "all", libpaths = libs, prompt = FALSE, force = TRUE,
                          project = tempdir()),
    file = lock
  )
))

# A record must match the copy R loads first, and come from a source renv can install again
rec <- renv::lockfile_read(lock)$Packages
ip <- installed.packages(lib.loc = libs)
ip <- ip[!duplicated(ip[, "Package"]), , drop = FALSE]
loaded <- setNames(ip[, "Version"], ip[, "Package"])
for (r in rec) {
  if (!identical(r$Version, unname(loaded[r$Package])))
    message("version differs from the loaded copy: ", r$Package, " ", r$Version,
            " vs ", loaded[r$Package])
  if (!r$Source %in% c("Repository", "Bioconductor", "GitHub", "URL"))
    message("source renv cannot install again: ", r$Package, " (", r$Source, ")")
}
message(length(rec), " packages recorded from ", paste(libs, collapse = ", "),
        "; R ", getRversion())

cat(readLines(lock), sep = "\n")
