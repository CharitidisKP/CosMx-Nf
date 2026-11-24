## Gene selection for the PCA, sourced by integrate.R and cluster_checks.R ##

HVG_METHODS <- c("variance", "variance_batch", "all", "file")

## Genes for the PCA from a genes x cells normalised matrix, by one of HVG_METHODS ##
select_genes <- function(nm, method = "variance", n_hvgs = 4000L, block = NULL, gene_file = "") {
  if (!method %in% HVG_METHODS) {
    stop("hvg_method must be one of: ", paste(HVG_METHODS, collapse = ", "), call. = FALSE)
  }
  if (method == "variance_batch" && length(unique(block)) < 2L) {
    message("gene_selection: no batch with 2 or more levels, variance_batch falls back to variance")
    method <- "variance"
  }
  switch(
    method,
    variance = suppressWarnings(
      scran::getTopHVGs(suppressWarnings(scran::modelGeneVar(nm)), n = n_hvgs)
    ),
    variance_batch = suppressWarnings(
      scran::getTopHVGs(suppressWarnings(scran::modelGeneVar(nm, block = block)), n = n_hvgs)
    ),
    all = rownames(nm),
    file = file_genes(gene_file, rownames(nm))
  )
}

## The genes listed in a file, one per line, that are on the panel ##
file_genes <- function(gene_file, panel) {
  if (!nzchar(gene_file)) stop("hvg_method file needs hvg_file", call. = FALSE)
  genes <- unique(trimws(readLines(gene_file, warn = FALSE)))
  genes <- genes[nzchar(genes)]
  miss <- setdiff(genes, panel)
  if (length(miss)) {
    message("gene_selection: ", length(miss), " genes in ", basename(gene_file), " are not on the panel")
  }
  keep <- intersect(genes, panel)
  if (!length(keep)) stop("none of the genes in ", gene_file, " are on the panel", call. = FALSE)
  keep
}

## NanoString's subclustering genes: variance above a loess fit on the mean, from raw counts ##
loess_genes <- function(counts, var_ratio = 1, min_mean = 0.2, span = 0.3) {
  n <- ncol(counts)
  mu <- Matrix::rowMeans(counts)
  v <- (Matrix::rowSums(counts^2) - n * mu^2) / (n - 1)
  ok <- mu > 0 & v > 0
  if (sum(ok) < 10L) return(character())
  x <- log10(mu[ok])
  y <- log10(v[ok])
  fit <- stats::loess(y ~ x, span = span)
  ratio <- 10^(y - stats::predict(fit))
  names(ratio)[ratio > var_ratio & mu[ok] > min_mean]
}

## The loess genes of every group with at least min_cells cells, combined ##
per_group_genes <- function(counts, groups, min_cells = 500L, ...) {
  groups <- as.character(groups)
  big <- names(which(table(groups) >= min_cells))
  as.character(sort(unique(unlist(lapply(big, function(g) {
    loess_genes(counts[, groups == g, drop = FALSE], ...)
  })))))
}
