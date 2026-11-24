# Builds the KPMP profile: mean expression per author_cell_type with at least 200 cells
#   Rscript kpmp_profile.R <KPMP atlas .h5ad> <KPMP.rds>
suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
H5 <- args[[1]]
out <- args[[2]]
LABEL <- "author_cell_type"
MIN_CELLS <- 200L

# anndata stores a categorical obs column as codes plus categories
h5col <- function(h5, path) {
  x <- h5read(h5, path)
  if (is.list(x) && all(c("codes", "categories") %in% names(x)))
    factor(x$categories[x$codes + 1L], levels = x$categories) else as.vector(x)
}

GRP <- if (any(h5ls(H5)$group == "/raw")) "raw/X" else "X"
VAR <- if (GRP == "raw/X") "raw/var" else "var"

sym    <- as.character(h5col(H5, paste0(VAR, "/feature_name")))
indptr <- as.numeric(h5read(H5, paste0(GRP, "/indptr")))
n_obs  <- length(indptr) - 1L
n_var  <- length(sym)
f      <- h5col(H5, paste0("obs/", LABEL))
S      <- matrix(0, n_var, nlevels(f), dimnames = list(sym, levels(f)))

# Counts summed per cell type, 20,000 cells at a time
CH <- 20000L
for (s in seq(1L, n_obs, by = CH)) {
  e  <- min(s + CH - 1L, n_obs)
  lo <- indptr[s] + 1; hi <- indptr[e + 1L]
  if (hi < lo) next
  d  <- as.numeric(h5read(H5, paste0(GRP, "/data"),    index = list(lo:hi)))
  j  <- as.integer(h5read(H5, paste0(GRP, "/indices"), index = list(lo:hi)))
  Xt <- sparseMatrix(i = j, p = indptr[s:(e + 1L)] - indptr[s], x = d,
                     dims = c(n_var, e - s + 1L), index1 = FALSE)
  ff <- f[s:e]; ok <- !is.na(ff)
  I  <- sparseMatrix(i = which(ok), j = as.integer(ff[ok]), x = 1,
                     dims = c(e - s + 1L, nlevels(f)))
  S  <- S + as.matrix(Xt %*% I)
}

# Mean per cell type, summed over duplicate genes, small cell types dropped
n <- as.integer(table(f))
p <- rowsum(sweep(S, 2, n, "/"), group = sym)
p <- p[, n >= MIN_CELLS, drop = FALSE]
saveRDS(p, out)
message(nrow(p), " genes x ", ncol(p), " cell types written to ", out)
