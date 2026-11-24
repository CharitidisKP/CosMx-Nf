## Scores for the discovery checks, sourced by cluster_checks.R ##

## Each cluster's best Jaccard overlap with a cluster from re-clustering the same cells ##
best_jaccard <- function(ref, sub) {
  inter <- unclass(table(as.character(ref), as.character(sub)))
  uni <- outer(rowSums(inter), colSums(inter), "+") - inter
  apply(inter / uni, 1, max)
}

## Median best Jaccard, stable clusters and the share of cells in them, for one resolution ##
stability_row <- function(j, n_cells, cut = 0.75) {
  n_cells <- n_cells[names(j)]
  data.frame(
    median_jaccard = round(stats::median(j), 3),
    stable_clusters = sum(j >= cut),
    cells_in_stable_pct = round(100 * sum(n_cells[j >= cut]) / sum(n_cells), 1)
  )
}

## The resolution with the highest median Jaccard, the coarser one on a tie ##
pick_resolution <- function(resolution, median_jaccard) {
  if (length(resolution) != length(median_jaccard)) {
    stop("pick_resolution: needs one score per resolution", call. = FALSE)
  }
  ok <- is.finite(median_jaccard)
  if (!any(ok)) {
    return(NA_character_)
  }
  best <- resolution[ok][median_jaccard[ok] == max(median_jaccard[ok])]
  as.character(best[which.min(as.numeric(best))])
}

## Cosines of the principal angles between two embeddings of the same cells ##
subspace_cosines <- function(a, b, k = min(ncol(a), ncol(b))) {
  qa <- qr.Q(qr(a[, seq_len(k), drop = FALSE]))
  qb <- qr.Q(qr(b[, seq_len(k), drop = FALSE]))
  pmin(1, svd(crossprod(qa, qb))$d)
}

## One row per gene set: its size, overlap with the run's set and the markers it keeps ##
gene_set_table <- function(sets, ref, markers) {
  data.frame(
    method = names(sets),
    n_genes = vapply(sets, length, integer(1)),
    jaccard_vs_run = vapply(
      sets,
      function(s) round(length(intersect(s, ref)) / length(union(s, ref)), 3),
      numeric(1)
    ),
    markers_kept_pct = vapply(sets, function(s) round(100 * mean(markers %in% s), 1), numeric(1)),
    markers_dropped = vapply(sets, function(s) paste(setdiff(markers, s), collapse = " "), character(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

## Noisy ML dimension of each cell in ix for every k and noise level, from one search of its nearest cells ##
cell_dims <- function(x, ix, ks, sigmas) {
  kmax <- max(ks)
  grid <- expand.grid(k = ks, s = sigmas)
  d <- intrinsicDimension:::yaIkNN(x, ix, kmax)[, kmax + seq_len(kmax), drop = FALSE]
  est <- vapply(seq_len(nrow(d)), function(i) {
    as.numeric(mapply(function(k, s) {
      intrinsicDimension:::maxLikDimEstFromR(
        d[i, seq_len(k)], intrinsicDimension::dnoiseGaussH, s, ncol(x),
        unbiased = TRUE
      )
    }, grid$k, grid$s))
  }, numeric(nrow(grid)))
  matrix(est, ncol = nrow(grid), byrow = TRUE)
}

## Median of the cell estimates for each noise level and k, with the cells split across cores ##
dims_estimate <- function(x, at, ks, sigmas, cores = 1L, per_cell = cell_dims) {
  chunks <- lapply(parallel::splitIndices(length(at), min(cores, length(at))), function(i) at[i])
  res <- parallel::mclapply(chunks, function(ix) per_cell(x, ix, ks, sigmas), mc.cores = cores)
  ok <- vapply(res, is.matrix, logical(1))
  if (!all(ok)) stop("dims_estimate: ", sum(!ok), " of ", length(res), " chunks failed", call. = FALSE)
  data.frame(
    noise = rep(names(sigmas), each = length(ks)),
    k = rep(ks, length(sigmas)),
    estimate = apply(do.call(rbind, res), 2, stats::median, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
