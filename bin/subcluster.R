#!/usr/bin/env Rscript
## Re-cluster each cluster of one basis column on its own embedding. Adds one column and
## changes nothing else, so annotate can run off the subclusters instead ##
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON")),
  basis = arg("basis"),
  clusters = arg("clusters"),
  resolution = as.numeric(arg("resolution", "0.5")),
  knn_k = as.integer(arg("knn_k", "20")),
  n_iterations = as.integer(arg("n_iterations", "1000")),
  n_pcs = as.integer(arg("n_pcs", "30")),
  dims_use = arg("dims_use", "1:20"),
  n_hvgs = as.integer(arg("n_hvgs", "2000")),
  markers = arg("markers", ""),
  batch_column = arg("batch_column", ""),
  batch_correct = arg("batch_correct", "auto"),
  min_cells = as.integer(arg("min_cells", "200")),
  seed = as.integer(arg("seed", "42"))
)
if (!nzchar(opt$basis)) {
  stop("subcluster.R: --basis is required (e.g. leiden_clus_res0.3)")
}
if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
parse_dims <- function(s) {
  p <- strsplit(s, ":")[[1]]
  if (length(p) == 2) {
    as.integer(p[1]):as.integer(p[2])
  } else {
    as.integer(strsplit(s, ",")[[1]])
  }
}

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
if (!opt$basis %in% names(pd)) {
  stop(
    "subcluster.R: --basis '",
    opt$basis,
    "' not on the object. Available cluster columns: ",
    paste(grep("^leiden_clus_res", names(pd), value = TRUE), collapse = ", ")
  )
}

## The new column carries its basis, so two runs never collide ##
out_col <- paste0("sub_", sub("^leiden_clus_", "", opt$basis))
basis_v <- as.character(pd[[opt$basis]])
all_cl <- unique(basis_v)
all_cl <- all_cl[order(
  suppressWarnings(as.numeric(all_cl)),
  all_cl,
  na.last = TRUE
)]

want <- trimws(strsplit(opt$clusters, ",")[[1]])
want <- want[nzchar(want)]
if (!length(want)) {
  want <- all_cl
}
unknown <- setdiff(want, all_cl)
if (length(unknown)) {
  stop(
    "subcluster.R: --clusters not in ",
    opt$basis,
    ": ",
    paste(unknown, collapse = ", "),
    ". Available: ",
    paste(all_cl, collapse = ", ")
  )
}

canon <- if (nzchar(opt$markers) && file.exists(opt$markers)) {
  unique(read.csv(opt$markers)$gene)
} else {
  character()
}
bvars <- trimws(strsplit(opt$batch_column, ",")[[1]])
bvars <- bvars[nzchar(bvars)]

message(
  "subcluster.R: ",
  opt$basis,
  " -> ",
  out_col,
  "; ",
  length(want),
  " of ",
  length(all_cl),
  " cluster(s) at resolution ",
  opt$resolution
)

## Own HVGs, PCA and Harmony per cluster. The global embedding already made this split ##
sub_one <- function(cl) {
  ids <- pd$cell_ID[basis_v == cl]
  if (length(ids) < opt$min_cells) {
    message("  ", cl, ": ", length(ids), " cells < min_cells, left whole")
    return(data.frame(cell_ID = ids, sub = "1"))
  }
  g <- subsetGiotto(merged, cell_ids = ids)

  nm <- as.matrix(getExpression(
    g,
    feat_type = "rna",
    values = "normalized",
    output = "matrix"
  ))
  dec <- suppressWarnings(scran::modelGeneVar(nm))
  hvf <- suppressWarnings(scran::getTopHVGs(dec, n = opt$n_hvgs))
  feats <- union(hvf, intersect(canon, rownames(nm)))
  rm(nm)

  ncp <- min(opt$n_pcs, length(feats) - 1L, length(ids) - 1L)
  dims <- parse_dims(opt$dims_use)
  dims <- dims[dims <= ncp]
  g <- runPCA(
    g,
    feat_type = "rna",
    feats_to_use = feats,
    ncp = ncp,
    name = "pca_sub"
  )

  ## Same switch as integrate.R. A cluster can be single sample, so vars are re-tested ##
  sd <- pDataDT(g)
  have <- Filter(function(v) v %in% names(sd), bvars)
  use <- Filter(function(v) length(unique(sd[[v]])) >= 2L, have)
  drop <- setdiff(have, use)
  do_bc <- switch(
    tolower(opt$batch_correct),
    "true" = TRUE,
    "false" = FALSE,
    length(use) > 0L
  )
  if (do_bc && !length(use)) {
    do_bc <- FALSE
  }
  has_red <- function(obj, ty, nm) {
    !is.null(tryCatch(
      getDimReduction(
        obj,
        reduction = "cells",
        reduction_method = ty,
        name = nm,
        output = "matrix"
      ),
      error = function(e) NULL
    ))
  }

  red <- c(type = "pca", name = "pca_sub")
  if (do_bc) {
    h <- tryCatch(
      runGiottoHarmony(
        g,
        vars_use = use,
        dim_reduction_to_use = "pca",
        dim_reduction_name = "pca_sub",
        dimensions_to_use = dims,
        name = "harmony_sub",
        seed_number = opt$seed
      ),
      error = function(e) {
        message(
          "  ",
          cl,
          ": harmony failed (",
          conditionMessage(e),
          "), using PCA"
        )
        NULL
      }
    )
    hit <- if (is.null(h)) {
      NULL
    } else {
      Filter(
        function(nm) has_red(h, "harmony", nm),
        c("harmony_sub", "harmony")
      )
    }
    if (length(hit)) {
      g <- h
      red <- c(type = "harmony", name = hit[[1]])
    } else if (!is.null(h)) {
      message("  ", cl, ": harmony ran but no reduction found, using PCA")
    }
  }

  g <- createNearestNetwork(
    g,
    feat_type = "rna",
    type = "sNN",
    dim_reduction_to_use = red[["type"]],
    dim_reduction_name = red[["name"]],
    dimensions_to_use = dims,
    k = min(opt$knn_k, length(ids) - 1L),
    name = "NN.sub"
  )
  g <- doLeidenCluster(
    g,
    feat_type = "rna",
    nn_network_to_use = "sNN",
    network_name = "NN.sub",
    resolution = opt$resolution,
    n_iterations = opt$n_iterations,
    name = "sub",
    seed_number = opt$seed
  )
  g <- runUMAP(
    g,
    feat_type = "rna",
    dim_reduction_to_use = red[["type"]],
    dim_reduction_name = red[["name"]],
    dimensions_to_use = dims,
    name = "umap_sub",
    seed_number = opt$seed
  )
  u <- getDimReduction(
    g,
    reduction = "cells",
    reduction_method = "umap",
    name = "umap_sub",
    output = "matrix"
  )

  m <- pDataDT(g)
  message(
    "  ",
    cl,
    ": ",
    length(ids),
    " cells -> ",
    length(unique(m$sub)),
    " subclusters on ",
    red[["name"]],
    if (length(use)) paste0(" (", paste(use, collapse = ", "), ")") else "",
    if (length(drop)) {
      paste0("; single level here: ", paste(drop, collapse = ", "))
    } else {
      ""
    }
  )
  data.frame(
    cell_ID = m$cell_ID,
    sub = as.character(m$sub),
    U1 = u[match(m$cell_ID, rownames(u)), 1],
    U2 = u[match(m$cell_ID, rownames(u)), 2]
  )
}

## Build on an existing column so a second run adds its cluster instead of dropping the first ##
sub_v <- if (out_col %in% names(pd)) as.character(pd[[out_col]]) else basis_v
umaps <- list()
for (cl in want) {
  r <- sub_one(cl)
  if (is.null(r$U1)) {
    sub_v[match(r$cell_ID, pd$cell_ID)] <- paste0(cl, ".", r$sub)
    next
  }
  sub_v[match(r$cell_ID, pd$cell_ID)] <- paste0(cl, ".", r$sub)
  ## One embedding per parent: cells outside it are NA, so each gets its own view ##
  nm <- paste0(out_col, "_c", cl, "_UMAP")
  for (k in 1:2) {
    v <- rep(NA_real_, nrow(pd))
    v[match(r$cell_ID, pd$cell_ID)] <- r[[paste0("U", k)]]
    umaps[[paste0(nm, k)]] <- v
  }
  r$subcluster <- factor(
    paste0(cl, ".", r$sub),
    levels = paste0(cl, ".", sort(unique(as.numeric(r$sub))))
  )
  p <- ggplot(r, aes(U1, U2, colour = subcluster)) +
    geom_point(size = 0.3) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(title = paste0(out_col, ", cluster ", cl), colour = NULL) +
    theme_minimal(base_size = 11)
  ggsave(
    file.path(opt$outdir, paste0("umap_", out_col, "_c", cl, ".png")),
    p,
    width = 7,
    height = 6,
    dpi = 150,
    bg = "white"
  )
}

## "3.10" sorts after "3.2" ##
key <- do.call(
  rbind,
  lapply(strsplit(unique(sub_v), ".", fixed = TRUE), function(z) {
    suppressWarnings(as.numeric(c(z, NA_character_))[1:2])
  })
)
lev <- unique(sub_v)[order(key[, 1], key[, 2], na.last = FALSE)]

meta <- stats::setNames(
  data.frame(cell_ID = pd$cell_ID, factor(sub_v, levels = lev)),
  c("cell_ID", out_col)
)
for (nm in names(umaps)) {
  meta[[nm]] <- umaps[[nm]]
}

merged <- addCellMetadata(
  merged,
  by_column = TRUE,
  column_cell_ID = "cell_ID",
  new_metadata = meta
)

## Write ##
out <- stats::setNames(
  data.frame(pd$cell_ID, basis_v, sub_v),
  c("cell_ID", opt$basis, out_col)
)
for (nm in names(umaps)) {
  out[[nm]] <- umaps[[nm]]
}
write.csv(
  out,
  file.path(opt$outdir, paste0(out_col, "_cells.csv")),
  row.names = FALSE
)

sd <- pDataDT(merged)
summ <- as.data.frame(
  table(sd[[out_col]], dnn = "subcluster"),
  responseName = "n_cells"
)
summ$parent <- sub("\\..*$", "", as.character(summ$subcluster))
if ("sample_id" %in% names(sd)) {
  mix <- table(sd[[out_col]], sd$sample_id)
  top <- colnames(mix)[max.col(mix, "first")]
  summ$top_sample <- top[match(as.character(summ$subcluster), rownames(mix))]
  summ$top_sample_pct <- round(100 * apply(mix, 1, max) / rowSums(mix))[
    match(as.character(summ$subcluster), rownames(mix))
  ]
}
write.csv(
  summ,
  file.path(opt$outdir, paste0(out_col, "_summary.csv")),
  row.names = FALSE
)

if (
  "umap" %in%
    tryCatch(
      as.character(GiottoClass::list_dim_reductions(merged)$name),
      error = function(e) character()
    )
) {
  p <- plotUMAP(
    merged,
    feat_type = "rna",
    dim_reduction_name = "umap",
    cell_color = out_col,
    point_size = 0.4,
    show_plot = FALSE,
    return_plot = TRUE
  )
  ggsave(
    file.path(opt$outdir, paste0("umap_", out_col, ".png")),
    p,
    width = 9,
    height = 7,
    dpi = 150,
    bg = "white"
  )
}

saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
message("subcluster.R: ", out_col, " has ", length(lev), " levels")
writeLines(out_col, file.path(opt$outdir, "subcluster_column.txt"))
