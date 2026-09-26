#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON")),
  auto_lo = as.integer(arg("auto_lo", "2")),
  auto_hi = as.integer(arg("auto_hi", "20")),
  n_starts = as.integer(arg("n_starts", "5")),
  seed = as.integer(arg("seed", "42")),
  cohort = tolower(arg("cohort", "false")) %in% c("true", "1", "yes"),
  cohort_vars = arg("cohort_vars", ""),
  cohort_column = arg("cohort_column", ""),
  n_cohorts = arg("n_cohorts", ""),
  nb_k = as.integer(arg("nb_k", "50")), # neighbours per cell for the spatial component
  nb_pcs = as.integer(arg("nb_pcs", "10")), # 0 disables it
  reference = arg("reference", ""), # genes x cell_types profile matrix (.csv or .rds)
  semi_basis = arg("semi_basis", ""), # cluster column whose level count sets n_clusts
  semi_n = arg("semi_n", "") # explicit override for that count
)
if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(Matrix)
  library(InSituType)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

## n_clusts may be a range. InSituType picks the best count inside it ##
parse_n <- function(s) {
  p <- strsplit(s, ":")[[1]]
  if (length(p) == 2) {
    as.integer(p[1]):as.integer(p[2])
  } else {
    as.integer(strsplit(s, ",")[[1]])
  }
}

merged <- loadGiotto(opt$input)
counts <- t(getExpression(
  merged,
  feat_type = "rna",
  values = "raw",
  output = "matrix"
))
md <- pDataDT(merged, feat_type = "rna")
neg <- md$neg_mean[match(rownames(counts), md$cell_ID)]
neg[is.na(neg)] <- stats::median(neg, na.rm = TRUE)

## Cohorting splits cells on the CosMx immunofluorescence / morphology channels ##
cohort <- NULL
if (opt$cohort && nzchar(opt$cohort_column)) {
  if (opt$cohort_column %in% names(md)) {
    cohort <- as.character(md[[opt$cohort_column]])[match(
      rownames(counts),
      md$cell_ID
    )]
    cohort[is.na(cohort)] <- "unassigned"
    message(
      "InSituType: cohorting on metadata column '",
      opt$cohort_column,
      "' -> ",
      length(unique(cohort)),
      " cohorts"
    )
  } else {
    message(
      "InSituType: --cohort_column '",
      opt$cohort_column,
      "' not in metadata. Falling back to fastCohorting on --cohort_vars"
    )
  }
}
## Neighbourhood expression, so the cohort carries spatial context and not morphology alone ##
nb_pcs <- NULL
if (opt$cohort && is.null(cohort) && opt$nb_pcs > 0L) {
  loc <- getSpatialLocations(merged, output = "data.table")
  i <- match(rownames(counts), loc$cell_ID)
  nb <- InSituCor:::nearestNeighborGraph(
    x = loc$sdimx[i],
    y = loc$sdimy[i],
    N = opt$nb_k,
    subset = md$sample_id[match(rownames(counts), md$cell_ID)]
  )
  nb_expr <- InSituCor:::get_neighborhood_expression(
    counts = counts,
    neighbors = nb
  )
  nb_pcs <- irlba::prcomp_irlba(nb_expr, n = opt$nb_pcs)$x
  colnames(nb_pcs) <- paste0("nbPC", seq_len(ncol(nb_pcs)))
  message(
    "InSituType: ",
    ncol(nb_pcs),
    " neighbourhood PCs from a ",
    opt$nb_k,
    " nearest neighbour graph"
  )
}

if (opt$cohort && is.null(cohort)) {
  cv <- trimws(strsplit(opt$cohort_vars, ",")[[1]])
  have <- intersect(cv[nzchar(cv)], names(md))
  miss <- setdiff(cv[nzchar(cv)], have)
  if (length(miss)) {
    message(
      "InSituType: cohort_vars not in metadata, ignored: ",
      paste(miss, collapse = ", ")
    )
  }
  if (length(have) >= 2L) {
    cm <- as.matrix(md[match(rownames(counts), md$cell_ID), ..have])
    cm <- apply(cm, 2, function(x) {
      x <- suppressWarnings(as.numeric(x))
      x[!is.finite(x)] <- stats::median(x[is.finite(x)], na.rm = TRUE)
      x
    })
    keep <- apply(cm, 2, function(x) stats::sd(x) > 0) # a constant channel breaks the transform
    if (sum(keep) >= 2L) {
      cm <- cbind(cm[, keep, drop = FALSE], nb_pcs) # cbind with NULL is a no-op
      cohort <- InSituType::fastCohorting(
        cm,
        gaussian_transform = TRUE,
        n_cohorts = if (nzchar(opt$n_cohorts)) {
          as.integer(opt$n_cohorts)
        } else {
          NULL
        }
      )
      ctab <- table(cohort)
      message(
        "InSituType: cohorting on ",
        paste(colnames(cm), collapse = ", "),
        " -> ",
        length(ctab),
        " cohorts (",
        paste(ctab, collapse = "/"),
        ")"
      )
      write.csv(
        data.frame(cohort = names(ctab), n_cells = as.integer(ctab)),
        file.path(opt$outdir, "insitutype_cohorts.csv"),
        row.names = FALSE
      )
    } else {
      message(
        "InSituType: <2 non-constant cohort vars. Clustering without cohorts"
      )
    }
  } else {
    message("InSituType: <2 cohort vars available. Clustering without cohorts")
  }
}

set.seed(opt$seed)
ist <- InSituType::insitutype(
  x = counts,
  neg = neg,
  cohort = cohort,
  n_clusts = opt$auto_lo:opt$auto_hi,
  n_starts = opt$n_starts
)
n_pick <- length(unique(ist$clust))
stopifnot(n_pick > 1)

newmd <- data.frame(
  cell_ID = names(ist$clust),
  insitutype_unsup = as.character(ist$clust)
)
# Keep the cohort assignment: it tracks segmentation quality and IF intensity, and the
# clustering is not reproducible without it.
if (!is.null(cohort)) {
  newmd$insitutype_cohort <- as.character(cohort)[match(
    newmd$cell_ID,
    rownames(counts)
  )]
}
merged <- addCellMetadata(
  merged,
  by_column = TRUE,
  column_cell_ID = "cell_ID",
  new_metadata = newmd
)

tab <- table(ist$clust)
write.csv(
  data.frame(cluster = names(tab), n_cells = as.integer(tab)),
  file.path(opt$outdir, "insitutype_summary.csv"),
  row.names = FALSE
)

## --- semi-supervised pass --- Runs only when a reference profile matrix is given ##
if (nzchar(opt$reference)) {
  ref <- if (grepl("\\.rds$", opt$reference, ignore.case = TRUE)) {
    readRDS(opt$reference)
  } else {
    as.matrix(read.csv(opt$reference, row.names = 1, check.names = FALSE))
  }
  shared <- intersect(rownames(ref), colnames(counts))
  message(
    "insitutype: reference ",
    ncol(ref),
    " types, ",
    nrow(ref),
    " genes; ",
    length(shared),
    " shared with the panel"
  )
  if (length(shared) < 50L) {
    message("insitutype: <50 shared gene. Skipping semi-supervised pass.")
  } else {
    n_semi <- if (nzchar(opt$semi_n)) {
      parse_n(opt$semi_n)
    } else if (
      nzchar(opt$semi_basis) &&
        opt$semi_basis %in% names(md)
    ) {
      length(unique(md[[opt$semi_basis]]))
    } else {
      NA_integer_
    }
    if (anyNA(n_semi)) {
      stop(
        "--reference given but n_clusts unresolved: pass --semi_n or a --semi_basis column present in the object"
      )
    }
    message(
      "insitutype: semi-supervised with n_clusts = ",
      paste0(min(n_semi), "-", max(n_semi)),
      if (nzchar(opt$semi_n)) {
        " (explicit)"
      } else {
        paste0(" (levels of ", opt$semi_basis, ")")
      }
    )
    fit <- function(...) {
      InSituType::insitutype(
        x = counts,
        neg = neg,
        cohort = cohort,
        reference_profiles = ref,
        n_clusts = n_semi,
        align_genes = TRUE,
        n_starts = opt$n_starts,
        ...
      )
    }
    set.seed(opt$seed)

    # Anchor selection can fail when confidence scores are too uniform.
    # update_reference_profiles = FALSE skips anchoring and uses the reference as given.
    # refinement = FALSE does NOT help here despite what the error message suggests.
    semi <- tryCatch(fit(), error = function(e) {
      if (!grepl("anchor", conditionMessage(e), ignore.case = TRUE)) {
        stop(e)
      }
      message(
        "insitutype: anchor selection failed. Retrying with update_reference_profiles = FALSE"
      )
      set.seed(opt$seed)
      fit(update_reference_profiles = FALSE)
    })
    merged <- addCellMetadata(
      merged,
      by_column = TRUE,
      column_cell_ID = "cell_ID",
      new_metadata = data.frame(
        cell_ID = names(semi$clust),
        insitutype_semisup = as.character(semi$clust)
      )
    )
    st <- table(semi$clust)
    write.csv(
      data.frame(
        cluster = names(st),
        n_cells = as.integer(st),
        novel = !names(st) %in% colnames(ref)
      ),
      file.path(opt$outdir, "insitutype_semisup_summary.csv"),
      row.names = FALSE
    )
  }
}

saveGiotto(merged, dir = opt$outdir, foldername = "giotto", overwrite = TRUE)
