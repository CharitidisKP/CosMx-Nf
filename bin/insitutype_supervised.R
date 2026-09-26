#!/usr/bin/env Rscript
## Supervised InSituType against one reference profile matrix. insitutype() with n_clusts = 0.
## Keeps anchor selection and platform rescaling. insitutypeML would skip both. ##
## Writes per-cell CSVs. attach_supervised.R merges them onto the Giotto object.
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  input = arg("input"),
  outdir = arg("outdir", "."),
  python = arg("python", Sys.getenv("RETICULATE_PYTHON")),
  reference = arg("reference"),
  name = arg("name", "ref"), # column suffix: insitutype_sup_<name>
  refine = tolower(arg("refine", "false")) %in% c("true", "1", "yes"),
  conf = as.numeric(arg("conf_threshold", "0.8")),
  min_genes = as.integer(arg("min_gene_overlap", "100")),
  cohort = tolower(arg("cohort", "false")) %in% c("true", "1", "yes"),
  cohort_vars = arg("cohort_vars", ""),
  n_cohorts = arg("n_cohorts", ""),
  seed = as.integer(arg("seed", "42"))
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
if (!nzchar(opt$reference)) {
  stop("--reference is required")
}

merged <- loadGiotto(opt$input)
counts <- t(getExpression(merged, values = "raw", output = "matrix")) # cells x genes
md <- pDataDT(merged)
neg <- md$neg_mean[match(rownames(counts), md$cell_ID)]
neg[is.na(neg)] <- stats::median(neg, na.rm = TRUE)

## Same IF/morphology cohorting as the unsupervised pass. Without it only the RNA panel decides ##
cohort <- NULL
if (opt$cohort) {
  cv <- trimws(strsplit(opt$cohort_vars, ",")[[1]])
  have <- intersect(cv[nzchar(cv)], names(md))
  if (length(have) >= 2L) {
    cm <- as.matrix(md[match(rownames(counts), md$cell_ID), ..have])
    cm <- apply(cm, 2, function(x) {
      x <- suppressWarnings(as.numeric(x))
      x[!is.finite(x)] <- stats::median(x[is.finite(x)], na.rm = TRUE)
      x
    })
    keep <- apply(cm, 2, function(x) stats::sd(x) > 0)
    if (sum(keep) >= 2L) {
      cohort <- InSituType::fastCohorting(
        cm[, keep, drop = FALSE],
        gaussian_transform = TRUE,
        n_cohorts = if (nzchar(opt$n_cohorts)) {
          as.integer(opt$n_cohorts)
        } else {
          NULL
        }
      )
      message(
        "supervised[",
        opt$name,
        "]: ",
        length(unique(cohort)),
        " cohorts from ",
        paste(colnames(cm)[keep], collapse = ", ")
      )
    } else {
      message(
        "supervised[",
        opt$name,
        "]: <2 non-constant cohort vars; no cohorting"
      )
    }
  } else {
    message(
      "supervised[",
      opt$name,
      "]: <2 cohort vars available; no cohorting"
    )
  }
}

ref <- if (grepl("\\.rds$", opt$reference, ignore.case = TRUE)) {
  as.matrix(readRDS(opt$reference))
} else {
  as.matrix(read.csv(opt$reference, row.names = 1, check.names = FALSE))
}
shared <- intersect(rownames(ref), colnames(counts))
message(
  "supervised[",
  opt$name,
  "]: reference ",
  ncol(ref),
  " types x ",
  nrow(ref),
  " genes; ",
  length(shared),
  " shared with the panel"
)
if (length(shared) < opt$min_genes) {
  stop(
    "supervised[",
    opt$name,
    "]: only ",
    length(shared),
    " genes shared with the panel (need >= ",
    opt$min_genes,
    "). Wrong profile for this panel, or gene symbols use a different convention."
  )
}

# Label normalisation: HCA names carry dots ("T.cell"), CosMx profiles use spaces. Downstream
# joins key on the label string, so settle on one convention here rather than in every consumer.
norm_lab <- function(x) gsub("[._[:space:]]+", "_", trimws(as.character(x)))

set.seed(opt$seed)
fit <- function(update = TRUE) {
  InSituType::insitutype(
    x = counts,
    neg = neg,
    cohort = cohort,
    reference_profiles = ref,
    n_clusts = 0,
    align_genes = TRUE,
    update_reference_profiles = update,
    rescale = TRUE,
    refit = FALSE
  )
}
sup <- tryCatch(fit(TRUE), error = function(e) {
  if (!grepl("anchor", conditionMessage(e), ignore.case = TRUE)) {
    stop(e)
  }
  message(
    "supervised[",
    opt$name,
    "]: anchor selection failed. Falling back to the reference uncorrected"
  )
  set.seed(opt$seed)
  fit(FALSE)
})

sup$clust <- norm_lab(sup$clust)
if (!is.null(sup$logliks)) {
  colnames(sup$logliks) <- norm_lab(colnames(sup$logliks))
}
if (!is.null(sup$profiles)) {
  colnames(sup$profiles) <- norm_lab(colnames(sup$profiles))
}
message(
  "supervised[",
  opt$name,
  "]: ",
  length(unique(sup$clust)),
  " cell types assigned"
)

## Anchors drove the platform correction, so record which cells and which types ##
out_anchor <- rep(NA_character_, nrow(counts))
if (!is.null(sup$anchors)) {
  an <- sup$anchors
  an[] <- norm_lab(an)
  out_anchor <- unname(an[match(rownames(counts), names(an))])
  at <- table(out_anchor[!is.na(out_anchor)])
  message(
    "supervised[",
    opt$name,
    "]: ",
    sum(at),
    " anchor cells across ",
    length(at),
    " types"
  )
}

## insitutype() returns $prob as a cells x types matrix. insitutypeML returns a per-cell scalar.
## Both branches kept so the function can be swapped back without touching anything below.
score_of <- function(res) {
  p <- res$prob
  if (is.matrix(p)) {
    vapply(
      seq_along(res$clust),
      function(i) {
        ct <- res$clust[i]
        if (is.na(ct) || !ct %in% colnames(p)) NA_real_ else p[i, ct]
      },
      numeric(1)
    )
  } else {
    as.numeric(p)
  }
}

# Row-wise softmax of the loglik matrix -> per-cell posterior over all types.
posterior <- function(res) {
  ll <- res$logliks
  if (is.null(ll) || !is.matrix(ll)) {
    return(NULL)
  }
  e <- exp(ll - apply(ll, 1, max, na.rm = TRUE))
  e / rowSums(e, na.rm = TRUE)
}

# Per-cell decision statistics. score / margin / entropy are the probability-scale confidence
# measures (same softmax as InSituType:::logliks2probs). loglik_margin is the best-vs-runner-up
# gap in nats, which stays informative when two correlated profiles compress the probability margin.
cell_stats <- function(res, post) {
  n <- length(res$clust)
  ll <- res$logliks
  na_df <- data.frame(
    second_type = rep(NA_character_, n),
    second_score = NA_real_,
    margin = NA_real_,
    entropy = NA_real_,
    loglik_top = NA_real_,
    loglik_margin = NA_real_,
    stringsAsFactors = FALSE
  )
  if (is.null(post) || is.null(ll) || !is.matrix(ll)) {
    return(na_df)
  }
  ord <- t(apply(ll, 1, function(r) order(r, decreasing = TRUE)[1:2])) # rank on logliks
  cn <- colnames(post)
  p1 <- post[cbind(seq_len(n), ord[, 1])]
  p2 <- post[cbind(seq_len(n), ord[, 2])]
  l1 <- ll[cbind(seq_len(n), ord[, 1])]
  l2 <- ll[cbind(seq_len(n), ord[, 2])]
  ent <- -rowSums(post * log(pmax(post, .Machine$double.xmin)), na.rm = TRUE)
  data.frame(
    second_type = cn[ord[, 2]],
    second_score = p2,
    margin = p1 - p2,
    entropy = ent,
    loglik_top = l1,
    loglik_margin = l1 - l2,
    stringsAsFactors = FALSE
  )
}

post <- posterior(sup)
out <- data.frame(
  cell_ID = rownames(counts),
  celltype = as.character(sup$clust),
  anchor = out_anchor,
  score = score_of(sup),
  stringsAsFactors = FALSE
)
out <- cbind(out, cell_stats(sup, post))

message(sprintf(
  "supervised[%s]: score quantiles (0/25/50/75/100%%): %s",
  opt$name,
  paste(signif(stats::quantile(out$score, na.rm = TRUE), 4), collapse = " / ")
))
message(sprintf(
  "supervised[%s]: loglik_margin quantiles (0/25/50/75/100%%): %s",
  opt$name,
  paste(
    signif(stats::quantile(out$loglik_margin, na.rm = TRUE), 4),
    collapse = " / "
  )
))
# Flags the degenerate case where every posterior is one-hot, so no conf threshold below 1
# could select a type.
if (isTRUE(mean(out$score > 1 - 1e-9, na.rm = TRUE) > 0.999)) {
  message(
    "supervised[",
    opt$name,
    "]: every posterior is 1.000 -- score/margin/entropy carry no",
    " information for this run; rank ambiguity on loglik_margin instead."
  )
}

# Per-type mean posterior (the quantity refineClusters thresholds on), plus the score
# distribution per type so a weak call is visible without reloading the logliks.
conf <- if (is.null(post)) {
  NULL
} else {
  vapply(
    colnames(post),
    function(ct) {
      i <- which(sup$clust == ct)
      if (!length(i)) NA_real_ else mean(post[i, ct], na.rm = TRUE)
    },
    numeric(1)
  )
}
if (!is.null(conf)) {
  tb <- table(sup$clust)
  cf <- data.frame(
    cell_type = names(conf),
    mean_conf = round(unname(conf), 4),
    n_cells = as.integer(tb[names(conf)]),
    stringsAsFactors = FALSE
  )
  agg <- function(f) {
    vapply(
      cf$cell_type,
      function(ct) {
        v <- out$score[out$celltype == ct]
        if (!length(v)) NA_real_ else f(v)
      },
      numeric(1)
    )
  }
  cf$mean_score <- round(agg(function(v) mean(v, na.rm = TRUE)), 4)
  cf$median_score <- round(agg(function(v) stats::median(v, na.rm = TRUE)), 4)
  cf$frac_score_0.8 <- round(agg(function(v) mean(v >= 0.8, na.rm = TRUE)), 4)
  cf$mean_margin <- round(
    vapply(
      cf$cell_type,
      function(ct) {
        v <- out$margin[out$celltype == ct]
        if (!length(v)) NA_real_ else mean(v, na.rm = TRUE)
      },
      numeric(1)
    ),
    4
  )
  # loglik-scale equivalents, for picking types to drop when the probability margin between two
  # correlated reference profiles is too compressed to separate them.
  cf$median_loglik_margin <- round(
    vapply(
      cf$cell_type,
      function(ct) {
        v <- out$loglik_margin[out$celltype == ct]
        if (!length(v)) NA_real_ else stats::median(v, na.rm = TRUE)
      },
      numeric(1)
    ),
    3
  )
  cf$median_loglik_top <- round(
    vapply(
      cf$cell_type,
      function(ct) {
        v <- out$loglik_top[out$celltype == ct]
        if (!length(v)) NA_real_ else stats::median(v, na.rm = TRUE)
      },
      numeric(1)
    ),
    1
  )
  write.csv(
    cf[order(cf$mean_conf), ],
    file.path(
      opt$outdir,
      paste0("insitutype_sup_", opt$name, "_confidence.csv")
    ),
    row.names = FALSE
  )
}

# The updated reference profiles. Every per-cell statistic is in the _cells.csv
if (!is.null(sup$profiles)) {
  write.csv(
    as.matrix(sup$profiles),
    file.path(opt$outdir, paste0("insitutype_sup_", opt$name, "_profiles.csv"))
  )
}

# --- optional refinement --------------------------------------------------------------------
# refineClusters() edits the logliks matrix: types whose mean posterior falls below conf are
# dropped and their cells reassigned to the next-best type. It needs >= 2 types to survive.
if (opt$refine) {
  if (is.null(conf)) {
    message("supervised[", opt$name, "]: no logliks matrix; refinement skipped")
  } else {
    to_delete <- names(conf)[!is.na(conf) & conf < opt$conf]
    remaining <- setdiff(colnames(sup$logliks), to_delete)
    if (!length(to_delete)) {
      # Distinguishes "every type passed the threshold" from the degenerate case where mean_conf
      # is exactly 1.0 everywhere and no threshold below 1 could fire.
      message(
        "supervised[",
        opt$name,
        "]: no type below conf ",
        opt$conf,
        "; nothing to refine",
        if (all(conf > 1 - 1e-9, na.rm = TRUE)) {
          paste0(
            " -- NOTE: mean_conf is exactly 1.000 for all ",
            length(conf),
            " types, so no threshold below 1 could select one. Drop types by",
            " mean_margin or median_loglik_margin in the confidence CSV instead."
          )
        } else {
          ""
        }
      )
    } else if (length(remaining) < 2L) {
      message(
        "supervised[",
        opt$name,
        "]: conf ",
        opt$conf,
        " would leave ",
        length(remaining),
        " type(s); refinement skipped (needs >= 2)"
      )
    } else {
      message(
        "supervised[",
        opt$name,
        "]: refining, dropping ",
        length(to_delete),
        " low-confidence type(s): ",
        paste(to_delete, collapse = ", ")
      )
      ref_res <- tryCatch(
        suppressMessages(InSituType::refineClusters(
          logliks = sup$logliks,
          to_delete = to_delete,
          counts = counts,
          neg = neg,
          cohort = cohort
        )),
        error = function(e) {
          message(
            "supervised[",
            opt$name,
            "]: refineClusters failed: ",
            conditionMessage(e)
          )
          NULL
        }
      )
      if (!is.null(ref_res)) {
        ref_res$clust <- norm_lab(ref_res$clust)
        if (!is.null(ref_res$prob) && is.matrix(ref_res$prob)) {
          colnames(ref_res$prob) <- norm_lab(colnames(ref_res$prob))
        }
        if (!is.null(ref_res$logliks)) {
          colnames(ref_res$logliks) <- norm_lab(colnames(ref_res$logliks))
        }
        rpost <- posterior(ref_res)
        rs <- cell_stats(ref_res, rpost)
        out$celltype_refined <- as.character(ref_res$clust)
        out$score_refined <- score_of(ref_res)
        out$second_type_refined <- rs$second_type
        out$second_score_refined <- rs$second_score
        out$margin_refined <- rs$margin
        out$entropy_refined <- rs$entropy
        out$loglik_top_refined <- rs$loglik_top
        out$loglik_margin_refined <- rs$loglik_margin
        message(
          "supervised[",
          opt$name,
          "]: ",
          length(unique(na.omit(ref_res$clust))),
          " cell types after refinement"
        )
      }
    }
  }
}

write.csv(
  out,
  file.path(opt$outdir, paste0("insitutype_sup_", opt$name, "_cells.csv")),
  row.names = FALSE
)
st <- table(out$celltype)
summ <- data.frame(cell_type = names(st), n_cells = as.integer(st))
if (!is.null(out$celltype_refined)) {
  rt <- table(out$celltype_refined)
  summ$n_cells_refined <- as.integer(rt[match(summ$cell_type, names(rt))])
}
write.csv(
  summ,
  file.path(opt$outdir, paste0("insitutype_sup_", opt$name, "_summary.csv")),
  row.names = FALSE
)
