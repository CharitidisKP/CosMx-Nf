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
  samplesheet = arg("samplesheet"),
  compare_by = arg("compare_by", "timepoint"),
  compare_baseline = arg("compare_baseline", "T0"),
  compare_within = arg("compare_within", "treatment"),
  pair_by = arg("pair_by", "subject_id"),
  celltype_column = arg("celltype_column", "cell_type"),
  exclude = arg("exclude"),
  sample_col = arg("sample_col", "sample_id"),
  fov_col = arg("fov_col", "fov"),
  min_cells = as.integer(arg("min_cells", "20")),
  fov_min_cells = as.integer(arg("fov_min_cells", "10")),
  fdr = as.numeric(arg("fdr", "0.05")),
  interaction = tolower(arg("interaction", "false")) %in% c("true", "1"),
  radius = as.numeric(arg("radius", "415")),
  overlap_flag = as.numeric(arg("overlap_flag", "2"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}
## edgeR 4 renamed calcNormFactors ##
norm_lib <- if (exists("normLibSizes", asNamespace("edgeR"))) {
  edgeR::normLibSizes
} else {
  edgeR::calcNormFactors
}
slug <- function(s) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", s))

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
raw <- getExpression(merged, values = "raw", output = "matrix")
locs <- getSpatialLocations(merged, output = "data.table")

## Only these tables are used below. The object would hold memory for the whole run ##
rm(merged)
invisible(gc())

## Sample facts come from the samplesheet, as in composition.R ##
ss <- fread(opt$samplesheet, colClasses = "character")
ss <- ss[sample_id %in% pd[[opt$sample_col]]]
within <- split_arg(opt$compare_within)
miss <- setdiff(c(opt$compare_by, within), names(ss))
if (length(miss)) {
  stop(
    "pseudobulk_de.R: samplesheet has no column ",
    paste(miss, collapse = ", ")
  )
}
ss[, grp := "all"]
if (length(within)) {
  ss[, grp := do.call(paste, c(.SD, sep = "/")), .SDcols = within]
}
ss[, level := get(opt$compare_by)]
ss[,
  subject := if (opt$pair_by %in% names(ss)) get(opt$pair_by) else NA_character_
]
ss[!nzchar(subject), subject := NA_character_]

## Typed cells only ##
has_fov <- opt$fov_col %in% names(pd)
if (!has_fov) {
  message("pseudobulk_de.R: no column '", opt$fov_col, "', FOV check skipped")
}
cells <- pd[,
  .(
    cell_ID,
    sample_id = get(opt$sample_col),
    cell_type = get(opt$celltype_column),
    fov = if (has_fov) as.character(get(opt$fov_col)) else NA_character_
  )
]
cells <- cells[
  !is.na(cell_type) &
    nzchar(cell_type) &
    !cell_type %in% split_arg(opt$exclude) &
    sample_id %in% ss$sample_id
]

## smiDE's overlap ratio per sample, cell type and gene: the mean in neighbouring cells of
## other types over the mean in the cell type itself. Every labelled cell is a neighbour ##
lab <- pd[
  !is.na(get(opt$celltype_column)) & get(opt$sample_col) %in% ss$sample_id,
  .(
    cell_ID,
    sample_id = get(opt$sample_col),
    cell_type = get(opt$celltype_column)
  )
]
lab[locs, `:=`(sdimx = i.sdimx, sdimy = i.sdimy), on = "cell_ID"]
ovl <- rbindlist(lapply(split(lab, by = "sample_id"), function(s) {
  o <- smiDE::overlap_ratio_metric(
    assay_matrix = raw[, s$cell_ID],
    metadata = as.data.frame(s),
    cluster_col = "cell_type",
    cellid_col = "cell_ID",
    sdimx_col = "sdimx",
    sdimy_col = "sdimy",
    radius = opt$radius
  )
  as.data.table(o)[, .(
    sample_id = s$sample_id[1],
    cell_type,
    gene = target,
    ratio
  )]
}))
raw <- raw[, cells$cell_ID]

## Summed raw counts per key, genes x keys ##
pseudobulk <- function(key) {
  f <- factor(key)
  m <- raw %*% Matrix::t(Matrix::fac2sparse(f))
  colnames(m) <- levels(f)
  m
}
cells[, key := paste(sample_id, cell_type, sep = "||")]
dsg <- cells[, .(n_cells = .N), by = .(key, sample_id, cell_type)]
dsg[, kept := n_cells >= opt$min_cells]
pb <- pseudobulk(cells$key)
if (has_fov) {
  cells[, fkey := paste(sample_id, fov, cell_type, sep = "||")]
  dfov <- cells[, .(n_cells = .N), by = .(fkey, sample_id, fov, cell_type)]
  dfov <- dfov[n_cells >= opt$fov_min_cells]
  pbf <- pseudobulk(cells$fkey)
}

## Comparisons: each later level against the baseline, inside each group ##
cmp <- ss[,
  if (opt$compare_baseline %in% level) {
    .(later = setdiff(unique(level), opt$compare_baseline))
  },
  by = grp
]
cmp[, name := paste0(grp, ": ", later, " vs ", opt$compare_baseline)]
side_of <- function(i) {
  s <- ss[grp == cmp$grp[i] & level %in% c(cmp$later[i], opt$compare_baseline)]
  setNames(
    fifelse(s$level == opt$compare_baseline, "reference", "case"),
    s$sample_id
  )
}

## logCPM on plain library sizes, for the FOV check ##
lcpm <- function(m) cpm(as.matrix(m), log = TRUE, prior.count = 2)

## One test per comparison and cell type. With residual df the paired QL test runs
## (~ subject + cond, or ~ subject + arm:post for the DiD), otherwise effect sizes only ##
run_test <- function(name, d, ct, did = FALSE) {
  status <- function(msg) {
    message("pseudobulk_de.R: ", name, ", ", ct, ": ", msg)
    list(summary = data.table(comparison = name, cell_type = ct, status = msg))
  }
  d <- merge(
    d,
    dsg[cell_type == ct & (kept)],
    by = "sample_id"
  )
  d[, subject := ss$subject[match(sample_id, ss$sample_id)]]
  if (uniqueN(d$cond) < 2) {
    return(status(sprintf("under %d cells on one side", opt$min_cells)))
  }
  both <- d[!is.na(subject), uniqueN(cond), by = subject][V1 == 2, subject]
  paired <- length(both) > 0L
  if (paired) {
    d <- d[subject %in% both]
  }
  if (did && d[, uniqueN(arm)] < 2) {
    return(status("one arm has no paired subject"))
  }

  y <- DGEList(as.matrix(pb[, d$key, drop = FALSE]))
  keep <- filterByExpr(y, group = d$cond)
  if (sum(keep) < 2L) {
    return(status("fewer than 2 expressed genes"))
  }
  y <- norm_lib(y[keep, , keep.lib.sizes = FALSE])
  d[, cond := factor(cond, levels = c("reference", "case"))]
  design <- if (did) {
    d[, post := as.numeric(cond == "case")]
    model.matrix(~ subject + arm:post, d)
  } else if (length(both) >= 2L) {
    model.matrix(~ subject + cond, d)
  } else {
    model.matrix(~cond, d)
  }
  df_res <- nrow(design) - qr(design)$rank
  lc <- cpm(y, log = TRUE, prior.count = 2)
  mean_lc <- function(keys) rowMeans(lc[, keys, drop = FALSE])

  ## Each subject's own change ##
  per <- rbindlist(lapply(both, function(s) {
    data.table(
      comparison = name,
      cell_type = ct,
      subject = s,
      gene = rownames(lc),
      log2_change = mean_lc(d[subject == s & cond == "case", key]) -
        mean_lc(d[subject == s & cond == "reference", key])
    )
  }))

  if (df_res >= 1L) {
    fit <- glmQLFit(y, design)
    test <- if (did) {
      arms <- grep(":post$", colnames(design), value = TRUE)
      con <- setNames(numeric(ncol(design)), colnames(design))
      con[arms] <- c(1, -1)
      glmQLFTest(fit, contrast = con)
    } else {
      glmQLFTest(fit, coef = "condcase")
    }
    tt <- as.data.table(
      topTags(test, n = Inf, sort.by = "none")$table,
      keep.rownames = "gene"
    )
  } else if (did) {
    return(status("df 0, no estimate"))
  } else {
    tt <- data.table(
      gene = rownames(lc),
      logFC = mean_lc(d[cond == "case", key]) -
        mean_lc(d[cond == "reference", key]),
      logCPM = aveLogCPM(y),
      F = NA_real_,
      PValue = NA_real_,
      FDR = NA_real_
    )
  }

  ## Every subject moves the way the pooled estimate does. Not defined for the DiD ##
  tt[, agree := NA]
  if (!did && length(both) >= 2L) {
    ag <- merge(per, tt[, .(gene, logFC)], by = "gene")[,
      .(agree = all(sign(log2_change) == sign(logFC) & log2_change != 0)),
      by = gene
    ]
    tt[ag, agree := i.agree, on = "gene"]
  }

  ## Share of case-side FOVs on the pooled side of their subject's baseline sample ##
  tt[, `:=`(fov_agree = NA_real_, n_fov = NA_integer_)]
  if (has_fov && !did && paired) {
    fv <- rbindlist(lapply(both, function(s) {
      case_fov <- dfov[
        cell_type == ct & sample_id %in% d[subject == s & cond == "case", sample_id]
      ]
      ref <- d[subject == s & cond == "reference", key]
      if (!nrow(case_fov)) {
        return(NULL)
      }
      g <- rownames(y)
      dlt <- lcpm(pbf[g, case_fov$fkey, drop = FALSE]) -
        rowMeans(lcpm(pb[g, ref, drop = FALSE]))
      data.table(gene = rep(g, ncol(dlt)), d = as.vector(dlt))
    }))
    if (nrow(fv)) {
      fa <- merge(fv, tt[, .(gene, logFC)], by = "gene")[,
        .(fov_agree = mean(sign(d) == sign(logFC)), n_fov = .N),
        by = gene
      ]
      tt[fa, `:=`(fov_agree = i.fov_agree, n_fov = i.n_fov), on = "gene"]
    }
  }

  ## Bleed could fake a change only where neighbours carry far more of the gene than the
  ## cell type does at baseline, so the highest baseline ratio decides the flag ##
  ov <- ovl[
    cell_type == ct & sample_id %in% d[cond == "reference", sample_id],
    .(
      overlap_ratio = if (all(is.na(ratio))) {
        NA_real_
      } else {
        max(ratio, na.rm = TRUE)
      }
    ),
    by = gene
  ]
  tt[, overlap_ratio := NA_real_]
  tt[ov, overlap_ratio := i.overlap_ratio, on = "gene"]
  tt[, bleed_risk := overlap_ratio > opt$overlap_flag]

  tt[, `:=`(comparison = name, cell_type = ct, paired = paired, df = df_res)]
  if (did) {
    tt[, note := "treatment confounded with batch; exploratory"]
  }
  hit <- tt[, !is.na(FDR) & FDR < opt$fdr]
  list(
    summary = data.table(
      comparison = name,
      cell_type = ct,
      status = if (df_res >= 1L) "tested" else "df 0, effect sizes only",
      paired = paired,
      n_subjects = length(both),
      n_samples = nrow(d),
      df = df_res,
      n_genes = nrow(tt),
      hits = sum(hit),
      hits_agree = sum(hit & tt$agree %in% TRUE),
      hits_flagged = sum(hit & tt$bleed_risk %in% TRUE)
    ),
    de = tt,
    per = per
  )
}

cts <- sort(unique(dsg$cell_type))
out <- list()
for (i in seq_len(nrow(cmp))) {
  side <- side_of(i)
  d <- data.table(sample_id = names(side), cond = unname(side))
  for (ct in cts) {
    out[[length(out) + 1L]] <- run_test(cmp$name[i], copy(d), ct)
  }
}

## Difference in differences between groups. Each subject's own baseline cancels its batch,
## but treatment still travels with batch, so the result is exploratory ##
if (opt$interaction && nrow(cmp) >= 2L) {
  for (pr in combn(seq_len(nrow(cmp)), 2L, simplify = FALSE)) {
    if (cmp$grp[pr[1]] == cmp$grp[pr[2]]) {
      next
    }
    name <- sprintf("DiD: (%s) - (%s)", cmp$name[pr[1]], cmp$name[pr[2]])
    d <- rbindlist(lapply(pr, function(i) {
      side <- side_of(i)
      data.table(
        sample_id = names(side),
        cond = unname(side),
        arm = paste0("arm", match(i, pr))
      )
    }))
    for (ct in cts) {
      out[[length(out) + 1L]] <- run_test(name, copy(d), ct, did = TRUE)
    }
  }
}
pick <- function(el) rbindlist(lapply(out, `[[`, el), fill = TRUE)
smry <- pick("summary")
de <- pick("de")
per <- pick("per")

## Write ##
long <- as.data.table(Matrix::summary(as(pb, "CsparseMatrix")))
long <- long[, .(
  key = colnames(pb)[j],
  gene = rownames(pb)[i],
  count = x
)]
long[dsg, `:=`(sample_id = i.sample_id, cell_type = i.cell_type), on = "key"]
fwrite(
  long[, .(sample_id, cell_type, gene, count)],
  file.path(opt$outdir, "pb_counts.csv.gz")
)
fwrite(
  merge(
    dsg[, .(sample_id, cell_type, n_cells, kept)],
    ss[, .(sample_id, group = grp, level, subject)],
    by = "sample_id"
  ),
  file.path(opt$outdir, "pb_design.csv")
)
fwrite(smry, file.path(opt$outdir, "pb_summary.csv"))
if (nrow(de)) {
  setcolorder(de, c("comparison", "cell_type", "gene"))
  setorder(de, comparison, cell_type, PValue, na.last = TRUE)
}
fwrite(de, file.path(opt$outdir, "pb_de.csv"))
fwrite(per, file.path(opt$outdir, "pb_per_subject_lfc.csv"))

tested <- if (nrow(de)) de[!is.na(FDR)] else de
for (nm in unique(tested$comparison)) {
  v <- tested[comparison == nm]
  v[,
    call := fcase(
      FDR < opt$fdr & bleed_risk %in% TRUE,
      "hit, possible bleed",
      FDR < opt$fdr & agree %in% TRUE,
      "hit, subjects agree",
      FDR < opt$fdr,
      "hit",
      default = "not significant"
    )
  ]
  p <- ggplot(v, aes(logFC, -log10(FDR), colour = call)) +
    geom_point(size = 0.4) +
    scale_colour_manual(
      values = c(
        `hit, subjects agree` = "firebrick",
        hit = "orange",
        `hit, possible bleed` = "steelblue",
        `not significant` = "grey70"
      )
    ) +
    facet_wrap(~cell_type) +
    labs(x = "log2 fold change", y = "-log10 FDR", colour = NULL, title = nm) +
    theme_minimal(base_size = 10)
  ggsave(
    file.path(opt$outdir, paste0("pb_volcano_", slug(nm), ".png")),
    p,
    width = 10,
    height = 3 + 2.5 * ceiling(uniqueN(v$cell_type) / 3),
    dpi = 150,
    bg = "white"
  )
}

## The first two subjects' changes against each other, per cell type ##
for (nm in unique(per$comparison)) {
  pp <- per[comparison == nm]
  subs <- head(sort(unique(pp$subject)), 2L)
  if (length(subs) < 2L || startsWith(nm, "DiD")) {
    next
  }
  w <- dcast(
    pp[subject %in% subs],
    cell_type + gene ~ subject,
    value.var = "log2_change"
  )
  setnames(w, subs, c("s1", "s2"))
  w[de[comparison == nm], hit := i.FDR < opt$fdr, on = .(cell_type, gene)]
  w[is.na(hit), hit := FALSE]
  p <- ggplot(w, aes(s1, s2, colour = hit)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(size = 0.4) +
    scale_colour_manual(
      values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
      labels = c(`FALSE` = "not significant", `TRUE` = "FDR hit")
    ) +
    facet_wrap(~cell_type) +
    labs(
      x = paste(subs[1], "log2 change"),
      y = paste(subs[2], "log2 change"),
      colour = NULL,
      title = nm
    ) +
    theme_minimal(base_size = 10)
  ggsave(
    file.path(opt$outdir, paste0("pb_lfc_agreement_", slug(nm), ".png")),
    p,
    width = 10,
    height = 3 + 2.5 * ceiling(uniqueN(w$cell_type) / 3),
    dpi = 150,
    bg = "white"
  )
}
message(sprintf(
  "pseudobulk_de.R: %d tests, %d hits at FDR < %.2f",
  if (nrow(smry)) nrow(smry[status == "tested"]) else 0L,
  if (nrow(smry)) sum(smry$hits, na.rm = TRUE) else 0L,
  opt$fdr
))
