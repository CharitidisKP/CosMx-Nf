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
  model_dir = arg("model_dir"),
  samplesheet = arg("samplesheet"),
  celltype_column = arg("celltype_column", "cell_type"),
  exclude = arg("exclude"),
  receivers = arg("receivers"),
  compare_by = arg("compare_by", "timepoint"),
  compare_baseline = arg("compare_baseline", "T0"),
  compare_within = arg("compare_within", "treatment"),
  pair_by = arg("pair_by", "subject_id"),
  sample_col = arg("sample_col", "sample_id"),
  expr_frac = as.numeric(arg("expr_frac", "0.10")),
  top_ligands = as.integer(arg("top_ligands", "20")),
  n_targets = as.integer(arg("n_targets", "200")),
  fdr = as.numeric(arg("fdr", "0.05"))
)

if (nzchar(opt$python)) {
  Sys.setenv(RETICULATE_PYTHON = opt$python)
}
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
suppressPackageStartupMessages({
  library(Giotto)
  library(data.table)
  library(Matrix)
  library(nichenetr)
  library(ggplot2)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}
model_file <- function(p) {
  f <- list.files(opt$model_dir, p, full.names = TRUE)
  if (!length(f)) {
    stop("cci_nichenet.R: no ", p, " in ", opt$model_dir)
  }
  f[1]
}

ligand_target_matrix <- readRDS(model_file("^ligand_target_matrix.*\\.rds$"))
lr_network <- as.data.table(readRDS(model_file("^lr_network.*\\.rds$")))
lr_network <- unique(lr_network[, .(from, to)])
weighted_networks <- readRDS(model_file("^weighted_networks.*\\.rds$"))

merged <- loadGiotto(opt$input)
pd <- pDataDT(merged)
norm <- getExpression(merged, values = "normalized", output = "matrix")
net <- getSpatialNetwork(
  merged,
  name = "Delaunay_network",
  output = "networkDT"
)
detection <- function(ids) Matrix::rowMeans(norm[, ids, drop = FALSE] > 0)

## Comparisons come from the samplesheet: each level of compare_by against the
## baseline, inside each compare_within group ##
ss <- fread(opt$samplesheet, colClasses = "character")
ss <- ss[sample_id %in% pd[[opt$sample_col]]]
within <- split_arg(opt$compare_within)
miss <- setdiff(c(opt$compare_by, within), names(ss))
if (length(miss)) {
  stop(
    "cci_nichenet.R: samplesheet has no column ",
    paste(miss, collapse = ", ")
  )
}
ss[, grp := "all"]
if (length(within)) {
  ss[, grp := do.call(paste, c(.SD, sep = "/")), .SDcols = within]
}
ss[,
  subject := if (opt$pair_by %in% names(ss)) get(opt$pair_by) else NA_character_
]
ss[!nzchar(subject), subject := NA_character_]
cmp <- ss[,
  if (opt$compare_baseline %in% get(opt$compare_by)) {
    .(level = setdiff(unique(get(opt$compare_by)), opt$compare_baseline))
  },
  by = grp
]
if (!nrow(cmp)) {
  stop(
    "cci_nichenet.R: no group has ",
    opt$compare_baseline,
    " and a later level"
  )
}

## The gene set is what changes in the receiver between the two sides. Blocked on
## subject when every subject has both, so pairs compare within the person ##
run_one <- function(g, l, r) {
  name <- paste0(g, ": ", l, " vs ", opt$compare_baseline)
  side <- ss[grp == g, setNames(get(opt$compare_by), sample_id)]
  rc <- pd[get(opt$celltype_column) == r & get(opt$sample_col) %in% names(side)]
  rc[, cond := side[get(opt$sample_col)]]
  rc <- rc[cond %in% c(l, opt$compare_baseline)]
  rc[, cond := fifelse(cond == l, "case", "reference")]
  rc[, subject := ss$subject[match(get(opt$sample_col), ss$sample_id)]]
  status <- function(msg) {
    message("cci_nichenet.R: ", name, ", ", r, ": ", msg)
    list(summary = data.table(comparison = name, receiver = r, status = msg))
  }
  if (uniqueN(rc$cond) < 2) {
    return(status("absent from one side"))
  }
  paired <- !anyNA(rc$subject) &&
    rc[, uniqueN(cond), by = subject][, all(V1 == 2)]
  mk <- scran::findMarkers(
    norm[, rc$cell_ID],
    groups = rc$cond,
    block = if (paired) rc$subject else NULL,
    direction = "up"
  )[["case"]]
  geneset <- intersect(
    rownames(mk)[which(mk$FDR < opt$fdr)],
    rownames(ligand_target_matrix)
  )
  expr_r <- names(which(detection(rc$cell_ID) >= opt$expr_frac))

  ## Senders are whatever touches the receiver on the case side, minus excluded labels ##
  rcase <- rc[cond == "case", cell_ID]
  nb <- setdiff(
    unique(c(net$to[net$from %in% rcase], net$from[net$to %in% rcase])),
    rcase
  )
  nb <- nb[
    !pd[match(nb, cell_ID), get(opt$celltype_column)] %in%
      split_arg(opt$exclude)
  ]
  expr_s <- names(which(detection(nb) >= opt$expr_frac))
  potential <- intersect(
    lr_network[from %in% expr_s & to %in% expr_r, unique(from)],
    colnames(ligand_target_matrix)
  )
  if (length(geneset) < 5 || !length(potential)) {
    return(status(sprintf(
      "%d target genes, %d ligands: too few",
      length(geneset),
      length(potential)
    )))
  }

  acts <- as.data.table(predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = intersect(
      expr_r,
      rownames(ligand_target_matrix)
    ),
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential
  ))
  setorder(acts, -aupr_corrected)
  best <- head(acts$test_ligand, opt$top_ligands)

  ## Which neighbouring cell types express each top ligand ##
  nb_ct <- pd[match(nb, cell_ID), get(opt$celltype_column)]
  senders <- rbindlist(
    lapply(split(nb, nb_ct), function(ids) {
      data.table(
        ligand = best,
        n_cells = length(ids),
        detection = detection(ids)[best]
      )
    }),
    idcol = "sender"
  )
  tag <- function(d) d[, `:=`(comparison = name, receiver = r)][]

  list(
    summary = data.table(
      comparison = name,
      receiver = r,
      status = "ok",
      paired = paired,
      n_case = sum(rc$cond == "case"),
      n_reference = sum(rc$cond == "reference"),
      n_neighbours = length(nb),
      n_geneset = length(geneset),
      n_ligands = length(potential)
    ),
    activities = tag(acts),
    targets = tag(rbindlist(lapply(best, function(x) {
      as.data.table(get_weighted_ligand_target_links(
        x,
        geneset,
        ligand_target_matrix,
        n = opt$n_targets
      ))
    }))),
    receptors = tag(as.data.table(get_weighted_ligand_receptor_links(
      best,
      expr_r,
      lr_network,
      weighted_networks$lr_sig
    ))),
    senders = tag(senders)
  )
}
out <- unlist(
  lapply(seq_len(nrow(cmp)), function(i) {
    lapply(split_arg(opt$receivers), function(r) {
      run_one(cmp$grp[i], cmp$level[i], r)
    })
  }),
  recursive = FALSE
)
pick <- function(el) rbindlist(lapply(out, `[[`, el), fill = TRUE)

act <- pick("activities")

## Write ##
if (nrow(act)) {
  top <- act[
    order(-aupr_corrected),
    head(.SD, 10),
    by = .(comparison, receiver)
  ]
  p <- ggplot(top, aes(receiver, test_ligand, fill = aupr_corrected)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "firebrick") +
    facet_wrap(~comparison, scales = "free_y") +
    labs(x = "Receiver", y = "Ligand", fill = "AUPR corrected") +
    theme_minimal(base_size = 10)
  ggsave(
    file.path(opt$outdir, "nichenet_ligand_activity.png"),
    p,
    width = max(6, 2 + 3 * uniqueN(top$comparison)),
    height = max(5, 1.5 + 0.2 * uniqueN(top$test_ligand)),
    dpi = 150,
    bg = "white"
  )
}
fwrite(pick("summary"), file.path(opt$outdir, "nichenet_summary.csv"))
fwrite(act, file.path(opt$outdir, "nichenet_ligand_activities.csv"))
fwrite(
  pick("targets"),
  file.path(opt$outdir, "nichenet_ligand_target_links.csv")
)
fwrite(
  pick("receptors"),
  file.path(opt$outdir, "nichenet_ligand_receptor_links.csv")
)
fwrite(
  pick("senders"),
  file.path(opt$outdir, "nichenet_ligand_senders.csv")
)
