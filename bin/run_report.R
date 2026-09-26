#!/usr/bin/env Rscript
## Writes run_report.md from every table published so far. The stage that calls it hands
## its own tables over in --current, because publishing into --results is asynchronous ##
a <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = "") {
  i <- which(a == paste0("--", k))
  if (length(i)) a[[i + 1L]] else d
}
opt <- list(
  results = arg("results"),
  current = arg("current", "current"),
  samplesheet = arg("samplesheet"),
  run_id = arg("run_id", "default"),
  stage = arg("stage", "discover"),
  commit = arg("commit", "unknown"),
  basis = arg("basis"),
  untyped = arg("untyped"),
  chosen_resolutions = arg("chosen_resolutions"),
  updated = arg("updated", format(Sys.time(), "%Y-%m-%d %H:%M")),
  outdir = arg("outdir", ".")
)

suppressPackageStartupMessages({
  library(data.table)
})
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
split_arg <- function(s) {
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}

## Tables by file name: this stage's own first (one numbered folder per file), then the
## published tree ##
find_tables <- function(pattern) {
  cur <- Sys.glob(paste0(opt$current, "*"))
  hits <- c(
    if (length(cur)) list.files(cur, pattern, full.names = TRUE),
    list.files(
      file.path(opt$results, "tables"),
      pattern,
      recursive = TRUE,
      full.names = TRUE
    )
  )
  hits[!duplicated(basename(hits))]
}
read_tables <- function(pattern, ...) {
  f <- find_tables(pattern)
  if (!length(f)) {
    return(NULL)
  }
  d <- rbindlist(lapply(f, fread, ...), fill = TRUE)
  if (nrow(d)) d else NULL
}
## A section whose tables cannot be read says so, and the others still appear ##
guard <- function(expr) {
  tryCatch(
    expr,
    error = function(e) paste0("_Not summarised: ", conditionMessage(e), "_")
  )
}
has <- function(d, cols) !is.null(d) && all(cols %in% names(d))
num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

## Figures belong to the step that drew them. The calling stage's own may still be
## copying, so they are linked without checking ##
stage_steps <- list(
  discover = c(
    "01_qc",
    "02_integration",
    "03_clustering",
    "04_reference_typing"
  ),
  subcluster = "05_subclustering",
  annotate = "06_annotation",
  spatial = "07_spatial",
  cci = "08_cci",
  de = "09_de",
  diagnose = c("03_clustering", "04_reference_typing", "05_subclustering"),
  supervised = "04_reference_typing"
)
fig <- function(path, alt = basename(path)) {
  step <- strsplit(path, "/", fixed = TRUE)[[1]][2]
  mine <- step %in% stage_steps[[opt$stage]]
  if (mine || file.exists(file.path(opt$results, path))) {
    c(paste0("![", alt, "](", path, ")"), "")
  }
}
fig_list <- function(paths) {
  paths <- paths[
    vapply(
      paths,
      function(p) {
        step <- strsplit(p, "/", fixed = TRUE)[[1]][2]
        step %in%
          stage_steps[[opt$stage]] ||
          file.exists(file.path(opt$results, p))
      },
      logical(1L)
    )
  ]
  if (length(paths)) paste0("- [", basename(paths), "](", paths, ")")
}

md_table <- function(d, digits = 2L) {
  d <- as.data.frame(d)
  for (j in seq_along(d)) {
    if (is.numeric(d[[j]])) {
      d[[j]] <- format(round(d[[j]], digits), big.mark = ",", trim = TRUE)
    }
    d[[j]] <- ifelse(is.na(d[[j]]) | d[[j]] == "NA", "", as.character(d[[j]]))
  }
  c(
    paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")),
    ""
  )
}
section <- function(title, body) {
  body <- unlist(body)
  if (length(body)) c(paste("##", title), "", body, "") else NULL
}

## Header and stages run ##
basis <- opt$basis
pj <- file.path(opt$results, "pipeline_info", "annotate_params.json")
if (
  !nzchar(basis) &&
    file.exists(pj) &&
    requireNamespace("jsonlite", quietly = TRUE)
) {
  b <- tryCatch(jsonlite::fromJSON(pj)$cluster_basis, error = function(e) NULL)
  if (is.character(b) && length(b) == 1L) basis <- b
}
header <- c(
  paste("# Run report:", opt$run_id),
  "",
  sprintf(
    "Updated %s, after the **%s** stage. Pipeline commit `%s`.",
    opt$updated,
    opt$stage,
    opt$commit
  ),
  ""
)
log_f <- file.path(opt$results, "pipeline_info", "stage_log.tsv")
stages <- if (file.exists(log_f)) {
  lg <- fread(log_f, sep = "\t", quote = "")
  lg[,
    .(launches = .N, last_started = max(started), commit = commit[.N]),
    by = stage
  ]
}

## 1. Samples and QC ##
ss <- fread(opt$samplesheet, colClasses = "character")
included <- if ("include" %in% names(ss)) {
  ss$sample_id[toupper(trimws(ss$include)) != "FALSE"]
} else {
  ss$sample_id
}
s1 <- guard({
  qc <- read_tables("_qc_summary\\.csv$")
  ld <- read_tables("_load_summary\\.csv$")
  if (has(qc, c("sample_id", "cells_before", "cells_after"))) {
    t1 <- merge(
      ss[,
        intersect(
          c("sample_id", "treatment", "timepoint", "subject_id"),
          names(ss)
        ),
        with = FALSE
      ],
      qc,
      by = "sample_id"
    )
    t1[, pct_kept := 100 * cells_after / cells_before]
    keep_cols <- intersect(
      c(
        "sample_id",
        "treatment",
        "timepoint",
        "subject_id",
        "cells_before",
        "cells_after",
        "pct_kept",
        "min_count",
        "median_counts_kept",
        "median_genes_kept"
      ),
      names(t1)
    )
    c(
      md_table(t1[, ..keep_cols], 1L),
      sprintf(
        "%s cells kept of %s loaded; %s genes after QC.",
        format(sum(t1$cells_after), big.mark = ","),
        format(sum(t1$cells_before), big.mark = ","),
        if (has(qc, "genes_after")) num(max(qc$genes_after)) else "?"
      ),
      "",
      fig_list(paste0("figures/01_qc/qc_", t1$sample_id, ".png"))
    )
  } else if (!is.null(ld)) {
    md_table(ld)
  }
})

## 2. Integration ##
s2 <- guard({
  ig <- read_tables("^integrate_summary\\.csv$")
  if (!is.null(ig)) {
    c(
      sprintf(
        "- Scale factor %s (median library size %s); %s variable genes, %s PCA features.",
        ig$scalefactor[1],
        ig$lib_median[1],
        num(ig$n_hvg[1]),
        num(ig$n_pca_feats[1])
      ),
      sprintf(
        "- Batch correction: %s on %s (levels %s); clustering on %s, dimensions %s.",
        if (isTRUE(as.logical(ig$batch_corrected[1]))) "Harmony" else "none",
        ig$batch_used[1],
        ig$batch_levels[1],
        ig$reduction[1],
        ig$dims[1]
      ),
      "",
      fig("figures/02_integration/scree.png", "PCA scree"),
      fig(
        "figures/03_clustering/umap_pca_by_batch.png",
        "Before Harmony, by batch"
      ),
      fig("figures/03_clustering/umap_by_batch.png", "After Harmony, by batch")
    )
  }
})

## 3. Clustering ##
chosen_res <- split_arg(opt$chosen_resolutions)
s3 <- guard({
  cs <- read_tables(
    "^cluster_summary\\.csv$",
    colClasses = c(resolution = "character")
  )
  if (has(cs, c("resolution", "n_clusters"))) {
    cs[,
      max_sample_share := vapply(
        resolution,
        function(r) {
          mx <- read_tables(paste0("^sample_mix_res", r, "\\.csv$"))
          if (is.null(mx)) {
            return(NA_real_)
          }
          max(as.matrix(mx[, -1]), na.rm = TRUE)
        },
        numeric(1L)
      )
    ]
    cs[,
      mostly_one_sample := vapply(
        resolution,
        function(r) {
          mx <- read_tables(paste0("^sample_mix_res", r, "\\.csv$"))
          if (is.null(mx)) {
            return(NA_integer_)
          }
          sum(apply(as.matrix(mx[, -1]), 1, max) > 0.8)
        },
        integer(1L)
      )
    ]
    cs[, chosen := ifelse(resolution %in% chosen_res, "yes", "")]
    c(
      md_table(cs[, .(
        resolution,
        n_clusters,
        max_sample_share,
        mostly_one_sample,
        chosen
      )]),
      paste(
        "Max sample share: the largest share of one cluster that comes from a single sample.",
        "Mostly one sample: clusters with more than 80% of their cells from one sample.",
        "A small fragment can reach 1.00 on its own."
      ),
      "",
      unlist(lapply(chosen_res, function(r) {
        fig(
          paste0("figures/03_clustering/umap_leiden_res", r, ".png"),
          paste("Leiden", r)
        )
      }))
    )
  }
})

## 4. Reference typing ##
s4 <- guard({
  ht <- read_tables("^hieratype_summary\\.csv$")
  iu <- read_tables("^insitutype_summary\\.csv$")
  conf_files <- find_tables("^insitutype_sup_.+_confidence\\.csv$")
  c(
    if (has(ht, c("ht_call", "pct"))) {
      sprintf(
        "- HieraType called %.1f%% of cells (the rest `unknown`).",
        100 - sum(ht[ht_call == "unknown", pct])
      )
    },
    if (has(iu, c("cluster", "n_cells"))) {
      sprintf("- InSituType unsupervised: %d clusters.", nrow(iu))
    },
    if (length(conf_files)) {
      sup <- rbindlist(lapply(conf_files, function(f) {
        d <- fread(f)
        nm <- sub("^insitutype_sup_(.+)_confidence\\.csv$", "\\1", basename(f))
        if (!all(c("n_cells", "frac_score_0.8", "cell_type") %in% names(d))) {
          return(NULL)
        }
        d <- d[order(-n_cells)]
        data.table(
          profile = nm,
          types_used = sum(d$n_cells > 0),
          pct_score_0.8 = 100 *
            sum(d$n_cells * d$frac_score_0.8) /
            sum(d$n_cells),
          largest_types = paste(head(d$cell_type, 4), collapse = ", ")
        )
      }))
      if (nrow(sup)) c("", md_table(sup, 1L))
    },
    unlist(lapply(
      c(
        "ht_call",
        "insitutype_unsup",
        sub("_confidence\\.csv$", "", basename(conf_files))
      ),
      function(cl) {
        fig(paste0("figures/04_reference_typing/umap_", cl, ".png"), cl)
      }
    ))
  )
})

## 5. Subclustering ##
s5 <- guard({
  sub_files <- find_tables("^sub_.+_summary\\.csv$")
  if (length(sub_files)) {
    unlist(lapply(sub_files, function(f) {
      d <- fread(
        f,
        colClasses = c(subcluster = "character", parent = "character")
      )
      sc <- sub("_summary\\.csv$", "", basename(f))
      per <- d[, .(subclusters = .N, cells = sum(n_cells)), by = parent]
      multi <- per[subclusters > 1]
      c(
        sprintf(
          "- `%s`: %d subclusters from %d parent clusters. Split: %s.",
          sc,
          nrow(d),
          nrow(per),
          if (nrow(multi)) {
            paste0(multi$parent, " into ", multi$subclusters, collapse = ", ")
          } else {
            "none"
          }
        ),
        "",
        fig(paste0("figures/05_subclustering/umap_", sc, ".png"), sc)
      )
    }))
  }
})

## 6. Annotation ##
untyped <- split_arg(opt$untyped)
s6 <- guard({
  ct <- read_tables("^cell_type_counts\\.csv$")
  ## Runs annotated before the refactor wrote the count as Freq ##
  if (has(ct, "Freq") && !has(ct, "n_cells")) {
    setnames(ct, "Freq", "n_cells")
  }
  cb <- read_tables("^composition_by_sample\\.csv$")
  pc <- read_tables("^composition_paired_change\\.csv$")
  if (has(ct, c("cell_type", "n_cells"))) {
    ct[, pct := 100 * n_cells / sum(n_cells)]
    setorder(ct, -n_cells)
    c(
      if (nzchar(basis)) sprintf("Labels on `%s`.", basis),
      "",
      md_table(ct, 1L),
      if (
        has(cb, c("sample_id", "cell_type", "share_of_all")) && length(untyped)
      ) {
        u <- cb[
          cell_type %in% untyped,
          .(untyped_pct = 100 * sum(share_of_all)),
          by = sample_id
        ]
        c("Untyped share per sample:", "", md_table(u, 1L))
      },
      if (has(pc, c("comparison", "cell_type", "log2_change"))) {
        chg <- pc[,
          .(
            subjects = .N,
            mean_log2_change = mean(log2_change),
            same_direction = abs(sum(sign(log2_change))) == .N
          ),
          by = .(comparison, cell_type)
        ]
        chg <- chg[order(comparison, -abs(mean_log2_change))]
        chg <- chg[, head(.SD, 8), by = comparison]
        c(
          "Largest paired changes in the share of typed cells (log2):",
          "",
          md_table(chg, 2L)
        )
      },
      fig("figures/06_annotation/umap_cell_type.png", "Cell types"),
      fig("figures/06_annotation/marker_dotplot.png", "Marker dotplot"),
      fig("figures/06_annotation/composition_barplot.png", "Composition"),
      if (!is.null(pc)) {
        fig(
          "figures/06_annotation/composition_paired_change.png",
          "Paired change"
        )
      }
    )
  }
})

## 7. Spatial ##
s7 <- guard({
  sn <- read_tables("^spatial_net_summary\\.csv$")
  pe <- read_tables("^proximity_enrichment\\.csv$")
  bk <- read_tables(
    "^banksy_summary\\.csv$",
    colClasses = c(lambda = "character")
  )
  bc <- read_tables("^bcell_summary\\.csv$")
  c(
    if (has(sn, c("n_edges", "cross_sample_edges"))) {
      sprintf(
        "- Delaunay network: %s edges, %s across samples.",
        num(sn$n_edges[1]),
        sn$cross_sample_edges[1]
      )
    },
    if (
      has(pe, c("sample", "unified_int", "type_int", "enrichm", "p.adj_higher"))
    ) {
      top <- pe[type_int == "hetero" & p.adj_higher < 0.05 & enrichm > 0]
      top <- top[,
        .(samples = .N, mean_log2_enrichment = mean(enrichm)),
        by = unified_int
      ][order(-samples, -mean_log2_enrichment)]
      c(
        "- Cell-type pairs enriched as neighbours (FDR < 0.05), most samples first:",
        "",
        md_table(head(top, 10), 2L)
      )
    },
    if (has(bk, c("lambda", "n_niches"))) {
      sprintf(
        "- BANKSY: %s.",
        paste0(bk$n_niches, " niches at lambda ", bk$lambda, collapse = ", ")
      )
    },
    if (!is.null(bc)) c("- B-lineage subclustering:", "", md_table(bc)),
    fig(
      "figures/07_spatial/network/proximity_heatmap.png",
      "Proximity enrichment"
    ),
    if (has(bk, "lambda")) {
      unlist(lapply(bk$lambda, function(l) {
        fig(
          paste0(
            "figures/07_spatial/banksy/niche_celltype_heatmap_lam",
            l,
            ".png"
          ),
          paste("Niches, lambda", l)
        )
      }))
    },
    if (!is.null(sn)) {
      c(
        "",
        "Cell-type maps per sample:",
        "",
        fig_list(paste0(
          "figures/07_spatial/network/spatial_celltype_",
          included,
          ".png"
        ))
      )
    }
  )
})

## 8. Cell-cell communication ##
s8 <- guard({
  li <- read_tables("^liana_ligand_receptor\\.csv$")
  nn_sum <- read_tables("^nichenet_summary\\.csv$")
  nn_act <- read_tables("^nichenet_ligand_activities\\.csv$")
  mi <- read_tables("^misty_improvements\\.csv$")
  sx <- read_tables("^sparkx_svgs\\.csv$")
  nv <- read_tables("^nnsvg_spatially_variable_genes\\.csv$")
  ic <- read_tables("^insitucor_modules\\.csv$")
  c(
    if (
      has(
        li,
        c(
          "source",
          "target",
          "ligand.complex",
          "receptor.complex",
          "aggregate_rank"
        )
      )
    ) {
      c(
        "LIANA, top 10 by aggregate rank (contact-weighted):",
        "",
        md_table(
          head(
            li[
              order(aggregate_rank),
              .(
                source,
                target,
                ligand.complex,
                receptor.complex,
                aggregate_rank
              )
            ],
            10
          ),
          4L
        ),
        fig("figures/08_cci/liana/liana_top_interactions.png", "LIANA")
      )
    },
    if (
      has(nn_act, c("comparison", "receiver", "test_ligand", "aupr_corrected"))
    ) {
      top <- nn_act[
        order(-aupr_corrected),
        .(top_ligands = paste(head(test_ligand, 5), collapse = ", ")),
        by = .(comparison, receiver)
      ]
      c(
        "NicheNet, top ligands per comparison and receiver:",
        "",
        md_table(top),
        fig("figures/08_cci/nichenet/nichenet_ligand_activity.png", "NicheNet")
      )
    } else if (!is.null(nn_sum)) {
      c("NicheNet:", "", md_table(nn_sum))
    },
    if (has(mi, c("target", "measure", "value"))) {
      g <- mi[
        measure == "gain.R2",
        .(mean_gain_r2 = mean(value, na.rm = TRUE)),
        by = target
      ]
      if (nrow(g)) {
        c(
          "MISTy, targets best explained by their surroundings (mean gain in R2):",
          "",
          md_table(head(g[order(-mean_gain_r2)], 10), 2L),
          fig("figures/08_cci/misty/misty_gain_r2.png", "MISTy")
        )
      }
    },
    if (has(sx, c("sample", "significant"))) {
      c(
        sprintf(
          "- SPARK-X spatially variable genes per sample: %s.",
          paste0(
            sx[, .(n = sum(significant)), by = sample][, paste0(
              sample,
              " ",
              n
            )],
            collapse = ", "
          )
        ),
        fig("figures/08_cci/sparkx/sparkx_svg_counts.png", "SPARK-X")
      )
    },
    if (has(nv, c("sample", "padj"))) {
      sprintf(
        "- nnSVG genes at adjusted p < 0.05: %s.",
        paste0(
          nv[, .(n = sum(padj < 0.05, na.rm = TRUE)), by = sample][,
            paste0(sample, " ", n)
          ],
          collapse = ", "
        )
      )
    },
    if (has(ic, "module")) {
      sprintf("- InSituCor: %d modules.", uniqueN(ic$module))
    }
  )
})

## 9. Differential expression ##
s9 <- guard({
  sm <- read_tables("^smide_meta_summary\\.csv$")
  smd <- read_tables("^smide_meta\\.csv$")
  c(
    if (!is.null(sm)) {
      c("smiDE meta-analysis, genes at FDR below the cutoff:", "", md_table(sm))
    },
    if (
      has(smd, c("cell_type", "result_component", "target", "log_fc", "fdr"))
    ) {
      top <- smd[
        result_component == "one.vs.rest" & fdr < 0.05,
        .(top_genes = paste(head(target[order(fdr)], 8), collapse = ", ")),
        by = cell_type
      ]
      if (nrow(top)) {
        c(
          "Top genes per cell type (niche against the rest):",
          "",
          md_table(top)
        )
      }
    },
    if (!is.null(sm)) {
      fig("figures/09_de/smide_meta_volcano.png", "smiDE volcano")
    }
  )
})

## Write ##
out <- c(
  header,
  if (!is.null(stages)) c("## Stages run", "", md_table(stages)),
  section("1. Samples and QC", s1),
  section("2. Integration", s2),
  section("3. Clustering", s3),
  section("4. Reference typing", s4),
  section("5. Subclustering", s5),
  section("6. Annotation", s6),
  section("7. Spatial", s7),
  section("8. Cell-cell communication", s8),
  section("9. Differential expression", s9)
)
writeLines(out, file.path(opt$outdir, "run_report.md"))
message(
  "run_report.R: ",
  length(out),
  " lines written after the ",
  opt$stage,
  " stage"
)
