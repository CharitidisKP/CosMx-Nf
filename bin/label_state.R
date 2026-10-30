## Label state: the cluster labels of one clustering and the per-cell edits laid over them.
## Base R only. bin/apply_labels.R applies it to the object and shiny/app.R shows it live, so
## the app draws what annotate will write ##

## An edit with this label hands the cell back to its cluster's label ##
KEEP_LABEL <- "(cluster label)"

## A cluster labels CSV (cluster_id,label) as a named vector, cluster_id -> label. Read as text:
## as numbers "3.1" and "3.10" are both 3.1 ##
read_cluster_labels <- function(file) {
  d <- utils::read.csv(file, colClasses = "character")
  if (!all(c("cluster_id", "label") %in% names(d))) {
    stop(file, " needs columns cluster_id and label")
  }
  dup <- unique(d$cluster_id[duplicated(d$cluster_id)])
  if (length(dup)) {
    stop("duplicate cluster_id in ", file, ": ", paste(dup, collapse = ", "))
  }
  d$label[is.na(d$label)] <- ""
  stats::setNames(trimws(d$label), d$cluster_id)
}

empty_overrides <- function() {
  data.frame(
    cell_ID = character(),
    label = character(),
    confidence = character(),
    note = character(),
    stringsAsFactors = FALSE
  )
}

## Label edits as cell_ID, label, confidence (high, the default, or low) and note. Rows keep
## their order: of a cell's edits, the later one wins ##
check_overrides <- function(d, what = "label overrides") {
  if (!all(c("cell_ID", "label") %in% names(d))) {
    stop(what, " needs columns cell_ID and label")
  }
  for (k in c("confidence", "note")) {
    if (!k %in% names(d)) d[[k]] <- rep("", nrow(d))
  }
  for (k in c("cell_ID", "label", "confidence", "note")) {
    v <- as.character(d[[k]])
    v[is.na(v)] <- ""
    d[[k]] <- trimws(v)
  }
  d$confidence <- tolower(d$confidence)
  d$confidence[!nzchar(d$confidence)] <- "high"
  bad <- setdiff(unique(d$confidence), c("high", "low"))
  if (length(bad)) {
    stop(
      "confidence must be high or low in ",
      what,
      ", found: ",
      paste(bad, collapse = ", ")
    )
  }
  d <- d[
    nzchar(d$cell_ID) & nzchar(d$label),
    c("cell_ID", "label", "confidence", "note"),
    drop = FALSE
  ]
  rownames(d) <- NULL
  d
}

read_label_overrides <- function(file) {
  if (is.null(file) || !nzchar(file) || !file.exists(file)) {
    return(empty_overrides())
  }
  check_overrides(utils::read.csv(file, colClasses = "character"), file)
}

## Each cell's latest high-confidence edit, and its latest low-confidence edit when that is newer
## than the high one. These are the only edits that count ##
current_edits <- function(overrides) {
  ov <- overrides
  ov$row <- seq_len(nrow(ov))
  latest <- function(conf) {
    x <- ov[ov$confidence == conf, , drop = FALSE]
    x[!duplicated(x$cell_ID, fromLast = TRUE), , drop = FALSE]
  }
  hi <- latest("high")
  lo <- latest("low")
  newer <- lo$row > hi$row[match(lo$cell_ID, hi$cell_ID)]
  newer[is.na(newer)] <- TRUE
  list(high = hi, low = lo[newer, , drop = FALSE])
}

## Each cell's labels:
##   cell_type          its cluster's label, or its latest high-confidence edit. Every analysis
##                      reads it
##   cell_type_display  cell_type, or its latest low-confidence edit when that is newer than
##                      its latest high one. Figures draw it
##   label_confidence   low where cell_type_display is a low-confidence edit, else high
##   cell_type_source   manual where a high-confidence edit set cell_type, else cluster
## A cluster without a label leaves its cells NA ##
apply_label_state <- function(cell_ids, clusters, cluster_labels, overrides) {
  cell_ids <- as.character(cell_ids)
  ct <- unname(cluster_labels[as.character(clusters)])
  ed <- current_edits(overrides[overrides$cell_ID %in% cell_ids, , drop = FALSE])
  hi <- ed$high[ed$high$label != KEEP_LABEL, , drop = FALSE]
  lo <- ed$low[ed$low$label != KEEP_LABEL, , drop = FALSE]
  i <- match(hi$cell_ID, cell_ids)
  ct[i] <- hi$label
  src <- rep("cluster", length(cell_ids))
  src[i] <- "manual"
  disp <- ct
  j <- match(lo$cell_ID, cell_ids)
  disp[j] <- lo$label
  conf <- rep("high", length(cell_ids))
  conf[j] <- "low"
  data.frame(
    cell_ID = cell_ids,
    cell_type = ct,
    cell_type_display = disp,
    label_confidence = conf,
    cell_type_source = src,
    stringsAsFactors = FALSE
  )
}

## The same edits cut to the rows that count, in their order. Applying the cut gives the same
## labels, so a file of edits stays small however often a cell is relabelled ##
compact_overrides <- function(overrides) {
  ed <- current_edits(overrides)
  rows <- sort(c(ed$high$row, ed$low$row))
  keep <- seq_len(nrow(overrides)) %in% rows
  out <- overrides[keep, c("cell_ID", "label", "confidence", "note"), drop = FALSE]
  out <- out[out$label != KEEP_LABEL, , drop = FALSE]
  rownames(out) <- NULL
  out
}

## Writes through a temporary file in the same folder, so a reader never sees half a file. The
## folder is made when missing ##
write_csv_atomic <- function(d, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".tmp_", tmpdir = dirname(file), fileext = ".csv")
  utils::write.csv(d, tmp, row.names = FALSE)
  if (!file.rename(tmp, file)) {
    unlink(tmp)
    stop("could not write ", file)
  }
  invisible(file)
}
