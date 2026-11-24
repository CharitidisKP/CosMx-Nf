## Helpers of the run report that need no results folder, sourced by _helpers.R and the tests ##

## A value from nextflow.config without its quotes ##
unquote <- function(x) {
  x <- trimws(x)
  sub("^'(.*)'$", "\\1", sub('^"(.*)"$', "\\1", x))
}

## A setting's comment, without the group header or tool name above it ##
param_doc <- function(raw, at_start) {
  if (!length(raw)) {
    return("")
  }
  txt <- trimws(raw)
  title <- function(x, n) lengths(strsplit(x, "\\s+")) <= n & !grepl("[.,;:()=]", x)
  if (at_start && (length(txt) > 1 || title(txt[1], 4))) {
    raw <- raw[-1]
    txt <- txt[-1]
  }
  ## Indented lines list a file's columns, which stay in the config ##
  ind <- which(grepl("^\\s{2,}", raw))
  if (length(ind)) {
    txt <- txt[seq_len(ind[1] - 1)]
  }
  txt <- txt[!(title(txt, 2) & seq_along(txt) < length(txt))]
  if (!length(txt)) {
    return("")
  }
  head <- seq_len(length(txt) - 1)
  txt[head] <- ifelse(grepl("[.!?:]$", txt[head]), txt[head], paste0(txt[head], "."))
  paste(txt, collapse = " ")
}

## Every setting in the params block of nextflow.config, with its default and its comment ##
config_params <- function(lines) {
  none <- data.frame(name = character(), default = character(), doc = character())
  start <- grep("^params\\s*\\{", lines)[1]
  if (is.na(start)) {
    return(none)
  }
  ends <- grep("^\\}", lines)
  end <- ends[ends > start][1]
  if (is.na(end) || end <= start + 1L) {
    return(none)
  }
  out <- list()
  pending <- character()
  fresh <- TRUE
  at_start <- TRUE
  for (l in lines[(start + 1L):(end - 1L)]) {
    t <- trimws(l)
    if (!nzchar(t)) {
      pending <- character()
      fresh <- TRUE
    } else if (startsWith(t, "//")) {
      if (!length(pending)) at_start <- fresh
      pending <- c(pending, sub("^\\s*//", "", l))
      fresh <- FALSE
    } else {
      m <- regmatches(t, regexec("^([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.*)$", t))[[1]]
      if (length(m)) {
        out[[length(out) + 1L]] <- data.frame(
          name = m[2], default = unquote(m[3]), doc = param_doc(pending, at_start)
        )
      }
      pending <- character()
      fresh <- FALSE
    }
  }
  if (length(out)) do.call(rbind, out) else none
}

## The first rows of a table, ordered by one column, ascending or with a minus sign descending ##
cap_rows <- function(df, cap = 5000L, order = NULL) {
  n <- nrow(df)
  if (n <= cap) {
    return(list(df = df, total = n, by = NULL))
  }
  col <- if (length(order)) sub("^-", "", order) else ""
  if (nzchar(col) && col %in% names(df)) {
    o <- base::order(df[[col]], decreasing = startsWith(order, "-"), na.last = TRUE)
    return(list(df = df[o[seq_len(cap)], , drop = FALSE], total = n, by = col))
  }
  list(df = df[seq_len(cap), , drop = FALSE], total = n, by = NULL)
}

## Columns holding p values or adjusted p values ##
is_p_col <- function(n) {
  x <- tolower(n)
  grepl("(^|[._])p([._]|$)", x) |
    grepl("p\\.?val", x) |
    grepl("(^|[._])(padj|fdr|qvalue)([._]|$)", x)
}

## Text columns holding only decimal numbers, as numbers, except the columns asked for as text ##
as_numbers <- function(d, colClasses = NULL) {
  if (length(colClasses) && is.null(names(colClasses))) {
    return(d)
  }
  keep <- if (is.list(colClasses)) unlist(colClasses[names(colClasses) == "character"]) else names(colClasses)[colClasses == "character"]
  if (is.numeric(keep)) keep <- names(d)[keep]
  for (col in setdiff(names(d)[vapply(d, is.character, logical(1))], keep)) {
    x <- d[[col]]
    has <- !is.na(x) & nzchar(x)
    first <- match(TRUE, has)
    if (is.na(first) || is.na(suppressWarnings(as.numeric(x[first])))) next
    num <- suppressWarnings(as.numeric(x))
    if (anyNA(num[has]) || !any(grepl(".", x[has], fixed = TRUE))) next
    if (data.table::is.data.table(d)) data.table::set(d, j = col, value = num) else d[[col]] <- num
  }
  d
}

## Cluster IDs in order by number at every level, 6.2 before 6.10 and 5.5.2 before 5.5.10, then any names ##
cluster_levels <- function(x) {
  u <- unique(as.character(x[!is.na(x)]))
  parts <- strsplit(u, ".", fixed = TRUE)
  part <- function(i) suppressWarnings(as.numeric(vapply(parts, function(p) if (length(p) >= i) p[i] else "0", "")))
  u[do.call(order, c(lapply(seq_len(max(lengths(parts), 1L)), part), list(u)))]
}

## Legend columns for n entries, 22 to a column, and none past 66 ##
legend_columns <- function(n) if (n > 66) 0L else as.integer(max(1, ceiling(n / 22)))

## A results file as a link from a page one folder below the results root ##
rel_href <- function(path, root) {
  rel <- substring(path, nchar(root) + 2L)
  parts <- strsplit(rel, "/", fixed = TRUE)
  paste0("../", vapply(parts, function(p) {
    paste(vapply(p, utils::URLencode, "", reserved = TRUE), collapse = "/")
  }, ""))
}

## Width and height of a PNG from its header ##
png_size <- function(f) {
  b <- tryCatch(readBin(f, "raw", 24L), error = function(e) raw())
  sig <- as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
  if (length(b) < 24L || !identical(b[1:8], sig)) {
    return(c(NA_integer_, NA_integer_))
  }
  be <- function(x) as.integer(sum(as.numeric(x) * 256^(3:0)))
  c(be(b[17:20]), be(b[21:24]))
}

## The run's launches in Nextflow's history, in order, with the stage each ran ##
history_launches <- function(lines, run_id) {
  none <- data.frame(
    time = character(), duration = character(), name = character(), status = character(),
    session = character(), stage = character(), command = character()
  )
  f <- strsplit(lines, "\t", fixed = TRUE)
  f <- f[lengths(f) >= 7L]
  if (!length(f)) {
    return(none)
  }
  cmd <- vapply(f, function(x) x[length(x)], "")
  keep <- grepl(paste0("--run_id ", run_id, " "), paste0(cmd, " "), fixed = TRUE)
  f <- f[keep]
  cmd <- cmd[keep]
  if (!length(f)) {
    return(none)
  }
  hit <- regexpr("--stage [^ ]+", cmd)
  stage <- rep("discover", length(cmd))
  stage[hit > 0] <- sub("^--stage ", "", regmatches(cmd, hit))
  field <- function(i) vapply(f, function(x) x[i], "")
  data.frame(
    time = field(1L), duration = field(2L), name = field(3L), status = field(4L),
    session = field(6L), stage = stage, command = cmd
  )
}

## A stage's state from its last launch: done, failed, running, stopped or none ##
launch_state <- function(stage, launches, made = FALSE, current = "") {
  st <- launches$status[launches$stage == stage]
  if (!length(st)) {
    return(if (made) "done" else "none")
  }
  last <- st[length(st)]
  if (last == "OK") {
    "done"
  } else if (last == "ERR") {
    "failed"
  } else if (identical(stage, current)) {
    "running"
  } else {
    "stopped"
  }
}

## True when every file a stage wrote is older than the object it read ##
outdated_by <- function(obj_time, out_times) {
  !is.na(obj_time) && length(out_times) > 0L && !any(out_times > obj_time, na.rm = TRUE)
}

## The digest main.nf and cosmx-nf take of a clustering, its labels and its per-cell edits ##
label_digest <- function(basis, labels, overrides) {
  rd <- function(f) readBin(f, "raw", file.size(f))
  tmp <- tempfile()
  on.exit(unlink(tmp))
  writeBin(c(charToRaw(basis), as.raw(10L), rd(labels), as.raw(10L), rd(overrides)), tmp)
  unname(tools::md5sum(tmp))
}

## A launch's own options, without the paths and names cosmx-nf adds ##
stage_options <- function(cmd) {
  tok <- strsplit(trimws(cmd), "\\s+")[[1]]
  m <- match("main.nf", tok)
  if (!is.na(m)) tok <- tok[-seq_len(m)]
  pairs <- c("-log", "-work-dir", "-w", "--outdir", "--run_id", "--stage", "-params-file", "-c", "-config")
  out <- character()
  i <- 1L
  while (i <= length(tok)) {
    t <- tok[i]
    if (t %in% pairs) {
      i <- i + 2L
    } else if (t == "-resume") {
      i <- i + 1L + as.integer(i < length(tok) && !startsWith(tok[i + 1L], "-"))
    } else {
      out <- c(out, t)
      i <- i + 1L
    }
  }
  paste(out, collapse = " ")
}

## The md5 of a text ##
md5_text <- function(x) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  writeLines(as.character(x), tmp, useBytes = TRUE)
  unname(tools::md5sum(tmp))
}

## A figure's key from its name, the size and time of its input files, its code and its size ##
fig_key <- function(name, inputs, code, w, h) {
  info <- file.info(inputs)
  md5_text(c(name, inputs, info$size, format(info$mtime, "%Y%m%d%H%M%OS3"), code, w, h))
}
