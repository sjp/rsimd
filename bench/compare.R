# Compares two benchmark runs written by bench/run.R: for every op x type x
# size x implementation present in both, prints the old and new median times
# and their ratio (old / new, so above 1 means the new run is faster), marking
# changes beyond the threshold. Informational: always exits with status 0.
#
# Usage:
#   Rscript bench/compare.R old.rds new.rds [--threshold 15] [--md] [--force]
#
#   --threshold  percentage change to mark (default 15)
#   --md         print a markdown table (for a pull request or step summary)
#   --force      compare runs from different machines anyway

args <- commandArgs(TRUE)
threshold <- 15
markdown <- FALSE
force <- FALSE
files <- character()
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (a == "--md") {
    markdown <- TRUE
  } else if (a == "--force") {
    force <- TRUE
  } else if (a == "--threshold" || startsWith(a, "--threshold=")) {
    if (a == "--threshold") {
      i <- i + 1L
      a <- if (i <= length(args)) args[[i]] else ""
    }
    threshold <- suppressWarnings(as.numeric(sub("^--threshold=", "", a)))
    if (is.na(threshold) || threshold <= 0) {
      stop("--threshold must be a positive number", call. = FALSE)
    }
  } else if (startsWith(a, "--")) {
    stop("unknown argument '", a, "'", call. = FALSE)
  } else {
    files <- c(files, a)
  }
  i <- i + 1L
}
if (length(files) != 2L) {
  cat("Usage: Rscript bench/compare.R old.rds new.rds [--threshold 15] [--md] [--force]\n")
  quit(status = 0L)
}

old <- readRDS(files[[1L]])
new <- readRDS(files[[2L]])

if (!identical(old$metadata$machine_key, new$metadata$machine_key)) {
  msg <- sprintf(
    "the runs are from different machines:\n  old: %s\n  new: %s\n",
    old$metadata$machine_key, new$metadata$machine_key
  )
  if (!force) {
    cat("Not comparing: ", msg, "Use --force to compare them anyway.\n", sep = "")
    quit(status = 0L)
  }
  cat("Warning: ", msg, sep = "")
}

key_cols <- c("table", "op", "type", "mode", "n", "impl")
both <- merge(
  old$results[c(key_cols, "median_us")], new$results[c(key_cols, "median_us")],
  by = key_cols, suffixes = c("_old", "_new")
)
if (nrow(both) == 0L) {
  cat("The runs have no op x type x size x implementation in common.\n")
  quit(status = 0L)
}
table_order <- c("main", "na", "precision", "math")
both <- both[order(
  match(both$table, table_order), both$op, both$type, both$mode, both$n, both$impl
), ]
both$ratio <- both$median_us_old / both$median_us_new
limit <- threshold / 100
both$mark <- ifelse(both$ratio > 1 + limit, "▲ faster",
  ifelse(both$ratio < 1 / (1 + limit), "▼ slower", "")
)

fmt_time <- function(us) {
  ifelse(us < 1e3, sprintf("%s µs", signif(us, 3)),
    ifelse(us < 1e6, sprintf("%s ms", signif(us / 1e3, 3)), sprintf("%s s", signif(us / 1e6, 3)))
  )
}
out <- data.frame(
  table = both$table, op = both$op, type = both$type,
  mode = both$mode, n = sub("e\\+?0*", "e", format(both$n, scientific = TRUE)),
  impl = both$impl, old = fmt_time(both$median_us_old), new = fmt_time(both$median_us_new),
  ratio = sprintf("%.2f", both$ratio), change = both$mark, stringsAsFactors = FALSE
)

describe <- function(m) {
  sprintf(
    "%s, rsimd %s, git %s", m$timestamp, m$rsimd_version,
    if (is.na(m$git_sha)) "unknown" else substr(m$git_sha, 1L, 12L)
  )
}
summary_line <- sprintf(
  "%d of %d timings changed by more than %g%%: %d faster, %d slower (ratio = old / new).",
  sum(nzchar(both$mark)), nrow(both), threshold,
  sum(both$ratio > 1 + limit), sum(both$ratio < 1 / (1 + limit))
)

if (markdown) {
  cat(sprintf("Old: %s  \nNew: %s\n\n", describe(old$metadata), describe(new$metadata)))
  cat("| ", paste(names(out), collapse = " | "), " |\n", sep = "")
  cat("|", paste(rep("---", ncol(out)), collapse = "|"), "|\n", sep = "")
  for (r in seq_len(nrow(out))) {
    cat("| ", paste(unlist(out[r, ]), collapse = " | "), " |\n", sep = "")
  }
  cat("\n", summary_line, "\n", sep = "")
} else {
  cat(sprintf("Old: %s\nNew: %s\n\n", describe(old$metadata), describe(new$metadata)))
  print(out, row.names = FALSE, right = FALSE)
  cat("\n", summary_line, "\n", sep = "")
}
