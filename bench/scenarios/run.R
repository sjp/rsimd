# Statistical scenarios: calculations base R has no function for, timed as
# idiomatic base R code and as the same algebra with rsimd primitives, on
# every implementation tier, with memory and the error against base (or an
# exact reference). See README.md in this folder.
#
# Usage:
#   Rscript bench/scenarios/run.R [--quick] [--only id,id] [--family f]
#                                 [--flagship] [--out dir]
#
#   --quick     each scenario's quick sizes and shorter timing budgets
#               (about three minutes)
#   --only      comma-separated scenario ids
#   --family    one family: descriptive, likelihood, ml, montecarlo,
#               kernel, distance, rolling, accuracy
#   --flagship  the flagship scenarios only
#   --out       output directory (default: bench/results/)
#
# Writes scenarios-<time>-<machine>-<sha>.{rds,csv,md} and
# scenarios-latest.md. Exits with status 1 if a column's error exceeds its
# scenario's tolerance (the numbers are then not to be trusted).

suppressPackageStartupMessages(library(rsimd))
if (!requireNamespace("bench", quietly = TRUE)) {
  stop("the bench package is required: install.packages(\"bench\")", call. = FALSE)
}

script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    return(file.path(getwd(), "bench", "scenarios"))
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
}
here <- script_dir()
source(file.path(here, "common.R"))
source(file.path(dirname(here), "meta.R"))

# ---------------------------------------------------------------------------
# Command line

parse_args <- function(args) {
  opts <- list(quick = FALSE, only = NULL, family = NULL, flagship = FALSE, out = NULL)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    key <- sub("=.*$", "", a)
    value <- function() {
      if (grepl("=", a, fixed = TRUE)) {
        return(sub("^[^=]*=", "", a))
      }
      i <<- i + 1L
      if (i > length(args)) stop("missing value for ", a, call. = FALSE)
      args[[i]]
    }
    if (key == "--quick") {
      opts$quick <- TRUE
    } else if (key == "--flagship") {
      opts$flagship <- TRUE
    } else if (key == "--only") {
      opts$only <- strsplit(value(), ",", fixed = TRUE)[[1L]]
    } else if (key == "--family") {
      opts$family <- value()
    } else if (key == "--out") {
      opts$out <- value()
    } else if (key %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript bench/scenarios/run.R [--quick] [--only a,b] [--family f]",
        "[--flagship] [--out dir]\n"
      )
      quit(status = 0L)
    } else {
      stop("unknown argument '", a, "'", call. = FALSE)
    }
    i <- i + 1L
  }
  opts
}

opts <- parse_args(commandArgs(TRUE))
scenarios <- load_scenarios(here)
if (!is.null(opts$only)) {
  unknown <- setdiff(opts$only, names(scenarios))
  if (length(unknown)) {
    stop("unknown scenario(s): ", paste(unknown, collapse = ", "), "; available: ",
      paste(names(scenarios), collapse = ", "),
      call. = FALSE
    )
  }
  scenarios <- scenarios[opts$only]
}
if (!is.null(opts$family)) scenarios <- Filter(function(s) s$family == opts$family, scenarios)
if (opts$flagship) scenarios <- Filter(function(s) isTRUE(s$flagship), scenarios)
if (length(scenarios) == 0L) stop("no scenarios selected", call. = FALSE)

tiers <- simd_available()
seed <- 20261006L

# ---------------------------------------------------------------------------
# Running

# The tolerance of a column: a scalar applies to base, idiomatic and rsimd;
# a named vector gives it by precision mode (and "base", "idiomatic"), a
# missing name meaning the error is shown but not checked.
col_tolerance <- function(s, kind, precision) {
  tol <- s$tolerance
  if (is.null(names(tol))) {
    return(if (kind == "ref") Inf else tol[[1L]])
  }
  key <- switch(kind,
    rsimd = precision,
    base = "base",
    idiomatic = "idiomatic",
    ref = ""
  )
  if (key %in% names(tol)) tol[[key]] else Inf
}

rows <- list()
failures <- character()
gaps <- list()
started <- Sys.time()

add_row <- function(s, n, column, kind, tier, precision, t, err, tol) {
  rows[[length(rows) + 1L]] <<- data.frame(
    scenario = s$id, family = s$family, n = n, column = column, kind = kind,
    tier = tier, precision = precision, median_us = t$median_us, mem_bytes = t$mem_bytes,
    error = err, tolerance = tol, stringsAsFactors = FALSE
  )
  if (!is.na(err) && err > tol) {
    failures <<- c(failures, sprintf(
      "%s n=%s %s: error %.2e exceeds tolerance %.1e", s$id, fmt_n(n), column, err, tol
    ))
  }
}

measure <- function(s, f, d, want, part = NULL) {
  t <- time_cell(f, d, opts$quick)
  w <- if (is.null(part)) want else part(want, d)
  err_fun <- if (is.null(s$error)) rel_error else s$error
  t$error <- err_fun(t$value, w)
  t
}

for (s in scenarios) {
  sizes <- if (opts$quick) s$quick_sizes else s$sizes
  message(sprintf("%s (%s)", s$id, paste(fmt_n(sizes), collapse = ", ")))
  for (n in sizes) {
    set.seed(seed)
    d <- s$setup(n)
    base_value <- s$base(d)
    want <- if (is.null(s$reference)) base_value else s$reference(d)
    t <- measure(s, s$base, d, want)
    add_row(s, n, "base", "base", NA, NA, t, t$error, col_tolerance(s, "base", NA))
    if (!is.null(s$base_idiomatic)) {
      t <- measure(s, s$base_idiomatic, d, want)
      add_row(
        s, n, "idiomatic", "idiomatic", NA, NA, t, t$error,
        col_tolerance(s, "idiomatic", NA)
      )
    }
    for (r in names(s$refs)) {
      ref <- s$refs[[r]]
      if (!requireNamespace(ref$pkg, quietly = TRUE)) next
      t <- measure(s, ref$fun, d, base_value, ref$part)
      add_row(s, n, r, "ref", NA, NA, t, t$error, Inf)
    }
    for (p in s$precision) {
      old <- simd_precision(p)
      for (tier in tiers) {
        t <- simd_with_impl(tier, measure(s, s$rsimd, d, want))
        column <- if (length(s$precision) > 1L) paste0(tier, "/", p) else tier
        add_row(s, n, column, "rsimd", tier, p, t, t$error, col_tolerance(s, "rsimd", p))
      }
      simd_precision(old)
    }
  }
  # The cost of composing a missing primitive: the rsimd composition on the
  # auto tier against one pass standing in for it, at the largest n.
  if (!is.null(s$gap)) {
    n <- max(sizes)
    set.seed(seed)
    d <- s$setup(n)
    comp <- time_cell(s$rsimd, d, opts$quick)
    ideal <- time_cell(s$gap$ideal, d, opts$quick)
    gaps[[length(gaps) + 1L]] <- data.frame(
      gap = s$gap$id, scenario = s$id, n = n, composed_us = comp$median_us,
      ideal_us = ideal$median_us, composed_mem = comp$mem_bytes, ideal_mem = ideal$mem_bytes,
      stringsAsFactors = FALSE
    )
  }
}
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

res <- do.call(rbind, rows)
key <- paste(res$scenario, res$n)
base_t <- res$median_us[res$kind == "base"][match(key, key[res$kind == "base"])]
res$speedup <- base_t / res$median_us
gap_tab <- if (length(gaps)) do.call(rbind, gaps) else NULL

# ---------------------------------------------------------------------------
# Output

mi <- machine_info(dirname(here))
out_dir <- if (is.null(opts$out)) file.path(dirname(here), "results") else opts$out
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
slug <- tolower(gsub("[^A-Za-z0-9]+", "-", paste(mi$sysinfo[["sysname"]], mi$sysinfo[["machine"]])))
stem <- sprintf(
  "scenarios-%s-%s-%s", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), slug,
  if (is.na(mi$sha)) "nogit" else substr(mi$sha, 1L, 7L)
)
metadata <- c(mi, list(
  rsimd_version = as.character(utils::packageVersion("rsimd")),
  bench_version = as.character(utils::packageVersion("bench")),
  r_version = R.version.string, available = tiers,
  auto = simd_with_impl("auto", as.character(simd_current())),
  quick = opts$quick, seed = seed, elapsed_sec = elapsed
))
saveRDS(
  list(results = res, gaps = gap_tab, metadata = metadata),
  file.path(out_dir, paste0(stem, ".rds"))
)
utils::write.csv(res, file.path(out_dir, paste0(stem, ".csv")), row.names = FALSE)

auto <- metadata$auto
family_titles <- c(
  descriptive = "Descriptive statistics", likelihood = "Likelihoods, MCMC and MLE",
  ml = "Losses and log-sum-exp", montecarlo = "Monte Carlo weights", kernel = "Kernel methods",
  distance = "Distance and similarity search", rolling = "Rolling and cumulative statistics",
  accuracy = "Accuracy at scale"
)
fmt_x <- function(r) ifelse(is.na(r), "", ifelse(r < 0.95, sprintf("%.2f×", r), sprintf("%.1f×", r)))

cell <- function(r) {
  parts <- fmt_time(r$median_us)
  if (r$kind != "base") parts <- paste0(parts, " (", fmt_x(r$speedup), ")")
  paste0(
    parts, "<br>", fmt_mem(r$mem_bytes), ", err ", fmt_err(r$error),
    if (!is.na(r$error) && r$error > r$tolerance) " **FAIL**" else ""
  )
}

scenario_md <- function(s) {
  sub <- res[res$scenario == s$id, , drop = FALSE]
  cols <- unique(sub$column)
  head <- c(
    paste0("### ", s$title, if (isTRUE(s$flagship)) " ★" else "", if (isTRUE(s$control)) " (control)" else ""),
    "", paste0("`", s$id, "` (", s$family, "). ", s$why), ""
  )
  tab <- c(
    paste0("| n | ", paste(cols, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(cols) + 1L), collapse = "|"), "|")
  )
  for (n in unique(sub$n)) {
    cells <- vapply(cols, function(cc) {
      r <- sub[sub$n == n & sub$column == cc, , drop = FALSE]
      if (nrow(r) == 0L) "" else cell(r[1L, ])
    }, "")
    tab <- c(tab, paste0("| ", fmt_n(n), " | ", paste(cells, collapse = " | "), " |"))
  }
  # Crossover: the smallest n at which the best tier (first precision mode)
  # beats base.
  rs <- sub[sub$kind == "rsimd" & (is.na(sub$precision) | sub$precision == s$precision[1L]), ]
  best <- tapply(rs$speedup, rs$n, max)
  wins <- as.numeric(names(best))[best > 1]
  cross <- if (length(wins)) {
    sprintf("Crossover: rsimd beats base from n = %s on.", fmt_n(min(wins)))
  } else {
    "Crossover: rsimd does not beat base at the sizes run."
  }
  c(head, tab, "", cross, "")
}

caveat <- if (mi$sysinfo[["machine"]] %in% c("aarch64", "arm64") && mi$sysinfo[["sysname"]] == "Linux") {
  paste(
    "**Platform caveat:** on arm64 Linux base R's `long double` is 128-bit IEEE, computed in",
    "software, so base `sum()`, `mean()` and `var()` are far slower than on x86-64 (80-bit x87).",
    "Speedups of reduction-heavy scenarios are inflated here; quote x86-64 numbers."
  )
} else {
  ""
}

flag <- Filter(function(s) isTRUE(s$flagship), scenarios)
flag_md <- character()
if (length(flag)) {
  flag_md <- c(
    "## Flagships (largest n, auto tier)", "",
    "| Scenario | n | base | rsimd | speedup | memory base → rsimd | error |",
    "|---|---|---|---|---|---|---|"
  )
  for (s in flag) {
    sub <- res[res$scenario == s$id, ]
    n <- max(sub$n)
    b <- sub[sub$n == n & sub$kind == "base", ]
    r <- sub[sub$n == n & sub$kind == "rsimd" & sub$tier == auto &
      (sub$precision == s$precision[1L]), ][1L, ]
    flag_md <- c(flag_md, sprintf(
      "| %s | %s | %s | %s | %s | %s → %s | %s |", s$title, fmt_n(n), fmt_time(b$median_us),
      fmt_time(r$median_us), fmt_x(r$speedup), fmt_mem(b$mem_bytes), fmt_mem(r$mem_bytes),
      fmt_err(r$error)
    ))
  }
  flag_md <- c(flag_md, "")
}

gap_md <- character()
if (!is.null(gap_tab)) {
  gap_md <- c(
    "## Cost of composing missing primitives", "",
    paste(
      "Each gap's scenario on the auto tier at its largest n, against the passes over the data",
      "that a fused primitive would make, timed as simd_sum() (memory traffic only, no",
      "arithmetic). For compositions dominated by exp or log the ratio overstates what fusing",
      "could save; the memory column is what it would save. README.md describes each gap."
    ), "",
    "| Gap | Scenario | n | composed | one pass | ratio | memory composed → one pass |",
    "|---|---|---|---|---|---|---|"
  )
  for (k in seq_len(nrow(gap_tab))) {
    g <- gap_tab[k, ]
    gap_md <- c(gap_md, sprintf(
      "| %s | %s | %s | %s | %s | %.1f× | %s → %s |", g$gap, g$scenario, fmt_n(g$n),
      fmt_time(g$composed_us), fmt_time(g$ideal_us), g$composed_us / g$ideal_us,
      fmt_mem(g$composed_mem), fmt_mem(g$ideal_mem)
    ))
  }
  gap_md <- c(gap_md, "")
}

ref_pkgs <- unique(unlist(lapply(scenarios, function(s) vapply(s$refs, function(r) r$pkg, ""))))
ref_pkgs <- ref_pkgs[vapply(ref_pkgs, requireNamespace, NA, quietly = TRUE)]

md <- c(
  "# rsimd statistical scenarios", "",
  "| | |", "|---|---|",
  sprintf("| Timestamp (UTC) | %s |", mi$timestamp),
  sprintf("| rsimd | %s (git %s) |", metadata$rsimd_version, if (is.na(mi$sha)) "unknown" else mi$sha),
  sprintf("| R | %s |", R.version.string),
  sprintf("| OS | %s %s (%s) |", mi$sysinfo[["sysname"]], mi$sysinfo[["release"]], mi$sysinfo[["machine"]]),
  sprintf("| CPU | %s, %s cores |", if (is.na(mi$model)) "unknown" else mi$model, mi$cores),
  sprintf("| CPU features | %s |", paste(mi$features, collapse = ", ")),
  sprintf("| Tiers | %s (auto: %s) |", paste(tiers, collapse = ", "), auto),
  sprintf("| Reference packages | %s |", paste(ref_pkgs, collapse = ", ")),
  sprintf("| Run | %s, %.0f s |", if (opts$quick) "quick" else "full", elapsed), "",
  paste(
    "Cells: median time per call, then (for all but base) the speedup over the fair base R",
    "column (above 1× is faster); memory allocated per call; and the error against base R (or",
    "the exact answer where a scenario has one; normwise relative unless the scenario defines",
    "its own). `idiomatic` is the usual base R form where the fair base uses rsimd's cheaper",
    "algebra. Reference columns come from CRAN packages when installed; their errors are shown,",
    "not checked. Timings are informational and vary between machines and runs."
  ),
  "", caveat, "",
  flag_md,
  unlist(lapply(unique(vapply(scenarios, function(s) s$family, "")), function(f) {
    fs <- Filter(function(s) s$family == f, scenarios)
    c(paste("##", family_titles[[f]]), "", unlist(lapply(fs, scenario_md)))
  })),
  gap_md
)
writeLines(md, file.path(out_dir, paste0(stem, ".md")), useBytes = TRUE)
invisible(file.copy(file.path(out_dir, paste0(stem, ".md")),
  file.path(out_dir, "scenarios-latest.md"),
  overwrite = TRUE
))

message(sprintf(
  "Wrote %s.{rds,csv,md} and scenarios-latest.md to %s (%.0f s)",
  stem, normalizePath(out_dir), elapsed
))
if (length(failures)) {
  message("Errors beyond tolerance:\n  ", paste(failures, collapse = "\n  "))
  quit(status = 1L)
}
