# Times the elementary functions whose implementation depends on the tier
# (issue 036): the real functions that the SIMD tiers may hand to the C
# math library, and the complex inverse functions on base R's branch-cut
# lines, where every lane is computed by the scalar code. Each case runs on
# every implementation in simd_available(), in both accuracy modes, and in
# base R, so the table shows where a SIMD tier is slower than `none`.
#
# Usage: Rscript bench/math_paths.R [--n 1e6] [--rounds 3] [--iterations 7]
#                                    [--ops log,asin] [--out dir]
#        (rsimd and bench installed; --ops keeps the cases of those
#        functions, real and complex)
#
# Each timing is the median over `rounds` rounds of the median of
# `iterations` calls; the implementations take turns within a round, so
# slow drift of the machine's speed affects them alike. The table goes to
# stdout and, with a csv of every timing, to <out>/math-paths-<stamp>.md.

library(rsimd)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
n <- as.numeric(opt("n", "1e6"))
rounds <- as.integer(opt("rounds", "3"))
iterations <- as.integer(opt("iterations", "7"))
out_dir <- opt("out", file.path("bench", "results"))
only <- opt("ops", NULL)

set.seed(20261006)
x <- stats::runif(n, -10, 10)
pos <- 10^stats::runif(n, -3, 3)
y <- stats::runif(n, -3, 3)
ch <- stats::runif(n, -20, 20)
ge1 <- 1 + 10^stats::runif(n, -3, 3)
# Base R's cut lines: real |x| > 1 for asin, acos and atanh, imaginary
# |y| > 1 for asinh and atan.
sgn <- sample(c(-1, 1), n, TRUE)
cut_re <- as.complex(sgn * stats::runif(n, 1.01, 10))
cut_im <- complex(real = 0, imaginary = sgn * stats::runif(n, 1.01, 10))
off <- complex(real = stats::runif(n, -3, 3), imaginary = stats::runif(n, -3, 3))
# Mixed: a quarter of the elements on the cut, at random positions.
on_cut <- stats::runif(n) < 0.25
mix_re <- ifelse(on_cut, cut_re, off)
mix_im <- ifelse(on_cut, cut_im, off)

cases <- list(
  list(op = "exp", input = "(-10, 10)", simd = quote(simd_exp(x)), base = quote(exp(x))),
  list(op = "log", input = "10^(-3, 3)", simd = quote(simd_log(pos)), base = quote(log(pos))),
  list(op = "log2", input = "10^(-3, 3)", simd = quote(simd_log2(pos)), base = quote(log2(pos))),
  list(op = "asinh", input = "(-10, 10)", simd = quote(simd_asinh(x)), base = quote(asinh(x))),
  list(
    op = "acosh", input = "1 + 10^(-3, 3)", simd = quote(simd_acosh(ge1)),
    base = quote(acosh(ge1))
  ),
  list(op = "cosh", input = "(-20, 20)", simd = quote(simd_cosh(ch)), base = quote(cosh(ch))),
  # Not on any tier's list; timed so a change in the balance shows.
  list(
    op = "cbrt", input = "(-10, 10)", simd = quote(simd_cbrt(x)),
    base = quote(sign(x) * abs(x)^(1 / 3))
  ),
  list(
    op = "pow", input = "10^(-3, 3) ^ (-3, 3)", simd = quote(simd_pow(pos, y)),
    base = quote(pos^y)
  )
)
for (f in c("asin", "acos", "atanh", "asinh", "atan")) {
  real_cut <- f %in% c("asin", "acos", "atanh")
  for (kind in c("cut", "mixed", "off")) {
    arg <- switch(kind,
      cut = if (real_cut) "cut_re" else "cut_im",
      mixed = if (real_cut) "mix_re" else "mix_im",
      off = "off"
    )
    cases[[length(cases) + 1L]] <- list(
      op = paste0(f, " (complex)"),
      input = switch(kind,
        cut = if (real_cut) "real |x| > 1" else "imaginary |y| > 1",
        mixed = "25% cut",
        off = "off the axes"
      ),
      simd = call(paste0("simd_", f), as.name(arg)),
      base = call(f, as.name(arg))
    )
  }
}

if (!is.null(only)) {
  keep <- strsplit(only, ",", fixed = TRUE)[[1L]]
  cases <- Filter(function(cs) sub(" .*", "", cs$op) %in% keep, cases)
  if (length(cases) == 0L) stop("no case matches --ops ", only, call. = FALSE)
}

tiers <- simd_available()
impls <- c("base", as.vector(t(outer(tiers, c("accurate", "fast"), paste))))
impls <- impls[!grepl("^none fast$", impls)]

# The CPU, as far as the OS says (arm64 Linux gives only the implementer
# and part codes: 0x61 is Apple, 0x41 part 0xd49 Neoverse N2).
cpu_model <- function() {
  if (file.exists("/proc/cpuinfo")) {
    info <- readLines("/proc/cpuinfo", warn = FALSE)
    field <- function(f) {
      hit <- grep(paste0("^", f, "\\s*:"), info, value = TRUE)
      if (length(hit)) trimws(sub("^[^:]*:", "", hit[1L])) else NA_character_
    }
    if (!is.na(field("model name"))) {
      return(field("model name"))
    }
    if (!is.na(field("CPU implementer"))) {
      return(sprintf("implementer %s part %s", field("CPU implementer"), field("CPU part")))
    }
  }
  if (Sys.info()[["sysname"]] == "Darwin") {
    return(system2("sysctl", c("-n", "machdep.cpu.brand_string"), stdout = TRUE))
  }
  "unknown"
}

# Median time of one call of `e` in milliseconds.
time_one <- function(e) {
  m <- bench::mark(eval(e), iterations = iterations, check = FALSE, memory = FALSE)
  as.numeric(m$median) * 1e3
}
time_impl <- function(case, impl) {
  if (impl == "base") {
    return(time_one(case$base))
  }
  parts <- strsplit(impl, " ", fixed = TRUE)[[1L]]
  old <- simd_math_accuracy(parts[2L])
  on.exit(simd_math_accuracy(old))
  simd_with_impl(parts[1L], time_one(case$simd))
}

rows <- list()
for (case in cases) {
  t <- matrix(NA_real_, rounds, length(impls), dimnames = list(NULL, impls))
  for (r in seq_len(rounds)) for (im in impls) t[r, im] <- time_impl(case, im)
  med <- apply(t, 2L, stats::median)
  rows[[length(rows) + 1L]] <- data.frame(
    op = case$op, input = case$input, impl = impls, median_ms = unname(med),
    stringsAsFactors = FALSE
  )
  message(sprintf(
    "%-18s %-22s %s", case$op, case$input,
    paste(sprintf("%s %.2f", impls, med), collapse = ", ")
  ))
}
res <- do.call(rbind, rows)

# One row per case, one column per implementation, then each SIMD tier's
# accurate time over none's (above 1 is slower than none).
wide <- reshape(res, idvar = c("op", "input"), timevar = "impl", direction = "wide")
names(wide) <- sub("^median_ms\\.", "", names(wide))
simd <- setdiff(tiers, "none")
for (tr in simd) {
  wide[[paste(tr, "/ none")]] <- wide[[paste(tr, "accurate")]] / wide[["none accurate"]]
}
fmt <- function(v) formatC(v, format = "f", digits = 2)
md <- c(
  sprintf("# Math paths (n = %s, median of %d x %d calls)", format(n), rounds, iterations),
  "",
  sprintf(
    "| | |\n|---|---|\n| rsimd | %s |\n| R | %s |\n| Machine | %s |\n| CPU | %s |",
    utils::packageVersion("rsimd"), R.version.string,
    paste(Sys.info()[c("sysname", "machine")], collapse = " "),
    cpu_model()
  ),
  "",
  paste0("| ", paste(names(wide), collapse = " | "), " |"),
  paste0("|", strrep("---|", ncol(wide))),
  apply(wide, 1L, function(r) {
    num <- suppressWarnings(as.numeric(r))
    paste0("| ", paste(ifelse(is.na(num), r, fmt(num)), collapse = " | "), " |")
  })
)
writeLines(md)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
writeLines(md, file.path(out_dir, sprintf("math-paths-%s.md", stamp)))
utils::write.csv(res, file.path(out_dir, sprintf("math-paths-%s.csv", stamp)), row.names = FALSE)
