# A1 and A2: accuracy at scale. The precision modes buy base R's accuracy
# for a fraction of its time; "fast" is honestly less accurate, and the
# tables show it. Errors are against the exact answer (A1) or base R's
# two-pass var() (A2).
#
# Sample:
#   s <- source("bench/scenarios/16-accuracy.R")$value
#   d <- s[[2]]$setup(1e6); simd_precision("pairwise"); s[[2]]$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

invisible(list(
  scenario(
    id = "sum_cancel", family = "accuracy",
    title = "sum(c(1e16, rep(1, n), -1e16)); the exact answer is n",
    why = "Large terms that cancel: the fast vector sum loses the small ones.",
    sizes = c(1e5, 1e6, 1e7),
    quick_sizes = c(1e5, 1e6),
    setup = function(n) list(x = c(1e16, rep(1, n), -1e16), n = n),
    base = function(d) sum(d$x),
    rsimd = function(d) simd_sum(d$x),
    reference = function(d) d$n,
    precision = c("fast", "pairwise", "compensated"),
    # fast and pairwise are shown, not checked; base R (long double) and
    # compensated are exact here.
    tolerance = c(base = 0, fast = Inf, pairwise = Inf, compensated = 0)
  ),
  scenario(
    id = "var_offset", family = "accuracy", flagship = TRUE,
    title = "Variance of rnorm(n) + 1e9",
    why = "A large offset: one-pass formulas lose everything, two-pass sums lose digits in fast mode.",
    sizes = c(1e5, 1e6),
    setup = function(n) list(x = rnorm(n) + 1e9),
    base = function(d) var(d$x),
    # The "textbook" one-pass formula many hand-written snippets use.
    base_idiomatic = function(d) {
      n <- length(d$x)
      (sum(d$x^2) - n * mean(d$x)^2) / (n - 1)
    },
    rsimd = function(d) simd_var(d$x),
    # Errors are against base R's two-pass var(), accurate here (its sums
    # run in long double).
    precision = c("fast", "pairwise", "compensated"),
    tolerance = c(fast = Inf, pairwise = 1e-12, compensated = 1e-12)
  )
))
