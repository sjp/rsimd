# S15: rolling mean and standard deviation (window 50) from differences of
# cumulative sums; S15b: the trapezoid area under a curve.
#
# The cumulative-sum method cancels badly when the data sit far from zero
# (its sd loses digits on x + 1e6): a known caveat of this formulation, on
# both sides alike.
#
# Sample:
#   s <- source("bench/scenarios/15-rolling.R")$value
#   d <- s[[1]]$setup(1e5); str(s[[1]]$rsimd(d))
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

# The rolling means: the first part of the result.
roll_mean_part <- function(v, d) v[seq_len(length(d$x) - d$k + 1L)]

invisible(list(
  scenario(
    id = "rolling", family = "rolling",
    title = "Rolling mean and sd, window 50 (cumulative-sum differences)",
    why = "Base R has stats::filter() for a rolling mean, nothing for a rolling sd.",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = function(n) list(x = cumsum(rnorm(n)) / 10, k = 50L),
    base = function(d) {
      x <- d$x
      k <- d$k
      n <- length(x)
      c1 <- c(0, cumsum(x))
      c2 <- c(0, cumsum(x * x))
      s1 <- c1[(k + 1):(n + 1)] - c1[1:(n - k + 1)]
      s2 <- c2[(k + 1):(n + 1)] - c2[1:(n - k + 1)]
      c(s1 / k, sqrt(pmax((s2 - s1 * s1 / k) / (k - 1), 0)))
    },
    rsimd = function(d) {
      x <- d$x
      k <- d$k
      n <- length(x)
      c1 <- c(0, simd_cumsum(x))
      c2 <- c(0, simd_cumsum(simd_mul(x, x)))
      s1 <- simd_sub(c1[(k + 1):(n + 1)], c1[1:(n - k + 1)])
      s2 <- simd_sub(c2[(k + 1):(n + 1)], c2[1:(n - k + 1)])
      c(
        simd_mul(s1, 1 / k),
        simd_sqrt(simd_pmax(simd_mul(simd_mul_add(simd_mul(s1, s1), -1 / k, s2), 1 / (k - 1)), 0))
      )
    },
    # stats::filter() gives the rolling mean only (sides = 1 aligns it to
    # the window's end); timed, error shown but not checked.
    refs = list(
      `stats::filter` = list(pkg = "stats", part = roll_mean_part, fun = function(d) {
        as.numeric(stats::filter(d$x, rep(1 / d$k, d$k), sides = 1))[-seq_len(d$k - 1)]
      }),
      `data.table::frollmean` = list(pkg = "data.table", part = roll_mean_part, fun = function(d) {
        data.table::frollmean(d$x, d$k)[-seq_len(d$k - 1)]
      })
    ),
    # Base R's cumsum() accumulates in long double, rsimd's in double (and
    # reassociated in blocks); the sd, a difference of large partial sums,
    # amplifies that by the size of the sums over the window's, which grows
    # with n (about 1e-8 at 1e6, on every tier including none).
    tolerance = 1e-7,
    gap = list(id = "G8", ideal = function(d) simd_cumsum(d$x))
  ),
  scenario(
    id = "auc", family = "rolling",
    title = "Trapezoid area under a curve",
    why = "Pharmacokinetic AUC and ROC AUC; base R has no integrate-by-trapezoid for data.",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = function(n) {
      t <- cumsum(rexp(n))
      list(t = t, y = exp(-t / max(t) * 3) + runif(n, 0, 0.01))
    },
    base = function(d) {
      n <- length(d$t)
      sum((d$t[-1] - d$t[-n]) * (d$y[-1] + d$y[-n])) / 2
    },
    base_idiomatic = function(d) sum(diff(d$t) * (head(d$y, -1) + tail(d$y, -1))) / 2,
    rsimd = function(d) {
      n <- length(d$t)
      simd_dot(simd_sub(d$t[-1], d$t[-n]), simd_add(d$y[-1], d$y[-n])) / 2
    },
    tolerance = 1e-12,
    gap = list(id = "G9", ideal = function(d) simd_dot(d$t, d$y))
  )
))
