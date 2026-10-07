# S3: weighted mean and variance, and S3b: the Gini coefficient of a
# sorted sample.
#
# Sample:
#   s <- source("bench/scenarios/03-weighted.R")$value
#   d <- s[[1]]$setup(1e5); s[[1]]$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

invisible(list(
  # Frequency weights: the variance divides by sum(w) - 1.
  scenario(
    id = "weighted_var", family = "descriptive",
    title = "Weighted mean and variance",
    why = "Base R has weighted.mean() but no weighted variance (cov.wt() is matrix-only).",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = function(n) list(x = rnorm(n, 10), w = runif(n, 0.5, 3)),
    base = function(d) {
      sw <- sum(d$w)
      m <- sum(d$w * d$x) / sw
      dev <- d$x - m
      c(mean = m, var = sum(d$w * dev * dev) / (sw - 1))
    },
    base_idiomatic = function(d) {
      m <- weighted.mean(d$x, d$w)
      c(mean = m, var = sum(d$w * (d$x - m)^2) / (sum(d$w) - 1))
    },
    rsimd = function(d) {
      sw <- simd_sum(d$w)
      m <- simd_dot(d$w, d$x) / sw
      dev <- simd_sub(d$x, m)
      c(mean = m, var = simd_dot(d$w, simd_mul(dev, dev)) / (sw - 1))
    },
    refs = list(matrixStats = list(pkg = "matrixStats", fun = function(d) {
      c(mean = matrixStats::weightedMean(d$x, d$w), var = matrixStats::weightedVar(d$x, d$w))
    })),
    tolerance = 1e-12,
    gap = list(id = "G5", ideal = function(d) simd_dot(d$w, d$x))
  ),
  # G = sum((2i - n - 1) x_(i)) / (n sum(x)) for sorted x.
  scenario(
    id = "gini", family = "descriptive",
    title = "Gini coefficient (sample already sorted)",
    why = "Base R has no Gini coefficient.",
    sizes = c(1e3, 1e4, 1e5, 1e6),
    setup = function(n) list(xs = sort(rlnorm(n)), i = seq_len(n)),
    base = function(d) {
      n <- length(d$xs)
      sum((2 * d$i - n - 1) * d$xs) / (n * sum(d$xs))
    },
    rsimd = function(d) {
      n <- length(d$xs)
      simd_dot(simd_mul_add(d$i, 2, -n - 1), d$xs) / (n * simd_sum(d$xs))
    },
    tolerance = 1e-12
  ),
  scenario(
    id = "gini_sort", family = "descriptive",
    title = "Gini coefficient including the sort",
    why = "As gini, with the sort() both sides need, which dominates.",
    sizes = c(1e3, 1e5, 1e6),
    setup = function(n) list(x = rlnorm(n)),
    base = function(d) {
      xs <- sort(d$x)
      n <- length(xs)
      sum((2 * seq_len(n) - n - 1) * xs) / (n * sum(xs))
    },
    rsimd = function(d) {
      xs <- sort(d$x)
      n <- length(xs)
      simd_dot(simd_mul_add(seq_len(n), 2, -n - 1), xs) / (n * simd_sum(xs))
    },
    tolerance = 1e-12
  )
))
