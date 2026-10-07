# S8: log-sum-exp, the normaliser of a log-probability vector, computed
# stably (shifted by the maximum).
#
# Sample:
#   s <- source("bench/scenarios/08-logsumexp.R")$value
#   d <- s$setup(1e5); s$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

invisible(scenario(
  id = "logsumexp", family = "ml",
  title = "log-sum-exp",
  why = "Base R has no logSumExp(); mixture models and HMMs need it everywhere.",
  sizes = c(1e3, 1e4, 1e5, 1e6),
  setup = function(n) list(z = rnorm(n, -50, 10)),
  base = function(d) {
    m <- max(d$z)
    m + log(sum(exp(d$z - m)))
  },
  rsimd = function(d) {
    m <- simd_max(d$z)
    m + log(simd_sum(simd_exp(simd_sub(d$z, m))))
  },
  refs = list(matrixStats = list(pkg = "matrixStats", fun = function(d) {
    matrixStats::logSumExp(d$z)
  })),
  # Base R's sum() adds in long double; on macOS arm64 that is double, and
  # its sequential sum is itself off by up to about 3e-13 at n = 1e6.
  tolerance = if (isTRUE(.Machine$longdouble.digits >= 64L)) 1e-13 else 1e-12,
  gap = list(id = "G3", ideal = function(d) simd_sum(d$z))
))
