# S1: skewness and kurtosis of a sample. Base R has var() and sd() but no
# higher moments.
#
# Sample (rsimd and base R only):
#   s <- source("bench/scenarios/01-moments.R")$value
#   d <- s$setup(1e5); s$rsimd(d); s$base(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

invisible(scenario(
  id = "moments", family = "descriptive", flagship = TRUE,
  title = "Skewness and kurtosis",
  why = "Base R has var() and sd() but no skewness or kurtosis.",
  sizes = c(1e3, 1e4, 1e5, 1e6),
  setup = function(n) list(x = rgamma(n, 2)),
  # Both sides: the deviations d, their squares d2, then sum(d2 * d) and
  # sum(d2^2) (fairness rule 1: the cubes are products on both sides).
  base = function(d) {
    x <- d$x
    n <- length(x)
    dev <- x - mean(x)
    d2 <- dev * dev
    s2 <- mean(d2)
    c(skewness = sum(d2 * dev) / n / s2^1.5, kurtosis = sum(d2 * d2) / n / s2^2)
  },
  base_idiomatic = function(d) {
    x <- d$x
    m <- mean(x)
    s2 <- mean((x - m)^2)
    c(skewness = mean((x - m)^3) / s2^1.5, kurtosis = mean((x - m)^4) / s2^2)
  },
  rsimd = function(d) {
    x <- d$x
    n <- length(x)
    dev <- simd_sub(x, simd_mean(x))
    d2 <- simd_mul(dev, dev)
    s2 <- simd_mean(d2)
    c(skewness = simd_dot(d2, dev) / n / s2^1.5, kurtosis = simd_sum_sq(d2) / n / s2^2)
  },
  refs = list(moments = list(pkg = "moments", fun = function(d) {
    c(skewness = moments::skewness(d$x), kurtosis = moments::kurtosis(d$x))
  })),
  tolerance = 1e-12,
  gap = list(id = "G4", ideal = function(d) simd_sum(d$x))
))
