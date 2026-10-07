# S4: geometric and harmonic means and the coefficient of variation.
#
# Sample:
#   s <- source("bench/scenarios/04-means.R")$value
#   d <- s$setup(1e5); s$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

invisible(scenario(
  id = "means", family = "descriptive",
  title = "Geometric and harmonic means, coefficient of variation",
  why = "Base R has mean() but no geometric or harmonic mean.",
  sizes = c(1e3, 1e4, 1e5, 1e6),
  setup = function(n) list(x = rlnorm(n)),
  base = function(d) {
    x <- d$x
    c(geometric = exp(mean(log(x))), harmonic = length(x) / sum(1 / x), cv = sd(x) / mean(x))
  },
  rsimd = function(d) {
    x <- d$x
    c(
      geometric = exp(simd_mean(simd_log(x))), harmonic = length(x) / simd_sum(simd_recip(x)),
      cv = simd_sd(x) / simd_mean(x)
    )
  },
  tolerance = 1e-12,
  # One pass per statistic (sum of logs, of reciprocals, the sd).
  gap = list(id = "G6", ideal = function(d) for (k in 1:3) simd_sum(d$x))
))
