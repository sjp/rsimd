# S13: the energy distance between two samples, by a loop over points
# (O(n^2)); the n x n form is shown where it fits.
#
# Sample:
#   s <- source("bench/scenarios/12-energy.R")$value
#   d <- s$setup(2000); s$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

# E = 2 mean|a - b| - mean|a - a'| - mean|b - b'|.
energy_loop <- function(d, sumabs) {
  a <- d$a
  b <- d$b
  sab <- 0
  saa <- 0
  sbb <- 0
  for (i in seq_along(a)) {
    sab <- sab + sumabs(b, a[i])
    saa <- saa + sumabs(a, a[i])
  }
  for (i in seq_along(b)) sbb <- sbb + sumabs(b, b[i])
  2 * sab / (length(a) * length(b)) - saa / length(a)^2 - sbb / length(b)^2
}

invisible(scenario(
  id = "energy", family = "distance",
  title = "Energy distance between two samples (O(n^2) loop; n per sample)",
  why = "Base R has no energy statistic; a permutation test repeats it many times.",
  sizes = c(1e3, 2e3, 5e3),
  quick_sizes = c(1e3, 2e3),
  setup = function(n) list(a = rnorm(n), b = rnorm(n, 0.2, 1.1)),
  base = function(d) energy_loop(d, function(v, c) sum(abs(v - c))),
  base_idiomatic = function(d) {
    2 * mean(abs(outer(d$a, d$b, "-"))) - mean(abs(outer(d$a, d$a, "-"))) -
      mean(abs(outer(d$b, d$b, "-")))
  },
  rsimd = function(d) energy_loop(d, function(v, c) simd_sum_abs(simd_sub(v, c))),
  tolerance = 1e-12
))
