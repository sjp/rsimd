# S7: the Box-Cox profile log-likelihood over a grid of lambda, for a
# positive sample. Full demo: Rscript bench/scenarios/07-boxcox.R
#
# Sample:
#   s <- source("bench/scenarios/07-boxcox.R")$value
#   d <- s$setup(1e5); s$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

boxcox_grid <- seq(-2, 2, by = 0.1)

# l(lambda) = -n/2 log(var(y_lambda)) + (lambda - 1) sum(log x), with
# var(y_lambda) = var(x^lambda) / lambda^2 (the "- 1" of the transform does
# not change the variance), and var(log x) at lambda = 0.
boxcox_base <- function(d) {
  x <- d$x
  n <- length(x)
  slog <- sum(log(x))
  vapply(d$grid, function(l) {
    v <- if (l == 0) var(log(x)) else var(x^l) / l^2
    -n / 2 * log(v) + (l - 1) * slog
  }, 0)
}
boxcox_idiomatic <- function(d) {
  x <- d$x
  n <- length(x)
  slog <- sum(log(x))
  vapply(d$grid, function(l) {
    y <- if (l == 0) log(x) else (x^l - 1) / l
    -n / 2 * log(var(y)) + (l - 1) * slog
  }, 0)
}
boxcox_rsimd <- function(d) {
  x <- d$x
  n <- length(x)
  lx <- simd_log(x)
  slog <- simd_sum(lx)
  vapply(d$grid, function(l) {
    v <- if (l == 0) simd_var(lx) else simd_var(simd_pow(x, l)) / l^2
    -n / 2 * log(v) + (l - 1) * slog
  }, 0)
}

bc <- scenario(
  id = "boxcox", family = "likelihood", flagship = TRUE,
  title = "Box-Cox profile log-likelihood over 41 values of lambda",
  why = "Choosing lambda means one transform and one variance per grid point.",
  sizes = c(1e3, 1e4, 1e5),
  setup = function(n) list(x = rgamma(n, 3), grid = boxcox_grid),
  base = boxcox_base,
  base_idiomatic = boxcox_idiomatic,
  rsimd = boxcox_rsimd,
  # MASS::boxcox() works through lm() and scales by the geometric mean, so
  # its values differ from these by a constant per n; it is timed, and its
  # error is shown but not checked.
  refs = list(MASS = list(pkg = "MASS", fun = function(d) {
    MASS::boxcox(d$x ~ 1, lambda = d$grid, plotit = FALSE)$y
  })),
  tolerance = 1e-10
)

# Full demo: the grid at n = 1e6, with the chosen lambda.
if (sys.nframe() == 0L) {
  set.seed(20261006)
  d <- bc$setup(1e6)
  tb <- system.time(lb <- bc$base(d))[["elapsed"]]
  tr <- system.time(lr <- bc$rsimd(d))[["elapsed"]]
  cat(sprintf(
    "Box-Cox, n = 1e6, 41 lambdas: base %.2f s, rsimd (%s) %.2f s, %.1fx; lambda %.1f / %.1f\n",
    tb, simd_current(), tr, tb / tr, d$grid[which.max(lb)], d$grid[which.max(lr)]
  ))
}

invisible(bc)
