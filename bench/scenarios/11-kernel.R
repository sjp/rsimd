# S12: leave-one-out likelihood cross-validation of a Gaussian kernel
# density estimate's bandwidth; C1 (control): an RBF kernel matrix.
#
# Sample:
#   s <- source("bench/scenarios/11-kernel.R")$value
#   d <- s[[1]]$setup(2000); s[[1]]$rsimd(d)
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

# sum_i log( sum_{j != i} K((x_i - x_j) / h) / ((n - 1) h) ), K Gaussian:
# per point, sum(exp(c0 * d^2)) - 1 removes the j = i term.
kde_setup <- function(n) list(x = c(rnorm(n / 2), rnorm(n / 2, 3, 0.5)), h = 0.3)
kde_const <- function(d) c(c0 = -0.5 / d$h^2, lc = log((length(d$x) - 1) * d$h * sqrt(2 * pi)))

invisible(list(
  scenario(
    id = "kde_cv", family = "kernel",
    title = "KDE leave-one-out likelihood CV, one bandwidth (O(n^2), loop over points)",
    why = "Base R selects bandwidths by bw.ucv()/bw.SJ(), not by likelihood CV.",
    sizes = c(1e3, 2e3, 4e3),
    quick_sizes = c(1e3, 2e3),
    setup = kde_setup,
    base = function(d) {
      k <- kde_const(d)
      x <- d$x
      s <- 0
      for (i in seq_along(x)) {
        dd <- x - x[i]
        s <- s + log(sum(exp(dd * dd * k[["c0"]])) - 1)
      }
      s - length(x) * k[["lc"]]
    },
    # The n x n form, which fits in memory at these sizes.
    base_idiomatic = function(d) {
      k <- kde_const(d)
      D <- outer(d$x, d$x, "-")
      sum(log(rowSums(exp(D * D * k[["c0"]])) - 1)) - length(d$x) * k[["lc"]]
    },
    rsimd = function(d) {
      k <- kde_const(d)
      x <- d$x
      s <- 0
      for (i in seq_along(x)) {
        dd <- simd_sub(x, x[i])
        s <- s + log(simd_sum(simd_exp(simd_mul(simd_mul(dd, dd), k[["c0"]]))) - 1)
      }
      s - length(x) * k[["lc"]]
    },
    tolerance = 1e-12,
    gap = list(id = "G6", ideal = function(d) for (i in seq_along(d$x)) simd_sum(d$x))
  ),
  # C1: exp(-gamma D^2) of a matrix of squared distances. Kept as a control:
  # where the C library's exp is about as fast as SLEEF's (arm64), the
  # extra multiply pass and its allocation make rsimd slower. n is the side.
  scenario(
    id = "rbf_matrix", family = "kernel", control = TRUE,
    title = "RBF kernel matrix exp(-gamma D^2) (control; n = matrix side)",
    why = "One exp per element: the gain is SLEEF's exp against the C library's, nothing else.",
    sizes = c(500, 1000, 2000),
    setup = function(n) {
      x <- matrix(rnorm(n * 3), n, 3)
      list(D2 = as.matrix(dist(x))^2, gamma = 0.5)
    },
    base = function(d) exp(-d$gamma * d$D2),
    rsimd = function(d) {
      k <- simd_exp(simd_mul(d$D2, -d$gamma))
      dim(k) <- dim(d$D2)
      k
    },
    tolerance = 1e-14,
    gap = list(id = "G2", ideal = function(d) simd_exp(d$D2))
  )
))
