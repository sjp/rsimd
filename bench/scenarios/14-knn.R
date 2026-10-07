# C2 and C3 (controls): Euclidean distances from a query to every column
# of a 64 x n matrix, the first step of k-nearest neighbours. rsimd has no
# vector-against-columns form, so C2 slices columns and loses; C3 uses the
# BLAS identity |x|^2 - 2 x'q + |q|^2, where neither side gains.
#
# Sample:
#   s <- source("bench/scenarios/14-knn.R")$value
#   d <- s[[1]]$setup(2000); head(s[[1]]$rsimd(d))
suppressPackageStartupMessages(library(rsimd))
if (!exists("scenario", mode = "function")) scenario <- function(...) list(...)

knn_setup <- function(n) {
  X <- matrix(rnorm(64 * n), 64)
  list(X = X, q = rnorm(64), nx = colSums(X^2))
}

invisible(list(
  scenario(
    id = "knn_columns", family = "distance", control = TRUE,
    title = "kNN distances by a loop of simd_dist() over columns (control)",
    why = "Shows the missing vector-against-matrix-columns form: slicing X[, j] copies.",
    sizes = c(2e3, 2e4),
    setup = knn_setup,
    base = function(d) sqrt(colSums((d$X - d$q)^2)),
    rsimd = function(d) vapply(seq_len(ncol(d$X)), function(j) simd_dist(d$X[, j], d$q), 0),
    tolerance = 1e-12,
    gap = list(id = "G1", ideal = function(d) simd_sum(d$X))
  ),
  scenario(
    id = "knn_blas", family = "distance", control = TRUE,
    title = "kNN distances by the BLAS identity (control)",
    why = "crossprod() does the work on both sides; rsimd only does the last elementwise steps.",
    sizes = c(2e3, 2e4),
    setup = knn_setup,
    base = function(d) {
      sqrt(pmax(d$nx - 2 * drop(crossprod(d$X, d$q)) + sum(d$q^2), 0))
    },
    rsimd = function(d) {
      xq <- drop(crossprod(d$X, d$q))
      simd_sqrt(simd_pmax(simd_add(simd_mul_add(xq, -2, d$nx), sum(d$q^2)), 0))
    },
    tolerance = 1e-12
  )
))
