# Generators of edge-case inputs for kernel tests.

# Lanes per vector over all tiers: 2 (128-bit, 64-bit lanes) up to 32
# (2048-bit SVE); 32-bit lanes and the 64-bit lanes of integer sums fall in
# the same range.
lane_widths <- c(2, 4, 8, 16, 32)

# Lengths around the vector width W: empty, single, the tail cases, the
# fast-mode block of four vectors (4 * W * 4 + 5 crosses it), the ALTREP
# region size (4096) and a length above the interrupt stride. For several
# W, the sorted union. Reduced runs (helper-subset.R) keep the lengths up
# to 1000, and 1000 itself.
edge_lengths <- function(W = lane_widths) {
  one <- function(W) c(0, 1, W - 1, W, W + 1, 2 * W + 3, 4 * W * 4 + 5, 4095, 4096, 4097, 1e6)
  n <- sort(unique(unlist(lapply(W, one))))
  if (reduced_lengths()) c(n[n < 1000], 1000) else n
}

# Special doubles. Sums that mix double.xmax with other large values
# overflow or not depending on the order of additions, which differs
# between tiers, so oracle sweeps use xmax = FALSE.
edge_doubles <- function(xmax = TRUE) {
  x <- c(
    0, -0, 1, -1, .Machine$double.xmin, 5e-324, .Machine$double.xmax,
    Inf, -Inf, NA, NaN, 1e16, -1e16, pi
  )
  if (xmax) x else x[x != .Machine$double.xmax | is.na(x)]
}
edge_ints <- function() c(0L, 1L, -1L, .Machine$integer.max, -.Machine$integer.max, NA)
edge_lgl <- function() c(TRUE, FALSE, NA)

# Compact (ALTREP) sequences of length n >= 2, which the access layer reads
# in regions without expanding them. (Length 1 and rev() results are
# ordinary vectors.)
altrep_inputs <- function(n) {
  stopifnot(n >= 2)
  list(colon = 1:n, seq_len = seq_len(n), double = as.double(1:n))
}

# Is x read through the ALTREP region path?
takes_region_path <- function(x) identical(.debug_regions(x)$path, "regions")

# Evaluates `code` with the random seed set to `seed`, restoring the
# previous state of the generator afterwards.
with_seed <- function(seed, code) {
  genv <- globalenv()
  old <- if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
    get(".Random.seed", envir = genv, inherits = FALSE)
  }
  on.exit(
    if (is.null(old)) {
      rm(".Random.seed", envir = genv)
    } else {
      assign(".Random.seed", old, envir = genv)
    }
  )
  set.seed(seed)
  code
}

# A random vector of `type` ("double", "integer" or "logical") and length
# n from seed `seed`, with about na_frac NA; doubles also get about
# nan_frac NaN and inf_frac +-Inf. Doubles span several magnitudes,
# integers [-1e6, 1e6].
rand_vec <- function(type, n, na_frac = 0.05, nan_frac = 0.02, inf_frac = 0.01, seed = 1L) {
  with_seed(seed, {
    x <- switch(type,
      double = stats::rnorm(n) * 10^stats::runif(n, -3, 3),
      integer = sample.int(2000001L, n, replace = TRUE) - 1000001L,
      logical = stats::runif(n) < 0.5,
      stop("unknown type ", type)
    )
    put <- function(x, frac, value) {
      k <- stats::rbinom(1, n, frac)
      x[sample.int(n, k)] <- value
      x
    }
    if (type == "double") {
      x <- put(x, inf_frac, Inf)
      x <- put(x, inf_frac / 2, -Inf)
      x <- put(x, nan_frac, NaN)
    }
    put(x, na_frac, NA)
  })
}
