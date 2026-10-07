# Helpers for the NA, precision and overflow tests, which drive the internal
# self-test kernels (.debug_fold, .debug_lgl, .debug_arith).

# Unit in the last place of a finite, non-zero double.
ulp <- function(x) 2^(floor(log2(abs(x))) - 52)

# Expects f() to give identical results on every available tier, comparing
# each with the none tier.
expect_tiers_identical <- function(f, label = "") {
  res <- for_each_tier(function(tier) f())
  for (tier in setdiff(names(res), "none")) {
    check_identical(res[[tier]], res[["none"]], info = paste(label, "tier", tier))
  }
  invisible(res)
}

# The value of a self-test fold.
fold <- function(x, y = NULL, term = "x", na_rm = FALSE, na_check = TRUE, precision = 0L) {
  .debug_fold(x, y, term, na_rm, na_check, precision)$value
}

# Precision codes, as the C side numbers them.
prec_codes <- c(fast = 0L, pairwise = 1L, compensated = 2L)

# Has base::sum() an accumulator wider than double? Not under a noLD build,
# nor where long double is double (arm64 macOS), nor under valgrind, which
# computes x87 long double in double precision (R measures the digits at
# startup, so they show 53 there although sizeof is 16).
has_wide_long_double <- function() {
  isTRUE(capabilities("long.double")) && .Machine$sizeof.longdouble > 8 &&
    isTRUE(.Machine$longdouble.digits >= 64)
}

# 257 doubles with 1e16, 1, -1e16 at positions 1, 129 and 257: in fast mode
# all three land in accumulator 0, lane 0 for every lane width up to 32, so
# the fast sum is 0 on every tier.
tier_proof_cancel <- function() {
  x <- numeric(257)
  x[c(1, 129, 257)] <- c(1e16, 1, -1e16)
  x
}
