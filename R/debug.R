# Internal drivers of the C access layer (src/rvec.c), used by the tests.
# Not part of the package interface.

# Reads x in chunks: list(path ("contiguous" or "regions"), regions, sum,
# na, no_na_hint, stride, type). Allocates nothing proportional to x.
.debug_regions <- function(x) .Call(C_simd_debug_regions, x)

# A copy of x read through the chunk loop, with the result attribute policy
# applied; no_na: the op guarantees an NA-free result.
.debug_copy <- function(x, no_na = FALSE) .Call(C_simd_debug_copy, x, no_na)

# The tier a call with operands x and y runs its kernels on (the selected
# one, or the simd_vec operands' pin).
.debug_active <- function(x, y = NULL) .Call(C_simd_debug_active, x, y)

# x + y as double through the binary chunk loop and broadcast rule:
# list(value, n, x_scalar, y_scalar, regions).
.debug_bin <- function(x, y) .Call(C_simd_debug_bin, x, y)

# Reduction `op` started from its identity, with the result fields in
# `fields` overridden, finished over n elements of `type`:
# list(init, value).
.debug_finish <- function(op, type, n, fields = list(), precision = 0L, na_rm = FALSE) {
  .Call(C_simd_debug_finish, op, type, n, fields, precision, na_rm)
}

# The options an entry point would derive from these arguments (NULL
# na_check and precision: the options).
.debug_opts <- function(x, na_rm = NULL, na_check = NULL, precision = NULL) {
  .Call(C_simd_debug_opts, x, na_rm, na_check, precision)
}

# Records n deferred warnings, "<i>:" followed by len x's, and returns n:
# at most the first 8 are issued, each cut to 255 bytes.
.debug_warn <- function(n, len) .Call(C_simd_debug_warn, n, len)

# Self-test kernels: the NA, precision and overflow rules run through the
# active implementation.

# A sum-like reduction folded chunk by chunk with term "x" (sum), "sq"
# (sum_sq), "abs" (sum_abs) or "xy" (dot, with y): list(value, count,
# saw_na, saw_nan, any_true, any_false). Precision defaults to the current
# mode (NULL: option rsimd.precision).
.debug_fold <- function(x, y = NULL, term = "x", na_rm = FALSE, na_check = TRUE,
                        precision = NULL) {
  .Call(C_simd_debug_fold, x, y, term, na_rm, na_check, precision)
}

# any(x) ("any"), all(x) ("all") or a full scan finished as any ("scan"),
# with the same list as .debug_fold().
.debug_lgl <- function(x, op = "any", na_rm = FALSE, na_check = TRUE) {
  .Call(C_simd_debug_lgl, x, op, na_rm, na_check)
}

# Elementwise op on x and y (unary ops ignore y), warning once on integer
# overflow.
.debug_arith <- function(op, x, y = x, na_check = TRUE) {
  .Call(C_simd_debug_arith, x, y, op, na_check)
}

# The distance between the doubles a and b in units in the last place, for
# the accuracy tests: 0 for two NAs, two other NaNs or equal values (+0 and
# -0 included), Inf when only one is missing or they are NA and NaN.
.ulp_dist <- function(a, b) .Call(C_simd_ulp_dist, a, b)

# Sets the complex arithmetic variants by code: c(mul_re, mul_im, div,
# cp_re, cp_im), -1 for base R's operator alone (kernel_types.h), and
# returns the previous codes. For testing the kernels of variants other
# than this R build's; restore the returned codes afterwards.
.debug_c128_variants <- function(codes) .Call(C_simd_c128_set_variants, as.integer(codes))
