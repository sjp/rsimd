# Internal drivers of the C access layer (src/rvec.c), used by the tests.
# Not part of the package interface.

# Reads x in chunks: list(path ("contiguous" or "regions"), regions, sum,
# na, no_na_hint, stride, type). Allocates nothing proportional to x.
.debug_regions <- function(x) .Call(C_simd_debug_regions, x)

# A copy of x read through the chunk loop, with the result attribute policy
# applied; no_na: the op guarantees an NA-free result.
.debug_copy <- function(x, no_na = FALSE) .Call(C_simd_debug_copy, x, no_na)

# x + y as double through the binary chunk loop and broadcast rule:
# list(value, n, x_scalar, y_scalar, regions).
.debug_bin <- function(x, y) .Call(C_simd_debug_bin, x, y)

# Reduction `op` started from its identity, with the result fields in
# `fields` overridden, finished over n elements of `type`:
# list(init, value).
.debug_finish <- function(op, type, n, fields = list(), precision = 0L) {
  .Call(C_simd_debug_finish, op, type, n, fields, precision)
}

# The options an entry point would derive from these arguments.
.debug_opts <- function(x, na_rm = NULL, na_check = NULL, precision = NULL) {
  .Call(C_simd_debug_opts, x, na_rm, na_check, precision)
}
