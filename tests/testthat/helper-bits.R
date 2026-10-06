# Helpers for the predicate, bitwise and conversion tests.

# Doubles with the given bit patterns (16 hex digits each, most significant
# first), e.g. "7ff00000000007a2" for NA_real_.
from_bits <- function(hex) {
  vapply(hex, function(h) {
    bytes <- strtoi(substring(h, seq(15, 1, -2), seq(16, 2, -2)), 16L)
    readBin(as.raw(bytes), "double", size = 8, endian = "little")
  }, double(1), USE.NAMES = FALSE)
}

# NaNs with assorted payloads and signs; the low word 1954 (0x7a2) makes NA.
nan_payloads <- function() {
  from_bits(c(
    "7ff8000000000000", "fff8000000000000", "7ff00000000007a2", "7ff80000000007a2",
    "fff80000000007a2", "fff00000000007a2", "7ff0000000000001", "7ff00000000007a3",
    "7ffc0000000007a2", "7fffffffffffffff", "7ff40000deadbeef", "7ff00001000007a2"
  ))
}

# Doubles to test predicates and conversions on: signed zeros, denormals,
# infinities, every NaN payload above and numbers around the integer and
# raw limits.
pred_doubles <- function() {
  c(
    0, -0, 1, -1, 0.5, -0.5, 5e-324, -5e-324, .Machine$double.xmin, .Machine$double.xmax,
    -.Machine$double.xmax, Inf, -Inf, NA, NaN, 2^31, -2^31, 2^53, 255, 256, pi,
    nan_payloads()
  )
}

# x as unsigned 32-bit values in doubles (NA stays NA).
u32 <- function(x) as.double(x) %% 2^32

# Reference bit counts of integers (NA for NA). popcount_ref() adds the
# counts of the four bytes, looked up in byte_pop; the others work from the
# unsigned values in doubles.
byte_pop <- vapply(0:255, function(b) sum(bitwAnd(b, 2L^(0:7)) != 0L), integer(1))
popcount_ref <- function(x) {
  stopifnot(is.integer(x))
  r <- integer(length(x))
  for (k in 0:3) r <- r + byte_pop[bitwAnd(bitwShiftR(x, 8L * k), 255L) + 1L]
  r
}
lzcnt_ref <- function(x) {
  u <- u32(x)
  r <- ifelse(u == 0, 32, 31 - floor(log2(pmax(u, 1))))
  as.integer(r)
}
tzcnt_ref <- function(x) {
  u <- u32(x)
  r <- numeric(length(u))
  for (k in 1:32) r <- r + (u %% 2^k == 0)
  as.integer(r)
}

# Unsigned 32-bit values (doubles in [0, 2^32)) back to integers, with
# 0x80000000 reading as NA.
from_u32 <- function(u) {
  u[!is.na(u) & u >= 2^31] <- u[!is.na(u) & u >= 2^31] - 2^32
  suppressWarnings(as.integer(ifelse(u == -2^31, NA, u)))
}

# The value and the warning messages of evaluating `expr`.
value_and_warnings <- function(expr) {
  value <- NULL
  w <- warnings_of(value <- expr)
  list(value = value, warnings = w)
}

# Expects f(...) on every tier to give base R's value and warnings for the
# same call of base_f (identical value, same messages in the same order).
expect_tiers_like_base <- function(base_f, f, ...) {
  want <- value_and_warnings(base_f(...))
  got <- expect_tier_warnings(want$warnings, f, ...)
  check_identical(got, want$value)
  expect_tiers_give(want$value, function(...) suppressWarnings(f(...)), ...)
}
