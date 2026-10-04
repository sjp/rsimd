# Three-valued logic, bitwise ops, shifts, rotates and bit counts. Each
# function checks the types the C side cannot name the function for and
# calls the entry point with the op's name; shift and rotate counts are
# checked here.

.logic <- function(x, y, op, fun) {
  .sync_impl()
  .check_supported(x, fun, c("integer64", "complex"), character(0))
  if (!is.null(y)) .check_supported(y, fun, c("integer64", "complex"), character(0), "y")
  .Call(C_simd_logic, x, y, op)
}

simd_and <- function(x, y) .logic(x, y, "and", "simd_and")
simd_or <- function(x, y) .logic(x, y, "or", "simd_or")
simd_xor <- function(x, y) .logic(x, y, "xor", "simd_xor")
simd_not <- function(x) .logic(x, NULL, "not", "simd_not")

# The bit ops take integer, logical and raw vectors; integer64 is planned.
.bit <- function(x, y, op, fun, na_check, k = NULL, raw = TRUE) {
  .sync_impl()
  unsupported <- c("integer64", "complex", "double", if (!raw) "raw")
  .check_supported(x, fun, unsupported, "integer64")
  if (!is.null(y)) .check_supported(y, fun, unsupported, "integer64", "y")
  .Call(C_simd_bit, x, y, op, k, na_check)
}

# n as a single number truncated toward zero, or an error.
.count_arg <- function(n) {
  if (!(is.numeric(n) || is.logical(n)) || length(n) != 1L || inherits(n, "integer64")) {
    stop("'n' must be a single number", call. = FALSE)
  }
  trunc(as.double(n))
}

# The count of a shift of x: for integers NA_integer_ (every element NA, as
# in bitwShiftL()) unless it is in 0..31; for raw it must be in 0..8.
.shift_count <- function(n, x) {
  n <- .count_arg(n)
  if (is.raw(x)) {
    if (is.na(n) || n < 0 || n > 8) {
      stop("argument 'n' must be a small integer", call. = FALSE)
    }
    return(as.integer(n))
  }
  if (is.na(n) || n < 0 || n > 31) NA_integer_ else as.integer(n)
}

# The count of a rotate of x, modulo the width (32, or 8 for raw); a
# negative count rotates the other way.
.rotate_count <- function(n, x) {
  n <- .count_arg(n)
  if (!is.finite(n)) stop("'n' must be a finite number", call. = FALSE)
  w <- if (is.raw(x)) 8 else 32
  as.integer(n - w * floor(n / w))
}

simd_bit_and <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, y, "and", "simd_bit_and", na_check)
}

simd_bit_or <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, y, "or", "simd_bit_or", na_check)
}

simd_bit_xor <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, y, "xor", "simd_bit_xor", na_check)
}

simd_bit_not <- function(x, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, NULL, "not", "simd_bit_not", na_check)
}

simd_shl <- function(x, n, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, NULL, "shl", "simd_shl", na_check, .shift_count(n, x))
}

simd_shr <- function(x, n, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, NULL, "shr", "simd_shr", na_check, .shift_count(n, x))
}

simd_sar <- function(x, n, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, NULL, "sar", "simd_sar", na_check, .shift_count(n, x), raw = FALSE)
}

simd_rotl <- function(x, n, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, NULL, "rotl", "simd_rotl", na_check, .rotate_count(n, x))
}

simd_rotr <- function(x, n, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, NULL, "rotr", "simd_rotr", na_check, .rotate_count(n, x))
}

simd_popcount <- function(x, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, NULL, "popcount", "simd_popcount", na_check)
}

simd_lzcnt <- function(x, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, NULL, "lzcnt", "simd_lzcnt", na_check)
}

simd_tzcnt <- function(x, na_check = getOption("rsimd.na_check", TRUE)) {
  .bit(x, NULL, "tzcnt", "simd_tzcnt", na_check)
}

simd_popcount_total <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .sync_impl()
  .check_supported(x, "simd_popcount_total", c("integer64", "complex", "double"), character(0))
  .Call(C_simd_popcount_total, x, na.rm, na_check)
}
