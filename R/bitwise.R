# Three-valued logic, bitwise ops, shifts, rotates and bit counts. Each
# function checks the types the C side cannot name the function for and
# calls the entry point with the op's name; shift and rotate counts are
# checked here.

# The logic ops take neither integer64 nor complex operands.
.logic_check <- function(fun, x, y = NULL) {
  .check_supported(x, fun, c("integer64", "complex"), character(0))
  if (!is.null(y)) .check_supported(y, fun, c("integer64", "complex"), character(0), "y")
  invisible()
}

simd_and <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) .logic_check("simd_and", x, y)
  .Call(C_simd_logic, x, y, "and")
}

simd_or <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) .logic_check("simd_or", x, y)
  .Call(C_simd_logic, x, y, "or")
}

simd_xor <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) .logic_check("simd_xor", x, y)
  .Call(C_simd_logic, x, y, "xor")
}

simd_not <- function(x) {
  if (is.object(x) || is.complex(x)) .logic_check("simd_not", x)
  .Call(C_simd_logic, x, NULL, "not")
}

# The bit ops take integer, logical, integer64 and raw vectors (sar not
# raw); operands with a class or of another type are checked here.
.bit_check <- function(fun, x, y = NULL, raw = TRUE) {
  unsupported <- c("complex", "double", if (!raw) "raw")
  .check_supported(x, fun, unsupported, character(0))
  if (!is.null(y)) .check_supported(y, fun, unsupported, character(0), "y")
  invisible()
}

# n as a single number truncated toward zero, or an error.
.count_arg <- function(n) {
  if (!(is.numeric(n) || is.logical(n)) || length(n) != 1L || inherits(n, "integer64")) {
    .stop("'n' must be a single number")
  }
  trunc(as.double(n))
}

# The count of a shift of x: for integers NA_integer_ (every element NA, as
# in bitwShiftL()) unless it is in 0..31 (0..63 for integer64); for raw it
# must be in 0..8.
.shift_count <- function(n, x) {
  # A plain number in 0..8 (0..31 unless x is raw) is a valid count.
  if (!is.object(n) && (is.integer(n) || is.double(n)) && length(n) == 1L && !is.na(n) &&
    n >= 0 && n < (if (is.raw(x)) 9 else 32)) {
    return(as.integer(n))
  }
  n <- .count_arg(n)
  if (is.raw(x)) {
    if (is.na(n) || n < 0 || n > 8) {
      .stop("argument 'n' must be a small integer")
    }
    return(as.integer(n))
  }
  w <- if (inherits(x, "integer64")) 63 else 31
  if (is.na(n) || n < 0 || n > w) NA_integer_ else as.integer(n)
}

# The count of a rotate of x, modulo the width (32, 64 for integer64 or 8
# for raw); a negative count rotates the other way.
.rotate_count <- function(n, x) {
  # A plain number in 0..7 (0..31 unless x is raw) is its own count.
  if (!is.object(n) && (is.integer(n) || is.double(n)) && length(n) == 1L && !is.na(n) &&
    n >= 0 && n < (if (is.raw(x)) 8 else 32)) {
    return(as.integer(n))
  }
  n <- .count_arg(n)
  if (!is.finite(n)) .stop("'n' must be a finite number")
  w <- if (is.raw(x)) 8 else if (inherits(x, "integer64")) 64 else 32
  as.integer(n - w * floor(n / w))
}

simd_bit_and <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.double(x) || is.double(y) || is.complex(x) ||
    is.complex(y)) {
    .bit_check("simd_bit_and", x, y)
  }
  .Call(C_simd_bit, x, y, "and", NULL, na_check)
}

simd_bit_or <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.double(x) || is.double(y) || is.complex(x) ||
    is.complex(y)) {
    .bit_check("simd_bit_or", x, y)
  }
  .Call(C_simd_bit, x, y, "or", NULL, na_check)
}

simd_bit_xor <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.double(x) || is.double(y) || is.complex(x) ||
    is.complex(y)) {
    .bit_check("simd_bit_xor", x, y)
  }
  .Call(C_simd_bit, x, y, "xor", NULL, na_check)
}

simd_bit_not <- function(x, na_check = NULL) {
  if (is.object(x) || is.double(x) || is.complex(x)) .bit_check("simd_bit_not", x)
  .Call(C_simd_bit, x, NULL, "not", NULL, na_check)
}

simd_shl <- function(x, n, na_check = NULL) {
  if (is.object(x) || is.double(x) || is.complex(x)) .bit_check("simd_shl", x)
  .Call(C_simd_bit, x, NULL, "shl", .shift_count(n, x), na_check)
}

simd_shr <- function(x, n, na_check = NULL) {
  if (is.object(x) || is.double(x) || is.complex(x)) .bit_check("simd_shr", x)
  .Call(C_simd_bit, x, NULL, "shr", .shift_count(n, x), na_check)
}

simd_sar <- function(x, n, na_check = NULL) {
  if (is.object(x) || is.double(x) || is.complex(x) || is.raw(x)) {
    .bit_check("simd_sar", x, raw = FALSE)
  }
  .Call(C_simd_bit, x, NULL, "sar", .shift_count(n, x), na_check)
}

simd_rotl <- function(x, n, na_check = NULL) {
  if (is.object(x) || is.double(x) || is.complex(x)) .bit_check("simd_rotl", x)
  .Call(C_simd_bit, x, NULL, "rotl", .rotate_count(n, x), na_check)
}

simd_rotr <- function(x, n, na_check = NULL) {
  if (is.object(x) || is.double(x) || is.complex(x)) .bit_check("simd_rotr", x)
  .Call(C_simd_bit, x, NULL, "rotr", .rotate_count(n, x), na_check)
}

simd_popcount <- function(x, na_check = NULL) {
  if (is.object(x) || is.double(x) || is.complex(x)) .bit_check("simd_popcount", x)
  .Call(C_simd_bit, x, NULL, "popcount", NULL, na_check)
}

simd_lzcnt <- function(x, na_check = NULL) {
  if (is.object(x) || is.double(x) || is.complex(x)) .bit_check("simd_lzcnt", x)
  .Call(C_simd_bit, x, NULL, "lzcnt", NULL, na_check)
}

simd_tzcnt <- function(x, na_check = NULL) {
  if (is.object(x) || is.double(x) || is.complex(x)) .bit_check("simd_tzcnt", x)
  .Call(C_simd_bit, x, NULL, "tzcnt", NULL, na_check)
}

simd_popcount_total <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_popcount_total")
  if (is.object(x) || is.double(x) || is.complex(x)) {
    .check_supported(x, "simd_popcount_total", c("complex", "double"), character(0))
  }
  .Call(C_simd_popcount_total, x, na.rm, na_check)
}
