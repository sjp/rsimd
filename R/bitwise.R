# Three-valued logic, bitwise ops, shifts, rotates and bit counts. Each
# function checks the types the C side takes but the function does not, and
# calls the entry point with the op's name; shift and rotate counts are
# checked here. The regular ones are made by .logic2(), .bit2(), .bit1() and
# .bit_count() (see .wrapper()).

# The logic ops do not take complex operands.
.logic_check <- function(fun, x, y = NULL) {
  .check_supported(x, fun, c("complex"), character(0))
  if (!is.null(y)) .check_supported(y, fun, c("complex"), character(0), "y")
  invisible()
}

# The two-operand logic op for `op`, named simd_<op>.
.logic2 <- function(op) {
  fun <- paste0("simd_", op)
  .wrapper(alist(x = , y = ), bquote({
    if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) .logic_check(.(fun), x, y)
    .Call(C_simd_logic, x, y, .(op))
  }))
}

simd_and <- .logic2("and")
simd_or <- .logic2("or")
simd_xor <- .logic2("xor")

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

# The two-operand bit op for `op`, named simd_bit_<op>.
.bit2 <- function(op) {
  .wrapper(alist(x = , y = , na_check = NULL), bquote({
    if (is.object(x) || is.object(y) || is.double(x) || is.double(y) || is.complex(x) ||
      is.complex(y)) {
      .bit_check(.(paste0("simd_bit_", op)), x, y)
    }
    .Call(C_simd_bit, x, y, .(op), NULL, na_check)
  }))
}

simd_bit_and <- .bit2("and")
simd_bit_or <- .bit2("or")
simd_bit_xor <- .bit2("xor")

# The one-operand bit op named `fun` for `op`.
.bit1 <- function(fun, op) {
  .wrapper(alist(x = , na_check = NULL), bquote({
    if (is.object(x) || is.double(x) || is.complex(x)) .bit_check(.(fun), x)
    .Call(C_simd_bit, x, NULL, .(op), NULL, na_check)
  }))
}

simd_bit_not <- .bit1("simd_bit_not", "not")

# The shift or rotate for `op`, named simd_<op>, whose count `count`
# (.shift_count or .rotate_count) gives.
.bit_count <- function(op, count) {
  .wrapper(alist(x = , n = , na_check = NULL), bquote({
    if (is.object(x) || is.double(x) || is.complex(x)) .bit_check(.(paste0("simd_", op)), x)
    .Call(C_simd_bit, x, NULL, .(op), .(as.name(count))(n, x), na_check)
  }))
}

simd_shl <- .bit_count("shl", ".shift_count")
simd_shr <- .bit_count("shr", ".shift_count")

simd_sar <- function(x, n, na_check = NULL) {
  if (is.object(x) || is.double(x) || is.complex(x) || is.raw(x)) {
    .bit_check("simd_sar", x, raw = FALSE)
  }
  .Call(C_simd_bit, x, NULL, "sar", .shift_count(n, x), na_check)
}

simd_rotl <- .bit_count("rotl", ".rotate_count")
simd_rotr <- .bit_count("rotr", ".rotate_count")
simd_popcount <- .bit1("simd_popcount", "popcount")
simd_lzcnt <- .bit1("simd_lzcnt", "lzcnt")
simd_tzcnt <- .bit1("simd_tzcnt", "tzcnt")

simd_popcount_total <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_popcount_total")
  if (is.object(x) || is.double(x) || is.complex(x)) {
    .check_supported(x, "simd_popcount_total", c("complex", "double"), character(0))
  }
  .Call(C_simd_popcount_total, x, na.rm, na_check)
}
