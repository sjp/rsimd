# Elementary functions. Each function rejects complex and integer64 input,
# which the C side cannot name the function for, and calls one of the
# entry points with the function's name and the code of option
# rsimd.math_accuracy; the C side reads integer and logical input as
# doubles, applies the broadcast rule and warns.

.math_check <- function(fun, args) {
  for (arg in names(args)) {
    .check_supported(args[[arg]], fun, c("integer64", "complex"), character(), arg)
  }
  invisible()
}

# `accuracy` is 0L for the functions that rsimd.math_accuracy does not
# affect, which then do not read it.
.math1 <- function(x, op, accuracy = .math_accuracy_code()) {
  .sync_impl()
  .math_check(paste0("simd_", op), list(x = x))
  .Call(C_simd_math1, x, op, accuracy)
}

.math2 <- function(x, y, op, names = c("x", "y"), accuracy = .math_accuracy_code()) {
  .sync_impl()
  args <- list(x, y)
  names(args) <- names
  .math_check(paste0("simd_", op), args)
  .Call(C_simd_math2, x, y, op, accuracy)
}

# The whole-number argument n of simd_scaleb() and simd_rootn() as an
# integer vector: integer and logical as they are, doubles only if whole
# and in the integer range (NA and NaN become NA). Attributes are kept, so
# a simd_vec keeps its pin.
.whole_arg <- function(n) {
  if (is.double(n)) {
    bad <- !is.na(n) & (abs(n) > .Machine$integer.max | n != trunc(n))
    if (any(bad)) stop("'n' must be whole numbers in the integer range", call. = FALSE)
  }
  if (!is.integer(n)) storage.mode(n) <- "integer"
  n
}

.math_n <- function(x, n, op) {
  .sync_impl()
  fun <- paste0("simd_", op)
  .math_check(fun, list(x = x, n = n))
  if (!(is.numeric(n) || is.logical(n))) {
    stop("non-numeric argument to mathematical function", call. = FALSE)
  }
  .Call(C_simd_math2, x, .whole_arg(n), op, 0L)
}

simd_exp <- function(x) .math1(x, "exp")

simd_exp2 <- function(x) .math1(x, "exp2")

simd_exp10 <- function(x) .math1(x, "exp10")

simd_expm1 <- function(x) .math1(x, "expm1")

simd_log <- function(x, base = exp(1)) {
  if (missing(base)) {
    return(.math1(x, "log"))
  }
  .sync_impl()
  .math_check("simd_log", list(x = x))
  .Call(C_simd_log, x, base, .math_accuracy_code())
}

simd_log2 <- function(x) .math1(x, "log2")

simd_log10 <- function(x) .math1(x, "log10")

simd_log1p <- function(x) .math1(x, "log1p")

simd_pow <- function(x, y) .math2(x, y, "pow")

simd_cbrt <- function(x) .math1(x, "cbrt")

simd_hypot <- function(x, y) .math2(x, y, "hypot")

simd_sin <- function(x) .math1(x, "sin")

simd_cos <- function(x) .math1(x, "cos")

simd_tan <- function(x) .math1(x, "tan")

simd_sincos <- function(x) {
  .sync_impl()
  .math_check("simd_sincos", list(x = x))
  .Call(C_simd_sincos, x, .math_accuracy_code())
}

simd_asin <- function(x) .math1(x, "asin")

simd_acos <- function(x) .math1(x, "acos")

simd_atan <- function(x) .math1(x, "atan")

simd_atan2 <- function(y, x) .math2(y, x, "atan2", names = c("y", "x"))

simd_sinpi <- function(x) .math1(x, "sinpi")

simd_cospi <- function(x) .math1(x, "cospi")

simd_tanpi <- function(x) .math1(x, "tanpi")

simd_sinh <- function(x) .math1(x, "sinh")

simd_cosh <- function(x) .math1(x, "cosh")

simd_tanh <- function(x) .math1(x, "tanh")

simd_asinh <- function(x) .math1(x, "asinh")

simd_acosh <- function(x) .math1(x, "acosh")

simd_atanh <- function(x) .math1(x, "atanh")

simd_ilogb <- function(x) {
  .sync_impl()
  .math_check("simd_ilogb", list(x = x))
  .Call(C_simd_ilogb, x)
}

simd_scaleb <- function(x, n) .math_n(x, n, "scaleb")

simd_rootn <- function(x, n) .math_n(x, n, "rootn")

simd_nextafter <- function(x, y) .math2(x, y, "nextafter", accuracy = 0L)

simd_next_up <- function(x) .math1(x, "next_up", accuracy = 0L)

simd_next_down <- function(x) .math1(x, "next_down", accuracy = 0L)

simd_remainder <- function(x, y) .math2(x, y, "remainder", accuracy = 0L)

simd_rsqrt <- function(x) .math1(x, "rsqrt", accuracy = 0L)

simd_recip_approx <- function(x) .math1(x, "recip_approx", accuracy = 0L)

simd_rsqrt_approx <- function(x) .math1(x, "rsqrt_approx", accuracy = 0L)
