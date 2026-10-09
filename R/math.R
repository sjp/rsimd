# Elementary functions. Each function rejects integer64 input, and complex
# input except for the functions base R defines on complex numbers, which
# go to the complex entry points; the others call one of the entry points
# with the function's name and NULL for the accuracy (the C side reads
# option rsimd.math_accuracy) or 0L for the functions it does not affect.
# The C side reads integer and logical input as doubles, applies the
# broadcast rule and warns. Each function makes its own .Call(), so that
# the warnings name the user's call; the regular ones are made by .math1(),
# .cmath1() and .math2() (see .wrapper()).

# Checks the operands in `args` (a named list) of function `fun`.
.math_check <- function(fun, args) {
  for (arg in names(args)) {
    .check_supported(args[[arg]], fun, c("integer64", "complex"), character(), arg)
  }
  invisible()
}

# The operands of x^y ("pow"), log(x, base = y) ("logb") or atan2(x, y)
# ("atan2") when either is complex, as a list: checked, and the other
# converted to complex, as in base R. The integer64 warning names `call`,
# the user's call.
.cmath_args <- function(x, y, op, call = .user_call()) {
  # Unclassed complex, double, integer or logical operands need no checks.
  if (!is.object(x) && !is.object(y) && .cmath_plain(x) && .cmath_plain(y)) {
    if (!is.complex(x)) x <- as.complex(x)
    if (!is.complex(y)) y <- as.complex(y)
    return(list(x, y))
  }
  fun <- paste0("simd_", if (op == "logb") "log" else op)
  .check_complex_i64(fun, x, y)
  args <- list(x, y)
  names(args) <- switch(op, logb = c("x", "base"), atan2 = c("y", "x"), c("x", "y"))
  for (arg in names(args)) {
    a <- args[[arg]]
    if (!(is.numeric(a) || is.logical(a) || is.complex(a))) {
      .check_data(a, arg)
      .stop(fun, "() does not support '", arg, "' of type ", typeof(a))
    }
  }
  p <- .promote_pair(x, y, call)
  list(p$x, p$y)
}

.cmath_plain <- function(x) is.complex(x) || is.double(x) || is.integer(x) || is.logical(x)

# The whole-number argument n of simd_scaleb() and simd_rootn() as an
# integer vector: integer and logical as they are, doubles only if whole
# and in the integer range (NA and NaN become NA). Attributes are kept, so
# a simd_vec keeps its pin.
.whole_arg <- function(n) {
  if (is.double(n)) {
    bad <- !is.na(n) & (abs(n) > .Machine$integer.max | n != trunc(n))
    if (any(bad)) .stop("'n' must be whole numbers in the integer range")
  }
  if (!is.integer(n)) storage.mode(n) <- "integer"
  n
}

# Checks x and n of simd_scaleb() or simd_rootn() (`fun`) and returns n as
# .whole_arg() does.
.n_arg <- function(fun, x, n) {
  .math_check(fun, list(x = x, n = n))
  if (!(is.numeric(n) || is.logical(n))) {
    .check_data(n, "n")
    .stop(fun, "() does not support 'n' of type ", typeof(n))
  }
  .whole_arg(n)
}

# The real-only function for `op`, named simd_<op>. `acc` is the accuracy
# passed to the C side: NULL, or 0L for the functions it does not affect.
.math1 <- function(op, acc) {
  .wrapper(alist(x = ), bquote({
    if (is.object(x) || is.complex(x)) .math_check(.(paste0("simd_", op)), list(x = x))
    .Call(C_simd_math1, x, .(op), .(acc))
  }))
}

# The function for `op`, named simd_<op>, which base R also defines on
# complex numbers.
.cmath1 <- function(op) {
  .wrapper(alist(x = ), bquote({
    if (is.complex(x)) {
      return(.Call(C_simd_cmath1, x, .(op), NULL))
    }
    if (is.object(x)) .math_check(.(paste0("simd_", op)), list(x = x))
    .Call(C_simd_math1, x, .(op), NULL)
  }))
}

# The real-only two-operand function for `op`, named simd_<op>; `acc` as
# for .math1().
.math2 <- function(op, acc) {
  .wrapper(alist(x = , y = ), bquote({
    if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
      .math_check(.(paste0("simd_", op)), list(x = x, y = y))
    }
    .Call(C_simd_math2, x, y, .(op), .(acc))
  }))
}

simd_exp <- .cmath1("exp")
simd_exp2 <- .math1("exp2", NULL)
simd_exp10 <- .math1("exp10", NULL)
simd_expm1 <- .math1("expm1", NULL)

simd_log <- function(x, base = exp(1)) {
  if (missing(base)) {
    if (is.complex(x)) {
      return(.Call(C_simd_cmath1, x, "log", NULL))
    }
    if (is.object(x)) .math_check("simd_log", list(x = x))
    return(.Call(C_simd_math1, x, "log", NULL))
  }
  if (is.complex(x) || is.complex(base)) {
    p <- .cmath_args(x, base, "logb")
    return(.Call(C_simd_cmath2, p[[1L]], p[[2L]], "logb", NULL, "log"))
  }
  if (is.object(x)) .math_check("simd_log", list(x = x))
  .Call(C_simd_log, x, base, NULL)
}

simd_log2 <- function(x) {
  if (is.complex(x)) {
    p <- .cmath_args(x, 2, "logb")
    return(.Call(C_simd_cmath2, p[[1L]], p[[2L]], "logb", NULL, "log2"))
  }
  if (is.object(x)) .math_check("simd_log2", list(x = x))
  .Call(C_simd_math1, x, "log2", NULL)
}

simd_log10 <- function(x) {
  if (is.complex(x)) {
    p <- .cmath_args(x, 10, "logb")
    return(.Call(C_simd_cmath2, p[[1L]], p[[2L]], "logb", NULL, "log10"))
  }
  if (is.object(x)) .math_check("simd_log10", list(x = x))
  .Call(C_simd_math1, x, "log10", NULL)
}

simd_log1p <- .math1("log1p", NULL)

# None of these depends on rsimd.math_accuracy.
simd_exp2m1 <- .math1("exp2m1", 0L)
simd_exp10m1 <- .math1("exp10m1", 0L)
simd_log2p1 <- .math1("log2p1", 0L)
simd_log10p1 <- .math1("log10p1", 0L)

simd_pow <- function(x, y) {
  if (is.complex(x) || is.complex(y)) {
    p <- .cmath_args(x, y, "pow")
    return(.Call(C_simd_cmath2, p[[1L]], p[[2L]], "pow", NULL, "^"))
  }
  if (is.object(x) || is.object(y)) .math_check("simd_pow", list(x = x, y = y))
  .Call(C_simd_math2, x, y, "pow", NULL)
}

simd_cbrt <- .math1("cbrt", NULL)
simd_hypot <- .math2("hypot", NULL)
simd_sin <- .cmath1("sin")
simd_cos <- .cmath1("cos")
simd_tan <- .cmath1("tan")

simd_sincos <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_sincos", list(x = x))
  .Call(C_simd_sincos, x, NULL, FALSE)
}

simd_asin <- .cmath1("asin")
simd_acos <- .cmath1("acos")
simd_atan <- .cmath1("atan")

simd_atan2 <- function(y, x) {
  if (is.complex(y) || is.complex(x)) {
    p <- .cmath_args(y, x, "atan2")
    return(.Call(C_simd_cmath2, p[[1L]], p[[2L]], "atan2", NULL, "atan2"))
  }
  if (is.object(y) || is.object(x)) .math_check("simd_atan2", list(y = y, x = x))
  .Call(C_simd_math2, y, x, "atan2", NULL)
}

simd_sinpi <- .math1("sinpi", NULL)
simd_cospi <- .math1("cospi", NULL)
simd_tanpi <- .math1("tanpi", NULL)

simd_sincospi <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_sincospi", list(x = x))
  .Call(C_simd_sincos, x, NULL, TRUE)
}

simd_sinh <- .cmath1("sinh")
simd_cosh <- .cmath1("cosh")
simd_tanh <- .cmath1("tanh")
simd_asinh <- .cmath1("asinh")
simd_acosh <- .cmath1("acosh")
simd_atanh <- .cmath1("atanh")

simd_ilogb <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_ilogb", list(x = x))
  .Call(C_simd_ilogb, x)
}

simd_scaleb <- function(x, n) {
  .Call(C_simd_math2, x, .n_arg("simd_scaleb", x, n), "scaleb", 0L)
}

simd_rootn <- function(x, n) {
  .Call(C_simd_math2, x, .n_arg("simd_rootn", x, n), "rootn", 0L)
}

simd_nextafter <- .math2("nextafter", 0L)
simd_next_up <- .math1("next_up", 0L)
simd_next_down <- .math1("next_down", 0L)
simd_remainder <- .math2("remainder", 0L)
simd_rsqrt <- .math1("rsqrt", 0L)
simd_recip_approx <- .math1("recip_approx", 0L)
simd_rsqrt_approx <- .math1("rsqrt_approx", 0L)
