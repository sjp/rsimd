# Elementary functions. Each function rejects integer64 input, and complex
# input except for the functions base R defines on complex numbers, which
# go to the complex entry points; the others call one of the entry points
# with the function's name and NULL for the accuracy (the C side reads
# option rsimd.math_accuracy) or 0L for the functions it does not affect.
# The C side reads integer and logical input as doubles, applies the
# broadcast rule and warns. Each function makes its own .Call(), so that
# the warnings name the user's call.

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
.cmath_args <- function(x, y, op, call = sys.call(-1L)) {
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

simd_exp <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "exp", NULL))
  }
  if (is.object(x)) .math_check("simd_exp", list(x = x))
  .Call(C_simd_math1, x, "exp", NULL)
}

simd_exp2 <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_exp2", list(x = x))
  .Call(C_simd_math1, x, "exp2", NULL)
}

simd_exp10 <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_exp10", list(x = x))
  .Call(C_simd_math1, x, "exp10", NULL)
}

simd_expm1 <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_expm1", list(x = x))
  .Call(C_simd_math1, x, "expm1", NULL)
}

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

simd_log1p <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_log1p", list(x = x))
  .Call(C_simd_math1, x, "log1p", NULL)
}

# None of these depends on rsimd.math_accuracy.
simd_exp2m1 <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_exp2m1", list(x = x))
  .Call(C_simd_math1, x, "exp2m1", 0L)
}

simd_exp10m1 <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_exp10m1", list(x = x))
  .Call(C_simd_math1, x, "exp10m1", 0L)
}

simd_log2p1 <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_log2p1", list(x = x))
  .Call(C_simd_math1, x, "log2p1", 0L)
}

simd_log10p1 <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_log10p1", list(x = x))
  .Call(C_simd_math1, x, "log10p1", 0L)
}

simd_pow <- function(x, y) {
  if (is.complex(x) || is.complex(y)) {
    p <- .cmath_args(x, y, "pow")
    return(.Call(C_simd_cmath2, p[[1L]], p[[2L]], "pow", NULL, "^"))
  }
  if (is.object(x) || is.object(y)) .math_check("simd_pow", list(x = x, y = y))
  .Call(C_simd_math2, x, y, "pow", NULL)
}

simd_cbrt <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_cbrt", list(x = x))
  .Call(C_simd_math1, x, "cbrt", NULL)
}

simd_hypot <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    .math_check("simd_hypot", list(x = x, y = y))
  }
  .Call(C_simd_math2, x, y, "hypot", NULL)
}

simd_sin <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "sin", NULL))
  }
  if (is.object(x)) .math_check("simd_sin", list(x = x))
  .Call(C_simd_math1, x, "sin", NULL)
}

simd_cos <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "cos", NULL))
  }
  if (is.object(x)) .math_check("simd_cos", list(x = x))
  .Call(C_simd_math1, x, "cos", NULL)
}

simd_tan <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "tan", NULL))
  }
  if (is.object(x)) .math_check("simd_tan", list(x = x))
  .Call(C_simd_math1, x, "tan", NULL)
}

simd_sincos <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_sincos", list(x = x))
  .Call(C_simd_sincos, x, NULL, FALSE)
}

simd_asin <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "asin", NULL))
  }
  if (is.object(x)) .math_check("simd_asin", list(x = x))
  .Call(C_simd_math1, x, "asin", NULL)
}

simd_acos <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "acos", NULL))
  }
  if (is.object(x)) .math_check("simd_acos", list(x = x))
  .Call(C_simd_math1, x, "acos", NULL)
}

simd_atan <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "atan", NULL))
  }
  if (is.object(x)) .math_check("simd_atan", list(x = x))
  .Call(C_simd_math1, x, "atan", NULL)
}

simd_atan2 <- function(y, x) {
  if (is.complex(y) || is.complex(x)) {
    p <- .cmath_args(y, x, "atan2")
    return(.Call(C_simd_cmath2, p[[1L]], p[[2L]], "atan2", NULL, "atan2"))
  }
  if (is.object(y) || is.object(x)) .math_check("simd_atan2", list(y = y, x = x))
  .Call(C_simd_math2, y, x, "atan2", NULL)
}

simd_sinpi <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_sinpi", list(x = x))
  .Call(C_simd_math1, x, "sinpi", NULL)
}

simd_cospi <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_cospi", list(x = x))
  .Call(C_simd_math1, x, "cospi", NULL)
}

simd_tanpi <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_tanpi", list(x = x))
  .Call(C_simd_math1, x, "tanpi", NULL)
}

simd_sincospi <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_sincospi", list(x = x))
  .Call(C_simd_sincos, x, NULL, TRUE)
}

simd_sinh <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "sinh", NULL))
  }
  if (is.object(x)) .math_check("simd_sinh", list(x = x))
  .Call(C_simd_math1, x, "sinh", NULL)
}

simd_cosh <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "cosh", NULL))
  }
  if (is.object(x)) .math_check("simd_cosh", list(x = x))
  .Call(C_simd_math1, x, "cosh", NULL)
}

simd_tanh <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "tanh", NULL))
  }
  if (is.object(x)) .math_check("simd_tanh", list(x = x))
  .Call(C_simd_math1, x, "tanh", NULL)
}

simd_asinh <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "asinh", NULL))
  }
  if (is.object(x)) .math_check("simd_asinh", list(x = x))
  .Call(C_simd_math1, x, "asinh", NULL)
}

simd_acosh <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "acosh", NULL))
  }
  if (is.object(x)) .math_check("simd_acosh", list(x = x))
  .Call(C_simd_math1, x, "acosh", NULL)
}

simd_atanh <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "atanh", NULL))
  }
  if (is.object(x)) .math_check("simd_atanh", list(x = x))
  .Call(C_simd_math1, x, "atanh", NULL)
}

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

simd_nextafter <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    .math_check("simd_nextafter", list(x = x, y = y))
  }
  .Call(C_simd_math2, x, y, "nextafter", 0L)
}

simd_next_up <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_next_up", list(x = x))
  .Call(C_simd_math1, x, "next_up", 0L)
}

simd_next_down <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_next_down", list(x = x))
  .Call(C_simd_math1, x, "next_down", 0L)
}

simd_remainder <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    .math_check("simd_remainder", list(x = x, y = y))
  }
  .Call(C_simd_math2, x, y, "remainder", 0L)
}

simd_rsqrt <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_rsqrt", list(x = x))
  .Call(C_simd_math1, x, "rsqrt", 0L)
}

simd_recip_approx <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_recip_approx", list(x = x))
  .Call(C_simd_math1, x, "recip_approx", 0L)
}

simd_rsqrt_approx <- function(x) {
  if (is.object(x) || is.complex(x)) .math_check("simd_rsqrt_approx", list(x = x))
  .Call(C_simd_math1, x, "rsqrt_approx", 0L)
}
