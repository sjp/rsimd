# Elementwise arithmetic. Each function calls one of four entry points with
# the op's name, and the C side applies the broadcast rule, picks the kernel
# and warns. Operands with a class, complex operands and, for the _wrap ops,
# double operands first go through .ew_args(), which checks the types the
# C side cannot name the function for and converts integer64 operands
# mixed with doubles to double, with a warning. Each function makes its own
# .Call(), so that warnings and errors from the C side name the user's call.

# Functions that take complex operands. A double, integer or logical
# operand of a binary one is converted to complex when the other operand
# is complex, as in base R.
.ew_complex <- c("simd_add", "simd_sub", "simd_mul", "simd_div", "simd_neg")

# Functions that take integer64 operands.
.ew_i64 <- c(
  "simd_add", "simd_sub", "simd_mul", "simd_div", "simd_add_wrap", "simd_sub_wrap",
  "simd_mul_wrap", "simd_neg", "simd_abs", "simd_neg_wrap", "simd_abs_wrap", "simd_idiv",
  "simd_mod", "simd_pmin", "simd_pmax", "simd_pmin_num", "simd_pmax_num", "simd_clamp",
  "simd_sign", "simd_mul_add", "simd_add_mul", "simd_mul_add_approx"
)

# Checks every operand in `args` (a named list) for function `fun`; the
# _wrap ops (`wrap = TRUE`) also reject doubles.
.ew_check <- function(fun, args, wrap = FALSE) {
  if (any(fun == .ew_i64) && all(vapply(args, .is_i64_data, NA))) {
    return(invisible())
  }
  cplx <- fun %in% .ew_complex
  unsupported <- c(if (!(fun %in% .ew_i64)) "integer64", if (!cplx) "complex", if (wrap) "double")
  if (cplx && length(args) == 2L) .check_complex_i64(fun, args[[1L]], args[[2L]])
  for (arg in names(args)) .check_supported(args[[arg]], fun, unsupported, character(), arg)
  invisible()
}

# Whether x is an integer64, plain or as a simd_vec (no other class), which
# needs no further checks for a function that takes integer64.
.is_i64_data <- function(x) {
  cls <- oldClass(x)
  n <- length(cls)
  n > 0L && n <= 2L && cls[[n]] == "integer64" && (n == 1L || cls[[1L]] == "simd_vec")
}

# Complex operands do not combine with integer64 ones.
.check_complex_i64 <- function(fun, x, y) {
  i64 <- inherits(x, "integer64") || inherits(y, "integer64")
  if (i64 && (is.complex(x) || is.complex(y))) {
    .stop(fun, "() cannot combine complex and integer64 operands")
  }
  invisible()
}

# The operands in `args` (a named list) of elementwise function `fun`,
# returned as they are when none has a class other than a simd_vec or
# integer64 (.is_i64_data()) and they are all complex (for a function that
# takes complex) or none is (nor, for the _wrap ops, double; with an
# integer64 operand, for a function that takes integer64 other than
# simd_div, nor double); otherwise
# checked (.ew_check()) and converted for the C side: with a complex
# operand the others become complex, otherwise integer64 operands become
# double when another operand is a double (always for simd_div, whose
# result is double). The integer64 warning names `call`, the user's call.
.ew_args <- function(fun, args, wrap = FALSE, call = sys.call(-1L)) {
  plain <- TRUE
  ncplx <- 0L
  i64 <- dbl <- raw <- FALSE
  for (a in args) {
    if (is.object(a)) {
      # .is_i64_data() inline: this loop is the per-call cost.
      cls <- oldClass(a)
      n <- length(cls)
      if (cls[[n]] == "integer64" && (n == 1L || (n == 2L && cls[[1L]] == "simd_vec"))) {
        i64 <- TRUE
        next
      }
      if (n != 1L || cls != "simd_vec") plain <- FALSE
    }
    if (is.complex(a)) ncplx <- ncplx + 1L
    if (is.double(a)) dbl <- TRUE
    if (is.raw(a)) raw <- TRUE
  }
  if (plain && (if (ncplx) {
    ncplx == length(args) && !i64 && fun %in% .ew_complex
  } else if (i64) {
    !dbl && fun != "simd_div" && any(fun == .ew_i64)
  } else {
    !(wrap && dbl)
  })) {
    return(args)
  }
  .ew_check(fun, args, wrap)
  # A raw operand is left for the C side, which gives base R's message.
  if (length(args) == 2L && ncplx && !raw) {
    p <- .promote_pair(args[[1L]], args[[2L]], call)
    return(list(p$x, p$y))
  }
  .i64_to_double(args, call, always = fun == "simd_div")
}

# na.rm of pmin/pmax: TRUE or FALSE.
.arg_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .stop("'", name, "' must be TRUE or FALSE")
  }
  x
}

simd_add <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_add", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "add", na_check)
}

simd_sub <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_sub", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "sub", na_check)
}

simd_mul <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_mul", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "mul", na_check)
}

simd_div <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_div", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "div", na_check)
}

simd_idiv <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_idiv", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "idiv", na_check)
}

simd_mod <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_mod", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "mod", na_check)
}

simd_add_wrap <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y) ||
    is.double(x) || is.double(y)) {
    p <- .ew_args("simd_add_wrap", list(x = x, y = y), wrap = TRUE)
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "add_wrap", na_check)
}

simd_sub_wrap <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y) ||
    is.double(x) || is.double(y)) {
    p <- .ew_args("simd_sub_wrap", list(x = x, y = y), wrap = TRUE)
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "sub_wrap", na_check)
}

simd_mul_wrap <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y) ||
    is.double(x) || is.double(y)) {
    p <- .ew_args("simd_mul_wrap", list(x = x, y = y), wrap = TRUE)
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "mul_wrap", na_check)
}

simd_neg <- function(x) {
  if (is.object(x) || is.complex(x)) .ew_check("simd_neg", list(x = x))
  .Call(C_simd_ew1, x, "neg")
}

# abs of complex z is its modulus, as in base R.
simd_abs <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath, x, "mod", NULL))
  }
  if (is.object(x)) .ew_check("simd_abs", list(x = x))
  .Call(C_simd_ew1, x, "abs")
}

simd_neg_wrap <- function(x) {
  if (is.object(x) || is.complex(x) || is.double(x)) .ew_check("simd_neg_wrap", list(x = x), TRUE)
  .Call(C_simd_ew1, x, "neg")
}

simd_abs_wrap <- function(x) {
  if (is.object(x) || is.complex(x) || is.double(x)) .ew_check("simd_abs_wrap", list(x = x), TRUE)
  .Call(C_simd_ew1, x, "abs")
}

simd_sign <- function(x) {
  if (is.object(x) || is.complex(x)) .ew_check("simd_sign", list(x = x))
  .Call(C_simd_ew1, x, "sign")
}

simd_copysign <- function(x, sign, na_check = NULL) {
  if (is.object(x) || is.object(sign) || is.complex(x) || is.complex(sign)) {
    p <- .ew_args("simd_copysign", list(x = x, sign = sign))
    x <- p[[1L]]
    sign <- p[[2L]]
  }
  .Call(C_simd_ew2, x, sign, "copysign", na_check)
}

simd_recip <- function(x) {
  if (is.object(x) || is.complex(x)) .ew_check("simd_recip", list(x = x))
  .Call(C_simd_ew1, x, "recip")
}

simd_sqrt <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "sqrt", NULL))
  }
  if (is.object(x)) .ew_check("simd_sqrt", list(x = x))
  .Call(C_simd_ew1, x, "sqrt")
}

simd_fma <- function(x, y, z, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.object(z) ||
    is.complex(x) || is.complex(y) || is.complex(z)) {
    p <- .ew_args("simd_fma", list(x = x, y = y, z = z))
    x <- p[[1L]]
    y <- p[[2L]]
    z <- p[[3L]]
  }
  .Call(C_simd_ew3, x, y, z, "fma", na_check)
}

simd_mul_add <- function(x, y, z, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.object(z) ||
    is.complex(x) || is.complex(y) || is.complex(z)) {
    p <- .ew_args("simd_mul_add", list(x = x, y = y, z = z))
    x <- p[[1L]]
    y <- p[[2L]]
    z <- p[[3L]]
  }
  .Call(C_simd_ew3, x, y, z, "mul_add", na_check)
}

simd_mul_add_approx <- function(x, y, z, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.object(z) ||
    is.complex(x) || is.complex(y) || is.complex(z)) {
    p <- .ew_args("simd_mul_add_approx", list(x = x, y = y, z = z))
    x <- p[[1L]]
    y <- p[[2L]]
    z <- p[[3L]]
  }
  .Call(C_simd_ew3, x, y, z, "mul_add_approx", na_check)
}

simd_add_mul <- function(x, y, z, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.object(z) ||
    is.complex(x) || is.complex(y) || is.complex(z)) {
    p <- .ew_args("simd_add_mul", list(x = x, y = y, z = z))
    x <- p[[1L]]
    y <- p[[2L]]
    z <- p[[3L]]
  }
  .Call(C_simd_ew3, x, y, z, "add_mul", na_check)
}

simd_lerp <- function(x, y, t, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.object(t) ||
    is.complex(x) || is.complex(y) || is.complex(t)) {
    p <- .ew_args("simd_lerp", list(x = x, y = y, t = t))
    x <- p[[1L]]
    y <- p[[2L]]
    t <- p[[3L]]
  }
  .Call(C_simd_ew3, x, y, t, "lerp", na_check)
}

# pmin, pmax and clamp have no NA check to skip: na_check is TRUE.
simd_pmin <- function(x, y, ..., na.rm = FALSE) {
  if (...length()) .dots_error("simd_pmin", "two vectors; nest calls for more", named = "na.rm")
  op <- if (.arg_flag(na.rm, "na.rm")) "pmin_num" else "pmin"
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_pmin", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, op, TRUE)
}

simd_pmax <- function(x, y, ..., na.rm = FALSE) {
  if (...length()) .dots_error("simd_pmax", "two vectors; nest calls for more", named = "na.rm")
  op <- if (.arg_flag(na.rm, "na.rm")) "pmax_num" else "pmax"
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_pmax", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, op, TRUE)
}

simd_pmin_num <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_pmin_num", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "pmin_num", TRUE)
}

simd_pmax_num <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_pmax_num", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, "pmax_num", TRUE)
}

simd_clamp <- function(x, lo, hi) {
  if (is.object(x) || is.object(lo) || is.object(hi) ||
    is.complex(x) || is.complex(lo) || is.complex(hi)) {
    p <- .ew_args("simd_clamp", list(x = x, lo = lo, hi = hi))
    x <- p[[1L]]
    lo <- p[[2L]]
    hi <- p[[3L]]
  }
  .Call(C_simd_ew3, x, lo, hi, "clamp", TRUE)
}

simd_floor <- function(x) {
  if (is.object(x) || is.complex(x)) .ew_check("simd_floor", list(x = x))
  .Call(C_simd_ew1, x, "floor")
}

simd_ceiling <- function(x) {
  if (is.object(x) || is.complex(x)) .ew_check("simd_ceiling", list(x = x))
  .Call(C_simd_ew1, x, "ceiling")
}

simd_trunc <- function(x) {
  if (is.object(x) || is.complex(x)) .ew_check("simd_trunc", list(x = x))
  .Call(C_simd_ew1, x, "trunc")
}

simd_round <- function(x, digits = 0) {
  if (!(is.numeric(digits) || is.logical(digits)) || length(digits) != 1L) {
    .stop("'digits' must be a single number")
  }
  if (is.object(x) || is.complex(x)) .ew_check("simd_round", list(x = x))
  if (!is.na(digits) && digits == 0) {
    return(.Call(C_simd_ew1, x, "round"))
  }
  .Call(C_simd_round_digits, x, digits)
}
