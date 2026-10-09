# Elementwise arithmetic. Each function calls one of four entry points with
# the op's name, and the C side applies the broadcast rule, picks the kernel
# and warns. Operands with a class, complex operands (for simd_add, _sub,
# _mul and _div only beside a non-complex one: the C side takes two complex
# operands) and, for the _wrap ops, double operands first go through
# .ew_args(), which checks the types the C side takes but the function does
# not, and converts integer64 operands mixed with doubles to double, with a
# warning. Each function makes its own .Call(), so that warnings and errors
# from the C side name the user's call; the regular ones are made by .ew1(),
# .ew2(), .ew3() and .pmin_num() (see .wrapper()).

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
.ew_args <- function(fun, args, wrap = FALSE, call = .user_call()) {
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
  # A raw operand is left for the C side, which rejects it.
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

# The two-operand function for `op`, named simd_<op>. `guard` is the test
# for operands that go through .ew_args(): "mixed" for a complex operand
# beside a non-complex one, "complex" for any complex operand and "wrap" for
# any complex or double operand.
.ew2 <- function(op, guard) {
  fun <- paste0("simd_", op)
  test <- switch(guard,
    mixed = quote(is.object(x) || is.object(y) || is.complex(x) != is.complex(y)),
    complex = quote(is.object(x) || is.object(y) || is.complex(x) || is.complex(y)),
    wrap = quote(is.object(x) || is.object(y) || is.complex(x) || is.complex(y) ||
      is.double(x) || is.double(y))
  )
  args <- if (guard == "wrap") {
    bquote(.ew_args(.(fun), list(x = x, y = y), wrap = TRUE))
  } else {
    bquote(.ew_args(.(fun), list(x = x, y = y)))
  }
  .wrapper(alist(x = , y = , na_check = NULL), bquote({
    if (.(test)) {
      p <- .(args)
      x <- p[[1L]]
      y <- p[[2L]]
    }
    .Call(C_simd_ew2, x, y, .(op), na_check)
  }))
}

simd_add <- .ew2("add", "mixed")
simd_sub <- .ew2("sub", "mixed")
simd_mul <- .ew2("mul", "mixed")
simd_div <- .ew2("div", "mixed")
simd_idiv <- .ew2("idiv", "complex")
simd_mod <- .ew2("mod", "complex")
simd_add_wrap <- .ew2("add_wrap", "wrap")
simd_sub_wrap <- .ew2("sub_wrap", "wrap")
simd_mul_wrap <- .ew2("mul_wrap", "wrap")

# The one-operand function named `fun` for `op`; with `wrap = TRUE` it also
# rejects doubles.
.ew1 <- function(fun, op, wrap = FALSE) {
  .wrapper(alist(x = ), if (wrap) {
    bquote({
      if (is.object(x) || is.complex(x) || is.double(x)) .ew_check(.(fun), list(x = x), TRUE)
      .Call(C_simd_ew1, x, .(op))
    })
  } else {
    bquote({
      if (is.object(x) || is.complex(x)) .ew_check(.(fun), list(x = x))
      .Call(C_simd_ew1, x, .(op))
    })
  })
}

simd_neg <- .ew1("simd_neg", "neg")

# abs of complex z is its modulus, as in base R.
simd_abs <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath, x, "mod", NULL))
  }
  if (is.object(x)) .ew_check("simd_abs", list(x = x))
  .Call(C_simd_ew1, x, "abs")
}

simd_neg_wrap <- .ew1("simd_neg_wrap", "neg", wrap = TRUE)
simd_abs_wrap <- .ew1("simd_abs_wrap", "abs", wrap = TRUE)
simd_sign <- .ew1("simd_sign", "sign")

simd_copysign <- function(x, sign, na_check = NULL) {
  if (is.object(x) || is.object(sign) || is.complex(x) || is.complex(sign)) {
    p <- .ew_args("simd_copysign", list(x = x, sign = sign))
    x <- p[[1L]]
    sign <- p[[2L]]
  }
  .Call(C_simd_ew2, x, sign, "copysign", na_check)
}

simd_recip <- .ew1("simd_recip", "recip")

simd_sqrt <- function(x) {
  if (is.complex(x)) {
    return(.Call(C_simd_cmath1, x, "sqrt", NULL))
  }
  if (is.object(x)) .ew_check("simd_sqrt", list(x = x))
  .Call(C_simd_ew1, x, "sqrt")
}

# The three-operand function for `op`, named simd_<op>.
.ew3 <- function(op) {
  .wrapper(alist(x = , y = , z = , na_check = NULL), bquote({
    if (is.object(x) || is.object(y) || is.object(z) ||
      is.complex(x) || is.complex(y) || is.complex(z)) {
      p <- .ew_args(.(paste0("simd_", op)), list(x = x, y = y, z = z))
      x <- p[[1L]]
      y <- p[[2L]]
      z <- p[[3L]]
    }
    .Call(C_simd_ew3, x, y, z, .(op), na_check)
  }))
}

simd_fma <- .ew3("fma")
simd_mul_add <- .ew3("mul_add")
simd_mul_add_approx <- .ew3("mul_add_approx")
simd_add_mul <- .ew3("add_mul")

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
  # .arg_flag() inline: a closure call costs as much as the rest.
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    .stop("'na.rm' must be TRUE or FALSE")
  }
  op <- if (na.rm) "pmin_num" else "pmin"
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_pmin", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, op, TRUE)
}

simd_pmax <- function(x, y, ..., na.rm = FALSE) {
  if (...length()) .dots_error("simd_pmax", "two vectors; nest calls for more", named = "na.rm")
  # .arg_flag() inline: a closure call costs as much as the rest.
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    .stop("'na.rm' must be TRUE or FALSE")
  }
  op <- if (na.rm) "pmax_num" else "pmax"
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    p <- .ew_args("simd_pmax", list(x = x, y = y))
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, op, TRUE)
}

# simd_pmin_num() and simd_pmax_num() (`op`).
.pmin_num <- function(op) {
  .wrapper(alist(x = , y = ), bquote({
    if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
      p <- .ew_args(.(paste0("simd_", op)), list(x = x, y = y))
      x <- p[[1L]]
      y <- p[[2L]]
    }
    .Call(C_simd_ew2, x, y, .(op), TRUE)
  }))
}

simd_pmin_num <- .pmin_num("pmin_num")
simd_pmax_num <- .pmin_num("pmax_num")

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

simd_floor <- .ew1("simd_floor", "floor")
simd_ceiling <- .ew1("simd_ceiling", "ceiling")
simd_trunc <- .ew1("simd_trunc", "trunc")

simd_round <- function(x, digits = 0) {
  if (!(is.numeric(digits) || is.logical(digits)) || length(digits) != 1L ||
    inherits(digits, "integer64")) {
    .stop("'digits' must be a single number")
  }
  digits <- as.double(unclass(digits))
  if (is.object(x) || is.complex(x)) .ew_check("simd_round", list(x = x))
  if (!is.na(digits) && digits == 0) {
    return(.Call(C_simd_ew1, x, "round"))
  }
  .Call(C_simd_round_digits, x, digits)
}
