# Elementwise arithmetic. Each function checks the types the C side cannot
# name the function for (complex, integer64, double for the _wrap ops) and
# calls one of four entry points with the op's name; the C side applies the
# broadcast rule, picks the kernel and warns. integer64 operands mixed with
# doubles are converted to double here, with a warning.

# Functions that take complex operands. A double, integer or logical
# operand of a binary one is converted to complex when the other operand
# is complex, as in base R.
.ew_complex <- c("simd_add", "simd_sub", "simd_neg")

# Functions that take integer64 operands.
.ew_i64 <- c(
  "simd_add", "simd_sub", "simd_mul", "simd_div", "simd_add_wrap", "simd_sub_wrap",
  "simd_mul_wrap", "simd_neg", "simd_abs", "simd_neg_wrap", "simd_abs_wrap", "simd_idiv",
  "simd_mod", "simd_pmin", "simd_pmax", "simd_pmin_num", "simd_pmax_num", "simd_clamp",
  "simd_sign", "simd_mul_add", "simd_add_mul"
)

# Checks every operand in `args` (a named list) for function `fun`; the
# _wrap ops (`wrap = TRUE`) also reject doubles.
.ew_check <- function(fun, args, wrap = FALSE) {
  cplx <- fun %in% .ew_complex
  unsupported <- c(if (!(fun %in% .ew_i64)) "integer64", if (!cplx) "complex", if (wrap) "double")
  if (cplx && length(args) == 2L) .check_complex_i64(fun, args[[1L]], args[[2L]])
  for (arg in names(args)) .check_supported(args[[arg]], fun, unsupported, character(), arg)
  invisible()
}

# Complex operands do not combine with integer64 ones.
.check_complex_i64 <- function(fun, x, y) {
  i64 <- inherits(x, "integer64") || inherits(y, "integer64")
  if (i64 && (is.complex(x) || is.complex(y))) {
    stop(fun, "() cannot combine complex and integer64 operands", call. = FALSE)
  }
  invisible()
}

# TRUE for an operand none of the checks below can reject or convert: no
# class (so not integer64), not complex and, for the _wrap ops, not double.
# Operands that are all plain go straight to the C side, as the checks cost
# several microseconds per call.
.ew_plain <- function(x, wrap) {
  !is.object(x) && !is.complex(x) && !(wrap && is.double(x))
}

.ew1 <- function(x, op, fun, wrap = FALSE) {
  .sync_impl()
  if (.ew_plain(x, wrap)) {
    return(.Call(C_simd_ew1, x, op))
  }
  .ew_check(fun, list(x = x), wrap)
  .Call(C_simd_ew1, x, op)
}

.ew2 <- function(x, y, op, fun, na_check, wrap = FALSE, names = c("x", "y")) {
  .sync_impl()
  if (.ew_plain(x, wrap) && .ew_plain(y, wrap)) {
    return(.Call(C_simd_ew2, x, y, op, na_check))
  }
  args <- list(x, y)
  names(args) <- names
  .ew_check(fun, args, wrap)
  # A raw operand is left for the C side, which gives base R's message.
  if ((is.complex(x) || is.complex(y)) && !is.raw(x) && !is.raw(y)) {
    p <- .promote_pair(x, y, sys.call(-1L))
    x <- p$x
    y <- p$y
  } else {
    p <- .i64_to_double(list(x, y), sys.call(-1L), always = op == "div")
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_ew2, x, y, op, na_check)
}

.ew3 <- function(x, y, z, op, fun, na_check, names = c("x", "y", "z")) {
  .sync_impl()
  if (.ew_plain(x, FALSE) && .ew_plain(y, FALSE) && .ew_plain(z, FALSE)) {
    return(.Call(C_simd_ew3, x, y, z, op, na_check))
  }
  args <- list(x, y, z)
  names(args) <- names
  .ew_check(fun, args)
  p <- .i64_to_double(list(x, y, z), sys.call(-1L))
  .Call(C_simd_ew3, p[[1L]], p[[2L]], p[[3L]], op, na_check)
}

# na.rm of pmin/pmax: TRUE or FALSE.
.arg_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be TRUE or FALSE", call. = FALSE)
  }
  x
}

simd_add <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "add", "simd_add", na_check)
}

simd_sub <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "sub", "simd_sub", na_check)
}

simd_mul <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "mul", "simd_mul", na_check)
}

simd_div <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "div", "simd_div", na_check)
}

simd_idiv <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "idiv", "simd_idiv", na_check)
}

simd_mod <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "mod", "simd_mod", na_check)
}

simd_add_wrap <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "add_wrap", "simd_add_wrap", na_check, wrap = TRUE)
}

simd_sub_wrap <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "sub_wrap", "simd_sub_wrap", na_check, wrap = TRUE)
}

simd_mul_wrap <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, y, "mul_wrap", "simd_mul_wrap", na_check, wrap = TRUE)
}

simd_neg <- function(x) .ew1(x, "neg", "simd_neg")

simd_abs <- function(x) .ew1(x, "abs", "simd_abs")

simd_neg_wrap <- function(x) .ew1(x, "neg", "simd_neg_wrap", wrap = TRUE)

simd_abs_wrap <- function(x) .ew1(x, "abs", "simd_abs_wrap", wrap = TRUE)

simd_sign <- function(x) .ew1(x, "sign", "simd_sign")

simd_copysign <- function(x, sign, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew2(x, sign, "copysign", "simd_copysign", na_check, names = c("x", "sign"))
}

simd_recip <- function(x) .ew1(x, "recip", "simd_recip")

simd_sqrt <- function(x) .ew1(x, "sqrt", "simd_sqrt")

simd_fma <- function(x, y, z, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew3(x, y, z, "fma", "simd_fma", na_check)
}

simd_mul_add <- function(x, y, z, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew3(x, y, z, "mul_add", "simd_mul_add", na_check)
}

simd_add_mul <- function(x, y, z, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew3(x, y, z, "add_mul", "simd_add_mul", na_check)
}

simd_lerp <- function(x, y, t, na_check = getOption("rsimd.na_check", TRUE)) {
  .ew3(x, y, t, "lerp", "simd_lerp", na_check, names = c("x", "y", "t"))
}

simd_pmin <- function(x, y, na.rm = FALSE) {
  op <- if (.arg_flag(na.rm, "na.rm")) "pmin_num" else "pmin"
  .ew2(x, y, op, "simd_pmin", NULL)
}

simd_pmax <- function(x, y, na.rm = FALSE) {
  op <- if (.arg_flag(na.rm, "na.rm")) "pmax_num" else "pmax"
  .ew2(x, y, op, "simd_pmax", NULL)
}

simd_pmin_num <- function(x, y) .ew2(x, y, "pmin_num", "simd_pmin_num", NULL)

simd_pmax_num <- function(x, y) .ew2(x, y, "pmax_num", "simd_pmax_num", NULL)

simd_clamp <- function(x, lo, hi) {
  .ew3(x, lo, hi, "clamp", "simd_clamp", NULL, names = c("x", "lo", "hi"))
}

simd_floor <- function(x) .ew1(x, "floor", "simd_floor")

simd_ceiling <- function(x) .ew1(x, "ceiling", "simd_ceiling")

simd_trunc <- function(x) .ew1(x, "trunc", "simd_trunc")

simd_round <- function(x, digits = 0) {
  if (!(is.numeric(digits) || is.logical(digits)) || length(digits) != 1L) {
    stop("'digits' must be a single number", call. = FALSE)
  }
  if (!is.na(digits) && digits == 0) {
    return(.ew1(x, "round", "simd_round"))
  }
  .sync_impl()
  .ew_check("simd_round", list(x = x))
  .Call(C_simd_round_digits, x, digits)
}
