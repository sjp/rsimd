# Predicates (elementwise, any, all) and elementwise comparisons. Each
# function checks the types the C side cannot name the function for
# (complex, and raw where it is not taken) and calls the entry point with
# the op's name.

# Predicates. mode: 0 elementwise, 1 any, 2 all. Only na, nan, finite and
# infinite take raw and complex input; the others check for them here.
.pred_check <- function(x, fun) .check_supported(x, fun, c("complex", "raw"), character())

simd_is_na <- function(x) .Call(C_simd_pred, x, "na", 0L)
simd_is_na_any <- function(x) .Call(C_simd_pred, x, "na", 1L)
simd_is_na_all <- function(x) .Call(C_simd_pred, x, "na", 2L)

simd_is_nan <- function(x) .Call(C_simd_pred, x, "nan", 0L)
simd_is_nan_any <- function(x) .Call(C_simd_pred, x, "nan", 1L)
simd_is_nan_all <- function(x) .Call(C_simd_pred, x, "nan", 2L)

simd_is_finite <- function(x) .Call(C_simd_pred, x, "finite", 0L)
simd_is_finite_any <- function(x) .Call(C_simd_pred, x, "finite", 1L)
simd_is_finite_all <- function(x) .Call(C_simd_pred, x, "finite", 2L)

simd_is_infinite <- function(x) .Call(C_simd_pred, x, "infinite", 0L)
simd_is_infinite_any <- function(x) .Call(C_simd_pred, x, "infinite", 1L)
simd_is_infinite_all <- function(x) .Call(C_simd_pred, x, "infinite", 2L)

simd_is_negative <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_negative")
  .Call(C_simd_pred, x, "negative", 0L)
}

simd_is_negative_any <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_negative_any")
  .Call(C_simd_pred, x, "negative", 1L)
}

simd_is_negative_all <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_negative_all")
  .Call(C_simd_pred, x, "negative", 2L)
}

simd_is_zero <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_zero")
  .Call(C_simd_pred, x, "zero", 0L)
}

simd_is_zero_any <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_zero_any")
  .Call(C_simd_pred, x, "zero", 1L)
}

simd_is_zero_all <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_zero_all")
  .Call(C_simd_pred, x, "zero", 2L)
}

simd_is_normal <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_normal")
  .Call(C_simd_pred, x, "normal", 0L)
}

simd_is_normal_any <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_normal_any")
  .Call(C_simd_pred, x, "normal", 1L)
}

simd_is_normal_all <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_normal_all")
  .Call(C_simd_pred, x, "normal", 2L)
}

simd_is_subnormal <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_subnormal")
  .Call(C_simd_pred, x, "subnormal", 0L)
}

simd_is_subnormal_any <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_subnormal_any")
  .Call(C_simd_pred, x, "subnormal", 1L)
}

simd_is_subnormal_all <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_subnormal_all")
  .Call(C_simd_pred, x, "subnormal", 2L)
}

simd_is_whole <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_whole")
  .Call(C_simd_pred, x, "whole", 0L)
}

simd_is_whole_any <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_whole_any")
  .Call(C_simd_pred, x, "whole", 1L)
}

simd_is_whole_all <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_whole_all")
  .Call(C_simd_pred, x, "whole", 2L)
}

simd_is_even <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_even")
  .Call(C_simd_pred, x, "even", 0L)
}

simd_is_even_any <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_even_any")
  .Call(C_simd_pred, x, "even", 1L)
}

simd_is_even_all <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_even_all")
  .Call(C_simd_pred, x, "even", 2L)
}

simd_is_odd <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_odd")
  .Call(C_simd_pred, x, "odd", 0L)
}

simd_is_odd_any <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_odd_any")
  .Call(C_simd_pred, x, "odd", 1L)
}

simd_is_odd_all <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_odd_all")
  .Call(C_simd_pred, x, "odd", 2L)
}

simd_is_pow2 <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_pow2")
  .Call(C_simd_pred, x, "pow2", 0L)
}

simd_is_pow2_any <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_pow2_any")
  .Call(C_simd_pred, x, "pow2", 1L)
}

simd_is_pow2_all <- function(x) {
  if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, "simd_is_pow2_all")
  .Call(C_simd_pred, x, "pow2", 2L)
}

# A complex comparison takes only == and !=; the others give base R's
# error. Operands with a class, complex operands and raw with non-raw ones
# are converted by .cmp_args() first; the integer64 warning names the
# user's call.
.cmp_args <- function(x, y, op, call = sys.call(-1L)) {
  if ((is.complex(x) || is.complex(y)) && !(op %in% c("eq", "ne"))) {
    stop("invalid comparison with complex values", call. = FALSE)
  }
  .cmp_operands(x, y, call)
}

# The operands of a comparison or of simd_hamming(), as a list, converted
# as base R's comparisons do. A raw operand with a non-raw one becomes
# logical when the other is logical, else integer; raw with raw compares
# bytes. An integer64 operand with a double becomes double, with a warning
# naming `call`. A complex operand makes the other complex (integer64
# through double, with the same warning).
.cmp_operands <- function(x, y, call) {
  if (!inherits(x, "integer64") && !inherits(y, "integer64") && is.complex(x) == is.complex(y) &&
    is.raw(x) == is.raw(y)) {
    return(list(x, y))
  }
  if (is.complex(x) || is.complex(y)) {
    return(list(.as_complex_op(x, call), .as_complex_op(y, call)))
  }
  if (is.raw(x) != is.raw(y)) {
    to <- if (is.logical(x) || is.logical(y)) as.logical else as.integer
    if (is.raw(x)) x <- .sv_like(to(x), x) else y <- .sv_like(to(y), y)
  }
  .i64_to_double(list(x, y), call)
}

# An operand of a complex comparison as complex.
.as_complex_op <- function(a, call) {
  if (is.complex(a)) {
    return(a)
  }
  if (inherits(a, "integer64")) {
    warning(simpleWarning("integer64 coerced to double", call))
    a <- .as_etype(a, "double")
  }
  .sv_like(as.complex(a), a)
}

simd_eq <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y) || is.raw(x) != is.raw(y)) {
    p <- .cmp_args(x, y, "eq")
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_cmp, x, y, "eq")
}

simd_ne <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y) || is.raw(x) != is.raw(y)) {
    p <- .cmp_args(x, y, "ne")
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_cmp, x, y, "ne")
}

simd_lt <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y) || is.raw(x) != is.raw(y)) {
    p <- .cmp_args(x, y, "lt")
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_cmp, x, y, "lt")
}

simd_le <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y) || is.raw(x) != is.raw(y)) {
    p <- .cmp_args(x, y, "le")
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_cmp, x, y, "le")
}

simd_gt <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y) || is.raw(x) != is.raw(y)) {
    p <- .cmp_args(x, y, "gt")
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_cmp, x, y, "gt")
}

simd_ge <- function(x, y) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y) || is.raw(x) != is.raw(y)) {
    p <- .cmp_args(x, y, "ge")
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_cmp, x, y, "ge")
}
