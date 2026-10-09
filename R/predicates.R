# Predicates (elementwise, any, all) and elementwise comparisons. Each
# function checks the types the C side takes but the function does not
# (complex, and raw where it is not taken) and calls the entry point with
# the op's name. They are made by .pred() and .cmp() (see .wrapper()).

# Predicates. mode: 0 elementwise, 1 any, 2 all. Only na, nan, finite and
# infinite take raw and complex input; the others check for them here.
.pred_check <- function(x, fun) .check_supported(x, fun, c("complex", "raw"), character())

# The predicate for `op` in `mode`, named simd_is_<op>, with _any or _all
# for modes 1 and 2.
.pred <- function(op, mode) {
  fun <- paste0("simd_is_", op, c("", "_any", "_all")[[mode + 1L]])
  .wrapper(alist(x = ), if (op %in% c("na", "nan", "finite", "infinite")) {
    bquote(.Call(C_simd_pred, x, .(op), .(mode)))
  } else {
    bquote({
      if (is.object(x) || is.complex(x) || is.raw(x)) .pred_check(x, .(fun))
      .Call(C_simd_pred, x, .(op), .(mode))
    })
  })
}

simd_is_na <- .pred("na", 0L)
simd_is_na_any <- .pred("na", 1L)
simd_is_na_all <- .pred("na", 2L)

simd_is_nan <- .pred("nan", 0L)
simd_is_nan_any <- .pred("nan", 1L)
simd_is_nan_all <- .pred("nan", 2L)

simd_is_finite <- .pred("finite", 0L)
simd_is_finite_any <- .pred("finite", 1L)
simd_is_finite_all <- .pred("finite", 2L)

simd_is_infinite <- .pred("infinite", 0L)
simd_is_infinite_any <- .pred("infinite", 1L)
simd_is_infinite_all <- .pred("infinite", 2L)

simd_is_negative <- .pred("negative", 0L)
simd_is_negative_any <- .pred("negative", 1L)
simd_is_negative_all <- .pred("negative", 2L)

simd_is_zero <- .pred("zero", 0L)
simd_is_zero_any <- .pred("zero", 1L)
simd_is_zero_all <- .pred("zero", 2L)

simd_is_normal <- .pred("normal", 0L)
simd_is_normal_any <- .pred("normal", 1L)
simd_is_normal_all <- .pred("normal", 2L)

simd_is_subnormal <- .pred("subnormal", 0L)
simd_is_subnormal_any <- .pred("subnormal", 1L)
simd_is_subnormal_all <- .pred("subnormal", 2L)

simd_is_whole <- .pred("whole", 0L)
simd_is_whole_any <- .pred("whole", 1L)
simd_is_whole_all <- .pred("whole", 2L)

simd_is_even <- .pred("even", 0L)
simd_is_even_any <- .pred("even", 1L)
simd_is_even_all <- .pred("even", 2L)

simd_is_odd <- .pred("odd", 0L)
simd_is_odd_any <- .pred("odd", 1L)
simd_is_odd_all <- .pred("odd", 2L)

simd_is_pow2 <- .pred("pow2", 0L)
simd_is_pow2_any <- .pred("pow2", 1L)
simd_is_pow2_all <- .pred("pow2", 2L)
# A complex comparison takes only == and !=; the others give base R's
# error. Operands with a class, complex operands and raw with non-raw ones
# are converted by .cmp_args() first; the integer64 warning names the
# user's call.
.cmp_args <- function(x, y, op, call = .user_call()) {
  if ((is.complex(x) || is.complex(y)) && !(op %in% c("eq", "ne"))) {
    .stop("invalid comparison with complex values")
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

# The comparison for `op`, named simd_<op>.
.cmp <- function(op) {
  .wrapper(alist(x = , y = ), bquote({
    if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y) || is.raw(x) != is.raw(y)) {
      p <- .cmp_args(x, y, .(op))
      x <- p[[1L]]
      y <- p[[2L]]
    }
    .Call(C_simd_cmp, x, y, .(op))
  }))
}

simd_eq <- .cmp("eq")
simd_ne <- .cmp("ne")
simd_lt <- .cmp("lt")
simd_le <- .cmp("le")
simd_gt <- .cmp("gt")
simd_ge <- .cmp("ge")
