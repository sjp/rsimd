# Predicates (elementwise, any, all) and elementwise comparisons. Each
# function checks the types the C side cannot name the function for
# (complex, and raw where it is not taken) and calls the entry point with
# the op's name.

# mode: 0 elementwise, 1 any, 2 all. negative and zero take neither raw
# nor complex input.
.pred <- function(x, op, mode, fun) {
  .sync_impl()
  if (op %in% c("negative", "zero")) .check_supported(x, fun, c("complex", "raw"), character())
  .Call(C_simd_pred, x, op, mode)
}

simd_is_na <- function(x) .pred(x, "na", 0L, "simd_is_na")
simd_is_na_any <- function(x) .pred(x, "na", 1L, "simd_is_na_any")
simd_is_na_all <- function(x) .pred(x, "na", 2L, "simd_is_na_all")

simd_is_nan <- function(x) .pred(x, "nan", 0L, "simd_is_nan")
simd_is_nan_any <- function(x) .pred(x, "nan", 1L, "simd_is_nan_any")
simd_is_nan_all <- function(x) .pred(x, "nan", 2L, "simd_is_nan_all")

simd_is_finite <- function(x) .pred(x, "finite", 0L, "simd_is_finite")
simd_is_finite_any <- function(x) .pred(x, "finite", 1L, "simd_is_finite_any")
simd_is_finite_all <- function(x) .pred(x, "finite", 2L, "simd_is_finite_all")

simd_is_infinite <- function(x) .pred(x, "infinite", 0L, "simd_is_infinite")
simd_is_infinite_any <- function(x) .pred(x, "infinite", 1L, "simd_is_infinite_any")
simd_is_infinite_all <- function(x) .pred(x, "infinite", 2L, "simd_is_infinite_all")

simd_is_negative <- function(x) .pred(x, "negative", 0L, "simd_is_negative")
simd_is_negative_any <- function(x) .pred(x, "negative", 1L, "simd_is_negative_any")
simd_is_negative_all <- function(x) .pred(x, "negative", 2L, "simd_is_negative_all")

simd_is_zero <- function(x) .pred(x, "zero", 0L, "simd_is_zero")
simd_is_zero_any <- function(x) .pred(x, "zero", 1L, "simd_is_zero_any")
simd_is_zero_all <- function(x) .pred(x, "zero", 2L, "simd_is_zero_all")

# A raw operand compared with a non-raw one is converted as base R does: to
# logical when the other is logical, else to integer. Raw with raw compares
# bytes. An integer64 operand compared with a double is converted to double,
# with a warning.
.cmp <- function(x, y, op, fun) {
  .sync_impl()
  .check_supported(x, fun, "complex", character(0))
  .check_supported(y, fun, "complex", character(0), "y")
  if (is.raw(x) != is.raw(y)) {
    to <- if (is.logical(x) || is.logical(y)) as.logical else as.integer
    if (is.raw(x)) x <- to(x) else y <- to(y)
  }
  p <- .i64_to_double(list(x, y), sys.call(-1L))
  .Call(C_simd_cmp, p[[1L]], p[[2L]], op)
}

simd_eq <- function(x, y) .cmp(x, y, "eq", "simd_eq")
simd_ne <- function(x, y) .cmp(x, y, "ne", "simd_ne")
simd_lt <- function(x, y) .cmp(x, y, "lt", "simd_lt")
simd_le <- function(x, y) .cmp(x, y, "le", "simd_le")
simd_gt <- function(x, y) .cmp(x, y, "gt", "simd_gt")
simd_ge <- function(x, y) .cmp(x, y, "ge", "simd_ge")
