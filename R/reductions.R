# Reductions. Each function makes its own .Call(), so that warnings and
# errors from the C side name the user's call; NULL for na_check and the
# precision tell the C side to read options rsimd.na_check and
# rsimd.precision. Types the C side cannot name the function for are
# checked here first, for operands with a class or of a rejected type.

simd_sum <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_sum")
  .Call(C_simd_sum, x, na.rm, na_check, NULL)
}

simd_prod <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_prod")
  if (is.object(x)) .check_supported(x, "simd_prod", "integer64")
  .Call(C_simd_prod, x, na.rm, na_check, NULL)
}

# The operands of prod(x + y) or prod(x - y) (function `fun`), as a list:
# a complex operand makes the other complex, as in base R.
.prod2_args <- function(fun, x, y) {
  .check_supported(x, fun, "integer64")
  .check_supported(y, fun, "integer64", arg = "y")
  if (is.complex(x) != is.complex(y)) {
    if (is.numeric(x) || is.logical(x)) x <- .as_etype(x, "complex")
    if (is.numeric(y) || is.logical(y)) y <- .as_etype(y, "complex")
  }
  list(x, y)
}

# prod(x + y) (op 0) and prod(x - y) (op 1) without the vector of sums.
simd_prod_sums <- function(x, y, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_prod_sums", "two vectors")
  if (is.object(x) || is.object(y) || is.complex(x) != is.complex(y)) {
    p <- .prod2_args("simd_prod_sums", x, y)
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_prod2, x, y, 0L, na.rm, na_check, NULL)
}

simd_prod_diffs <- function(x, y, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_prod_diffs", "two vectors")
  if (is.object(x) || is.object(y) || is.complex(x) != is.complex(y)) {
    p <- .prod2_args("simd_prod_diffs", x, y)
    x <- p[[1L]]
    y <- p[[2L]]
  }
  .Call(C_simd_prod2, x, y, 1L, na.rm, na_check, NULL)
}

simd_mean <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_mean")
  if (is.object(x)) .check_supported(x, "simd_mean", "integer64")
  .Call(C_simd_mean, x, na.rm, na_check, NULL)
}

# min, max and range share one kernel; op is 0, 1 or 2. Complex input is
# left to the C side, which rejects it with base R's message. With absval
# (max_abs and min_abs) the kernel reads abs(x), Mod(x) for complex x (in
# the math accuracy mode), and raw input is rejected here.
simd_min <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_min")
  .Call(C_simd_minmax, x, 0L, na.rm, na_check, FALSE, 0L)
}

simd_max <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_max")
  .Call(C_simd_minmax, x, 1L, na.rm, na_check, FALSE, 0L)
}

simd_range <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_range")
  .Call(C_simd_minmax, x, 2L, na.rm, na_check, FALSE, 0L)
}

simd_max_abs <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_max_abs")
  if (is.object(x) || is.raw(x)) .check_supported(x, "simd_max_abs", "raw", character())
  .Call(C_simd_minmax, x, 1L, na.rm, na_check, TRUE, NULL)
}

simd_min_abs <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_min_abs")
  if (is.object(x) || is.raw(x)) .check_supported(x, "simd_min_abs", "raw", character())
  .Call(C_simd_minmax, x, 0L, na.rm, na_check, TRUE, NULL)
}

simd_which_min <- function(x) {
  if (is.object(x) || is.complex(x)) {
    .check_supported(x, "simd_which_min", "complex", character())
  }
  .Call(C_simd_which, x, FALSE, FALSE, 0L)
}

simd_which_max <- function(x) {
  if (is.object(x) || is.complex(x)) {
    .check_supported(x, "simd_which_max", "complex", character())
  }
  .Call(C_simd_which, x, TRUE, FALSE, 0L)
}

simd_which_min_abs <- function(x) {
  if (is.object(x) || is.raw(x)) .check_supported(x, "simd_which_min_abs", "raw", character())
  .Call(C_simd_which, x, FALSE, TRUE, NULL)
}

simd_which_max_abs <- function(x) {
  if (is.object(x) || is.raw(x)) .check_supported(x, "simd_which_max_abs", "raw", character())
  .Call(C_simd_which, x, TRUE, TRUE, NULL)
}

simd_any <- function(x, ..., na.rm = FALSE) {
  if (...length()) .dots_error("simd_any", named = "na.rm")
  if (is.object(x) || is.complex(x)) .check_supported(x, "simd_any", "complex", character())
  .Call(C_simd_anyall, x, FALSE, na.rm)
}

simd_all <- function(x, ..., na.rm = FALSE) {
  if (...length()) .dots_error("simd_all", named = "na.rm")
  if (is.object(x) || is.complex(x)) .check_supported(x, "simd_all", "complex", character())
  .Call(C_simd_anyall, x, TRUE, na.rm)
}

simd_any_na <- function(x) .Call(C_simd_na, x, 0L)

simd_count_na <- function(x) .Call(C_simd_na, x, 1L)

simd_which_na <- function(x) .Call(C_simd_na, x, 2L)

# which(x) and sum(x) of a logical x; other types are errors.
.not_logical <- c("double", "integer", "raw", "complex", "integer64")

simd_which <- function(x) {
  if (is.object(x) || !is.logical(x)) .check_supported(x, "simd_which", .not_logical, character())
  .Call(C_simd_true, x, 1L, FALSE)
}

simd_count <- function(x, ..., na.rm = FALSE) {
  if (...length()) .dots_error("simd_count", named = "na.rm")
  if (is.object(x) || !is.logical(x)) .check_supported(x, "simd_count", .not_logical, character())
  .Call(C_simd_true, x, 0L, na.rm)
}

# sum_sq, norm and sum_abs share one entry point; op is 0, 1 or 2. They
# take neither integer64 nor complex input.
simd_sum_sq <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_sum_sq")
  if (is.object(x) || is.complex(x)) .check_real(x, "simd_sum_sq")
  .Call(C_simd_sum_sq, x, 0L, na.rm, na_check, NULL)
}

simd_norm <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_norm")
  if (is.object(x) || is.complex(x)) .check_real(x, "simd_norm")
  .Call(C_simd_sum_sq, x, 1L, na.rm, na_check, NULL)
}

simd_sum_abs <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_sum_abs")
  if (is.object(x) || is.complex(x)) .check_real(x, "simd_sum_abs")
  .Call(C_simd_sum_sq, x, 2L, na.rm, na_check, NULL)
}

# Errors if x (argument `arg` of `fun`) is integer64 or complex.
.check_real <- function(x, fun, arg = "x") {
  .check_supported(x, fun, c("integer64", "complex"), later = character(), arg = arg)
}

# dot, dist and cosine share one entry point; op is 0, 1 or 2.
simd_dot <- function(x, y, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_dot", "two vectors")
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    .check_real(x, "simd_dot")
    .check_real(y, "simd_dot", "y")
  }
  .Call(C_simd_dot, x, y, 0L, na.rm, na_check, NULL)
}

simd_dist <- function(x, y, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_dist", "two vectors")
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    .check_real(x, "simd_dist")
    .check_real(y, "simd_dist", "y")
  }
  .Call(C_simd_dot, x, y, 1L, na.rm, na_check, NULL)
}

simd_cosine <- function(x, y, na_check = NULL) {
  if (is.object(x) || is.object(y) || is.complex(x) || is.complex(y)) {
    .check_real(x, "simd_cosine")
    .check_real(y, "simd_cosine", "y")
  }
  .Call(C_simd_dot, x, y, 2L, FALSE, na_check, NULL)
}

simd_var <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_var")
  if (is.object(x) || is.complex(x)) .check_real(x, "simd_var")
  .Call(C_simd_var, x, FALSE, na.rm, na_check, NULL)
}

simd_sd <- function(x, ..., na.rm = FALSE, na_check = NULL) {
  if (...length()) .dots_error("simd_sd")
  if (is.object(x) || is.complex(x)) .check_real(x, "simd_sd")
  .Call(C_simd_var, x, TRUE, na.rm, na_check, NULL)
}

# The scans share one entry point; op is 0 (cumsum), 1 (cumprod), 2
# (cummin) or 3 (cummax). Only cumprod rejects integer64, and only cummin
# and cummax complex.
simd_cumsum <- function(x) {
  if (is.object(x)) .check_data(x)
  .Call(C_simd_scan, x, 0L, NULL)
}

simd_cumprod <- function(x) {
  if (is.object(x)) .check_supported(x, "simd_cumprod", "integer64", later = character())
  .Call(C_simd_scan, x, 1L, NULL)
}

simd_cummin <- function(x) {
  if (is.object(x) || is.complex(x)) {
    .check_supported(x, "simd_cummin", "complex", later = character())
  }
  .Call(C_simd_scan, x, 2L, NULL)
}

simd_cummax <- function(x) {
  if (is.object(x) || is.complex(x)) {
    .check_supported(x, "simd_cummax", "complex", later = character())
  }
  .Call(C_simd_scan, x, 3L, NULL)
}
