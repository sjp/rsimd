simd_sum <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .sync_impl()
  .Call(C_simd_sum, x, na.rm, na_check, .precision_code())
}

simd_prod <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .sync_impl()
  .check_supported(x, "simd_prod", "integer64")
  .Call(C_simd_prod, x, na.rm, na_check, .precision_code())
}

# prod(x + y) (op 0) and prod(x - y) (op 1) without the vector of sums; a
# complex operand makes the other complex, as in base R.
.simd_prod2 <- function(x, y, op, fun, na.rm, na_check) {
  .sync_impl()
  .check_supported(x, fun, "integer64")
  .check_supported(y, fun, "integer64", arg = "y")
  if (is.complex(x) != is.complex(y)) {
    if (is.numeric(x) || is.logical(x)) x <- .as_etype(x, "complex")
    if (is.numeric(y) || is.logical(y)) y <- .as_etype(y, "complex")
  }
  .Call(C_simd_prod2, x, y, op, na.rm, na_check, .precision_code())
}

simd_prod_sums <- function(x, y, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_prod2(x, y, 0L, "simd_prod_sums", na.rm, na_check)
}

simd_prod_diffs <- function(x, y, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_prod2(x, y, 1L, "simd_prod_diffs", na.rm, na_check)
}

simd_mean <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .sync_impl()
  .check_supported(x, "simd_mean", "integer64")
  .Call(C_simd_mean, x, na.rm, na_check, .precision_code())
}

# min, max and range share one kernel; op is 0, 1 or 2. Complex input is
# left to the C side, which rejects it with base R's message. With absval
# (max_abs and min_abs) the kernel reads abs(x), Mod(x) for complex x, and
# raw input is rejected here.
.simd_minmax <- function(x, op, fun, na.rm, na_check, absval = FALSE) {
  .sync_impl()
  if (absval) .check_supported(x, fun, "raw", character())
  .Call(C_simd_minmax, x, op, na.rm, na_check, absval, if (absval) .math_accuracy_code() else 0L)
}

simd_min <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_minmax(x, 0L, "simd_min", na.rm, na_check)
}

simd_max <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_minmax(x, 1L, "simd_max", na.rm, na_check)
}

simd_range <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_minmax(x, 2L, "simd_range", na.rm, na_check)
}

simd_max_abs <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_minmax(x, 1L, "simd_max_abs", na.rm, na_check, absval = TRUE)
}

simd_min_abs <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_minmax(x, 0L, "simd_min_abs", na.rm, na_check, absval = TRUE)
}

simd_which_min <- function(x) {
  .sync_impl()
  .check_supported(x, "simd_which_min", "complex", character())
  .Call(C_simd_which, x, FALSE, FALSE, 0L)
}

simd_which_max <- function(x) {
  .sync_impl()
  .check_supported(x, "simd_which_max", "complex", character())
  .Call(C_simd_which, x, TRUE, FALSE, 0L)
}

simd_which_min_abs <- function(x) {
  .sync_impl()
  .check_supported(x, "simd_which_min_abs", "raw", character())
  .Call(C_simd_which, x, FALSE, TRUE, .math_accuracy_code())
}

simd_which_max_abs <- function(x) {
  .sync_impl()
  .check_supported(x, "simd_which_max_abs", "raw", character())
  .Call(C_simd_which, x, TRUE, TRUE, .math_accuracy_code())
}

simd_any <- function(x, na.rm = FALSE) {
  .sync_impl()
  .check_supported(x, "simd_any", "complex", character())
  .Call(C_simd_anyall, x, FALSE, na.rm)
}

simd_all <- function(x, na.rm = FALSE) {
  .sync_impl()
  .check_supported(x, "simd_all", "complex", character())
  .Call(C_simd_anyall, x, TRUE, na.rm)
}

simd_any_na <- function(x) {
  .sync_impl()
  .Call(C_simd_na, x, 0L)
}

simd_count_na <- function(x) {
  .sync_impl()
  .Call(C_simd_na, x, 1L)
}

simd_which_na <- function(x) {
  .sync_impl()
  .Call(C_simd_na, x, 2L)
}

# which(x) and sum(x) of a logical x; other types are errors.
.not_logical <- c("double", "integer", "raw", "complex", "integer64")

simd_which <- function(x) {
  .sync_impl()
  .check_supported(x, "simd_which", .not_logical, character())
  .Call(C_simd_true, x, 1L, FALSE)
}

simd_count <- function(x, na.rm = FALSE) {
  .sync_impl()
  .check_supported(x, "simd_count", .not_logical, character())
  .Call(C_simd_true, x, 0L, na.rm)
}

# sum_sq, norm and sum_abs share one entry point; op is 0, 1 or 2.
.simd_sum_sq <- function(x, op, fun, na.rm, na_check) {
  .sync_impl()
  .check_supported(x, fun, c("integer64", "complex"), later = character())
  .Call(C_simd_sum_sq, x, op, na.rm, na_check, .precision_code())
}

simd_sum_sq <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_sum_sq(x, 0L, "simd_sum_sq", na.rm, na_check)
}

simd_norm <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_sum_sq(x, 1L, "simd_norm", na.rm, na_check)
}

simd_sum_abs <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_sum_sq(x, 2L, "simd_sum_abs", na.rm, na_check)
}

# dot, dist and cosine share one entry point; op is 0, 1 or 2.
.simd_dot <- function(x, y, op, fun, na.rm, na_check) {
  .sync_impl()
  .check_supported(x, fun, c("integer64", "complex"), later = character())
  .check_supported(y, fun, c("integer64", "complex"), later = character(), arg = "y")
  .Call(C_simd_dot, x, y, op, na.rm, na_check, .precision_code())
}

simd_dot <- function(x, y, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_dot(x, y, 0L, "simd_dot", na.rm, na_check)
}

simd_dist <- function(x, y, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_dot(x, y, 1L, "simd_dist", na.rm, na_check)
}

simd_cosine <- function(x, y, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_dot(x, y, 2L, "simd_cosine", FALSE, na_check)
}

.simd_var <- function(x, sd, fun, na.rm, na_check) {
  .sync_impl()
  .check_supported(x, fun, c("integer64", "complex"), later = character())
  .Call(C_simd_var, x, sd, na.rm, na_check, .precision_code())
}

simd_var <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_var(x, FALSE, "simd_var", na.rm, na_check)
}

simd_sd <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .simd_var(x, TRUE, "simd_sd", na.rm, na_check)
}

# The scans share one entry point; op is 0 (cumsum), 1 (cumprod), 2
# (cummin) or 3 (cummax). Only cumprod rejects integer64, and only cummin
# and cummax complex.
.simd_scan <- function(x, op, fun) {
  .sync_impl()
  .check_supported(x, fun, c(if (op == 1L) "integer64", if (op >= 2L) "complex"),
    later = character()
  )
  .Call(C_simd_scan, x, op, .precision_code())
}

simd_cumsum <- function(x) .simd_scan(x, 0L, "simd_cumsum")

simd_cumprod <- function(x) .simd_scan(x, 1L, "simd_cumprod")

simd_cummin <- function(x) .simd_scan(x, 2L, "simd_cummin")

simd_cummax <- function(x) .simd_scan(x, 3L, "simd_cummax")
