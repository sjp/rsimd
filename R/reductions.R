simd_sum <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .sync_impl()
  .check_supported(x, "simd_sum", "integer64")
  .Call(C_simd_sum, x, na.rm, na_check, .precision_code())
}

simd_prod <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .sync_impl()
  .check_supported(x, "simd_prod", c("integer64", "complex"), later = "integer64")
  .Call(C_simd_prod, x, na.rm, na_check)
}

simd_mean <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .sync_impl()
  .check_supported(x, "simd_mean", c("integer64", "complex"), later = "integer64")
  .Call(C_simd_mean, x, na.rm, na_check, .precision_code())
}

# min, max and range share one kernel; op is 0, 1 or 2. Complex input is
# left to the C side, which rejects it with base R's message.
.simd_minmax <- function(x, op, fun, na.rm, na_check) {
  .sync_impl()
  .check_supported(x, fun, "integer64")
  .Call(C_simd_minmax, x, op, na.rm, na_check)
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

simd_which_min <- function(x) {
  .sync_impl()
  .check_supported(x, "simd_which_min", c("integer64", "complex"), later = "integer64")
  .Call(C_simd_which, x, FALSE)
}

simd_which_max <- function(x) {
  .sync_impl()
  .check_supported(x, "simd_which_max", c("integer64", "complex"), later = "integer64")
  .Call(C_simd_which, x, TRUE)
}

simd_any <- function(x, na.rm = FALSE) {
  .sync_impl()
  .check_supported(x, "simd_any", c("integer64", "complex"), later = "integer64")
  .Call(C_simd_anyall, x, FALSE, na.rm)
}

simd_all <- function(x, na.rm = FALSE) {
  .sync_impl()
  .check_supported(x, "simd_all", c("integer64", "complex"), later = "integer64")
  .Call(C_simd_anyall, x, TRUE, na.rm)
}

simd_any_na <- function(x) {
  .sync_impl()
  .check_supported(x, "simd_any_na", "integer64")
  .Call(C_simd_na, x, 0L)
}

simd_count_na <- function(x) {
  .sync_impl()
  .check_supported(x, "simd_count_na", "integer64")
  .Call(C_simd_na, x, 1L)
}

simd_which_na <- function(x) {
  .sync_impl()
  .check_supported(x, "simd_which_na", "integer64")
  .Call(C_simd_na, x, 2L)
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
# (cummin) or 3 (cummax).
.simd_scan <- function(x, op, fun, later = character()) {
  .sync_impl()
  .check_supported(x, fun, c("integer64", "complex"), later = later)
  .Call(C_simd_scan, x, op, .precision_code())
}

simd_cumsum <- function(x) .simd_scan(x, 0L, "simd_cumsum", later = "integer64")

simd_cumprod <- function(x) .simd_scan(x, 1L, "simd_cumprod")

simd_cummin <- function(x) .simd_scan(x, 2L, "simd_cummin")

simd_cummax <- function(x) .simd_scan(x, 3L, "simd_cummax")
