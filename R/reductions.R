simd_sum <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .sync_impl()
  .check_supported(x, "simd_sum", c("integer64", "complex"))
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
