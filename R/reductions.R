simd_sum <- function(x, na.rm = FALSE, na_check = getOption("rsimd.na_check", TRUE)) {
  .sync_impl()
  .check_supported(x, "simd_sum", c("integer64", "complex"))
  .Call(C_simd_sum, x, na.rm, na_check, .precision_code())
}
