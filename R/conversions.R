# Type conversions between double, integer, logical, raw and integer64
# vectors.

.cvt_modes <- c("checked", "saturating", "truncating")

.convert <- function(x, to, mode, fun) {
  .sync_impl()
  .check_supported(x, fun, "complex", character(0))
  .Call(C_simd_convert, x, to, match(mode, .cvt_modes) - 1L, FALSE)
}

simd_as_integer <- function(x, mode = c("checked", "saturating", "truncating")) {
  .convert(x, "integer", match.arg(mode), "simd_as_integer")
}

simd_as_double <- function(x) .convert(x, "double", "checked", "simd_as_double")

simd_as_logical <- function(x) .convert(x, "logical", "checked", "simd_as_logical")

simd_as_raw <- function(x, mode = c("checked", "saturating", "truncating")) {
  .convert(x, "raw", match.arg(mode), "simd_as_raw")
}

simd_as_integer64 <- function(x, mode = c("checked", "saturating", "truncating")) {
  .convert(x, "integer64", match.arg(mode), "simd_as_integer64")
}
