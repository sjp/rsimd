# Type conversions between double, integer, logical and raw vectors.

.cvt_modes <- c("checked", "saturating", "truncating")

# Conversions that will take integer64 input: their messages say "yet".
.cvt_later_i64 <- c("simd_as_integer", "simd_as_double", "simd_as_logical")

.convert <- function(x, to, mode, fun) {
  .sync_impl()
  later <- if (fun %in% .cvt_later_i64) "integer64" else character(0)
  .check_supported(x, fun, c("integer64", "complex"), later)
  .Call(C_simd_convert, x, to, match(mode, .cvt_modes) - 1L)
}

simd_as_integer <- function(x, mode = c("checked", "saturating", "truncating")) {
  .convert(x, "integer", match.arg(mode), "simd_as_integer")
}

simd_as_double <- function(x) .convert(x, "double", "checked", "simd_as_double")

simd_as_logical <- function(x) .convert(x, "logical", "checked", "simd_as_logical")

simd_as_raw <- function(x, mode = c("checked", "saturating", "truncating")) {
  .convert(x, "raw", match.arg(mode), "simd_as_raw")
}
