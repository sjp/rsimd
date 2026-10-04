# Complex conjugate, real and imaginary parts. Other complex support lives
# in the functions that take complex input (simd_add, simd_sub, simd_neg,
# simd_sum and the predicates is_na, is_nan, is_finite, is_infinite).

.cplx <- function(z, op, fun) {
  .sync_impl()
  .check_supported(z, fun, "integer64", character(), "z")
  .Call(C_simd_cplx, z, op)
}

simd_conj <- function(z) .cplx(z, "conj", "simd_conj")

simd_re <- function(z) .cplx(z, "re", "simd_re")

simd_im <- function(z) .cplx(z, "im", "simd_im")
