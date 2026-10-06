# Complex conjugate, real and imaginary parts, argument, and the load-time
# choice of the complex arithmetic variants. Other complex support lives in
# the functions that take complex input (simd_add, simd_sub, simd_mul,
# simd_div, simd_neg, simd_abs, simd_sum, simd_prod, simd_mean,
# simd_cumsum, simd_cumprod, simd_eq, simd_ne and the predicates).

simd_conj <- function(z) {
  if (is.object(z)) .check_supported(z, "simd_conj", "integer64", character(), "z")
  .Call(C_simd_cplx, z, "conj")
}

simd_re <- function(z) {
  if (is.object(z)) .check_supported(z, "simd_re", "integer64", character(), "z")
  .Call(C_simd_cplx, z, "re")
}

simd_im <- function(z) {
  if (is.object(z)) .check_supported(z, "simd_im", "integer64", character(), "z")
  .Call(C_simd_cplx, z, "im")
}

# Arg of complex z; double, integer and logical x are complex with a zero
# imaginary part, as in base R.
simd_arg <- function(z) {
  .check_supported(z, "simd_arg", "integer64", character(), "z")
  if (!is.complex(z)) {
    if (!is.numeric(z) && !is.logical(z)) stop("non-numeric argument to function", call. = FALSE)
    z <- .sv_like(as.complex(z), z)
  }
  .Call(C_simd_cmath, z, "arg", NULL)
}

# Base R rounds complex products, quotients and cumprod as the compiler
# that built R does. Base R computes them on a fixed set of inputs, and the
# C side picks the kernel variants that reproduce every result (see
# attr(simd_current(), "complex")).
init_complex <- function() {
  p <- .Call(C_simd_c128_probe_inputs)
  .Call(C_simd_c128_probe, p$z, p$w, p$z * p$w, p$z / p$w, p$x, cumprod(p$x))
  invisible()
}
