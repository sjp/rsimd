# Hamming distances: the number of positions where two vectors differ
# (simd_hamming) and the number of bits in which they differ
# (simd_hamming_bits).

# The operands are converted as for the comparisons (.cmp_operands()).
simd_hamming <- function(x, y, na.rm = FALSE) {
  .sync_impl()
  na.rm <- .arg_flag(na.rm, "na.rm")
  p <- .cmp_operands(x, y, sys.call())
  .Call(C_simd_hamming, p[[1L]], p[[2L]], na.rm)
}

# The width class of an operand of simd_hamming_bits: 32-bit integers
# (integer, logical), raw or integer64; anything else errors.
.bits_class <- function(x, arg) {
  if (is.object(x)) .check_data(x, arg)
  if (inherits(x, "integer64")) {
    return("integer64")
  }
  switch(typeof(x),
    integer = ,
    logical = "integer",
    raw = "raw",
    stop("simd_hamming_bits() does not support '", arg, "' of type ", typeof(x), call. = FALSE)
  )
}

simd_hamming_bits <- function(x, y) {
  .sync_impl()
  cx <- .bits_class(x, "x")
  cy <- .bits_class(y, "y")
  if (cx != cy) {
    stop("simd_hamming_bits() needs 'x' and 'y' of the same width: ", cx, " and ", cy,
      call. = FALSE
    )
  }
  .Call(C_simd_hamming_bits, x, y)
}
