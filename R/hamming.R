# Hamming distances: the number of positions where two vectors differ
# (simd_hamming) and the number of bits in which they differ
# (simd_hamming_bits).

# The operands are converted as for the comparisons (raw with non-raw to
# logical or integer, integer64 with double to double, with a warning),
# and a non-complex operand with a complex one to complex (integer64
# through double, with the same warning).
simd_hamming <- function(x, y, na.rm = FALSE) {
  .sync_impl()
  na.rm <- .arg_flag(na.rm, "na.rm")
  if (is.complex(x) != is.complex(y)) {
    call <- sys.call()
    to <- function(a) {
      if (is.complex(a)) {
        return(a)
      }
      if (inherits(a, "integer64")) {
        warning(simpleWarning("integer64 coerced to double", call))
        a <- .as_etype(a, "double")
      }
      .sv_like(as.complex(a), a)
    }
    return(.Call(C_simd_hamming, to(x), to(y), na.rm))
  }
  if (is.raw(x) != is.raw(y)) {
    to <- if (is.logical(x) || is.logical(y)) as.logical else as.integer
    if (is.raw(x)) x <- .sv_like(to(x), x) else y <- .sv_like(to(y), y)
  }
  p <- .i64_to_double(list(x, y), sys.call())
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
