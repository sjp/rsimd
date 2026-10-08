# Type conversions between double, integer, logical, raw and integer64
# vectors. mode is passed to the C side as a code (0 "checked", 1
# "saturating", 2 "truncating"), partially matched as match.arg() does but
# with errors that name 'mode'.

.cvt_modes <- c("checked", "saturating", "truncating")

# The code of mode (missing, or the whole choice vector as match.arg()
# allows: the default, "checked").
.cvt_code <- function(mode, missing) {
  if (missing || identical(mode, .cvt_modes)) {
    return(0L)
  }
  if (!is.character(mode) || length(mode) != 1L) {
    .stop("'mode' must be a single string")
  }
  code <- pmatch(mode, .cvt_modes)
  if (is.na(code)) {
    .stop("'mode' must be one of \"checked\", \"saturating\" or \"truncating\"")
  }
  code - 1L
}

simd_as_integer <- function(x, mode = c("checked", "saturating", "truncating")) {
  code <- .cvt_code(mode, missing(mode))
  if (is.object(x) || is.complex(x)) .check_supported(x, "simd_as_integer", "complex", character(0))
  .Call(C_simd_convert, x, "integer", code, FALSE)
}

simd_as_double <- function(x) {
  if (is.object(x) || is.complex(x)) .check_supported(x, "simd_as_double", "complex", character(0))
  .Call(C_simd_convert, x, "double", 0L, FALSE)
}

simd_as_logical <- function(x) {
  if (is.object(x) || is.complex(x)) .check_supported(x, "simd_as_logical", "complex", character(0))
  .Call(C_simd_convert, x, "logical", 0L, FALSE)
}

simd_as_raw <- function(x, mode = c("checked", "saturating", "truncating")) {
  code <- .cvt_code(mode, missing(mode))
  if (is.object(x) || is.complex(x)) .check_supported(x, "simd_as_raw", "complex", character(0))
  .Call(C_simd_convert, x, "raw", code, FALSE)
}

simd_as_integer64 <- function(x, mode = c("checked", "saturating", "truncating")) {
  code <- .cvt_code(mode, missing(mode))
  if (is.object(x) || is.complex(x)) {
    .check_supported(x, "simd_as_integer64", "complex", character(0))
  }
  .Call(C_simd_convert, x, "integer64", code, FALSE)
}
