# Operands of a binary operation converted to their common element type
# (rsimd_promote() in src/rvec.c), as list(x, y). Only an operand whose type
# differs from the common type is converted, so same-type operands are
# passed through without a copy; two logicals stay logical (their storage
# is integer). Converting integer64 to double warns, naming the caller.
.promote_pair <- function(x, y) {
  types <- .Call(C_simd_promote, x, y)
  to <- types[[1L]]
  if (to == "double" && "integer64" %in% types[2:3]) {
    warning(simpleWarning("integer64 coerced to double", sys.call(-1L)))
  }
  list(
    x = if (types[[2L]] == to) x else .as_etype(x, to),
    y = if (types[[3L]] == to) y else .as_etype(y, to)
  )
}

# x converted to the element type named by `to` (a type name from
# C_simd_promote).
.as_etype <- function(x, to) {
  switch(to,
    double = {
      if (inherits(x, "integer64")) .need_bit64()
      as.double(x)
    },
    integer = as.integer(x),
    complex = as.complex(x),
    integer64 = {
      .need_bit64()
      bit64::as.integer64(x)
    },
    stop("internal error: cannot convert to ", to, call. = FALSE)
  )
}

# Loads bit64, whose methods convert integer64 values; integer64 inputs
# exist only when it is installed.
.need_bit64 <- function() {
  if (!requireNamespace("bit64", quietly = TRUE)) {
    stop("package 'bit64' is needed to convert integer64 values", call. = FALSE)
  }
  invisible()
}

# Errors if x is of one of the element types in `unsupported` ("integer64",
# "complex"), which function `fun` does not take for its argument `arg`;
# the message says "yet" for the types in `later`, which it is planned to
# take.
.check_supported <- function(x, fun, unsupported, later = unsupported, arg = "x") {
  type <- if (inherits(x, "integer64")) "integer64" else typeof(x)
  if (type %in% unsupported) {
    stop(fun, "() does not support '", arg, "' of type ", type, if (type %in% later) " yet",
      call. = FALSE
    )
  }
  invisible()
}
