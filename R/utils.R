# Operands of a binary operation converted to their common element type
# (rsimd_promote() in src/rvec.c), as list(x, y). Only an operand whose type
# differs from the common type is converted, so same-type operands are
# passed through without a copy; two logicals stay logical (their storage
# is integer). Converting integer64 to double warns, naming the caller.
.promote_pair <- function(x, y, call = sys.call(-1L)) {
  types <- .Call(C_simd_promote, x, y)
  to <- types[[1L]]
  if (to == "double" && "integer64" %in% types[2:3]) {
    warning(simpleWarning("integer64 coerced to double", call))
  }
  list(
    x = if (types[[2L]] == to) x else .as_etype(x, to),
    y = if (types[[3L]] == to) y else .as_etype(y, to)
  )
}

# x converted to the element type named by `to` (a type name from
# C_simd_promote). integer64 is converted natively (bit64 is not needed),
# without the precision warning: the caller warns about the coercion. A
# simd_vec stays one (with its pin and NA-free flag).
.as_etype <- function(x, to) {
  out <- switch(to,
    double = {
      if (inherits(x, "integer64")) .Call(C_simd_convert, x, "double", 0L, TRUE) else as.double(x)
    },
    integer = as.integer(x),
    complex = as.complex(x),
    integer64 = .Call(C_simd_convert, x, "integer64", 0L, FALSE),
    stop("internal error: cannot convert to ", to, call. = FALSE)
  )
  .sv_like(out, x)
}

# The operands in `args` (a list) with every integer64 one converted to
# double when another operand is a double, warning "integer64 coerced to
# double" with `call`; with `always` (simd_div, whose result is double)
# integer64 operands are converted in any case, the warning still being
# given only when a double is present.
.i64_to_double <- function(args, call, always = FALSE) {
  i64 <- vapply(args, inherits, NA, what = "integer64")
  if (!any(i64)) {
    return(args)
  }
  dbl <- vapply(args, function(a) is.double(a) && !inherits(a, "integer64"), NA)
  if (!any(dbl) && !always) {
    return(args)
  }
  if (any(dbl)) warning(simpleWarning("integer64 coerced to double", call))
  args[i64] <- lapply(args[i64], .as_etype, to = "double")
  args
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
