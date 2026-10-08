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
# C_simd_promote). integer64 is converted natively (an integer64 result
# loads bit64), without the precision warning: the caller warns about the
# coercion. A
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
  i64 <- dbl <- FALSE
  for (a in args) {
    if (inherits(a, "integer64")) i64 <- TRUE else if (is.double(a)) dbl <- TRUE
  }
  if (!i64 || !(dbl || always)) {
    return(args)
  }
  if (dbl) warning(simpleWarning("integer64 coerced to double", call))
  for (i in seq_along(args)) {
    if (inherits(args[[i]], "integer64")) args[[i]] <- .as_etype(args[[i]], "double")
  }
  args
}

# Errors if x has a class rsimd does not take (.check_data()), or is of one
# of the element types in `unsupported` ("integer64", "complex"), which function
# `fun` does not take for its argument `arg`; the message says "yet" for
# the types in `later`, which it is planned to take.
.check_supported <- function(x, fun, unsupported, later = unsupported, arg = "x") {
  if (is.object(x)) .check_data(x, arg)
  type <- if (inherits(x, "integer64")) "integer64" else typeof(x)
  if (type %in% unsupported) {
    stop(fun, "() does not support '", arg, "' of type ", type, if (type %in% later) " yet",
      call. = FALSE
    )
  }
  invisible()
}

# Errors unless x is data rsimd takes: an atomic vector of a supported type
# without a class, an integer64 or a simd_vec. The message is the one
# rsimd_check_atomic() in src/rvec.c gives, naming the argument `arg`.
.check_data <- function(x, arg = "x") {
  types <- c("double", "integer", "logical", "raw", "complex")
  ok <- inherits(x, c("integer64", "simd_vec")) ||
    (is.atomic(x) && !is.object(x) && typeof(x) %in% types)
  if (!ok) {
    what <- if (is.object(x)) class(x)[[1L]] else if (is.list(x)) "list" else typeof(x)
    stop("'", arg, "' must be an atomic vector (double, integer, logical, raw, complex or ",
      "integer64), not ", what,
      call. = FALSE
    )
  }
  invisible()
}

# Loads bit64, so that its S3 methods are registered before rsimd hands out
# an integer64 object (without them it sorts, prints and compares as the
# doubles its bits make); errors if it is not installed. The C side does
# the same for the integer64 results of the simd_* functions.
.need_bit64 <- function() {
  if (!isNamespaceLoaded("bit64") && !requireNamespace("bit64", quietly = TRUE)) {
    stop("integer64 results need package 'bit64'; install it with install.packages(\"bit64\")",
      call. = FALSE
    )
  }
  invisible()
}
