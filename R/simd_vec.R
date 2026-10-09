# The simd_vec wrapper: an atomic vector with class "simd_vec" (or
# c("simd_vec", "integer64")), an optional pinned implementation (attribute
# rsimd_impl) and a known-NA-free flag (TRUE or FALSE, held by attribute
# rsimd_na_token; absent for unknown).
#
# Base R copies attributes onto new data (pmin(), storage.mode<-, ...), so
# the flag is only read through .sv_flag(), which trusts it only while its
# token is this object's own and unshared and the class is the one it was
# set with (see rsimd_sv_flag() in src/rvec.h), and only set through
# .Call(C_simd_sv_stamp, ...), which makes a fresh token. The methods here
# strip a simd_vec with .subset(), never by copying its attributes, so that
# they do not share the token of the original and void its flag.
#
# The simd_* functions themselves handle simd_vec operands on the C side
# (src/rvec.c): they run under the operands' pin, skip NA checks for
# operands flagged NA-free and wrap value-vector results again. The methods
# here map R's generics onto those functions, and do in R what has no
# kernel (base fallbacks, subsetting, c(), printing).

# ---- Construction and accessors -------------------------------------------

simd_vec <- function(x, impl = NULL, scan_na = FALSE) {
  if (!is.logical(scan_na) || length(scan_na) != 1L || is.na(scan_na)) {
    .stop("'scan_na' must be TRUE or FALSE")
  }
  if (is_simd_vec(x)) {
    if (missing(impl) && !scan_na) {
      return(x)
    }
    impl <- if (missing(impl)) attr(x, "rsimd_impl", exact = TRUE) else .sv_check_impl(impl)
    na_free <- if (scan_na) !simd_any_na(x) else .sv_flag(x)
    return(.sv_new(.sv_strip(x), impl, na_free))
  }
  data <- .sv_strip(x)
  impl <- .sv_check_impl(impl)
  na_free <- if (scan_na) !simd_any_na(data) else NULL
  .sv_new(data, impl, na_free)
}

as_simd_vec <- function(x, ...) UseMethod("as_simd_vec")

as_simd_vec.default <- function(x, ...) simd_vec(x, ...)

as_simd_vec.simd_vec <- function(x, ...) simd_vec(x, ...)

is_simd_vec <- function(x) inherits(x, "simd_vec")

simd_impl <- function(x) {
  .sv_assert(x)
  attr(x, "rsimd_impl", exact = TRUE)
}

`simd_impl<-` <- function(x, value) {
  .sv_assert(x)
  .sv_release(x)
  flag <- .sv_flag(x)
  .sv_new(.sv_data_own(x), .sv_check_impl(value), flag)
}

simd_na_free <- function(x) {
  .sv_assert(x)
  .sv_flag(x)
}

# ---- Internal helpers ------------------------------------------------------

# The NA-free flag of simd_vec x: TRUE, FALSE, or NULL when unknown or not
# valid for this object.
.sv_flag <- function(x) .Call(C_simd_sv_flag, x)

# Replacement methods start with .sv_release(x): for f(y) <- value on a
# shared y, R passes a copy of y, and the token it shares with y would void
# y's flag for good. Removes, in place, a token that is not x's own.
.sv_release <- function(x) invisible(.Call(C_simd_sv_release, x))

.sv_assert <- function(x, arg = "x") {
  if (!is_simd_vec(x)) .stop("'", arg, "' must be a simd_vec")
  invisible()
}

# A pin value checked: NULL, or the name of an available tier.
.sv_check_impl <- function(impl) {
  if (is.null(impl)) {
    return(NULL)
  }
  if (identical(impl, "auto")) {
    .stop("'impl' must name an implementation, not \"auto\"; use NULL to follow the ",
      "global setting")
  }
  problem <- .impl_problem(impl)
  if (!is.null(problem)) .stop(problem)
  impl
}

# The data of x without attributes, except an integer64 class. Errors for
# anything that is not a supported atomic vector.
.sv_strip <- function(x) {
  .check_data(x)
  i64 <- inherits(x, "integer64")
  if (is.object(x) || !is.null(attributes(x))) {
    # Not attributes(x) <- NULL: on a shared x that copies the attributes,
    # and so the NA-free token (see the top of this file).
    # Bound to x alone, so that class<- below changes it in place.
    orig <- x
    x <- .Call(C_simd_sv_bare, orig)
    if (is.null(x)) {
      x <- .subset(orig, seq_along(orig))
      names(x) <- NULL
    }
  }
  if (i64) class(x) <- "integer64"
  x
}

# data (bare, or of class integer64) as a simd_vec, with NA-free flag
# na_free or, with scan_na, the one a scan finds. The attributes are set by
# one call of the replacement function, not bound to a name in between: a
# second replacement, or C_simd_sv_stamp on a bound value, would see data as
# shared and copy it. A fresh result is so changed in place, and a shared
# one is copied by R (lazily, through an ALTREP wrapper, for a long vector).
.sv_new <- function(data, impl = NULL, na_free = NULL, scan_na = FALSE) {
  i64 <- inherits(data, "integer64")
  if (i64) .need_bit64()
  if (scan_na) na_free <- !simd_any_na(data)
  cls <- if (i64) c("simd_vec", "integer64") else "simd_vec"
  .Call(
    C_simd_sv_stamp, `attributes<-`(data, list(class = cls, rsimd_impl = impl)),
    if (is.raw(data)) TRUE else na_free
  )
}

# The data of a simd_vec (integer64 keeps its class), as .sv_strip() but
# for replacement methods, after .sv_release(x): x is then R's
# private copy (f(y) <- value on a shared y) or the object being replaced,
# so its attributes are dropped in place, saving a copy of the data.
.sv_data_own <- function(x) {
  i64 <- inherits(x, "integer64")
  attributes(x) <- NULL
  if (i64) class(x) <- "integer64"
  x
}

# out (a bare result computed from x by base R) as a simd_vec like x: same
# pin, and the given NA-free flag. Returned unchanged when x is not a
# simd_vec or out already is one.
.sv_like <- function(out, x, na_free = .sv_flag(x)) {
  if (!is_simd_vec(x) || is_simd_vec(out)) {
    return(out)
  }
  .sv_new(out, attr(x, "rsimd_impl", exact = TRUE), na_free)
}

# The error for operands pinned to different tiers a and b.
.sv_pin_clash <- function(a, b) {
  .stop("operands pinned to different implementations ('", a, "' vs '", b,
    "'); unpin one with simd_impl(x) <- NULL")
}

# The pin shared by the simd_vec arguments (NULL when none is pinned), with
# the same rules and messages as the C side (src/rvec.c): different pins
# and pins to unavailable tiers are errors.
.sv_resolve <- function(args) {
  pin <- NULL
  for (a in args) {
    if (!is_simd_vec(a)) next
    p <- attr(a, "rsimd_impl", exact = TRUE)
    if (is.null(p)) next
    if (!is.null(pin) && !identical(pin, p)) .sv_pin_clash(pin, p)
    if (!p %in% simd_available()) {
      .stop("'x' is pinned to implementation '", p,
        "', which is not available on this machine; unpin it with simd_impl(x) <- NULL")
    }
    pin <- p
  }
  pin
}

# A base R function applied to the data of simd_vec arguments, the result
# wrapped with the arguments' common pin and an unknown NA-free flag.
# integer64 data dispatches to bit64's methods, which must be available.
.sv_fallback <- function(f, args, wrap = TRUE) {
  pin <- .sv_resolve(args)
  data <- lapply(args, function(a) if (is_simd_vec(a)) .sv_strip(a) else a)
  if (any(vapply(data, inherits, NA, what = "integer64"))) .need_bit64()
  out <- .relay(do.call(f, data))
  if (wrap && ((is.atomic(out) && !is.object(out)) || inherits(out, "integer64"))) {
    out <- .sv_new(.sv_strip(out), pin)
  }
  out
}

.sv_type <- function(x) if (inherits(x, "integer64")) "integer64" else typeof(x)

# Length rule of the binary operators, for base fallbacks: a
# zero-length operand gives a zero-length result, as in base R.
.sv_check_lengths <- function(e1, e2) {
  n1 <- length(e1)
  n2 <- length(e2)
  if (n1 != n2 && n1 > 1L && n2 > 1L) {
    .stop("lengths of 'x' (", n1, ") and 'y' (", n2,
      ") must be equal or one of them must be 1")
  }
  invisible()
}

# ---- Group generics --------------------------------------------------------

# A simd_vec operand wins over another class's Ops method (R >= 4.3), so
# that Ops.simd_vec can reject an operand of another class (Date, difftime
# ...): were both sides to return FALSE, R would warn "Incompatible methods"
# and run the internal operator on the bare data. bit64 also always claims
# the operation, so a plain integer64 on the left of a simd_vec still
# dispatches to bit64.
chooseOpsMethod.simd_vec <- function(x, y, mx, my, cl, reverse) TRUE

Ops.simd_vec <- function(e1, e2) {
  unary <- nargs() == 1L
  gen <- .Generic
  if (unary) {
    return(switch(gen,
      "+" = .sv_uplus(e1),
      "-" = simd_neg(e1),
      "!" = if (is.raw(e1)) simd_bit_not(e1) else simd_not(e1),
      .stop("invalid unary operator")
    ))
  }
  if (!(.sv_operand_ok(e1) && .sv_operand_ok(e2))) {
    .check_data(e1, "x")
    .check_data(e2, "y")
  }
  switch(gen,
    "+" = simd_add(e1, e2),
    "-" = simd_sub(e1, e2),
    "*" = simd_mul(e1, e2),
    "/" = simd_div(e1, e2),
    "%%" = if (.sv_cplx(e1, e2)) .sv_ops_fallback(gen, e1, e2) else simd_mod(e1, e2),
    "%/%" = if (.sv_cplx(e1, e2)) .sv_ops_fallback(gen, e1, e2) else simd_idiv(e1, e2),
    "^" = if (inherits(e1, "integer64") || inherits(e2, "integer64")) {
      .sv_ops_fallback(gen, e1, e2)
    } else {
      simd_pow(e1, e2)
    },
    "==" = if (.sv_cplx(e1, e2)) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_eq(e1, e2),
    "!=" = if (.sv_cplx(e1, e2)) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_ne(e1, e2),
    "<" = if (.sv_cplx(e1, e2)) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_lt(e1, e2),
    ">" = if (.sv_cplx(e1, e2)) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_gt(e1, e2),
    "<=" = if (.sv_cplx(e1, e2)) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_le(e1, e2),
    ">=" = if (.sv_cplx(e1, e2)) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_ge(e1, e2),
    "&" = if (is.raw(e1) || is.raw(e2)) simd_bit_and(e1, e2) else simd_and(e1, e2),
    "|" = if (is.raw(e1) || is.raw(e2)) simd_bit_or(e1, e2) else simd_or(e1, e2),
    .stop("operator '", gen, "' is not supported for simd_vec")
  )
}

# TRUE for an operand .check_data() accepts, tested with a few primitives
# (.check_data() itself is the slower path that names a bad operand).
.sv_operand_ok <- function(x) {
  if (is.object(x)) {
    inherits(x, c("simd_vec", "integer64"))
  } else {
    is.numeric(x) || is.logical(x) || is.complex(x) || is.raw(x)
  }
}

.sv_cplx <- function(e1, e2) is.complex(e1) || is.complex(e2)

# Unary +: base R makes a logical integer and rejects raw.
.sv_uplus <- function(x) {
  if (is.raw(x)) .stop("invalid argument to unary operator")
  if (is.logical(x)) simd_as_integer(x) else x
}

.sv_ops_fallback <- function(gen, e1, e2, wrap = TRUE) {
  .sv_check_lengths(e1, e2)
  .sv_fallback(get(gen, envir = baseenv()), list(e1, e2), wrap)
}

# Math functions with a kernel, by the element types it takes besides
# double, integer and logical.
.sv_math_kernels <- list(
  abs = simd_abs, sign = simd_sign, sqrt = simd_sqrt, floor = simd_floor,
  ceiling = simd_ceiling, trunc = simd_trunc, exp = simd_exp, expm1 = simd_expm1,
  log2 = simd_log2, log10 = simd_log10, log1p = simd_log1p, cos = simd_cos, sin = simd_sin,
  tan = simd_tan, cospi = simd_cospi, sinpi = simd_sinpi, tanpi = simd_tanpi,
  acos = simd_acos, asin = simd_asin, atan = simd_atan, cosh = simd_cosh, sinh = simd_sinh,
  tanh = simd_tanh, acosh = simd_acosh, asinh = simd_asinh, atanh = simd_atanh,
  cumsum = simd_cumsum, cumprod = simd_cumprod, cummax = simd_cummax, cummin = simd_cummin
)

# Math functions that take integer64 input in rsimd.
.sv_math_i64 <- c("abs", "sign", "cumsum", "cummax", "cummin")

# Math functions that take complex input in rsimd (log too).
.sv_math_complex <- c(
  "abs", "cumsum", "cumprod", "sqrt", "exp", "log2", "log10", "cos", "sin", "tan", "acos",
  "asin", "atan", "cosh", "sinh", "tanh", "acosh", "asinh", "atanh"
)

Math.simd_vec <- function(x, ...) {
  gen <- .Generic
  type <- .sv_type(x)
  kernel <- type %in% c("double", "integer", "logical", "raw") ||
    (type == "integer64" && gen %in% .sv_math_i64) ||
    (type == "complex" && gen %in% c(.sv_math_complex, "log"))
  if (kernel && gen == "round") {
    return(simd_round(x, ...))
  }
  if (kernel && gen == "log") {
    return(simd_log(x, ...))
  }
  if (kernel && !is.null(f <- .sv_math_kernels[[gen]])) {
    return(f(x))
  }
  # Without a bit64 method the primitive would read the integer64 bits as
  # doubles, so such data goes to double first, as bit64's sqrt() and log().
  if (type == "integer64") {
    .need_bit64()
    if (is.null(getS3method(gen, "integer64", optional = TRUE))) {
      return(.relay(get(gen, envir = baseenv())(simd_as_double(x), ...)))
    }
  }
  f <- get(gen, envir = baseenv())
  .sv_fallback(function(x) f(x, ...), list(x))
}

# Summary functions with a kernel, and the element types it rejects.
.sv_summary_kernels <- list(
  sum = list(simd_sum, character()),
  prod = list(simd_prod, "integer64"),
  min = list(simd_min, "complex"),
  max = list(simd_max, "complex"),
  range = list(simd_range, "complex"),
  any = list(simd_any, "complex"),
  all = list(simd_all, "complex")
)

Summary.simd_vec <- function(..., na.rm = FALSE) {
  gen <- .Generic
  args <- list(...)
  # range() takes 'finite' (an argument of range.default) through '...'.
  finite <- FALSE
  if (gen == "range" && "finite" %in% names(args)) {
    finite <- args[["finite"]]
    args <- args[names(args) != "finite"]
  }
  x <- if (length(args) == 1L) args[[1L]] else .sv_combine(args)
  k <- .sv_summary_kernels[[gen]]
  if (isFALSE(finite) && is_simd_vec(x) && !(.sv_type(x) %in% k[[2L]])) {
    return(k[[1L]](x, na.rm = na.rm))
  }
  f <- get(gen, envir = baseenv())
  g <- if (isFALSE(finite)) {
    function(x) f(x, na.rm = na.rm)
  } else {
    function(x) f(x, na.rm = na.rm, finite = finite)
  }
  .sv_fallback(g, list(x), wrap = FALSE)
}

# Re, Im, Mod, Arg and Conj: kernels for complex data, base R otherwise
# (whose result types for numeric input differ from the kernels').
Complex.simd_vec <- function(z) {
  if (!is.complex(z)) {
    f <- get(.Generic, envir = baseenv())
    return(.sv_fallback(f, list(z)))
  }
  switch(.Generic,
    Re = simd_re(z),
    Im = simd_im(z),
    Mod = simd_abs(z),
    Arg = simd_arg(z),
    Conj = simd_conj(z)
  )
}

mean.simd_vec <- function(x, trim = 0, na.rm = FALSE, ...) {
  if (!identical(trim, 0) || .sv_type(x) == "integer64") {
    return(.sv_fallback(function(x) mean(x, trim = trim, na.rm = na.rm, ...), list(x),
      wrap = FALSE
    ))
  }
  simd_mean(x, na.rm = na.rm)
}

# Registered for stats' generic (NAMESPACE). median.default sorts with [,
# which keeps the class, so it would return a length-one simd_vec.
median.simd_vec <- function(x, na.rm = FALSE, ...) {
  .sv_fallback(function(x) stats::median(x, na.rm = na.rm, ...), list(x), wrap = FALSE)
}

is.na.simd_vec <- function(x) simd_is_na(x)

# The data, for order() and so sort(): xtfrm.default would unclass(x),
# copying the attributes and so voiding x's NA-free flag.
xtfrm.simd_vec <- function(x) xtfrm(.sv_strip(x))

# The data, for match() and so %in%: mtfrm.default would compare through
# as.character(x), i.e. 15 significant digits. integer64 data keeps bit64's
# transform.
mtfrm.simd_vec <- function(x) {
  data <- simd_unwrap(x)
  if (inherits(data, "integer64")) mtfrm(data) else data
}

anyNA.simd_vec <- function(x, recursive = FALSE) {
  flag <- .sv_flag(x)
  if (isTRUE(flag)) {
    return(FALSE)
  }
  if (isFALSE(flag)) {
    return(TRUE)
  }
  simd_any_na(x)
}

# ---- Subsetting, combining, length -----------------------------------------

`[.simd_vec` <- function(x, i, ..., drop = TRUE) {
  # ...length() counts empty arguments, so x[, 1] and x[, ] fail here too.
  if (...length() > 0L) .stop("incorrect number of dimensions")
  if (missing(i)) {
    return(x)
  }
  impl <- attr(x, "rsimd_impl", exact = TRUE)
  # The result of an NA-free x has an NA only where the index is NA or past
  # the end; one pass over the (fresh) result says which, so the index is
  # never looked at.
  check <- isTRUE(.sv_flag(x))
  if (!inherits(x, "integer64")) {
    return(.sv_new(.subset(x, i), impl, scan_na = check))
  }
  j <- seq_along(x)[i]
  .sv_new(.sv_take(x, j), impl, if (check) !anyNA(j))
}

# The data of x (a simd_vec or its data) at positions j (NA giving a
# missing element), integer64 too, without copying the rest.
.sv_take <- function(x, j) {
  out <- .subset(x, j)
  if (!inherits(x, "integer64")) {
    return(out)
  }
  if (anyNA(j)) out[is.na(j)] <- unclass(.sv_na_i64())
  class(out) <- "integer64"
  out
}

.sv_na_i64 <- function() simd_as_integer64(NA)

`[<-.simd_vec` <- function(x, i, ..., value) {
  if (...length() > 0L) .stop("incorrect number of subscripts on matrix")
  # Checked before R coerces the data to the type of value, which would
  # then be blamed on x. NULL keeps base R's "replacement has length zero".
  if (!is.null(value)) .check_data(value, "value")
  .sv_release(x)
  # value's pin must match a pin of x, as in c(); an unpinned x stays
  # unpinned, and no pin is checked for availability (no kernel runs).
  if (is_simd_vec(value)) {
    pin <- attr(x, "rsimd_impl", exact = TRUE)
    p <- attr(value, "rsimd_impl", exact = TRUE)
    if (!is.null(pin) && !is.null(p) && !identical(pin, p)) .sv_pin_clash(pin, p)
    value <- .sv_strip(value)
  }
  data <- .sv_data_own(x)
  i64 <- inherits(data, "integer64")
  if (i64) {
    if (!is.null(value) && !inherits(value, "integer64")) value <- simd_as_integer64(value)
    data <- unclass(data)
    value <- unclass(value)
    n <- length(data)
  } else if (inherits(value, "integer64")) {
    .stop("cannot assign integer64 values into a ", typeof(data), " simd_vec")
  }
  # Base R's errors and warnings (a bad subscript, a value whose length
  # does not fit) name the user's call; one value at a logical or single
  # subscript raises none, and so skips .relay()'s handlers.
  plain <- length(value) == 1L &&
    (missing(i) || (!is.object(i) && (is.logical(i) || (is.numeric(i) && length(i) == 1L))))
  if (missing(i)) {
    if (plain) data[] <- value else .relay(data[] <- value)
  } else if (plain) {
    data[i] <- value
  } else {
    .relay(data[i] <- value)
  }
  if (i64) {
    if (length(data) > n) {
      # The gap R filled with double NA becomes integer64 NA.
      gap <- logical(n)
      gap[i] <- TRUE
      data[is.na(gap)] <- unclass(.sv_na_i64())
    }
    class(data) <- "integer64"
  }
  .sv_new(.sv_strip(data), attr(x, "rsimd_impl", exact = TRUE))
}

# x[[i]] <- value is x[i] <- value for a single element, with base R's
# checks (the default method would keep the attributes, and so a stale
# NA-free flag).
`[[<-.simd_vec` <- function(x, i, ..., value) {
  if (...length() > 0L) .stop("[[ ]] improper number of subscripts")
  .sv_release(x)
  if (length(i) != 1L) {
    what <- if (length(i)) "more" else "less"
    where <- if (length(i)) "vectorIndex" else "OneIndex"
    .stop("attempt to select ", what, " than one element in ", where)
  }
  if (length(value) != 1L) {
    if (length(value)) .stop("more elements supplied than there are to replace")
    .stop("replacement has length zero")
  }
  `[<-.simd_vec`(x, i, value = value)
}

`[[.simd_vec` <- function(x, i, ...) {
  if (...length() > 0L) .stop("incorrect number of subscripts")
  # An index other than one in range goes through .relay(), for base R's
  # error to name the user's call; the handlers are too slow for every call.
  if (inherits(x, "integer64")) {
    plain <- is.numeric(i) && !is.object(i) && length(i) == 1L && !is.na(i) && i >= 1 &&
      i < length(x) + 1
    return(.sv_take(x, if (plain) seq_along(x)[[i]] else .relay(seq_along(x)[[i]])))
  }
  # .subset() gives NA rather than an error out of range, so a single
  # element that is not NA is the answer; anything else is left to
  # .subset2().
  if (is.numeric(i)) {
    out <- .subset(x, i)
    if (length(out) == 1L && !is.na(out)) return(out)
  }
  .relay(.subset2(x, i))
}

# recursive and use.names change nothing for atomic parts without names.
c.simd_vec <- function(..., recursive = FALSE, use.names = TRUE) .sv_combine(list(...))

# c() of simd_vec and plain vectors: base R's type promotion (integer64
# absorbing the integer, logical and double parts, as bit64's c() does),
# the parts' common pin, and NA-free when every part is known to be.
.sv_combine <- function(args) {
  args <- args[lengths(args) > 0L | vapply(args, is_simd_vec, NA)]
  pin <- .sv_resolve(args)
  is_sv <- vapply(args, is_simd_vec, NA)
  for (a in args[!is_sv]) .check_data(a)
  flags <- lapply(args, function(a) if (is_simd_vec(a)) .sv_flag(a))
  # Checked before the integer64 conversion, whose NA bits read as -0.
  plain_na <- any(vapply(args[!is_sv], simd_any_na, NA))
  na_free <- if (any(vapply(flags, isFALSE, NA))) {
    FALSE
  } else if (all(vapply(flags[is_sv], isTRUE, NA))) {
    if (!plain_na) TRUE
  }
  # The result goes to .sv_new() unbound, so that it is not copied there.
  if (!length(args)) {
    return(.sv_new(logical(), pin, na_free))
  }
  if (any(vapply(args, inherits, NA, what = "integer64"))) {
    data <- lapply(args, function(a) {
      d <- .sv_strip(a)
      unclass(if (inherits(d, "integer64")) d else simd_as_integer64(d))
    })
    return(.sv_new(`class<-`(unlist(data, use.names = FALSE), "integer64"), pin, na_free))
  }
  # unlist() reads the data of the parts and ignores their attributes, so
  # they need no stripping (a copy each); except a part that is an ALTREP
  # wrapper (simd_vec(x) of a shared x), which unlist() reads element by
  # element, several times slower than the copy.
  args[is_sv] <- lapply(args[is_sv], function(a) if (.Call(C_simd_sv_altrep, a)) .sv_strip(a) else a)
  .sv_new(unlist(args, use.names = FALSE), pin, na_free)
}

`length<-.simd_vec` <- function(x, value) {
  .sv_release(x)
  flag <- .sv_flag(x)
  data <- .sv_data_own(x)
  n <- length(data)
  i64 <- inherits(data, "integer64")
  if (i64) data <- unclass(data)
  # A value other than one count goes through .relay(), for base R's
  # errors and warnings to name the user's call.
  if (is.numeric(value) && !is.object(value) && length(value) == 1L && !is.na(value)) {
    length(data) <- value
  } else {
    .relay(length(data) <- value)
  }
  if (i64) {
    if (length(data) > n) data[(n + 1L):length(data)] <- unclass(.sv_na_i64())
    class(data) <- "integer64"
  }
  na_free <- if (isTRUE(flag) && length(data) <= n) TRUE else NULL
  .sv_new(data, attr(x, "rsimd_impl", exact = TRUE), na_free)
}

# A simd_vec has no names, dimensions or dimnames. Setting them unwraps:
# the result is the plain data with the value set, so that base functions
# that name their result (quantile(), summary(), setNames(), lm()) return
# what they do for plain data. Setting NULL leaves the simd_vec as it is.
`names<-.simd_vec` <- function(x, value) .sv_set_attr(x, value, `names<-`)

`dim<-.simd_vec` <- function(x, value) .sv_set_attr(x, value, `dim<-`)

`dimnames<-.simd_vec` <- function(x, value) .sv_set_attr(x, value, `dimnames<-`)

.sv_set_attr <- function(x, value, set) {
  .sv_release(x)
  if (is.null(value)) {
    return(x)
  }
  # Not .sv_data_own(): x may be shared with the object it was copied from,
  # whose token a copy of the attributes would void.
  data <- .sv_strip(x)
  if (inherits(data, "integer64")) .need_bit64()
  .relay(set(data, value))
}

# ---- Unwrapping ------------------------------------------------------------

simd_unwrap <- function(x) {
  if (!is_simd_vec(x)) {
    .check_data(x)
    return(x)
  }
  data <- .sv_strip(x)
  if (inherits(data, "integer64")) .need_bit64()
  data
}

# Registered for bit64's generic when bit64 is loaded (NAMESPACE).
as.integer64.simd_vec <- function(x, ...) {
  data <- .sv_strip(x)
  if (inherits(data, "integer64")) data else bit64::as.integer64(data, ...)
}

as.vector.simd_vec <- function(x, mode = "any") {
  data <- .sv_strip(x)
  if (inherits(data, "integer64")) data <- simd_as_double(data)
  .relay(as.vector(data, mode))
}

as.double.simd_vec <- function(x, ...) {
  data <- .sv_strip(x)
  if (inherits(data, "integer64")) {
    unclass(simd_as_double(data))
  } else if (is.complex(data)) {
    # Its warning (imaginary parts discarded) names the user's call.
    .relay(as.double(data))
  } else {
    as.double(data)
  }
}

as.integer.simd_vec <- function(x, ...) {
  data <- .sv_strip(x)
  if (inherits(data, "integer64")) {
    unclass(simd_as_integer(data))
  } else if (is.double(data)) {
    # The kernel converts as base R does, and its warning (NAs out of
    # range) names the user's call.
    simd_as_integer(data)
  } else if (is.complex(data)) {
    # Its warning (imaginary parts discarded) names the user's call.
    .relay(as.integer(data))
  } else {
    as.integer(data)
  }
}

as.logical.simd_vec <- function(x, ...) {
  data <- .sv_strip(x)
  if (inherits(data, "integer64")) unclass(simd_as_logical(data)) else as.logical(data)
}

as.data.frame.simd_vec <- function(x, row.names = NULL, optional = FALSE, ...,
                                   nm = deparse1(substitute(x))) {
  as.data.frame(.sv_printable(x), row.names = row.names, optional = optional, ..., nm = nm)
}

# ---- Printing --------------------------------------------------------------

.sv_header <- function(x) {
  impl <- attr(x, "rsimd_impl", exact = TRUE)
  flag <- .sv_flag(x)
  sprintf(
    "impl=%s na_free=%s", if (is.null(impl)) "auto" else impl,
    if (is.null(flag)) "unknown" else as.character(flag)
  )
}

# The data in a form base R prints correctly: integer64 needs bit64, else
# it is shown as double.
.sv_printable <- function(x) {
  data <- .sv_strip(x)
  if (inherits(data, "integer64") && !requireNamespace("bit64", quietly = TRUE)) {
    data <- unclass(simd_as_double(data))
  }
  data
}

print.simd_vec <- function(x, ...) {
  # %.0f, as %d rejects the double length of a long vector.
  cat(sprintf("<simd_vec[%.0f] %s %s>\n", length(x), .sv_type(x), .sv_header(x)))
  if (length(x) > 0L) print(.sv_printable(x), ...)
  invisible(x)
}

format.simd_vec <- function(x, ...) format(.sv_printable(x), ...)

str.simd_vec <- function(object, ...) {
  n <- length(object)
  if (n == 0) {
    cat(sprintf("simd_vec [0] %s %s\n", .sv_type(object), .sv_header(object)))
    return(invisible())
  }
  shown <- format(.sv_printable(object[seq_len(min(n, 10L))]))
  cat(sprintf(
    "simd_vec [1:%.0f] %s %s %s%s\n", n, .sv_type(object), .sv_header(object),
    paste(shown, collapse = " "), if (n > 10L) " ..." else ""
  ))
  invisible()
}

# ---- Comparison ------------------------------------------------------------

# Data, pin and NA-free flag; not the flag's per-object token, which makes
# two separately built flagged simd_vecs differ for identical().
all.equal.simd_vec <- function(target, current, ...) {
  if (!is_simd_vec(current)) {
    return(paste0("target is simd_vec, current is ", data.class(current)))
  }
  show <- function(v, unset) if (is.null(v)) unset else as.character(v)
  msg <- NULL
  p1 <- attr(target, "rsimd_impl", exact = TRUE)
  p2 <- attr(current, "rsimd_impl", exact = TRUE)
  if (!identical(p1, p2)) {
    msg <- c(msg, paste0("pins differ: ", show(p1, "auto"), " vs ", show(p2, "auto")))
  }
  f1 <- .sv_flag(target)
  f2 <- .sv_flag(current)
  if (!identical(f1, f2)) {
    msg <- c(msg, paste0(
      "NA-free flags differ: ", show(f1, "unknown"), " vs ",
      show(f2, "unknown")
    ))
  }
  d1 <- .sv_strip(target)
  d2 <- .sv_strip(current)
  if (inherits(d1, "integer64") || inherits(d2, "integer64")) .need_bit64()
  data <- all.equal(d1, d2, ...)
  if (!isTRUE(data)) msg <- c(msg, data)
  if (is.null(msg)) TRUE else msg
}
