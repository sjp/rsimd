# The simd_vec wrapper: an atomic vector with class "simd_vec" (or
# c("simd_vec", "integer64")), an optional pinned implementation (attribute
# rsimd_impl) and a known-NA-free flag (attribute rsimd_na_free: TRUE,
# FALSE or absent for unknown).
#
# The simd_* functions themselves handle simd_vec operands on the C side
# (src/rvec.c): they run under the operands' pin, skip NA checks for
# operands flagged NA-free and wrap value-vector results again. The methods
# here map R's generics onto those functions, and do in R what has no
# kernel (base fallbacks, subsetting, c(), printing).

# ---- Construction and accessors -------------------------------------------

simd_vec <- function(x, impl = NULL, check_na = FALSE) {
  if (!is.logical(check_na) || length(check_na) != 1L || is.na(check_na)) {
    stop("'check_na' must be TRUE or FALSE", call. = FALSE)
  }
  if (is_simd_vec(x)) {
    if (!missing(impl)) attr(x, "rsimd_impl") <- .sv_check_impl(impl)
    if (check_na) attr(x, "rsimd_na_free") <- !simd_any_na(x)
    return(x)
  }
  data <- .sv_strip(x)
  impl <- .sv_check_impl(impl)
  na_free <- if (check_na) !simd_any_na(data) else NULL
  .sv_new(data, impl, na_free)
}

as_simd_vec <- function(x, ...) UseMethod("as_simd_vec")

as_simd_vec.default <- function(x, ...) simd_vec(x, ...)

as_simd_vec.simd_vec <- function(x, ...) x

is_simd_vec <- function(x) inherits(x, "simd_vec")

simd_impl <- function(x) {
  .sv_assert(x)
  attr(x, "rsimd_impl", exact = TRUE)
}

`simd_impl<-` <- function(x, value) {
  .sv_assert(x)
  attr(x, "rsimd_impl") <- .sv_check_impl(value)
  x
}

simd_na_free <- function(x) {
  .sv_assert(x)
  attr(x, "rsimd_na_free", exact = TRUE)
}

# ---- Internal helpers ------------------------------------------------------

.sv_assert <- function(x, arg = "x") {
  if (!is_simd_vec(x)) stop("'", arg, "' must be a simd_vec", call. = FALSE)
  invisible()
}

# A pin value checked: NULL, or the name of an available tier.
.sv_check_impl <- function(impl) {
  if (is.null(impl)) {
    return(NULL)
  }
  if (identical(impl, "auto")) {
    stop("'impl' must name an implementation, not \"auto\"; use NULL to follow the ",
      "global setting",
      call. = FALSE
    )
  }
  problem <- .impl_problem(impl)
  if (!is.null(problem)) stop(problem, call. = FALSE)
  impl
}

# The data of x without attributes, except an integer64 class. Errors for
# anything that is not a supported atomic vector.
.sv_strip <- function(x) {
  i64 <- inherits(x, "integer64")
  types <- c("double", "integer", "logical", "raw", "complex")
  ok <- i64 || is_simd_vec(x) || (is.atomic(x) && !is.object(x) && typeof(x) %in% types)
  if (!ok || is.factor(x)) {
    stop("'x' must be an atomic vector (double, integer, logical, raw, complex or ",
      "integer64), not ", class(x)[[1L]],
      call. = FALSE
    )
  }
  attributes(x) <- NULL
  if (i64) class(x) <- "integer64"
  x
}

.sv_new <- function(data, impl = NULL, na_free = NULL) {
  class(data) <- if (inherits(data, "integer64")) c("simd_vec", "integer64") else "simd_vec"
  attr(data, "rsimd_impl") <- impl
  attr(data, "rsimd_na_free") <- if (is.raw(data)) TRUE else na_free
  data
}

# The data of a simd_vec (integer64 keeps its class).
.sv_data <- function(x) .sv_strip(x)

# out (a bare result computed from x by base R) as a simd_vec like x: same
# pin, and the given NA-free flag. Returned unchanged when x is not a
# simd_vec or out already is one.
.sv_like <- function(out, x, na_free = attr(x, "rsimd_na_free", exact = TRUE)) {
  if (!is_simd_vec(x) || is_simd_vec(out)) {
    return(out)
  }
  .sv_new(out, attr(x, "rsimd_impl", exact = TRUE), na_free)
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
    if (!is.null(pin) && !identical(pin, p)) {
      stop("operands pinned to different implementations ('", pin, "' vs '", p,
        "'); unpin one with simd_impl(x) <- NULL",
        call. = FALSE
      )
    }
    if (!p %in% simd_available()) {
      stop("'x' is pinned to implementation '", p,
        "', which is not available on this machine; unpin it with simd_impl(x) <- NULL",
        call. = FALSE
      )
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
  data <- lapply(args, function(a) if (is_simd_vec(a)) .sv_data(a) else a)
  if (any(vapply(data, inherits, NA, what = "integer64"))) .sv_need_bit64()
  out <- do.call(f, data)
  if (wrap && ((is.atomic(out) && !is.object(out)) || inherits(out, "integer64"))) {
    out <- .sv_new(.sv_strip(out), pin)
  }
  out
}

.sv_need_bit64 <- function() {
  if (!requireNamespace("bit64", quietly = TRUE)) {
    stop("this operation on an integer64 simd_vec needs package 'bit64'", call. = FALSE)
  }
  invisible()
}

.sv_type <- function(x) if (inherits(x, "integer64")) "integer64" else typeof(x)

# Length rule of the binary operators (D24), for base fallbacks.
.sv_check_lengths <- function(e1, e2) {
  n1 <- length(e1)
  n2 <- length(e2)
  if (n1 != n2 && n1 != 1L && n2 != 1L) {
    stop("lengths of 'x' (", n1, ") and 'y' (", n2,
      ") must be equal or one of them must be 1",
      call. = FALSE
    )
  }
  invisible()
}

# ---- Group generics --------------------------------------------------------

# A simd_vec operand wins over another class's Ops method (R >= 4.3). bit64
# also always claims the operation, so a plain integer64 on the left of a
# simd_vec still dispatches to bit64.
chooseOpsMethod.simd_vec <- function(x, y, mx, my, cl, reverse) TRUE

Ops.simd_vec <- function(e1, e2) {
  unary <- nargs() == 1L
  gen <- .Generic
  if (unary) {
    return(switch(gen,
      "+" = .sv_uplus(e1),
      "-" = simd_neg(e1),
      "!" = if (is.raw(e1)) simd_bit_not(e1) else simd_not(e1),
      stop("invalid unary operator", call. = FALSE)
    ))
  }
  cplx <- is.complex(e1) || is.complex(e2)
  i64 <- inherits(e1, "integer64") || inherits(e2, "integer64")
  raw <- is.raw(e1) || is.raw(e2)
  switch(gen,
    "+" = simd_add(e1, e2),
    "-" = simd_sub(e1, e2),
    "*" = simd_mul(e1, e2),
    "/" = simd_div(e1, e2),
    "%%" = if (cplx) .sv_ops_fallback(gen, e1, e2) else simd_mod(e1, e2),
    "%/%" = if (cplx) .sv_ops_fallback(gen, e1, e2) else simd_idiv(e1, e2),
    "^" = if (cplx || i64) .sv_ops_fallback(gen, e1, e2) else simd_pow(e1, e2),
    "==" = if (cplx) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_eq(e1, e2),
    "!=" = if (cplx) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_ne(e1, e2),
    "<" = if (cplx) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_lt(e1, e2),
    ">" = if (cplx) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_gt(e1, e2),
    "<=" = if (cplx) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_le(e1, e2),
    ">=" = if (cplx) .sv_ops_fallback(gen, e1, e2, FALSE) else simd_ge(e1, e2),
    "&" = if (raw) simd_bit_and(e1, e2) else simd_and(e1, e2),
    "|" = if (raw) simd_bit_or(e1, e2) else simd_or(e1, e2),
    stop("operator '", gen, "' is not supported for simd_vec", call. = FALSE)
  )
}

# Unary +: base R makes a logical integer and rejects raw.
.sv_uplus <- function(x) {
  if (is.raw(x)) stop("invalid argument to unary operator", call. = FALSE)
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

Math.simd_vec <- function(x, ...) {
  gen <- .Generic
  type <- .sv_type(x)
  kernel <- type %in% c("double", "integer", "logical", "raw") ||
    (type == "integer64" && gen %in% .sv_math_i64)
  if (kernel && gen == "round") {
    return(simd_round(x, ...))
  }
  if (kernel && gen == "log") {
    return(simd_log(x, ...))
  }
  if (kernel && !is.null(f <- .sv_math_kernels[[gen]])) {
    return(f(x))
  }
  f <- get(gen, envir = baseenv())
  .sv_fallback(function(x) f(x, ...), list(x))
}

# Summary functions with a kernel, and the element types it rejects.
.sv_summary_kernels <- list(
  sum = list(simd_sum, character()),
  prod = list(simd_prod, c("complex", "integer64")),
  min = list(simd_min, "complex"),
  max = list(simd_max, "complex"),
  range = list(simd_range, "complex"),
  any = list(simd_any, "complex"),
  all = list(simd_all, "complex")
)

Summary.simd_vec <- function(..., na.rm = FALSE) {
  gen <- .Generic
  args <- list(...)
  x <- if (length(args) == 1L) args[[1L]] else .sv_combine(args)
  k <- .sv_summary_kernels[[gen]]
  if (is_simd_vec(x) && !(.sv_type(x) %in% k[[2L]])) {
    return(k[[1L]](x, na.rm = na.rm))
  }
  f <- get(gen, envir = baseenv())
  .sv_fallback(function(x) f(x, na.rm = na.rm), list(x), wrap = FALSE)
}

mean.simd_vec <- function(x, trim = 0, na.rm = FALSE, ...) {
  if (!identical(trim, 0) || .sv_type(x) %in% c("complex", "integer64")) {
    return(.sv_fallback(function(x) mean(x, trim = trim, na.rm = na.rm, ...), list(x),
      wrap = FALSE
    ))
  }
  simd_mean(x, na.rm = na.rm)
}

anyNA.simd_vec <- function(x, recursive = FALSE) {
  flag <- attr(x, "rsimd_na_free", exact = TRUE)
  if (isTRUE(flag)) {
    return(FALSE)
  }
  if (isFALSE(flag)) {
    return(TRUE)
  }
  simd_any_na(x)
}

# ---- Subsetting, combining, length -----------------------------------------

`[.simd_vec` <- function(x, i, ...) {
  if (missing(i)) {
    return(x)
  }
  if (...length() > 0L) stop("incorrect number of dimensions", call. = FALSE)
  data <- .sv_data(x)
  j <- seq_along(data)[i]
  out <- .sv_take(data, j)
  flag <- attr(x, "rsimd_na_free", exact = TRUE)
  na_free <- if (isTRUE(flag) && !anyNA(j)) TRUE else NULL
  .sv_new(out, attr(x, "rsimd_impl", exact = TRUE), na_free)
}

# data[j] for positions j (NA giving a missing element), integer64 too.
.sv_take <- function(data, j) {
  if (!inherits(data, "integer64")) {
    return(data[j])
  }
  out <- unclass(data)[j]
  if (anyNA(j)) out[is.na(j)] <- unclass(.sv_na_i64())
  class(out) <- "integer64"
  out
}

.sv_na_i64 <- function() simd_as_integer64(NA)

`[<-.simd_vec` <- function(x, i, value) {
  data <- .sv_data(x)
  if (is_simd_vec(value)) value <- .sv_data(value)
  if (inherits(data, "integer64")) {
    if (!inherits(value, "integer64")) value <- simd_as_integer64(value)
    data <- unclass(data)
    if (missing(i)) data[] <- unclass(value) else data[i] <- unclass(value)
    class(data) <- "integer64"
  } else if (inherits(value, "integer64")) {
    stop("cannot assign integer64 values into a ", typeof(data), " simd_vec", call. = FALSE)
  } else if (missing(i)) {
    data[] <- value
  } else {
    data[i] <- value
  }
  .sv_new(.sv_strip(data), attr(x, "rsimd_impl", exact = TRUE))
}

`[[.simd_vec` <- function(x, i, ...) {
  data <- .sv_data(x)
  if (inherits(data, "integer64")) .sv_take(data, seq_along(data)[[i]]) else data[[i]]
}

c.simd_vec <- function(...) .sv_combine(list(...))

# c() of simd_vec and plain vectors: base R's type promotion (integer64
# absorbing the integer, logical and double parts, as bit64's c() does),
# the parts' common pin, and NA-free when every part is known to be.
.sv_combine <- function(args) {
  args <- args[lengths(args) > 0L | vapply(args, is_simd_vec, NA)]
  pin <- .sv_resolve(args)
  flags <- lapply(args, function(a) if (is_simd_vec(a)) attr(a, "rsimd_na_free", exact = TRUE))
  data <- lapply(args, function(a) .sv_strip(a))
  if (any(vapply(data, inherits, NA, what = "integer64"))) {
    data <- lapply(data, function(d) {
      unclass(if (inherits(d, "integer64")) d else simd_as_integer64(d))
    })
    out <- unlist(data, use.names = FALSE)
    if (is.null(out)) out <- double()
    class(out) <- "integer64"
  } else {
    out <- unlist(data, use.names = FALSE)
    if (is.null(out)) out <- logical()
  }
  is_sv <- vapply(args, is_simd_vec, NA)
  na_free <- if (any(vapply(flags, isFALSE, NA))) {
    FALSE
  } else if (all(vapply(flags[is_sv], isTRUE, NA))) {
    if (all(is_sv) || !any(vapply(data[!is_sv], simd_any_na, NA))) TRUE
  }
  .sv_new(out, pin, na_free)
}

`length<-.simd_vec` <- function(x, value) {
  data <- .sv_data(x)
  n <- length(data)
  if (inherits(data, "integer64")) {
    raw <- unclass(data)
    length(raw) <- value
    if (value > n) raw[(n + 1L):value] <- unclass(.sv_na_i64())
    class(raw) <- "integer64"
    data <- raw
  } else {
    length(data) <- value
  }
  flag <- attr(x, "rsimd_na_free", exact = TRUE)
  na_free <- if (isTRUE(flag) && value <= n) TRUE else NULL
  .sv_new(data, attr(x, "rsimd_impl", exact = TRUE), na_free)
}

# ---- Unwrapping ------------------------------------------------------------

as.vector.simd_vec <- function(x, mode = "any") {
  data <- .sv_data(x)
  if (inherits(data, "integer64")) data <- simd_as_double(data)
  as.vector(data, mode)
}

as.double.simd_vec <- function(x, ...) {
  data <- .sv_data(x)
  if (inherits(data, "integer64")) unclass(simd_as_double(data)) else as.double(data)
}

as.integer.simd_vec <- function(x, ...) {
  data <- .sv_data(x)
  if (inherits(data, "integer64")) unclass(simd_as_integer(data)) else as.integer(data)
}

as.logical.simd_vec <- function(x, ...) {
  data <- .sv_data(x)
  if (inherits(data, "integer64")) unclass(simd_as_logical(data)) else as.logical(data)
}

as.data.frame.simd_vec <- function(x, row.names = NULL, optional = FALSE, ...,
                                   nm = deparse1(substitute(x))) {
  as.data.frame(.sv_data(x), row.names = row.names, optional = optional, ..., nm = nm)
}

# ---- Printing --------------------------------------------------------------

.sv_header <- function(x) {
  impl <- attr(x, "rsimd_impl", exact = TRUE)
  flag <- attr(x, "rsimd_na_free", exact = TRUE)
  sprintf(
    "impl=%s na_free=%s", if (is.null(impl)) "auto" else impl,
    if (is.null(flag)) "unknown" else as.character(flag)
  )
}

# The data in a form base R prints correctly: integer64 needs bit64, else
# it is shown as double.
.sv_printable <- function(x) {
  data <- .sv_data(x)
  if (inherits(data, "integer64") && !requireNamespace("bit64", quietly = TRUE)) {
    data <- unclass(simd_as_double(data))
  }
  data
}

print.simd_vec <- function(x, ...) {
  cat(sprintf("<simd_vec[%d] %s %s>\n", length(x), .sv_type(x), .sv_header(x)))
  if (length(x) > 0L) print(.sv_printable(x), ...)
  invisible(x)
}

format.simd_vec <- function(x, ...) format(.sv_printable(x), ...)

str.simd_vec <- function(object, ...) {
  n <- length(object)
  shown <- format(.sv_printable(object[seq_len(min(n, 10L))]))
  cat(sprintf(
    "simd_vec [1:%d] %s %s %s%s\n", n, .sv_type(object), .sv_header(object),
    paste(shown, collapse = " "), if (n > 10L) " ..." else ""
  ))
  invisible()
}
