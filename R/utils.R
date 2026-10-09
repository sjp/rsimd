# Operands of a binary operation converted to their common element type
# (rsimd_promote() in src/rvec.c), as list(x, y). Only an operand whose type
# differs from the common type is converted, so same-type operands are
# passed through without a copy; two logicals stay logical (their storage
# is integer). Converting integer64 to double warns, naming the caller.
.promote_pair <- function(x, y, call = .user_call()) {
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
    .stop("internal error: cannot convert to ", to)
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
    .stop(fun, "() does not support '", arg, "' of type ", type, if (type %in% later) " yet")
  }
  invisible()
}

# Errors unless x is data rsimd takes: an atomic vector of a supported type
# without a class, an integer64 or a simd_vec (NULL is an error). The message
# is the one rsimd_check_atomic() in src/rvec.c gives, naming the argument
# `arg`.
.check_data <- function(x, arg = "x") {
  types <- c("double", "integer", "logical", "raw", "complex")
  ok <- inherits(x, c("integer64", "simd_vec")) ||
    (is.atomic(x) && !is.object(x) && typeof(x) %in% types)
  if (!ok) {
    if (is.null(x)) .stop("'", arg, "' is NULL")
    what <- if (is.object(x)) class(x)[[1L]] else if (is.list(x)) "list" else typeof(x)
    .stop("'", arg, "' must be an atomic vector (double, integer, logical, raw, complex or ",
      "integer64), not ", what)
  }
  invisible()
}

# Loads bit64, so that its S3 methods are registered before rsimd hands out
# an integer64 object (without them it sorts, prints and compares as the
# doubles its bits make); errors if it is not installed. The C side does
# the same for the integer64 results of the simd_* functions.
.need_bit64 <- function() {
  if (!isNamespaceLoaded("bit64") && !requireNamespace("bit64", quietly = TRUE)) {
    .stop("integer64 results need package 'bit64'; install it with install.packages(\"bit64\")")
  }
  invisible()
}

# The error for arguments that fall into the `...` of a function taking
# `what` operands, whose na.rm (and na_check) must be named: simd_sum(x, y)
# would otherwise read y as na.rm. Names the user's call.
.dots_error <- function(fun, what = "one vector; combine several with c()",
                        named = "na.rm and na_check") {
  .stop(fun, "() takes ", what, ", and ", named, " must be named")
}

# An error made of the pasted arguments, as stop() makes it, that names the
# user's call (.user_call()) rather than the rsimd helper raising it.
.stop <- function(...) {
  stop(simpleError(.makeMessage(...), .user_call()))
}

# The call the user made into rsimd: the innermost call of an rsimd
# function from outside the package (an exported function, or the simd_vec
# method a base generic dispatched to, as the generic's call). Walks the
# stack, so only for errors and warnings.
.user_call <- function(last = sys.nframe() - 1L) {
  i <- .user_frame(last)
  if (i > 0L) .generic_call(sys.call(i), sys.frame(i)) else NULL
}

# The value of expr, base R code run for the user's call into rsimd, with
# any error or warning it raises given that call (.user_call()) in place
# of the internal one. The handlers cost about a microsecond, so hot paths
# call this only when expr can raise a condition.
.relay <- function(expr) {
  # The handlers run in frames above the signal's, so the user's call is
  # sought below this frame.
  n <- sys.nframe()
  withCallingHandlers(expr,
    error = function(e) {
      e$call <- .user_call(n)
      stop(e)
    },
    warning = function(w) {
      w$call <- .user_call(n)
      warning(w)
      invokeRestart("muffleWarning")
    }
  )
}

# `call`, of the function whose frame is `frame`, as the user wrote it: a
# method that a generic dispatched to gets the generic's name (v + 1, not
# Ops.simd_vec(v, 1)), and a replacement method the assignment
# (`*tmp*`[i] <- value, `*tmp*` being R's name for the object changed). A
# group generic such as Summary passes the values of the arguments, not
# their expressions: a simd_vec among them is shown as its data.
.generic_call <- function(call, frame) {
  gen <- get0(".Generic", envir = frame, inherits = FALSE)
  head <- call[[1L]]
  if (!is.character(gen) || length(gen) != 1L || !is.symbol(head) ||
    identical(as.character(head), gen)) {
    return(call)
  }
  # Plain R, not simd_unwrap(): this may run inside a .Call that is raising
  # an error (rsimd_error()).
  args <- lapply(as.list(call)[-1L], function(a) {
    if (!is_simd_vec(a)) return(a)
    attr(a, "rsimd_impl") <- attr(a, "rsimd_na_token") <- NULL
    oldClass(a) <- setdiff(oldClass(a), "simd_vec")
    a
  })
  n <- length(args)
  if (endsWith(gen, "<-") && n >= 2L && identical(names(args)[[n]], "value")) {
    target <- as.call(c(as.name(substr(gen, 1L, nchar(gen) - 2L)), args[-n]))
    return(call("<-", target, args[[n]]))
  }
  as.call(c(as.name(gen), args))
}

# The frame number of the call .user_call() returns (0 for none), as seen
# from the caller of .user_frame(), among frames 1 to last.
.user_frame <- function(last = sys.nframe() - 1L) {
  ns <- topenv()
  parents <- sys.parents()
  frames <- sys.frames()
  for (i in rev(seq_len(min(last, length(parents))))) {
    if (!identical(topenv(environment(sys.function(i))), ns)) next
    p <- parents[[i]]
    if (!identical(topenv(if (p == 0L) globalenv() else frames[[p]]), ns)) return(i)
  }
  0L
}

# The error the C side raises (rsimd_type_error(), rsimd_mix_error()) for
# an operand of `type` that the op does not take (with one of type
# `other`, for a pair that it does not take together). In the user's call
# of an exported simd_* function it is rsimd's message, naming the first
# argument of that function holding a value of `type` (preferring `arg`,
# the C side's name for it); elsewhere (the simd_vec methods of base
# generics) it is base R's `msg`, naming the call that made the .Call(), as
# an error raised in C would.
.type_error <- function(type, arg, msg, other = NULL) {
  i <- .user_frame()
  call <- if (i > 0L) sys.call(i)
  fun <- if (i > 0L) .simd_fun_name(call, sys.function(i))
  if (is.null(fun)) stop(simpleError(msg, if (i > 0L) .generic_call(call, sys.frame(i))))
  if (!is.null(other)) {
    stop(simpleError(paste0(fun, "() cannot combine ", type, " and ", other, " operands"), call))
  }
  frame <- sys.frame(i)
  has_type <- function(a) {
    exists(a, envir = frame, inherits = FALSE) &&
      identical(tryCatch(.sv_type(get(a, envir = frame)), error = function(e) NULL), type)
  }
  formal <- setdiff(names(formals(sys.function(i))), "...")
  found <- if (!is.null(arg) && arg %in% formal && has_type(arg)) arg else Find(has_type, formal)
  if (!is.null(found)) arg <- found
  if (is.null(arg)) arg <- formal[[1L]]
  stop(simpleError(paste0(fun, "() does not support '", arg, "' of type ", type), call))
}

# The name of the exported simd_* function `f`, called by `call`: the name
# in the call (plain or rsimd::name) when it is that function, else the
# first export that is (do.call(simd_var, ...)); NULL if f is none.
.simd_fun_name <- function(call, f) {
  ns <- topenv()
  head <- call[[1L]]
  if (is.call(head) && length(head) == 3L && as.character(head[[1L]]) %in% c("::", ":::")) {
    head <- head[[3L]]
  }
  names <- getNamespaceExports(ns)
  names <- names[startsWith(names, "simd_")]
  if (is.symbol(head) && as.character(head) %in% names) names <- c(as.character(head), names)
  for (name in names) {
    if (identical(get(name, envir = ns), f)) return(name)
  }
  NULL
}
