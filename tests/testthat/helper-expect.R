# Expectations that compare every tier with the none oracle and with base R.
#
# Doubles are compared with a bound that depends on the precision mode,
# since tiers add in different orders:
#   compensated  2^-50 relative;
#   pairwise     2^-50 relative, or (128 + log2 n) * 2^-52 * sum(abs(x))
#                absolute (each 128-element leaf is summed in fast mode);
#   fast         n * 8 * 2^-52 * sum(abs(x)) absolute.
# x is the input of the reduction and sum(abs(x)) is over its finite
# elements. Missing and infinite values must agree exactly: is.na(),
# is.nan() and the infinities.

eps <- 2^-52

# Batching. Reporting a testthat expectation costs about a millisecond, more
# than most of the kernel calls a sweep checks, so sweeps over many lengths,
# inputs and modes run inside batch_expectations(). There every helper in
# these files, and check_identical(), records a failed check instead of
# signalling an expectation, and one expectation at the end reports how
# many checks failed and the first 50. Nested batches join the outermost.
batch_state <- new.env(parent = emptyenv())

batch_expectations <- function(code) {
  if (!is.null(batch_state$problems)) {
    return(invisible(code))
  }
  batch_state$problems <- character(0)
  batch_state$checks <- 0L
  on.exit(batch_state$problems <- NULL)
  code
  problems <- batch_state$problems
  testthat::expect(length(problems) == 0L, sprintf(
    "%d of %d checks failed:\n%s", length(problems), batch_state$checks,
    paste(utils::head(problems, 50), collapse = "\n")
  ))
}

# testthat::expect(ok, message), or inside batch_expectations() a check
# recorded in the batch. `message` may be a function returning it, so that
# it is only built on failure.
simd_expect <- function(ok, message) {
  if (is.function(message)) message <- if (ok) "" else message()
  if (is.null(batch_state$problems)) {
    return(testthat::expect(ok, message))
  }
  batch_state$checks <- batch_state$checks + 1L
  if (!ok) batch_state$problems <- c(batch_state$problems, message)
  invisible(ok)
}

# expect_identical() and expect_true() that are batched inside
# batch_expectations(). num.eq = FALSE also tells 0 from -0 and NA from NaN.
check_identical <- function(object, expected, info = NULL, num.eq = TRUE) {
  simd_expect(identical(object, expected, num.eq = num.eq), function() {
    sprintf(
      "%s: %s not identical to %s", paste(info, collapse = " "),
      deparse1(object), deparse1(expected)
    )
  })
}

# expect_error(expr, regexp) that is batched inside batch_expectations().
check_error <- function(expr, regexp, info = NULL) {
  msg <- tryCatch({
    expr
    NULL
  }, error = conditionMessage)
  simd_expect(!is.null(msg) && grepl(regexp, msg), function() {
    sprintf(
      "%s: %s", paste(info, collapse = " "),
      if (is.null(msg)) "no error" else sprintf("error \"%s\" does not match %s", msg, regexp)
    )
  })
}

check_true <- function(object, info = NULL) {
  simd_expect(isTRUE(object), function() {
    sprintf("%s: %s is not TRUE", paste(info, collapse = " "), deparse1(object))
  })
}

# sum(abs(x)) over the finite elements of x, as a double.
abs_mass <- function(x) {
  if (!is.numeric(x) && !is.logical(x)) {
    return(0)
  }
  x <- as.double(x)
  sum(abs(x[is.finite(x)]))
}

# Absolute difference allowed between two tiers' results for input x.
mode_bound <- function(x, precision) {
  n <- max(length(x), 1)
  switch(precision,
    fast = n * 8 * eps * abs_mass(x),
    pairwise = (128 + log2(n)) * eps * abs_mass(x),
    compensated = 0
  )
}

# NULL if a and b (double vectors of equal length) agree within `bound`
# (absolute, per element) or within relative tolerance `rel`, else a
# description of the first difference. Missing values and infinities must
# match exactly; with nan = FALSE only is.na() is compared for them.
double_mismatch <- function(a, b, bound, rel, nan = TRUE) {
  if (!is.double(a) || !is.double(b)) {
    return(sprintf("types differ: %s vs %s", typeof(a), typeof(b)))
  }
  if (length(a) != length(b)) {
    return(sprintf("lengths differ: %d vs %d", length(a), length(b)))
  }
  if (!identical(is.na(a), is.na(b))) {
    return("is.na() differs")
  }
  if (nan && !identical(is.nan(a), is.nan(b))) {
    return("is.nan() differs")
  }
  ok <- !is.na(a)
  a <- a[ok]
  b <- b[ok]
  inf <- is.infinite(a) | is.infinite(b)
  if (!identical(a[inf], b[inf])) {
    return("infinite values differ")
  }
  a <- a[!inf]
  b <- b[!inf]
  diff <- abs(a - b)
  bad <- !(diff <= bound | diff <= rel * abs(b))
  if (any(bad)) {
    i <- which(bad)[1L]
    return(sprintf(
      "%.17g vs %.17g (difference %.3g, bound %.3g absolute or %.3g relative)",
      a[i], b[i], diff[i], bound, rel
    ))
  }
  NULL
}

# Evaluates f(...) with each available tier selected and expects every
# tier's result to equal the none tier's: identical for results that are
# not double, within the precision-mode bound above for doubles (or
# relative `tolerance` if given). The first argument in ... is the input x
# the bound is computed from. Evaluated in mode `precision`. One
# expectation, whose failure message names every tier that differs.
# Returns the results by tier, invisibly.
expect_simd_equal <- function(f, ..., tolerance = NULL, precision = simd_precision()) {
  old <- simd_precision(precision)
  on.exit(simd_precision(old))
  x <- ..1
  res <- with_each_tier(function() f(...))
  oracle <- res[["none"]]
  problems <- character(0)
  for (tier in setdiff(names(res), "none")) {
    got <- res[[tier]]
    problem <- if (!is.double(got) || !is.double(oracle)) {
      if (!identical(got, oracle)) {
        sprintf("%s vs %s", deparse1(got), deparse1(oracle))
      }
    } else if (!is.null(tolerance)) {
      double_mismatch(got, oracle, 0, tolerance)
    } else {
      double_mismatch(got, oracle, mode_bound(x, precision), 2^-50)
    }
    if (!is.null(problem)) {
      problems <- c(problems, sprintf(
        "tier %s differs from none (%s mode, n = %.0f): %s",
        tier, precision, length(x), problem
      ))
    }
  }
  simd_expect(length(problems) == 0L, paste(problems, collapse = "\n"))
  invisible(res)
}

# Evaluates simd_fn(...) with each available tier selected and expects it
# to equal base_fn(...): identical results unless they are double; doubles
# within the precision-mode bound plus base R's own rounding (its sum uses
# long double where available), and never tighter than 1e-12 relative, or
# within relative `tolerance` if given (plus n * eps, base R's own bound for
# a sum of positive terms, where it has no long double). When x holds NA and also NaN, or
# both Inf and -Inf (whose sum is NaN), only is.na() is compared, since base
# R's choice between NA and NaN depends on the order of the elements.
expect_simd_matches_base <- function(simd_fn, base_fn, ..., tolerance = NULL,
                                     precision = simd_precision()) {
  old <- simd_precision(precision)
  on.exit(simd_precision(old))
  x <- ..1
  expected <- base_fn(...)
  # Raw vectors have no is.nan() method, and no missing values.
  fp <- is.double(x) || is.complex(x)
  has_na <- if (fp) any(is.na(x) & !is.nan(x)) else anyNA(x)
  makes_nan <- fp &&
    (any(is.nan(x)) || (any(x == Inf, na.rm = TRUE) && any(x == -Inf, na.rm = TRUE)))
  nan <- !(has_na && makes_nan)
  base_eps <- if (has_wide_long_double()) 2^-63 else eps
  bound <- mode_bound(x, precision) + max(length(x), 1) * base_eps * abs_mass(x)
  res <- with_each_tier(function() simd_fn(...))
  problems <- character(0)
  for (tier in names(res)) {
    got <- res[[tier]]
    problem <- if (!is.double(got) || !is.double(expected)) {
      if (!identical(got, expected)) {
        sprintf("%s vs %s", deparse1(got), deparse1(expected))
      }
    } else if (!is.null(tolerance)) {
      base_tol <- if (has_wide_long_double()) 0 else max(length(x), 1) * eps
      double_mismatch(got, expected, 0, tolerance + base_tol, nan)
    } else {
      double_mismatch(got, expected, bound, max(1e-12, 2^-50), nan)
    }
    if (!is.null(problem)) {
      problems <- c(problems, sprintf(
        "tier %s differs from base R (%s mode, n = %.0f): %s",
        tier, precision, length(x), problem
      ))
    }
  }
  simd_expect(length(problems) == 0L, paste(problems, collapse = "\n"))
  invisible(res)
}

# Evaluates f(...) with each available tier selected and expects every
# tier's result to be bit-identical to the none tier's, including the sign
# of zero and NA versus NaN (identical(num.eq = FALSE)). For results that do
# not depend on the order of operations: min, max, which, any, all, counts.
# One expectation; returns the results by tier, invisibly.
expect_simd_identical <- function(f, ...) {
  x <- ..1
  res <- with_each_tier(function() f(...))
  oracle <- res[["none"]]
  problems <- character(0)
  for (tier in setdiff(names(res), "none")) {
    got <- res[[tier]]
    if (!identical(got, oracle, num.eq = FALSE)) {
      problems <- c(problems, sprintf(
        "tier %s differs from none (n = %.0f): %s vs %s",
        tier, length(x), deparse1(got), deparse1(oracle)
      ))
    }
  }
  simd_expect(length(problems) == 0L, paste(problems, collapse = "\n"))
  invisible(res)
}
