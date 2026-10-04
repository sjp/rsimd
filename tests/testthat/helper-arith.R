# Helpers for the elementwise arithmetic tests.

# Are a and b the same values: same type and length, missing values of the
# same kind (NA or NaN) in the same places, and every other element
# identical bit for bit (so 0 and -0 differ)? Attributes are compared too.
same_values <- function(a, b) {
  if (!identical(typeof(a), typeof(b)) || length(a) != length(b)) return(FALSE)
  if (!identical(attributes(a), attributes(b))) return(FALSE)
  if (!is.double(a)) return(identical(a, b))
  ok <- !is.na(a)
  identical(is.na(a), is.na(b)) && identical(is.nan(a), is.nan(b)) &&
    identical(a[ok], b[ok], num.eq = FALSE)
}

# A short description of the first element where a and b differ.
first_difference <- function(a, b) {
  if (!identical(typeof(a), typeof(b))) return(sprintf("types %s vs %s", typeof(a), typeof(b)))
  if (length(a) != length(b)) return(sprintf("lengths %d vs %d", length(a), length(b)))
  if (!identical(attributes(a), attributes(b))) return("attributes differ")
  kind <- function(v) {
    if (is.double(v)) {
      ifelse(is.na(v) & !is.nan(v), "NA", ifelse(is.nan(v), "NaN", sprintf("%a", v)))
    } else {
      ifelse(is.na(v), "NA", as.character(v))
    }
  }
  ka <- kind(a)
  kb <- kind(b)
  i <- which(ka != kb)[1L]
  sprintf("element %d: %s vs %s", i, ka[i], kb[i])
}

# Expects f(...) to give `expected` (same_values) with every available tier
# selected. One expectation naming each tier that differs.
expect_tiers_give <- function(expected, f, ...) {
  res <- with_each_tier(function() f(...))
  problems <- character(0)
  for (tier in names(res)) {
    if (!same_values(res[[tier]], expected)) {
      problems <- c(problems, sprintf(
        "tier %s: %s", tier, first_difference(res[[tier]], expected)
      ))
    }
  }
  testthat::expect(length(problems) == 0L, paste(problems, collapse = "\n"))
  invisible(res)
}

# The messages of the warnings that evaluating `expr` raises, muffled.
warnings_of <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  msgs
}

# Expects f(...) to raise exactly the warnings `msgs` (possibly none) on
# every tier, and returns the none tier's value.
expect_tier_warnings <- function(msgs, f, ...) {
  res <- with_each_tier(function() {
    value <- NULL
    w <- warnings_of(value <- f(...))
    list(value = value, warnings = w)
  })
  for (tier in names(res)) {
    expect_identical(res[[tier]]$warnings, msgs, info = paste("tier", tier))
  }
  invisible(res[["none"]]$value)
}

# Every pair of elements of a and b, as list(x, y) of equal length.
all_pairs <- function(a, b = a) {
  list(x = rep(a, each = length(b)), y = rep(b, times = length(a)))
}

# Is long double IEEE quadruple precision (arm64 Linux)? Then base R's
# %/% and %% are exact, as simd_idiv() and simd_mod() are.
has_quad_long_double <- function() {
  isTRUE(capabilities("long.double")) && isTRUE(.Machine$longdouble.digits >= 113)
}

# NA_real_ wherever x or y is NA: the rule of the elementwise double ops
# for missing values, applied to a base R result (base R gives NaN for
# NaN + NA on some CPUs).
na_merged <- function(r, ...) {
  for (v in list(...)) r[is.na(v) & !is.nan(v)] <- NA_real_
  r
}

# Integer and logical results of a two's complement operation on int32
# values, from exact doubles: NA where an operand is NA, and a wrapped
# result of -2^31 is NA too.
wrap_i32 <- function(v) {
  v <- v %% 2^32
  v[!is.na(v) & v >= 2^31] <- v[!is.na(v) & v >= 2^31] - 2^32
  as.integer(ifelse(v == -2^31, NA, v))
}
