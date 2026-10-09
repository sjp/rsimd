# The simd_vec wrapper: construction, implementation pins, the NA-free
# flag, the group generics and the other methods (R/simd_vec.R), and the
# C side that pins, flags and wraps the results of every simd_ function
# (src/rvec.c).

# The data of a simd_vec as a plain vector (integer64 keeps its class).
# .subset() copies no attributes: attributes(x) <- NULL on a shared x
# would copy them, sharing the NA-free token and voiding x's flag.
bare <- function(x) {
  i64 <- inherits(x, "integer64")
  x <- .subset(x, if (length(x)) TRUE else 0L)
  if (i64) class(x) <- "integer64"
  x
}

# A tier that is a known id but not available here, or NULL.
unavailable_tier <- function() {
  t <- setdiff(simd_tiers(), c("auto", simd_available()))
  if (length(t)) t[[1L]] else NULL
}

# x must be a simd_vec holding `value`, pinned to `impl`, with NA-free flag
# `na_free`.
expect_sv <- function(x, value, impl = NULL, na_free = NULL) {
  check_true(is_simd_vec(x), info = "is_simd_vec")
  check_identical(bare(x), value, info = c("value", impl))
  check_identical(simd_impl(x), impl, info = c("simd_impl", impl))
  check_identical(simd_na_free(x), na_free, info = c("simd_na_free", impl))
  # Last: attributes() marks the token shared, which voids the flag.
  check_identical(
    sort(names(attributes(x))),
    sort(c(
      "class", if (!is.null(impl)) "rsimd_impl",
      if (!is.null(na_free)) "rsimd_na_token"
    )),
    info = c("attributes", impl)
  )
}

# ---- Construction ----------------------------------------------------------

test_that("simd_vec() wraps the supported types and drops other attributes", {
  expect_sv(simd_vec(c(a = 1, b = 2)), c(1, 2))
  expect_sv(simd_vec(matrix(1:4, 2)), 1:4)
  expect_sv(simd_vec(c(TRUE, NA)), c(TRUE, NA))
  expect_sv(simd_vec(1i), 1i)
  expect_sv(simd_vec(numeric()), numeric())
  # Raw has no missing values.
  expect_sv(simd_vec(as.raw(1:2)), as.raw(1:2), na_free = TRUE)
  expect_identical(class(simd_vec(1)), "simd_vec")

  if (has_bit64()) {
    i64 <- simd_vec(simd_as_integer64(c(1, NA)))
    expect_identical(class(i64), c("simd_vec", "integer64"))
    expect_identical(bare(i64), simd_as_integer64(c(1, NA)))
  }
})

test_that("simd_vec() rejects unsupported input", {
  msg <- "'x' must be an atomic vector"
  expect_error(simd_vec("a"), msg)
  expect_error(simd_vec(factor("a")), msg)
  expect_error(simd_vec(list(1)), msg)
  expect_error(simd_vec(data.frame(a = 1)), msg)
  expect_error(simd_vec(NULL), msg)
  expect_error(simd_vec(Sys.Date()), msg)
  expect_error(simd_vec(1, check_na = NA), "'check_na' must be TRUE or FALSE")
})

test_that("check_na records the NA-free flag", {
  expect_identical(simd_na_free(simd_vec(c(1, 2), check_na = TRUE)), TRUE)
  expect_identical(simd_na_free(simd_vec(c(1, NA), check_na = TRUE)), FALSE)
  expect_identical(simd_na_free(simd_vec(c(1, NaN), check_na = TRUE)), FALSE)
  expect_identical(simd_na_free(simd_vec(c(1L, NA), check_na = TRUE)), FALSE)
  expect_identical(simd_na_free(simd_vec(c(1, NA))), NULL)
  skip_if_not(has_bit64(), "bit64 is not installed")
  expect_identical(simd_na_free(simd_vec(simd_as_integer64(c(1, NA)), check_na = TRUE)), FALSE)
})

test_that("pins are validated", {
  t <- simd_available()[[1L]]
  expect_identical(simd_impl(simd_vec(1, impl = t)), t)
  expect_identical(simd_impl(simd_vec(1)), NULL)
  expect_error(simd_vec(1, impl = "auto"), "not \"auto\"")
  expect_error(simd_vec(1, impl = "bogus"), "unknown implementation 'bogus'")
  expect_error(simd_vec(1, impl = 1), "'impl' must be a single string")
  t_na <- unavailable_tier()
  if (!is.null(t_na)) {
    expect_error(simd_vec(1, impl = t_na), "not available on this machine")
  }
})

test_that("simd_vec() of a simd_vec replaces only what is asked", {
  x <- simd_vec(c(1, 2), impl = "none", check_na = TRUE)
  expect_identical(simd_vec(x), x)
  expect_sv(simd_vec(x, impl = NULL), c(1, 2), NULL, TRUE)
  y <- simd_vec(simd_vec(c(1, NA)), check_na = TRUE)
  expect_identical(simd_na_free(y), FALSE)
})

test_that("accessors", {
  x <- simd_vec(1:3)
  expect_true(is_simd_vec(x))
  expect_false(is_simd_vec(1:3))
  expect_identical(as_simd_vec(x), x)
  expect_identical(as_simd_vec(1:3, impl = "none"), simd_vec(1:3, impl = "none"))
  simd_impl(x) <- "none"
  expect_identical(simd_impl(x), "none")
  simd_impl(x) <- NULL
  expect_identical(simd_impl(x), NULL)
  expect_error(simd_impl(x) <- "auto", "not \"auto\"")
  expect_error(simd_impl(1), "'x' must be a simd_vec")
  expect_error(simd_na_free(1), "'x' must be a simd_vec")
})

# ---- Implementation pins ---------------------------------------------------

test_that("a pinned operand selects its tier over the global selection", {
  for (t in tiers_to_test()) {
    x <- simd_vec(1, impl = t)
    for (g in tiers_to_test()) {
      expect_identical(simd_with_impl(g, .debug_active(x)), t)
      expect_identical(simd_with_impl(g, .debug_active(x, 2)), t)
      expect_identical(simd_with_impl(g, .debug_active(3, x)), t)
      # Unpinned simd_vec operands follow the global selection.
      expect_identical(simd_with_impl(g, .debug_active(simd_vec(1))), g)
      expect_identical(simd_with_impl(g, .debug_active(simd_vec(1), x)), t)
    }
    expect_identical(.debug_active(x, simd_vec(2, impl = t)), t)
  }
})

test_that("different pins are an error, and the next call is unaffected", {
  tiers <- simd_available()
  skip_if(length(tiers) < 2L, "only one tier available")
  a <- simd_vec(c(1, 2), impl = tiers[[1L]])
  b <- simd_vec(c(1, 2), impl = tiers[[2L]])
  msg <- sprintf(
    paste(
      "operands pinned to different implementations ('%s' vs '%s');",
      "unpin one with simd_impl(x) <- NULL"
    ),
    tiers[[1L]], tiers[[2L]]
  )
  expect_error(a + b, msg, fixed = TRUE)
  expect_error(simd_add(a, b), msg, fixed = TRUE)
  expect_error(simd_dot(a, b), msg, fixed = TRUE)
  expect_error(c(a, b), msg, fixed = TRUE)
  expect_error(sum(a, b), msg, fixed = TRUE)
  expect_error(simd_hamming(a, b), msg, fixed = TRUE)
  expect_error(simd_hamming(a, simd_vec(c(1i, 2i), impl = tiers[[2L]])), msg, fixed = TRUE)
  # The failed call switched the table to the first pin; the next unpinned
  # call runs on the selection again.
  simd_with_impl(tiers[[2L]], {
    expect_identical(.debug_active(1), tiers[[2L]])
  })
  simd_impl(b) <- NULL
  expect_sv(a + b, c(2, 4), tiers[[1L]])
})

test_that("[<- and [[<- reject a value pinned to another tier, as c() does", {
  tiers <- simd_available()
  skip_if(length(tiers) < 2L, "only one tier available")
  a <- simd_vec(c(1, 2), impl = tiers[[1L]])
  b <- simd_vec(5, impl = tiers[[2L]])
  msg <- sprintf("operands pinned to different implementations ('%s' vs '%s')", tiers[[1L]], tiers[[2L]])
  expect_error(a[2] <- b, msg, fixed = TRUE)
  expect_error(a[[2]] <- b, msg, fixed = TRUE)
  i <- simd_vec(c(1L, 2L), impl = tiers[[1L]])
  expect_error(i[2] <- simd_vec(5L, impl = tiers[[2L]]), msg, fixed = TRUE)
  expect_sv(a, c(1, 2), tiers[[1L]])
  # The same pin, or one side unpinned, keeps x's pin (or lack of one).
  a[2] <- simd_vec(5, impl = tiers[[1L]])
  expect_sv(a, c(1, 5), tiers[[1L]])
  a[[1]] <- simd_vec(4)
  expect_sv(a, c(4, 5), tiers[[1L]])
  u <- simd_vec(c(1, 2))
  u[2] <- b
  expect_sv(u, c(1, 5), NULL)
})

test_that("[<- does not check that the pins are available", {
  t_na <- unavailable_tier()
  skip_if(is.null(t_na), "every tier is available")
  x <- structure(c(1, 2), class = "simd_vec", rsimd_impl = t_na)
  x[2] <- structure(5, class = "simd_vec", rsimd_impl = t_na)
  expect_identical(attr(x, "rsimd_impl"), t_na)
  expect_identical(.subset(x, 1:2), c(1, 5))
})

test_that("a nested call from a warning handler leaves the outer call intact", {
  # Warnings are issued when the call has finished, so a handler that calls
  # rsimd cannot overwrite the outer call's simd_vec state or pin.
  nested <- function(expr, inner) {
    withCallingHandlers(expr, warning = function(w) {
      inner()
      invokeRestart("muffleWarning")
    })
  }
  x <- simd_vec(c(.Machine$integer.max, 1L))
  expect_sv(nested(simd_add(x, 1L), function() simd_add(1, 2)), c(NA, 2L))
  tiers <- simd_available()
  for (t in tiers) {
    other <- tiers[[if (t == tiers[[1L]]) length(tiers) else 1L]]
    xp <- simd_vec(c(.Machine$integer.max, 1L), impl = t)
    inner <- function() {
      simd_sum(1:3)
      simd_add(simd_vec(1, impl = other), 2)
      expect_identical(.debug_active(1), simd_current()[[1L]])
    }
    expect_sv(nested(simd_add(xp, 1L), inner), c(NA, 2L), t)
    expect_sv(nested(xp + 1L, inner), c(NA, 2L), t)
    expect_sv(nested(simd_sqrt(simd_vec(c(-1, 4), impl = t)), inner), c(NaN, 2), t)
    # Every warning of the call reaches the handler, in order.
    msgs <- character(0)
    withCallingHandlers(
      simd_as_integer(simd_vec(c(1e10, 2), impl = t)),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        inner()
        invokeRestart("muffleWarning")
      }
    )
    expect_identical(msgs, "NAs introduced by coercion to integer range")
  }
})

test_that("a call that failed with a pinned operand does not leave its tier active", {
  tiers <- simd_available()
  skip_if(length(tiers) < 2L, "only one tier available")
  bad <- simd_vec(c(1, 2), impl = tiers[[2L]])
  expect_error(simd_add(bad, c(1, 2, 3)), "lengths")
  expect_identical(simd_probe_slots()[["tier_name"]], simd_current()[[1L]])
})

test_that("a pin to an unavailable or invalid tier is an error when used", {
  t_na <- unavailable_tier()
  skip_if(is.null(t_na), "every tier is available")
  x <- structure(c(1, 2), class = "simd_vec", rsimd_impl = t_na)
  msg <- sprintf("'x' is pinned to implementation '%s', which is not available", t_na)
  expect_error(x + 1, msg, fixed = TRUE)
  expect_error(sum(x), msg, fixed = TRUE)
  expect_error(exp(x), msg, fixed = TRUE)
  expect_error(c(x, 1), msg, fixed = TRUE)
  expect_error(simd_is_even(x), msg, fixed = TRUE)
  expect_error(simd_hamming(x, 1), msg, fixed = TRUE)
  expect_error(simd_hamming(1i, x), sub("'x'", "'y'", msg), fixed = TRUE)
  expect_error(
    simd_hamming_bits(structure(1L, class = "simd_vec", rsimd_impl = t_na), 1L), msg,
    fixed = TRUE
  )
  simd_impl(x) <- NULL
  expect_identical(bare(x + 1), c(2, 3))
})

test_that("an invalid rsimd_impl attribute is an error", {
  x <- structure(c(1, 2), class = "simd_vec", rsimd_impl = "bogus")
  expect_error(x + 1, "'x' has an invalid 'rsimd_impl' attribute")
  y <- structure(c(1, 2), class = "simd_vec", rsimd_impl = 3)
  expect_error(simd_sum(y), "invalid 'rsimd_impl' attribute")
})

test_that("results carry the resolved pin", {
  t <- simd_available()[[1L]]
  x <- simd_vec(c(1, 2), impl = t)
  expect_identical(simd_impl(x + 1), t)
  expect_identical(simd_impl(1 + x), t)
  expect_identical(simd_impl(x + simd_vec(1)), t)
  expect_identical(simd_impl(simd_vec(1) + simd_vec(2)), NULL)
  expect_identical(simd_impl(exp(x)), t)
  expect_identical(simd_impl(x[1]), t)
  expect_identical(simd_impl(c(x, 3)), t)
  expect_identical(simd_impl(c(simd_vec(1), x)), t)
  expect_identical(simd_impl(gamma(x)), t)
})

# ---- simd_ functions given a simd_vec ---------------------------------------

test_that("value results of simd_ functions are simd_vecs, others are plain", {
  x <- simd_vec(c(0.5, -1, 2), impl = "none")
  expect_sv(simd_add(x, 1), c(1.5, 0, 3), "none")
  expect_sv(simd_sigmoid(x), simd_sigmoid(c(0.5, -1, 2)), "none")
  expect_sv(simd_softmax(x), simd_softmax(c(0.5, -1, 2)), "none")
  expect_sv(simd_as_integer(x, "truncating"), c(0L, -1L, 2L), "none")
  expect_sv(simd_cumsum(x), c(0.5, -0.5, 1.5), "none")
  sc <- simd_sincos(x)
  expect_sv(sc$sin, sin(c(0.5, -1, 2)), "none")
  expect_sv(sc$cos, cos(c(0.5, -1, 2)), "none")
  expect_sv(simd_conj(simd_vec(c(1 + 2i, NA))), c(1 - 2i, NA))
  expect_sv(simd_popcount(simd_vec(c(7L, NA))), c(3L, NA))
  # A length-1 operand is checked for NA directly.
  expect_sv(simd_popcount(simd_vec(7L)), 3L, na_free = TRUE)

  expect_identical(simd_sum(x), 1.5)
  expect_identical(simd_range(x), c(-1, 2))
  expect_identical(simd_which_max(x), 3L)
  expect_identical(simd_is_na(x), c(FALSE, FALSE, FALSE))
  expect_identical(simd_lt(x, 1), c(TRUE, TRUE, FALSE))
  expect_identical(simd_and(simd_vec(c(TRUE, NA)), TRUE), c(TRUE, NA))
  expect_identical(simd_var(x), var(c(0.5, -1, 2)))
  expect_identical(simd_which_na(simd_vec(c(1, NA))), 2L)
  expect_identical(simd_is_whole(x), c(FALSE, TRUE, TRUE))
  expect_identical(simd_is_pow2_any(x), TRUE)
  expect_identical(simd_hamming(x, c(0.5, 1, 2)), 1)
  expect_identical(simd_hamming(x, 1i), 3)
  expect_identical(simd_hamming(simd_vec(as.raw(1:3)), 2L), 2)
  expect_identical(simd_hamming_bits(simd_vec(c(1L, 2L)), 3L), 2)
})

test_that("the input is not modified", {
  x <- simd_vec(c(1, 2), check_na = TRUE)
  before <- x
  y <- simd_floor(x)
  simd_impl(y) <- "none"
  z <- simd_as_double(x)
  class(z) <- NULL
  expect_identical(x, before)
})

# ---- Ops -------------------------------------------------------------------

test_that("arithmetic operators equal the simd_ functions on every tier", {
  batch_expectations({
    x <- c(1.5, -2, 0, NA, 7.25)
    xi <- c(3L, -7L, 0L, NA, 100L)
    xl <- c(TRUE, FALSE, NA, TRUE, TRUE)
    ops <- list(
      `+` = simd_add, `-` = simd_sub, `*` = simd_mul, `/` = simd_div,
      `^` = simd_pow, `%%` = simd_mod, `%/%` = simd_idiv
    )
    for (t in tiers_to_test()) {
      for (op in names(ops)) {
        f <- get(op)
        for (a in list(x, xi, xl)) {
          for (b in list(a, rev(a), 2L, 0.5)) {
            if (op == "^" && is.logical(a) && is.logical(b)) next
            want <- simd_with_impl(t, ops[[op]](a, b))
            expect_sv(f(simd_vec(a, impl = t), b), want, t)
            expect_sv(f(a, simd_vec(b, impl = t)), want, t)
            expect_sv(f(simd_vec(a, impl = t), simd_vec(b)), want, t)
          }
        }
      }
      expect_sv(-simd_vec(x, impl = t), simd_with_impl(t, simd_neg(x)), t)
      expect_sv(-simd_vec(xi, impl = t), -xi, t)
      expect_sv(-simd_vec(xl, impl = t), -xl, t)
    }
  })
})

test_that("operators follow base R's types and warnings", {
  expect_identical(bare(simd_vec(1:3) + 1L), 2:4)
  expect_identical(bare(simd_vec(1:3) / 2L), c(0.5, 1, 1.5))
  expect_identical(bare(simd_vec(TRUE) + TRUE), 2L)
  expect_identical(bare(simd_vec(2L)^2L), 4)
  expect_warning(
    expect_identical(bare(simd_vec(.Machine$integer.max) + 1L), NA_integer_),
    "NAs produced by integer overflow"
  )
  expect_error(simd_vec(as.raw(1)) + as.raw(1), "non-numeric argument to binary operator")
  expect_error(simd_vec(1:3) + 1:2, "lengths of 'x' \\(3\\) and 'y' \\(2\\)")
})

test_that("unary + follows base R", {
  expect_sv(+simd_vec(c(1, 2)), c(1, 2))
  expect_sv(+simd_vec(TRUE, check_na = TRUE), 1L, na_free = TRUE)
  expect_error(+simd_vec(as.raw(1)), "invalid argument to unary operator")
})

test_that("comparisons return plain logical vectors", {
  x <- c(1.5, -2, NA, 4)
  ops <- list(
    `==` = simd_eq, `!=` = simd_ne, `<` = simd_lt, `>` = simd_gt,
    `<=` = simd_le, `>=` = simd_ge
  )
  for (t in tiers_to_test()) {
    for (op in names(ops)) {
      f <- get(op)
      want <- simd_with_impl(t, ops[[op]](x, 1.5))
      expect_identical(f(simd_vec(x, impl = t), 1.5), want)
      expect_identical(f(1.5, simd_vec(x, impl = t)), ops[[op]](1.5, x))
      expect_identical(f(simd_vec(1:3, impl = t), 2L), get(op)(1:3, 2L))
    }
    r <- as.raw(c(1, 5, 9))
    expect_identical(simd_vec(r, impl = t) == as.raw(5), r == as.raw(5))
    expect_identical(simd_vec(r, impl = t) < as.raw(5), r < as.raw(5))
  }
})

test_that("logical operators use three-valued logic, and bitwise ops on raw", {
  a <- c(TRUE, FALSE, NA, TRUE)
  b <- c(NA, TRUE, FALSE, TRUE)
  for (t in tiers_to_test()) {
    expect_identical(simd_vec(a, impl = t) & b, a & b)
    expect_identical(simd_vec(a, impl = t) | b, a | b)
    expect_identical(!simd_vec(a, impl = t), !a)
    expect_identical(simd_vec(c(0, 2, NA)) & TRUE, c(0, 2, NA) & TRUE)
    r <- as.raw(c(0x0f, 0xf0))
    expect_sv(simd_vec(r, impl = t) & as.raw(0x3c), r & as.raw(0x3c), t, TRUE)
    expect_sv(simd_vec(r, impl = t) | as.raw(0x3c), r | as.raw(0x3c), t, TRUE)
    expect_sv(!simd_vec(r, impl = t), !r, t, TRUE)
  }
})

test_that("complex operators: kernels for + and -, base for the rest", {
  z <- c(1 + 2i, -3i, NA)
  w <- c(2 - 1i, 1i, 1)
  expect_sv(simd_vec(z) + w, simd_add(z, w))
  expect_sv(simd_vec(z) - w, simd_sub(z, w))
  expect_sv(-simd_vec(z), simd_neg(z))
  expect_sv(simd_vec(z, impl = "none") * w, z * w, "none")
  expect_sv(simd_vec(z) / w, z / w)
  expect_sv(simd_vec(z)^2, z^2)
  expect_identical(simd_vec(z) == w, z == w)
  expect_identical(simd_vec(z) != z, z != z)
  expect_error(simd_vec(z) < w, "invalid comparison with complex values")
  expect_error(simd_vec(z) * c(w, 1), "lengths of 'x' \\(3\\) and 'y' \\(4\\)")
})

test_that("an operand of another class is rejected on either side", {
  msg <- function(arg, cls) {
    paste0(
      "'", arg, "' must be an atomic vector \\(double, integer, logical, raw, complex or ",
      "integer64\\), not ", cls
    )
  }
  sv <- simd_vec(1)
  # The simd_vec method wins over the other class's (chooseOpsMethod), so
  # there is no "Incompatible methods" warning.
  expect_no_warning(expect_error(Sys.Date() + sv, msg("x", "Date")))
  expect_no_warning(expect_error(sv + Sys.Date(), msg("y", "Date")))
  expect_no_warning(expect_error(sv == Sys.time(), msg("y", "POSIXct")))
  expect_no_warning(expect_error(as.difftime(1, units = "hours") * sv, msg("x", "difftime")))
  d <- structure(1, class = "foo")
  `Ops.foo` <- function(e1, e2) "foo"
  registerS3method("Ops", "foo", `Ops.foo`)
  expect_error(d + sv, msg("x", "foo"))
  expect_error(sv + d, msg("y", "foo"))
  # Paths that do not reach a kernel reject it too.
  expect_error(simd_vec(1i) == d, msg("y", "foo"))
  if (has_bit64()) expect_error(simd_vec(simd_as_integer64(2))^d, msg("y", "foo"))
  # Other types too, on either side.
  expect_error(sv + "a", msg("y", "character"))
  expect_error(list(1) * sv, msg("x", "list"))
  expect_error(sv == NULL, msg("y", "NULL"))
  expect_error(sv & sum, msg("y", "builtin"))
  # integer64 and simd_vec operands are taken.
  expect_sv(sv + simd_vec(2), 3)
})

# ---- Math ------------------------------------------------------------------

test_that("Math functions with a kernel equal the simd_ functions", {
  batch_expectations({
    x <- c(0.25, -0.75, 3.5, NA, NaN, Inf)
    pos <- c(0.25, 1, 3.5, NA, 1e300)
    fns <- c(
      "abs", "sign", "floor", "ceiling", "trunc", "exp", "expm1", "cos", "sin", "tan",
      "cospi", "sinpi", "tanpi", "atan", "cosh", "sinh", "tanh", "asinh",
      "cumsum", "cumprod", "cummax", "cummin"
    )
    on_pos <- c("sqrt", "log2", "log10", "log1p", "acosh")
    on_unit <- c("acos", "asin", "atanh")
    ns <- asNamespace("rsimd")
    for (t in tiers_to_test()) {
      for (f in c(fns, on_pos, on_unit)) {
        arg <- if (f %in% on_pos) pos else if (f %in% on_unit) c(-0.5, 0, 0.75, NA) else x
        if (f == "acosh") arg <- arg + 1
        kernel <- get(paste0("simd_", f), envir = ns)
        # Some give NaN with a warning, as base R does.
        want <- suppressWarnings(simd_with_impl(t, kernel(arg)))
        expect_sv(suppressWarnings(get(f)(simd_vec(arg, impl = t))), want, t)
      }
      expect_sv(round(simd_vec(x, impl = t)), simd_with_impl(t, simd_round(x)), t)
      expect_sv(round(simd_vec(x, impl = t), 2), simd_round(x, 2), t)
      expect_sv(round(simd_vec(x, impl = t), digits = -1), simd_round(x, -1), t)
      expect_sv(log(simd_vec(pos, impl = t)), simd_with_impl(t, simd_log(pos)), t)
      expect_sv(log(simd_vec(pos, impl = t), 2), simd_with_impl(t, simd_log(pos, 2)), t)
      expect_sv(log(simd_vec(pos, impl = t), base = 10), simd_with_impl(t, simd_log(pos, 10)), t)
    }
    # Integer and logical input.
    expect_sv(exp(simd_vec(0:2)), simd_exp(0:2))
    expect_sv(abs(simd_vec(c(-2L, NA))), c(2L, NA))
    expect_sv(cumsum(simd_vec(c(TRUE, TRUE))), 1:2)
    expect_error(exp(simd_vec(as.raw(1))), "non-numeric argument to mathematical function")
  })
})

test_that("Math functions without a kernel fall back to base R", {
  x <- c(0.5, 2.5, NA)
  for (f in c("gamma", "lgamma", "digamma", "trigamma")) {
    expect_sv(get(f)(simd_vec(x, impl = "none", check_na = TRUE)), get(f)(x), "none")
  }
  expect_sv(signif(simd_vec(pi * 100)), signif(pi * 100))
  expect_sv(signif(simd_vec(pi * 100), 2), signif(pi * 100, 2))
  z <- c(3 + 4i, -1i)
  expect_sv(abs(simd_vec(z)), abs(z))
  # exp and sqrt of complex use the kernels (test-complex-math.R).
  expect_sv(exp(simd_vec(z, impl = "none")), exp(z), "none")
  expect_sv(sqrt(simd_vec(z, impl = "none")), sqrt(z), "none")
  expect_sv(cumsum(simd_vec(z)), cumsum(z))
  expect_sv(round(simd_vec(z), 1), round(z, 1))
})

# ---- Summary, mean, anyNA -------------------------------------------------

test_that("Summary functions equal the simd_ functions", {
  x <- c(1.5, -2, NA, 4)
  fns <- list(
    sum = simd_sum, prod = simd_prod, min = simd_min, max = simd_max,
    range = simd_range
  )
  for (t in tiers_to_test()) {
    for (f in names(fns)) {
      for (na.rm in c(FALSE, TRUE)) {
        want <- simd_with_impl(t, fns[[f]](x, na.rm = na.rm))
        expect_identical(get(f)(simd_vec(x, impl = t), na.rm = na.rm), want)
      }
    }
    l <- c(TRUE, NA, FALSE)
    expect_identical(any(simd_vec(l, impl = t)), any(l))
    expect_identical(all(simd_vec(l, impl = t)), all(l))
    expect_identical(any(simd_vec(l, impl = t), na.rm = TRUE), any(l, na.rm = TRUE))
    expect_identical(all(simd_vec(l[-3], impl = t), na.rm = TRUE), TRUE)
    expect_identical(sum(simd_vec(1:10, impl = t)), 55L)
  }
  expect_identical(sum(simd_vec(1:3), 4:6, 0.5), 21.5)
  expect_identical(max(simd_vec(1:3), 10L), 10L)
  expect_identical(range(simd_vec(c(3, 1)), c(-1, 9)), c(-1, 9))
  expect_warning(expect_identical(min(simd_vec(numeric())), Inf), "no non-missing arguments")
  z <- c(1 + 1i, 2i)
  expect_identical(sum(simd_vec(z)), sum(z))
  expect_identical(prod(simd_vec(z)), prod(z))
  expect_error(min(simd_vec(z)), "invalid 'type' \\(complex\\)")
})

test_that("mean() uses simd_mean, with base R for trim", {
  x <- c(1, 2, 4, NA, 100)
  for (t in tiers_to_test()) {
    expect_identical(mean(simd_vec(x, impl = t)), simd_with_impl(t, simd_mean(x)))
    expect_identical(
      mean(simd_vec(x, impl = t), na.rm = TRUE),
      simd_with_impl(t, simd_mean(x, na.rm = TRUE))
    )
  }
  expect_identical(mean(simd_vec(x), trim = 0.25, na.rm = TRUE), mean(x, trim = 0.25, na.rm = TRUE))
  expect_identical(mean(simd_vec(1:4)), 2.5)
  expect_identical(mean(simd_vec(c(1i, 3i))), 2i)
})

test_that("anyNA() answers from a known flag without reading the data", {
  calls <- 0L
  local_mocked_bindings(simd_any_na = function(x) {
    calls <<- calls + 1L
    anyNA(unclass(x))
  })
  expect_false(anyNA(simd_vec(c(1, 2), check_na = TRUE)))
  expect_true(anyNA(simd_vec(c(1, NA), check_na = TRUE)))
  expect_identical(calls, 2L) # the two checks made by simd_vec()
  expect_true(anyNA(simd_vec(c(1, NA))))
  expect_false(anyNA(simd_vec(c(1, 2))))
  expect_identical(calls, 4L)
})

# ---- The NA-free flag -----------------------------------------------------

test_that("the flag is kept by ops that cannot make a missing value", {
  x <- simd_vec(c(-1.5, 2.25, Inf), check_na = TRUE)
  for (f in list(
    abs, sign, floor, ceiling, trunc, round, function(v) -v,
    function(v) round(v, 1), cummax, cummin, simd_conj, simd_re, simd_im
  )) {
    expect_identical(simd_na_free(f(x)), TRUE)
  }
  i <- simd_vec(c(-3L, 5L), check_na = TRUE)
  for (f in list(
    abs, function(v) -v, simd_abs_wrap, simd_neg_wrap, simd_popcount,
    simd_lzcnt, simd_tzcnt, cummax
  )) {
    expect_identical(simd_na_free(f(i)), TRUE)
  }
  expect_identical(simd_na_free(simd_pmin(x, simd_vec(0, check_na = TRUE))), TRUE)
  expect_identical(simd_na_free(simd_pmax(x, 0)), TRUE)
  expect_identical(simd_na_free(simd_pmin(x, c(0, 1, 2))), NULL)
  expect_identical(simd_na_free(simd_clamp(x, -1, 1)), TRUE)
  # 1/0 is Inf and 1/Inf is 0; copysign only moves a sign bit.
  z <- simd_vec(c(0, -0, Inf, 2), check_na = TRUE)
  expect_identical(simd_na_free(simd_recip(z)), TRUE)
  expect_identical(simd_na_free(simd_recip(i)), TRUE)
  expect_identical(simd_na_free(simd_copysign(z, -1)), TRUE)
  expect_identical(simd_na_free(simd_copysign(z, x[c(1, 2, 3, 1)])), TRUE)
  expect_identical(simd_na_free(simd_copysign(z, NaN)), NULL)
  expect_identical(simd_na_free(simd_copysign(z, c(1, NA, 1, 1))), NULL)
  expect_identical(simd_na_free(simd_recip(simd_vec(c(1, NaN), check_na = TRUE))), NULL)
  expect_identical(simd_na_free(simd_recip(simd_vec(c(1, 2)))), NULL)
  # FALSE and unknown are not kept.
  expect_identical(simd_na_free(abs(simd_vec(c(1, NA), check_na = TRUE))), NULL)
  expect_identical(simd_na_free(abs(simd_vec(c(1, 2)))), NULL)
})

test_that("the flag is dropped by ops that can make a missing value", {
  x <- simd_vec(c(-1.5, 2.25, Inf), check_na = TRUE)
  for (f in list(
    function(v) v + v, function(v) v - 1, function(v) v * 2, function(v) v / v,
    function(v) v^2, function(v) v %% 2, function(v) v %/% 2, exp, sqrt, log, sin,
    cumsum, cumprod, gamma, simd_sigmoid, simd_softmax
  )) {
    expect_identical(simd_na_free(suppressWarnings(f(x))), NULL)
  }
  i <- simd_vec(c(-3L, 5L), check_na = TRUE)
  for (f in list(
    function(v) v + 1L, function(v) v * v, function(v) v %/% 2L,
    function(v) v %% 2L, function(v) simd_add_wrap(v, 1L), function(v) simd_mul_wrap(v, v),
    function(v) simd_shl(v, 1), simd_bit_not, cumsum
  )) {
    expect_identical(simd_na_free(f(i)), NULL)
  }
  # A wrapped result of -2^31 reads as NA.
  big <- simd_vec(-.Machine$integer.max, check_na = TRUE)
  expect_identical(bare(simd_sub_wrap(big, 1L)), NA_integer_)
})

test_that("round with a missing digits drops the flag", {
  # digits is not an operand, yet NA or NaN digits make every element
  # missing, so the result must not be stamped NA-free.
  d <- c(1, 2.5, 3, -7.25, 1e300)
  for (t in simd_available()) {
    x <- simd_vec(d, check_na = TRUE, impl = t)
    i <- simd_vec(c(1L, -2L, 3L), check_na = TRUE, impl = t)
    cases <- list(
      list(r = simd_round(x, NA), n = 5),
      list(r = simd_round(x, NaN), n = 5),
      list(r = round(x, NA_real_), n = 5),
      list(r = simd_round(i, NA_integer_), n = 3)
    )
    for (cs in cases) {
      r <- cs$r
      expect_identical(simd_na_free(r), NULL, info = t)
      expect_true(anyNA(r), info = t)
      expect_true(simd_any_na(r), info = t)
      expect_identical(simd_count_na(r), cs$n, info = t)
      expect_true(is.na(simd_min(r)), info = t)
      expect_true(is.na(simd_sum(r)), info = t)
    }
    expect_identical(bare(simd_round(x, NA)), round(d, NA))
    expect_identical(simd_na_free(simd_round(x, 2)), TRUE, info = t)
  }
})

test_that("conversions keep the flag only when they cannot overflow", {
  d <- simd_vec(c(1.5, -3e9), check_na = TRUE)
  i <- simd_vec(c(1L, -5L), check_na = TRUE)
  l <- simd_vec(c(TRUE, FALSE), check_na = TRUE)
  expect_identical(simd_na_free(simd_as_double(i)), TRUE)
  expect_identical(simd_na_free(simd_as_logical(d)), TRUE)
  expect_identical(simd_na_free(simd_as_integer(l)), TRUE)
  if (has_bit64()) expect_identical(simd_na_free(simd_as_integer64(i)), TRUE)
  expect_identical(simd_na_free(simd_as_integer(d, "saturating")), TRUE)
  expect_identical(simd_na_free(suppressWarnings(simd_as_raw(d))), TRUE)
  expect_identical(simd_na_free(suppressWarnings(simd_as_integer(d))), NULL)
  expect_identical(simd_na_free(simd_as_integer(d, "truncating")), NULL)
})

test_that("subsetting keeps TRUE unless a missing element is selected, then FALSE", {
  x <- simd_vec(c(1, 2, 3), impl = "none", check_na = TRUE)
  expect_sv(x[2:3], c(2, 3), "none", TRUE)
  expect_sv(x[-1], c(2, 3), "none", TRUE)
  expect_sv(x[c(TRUE, FALSE, TRUE)], c(1, 3), "none", TRUE)
  expect_sv(x[0], numeric(), "none", TRUE)
  expect_sv(x[c(1, NA)], c(1, NA), "none", FALSE)
  expect_sv(x[5], NA_real_, "none", FALSE)
  expect_sv(x[c(TRUE, NA, FALSE)], c(1, NA), "none", FALSE)
  expect_identical(x[], x)
  y <- simd_vec(c(1, NA), check_na = TRUE)
  expect_sv(y[1], 1)
  expect_sv(simd_vec(c(1, 2))[1], 1)
  expect_error(x[1, 2], "incorrect number of dimensions")
})

test_that("c() is TRUE only when every part is known NA-free", {
  a <- simd_vec(c(1, 2), check_na = TRUE)
  b <- simd_vec(3, check_na = TRUE)
  f <- simd_vec(c(NA, 1), check_na = TRUE)
  expect_sv(c(a, b), c(1, 2, 3), na_free = TRUE)
  expect_sv(c(a, 4, 5L), c(1, 2, 4, 5), na_free = TRUE)
  expect_sv(c(a, c(4, NA)), c(1, 2, 4, NA))
  expect_sv(c(a, f), c(1, 2, NA, 1), na_free = FALSE)
  expect_sv(c(simd_vec(1), f), c(1, NA, 1), na_free = FALSE)
  expect_sv(c(a, simd_vec(3)), c(1, 2, 3))
  expect_sv(c(simd_vec(1:2), TRUE), c(1L, 2L, 1L))
  expect_sv(c(simd_vec(1:2), 0.5), c(1, 2, 0.5))
  expect_sv(c(simd_vec(1), c(a = 2)), c(1, 2))
  expect_sv(c(simd_vec(1), 1i), c(1 + 0i, 1i))
  # Long parts: a simd_vec of a shared vector (which R wraps) and computed
  # ones, the parts keeping their flags.
  d <- as.double(seq_len(1000))
  w <- simd_vec(d, check_na = TRUE)
  expect_sv(c(w, simd_vec(-d, check_na = TRUE), 1), c(d, -d, 1), na_free = TRUE)
  expect_sv(c(w, w + 1, simd_vec(1:2)), c(d, d + 1, 1, 2))
  expect_sv(c(simd_vec(1:3), d), c(1, 2, 3, d))
  expect_true(simd_na_free(w))
  expect_error(c(simd_vec(1), "a"), "'x' must be an atomic vector")
  expect_error(c(simd_vec(1), list(2)), "'x' must be an atomic vector")
})

test_that("assignment and length<- keep the class and pin", {
  x <- simd_vec(c(1, 2, 3), impl = "none", check_na = TRUE)
  x[2] <- 10
  expect_sv(x, c(1, 10, 3), "none")
  x[] <- 0
  expect_sv(x, c(0, 0, 0), "none")
  y <- simd_vec(1:3)
  y[2] <- 2.5
  expect_sv(y, c(1, 2.5, 3))
  y[4] <- simd_vec(9)
  expect_sv(y, c(1, 2.5, 3, 9))

  z <- simd_vec(c(1, 2, 3), impl = "none", check_na = TRUE)
  length(z) <- 2
  expect_sv(z, c(1, 2), "none", TRUE)
  length(z) <- 4
  expect_sv(z, c(1, 2, NA, NA), "none")
})

test_that("[ with a logical or non-negative index keeps base R's values and the flag", {
  d <- c(1.5, 2.5, 3.5, 4.5)
  x <- simd_vec(d, check_na = TRUE)
  for (i in list(
    c(TRUE, FALSE), c(TRUE, FALSE, TRUE, TRUE), logical(), c(3, 1), c(0, 2), 2.9,
    integer(), c(4L, 4L), matrix(c(3, 1)), c(1.999, 0.5)
  )) {
    out <- x[i]
    expect_identical(simd_unwrap(out), d[i], info = deparse(i))
    expect_true(simd_na_free(out), info = deparse(i))
  }
  # Past the end, NA or a name in the index: NA elements, found by the scan.
  for (i in list(c(TRUE, FALSE, TRUE, TRUE, TRUE), c(TRUE, NA), 5, c(1, NA), 2^40, Inf, "a")) {
    out <- x[i]
    expect_identical(simd_unwrap(out), d[i], info = deparse(i))
    expect_false(simd_na_free(out), info = deparse(i))
  }
  for (i in list(-1, c(-4, 0, -2), -5, simd_vec(c(TRUE, FALSE)), factor(c("b", "a")), x > 2)) {
    out <- x[i]
    expect_identical(simd_unwrap(out), d[if (is_simd_vec(i)) simd_unwrap(i) else i], info = deparse(i))
    expect_true(simd_na_free(out), info = deparse(i))
  }
  # An x of unknown flag is not scanned.
  expect_null(simd_na_free(simd_vec(d)[1:2]))
  y <- simd_vec(d, "none")
  expect_identical(simd_impl(y[c(TRUE, FALSE)]), "none")
})

test_that("[[ returns a plain element", {
  x <- simd_vec(c(1, 2), impl = "none")
  expect_identical(x[[2]], 2)
  expect_identical(simd_vec(c(TRUE, FALSE))[[1]], TRUE)
})

test_that("the flag lets the functions skip their missing-value checks", {
  x <- simd_vec(c(1, 2, 3), check_na = TRUE)
  expect_false(.debug_opts(x)$na_check)
  expect_true(.debug_opts(simd_vec(c(1, 2, 3)))$na_check)
  expect_true(.debug_opts(simd_vec(c(1, NA), check_na = TRUE))$na_check)
  expect_identical(sum(x), 6)
  expect_identical(bare(x + 1), c(2, 3, 4))
})

# ---- A stale flag can never be read ---------------------------------------

test_that("the flag belongs to the object it was set on", {
  x <- simd_vec(c(1, 2, 3), check_na = TRUE)
  expect_identical(simd_na_free(x), TRUE)
  # A base function that copies the attributes onto new data shares the
  # token: the copy and the original both read as unknown.
  y <- pmin(x, NA)
  expect_null(simd_na_free(y))
  expect_identical(max(y), NA_real_)
  expect_null(simd_na_free(x))
  expect_identical(sum(x), 6)

  # A forged flag (no token) or a token of another object is unknown.
  f <- structure(c(1, NA), class = "simd_vec", rsimd_na_free = TRUE)
  expect_null(simd_na_free(f))
  expect_identical(sum(f), NA_real_)
  x <- simd_vec(c(1, 2, 3), check_na = TRUE)
  g <- c(NA, 2, 3)
  attributes(g) <- attributes(x)
  expect_null(simd_na_free(g))
  expect_identical(sum(g), NA_real_)
  # Setting the old flag attribute on an object with a valid token of its
  # own changes nothing.
  f <- simd_vec(c(1L, NA, 3L), check_na = TRUE)
  attr(f, "rsimd_na_free") <- TRUE
  expect_false(simd_na_free(f))
  expect_identical(sum(f), NA_integer_)
  expect_true(anyNA(f))

  # Any change of class voids the flag: without one, [<- can change the
  # data in place.
  g <- simd_vec(1:3, check_na = TRUE)
  class(g) <- NULL
  g[2] <- NA
  class(g) <- "simd_vec"
  expect_null(simd_na_free(g))
  expect_true(anyNA(g))
  expect_identical(sum(g), NA_integer_)
  expect_identical(min(g), NA_integer_)
  expect_identical(mean(g), NA_real_)
  expect_identical(simd_count_na(g), 1)
  g <- simd_vec(1:3, check_na = TRUE)
  oldClass(g) <- c("foo", oldClass(g))
  expect_null(simd_na_free(g))
  # Nor does putting back the class vector of another simd_vec, made in R
  # or in C, flagged or not: a flagged object has a class vector of its own.
  k <- simd_vec(c(1, 2), check_na = TRUE)
  for (other in list(simd_vec(1), simd_vec(1) + 1, k, k + k, simd_vec(1L) + 1L)) {
    for (g in list(simd_vec(1:3, check_na = TRUE), simd_add(k, k))) {
      class(g) <- NULL
      g[2] <- NA
      class(g) <- class(other)
      expect_null(simd_na_free(g))
    }
  }
  expect_true(simd_na_free(k))

  # A reloaded flag is unknown (the token's pointer is not saved).
  x <- simd_vec(c(1, 2, 3), check_na = TRUE)
  expect_null(simd_na_free(unserialize(serialize(x, NULL))))
  path <- tempfile(fileext = ".rds")
  saveRDS(simd_vec(c(1, 2), check_na = TRUE), path)
  expect_null(simd_na_free(readRDS(path)))
  unlink(path)

  # Raw is NA-free whatever the attributes say.
  r <- simd_vec(as.raw(1:3))
  r2 <- pmin(r, as.raw(2))
  expect_identical(simd_na_free(r2), TRUE)
})

test_that("ordinary use keeps the flag", {
  x <- simd_vec(c(3, 1, 2), check_na = TRUE)
  f <- function(v) sum(v)
  f(x)
  l <- list(x)
  y <- x[1:2]
  z <- -x
  s <- sort(x)
  print_out <- capture.output(print(x), str(x))
  expect_identical(simd_na_free(x), TRUE)
  expect_identical(simd_na_free(l[[1L]]), TRUE)
  expect_identical(simd_na_free(y), TRUE)
  expect_identical(simd_na_free(z), TRUE)
  expect_identical(simd_na_free(s), TRUE)
  # Replacement on a copy: R hands the method a copy that shares the
  # token, which the method releases, so the original keeps its flag.
  w <- x
  w[1] <- NA
  expect_identical(simd_na_free(x), TRUE)
  w <- x
  w[[1]] <- NA
  expect_identical(simd_na_free(x), TRUE)
  w <- x
  length(w) <- 2
  expect_identical(simd_na_free(x), TRUE)
  w <- x
  simd_impl(w) <- "none"
  expect_identical(simd_na_free(x), TRUE)
  h <- function(v) {
    v[1] <- NA
    v
  }
  expect_identical(sum(h(x)), NA_real_)
  expect_identical(simd_na_free(x), TRUE)
  # A function that returns its input copies it, token and all.
  i <- simd_vec(1:3, check_na = TRUE)
  expect_sv(simd_as_integer(i), 1:3, NULL, TRUE)
  expect_identical(simd_na_free(i), TRUE)
  # simd_impl<- and simd_vec(impl =) keep the flag of an unshared object.
  u <- simd_vec(c(1, 2), check_na = TRUE)
  simd_impl(u) <- "none"
  expect_sv(u, c(1, 2), "none", TRUE)
  expect_sv(simd_vec(simd_vec(c(1, 2), check_na = TRUE), impl = "none"), c(1, 2), "none", TRUE)
})

test_that("[[<- assigns one element and drops the flag", {
  w <- simd_vec(c(1L, 2L, 3L), check_na = TRUE)
  w[[2]] <- NA
  expect_sv(w, c(1L, NA, 3L))
  expect_identical(sum(w), NA_integer_)
  expect_true(anyNA(w))
  expect_identical(simd_which_na(w), 2L)
  w2 <- simd_vec(c(1, 2, 3), impl = "none", check_na = TRUE)
  w2[[2]] <- NA
  expect_sv(w2, c(1, NA, 3), "none")
  expect_identical(max(w2), NA_real_)
  w2[[5]] <- 9
  expect_sv(w2, c(1, NA, 3, NA, 9), "none")
  expect_error(w2[[1:2]] <- 1, "more than one element")
  expect_error(w2[[integer()]] <- 1, "less than one element")
  expect_error(w2[[1]] <- 1:2, "more elements supplied")
  expect_error(w2[[1]] <- numeric(), "replacement has length zero")

  skip_if_not(has_bit64(), "bit64 is not installed")
  v <- simd_vec(simd_as_integer64(c(1, 2, 3)), check_na = TRUE)
  v[[2]] <- 5L
  expect_identical(bare(v), simd_as_integer64(c(1, 5, 3)))
  v[[5]] <- 7L
  expect_identical(bare(v), simd_as_integer64(c(1, 5, 3, NA, 7)))
  v[7] <- 1L
  expect_identical(bare(v), simd_as_integer64(c(1, 5, 3, NA, 7, NA, 1)))
  # A double value is truncated toward zero without a warning, as bit64's
  # [[<- does.
  expect_silent(v[[1]] <- 2.5)
  expect_silent(v[2] <- -2.5)
  expect_identical(bare(v), simd_as_integer64(c(2, -2, 3, NA, 7, NA, 1)))
  skip_if_not_installed("bit64")
  w <- bit64::as.integer64(c(1, 5, 3))
  w[[1]] <- 2.5
  w[2] <- -2.5
  expect_identical(bare(v)[1:3], w)
})

test_that("[<- and [[<- reject a value of a type rsimd does not take, naming 'value'", {
  data <- list(1:3, c(1, 2, 3))
  if (has_bit64()) data <- c(data, list(simd_as_integer64(1:3)))
  values <- list(character = "a", list = list(1), factor = factor("z"), Date = Sys.Date())
  for (d in data) {
    for (what in names(values)) {
      msg <- paste0("^'value' must be an atomic vector .*, not ", what, "$")
      x <- simd_vec(d, check_na = TRUE)
      expect_error(x[2] <- values[[what]], msg)
      expect_error(x[] <- values[[what]], msg)
      expect_error(x[[2]] <- values[[what]], msg)
      expect_sv(x, d)
    }
    x <- simd_vec(d)
    expect_error(x[2] <- NULL, "replacement has length zero")
    expect_error(x[[2]] <- NULL, "replacement has length zero")
  }
})

test_that("base functions on a flagged simd_vec agree with plain data", {
  batch_expectations({
    data <- list(
      double = c(1.5, -2, 3e10),
      integer = c(1L, -2L, 3L),
      logical = c(TRUE, FALSE, TRUE),
      complex = c(1 + 1i, -2i, 3),
      integer64 = if (has_bit64()) simd_as_integer64(c(1, -2, 3))
    )
    data <- Filter(Negate(is.null), data)
    fns <- list(
      pmin = function(x) pmin(x, NA),
      pmax = function(x) pmax(x, NA),
      ifelse = function(x) ifelse(c(TRUE, FALSE, TRUE), x, NA),
      replace = function(x) replace(x, 2, NA),
      storage_integer = function(x) {
        storage.mode(x) <- "integer"
        x
      },
      storage_double = function(x) {
        storage.mode(x) <- "double"
        x
      },
      mode_integer = function(x) {
        mode(x) <- "integer"
        x
      },
      Re = Re, Im = Im, Mod = Mod, Arg = Arg, Conj = Conj,
      rev = rev, sort = sort, unique = unique, head = function(x) head(x, 2),
      tail = function(x) tail(x, 2), round = function(x) round(x, 1), signif = signif,
      sqrt = sqrt, log = log, cumsum = cumsum, diff = diff,
      is_na_assign = function(x) {
        is.na(x) <- 2
        x
      },
      sub2_assign = function(x) {
        x[[2]] <- NA
        x
      },
      sub_assign = function(x) {
        x[2] <- NA
        x
      },
      append = function(x) append(x, NA),
      rep = function(x) rep(x, 2),
      rep_len = function(x) rep_len(x, 5),
      length_assign = function(x) `length<-`(x, 5),
      na_index = function(x) x[c(1, NA)],
      c = function(x) c(x, NA),
      unclass_reclass = function(x) {
        y <- unclass(x)
        y[2] <- NA
        class(y) <- class(x)
        y
      },
      attributes_copy = function(x) {
        y <- c(NA, 1, 2)
        attributes(y) <- attributes(x)
        y
      },
      attr_copy = function(x) {
        y <- x
        attr(y, "extra") <- 1
        y[2] <- NA
        y
      },
      structure = function(x) structure(rep(NA, length(x)), class = class(x)),
      times_na = function(x) x * NA,
      mod_zero = function(x) x %% 0,
      ifelse_cond = function(x) ifelse(x == 1, x, NA),
      split_unsplit = function(x) {
        f <- c(1, 2, 1)
        unsplit(split(x, f), f)
      },
      serialize = function(x) unserialize(serialize(x, NULL)),
      names_dim = function(x) {
        names(x) <- c("a", "b", "c")
        dim(x) <- c(3, 1)
        x[2] <- NA
        x
      }
    )
    # Data without a class, which rsimd takes (structure() above gives
    # plain data the class "numeric").
    plain <- function(r) if (is.object(r)) bare(r) else r
    plain_any_na <- function(r) simd_any_na(plain(r))
    for (type in names(data)) {
      for (f in names(fns)) {
        label <- paste(f, "on", type)
        want <- tryCatch(suppressWarnings(fns[[f]](data[[type]])), error = function(e) NULL)
        if (is.null(want) || !is.atomic(unclass(want))) next
        x <- simd_vec(data[[type]], check_na = TRUE)
        got <- tryCatch(suppressWarnings(fns[[f]](x)), error = function(e) NULL)
        if (is.null(got)) next
        # The core property: a TRUE flag is never on data with a missing value.
        if (is_simd_vec(got) && isTRUE(simd_na_free(got))) {
          check_identical(plain_any_na(got), FALSE, info = label)
        }
        # Base R on a plain integer64 works on its double bits without bit64;
        # bit64 has no length<- method (base pads with double NA), and an
        # integer64 with a double is double in rsimd (integer64 in bit64).
        i64_differs <- !isNamespaceLoaded("bit64") || f %in% c("length_assign", "mod_zero")
        if (type == "integer64" && i64_differs) next
        check_identical(anyNA(got), plain_any_na(want), info = label)
        if (is.numeric(got) || is.complex(got) || is.logical(got)) {
          check_identical(simd_which_na(got), simd_which_na(plain(want)), info = label)
          if (!inherits(got, "integer64") && !is.complex(got)) {
            check_identical(
              suppressWarnings(max(got)), suppressWarnings(max(.subset(want, seq_along(want)))),
              info = label
            )
          }
        }
      }
    }
  })
})

# ---- Other base generics ---------------------------------------------------

test_that("range(finite = TRUE) follows base R", {
  expect_identical(range(simd_vec(c(5, Inf)), finite = TRUE), c(5, 5))
  expect_identical(range(simd_vec(c(1, Inf, NA)), finite = TRUE), c(1, 1))
  expect_identical(range(simd_vec(c(1, NA)), finite = FALSE), c(NA_real_, NA_real_))
  expect_identical(range(simd_vec(c(2, -Inf)), 7, finite = TRUE), c(2, 7))
})

test_that("c(), [ and as_simd_vec() take base R's other arguments", {
  expect_sv(c(simd_vec(1L), recursive = TRUE), 1L)
  expect_sv(c(simd_vec(1L), 2L, use.names = FALSE), c(1L, 2L))
  w <- simd_vec(c(4, 5, 6), impl = "none")
  expect_sv(w[1, drop = FALSE], 4, "none")
  expect_sv(w[2:3, drop = TRUE], c(5, 6), "none")
  x <- simd_vec(c(1, NA))
  expect_sv(as_simd_vec(x, impl = "none", check_na = TRUE), c(1, NA), "none", FALSE)
  expect_identical(as_simd_vec(x), x)
})

test_that("setting names, dimensions or dimnames unwraps", {
  z <- simd_vec(c(1, 2, 3), impl = "none", check_na = TRUE)
  y <- z
  names(y) <- c("a", "b", "c")
  expect_identical(y, c(a = 1, b = 2, c = 3))
  y <- z
  dim(y) <- c(3L, 1L)
  expect_identical(y, matrix(c(1, 2, 3), 3))
  dimnames(y) <- list(NULL, "v")
  expect_identical(y, matrix(c(1, 2, 3), 3, dimnames = list(NULL, "v")))
  expect_error(dimnames(z) <- list("a"), "'dimnames' applied to non-array")
  # The original keeps its flag.
  expect_sv(z, c(1, 2, 3), "none", TRUE)
  # NULL leaves the simd_vec as it is.
  y <- simd_vec(c(1, 2, 3), impl = "none", check_na = TRUE)
  names(y) <- NULL
  dim(y) <- NULL
  dimnames(y) <- NULL
  expect_sv(y, c(1, 2, 3), "none", TRUE)
  expect_sv(unname(simd_vec(c(1, 2, 3), impl = "none")), c(1, 2, 3), "none")
  # integer64 data keeps its class.
  if (has_bit64()) {
    w <- simd_vec(simd_as_integer64(c(1, 2)))
    names(w) <- c("a", "b")
    expect_identical(class(w), "integer64")
    expect_identical(names(w), c("a", "b"))
  }
  # Base functions that name or shape their result give base's answer.
  v <- c(3, 1, 2, 5)
  expect_identical(quantile(simd_vec(v)), quantile(v))
  expect_identical(summary(simd_vec(v)), summary(v))
  expect_identical(stats::setNames(simd_vec(v), letters[1:4]), stats::setNames(v, letters[1:4]))
  expect_identical(table(simd_vec(v), dnn = NULL), table(v, dnn = NULL))
  expect_identical(matrix(simd_vec(v), 2), matrix(v, 2))
  expect_identical(outer(simd_vec(v), simd_vec(v)), outer(v, v))
  expect_identical(
    unname(coef(stats::lm(simd_vec(v) ~ seq_along(v)))),
    unname(coef(stats::lm(v ~ seq_along(v))))
  )
})

test_that("is.na() returns a plain logical vector for every type", {
  expect_identical(is.na(simd_vec(c(1, NA, NaN))), c(FALSE, TRUE, TRUE))
  expect_identical(is.na(simd_vec(c(1L, NA))), c(FALSE, TRUE))
  expect_identical(is.na(simd_vec(c(1i, NA))), c(FALSE, TRUE))
  expect_identical(is.na(simd_vec(as.raw(1))), FALSE)
  if (has_bit64()) expect_identical(is.na(simd_vec(simd_as_integer64(c(1, NA)))), c(FALSE, TRUE))
})

test_that("complex simd_vecs use the complex kernels", {
  z <- c(3 + 4i, -1i, NA, 2)
  for (t in tiers_to_test()) {
    x <- simd_vec(z, impl = t)
    expect_sv(Re(x), simd_with_impl(t, simd_re(z)), t)
    expect_sv(Im(x), simd_with_impl(t, simd_im(z)), t)
    expect_sv(Mod(x), simd_with_impl(t, simd_abs(z)), t)
    expect_sv(abs(x), simd_with_impl(t, simd_abs(z)), t)
    expect_sv(Arg(x), simd_with_impl(t, simd_arg(z)), t)
    expect_sv(Conj(x), simd_with_impl(t, simd_conj(z)), t)
    expect_sv(cumsum(x), simd_with_impl(t, simd_cumsum(z)), t)
    expect_sv(cumprod(x), simd_with_impl(t, simd_cumprod(z)), t)
    expect_identical(prod(x, na.rm = TRUE), simd_with_impl(t, simd_prod(z, na.rm = TRUE)))
    expect_identical(mean(x), simd_with_impl(t, simd_mean(z)))
  }
  # Numeric data: base R's result types.
  expect_sv(Re(simd_vec(1:2)), c(1, 2))
  expect_sv(Conj(simd_vec(1:2)), c(1, 2))
  expect_sv(Mod(simd_vec(-1L)), 1)
  expect_sv(Arg(simd_vec(-1)), pi)
  expect_identical(Re(simd_vec(complex(real = 1, imaginary = NA), check_na = TRUE))[[1]], 1)
  expect_false(anyNA(Re(simd_vec(complex(real = 1, imaginary = NA), check_na = TRUE))))
})

test_that("all.equal() compares data, pin and flag, not the token", {
  a <- simd_vec(c(1, 2), impl = "none", check_na = TRUE)
  b <- simd_vec(c(1, 2), impl = "none", check_na = TRUE)
  expect_false(identical(a, b))
  expect_true(all.equal(a, b))
  expect_true(all.equal(a, simd_vec(c(1, 2 + 1e-12), impl = "none", check_na = TRUE)))
  expect_match(all.equal(a, simd_vec(c(1, 3), impl = "none", check_na = TRUE)), "Mean relative")
  expect_match(all.equal(a, simd_vec(c(1, 2), check_na = TRUE)), "pins differ: none vs auto")
  expect_match(all.equal(a, simd_vec(c(1, 2), impl = "none")), "flags differ: TRUE vs unknown")
  expect_match(all.equal(a, c(1, 2)), "target is simd_vec, current is numeric")
})

# ---- Base functions, unwrapping, printing ---------------------------------

test_that("rev, sort, head and tail keep the class; rep and unique do not", {
  x <- simd_vec(c(3, 1, 2), impl = "none", check_na = TRUE)
  expect_sv(rev(x), c(2, 1, 3), "none", TRUE)
  expect_sv(sort(x), c(1, 2, 3), "none", TRUE)
  expect_sv(head(x, 2), c(3, 1), "none", TRUE)
  expect_sv(tail(x, 1), 2, "none", TRUE)
  expect_false(is_simd_vec(rep(x, 2)))
  expect_false(is_simd_vec(unique(x)))
})

test_that("match and %in% compare the data, not its 15-digit strings", {
  d <- c(0.1 + 0.2, 1/3, 1e15 + 1)
  tab <- c(0.3, 0.333333333333333, 1e15 + 2)
  x <- simd_vec(d, impl = "none", check_na = TRUE)
  expect_identical(match(x, tab), match(d, tab))
  expect_identical(match(tab, x), match(tab, d))
  expect_identical(match(x, x), 1:3)
  expect_identical(x %in% 0.3, d %in% 0.3)
  expect_identical(is.element(x, d), rep(TRUE, 3))
  expect_identical(setdiff(x, tab), setdiff(d, tab))
  expect_identical(match(simd_vec(c(2L, NA)), c(NA, 2L)), c(2L, 1L))
  expect_identical(mtfrm(x), d)
  skip_if_not_installed("bit64")
  b <- simd_vec(bit64::as.integer64(c("9007199254740993", NA)))
  expect_identical(match(b, bit64::as.integer64(c(NA, "9007199254740993"))), c(2L, 1L))
  expect_identical(match(b, bit64::as.integer64("9007199254740992")), c(NA_integer_, NA_integer_))
})

test_that("unwrapping returns plain vectors", {
  x <- simd_vec(c(1.5, -2), impl = "none")
  expect_identical(as.vector(x), c(1.5, -2))
  expect_identical(as.vector(x, "integer"), c(1L, -2L))
  expect_identical(as.double(x), c(1.5, -2))
  expect_identical(as.numeric(x), c(1.5, -2))
  expect_identical(as.integer(x), c(1L, -2L))
  expect_identical(as.logical(x), c(TRUE, TRUE))
  expect_identical(as.double(simd_vec(1:2)), c(1, 2))
  expect_identical(data.frame(v = x), data.frame(v = c(1.5, -2)))
  expect_identical(as.data.frame(x)[[1L]], c(1.5, -2))
  expect_identical(names(as.data.frame(x)), "x")
  expect_identical(attr(unclass(x), "rsimd_impl"), "none")
  expect_identical(var(x), var(c(1.5, -2)))
})

test_that("print, format and str", {
  x <- simd_vec(c(1, 2.5, NA), impl = "none")
  expect_output(
    print(x), "<simd_vec[3] double impl=none na_free=unknown>\n[1] 1.0 2.5  NA",
    fixed = TRUE
  )
  capture.output(res <- withVisible(print(x)))
  expect_false(res$visible)
  expect_identical(res$value, x)
  expect_output(
    print(simd_vec(1:2, check_na = TRUE)), "<simd_vec[2] integer impl=auto na_free=TRUE>\n[1] 1 2",
    fixed = TRUE
  )
  expect_output(
    print(simd_vec(c(TRUE, NA), check_na = TRUE)),
    "<simd_vec[2] logical impl=auto na_free=FALSE>",
    fixed = TRUE
  )
  expect_identical(
    capture.output(print(simd_vec(numeric()))), "<simd_vec[0] double impl=auto na_free=unknown>"
  )
  expect_identical(format(x), format(c(1, 2.5, NA)))
  expect_output(str(x), "simd_vec [1:3] double impl=none na_free=unknown 1.0 2.5  NA", fixed = TRUE)
  expect_output(
    str(simd_vec(1:20)), "simd_vec [1:20] integer impl=auto na_free=unknown  1  2",
    fixed = TRUE
  )
  expect_output(str(simd_vec(1:20)), "10 ...", fixed = TRUE)
  expect_output(str(list(a = x)), "simd_vec [1:3] double", fixed = TRUE)
  expect_identical(
    capture.output(str(simd_vec(numeric()))), "simd_vec [0] double impl=auto na_free=unknown"
  )
})

test_that("print and str give the full length of a long simd_vec", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("RSIMD_EXTENDED_TESTS"), "true"),
    "RSIMD_EXTENDED_TESTS is not true"
  )
  skip_if(.Machine$sizeof.pointer < 8, "no long vectors on 32-bit platforms")
  x <- simd_vec(raw(2^31 + 1), check_na = TRUE)
  out <- capture.output(print(x, max = 3))
  expect_identical(out[1L], "<simd_vec[2147483649] raw impl=auto na_free=TRUE>")
  expect_output(str(x), "simd_vec [1:2147483649] raw impl=auto na_free=TRUE 00 00", fixed = TRUE)
})

# ---- integer64 -------------------------------------------------------------

test_that("integer64 simd_vecs run the integer64 kernels", {
  skip_if_not(has_bit64(), "bit64 is not installed")
  v <- simd_as_integer64(c(1, -2, NA, 4))
  x <- simd_vec(v, impl = "none")
  expect_identical(class(x + x), c("simd_vec", "integer64"))
  expect_identical(bare(x + x), simd_add(v, v))
  expect_identical(bare(x - 1L), simd_sub(v, 1L))
  expect_identical(bare(x * x), simd_mul(v, v))
  expect_identical(bare(-x), simd_neg(v))
  expect_identical(bare(abs(x)), simd_abs(v))
  expect_identical(bare(cumsum(x)), simd_cumsum(v))
  expect_identical(simd_impl(x * x), "none")
  expect_identical(sum(x, na.rm = TRUE), simd_sum(v, na.rm = TRUE))
  expect_identical(min(x, na.rm = TRUE), simd_min(v, na.rm = TRUE))
  expect_identical(max(x, na.rm = TRUE), simd_max(v, na.rm = TRUE))
  expect_identical(x == 1L, simd_eq(v, 1L))
  expect_warning(
    expect_identical(bare(x + 0.5), simd_add(simd_as_double(v), 0.5)),
    "integer64 coerced"
  )

  expect_identical(bare(x[c(2, NA)]), simd_as_integer64(c(-2, NA)))
  expect_identical(x[[1]], simd_as_integer64(1))
  expect_identical(bare(c(x, 5L)), simd_as_integer64(c(1, -2, NA, 4, 5)))
  y <- x
  y[1] <- 10L
  expect_identical(bare(y), simd_as_integer64(c(10, -2, NA, 4)))
  length(y) <- 5
  expect_identical(bare(y), simd_as_integer64(c(10, -2, NA, 4, NA)))
  expect_identical(as.double(x), c(1, -2, NA, 4))
  expect_identical(as.integer(x), c(1L, -2L, NA, 4L))
  expect_output(print(x), "<simd_vec[4] integer64 impl=none na_free=unknown>", fixed = TRUE)
})

test_that("integer64 simd_vecs work with bit64 loaded", {
  skip_if_not_installed("bit64")
  v <- bit64::as.integer64(c(1, -2, NA, 4))
  x <- simd_vec(v, impl = "none")
  expect_identical(bare(x + x), v + v)
  expect_identical(bare(x - v), v - v)
  expect_identical(bare(x * 3L), v * 3L)
  expect_identical(sum(x, na.rm = TRUE), sum(v, na.rm = TRUE))
  expect_identical(min(x, na.rm = TRUE), min(v, na.rm = TRUE))
  expect_identical(max(x, na.rm = TRUE), max(v, na.rm = TRUE))
  # Fallbacks reach bit64's methods.
  expect_identical(bare(x^2L), bare(.sv_strip(v^2L)))
  expect_identical(prod(x, na.rm = TRUE), prod(v, na.rm = TRUE))
  expect_identical(simd_impl(x^2L), "none")
  # bit64 claims operations with a plain integer64 on the left.
  expect_identical(class(v + x), "integer64")
  expect_output(print(x), "<NA>", fixed = TRUE)
})
