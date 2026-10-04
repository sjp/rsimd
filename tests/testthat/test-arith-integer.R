# Elementwise arithmetic on integer and logical vectors: bit-identical to
# base R on every tier, NA positions and the single overflow warning
# included.

imax <- .Machine$integer.max
ovf_msg <- "NAs produced by integer overflow"

# Integers around the overflow boundaries and the products near +-2^31.
boundary_ints <- function() {
  c(
    0L, 1L, -1L, 2L, -2L, 3L, 46340L, -46340L, 46341L, -46341L, 65536L, -65536L,
    imax, -imax, imax - 1L, -imax + 1L, 1073741824L, -1073741824L, NA
  )
}

test_that("checked add, sub and mul match base R with one warning", {
  p <- all_pairs(boundary_ints())
  ops <- list(simd_add = `+`, simd_sub = `-`, simd_mul = `*`)
  for (name in names(ops)) {
    f <- get(name)
    expected <- suppressWarnings(ops[[name]](p$x, p$y))
    got <- expect_tier_warnings(ovf_msg, f, p$x, p$y)
    expect_identical(got, expected, info = name)
    expect_tiers_give(expected, function() suppressWarnings(f(p$x, p$y)))
  }
})

test_that("no overflow warning when only NA operands are involved", {
  expect_tier_warnings(character(0), simd_add, c(1L, NA, 3L), c(NA, 2L, 4L))
  expect_identical(simd_add(c(1L, NA, 3L), c(NA, 2L, 4L)), c(NA, NA, 7L))
  # -imax - 1 is in two's complement range but is NA in R.
  expect_tier_warnings(ovf_msg, simd_sub, -imax, 1L)
  expect_tier_warnings(ovf_msg, simd_add, imax, 1L)
  expect_tier_warnings(ovf_msg, simd_mul, 65536L, 32768L)
  expect_tier_warnings(character(0), simd_mul, 46340L, 46340L)
})

test_that("checked ops match base R on random vectors of edge lengths", {
  for (n in c(0, 1, 3, 4, 5, 7, 8, 9, 15, 16, 17, 33, 67, 4095, 4097)) {
    x <- rand_vec("integer", n, seed = n + 1L) * 1000L
    y <- rand_vec("integer", n, seed = n + 2L)
    expect_tiers_give(suppressWarnings(x * y), function() suppressWarnings(simd_mul(x, y)))
    expect_tiers_give(suppressWarnings(x + y), function() suppressWarnings(simd_add(x, y)))
    expect_tiers_give(suppressWarnings(x - y), function() suppressWarnings(simd_sub(x, y)))
  }
})

test_that("wrapping ops wrap without a warning and keep NA", {
  p <- all_pairs(boundary_ints())
  xd <- as.double(p$x)
  yd <- as.double(p$y)
  # Exact 64-bit products, wrapped by hand in two halves of 16 bits.
  mul_wrap_ref <- function(a, b) {
    lo <- (a %% 65536) * b
    hi <- ((a %/% 65536) %% 65536) * b
    wrap_i32((hi %% 65536) * 65536 + lo %% 2^32)
  }
  expect_tier_warnings(character(0), simd_add_wrap, p$x, p$y)
  expect_tiers_give(wrap_i32(xd + yd), simd_add_wrap, p$x, p$y)
  expect_tiers_give(wrap_i32(xd - yd), simd_sub_wrap, p$x, p$y)
  expect_tiers_give(mul_wrap_ref(xd, yd), simd_mul_wrap, p$x, p$y)
  # The wrapped result -2^31 is the bit pattern of NA_integer_.
  expect_identical(simd_add_wrap(imax, 1L), NA_integer_)
  expect_identical(simd_add_wrap(imax, 2L), -imax)
  expect_identical(simd_mul_wrap(65536L, 65537L), 65536L)
  expect_identical(simd_sub_wrap(NA, 1L), NA_integer_)
})

test_that("neg and abs keep integers integer and NA as NA", {
  x <- boundary_ints()
  expect_tiers_give(-x, simd_neg, x)
  expect_tiers_give(abs(x), simd_abs, x)
  expect_tiers_give(-x, simd_neg_wrap, x)
  expect_tiers_give(abs(x), simd_abs_wrap, x)
  expect_tiers_give(-c(TRUE, FALSE, NA), simd_neg, c(TRUE, FALSE, NA))
  expect_tiers_give(abs(c(TRUE, FALSE, NA)), simd_abs, c(TRUE, FALSE, NA))
})

test_that("integer %/% and %% match base R, zero divisor NA without warning", {
  vals <- c(0L, 1L, -1L, 2L, -2L, 3L, -3L, 5L, -5L, 7L, -7L, imax, -imax, NA)
  p <- all_pairs(vals)
  expect_tier_warnings(character(0), simd_idiv, p$x, p$y)
  expect_tiers_give(p$x %/% p$y, simd_idiv, p$x, p$y)
  expect_tiers_give(p$x %% p$y, simd_mod, p$x, p$y)
  expect_identical(simd_idiv(-7L, 2L), -4L)
  expect_identical(simd_mod(-7L, 2L), 1L)
  expect_identical(simd_idiv(5L, 0L), NA_integer_)
  expect_identical(simd_mod(5L, 0L), NA_integer_)
  x <- rand_vec("integer", 1003, seed = 3L)
  y <- rand_vec("integer", 1003, seed = 4L) %/% 1000L
  expect_tiers_give(x %/% y, simd_idiv, x, y)
  expect_tiers_give(x %% y, simd_mod, x, y)
})

test_that("integer %/% and %% by a scalar match base R for every kind of divisor", {
  # A scalar divisor takes a division-free path (multiplication by a magic
  # number), which differs by divisor: small, powers of two and their
  # neighbours, +-1, the extremes.
  p2 <- as.integer(2^(1:30))
  divisors <- c(
    -9:-1, 1:9, p2, -p2, p2 - 1L, 1L - p2, p2 + 1L, -p2 - 1L,
    imax, -imax, imax - 1L, 641L, -6700417L, 1000000007L
  )
  x <- c(
    0L, 1L, -1L, imax, -imax, imax - 1L, -imax + 1L, NA,
    rand_vec("integer", 997, seed = 5L), seq(-3000L, 3000L, by = 7L)
  )
  for (d in divisors) {
    info <- paste("divisor", d)
    expect_tiers_give(x %/% d, simd_idiv, x, d)
    expect_tiers_give(x %% d, simd_mod, x, d)
    ok <- x[!is.na(x)]
    expect_identical(simd_idiv(ok, d, na_check = FALSE), ok %/% d, info = info)
    expect_identical(simd_mod(ok, d, na_check = FALSE), ok %% d, info = info)
  }
  expect_identical(simd_idiv(x, NA_integer_), x %/% NA_integer_)
  expect_identical(simd_mod(x, 0L), x %% 0L)
})

test_that("mul_add and add_mul follow base R's composition", {
  vals <- c(0L, 1L, -1L, 2L, 46341L, -46341L, 65536L, imax, -imax, NA)
  g <- expand.grid(x = vals, y = vals, z = c(0L, 1L, -1L, imax, -imax, NA))
  expected <- suppressWarnings(g$x * g$y + g$z)
  got <- expect_tier_warnings(ovf_msg, simd_mul_add, g$x, g$y, g$z)
  expect_identical(got, expected)
  expected <- suppressWarnings((g$x + g$y) * g$z)
  got <- expect_tier_warnings(ovf_msg, simd_add_mul, g$x, g$y, g$z)
  expect_identical(got, expected)
  # An overflowing intermediate is NA even if the final value would fit.
  expect_identical(suppressWarnings(simd_mul_add(65536L, 32768L, -imax)), NA_integer_)
  expect_identical(suppressWarnings(simd_add_mul(imax, 1L, 0L)), NA_integer_)
  # It stays NA when the caller turns NA checks off.
  expect_identical(
    suppressWarnings(simd_mul_add(65536L, 32768L, -imax, na_check = FALSE)), NA_integer_
  )
  expect_tier_warnings(character(0), simd_mul_add, c(NA, 2L), 3L, 4L)
})

test_that("results of every lgl/int pair are integer like base R", {
  vals <- list(logical = c(TRUE, FALSE, NA, TRUE), integer = c(3L, -2L, NA, imax))
  for (a in names(vals)) {
    for (b in names(vals)) {
      x <- vals[[a]]
      y <- vals[[b]]
      info <- paste(a, b)
      expect_identical(suppressWarnings(simd_add(x, y)), suppressWarnings(x + y), info = info)
      expect_identical(suppressWarnings(simd_mul(x, y)), suppressWarnings(x * y), info = info)
      expect_identical(simd_idiv(x, y), x %/% y, info = info)
      expect_identical(simd_div(x, y), x / y, info = info)
      expect_identical(simd_pmax(x, y), pmax(x, y), info = info)
    }
  }
})

test_that("broadcast works in each position", {
  x <- c(1L, -2L, NA, 4L, 5L)
  expect_tiers_give(x + 3L, simd_add, x, 3L)
  expect_tiers_give(3L - x, simd_sub, 3L, x)
  expect_tiers_give(x * 2L + 1L, simd_mul_add, x, 2L, 1L)
  expect_tiers_give(2L * x + 1L, simd_mul_add, 2L, x, 1L)
  expect_tiers_give(2L * 3L + x, simd_mul_add, 2L, 3L, x)
  expect_tiers_give(integer(0), simd_add, integer(0), 1L)
  expect_tiers_give(integer(0), simd_mul_add, 1L, integer(0), 2L)
  expect_tiers_give(7L, simd_add, 3L, 4L)
})

test_that("wrapping ops reject doubles, raw gives base R's errors", {
  expect_error(simd_add_wrap(1, 2L), "simd_add_wrap() does not support 'x' of type double",
    fixed = TRUE
  )
  expect_error(simd_mul_wrap(1L, 2), "'y' of type double", fixed = TRUE)
  expect_error(simd_abs_wrap(-1.5), "simd_abs_wrap() does not support 'x' of type double",
    fixed = TRUE
  )
  expect_error(simd_add(as.raw(1), as.raw(2)), "non-numeric argument to binary operator")
  expect_error(simd_mul_add(1L, as.raw(2), 3L), "non-numeric argument to binary operator")
  expect_error(simd_neg(as.raw(1)), "invalid argument to unary operator")
  expect_error(simd_abs(as.raw(1)), "non-numeric argument to mathematical function")
})

test_that("na_check = FALSE skips NA masks without harm", {
  x <- c(1L, 2L, 3L)
  expect_identical(simd_add(x, 1L, na_check = FALSE), x + 1L)
  expect_identical(simd_add_wrap(x, 1L, na_check = FALSE), x + 1L)
  expect_error(simd_add(x, 1L, na_check = NA), "'na_check' must be TRUE or FALSE")
  # NA operands with the check off: unspecified results, but no crash.
  expect_length(simd_add(c(NA, 1L), 1L, na_check = FALSE), 2L)
})

test_that("chunk boundaries and ALTREP inputs give the oracle's result", {
  n <- 2^20 + 7
  x <- seq_len(n)
  expect_true(takes_region_path(x))
  y <- rev(x)
  expect_simd_identical(simd_add, x, y)
  expect_identical(simd_add(x, y), rep(as.integer(n + 1), n))
  expect_identical(simd_mul(x, 3L), as.integer(x) * 3L)
  expect_identical(simd_idiv(x, 7L), x %/% 7L)
  expect_identical(simd_mod(1:5000, 7L), (1:5000) %% 7L)
  expect_identical(simd_sub_wrap(x, 1L), x - 1L)
  w <- warnings_of(r <- simd_add(x, imax - 2L))
  expect_identical(w, ovf_msg)
  expect_identical(r, suppressWarnings(x + (imax - 2L)))
})
