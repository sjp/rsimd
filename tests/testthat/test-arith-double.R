# Elementwise arithmetic on doubles: bit-identical to base R and across
# tiers, with NA_real_ wherever an operand is NA.

edge <- function() {
  c(
    0, -0, 1, -1, 0.5, -0.5, 2, -3, 1.5, 5.5, -5.5, 7, 1e16, -1e16, 1e300, 5e-324,
    .Machine$double.xmin, .Machine$double.xmax, Inf, -Inf, NA, NaN, pi
  )
}

test_that("add, sub, mul and div match base R bit for bit", {
  p <- all_pairs(edge())
  ops <- list(simd_add = `+`, simd_sub = `-`, simd_mul = `*`, simd_div = `/`)
  for (name in names(ops)) {
    expected <- na_merged(ops[[name]](p$x, p$y), p$x, p$y)
    expect_tiers_give(expected, get(name), p$x, p$y)
  }
  expect_tiers_give(c(NaN, NaN), simd_sub, c(Inf, 0), c(Inf, NaN))
  expect_tiers_give(NaN, simd_mul, 0, Inf)
})

test_that("an NA operand gives NA, whatever the other operand", {
  x <- c(NA, NaN, NA, 1, NaN)
  y <- c(NaN, NA, NA, NA, 1)
  for (f in list(simd_add, simd_sub, simd_mul, simd_div, simd_idiv, simd_mod, simd_copysign)) {
    expect_tiers_give(c(NA, NA, NA, NA, NaN), f, x, y)
  }
  expect_tiers_give(c(NA, NA, NA, NA, NaN), simd_fma, x, y, 1)
  expect_tiers_give(c(NA_real_, NA_real_), simd_lerp, c(1, NaN), 2, c(NA, NA))
})

test_that("integer operands are read as doubles, NA as NA", {
  x <- c(1L, NA, -3L, .Machine$integer.max)
  y <- c(0.5, 2, NaN, -Inf)
  expect_tiers_give(na_merged(x + y, x, y), simd_add, x, y)
  expect_tiers_give(na_merged(y - x, x, y), simd_sub, y, x)
  expect_tiers_give(x / 3L, simd_div, x, 3L)
  expect_tiers_give(c(TRUE, NA) * 2.5, simd_mul, c(TRUE, NA), 2.5)
  expect_tiers_give(5L / 0L, simd_div, 5L, 0L)
  expect_tiers_give(NaN, simd_div, 0L, 0L)
})

test_that("%/% and %% reproduce base R's special cases", {
  expect_tiers_give(c(Inf, -Inf, NaN, NaN), simd_idiv, c(5, -5, 0, NaN), 0)
  expect_tiers_give(c(-4, -4, -4, -4), simd_idiv, c(-7, 7, -7.5, 7.5), c(2, -2, 2, -2))
  expect_tiers_give(c(1, -1, 0.5, 1.5), simd_mod, c(-7, 7, -5.5, 5.5), c(2, -2, 2, 2))
  expect_tiers_give(
    c(NaN, 2, Inf, -Inf, -2, NaN, NaN), simd_mod, c(Inf, 2, -2, 2, -2, -Inf, Inf),
    c(2, Inf, Inf, -Inf, -Inf, -Inf, Inf)
  )
  expect_tiers_give(
    c(-1, 0, 0, -1, 0, Inf), simd_idiv, c(-2, 2, -2, 2, 0, Inf),
    c(Inf, Inf, -Inf, -Inf, Inf, 2)
  )
  expect_tiers_give(1, simd_mod, 1e20, 3)
  expect_tiers_give(NaN, simd_mod, 5, 0)
  # A zero remainder is +0, except that an infinite divisor returns x.
  expect_tiers_give(c(0, 0, 0, 0, -0), simd_mod, c(-0, 0, -4, 4, -0), c(2, -2, 2, -2, Inf))
  expect_tiers_give(c(0, 0), simd_idiv, c(-0, 1), c(5, 5))
})

test_that("%/% and %% are the exact floor and floored remainder", {
  # x / y rounds up to 3 here, so floor(x / y) would be one too big.
  x <- 3 + 2^-51
  y <- 1 + 2^-52
  expect_identical(x / y, 3)
  expect_tiers_give(2, simd_idiv, x, y)
  expect_tiers_give(1, simd_mod, x, y)
  expect_tiers_give(-3, simd_idiv, -x, y)
  expect_tiers_give(y - 1, simd_mod, -x, y)
  # The quotient underflows to -0, but the floor is -1.
  expect_tiers_give(c(-1, -1, 0), simd_idiv, c(-1e-300, -1e-300, 1e-300), c(1e300, 1, 1e300))
  # Quotients of 2^52 and more go through the scalar fallback.
  big <- c(1e20, 2^60 + 2^8, -(2^55), .Machine$double.xmax, 1e300)
  d <- c(3, 7, 0.1, 3, 1e-10)
  res <- expect_simd_identical(simd_mod, big, d)
  r <- res[["none"]]
  expect_true(all(r >= 0 & r < d))
  expect_identical(r[1:2], c(1, 5)) # 2^60 = 8^20 and 2^8 are 1 and 4 mod 7
  expect_simd_identical(simd_idiv, big, d)
})

test_that("%/% and %% match base R on random doubles", {
  x <- rand_vec("double", 3001, seed = 11L)
  y <- rand_vec("double", 3001, seed = 12L)
  q <- abs(x / y)
  keep <- !(is.finite(q) & q >= 2^52)
  x <- x[keep]
  y <- y[keep]
  expect_simd_identical(simd_idiv, x, y)
  expect_simd_identical(simd_mod, x, y)
  if (has_quad_long_double()) {
    expect_tiers_give(na_merged(x %/% y, x, y), simd_idiv, x, y)
    expect_tiers_give(na_merged(x %% y, x, y), simd_mod, x, y)
  } else {
    # Base R's long double (or double) correction is inexact here: compare
    # away from integer quotients, and remainders within rounding.
    qk <- q[keep]
    away <- !is.finite(qk) | abs(qk - round(qk)) > 1e-6
    expect_identical(simd_idiv(x, y)[away], na_merged(x %/% y, x, y)[away])
    expect_equal(simd_mod(x, y), na_merged(x %% y, x, y), tolerance = 1e-9)
  }
  i <- rand_vec("integer", 500, seed = 13L)
  j <- c(-7, 2.5, 0, 3)
  expect_tiers_give(na_merged(i %% 2.5, i), simd_mod, i, 2.5)
  expect_tiers_give(
    na_merged(rep(i, 4) %/% rep(j, each = 500), i), simd_idiv, rep(i, 4),
    rep(j, each = 500)
  )
})

test_that("fma is fused on every tier; mul_add and add_mul are not", {
  expect_tiers_give(2^-54, simd_fma, 0.1, 10, -1)
  expect_tiers_give(0, simd_mul_add, 0.1, 10, -1)
  expect_tiers_give((0.1 + 0.2) * 3, simd_add_mul, 0.1, 0.2, 3)
  x <- rand_vec("double", 1031, seed = 21L)
  y <- rand_vec("double", 1031, seed = 22L)
  z <- rand_vec("double", 1031, seed = 23L)
  expect_simd_identical(simd_fma, x, y, z)
  expect_tiers_give(na_merged(x * y + z, x, y, z), simd_mul_add, x, y, z)
  expect_tiers_give(na_merged((x + y) * z, x, y, z), simd_add_mul, x, y, z)
  # Within one rounding of the unfused result.
  f <- simd_fma(x, y, z)
  u <- x * y + z
  ok <- is.finite(u) & is.finite(f) & u != 0
  expect_true(all(abs(f[ok] - u[ok]) <= 2 * 2^-52 * (abs(x * y)[ok] + abs(z[ok]))))
  expect_tiers_give(c(5, 7), simd_fma, 1:2, 2L, 3L)
})

test_that("lerp is exact at t = 0 and t = 1", {
  x <- c(1e16, -3, 0.1, 5e-324, 1e308)
  y <- c(1, 7, 0.3, -1e308, -1e308)
  expect_tiers_give(x, simd_lerp, x, y, 0)
  expect_tiers_give(y, simd_lerp, x, y, 1)
  expect_tiers_give(c(5, 6.25), simd_lerp, c(0, 5), 10, c(0.5, 0.25))
  t <- rand_vec("double", 517, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = 31L) / 1e3
  expect_simd_identical(simd_lerp, rep(2, 517), 3, t)
  expect_tiers_give(5.5, simd_lerp, 5L, 6L, 0.5)
})

test_that("ternary ops broadcast in every position", {
  v <- c(1, 2.5, -4, NA, 8)
  expect_tiers_give(na_merged(v * 2 + 3, v), simd_mul_add, v, 2, 3)
  expect_tiers_give(na_merged(2 * v + 3, v), simd_mul_add, 2, v, 3)
  expect_tiers_give(na_merged(2 * 3 + v, v), simd_mul_add, 2, 3, v)
  expect_tiers_give(na_merged((v + 2) * 3, v), simd_add_mul, v, 2, 3)
  expect_tiers_give(numeric(0), simd_fma, numeric(0), 2, 3)
  expect_error(simd_fma(1:3, 1:2, 1), "lengths of 'x' (3), 'y' (2) and 'z' (1) must be equal or 1",
    fixed = TRUE
  )
  expect_error(simd_lerp(1:3, 1, 1:2), "'t' (2)", fixed = TRUE)
  expect_error(simd_add(1:3, c(1, 2)), "lengths of 'x' (3) and 'y' (2)", fixed = TRUE)
})

test_that("neg, abs, recip and sign match base R", {
  x <- edge()
  expect_tiers_give(-x, simd_neg, x)
  expect_tiers_give(abs(x), simd_abs, x)
  expect_tiers_give(1 / x, simd_recip, x)
  expect_tiers_give(sign(x), simd_sign, x)
  expect_identical(1 / simd_sign(-0), Inf)
  expect_tiers_give(sign(c(-2L, 0L, NA, 5L)), simd_sign, c(-2L, 0L, NA, 5L))
  expect_tiers_give(1 / c(-2L, 0L, NA), simd_recip, c(-2L, 0L, NA))
  expect_identical(1 / simd_neg(0), -Inf)
  expect_identical(1 / simd_abs(-0), Inf)
})

test_that("sqrt warns once when a number gives NaN", {
  x <- c(-1, 4, NaN, NA, -Inf, Inf, -0, 2)
  expected <- suppressWarnings(sqrt(x))
  got <- expect_tier_warnings("NaNs produced", simd_sqrt, x)
  expect_true(same_values(got, expected))
  expect_tier_warnings(character(0), simd_sqrt, c(NaN, NA, 4, -0))
  expect_tiers_give(sqrt(c(4, 9, 2)), simd_sqrt, c(4L, 9L, 2L))
  expect_tiers_give(NA_real_, simd_sqrt, NA_integer_)
})

test_that("copysign takes the sign bit, and a missing sign propagates", {
  x <- c(1, -2, 0, -0, Inf, 3, NaN, 5)
  s <- c(-1, 1, -0, 1, -5, NaN, -1, NA)
  got <- simd_copysign(x, s)
  expect_identical(got[1:5], c(-1, 2, -0, 0, -Inf))
  expect_identical(1 / got[3:4], c(-Inf, Inf))
  expect_true(is.nan(got[6]))
  expect_true(is.nan(got[7]))
  expect_true(is.na(got[8]) && !is.nan(got[8]))
  expect_simd_identical(simd_copysign, x, s)
  expect_tiers_give(c(-2, 3), simd_copysign, c(2L, -3L), c(-1, 1))
  expect_tiers_give(c(NA_real_, NA_real_), simd_copysign, c(NA, 1), c(1, NA))
})

test_that("double results agree across tiers on random vectors of edge lengths", {
  for (n in c(1, 2, 3, 4, 5, 7, 8, 9, 16, 17, 33, 4097)) {
    x <- rand_vec("double", n, seed = n + 40L)
    y <- rand_vec("double", n, seed = n + 41L)
    for (f in list(simd_add, simd_mul, simd_div, simd_idiv, simd_mod, simd_copysign)) {
      expect_simd_identical(f, x, y)
    }
    for (f in list(simd_neg, simd_sign, simd_recip, simd_floor, simd_round)) {
      expect_simd_identical(f, x)
    }
    expect_simd_identical(simd_fma, x, y, 1.5)
    expect_simd_identical(simd_lerp, x, y, 0.25)
  }
})

test_that("chunk boundaries and ALTREP inputs give the oracle's result", {
  n <- 2^20 + 7
  x <- as.double(seq_len(n))
  expect_simd_identical(simd_div, x, 3)
  expect_simd_identical(simd_mod, seq_len(n), 2.5)
  expect_identical(simd_div(seq_len(n), 4L), seq_len(n) / 4L)
  expect_identical(simd_fma(seq_len(n), 2, 0.5), seq_len(n) * 2 + 0.5)
  expect_identical(simd_sqrt(seq_len(5000)), sqrt(seq_len(5000)))
})

test_that("complex and integer64 operands are rejected where not taken", {
  expect_error(simd_idiv(1i, 1), "simd_idiv() does not support 'x' of type complex",
    fixed = TRUE
  )
  expect_error(simd_mod(1, 1i), "simd_mod() does not support 'y' of type complex",
    fixed = TRUE
  )
  expect_error(simd_fma(1, 1, 1i), "'z' of type complex", fixed = TRUE)
  expect_error(simd_sqrt(1i), "simd_sqrt() does not support 'x' of type complex",
    fixed = TRUE
  )
  x64 <- structure(0, class = "integer64")
  expect_error(simd_fma(x64, 1, 1), "simd_fma() does not support 'x' of type integer64",
    fixed = TRUE
  )
  expect_error(simd_sqrt(x64), "integer64$")
  expect_error(simd_add("a", 1), "'x' must be an atomic vector")
})

test_that("results are bare vectors", {
  x <- c(a = 1, b = 2)
  expect_null(attributes(simd_add(x, 1)))
  m <- matrix(1:4, 2)
  expect_null(attributes(simd_neg(m)))
  expect_null(attributes(simd_fma(m, 2, 1)))
})
