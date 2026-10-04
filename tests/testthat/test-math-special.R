# Special values, the pi functions' exact cases and simd_pow()'s base R
# rules.

test_that("special values give the same kind of result as base R on every tier", {
  x <- math_specials()
  tab <- math1_table()
  refs <- list(cbrt = function(x) ifelse(is.na(x) | x == 0 | is.infinite(x), x, sign(x)))
  for (name in names(tab)) {
    spec <- tab[[name]]
    base_f <- if (is.null(spec[[2]])) refs[[name]] else spec[[2]]
    # Base R's sinpi() zero signs depend on the platform: glibc >= 2.41's
    # sinpi gives the sign of x (as simd_sinpi does, tested below), R's own
    # fallback gives +0 at every integer.
    zeros <- !(name %in% r46_plus_zero && getRversion() >= "4.6.0") && name != "sinpi"
    res <- with_each_tier(function() suppressWarnings(spec[[1]](x)))
    for (tier in names(res)) {
      got <- res[[tier]]
      info <- paste(name, "on", tier)
      # Every tier exactly as none for the kinds, infinities and zeros.
      bound <- if (name == "cbrt") 4 else 2
      expect_null(math_mismatch(got, res[["none"]], bound, x), label = info)
      if (is.null(base_f) || name == "tanpi") next
      want <- suppressWarnings(base_f(x))
      expect_identical(is.na(got), is.na(want), info = info)
      expect_identical(is.nan(got), is.nan(want), info = info)
      expect_identical(is.infinite(got), is.infinite(want), info = info)
      if (zeros) {
        expect_identical(zero_sign(got)[want %in% 0], zero_sign(want)[want %in% 0], info = info)
      }
    }
  }
})

test_that("the odd functions keep the sign of zero, as C does", {
  for (name in c(r46_plus_zero, "cbrt", "sinpi")) {
    f <- get(paste0("simd_", name))
    expect_tiers_give(c(-0, 0), f, c(-0, 0))
    expect_identical(zero_sign(f(c(-0, 0))), c(-1, 1), info = name)
  }
  expect_tiers_give(-pi, simd_atan2, -0, -0)
  expect_tiers_give(pi, simd_atan2, 0, -0)
  expect_tiers_give(c(-0, 0), simd_atan2, c(-0, 0), c(0, 0))
  expect_tiers_give(-Inf, simd_log, c(-0))
})

test_that("sinpi, cospi and tanpi are exact at the integers and half-integers", {
  h <- (-200:200) / 2
  n <- h[h == round(h)]
  half <- h[h != round(h)]
  # sinpi: a zero with the sign of x at the integers, +-1 at the halves.
  expect_tiers_give(ifelse(n < 0 | (n == 0 & 1 / n < 0), -0, 0), simd_sinpi, n)
  expect_tiers_give(ifelse((half - 0.5) %% 2 == 0, 1, -1), simd_sinpi, half)
  # cospi: +-1 at the integers, +0 at the halves.
  expect_tiers_give(ifelse(n %% 2 == 0, 1, -1), simd_cospi, n)
  expect_tiers_give(rep(0, length(half)), simd_cospi, half)
  # tanpi: +0 at every integer (also -0 and the negative ones), +-1 at the
  # quarters, NaN at the halves.
  expect_tiers_give(rep(0, length(n)), simd_tanpi, n)
  expect_tiers_give(0, simd_tanpi, -0)
  q <- (-40:40) / 4
  q <- q[q * 2 != round(q * 2)]
  expect_tiers_give(ifelse(q %% 1 == 0.25, 1, -1), simd_tanpi, q)
  expect_tier_warnings("NaNs produced", simd_tanpi, half)
  expect_true(all(is.nan(suppressWarnings(simd_tanpi(half)))))
  # The same values as base R where base R is exact too.
  expect_identical(simd_cospi(h), cospi(h))
  expect_identical(simd_tanpi(c(n, q)), tanpi(c(n, q)))
  if (base_has_platform_sinpi()) {
    expect_identical(1 / simd_sinpi(h), 1 / sinpi(h))
  }
})

test_that("the pi functions reduce huge arguments exactly", {
  # 2^53 + 1 rounds to 2^53; 2^52 + 1 and 1e9 + 1 are odd.
  x <- c(
    1e300, -1e300, 2^53 + 2, 2^53 + 1, 1e17, -1e17, 2^52 + 1, 2^51 + 0.5,
    2.5e8 + 0.5, 1e9 + 1, -(2^51 + 0.5)
  )
  expect_tiers_give(c(0, -0, 0, 0, 0, -0, 0, 1, 1, 0, -1), simd_sinpi, x)
  expect_tiers_give(c(1, 1, 1, 1, 1, 1, -1, 0, 0, -1, 0), simd_cospi, x)
  expect_tiers_give(
    c(0, 0, 0, 0, 0, 0, 0, NaN, NaN, 0, NaN), function(x) suppressWarnings(simd_tanpi(x)), x
  )
  # Beyond 2.5e8, where SLEEF's own sinpi and cospi give 0 and 1.
  q <- 2.5e8 + c(0.25, 0.75, 1.125, -0.375) + 2^40
  for (f in list(simd_sinpi, simd_cospi, simd_tanpi)) {
    expect_close_to(f(q), f(q - 2^40 - 2.5e8), 0, "reduction")
    expect_tiers_close(f, q, ulps = 2)
  }
  expect_close_to(simd_sinpi(2.5e8 + 0.25), sqrt(0.5), 1, "sinpi(1/4)")
  expect_tiers_give(c(0, 1, 0), simd_sinpi, c(1e300, 0.5, 2^53 + 1))
  expect_tiers_give(1, simd_cospi, 1e17)
  expect_identical(simd_sinpi(c(1e300, 2^53 + 1)), sinpi(c(1e300, 2^53 + 1)))
})

test_that("the pi functions give NaN with a warning for infinities", {
  for (f in list(simd_sinpi, simd_cospi, simd_tanpi)) {
    expect_tier_warnings("NaNs produced", f, c(Inf, 1))
    expect_tiers_give(c(NaN, NaN), function(x) suppressWarnings(f(x)), c(Inf, -Inf))
    expect_tier_warnings(character(), f, c(NaN, NA))
    expect_tiers_give(c(NaN, NA), f, c(NaN, NA))
  }
})

test_that("simd_pow() is bit-identical to base R's ^ at its special cases", {
  v <- c(
    0, -0, 1, -1, 2, -2, 3, 0.5, -0.5, 1 / 3, -8, 11, -11, 11.5, 12, 4, 10, 1e300, 1e-300,
    Inf, -Inf, NaN, NA
  )
  p <- all_pairs(v)
  # Missing values by the package's rule: 1 when x is 1 or y is 0,
  # otherwise NA when either operand is NA. (Base R warns about its own
  # fmod for (-Inf)^1e300.)
  want <- na_merged(suppressWarnings(p$x^p$y), p$x, p$y)
  want[p$x %in% 1 | p$y %in% 0] <- 1
  # Bit-identical where base R has a special case; elsewhere both call a
  # pow within 1 ULP.
  special <- is.na(p$x) | is.na(p$y) | !is.finite(p$x) | !is.finite(p$y) | p$x %in% c(0, 1) |
    p$y %in% c(0, 2) | (p$y %in% c(3, 4) & abs(p$x) <= 11)
  expect_tiers_give(want[special], simd_pow, p$x[special], p$y[special])
  expect_tiers_close(simd_pow, p$x, p$y, ulps = 2)
  for (tier in tiers_to_test()) {
    expect_close_to(simd_with_impl(tier, simd_pow(p$x, p$y)), want, 2, paste("pow on", tier))
  }
  # Signed zeros and infinities of R's rules.
  expect_tiers_give(
    c(0, Inf, 0, Inf, NaN, NaN, -Inf, Inf, 0, -8, 4), simd_pow,
    c(-0, -0, -0, -0, -Inf, -2, -Inf, -Inf, -Inf, -2, -2),
    c(3, -1, 0.5, -2, 0.5, Inf, 3, 2, -3, 3, 2)
  )
  # x^2, x^3 and x^4 for |x| <= 11 are products, as base R computes them.
  x <- math_random(500, -11, 11)
  for (k in 2:4) expect_tiers_give(x^k, simd_pow, x, k)
  y <- math_random(500, -11, 11, seed = 2L)
  expect_tiers_give(y * y, simd_pow, y, 2)
  expect_tiers_give(c(1, 1, 1, 1), simd_pow, c(NA, NaN, Inf, 1), c(0, 0, 0, NA))
})

test_that("simd_pow() never warns, like ^", {
  expect_tier_warnings(character(), simd_pow, c(-8, -1, -Inf), c(1 / 3, 0.5, 0.5))
  expect_true(all(is.nan(simd_pow(c(-8, -1), c(1 / 3, 0.5)))))
  expect_tiers_give(NaN, simd_pow, -8, 1 / 3)
})

test_that("simd_pow() takes integer and logical operands like ^", {
  x <- c(2L, -3L, NA, 0L, 1L, 10L)
  y <- c(3L, 2L, 0L, -1L, NA, 22L)
  expect_tiers_give(x^y, simd_pow, x, y)
  expect_tiers_give(c(TRUE, FALSE, NA)^2L, simd_pow, c(TRUE, FALSE, NA), 2L)
  expect_tiers_give(2^(0:30), simd_pow, 2L, 0:30)
  expect_identical(typeof(simd_pow(2L, 3L)), "double")
})
