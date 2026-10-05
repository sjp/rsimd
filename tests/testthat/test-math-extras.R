# The floating-point extras: ilogb, scaleb, rootn, nextafter, next_up,
# next_down, remainder and rsqrt (exact, or for rootn within 1 ULP), and the
# approximations recip_approx, rsqrt_approx and mul_add_approx.

# The exponent e with 2^e <= |x| < 2^(e + 1), computed exactly.
ilogb_ref <- function(x) {
  out <- rep(NA_integer_, length(x))
  ok <- !is.na(x) & is.finite(x) & x != 0
  a <- abs(x[ok])
  e <- floor(log2(a))
  e <- ifelse(2^e > a, e - 1, ifelse(2^(e + 1) <= a, e + 1, e))
  out[ok] <- as.integer(e)
  out
}

# Finite doubles over the whole range (subnormals too), both signs.
extras_doubles <- function(n = 2000, seed = 20261005L) {
  with_seed(seed, {
    e <- sample(-1074:1023, n, replace = TRUE)
    m <- 1 + stats::runif(n)
    v <- ifelse(e < -1022, 2^e, m * 2^e)
    v[!is.finite(v)] <- .Machine$double.xmax
    v * sample(c(-1, 1), n, replace = TRUE)
  })
}

test_that("simd_ilogb is the binary exponent, NA for zero, infinities and missing values", {
  x <- c(
    1, 2, 3, 0.5, 0.75, 1024, 1023, -8, 5e-324, 2^-1023, 2^-1022, 2.2e-308,
    .Machine$double.xmax, 0, -0, Inf, -Inf, NaN, NA, extras_doubles(500)
  )
  want <- ilogb_ref(x)
  expect_identical(want[1:19], c(
    0L, 1L, 1L, -1L, -1L, 10L, 9L, 3L, -1074L, -1023L, -1022L, -1023L, 1023L,
    rep(NA_integer_, 6)
  ))
  expect_tiers_give(want, simd_ilogb, x)
  expect_tiers_give(c(0L, 3L, NA, NA), simd_ilogb, c(1L, -10L, 0L, NA))
  expect_tiers_give(c(0L, NA, NA), simd_ilogb, c(TRUE, FALSE, NA))
  expect_tiers_give(integer(0), simd_ilogb, numeric(0))
  for (n in edge_lengths()) {
    v <- rep_len(c(3, -0.1, 0, 1e300), n)
    expect_tiers_give(ilogb_ref(v), simd_ilogb, v)
  }
})

test_that("simd_scaleb is x * 2^n exactly, rounded once, saturating", {
  x <- extras_doubles(1000)
  n <- with_seed(1L, sample(-60:60, 1000, replace = TRUE))
  # Within the normal range the product is exact.
  want <- x * 2^n
  ok <- abs(want) >= 2^-1022 & abs(want) < Inf
  expect_tiers_give(want[ok], simd_scaleb, x[ok], n[ok])
  cases <- list(
    list(1, 1024L, Inf), list(-1, 1024L, -Inf), list(1, 1023L, 2^1023),
    list(5e-324, 1074L, 1), list(1, -1074L, 5e-324), list(1, -1075L, 0),
    list(3 * 5e-324, -1L, 2 * 5e-324), list(5 * 5e-324, -1L, 2 * 5e-324),
    list(.Machine$double.xmax, -2098L, 5e-324), list(.Machine$double.xmax, -2099L, 0),
    list(1, .Machine$integer.max, Inf), list(1, -.Machine$integer.max, 0),
    list(-1, -.Machine$integer.max, -0), list(5e-324, .Machine$integer.max, Inf),
    list(0, 5000L, 0), list(-0, 5L, -0), list(Inf, -5000L, Inf), list(-Inf, 3L, -Inf),
    list(1.5, 0L, 1.5)
  )
  for (cs in cases) {
    expect_tiers_give(cs[[3]], simd_scaleb, cs[[1]], cs[[2]])
  }
  # Every tier as the none tier, over the whole range.
  big <- with_seed(2L, sample(-2200:2200, 1000, replace = TRUE))
  expect_tiers_identical(function() simd_scaleb(x, big), "scaleb")
  # Missing values.
  expect_tiers_give(c(NA, NaN, NA, NA), simd_scaleb, c(NA, NaN, 1, NaN), c(1L, 1L, NA, NA))
})

test_that("simd_scaleb and simd_rootn take whole numbers as n", {
  expect_identical(simd_scaleb(3, 2), 12)
  expect_identical(simd_scaleb(3, c(TRUE, FALSE)), c(6, 3))
  expect_identical(simd_rootn(8, 3), 2)
  expect_identical(simd_scaleb(1, NaN), NA_real_)
  expect_identical(simd_rootn(4, c(2, NA)), c(2, NA))
  expect_error(simd_scaleb(1, 1.5), "'n' must be whole numbers in the integer range")
  expect_error(simd_rootn(1, 2^31), "'n' must be whole numbers in the integer range")
  expect_error(simd_rootn(1, Inf), "'n' must be whole numbers")
  expect_error(simd_scaleb(1, "2"), "non-numeric argument")
  expect_error(simd_scaleb(1, 1i), "does not support 'n' of type complex")
  expect_error(simd_scaleb("1", 2L), "must be an atomic vector")
  expect_error(simd_scaleb(1:3, 1:2), "length")
  # Broadcast in either position.
  expect_identical(simd_scaleb(1:3, 1L), c(2, 4, 6))
  expect_identical(simd_scaleb(1, 1:3), c(2, 4, 8))
  if (requireNamespace("bit64", quietly = TRUE)) {
    expect_error(simd_rootn(8, bit64::as.integer64(3)), "does not support 'n' of type integer64")
  }
})

test_that("simd_rootn follows C23's special values", {
  x <- c(0, -0, 0, -0, 0, -0, 0, -0, Inf, -Inf, Inf, -Inf, -Inf, -Inf, -8, -8, 27, 2^60, 16, -16)
  n <- c(2L, 2L, 3L, 3L, -2L, -2L, -3L, -3L, 2L, 3L, -2L, -3L, 2L, -2L, 3L, -3L, 3L, 4L, 4L, 4L)
  want <- c(
    0, 0, 0, -0, Inf, Inf, Inf, -Inf, Inf, -Inf, 0, -0, NaN, NaN, -2, -0.5, 3, 2^15, 2, NaN
  )
  for_each_tier(function(tier) {
    expect_identical(zero_sign(suppressWarnings(simd_rootn(x, n))), zero_sign(want), info = tier)
  })
  expect_tiers_give(want, function(...) suppressWarnings(simd_rootn(...)), x, n)
  # n = 0 is NaN with a warning; n = 1 and -1 are x and 1 / x.
  quiet_rootn <- function(...) suppressWarnings(simd_rootn(...))
  expect_tiers_give(c(NaN, NaN, NaN), quiet_rootn, c(2, 0, Inf), 0L)
  expect_tiers_give(
    c(2, -0, Inf, 0.5, -Inf, 0), simd_rootn, c(2, -0, Inf, 2, -0, Inf),
    c(1L, 1L, 1L, -1L, -1L, -1L)
  )
  # Missing values come back unchanged, also with n = 0, and never warn.
  expect_tiers_give(
    c(NA, NaN, NA, NA, NaN), simd_rootn, c(NA, NaN, NA, 1, NaN), c(0L, 0L, 3L, NA, 2L)
  )
  for_each_tier(function(tier) {
    expect_warning(simd_rootn(-8, 2L), "NaNs produced", info = tier)
    expect_warning(simd_rootn(8, 0L), "NaNs produced", info = tier)
    expect_warning(simd_rootn(-Inf, 4L), "NaNs produced", info = tier)
    expect_no_warning(simd_rootn(c(NA, NaN, 8), c(2L, 0L, 3L)))
  })
})

test_that("simd_rootn is exact for exact powers", {
  b <- c(2:200, 1001, 3^10, 12345, 2^(8:30))
  for (n in c(3L, 4L, 5L, 6L, 7L, 10L, 13L, 20L, 31L, 60L)) {
    p <- b^n
    # Only exactly representable powers: those below 2^53, and powers of
    # two below 2^1000.
    keep <- p < 2^53 | (log2(b) == round(log2(b)) & p < 2^1000)
    expect_tiers_give(b[keep], simd_rootn, p[keep], n)
    if (n %% 2 == 1) expect_tiers_give(-b[keep], simd_rootn, -p[keep], n)
    # The reciprocal roots of powers of two are exact too.
    two <- keep & log2(b) == round(log2(b))
    expect_tiers_give(1 / b[two], simd_rootn, p[two], -n)
  }
  expect_tiers_give(c(3, 2^15, 10, 2), simd_rootn, c(27, 2^60, 1e10, 2^1000), c(3L, 4L, 10L, 1000L))
})

test_that("simd_rootn is within 2 ULP of the none tier on every tier, and close to x^(1/n)", {
  x <- abs(extras_doubles(1500))
  ns <- c(
    2L, 3L, 4L, 5L, 7L, 10L, 17L, 100L, 511L, 512L, 513L, 1000L, 123457L, .Machine$integer.max
  )
  for (n in c(ns, -ns)) {
    ref <- expect_tiers_close(simd_rootn, x, n, ulps = 2, label = paste("rootn n =", n))
    # pow(|x|, 1 / n) is off by up to about |log(x) / n| ULP.
    loose <- x^(1 / n)
    expect_close_to(ref, loose, ulps = 2 + 750 / abs(n), label = paste("rootn vs x^(1/n), n =", n))
  }
  # Odd n of negative x: minus the root of -x, on every tier.
  for (n in c(3L, 5L, -7L, 513L)) {
    neg <- with_each_tier(function() simd_rootn(-x, n))
    pos <- with_each_tier(function() simd_rootn(x, n))
    for (tier in names(neg)) expect_identical(neg[[tier]], -pos[[tier]], info = tier)
  }
  # n = 2 is sqrt exactly.
  expect_tiers_give(sqrt(x), simd_rootn, x, 2L)
})

test_that("simd_next_up, simd_next_down and simd_nextafter step one ULP", {
  x <- c(extras_doubles(1000), 1, -1, 2^-1022, -2^-1022, 2^-1022 - 5e-324)
  up <- simd_with_impl("none", simd_next_up(x))
  down <- simd_with_impl("none", simd_next_down(x))
  expect_true(all(up > x))
  expect_true(all(down < x))
  finite <- is.finite(up) & is.finite(down)
  expect_true(all(ulp_dist(x[finite], up[finite]) == 1))
  expect_true(all(ulp_dist(x[finite], down[finite]) == 1))
  expect_tiers_give(up, simd_next_up, x)
  expect_tiers_give(down, simd_next_down, x)
  big <- .Machine$double.xmax
  # The double below big, computed: R parses 1.7976931348623155e308 as Inf
  # where long double is double (arm64 macOS).
  below <- big - 2^971
  specials <- c(0, -0, 5e-324, -5e-324, big, -big, Inf, -Inf, NaN, NA)
  expect_tiers_give(
    c(5e-324, 5e-324, 1e-323, -0, Inf, -below, Inf, -big, NaN, NA),
    simd_next_up, specials
  )
  expect_tiers_give(
    c(-5e-324, -5e-324, 0, -1e-323, below, -Inf, big, -Inf, NaN, NA),
    simd_next_down, specials
  )
  expect_identical(zero_sign(simd_next_up(-5e-324)), -1)
  expect_identical(zero_sign(simd_next_down(5e-324)), 1)
  # Integer and logical input is read as double.
  expect_tiers_give(c(1 + 2^-52, NA, 5e-324), simd_next_up, c(1L, NA, 0L))
  expect_tiers_give(c(1 - 2^-53, -5e-324), simd_next_down, c(TRUE, FALSE))
  # nextafter: towards y, y itself when equal (so the sign of y's zero).
  y <- c(rep(Inf, 5), rep(-Inf, 5))
  xx <- c(x[1:5], x[1:5])
  expect_tiers_give(c(up[1:5], down[1:5]), simd_nextafter, xx, y)
  expect_tiers_give(
    c(-0, 0, 5e-324, -5e-324, 2, 1), simd_nextafter, c(0, -0, 0, -0, 2, 1), c(-0, 0, 1, -1, 2, 1)
  )
  expect_identical(zero_sign(simd_nextafter(c(0, -0), c(-0, 0))), c(-1, 1))
  expect_tiers_give(c(NA, NaN, NA, NA), simd_nextafter, c(NA, NaN, 1, NaN), c(1, 1, NA, NA))
  expect_tiers_identical(function() simd_nextafter(x, rev(x)), "nextafter")
})

test_that("simd_remainder is the IEEE remainder, exact on every tier", {
  expect_tiers_give(
    c(1, -1, -1, 1, -0, 0, 0.5, -0.5, 1, 5),
    simd_remainder, c(5, 7, -5, -7, -4, 4, 2.5, -2.5, 1, 5), c(2, 2, 2, -2, 2, 2, 1, 1, 2, Inf)
  )
  expect_identical(zero_sign(simd_remainder(c(-4, 4), 2)), c(-1, 1))
  # Unlike %%, the quotient is rounded to nearest (ties to even).
  expect_identical(simd_remainder(5, 3), -1)
  expect_identical(5 %% 3, 2)
  x <- extras_doubles(1500)
  y <- extras_doubles(1500, seed = 7L)
  r <- simd_with_impl("none", simd_remainder(x, y))
  expect_true(all(abs(r) <= abs(y) / 2))
  expect_tiers_identical(function() simd_remainder(x, y), "remainder")
  # Huge quotients, beyond the range of a double.
  big <- c(.Machine$double.xmax, 1e300, -2^1000, 3^600)
  tiny <- c(5e-324, 3e-310, 1e-300, -7e-200)
  expect_tiers_identical(
    function() simd_remainder(rep(big, each = 4), rep(tiny, 4)), "remainder huge"
  )
  expect_tiers_identical(
    function() simd_remainder(rep(tiny, each = 4), rep(big, 4)), "remainder tiny"
  )
  # NaN results warn; missing values do not.
  for_each_tier(function(tier) {
    expect_warning(simd_remainder(1, 0), "NaNs produced", info = tier)
    expect_warning(simd_remainder(Inf, 1), "NaNs produced", info = tier)
    got <- suppressWarnings(simd_remainder(c(1, Inf, -Inf), c(0, 1, 0)))
    expect_identical(got, c(NaN, NaN, NaN))
    expect_no_warning(simd_remainder(c(NA, NaN, 1), c(1, 1, NA)))
  })
  expect_tiers_give(c(NA, NaN, NA), simd_remainder, c(NA, NaN, 1), c(1, 1, NA))
})

test_that("simd_rsqrt is base R's 1 / sqrt(x) on every tier", {
  x <- c(abs(extras_doubles(1000)), 0, -0, Inf, -1, -Inf, NaN, NA, 4, 2, 1e-310)
  want <- suppressWarnings(1 / sqrt(x))
  expect_tiers_give(want, function(v) suppressWarnings(simd_rsqrt(v)), x)
  for_each_tier(function(tier) {
    expect_warning(simd_rsqrt(-1), "NaNs produced", info = tier)
    expect_no_warning(simd_rsqrt(c(NA, NaN, 0, -0)))
  })
  expect_identical(simd_rsqrt(c(4L, NA)), c(0.5, NA))
})

test_that("the approximations are within 2^-22 of the exact value on every tier", {
  x <- extras_doubles(4000, seed = 3L)
  x <- x[abs(x) >= 2^-1022 & abs(x) < 2^1022]
  bound <- 2^-22
  for_each_tier(function(tier) {
    r <- simd_recip_approx(x)
    expect_true(max(abs(r * x - 1)) <= bound, info = tier)
    s <- simd_rsqrt_approx(abs(x))
    expect_true(max(abs(s * sqrt(abs(x)) - 1)) <= bound * 1.01, info = tier)
  })
  # On none they are exact.
  expect_identical(simd_with_impl("none", simd_recip_approx(x)), 1 / x)
  expect_identical(simd_with_impl("none", simd_rsqrt_approx(abs(x))), 1 / sqrt(abs(x)))
})

test_that("the approximations are exact for special values and keep missing values", {
  sp <- c(0, -0, Inf, -Inf, NaN, NA, 5e-324, -1e-310, 2^1022, -.Machine$double.xmax, 2^-1023)
  expect_tiers_give(1 / sp, simd_recip_approx, sp)
  expect_identical(zero_sign(simd_recip_approx(c(Inf, -Inf))), c(1, -1))
  pos <- c(0, -0, Inf, NaN, NA, 5e-324, 1e-310, 2^-1023)
  expect_tiers_give(1 / sqrt(pos), simd_rsqrt_approx, pos)
  for_each_tier(function(tier) {
    expect_warning(simd_rsqrt_approx(c(4, -1)), "NaNs produced", info = tier)
    expect_identical(suppressWarnings(simd_rsqrt_approx(c(-1, -Inf, -5e-324))), c(NaN, NaN, NaN))
    expect_no_warning(simd_recip_approx(c(0, NA, NaN, Inf)))
  })
  # 2 is in range, so only within the bound: the hardware estimates (RCPPS,
  # FRECPE) need not be exact even for powers of two.
  for_each_tier(function(tier) {
    r <- simd_recip_approx(c(2L, NA, 0L))
    expect_identical(r[2:3], c(NA, Inf), info = tier)
    expect_true(abs(r[1] * 2 - 1) <= 2^-22, info = tier)
  })
  # Every length, so every tail.
  for (n in edge_lengths()) {
    v <- rep_len(c(3, 1e-308, 7, -0.25), n)
    for_each_tier(function(tier) {
      expect_true(all(abs(simd_recip_approx(v) * v - 1) <= 2^-22), info = paste(tier, n))
    })
  }
})

test_that("the extras ignore rsimd.math_accuracy", {
  x <- extras_doubles(500)
  for_each_tier(function(tier) {
    old <- simd_math_accuracy("fast")
    on.exit(simd_math_accuracy(old), add = TRUE)
    fast <- list(
      simd_recip_approx(x), simd_rsqrt_approx(abs(x)), simd_rsqrt(abs(x)), simd_rootn(abs(x), 3L),
      simd_remainder(x, rev(x)), simd_nextafter(x, 0)
    )
    simd_math_accuracy("accurate")
    acc <- list(
      simd_recip_approx(x), simd_rsqrt_approx(abs(x)), simd_rsqrt(abs(x)), simd_rootn(abs(x), 3L),
      simd_remainder(x, rev(x)), simd_nextafter(x, 0)
    )
    expect_identical(fast, acc, info = tier)
  })
})

test_that("simd_mul_add_approx is fused or unfused, and exact for integers", {
  x <- extras_doubles(500) / 1e300
  y <- extras_doubles(500, seed = 9L) / 1e300
  z <- -x * y
  fused <- simd_with_impl("none", simd_fma(x, y, z))
  plain <- simd_with_impl("none", simd_mul_add(x, y, z))
  expect_identical(simd_with_impl("none", simd_mul_add_approx(x, y, z)), plain)
  for_each_tier(function(tier) {
    got <- simd_mul_add_approx(x, y, z)
    same <- identical(got, fused) || identical(got, plain)
    expect_true(same, info = tier)
  })
  expect_identical(simd_mul_add_approx(0.1, 10, NA), NA_real_)
  # Integers: simd_mul_add's exact, checked result.
  for_each_tier(function(tier) {
    expect_identical(simd_mul_add_approx(1:3, 2L, 1L), c(3L, 5L, 7L), info = tier)
    expect_warning(
      r <- simd_mul_add_approx(.Machine$integer.max, 2L, -5L), "NAs produced by integer overflow"
    )
    expect_identical(r, NA_integer_)
    expect_identical(simd_mul_add_approx(TRUE, 3L, NA), NA_integer_)
  })
  if (requireNamespace("bit64", quietly = TRUE)) {
    a <- bit64::as.integer64(c(3, -4, NA))
    expect_identical(simd_mul_add_approx(a, 5L, 1L), simd_mul_add(a, 5L, 1L))
  }
})

test_that("the extras reject complex and integer64 input", {
  unary <- list(
    simd_ilogb, simd_next_up, simd_next_down, simd_rsqrt, simd_recip_approx, simd_rsqrt_approx
  )
  for (f in unary) {
    expect_error(f(1i), "does not support 'x' of type complex")
    expect_error(f("a"), "must be an atomic vector")
  }
  expect_error(simd_nextafter(1, 1i), "does not support 'y' of type complex")
  expect_error(simd_remainder(1i, 1), "does not support 'x' of type complex")
  if (requireNamespace("bit64", quietly = TRUE)) {
    expect_error(simd_ilogb(bit64::as.integer64(1)), "does not support 'x' of type integer64")
    expect_error(simd_scaleb(bit64::as.integer64(1), 1L), "does not support 'x' of type integer64")
  }
})

test_that("the extras honour a simd_vec's pin and re-wrap", {
  v <- simd_vec(c(8, 27, -0.5), impl = "none")
  for (f in list(
    simd_next_up, simd_next_down, simd_rsqrt, simd_recip_approx, simd_rsqrt_approx,
    function(x) simd_rootn(x, 3L), function(x) simd_scaleb(x, 2L),
    function(x) simd_nextafter(x, 0), function(x) simd_remainder(x, 2),
    function(x) simd_mul_add_approx(x, 2, 1)
  )) {
    r <- suppressWarnings(f(v))
    expect_true(is_simd_vec(r))
    expect_identical(simd_impl(r), "none")
  }
  i <- simd_ilogb(v)
  expect_true(is_simd_vec(i))
  expect_identical(as.vector(i), c(3L, 4L, -1L))
  # A pinned n counts too.
  expect_identical(simd_impl(simd_scaleb(1, simd_vec(2L, impl = "none"))), "none")
})
