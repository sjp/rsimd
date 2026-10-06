# Accuracy of the elementary functions: every tier against the none tier
# (libm, as base R) within 2 ULP, and against base R within 2 ULP, over a
# grid, seeded random values and subnormals in each function's domain. The
# extended run uses a fine grid and 1e5 random values per function.

check_math1 <- function(name, spec, extended) {
  x <- math_inputs(spec[[3]], extended)
  none_bound <- if (is.null(spec$none)) 2 else spec$none
  err <- trig_ref_err(x, name)
  expect_tiers_close(spec[[1]], x, ulps = none_bound, label = name, x = x, abs = err)
  base_f <- spec[[2]]
  if (is.null(base_f) || (isTRUE(spec$platform) && !base_has_platform_sinpi())) {
    return(invisible())
  }
  want <- suppressWarnings(base_f(x))
  base_bound <- if (is.null(spec$base)) 2 else spec$base
  # R 4.6 returns +0 for f(-0) in these functions; rsimd returns -0.
  zeros <- !(name %in% r46_plus_zero && getRversion() >= "4.6.0")
  for (tier in tiers_to_test()) {
    got <- simd_with_impl(tier, suppressWarnings(spec[[1]](x)))
    expect_close_to(got, want, base_bound, paste(name, "vs base R on", tier),
      x = x, zeros = zeros, abs = err, abs_below = Inf
    )
  }
}

test_that("unary functions match none and base R within their ULP bounds", {
  tab <- math1_table()
  for (name in names(tab)) check_math1(name, tab[[name]], extended = FALSE)
})

test_that("unary functions are accurate on 1e5 random values (extended)", {
  skip_unless_extended()
  tab <- math1_table()
  for (name in names(tab)) check_math1(name, tab[[name]], extended = TRUE)
})

test_that("neon hands log, log2, cosh, asinh, acosh and pow to the C library", {
  # issue 036: there the C library was measured faster than SLEEF, so these
  # are the none tier's results bit for bit (asinh, acosh and pow, which
  # SLEEF has no fast variant of, in fast mode too). Not on Windows.
  skip_if_not("neon" %in% tiers_to_test(), "no neon tier")
  skip_on_os("windows")
  tab <- math1_table()
  same <- function(f, ..., acc) {
    old <- simd_math_accuracy(acc)
    on.exit(simd_math_accuracy(old))
    identical(
      simd_with_impl("neon", suppressWarnings(f(...))),
      simd_with_impl("none", suppressWarnings(f(...))),
      num.eq = FALSE
    )
  }
  for (name in c("log", "log2", "cosh", "asinh", "acosh")) {
    x <- c(math_inputs(tab[[name]][[3]]), NA, NaN, Inf, -Inf, 0, -0)
    expect_true(same(tab[[name]][[1]], x, acc = "accurate"), label = name)
    if (name %in% c("asinh", "acosh")) {
      expect_true(same(tab[[name]][[1]], x, acc = "fast"), label = paste(name, "fast"))
    }
  }
  x <- c(math_specials(), math_random(1500, -300, 300, TRUE))
  expect_true(same(simd_log, abs(x), 3, acc = "accurate"), label = "log base 3")
  for (acc in c("accurate", "fast")) {
    expect_true(same(simd_pow, abs(x) %% 30, x %% 8, acc = acc), label = paste("pow", acc))
    expect_true(same(simd_pow, 2L, x, acc = acc), label = paste("pow int base", acc))
  }
})

test_that("cbrt and tanpi are close to their base R expressions", {
  x <- c(math_random(2000, -300, 300, log_scale = TRUE), -8, 27, -0.125)
  ref <- sign(x) * abs(x)^(1 / 3)
  for (tier in tiers_to_test()) {
    got <- simd_with_impl(tier, simd_cbrt(x))
    # x^(1/3) itself is off by up to (1/3 - fl(1/3)) * |log(x)|, ~1.3e-14.
    expect_lte(max(abs(got - ref) / abs(ref)), 2e-14, label = paste("cbrt on", tier))
  }
  # Negative input gives a negative root (C's cbrt), unlike (-8)^(1/3).
  # glibc's cbrt is not exact at every cube (cbrt(27) is 3 + 4.4e-16 on
  # aarch64), so only -8 is checked exactly.
  expect_tiers_give(-2, simd_cbrt, -8)
  # Base R's tanpi is tan(pi * x) after reduction, which is inaccurate near
  # the poles; away from them they agree closely.
  t <- math_random(2000, -0.4, 0.4)
  for (tier in tiers_to_test()) {
    got <- simd_with_impl(tier, simd_tanpi(t))
    expect_lte(max(ulp_dist(got, tanpi(t))), 8, label = paste("tanpi on", tier))
  }
})

test_that("binary functions match none and base R", {
  x <- c(math_specials(), math_random(1500, -300, 300, TRUE))
  y <- c(rev(math_specials()), math_random(1500, -300, 300, TRUE, seed = 7L))
  p <- all_pairs(math_specials())
  expect_tiers_close(simd_atan2, x, y, label = "atan2", ulps = 2)
  expect_tiers_close(simd_atan2, p$x, p$y, label = "atan2 specials", ulps = 2)
  for (tier in tiers_to_test()) {
    got <- simd_with_impl(tier, simd_atan2(x, y))
    expect_close_to(got, atan2(x, y), 2, paste("atan2 vs base R on", tier))
  }
  # hypot: none uses libm's hypot; base R has none, so compare with
  # sqrt(x^2 + y^2) where that neither overflows nor underflows.
  expect_tiers_close(simd_hypot, x, y, label = "hypot", ulps = 1)
  expect_tiers_close(simd_hypot, p$x, p$y, label = "hypot specials", ulps = 1)
  hx <- math_random(2000, -100, 100, TRUE)
  hy <- math_random(2000, -100, 100, TRUE, seed = 3L)
  for (tier in tiers_to_test()) {
    got <- simd_with_impl(tier, simd_hypot(hx, hy))
    expect_close_to(got, sqrt(hx^2 + hy^2), 2, paste("hypot vs sqrt on", tier))
  }
  # pow: SLEEF within 1 ULP of the exact value; base R's ^ is libm's pow.
  px <- c(math_random(2000, 0, 50), -math_random(500, 0, 50, seed = 5L), 1e-300, 1e300)
  py <- c(math_random(2000, -20, 20, seed = 9L), round(math_random(500, -30, 30)), 0.5, -0.5)
  expect_tiers_close(simd_pow, px, py, label = "pow", ulps = 2)
  for (tier in tiers_to_test()) {
    got <- simd_with_impl(tier, simd_pow(px, py))
    expect_close_to(got, px^py, 2, paste("pow vs base R on", tier))
  }
})

test_that("binary functions are accurate on 1e5 random pairs (extended)", {
  skip_unless_extended()
  x <- math_random(1e5, -300, 300, TRUE)
  y <- math_random(1e5, -300, 300, TRUE, seed = 11L)
  expect_tiers_close(simd_atan2, x, y, label = "atan2", ulps = 2)
  expect_tiers_close(simd_hypot, x, y, label = "hypot", ulps = 1)
  px <- math_random(1e5, -3, 3, TRUE)
  py <- math_random(1e5, -40, 40, seed = 12L)
  expect_tiers_close(simd_pow, px, py, label = "pow", ulps = 2)
  for (tier in tiers_to_test()) {
    got <- simd_with_impl(tier, simd_pow(px, py))
    expect_close_to(got, px^py, 2, paste("pow vs base R on", tier))
  }
})

test_that("sincos gives sin and cos", {
  x <- math_inputs(dom(-1e3, 1e3))
  res <- with_each_tier(function() suppressWarnings(simd_sincos(x)))
  # SLEEF's sincos is its own algorithm, so the results may differ from
  # simd_sin() and simd_cos() in the last bit; both are within 1 ULP.
  sin_err <- trig_ref_err(x, "sin")
  cos_err <- trig_ref_err(x, "cos")
  for (tier in names(res)) {
    expect_named(res[[tier]], c("sin", "cos"))
    expect_close_to(res[[tier]]$sin, res[["none"]]$sin, 2, paste("sin on", tier),
      x = x, abs = sin_err, abs_below = Inf
    )
    expect_close_to(res[[tier]]$cos, res[["none"]]$cos, 2, paste("cos on", tier),
      x = x, abs = cos_err, abs_below = Inf
    )
    expect_close_to(res[[tier]]$sin, suppressWarnings(sin(x)), 2, paste("base sin on", tier),
      zeros = getRversion() < "4.6.0", abs = sin_err, abs_below = Inf
    )
    expect_close_to(res[[tier]]$cos, suppressWarnings(cos(x)), 2, paste("base cos on", tier),
      abs = cos_err, abs_below = Inf
    )
  }
  expect_identical(res[["none"]]$sin, suppressWarnings(sin(x)))
})

test_that("log2 and log10 of exact powers are exact; exp2 and exp10 of integers too", {
  k <- -1074:1023
  expect_tiers_give(as.double(k), simd_log2, 2^k)
  expect_tiers_give(as.double(0:22), simd_log10, 10^(0:22))
  expect_tiers_give(2^k, simd_exp2, as.double(k))
  expect_tiers_give(10^(0:22), simd_exp10, as.double(0:22))
  expect_tiers_give(2^(0:60), simd_pow, 2, as.double(0:60))
})

test_that("exp and its relatives produce subnormals and overflow where base R does", {
  x <- c(709.78, 709.79, -745, -745.13, -746, -708.4, -720, 710, 1000, -1000)
  ref <- expect_tiers_close(simd_exp, x, ulps = 2)
  expect_identical(ref, exp(x))
  expect_identical(is.infinite(simd_exp(x)), is.infinite(exp(x)))
  expect_tiers_give(c(5e-324, 0, Inf), simd_exp, c(-745, -746, 709.79))
  expect_tiers_give(c(5e-324, 0, Inf, 2^-1074), simd_exp2, c(-1074, -1076, 1024, -1074))
  expect_tiers_give(c(0, Inf), simd_exp10, c(-324, 309))
  expect_tiers_give(c(-1, Inf, 0, -0), simd_expm1, c(-Inf, Inf, 0, -0))
  expect_tiers_close(simd_expm1, c(709.7, 709.78, 709.79, -40, 1e-300), ulps = 2)
})

test_that("sin, cos and tan reduce large arguments correctly", {
  x <- c(1e22, 1e300, 2^53, -2^1000, 1e15 + 0.3, 8.98846567431158e307)
  # Where the C library does not reduce these accurately (Windows), the
  # none tier and base R are no reference (trig_ref_err); the exact values
  # below still are.
  fns <- list(list(simd_sin, sin, "sin"), list(simd_cos, cos, "cos"), list(simd_tan, tan, "tan"))
  for (f in fns) {
    err <- trig_ref_err(x, f[[3]])
    expect_tiers_close(f[[1]], x, ulps = 2, abs = err)
    for (tier in tiers_to_test()) {
      expect_close_to(simd_with_impl(tier, f[[1]](x)), f[[2]](x), 2, paste("on", tier),
        abs = err, abs_below = Inf
      )
    }
  }
  # The active tier's own values are exact unless it calls a C library that
  # does not reduce accurately (the none tier on Windows).
  if (tier_reduces_trig()) {
    expect_equal(simd_sin(1e22), -0.8522008497671888, tolerance = 1e-15)
    # 1e300 written exactly: where long double is double (arm64 macOS), R
    # parses the decimal literal to a neighbouring double.
    expect_equal(simd_cos(0x1.7e43c8800759cp+996), -0.5753861119575491, tolerance = 1e-15)
  }
})

test_that("hyperbolic functions are accurate up to overflow", {
  x <- c(709, 709.5, 709.9, 710.4, 710.47, 710.48, 710.5, 711, 1e3, -710.4)
  for (f in list(list(simd_sinh, sinh), list(simd_cosh, cosh))) {
    expect_tiers_close(f[[1]], x, ulps = 2)
    expect_close_to(f[[1]](x), f[[2]](x), 2, "vs base R")
  }
  big <- c(1e154, 1.34e154, 1e200, 1e300, .Machine$double.xmax)
  expect_tiers_close(simd_asinh, c(big, -big), ulps = 2)
  expect_tiers_close(simd_acosh, big, ulps = 2)
  expect_true(all(is.finite(simd_asinh(c(big, -big)))))
  expect_close_to(simd_acosh(big), acosh(big), 2, "acosh vs base R")
})
