# The elementary functions' interface: the "NaNs produced" warning, missing
# values, argument types, lengths, broadcasting and simd_log()'s base.

math1_names <- function() names(math1_table())

test_that("one 'NaNs produced' warning per call, as base R", {
  expect_tier_warnings("NaNs produced", simd_log, c(-1, -2, -3))
  expect_tier_warnings("NaNs produced", simd_sqrt, c(-1, -2))
  cases <- list(
    log = -1, log2 = -1, log10 = -1, log1p = -2, asin = 2, acos = -2, acosh = 0.5,
    atanh = 2, sin = Inf, cos = -Inf, tan = Inf, sinpi = Inf, cospi = Inf, tanpi = 0.5
  )
  for (name in names(cases)) {
    f <- get(paste0("simd_", name))
    x <- rep(cases[[name]], 9)
    w <- expect_tier_warnings("NaNs produced", f, x)
    expect_true(all(is.nan(w) & !is.na(w) | is.nan(w)), info = name)
    expect_identical(suppressWarnings(f(x)), suppressWarnings(get(name)(x)), info = name)
  }
  expect_tier_warnings("NaNs produced", simd_sincos, c(1, Inf))
  expect_tier_warnings("NaNs produced", simd_log, -8, 2)
  expect_tier_warnings("NaNs produced", simd_log, 1, 1)
})

test_that("NA and NaN inputs propagate silently and keep their kind", {
  x <- c(NaN, NA, 1, NA, NaN)
  for (name in math1_names()) {
    f <- get(paste0("simd_", name))
    res <- expect_tier_warnings(character(), f, x)
    expect_identical(is.na(res[c(1, 2, 4, 5)]), rep(TRUE, 4), info = name)
    expect_identical(is.nan(res[c(1, 2, 4, 5)]), c(TRUE, FALSE, FALSE, TRUE), info = name)
  }
  sc <- expect_tier_warnings(character(), simd_sincos, x)
  expect_identical(is.nan(sc$sin), is.nan(x))
  expect_identical(is.nan(sc$cos), is.nan(x))
  # Two operands: NA if either is NA, else NaN if either is NaN (base R's
  # atan2), no warning.
  a <- c(NA, NaN, NA, 1, NaN, 1)
  b <- c(NaN, NA, NA, NA, 1, NaN)
  want <- c(NA, NA, NA, NA, NaN, NaN)
  expect_tiers_give(want, simd_atan2, a, b)
  expect_identical(atan2(a, b), want)
  expect_tiers_give(want, simd_hypot, a, b)
  expect_tier_warnings(character(), simd_hypot, a, b)
  expect_tiers_give(c(NaN, NaN, NA), simd_hypot, c(Inf, NaN, NA), c(NaN, Inf, Inf))
  expect_tiers_give(c(NA, NA, 1, 1, NaN), simd_pow, c(NA, 2, NA, 1, NaN), c(2, NA, 0, NA, 3))
})

test_that("an integer NA becomes NA_real_, and integer or logical input equals the double path", {
  x <- c(-3L, 0L, 1L, 2L, NA, 100L, .Machine$integer.max)
  for (name in math1_names()) {
    f <- get(paste0("simd_", name))
    res <- with_each_tier(function() suppressWarnings(f(x)))
    for (tier in names(res)) {
      expect_identical(res[[tier]], simd_with_impl(tier, suppressWarnings(f(as.double(x)))),
        info = paste(name, tier)
      )
    }
    expect_identical(suppressWarnings(f(c(TRUE, FALSE, NA))), suppressWarnings(f(c(1, 0, NA))),
      info = name
    )
  }
  # The double path, not base exp(): a tier's exp() may be 1 ulp from libm's.
  expect_identical(simd_exp(1L), simd_exp(1))
  expect_identical(simd_exp(TRUE), simd_exp(1))
  expect_identical(simd_exp(NA_integer_), NA_real_)
  expect_identical(simd_atan2(1L, c(2, NA)), simd_atan2(1, c(2, NA)))
  expect_identical(simd_hypot(3L, 4L), 5)
  expect_identical(simd_sincos(1:3), simd_sincos(c(1, 2, 3)))
})

test_that("unsupported argument types error", {
  for (name in c("exp", "log", "sinpi", "atanh")) {
    f <- get(paste0("simd_", name))
    if (name %in% c("sinpi")) {
      expect_error(f(1i), sprintf("simd_%s() does not support 'x' of type complex", name),
        fixed = TRUE
      )
    }
    # As every function of the package: the access layer names the type.
    expect_error(f("a"), "'x' must be an atomic vector", fixed = TRUE)
    expect_error(f(list(1)), "'x' must be an atomic vector", fixed = TRUE)
    expect_error(f(as.raw(1)), "non-numeric argument to mathematical function", fixed = TRUE)
  }
  # Complex input is taken where base R takes it (test-complex-math.R).
  expect_error(simd_hypot(2, 1i), "simd_hypot() does not support 'y' of type complex",
    fixed = TRUE
  )
  expect_error(simd_nextafter(1i, 1), "simd_nextafter() does not support 'x' of type complex",
    fixed = TRUE
  )
  expect_error(simd_hypot(1, as.raw(1)), "non-numeric argument to mathematical function",
    fixed = TRUE
  )
  expect_error(simd_pow(as.raw(1), 1), "non-numeric argument to binary operator", fixed = TRUE)
  expect_error(simd_pow("a", 1), "'x' must be an atomic vector", fixed = TRUE)
  expect_error(simd_sincos(1i), "simd_sincos() does not support 'x' of type complex",
    fixed = TRUE
  )
  expect_error(simd_log1p(1i), "simd_log1p() does not support 'x' of type complex", fixed = TRUE)
  skip_if_not_installed("bit64")
  expect_error(simd_exp(bit64::as.integer64(1)),
    "simd_exp() does not support 'x' of type integer64",
    fixed = TRUE
  )
})

test_that("results are bare double vectors of every length, chunked and ALTREP inputs too", {
  for (n in sweep_lengths()) {
    x <- math_random(n, -5, 5, seed = n %% 1000)
    for (f in list(simd_exp, simd_tanh, simd_sinpi)) {
      res <- expect_tiers_close(f, x, ulps = 2)
      expect_length(res, n)
    }
    expect_tiers_close(simd_atan2, x, 0.5, ulps = 2)
    expect_tiers_close(simd_pow, 1.5, x, ulps = 2)
  }
  for (x in altrep_inputs(5000)) {
    expect_true(takes_region_path(x))
    expect_identical(simd_log(x), simd_log(as.double(x)))
    expect_identical(simd_pow(x, 0.5), simd_pow(as.double(x), 0.5))
  }
  x <- c(a = 1, b = 2)
  dim(x) <- NULL
  named <- c(a = 0.5, b = 2)
  expect_null(attributes(simd_exp(named)))
  expect_null(attributes(simd_atan2(named, 1)))
  m <- matrix(1:4, 2)
  expect_null(attributes(simd_sin(m)))
  expect_identical(simd_exp(double()), double())
  expect_identical(simd_pow(double(), 2), double())
  expect_identical(simd_sincos(double()), list(sin = double(), cos = double()))
})

test_that("binary functions broadcast a length-1 operand and reject other mismatches", {
  x <- c(-2, -0.5, 0, 0.5, 2, 8)
  expect_tiers_give(x^3, simd_pow, x, 3)
  expect_tiers_close(simd_pow, 2, x, ulps = 2)
  expect_close_to(simd_pow(2, x), 2^x, 1, "pow vs base R")
  expect_tiers_close(simd_atan2, x, 1)
  expect_identical(simd_atan2(1, x), atan2(1, x))
  expect_identical(simd_hypot(x, 0), abs(x))
  expect_error(simd_pow(1:3, 1:2), "lengths of 'x' (3) and 'y' (2) must be equal or one of them must be 1",
    fixed = TRUE
  )
  expect_error(simd_atan2(1:3, 1:2), "lengths of 'y' (3) and 'x' (2) must be equal or one of them must be 1",
    fixed = TRUE
  )
  expect_identical(simd_hypot(1:2, double()), double())
})

test_that("simd_log() takes a base like base R's log()", {
  x <- c(8, 100, 1000, 0.1, 0, -1, Inf, NA, NaN, 3^(0:20))
  for (b in c(2, 10, 3, exp(1), 0.5, 1e-300)) {
    res <- with_each_tier(function() suppressWarnings(simd_log(x, b)))
    expect_tier_warnings("NaNs produced", simd_log, x, b)
    want <- suppressWarnings(log(x, b))
    for (tier in names(res)) {
      expect_close_to(res[[tier]], want, 2, paste("base", b, "on", tier))
    }
  }
  expect_identical(suppressWarnings(simd_log(x, 2)), suppressWarnings(simd_log2(x)))
  expect_identical(suppressWarnings(simd_log(x, 10)), suppressWarnings(simd_log10(x)))
  expect_identical(suppressWarnings(simd_log(x, 2L)), suppressWarnings(simd_log2(x)))
  expect_identical(suppressWarnings(simd_log(x)), suppressWarnings(log(x)))
  expect_identical(simd_log(c(1, NA, NaN), NA), log(c(1, NA, NaN), NA))
  expect_identical(simd_log(c(1, NA, NaN), NA_real_), rep(NA_real_, 3))
  expect_identical(simd_log(c(1, NA, NaN), NaN), c(NaN, NA, NaN))
  expect_identical(simd_log(c(1, 2), 0), log(c(1, 2), 0))
  expect_identical(suppressWarnings(simd_log(c(1, 2), -2)), c(NaN, NaN))
  expect_identical(simd_log(c(2, 0.5), 1), log(c(2, 0.5), 1))
  expect_identical(suppressWarnings(simd_log(c(2, 1), TRUE)), suppressWarnings(log(c(2, 1), TRUE)))
  # v1 takes a single base only (base R recycles a vector of bases).
  expect_error(simd_log(1, c(2, 3)), "'base' must be a single number", fixed = TRUE)
  expect_error(simd_log(1, double()), "'base' must be a single number", fixed = TRUE)
  expect_error(simd_log(1, "e"), "non-numeric argument to mathematical function", fixed = TRUE)
  expect_error(simd_log(1, NULL), "non-numeric argument to mathematical function", fixed = TRUE)
})
