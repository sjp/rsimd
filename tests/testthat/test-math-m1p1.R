# simd_exp2m1, simd_exp10m1, simd_log2p1 and simd_log10p1 (C23's exp2m1,
# exp10m1, log2p1 and log10p1), and simd_sincospi. Their accuracy against
# a long double reference is measured by tools/check_vector_layer.sh; here
# every tier is compared with the none tier and with base R where base R's
# expression is accurate, and the special values are checked exactly.

m1p1 <- list(
  exp2m1 = list(f = simd_exp2m1, base = function(x) 2^x - 1, inv = function(x) log2(1 + x)),
  exp10m1 = list(f = simd_exp10m1, base = function(x) 10^x - 1, inv = function(x) log10(1 + x)),
  log2p1 = list(f = simd_log2p1, base = function(x) log2(1 + x)),
  log10p1 = list(f = simd_log10p1, base = function(x) log10(1 + x))
)

# Arguments over the domain of each function: moderate, tiny, large and
# special values, and whole numbers.
m1p1_inputs <- function(name, n = 3000, seed = 20261006L) {
  with_seed(seed, {
    tiny <- 10^stats::runif(200, -320, -1) * sample(c(-1, 1), 200, replace = TRUE)
    whole <- as.double(-60:60)
    special <- c(0, -0, Inf, -Inf, NA, NaN, 5e-324, -5e-324, 1e-300)
    x <- switch(name,
      exp2m1 = c(stats::runif(n, -60, 60), stats::runif(200, 54, 1030), -1100.5, 53.5, 1024),
      exp10m1 = c(stats::runif(n, -18, 18), stats::runif(200, 17, 310), -330.5, 308.25, 309),
      c(
        10^stats::runif(n, -15, 300), -10^stats::runif(n, -15, -1e-9), -1, -2,
        2^(1:53) - 1, 10^(1:15) - 1
      )
    )
    c(special, tiny, whole, x)
  })
}

test_that("the m1/p1 functions agree with the none tier within 2 ULP on every tier", {
  for (name in names(m1p1)) {
    x <- m1p1_inputs(name)
    expect_tiers_close(m1p1[[name]]$f, x, ulps = 2, label = name, x = x)
    xi <- c(-60:60, NA)
    expect_tiers_close(m1p1[[name]]$f, xi, ulps = 2, label = paste(name, "integer"))
    for (n in sweep_lengths()) {
      v <- rep_len(c(0.3, -0.7, 1e-12, 25, NA, NaN, -0), n)
      expect_tiers_close(m1p1[[name]]$f, v, ulps = 2, label = paste(name, n))
    }
  }
})

test_that("the m1/p1 functions match base R where its expression is accurate", {
  with_seed(1L, {
    big <- c(stats::runif(2000, 1, 60), stats::runif(2000, -60, -1))
    pos <- c(stats::runif(2000, 1, 1e6), 10^stats::runif(2000, 0, 300))
    for (tier in tiers_to_test()) {
      simd_with_impl(tier, {
        # b^x - 1 loses nothing for |x| >= 1, nor log(1 + x) for x >= 1.
        expect_close_to(simd_exp2m1(big), 2^big - 1, 3, paste(tier, "exp2m1"))
        expect_close_to(simd_exp10m1(big / 4), 10^(big / 4) - 1, 3, paste(tier, "exp10m1"))
        expect_close_to(simd_log2p1(pos), log2(1 + pos), 3, paste(tier, "log2p1"))
        expect_close_to(simd_log10p1(pos), log10(1 + pos), 3, paste(tier, "log10p1"))
        # Near 0 the first terms of the series: x log(b) (1 + x log(b) / 2).
        small <- 10^stats::runif(500, -300, -9) * sample(c(-1, 1), 500, replace = TRUE)
        for (b in c(2, 10)) {
          t <- small * log(b)
          f <- if (b == 2) simd_exp2m1 else simd_exp10m1
          g <- if (b == 2) simd_log2p1 else simd_log10p1
          expect_close_to(f(small), t + t * t / 2, 3, paste(tier, "expm1 small", b),
            abs = 5e-324, abs_below = 1e-300
          )
          expect_close_to(g(small), (small - small * small / 2) / log(b), 3,
            paste(tier, "log1p small", b),
            abs = 5e-324, abs_below = 1e-300
          )
        }
      })
    }
  })
})

test_that("the m1/p1 functions are exact where b^x - 1 or log_b(1 + x) is", {
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_identical(simd_exp2m1(-53:53), 2^(-53:53) - 1, info = tier)
      expect_identical(simd_exp10m1(1:15), 10^(1:15) - 1, info = tier)
      expect_identical(simd_log2p1(2^(-52:53) - 1), as.double(-52:53), info = tier)
      expect_identical(simd_log10p1(10^(1:15) - 1), as.double(1:15), info = tier)
      for (f in list(simd_exp2m1, simd_exp10m1, simd_log2p1, simd_log10p1)) {
        expect_identical(1 / f(c(0, -0)), c(Inf, -Inf), info = tier)
        expect_identical(f(c(NA, NaN)), c(NA, NaN), info = tier)
        expect_identical(f(NA_integer_), NA_real_, info = tier)
        expect_identical(f(numeric(0)), numeric(0), info = tier)
        expect_identical(f(Inf), Inf, info = tier)
      }
      expect_identical(simd_exp2m1(c(-Inf, -1100, 1024)), c(-1, -1, Inf), info = tier)
      expect_identical(simd_exp10m1(c(-Inf, -400, 309)), c(-1, -1, Inf), info = tier)
      expect_identical(simd_log2p1(-1), -Inf, info = tier)
      expect_identical(simd_log10p1(-1), -Inf, info = tier)
      expect_identical(simd_log2p1(c(TRUE, FALSE)), c(1, 0), info = tier)
    })
  }
})

test_that("the m1/p1 functions warn as base R does, and only take real input", {
  expect_tier_warnings("NaNs produced", simd_log2p1, c(-2, 1, NA))
  expect_tier_warnings("NaNs produced", simd_log10p1, -Inf)
  expect_tier_warnings(character(0), simd_log2p1, c(-1, NaN, NA))
  expect_tier_warnings(character(0), simd_exp2m1, c(-Inf, Inf, 1e300, NaN))
  expect_identical(suppressWarnings(simd_log2p1(c(-2, -Inf))), c(NaN, NaN))
  for (f in c("simd_exp2m1", "simd_exp10m1", "simd_log2p1", "simd_log10p1")) {
    expect_error(get(f)(1i), paste0(f, "\\(\\) does not support 'x' of type complex"))
    expect_error(get(f)("a"), "not character")
  }
  if (has_bit64()) {
    expect_error(simd_exp2m1(simd_as_integer64(1)), "of type integer64")
  }
})

test_that("the m1/p1 functions do not depend on rsimd.math_accuracy", {
  for (name in names(m1p1)) {
    x <- m1p1_inputs(name, n = 500)
    for (tier in tiers_to_test()) {
      simd_with_impl(tier, {
        old <- simd_math_accuracy("fast")
        fast <- suppressWarnings(m1p1[[name]]$f(x))
        simd_math_accuracy(old)
        expect_identical(fast, suppressWarnings(m1p1[[name]]$f(x)), info = paste(tier, name))
      })
    }
  }
})

test_that("simd_sincospi is list(sin = simd_sinpi(x), cos = simd_cospi(x))", {
  with_seed(2L, {
    x <- c(
      stats::runif(2000, -4, 4), (-40:40) / 4, 1e300, -1e17, 2.5e8 + 0.5, 1e-310, -5e-324,
      0, -0, Inf, -Inf, NA, NaN
    )
  })
  for (mode in c("accurate", "fast")) {
    old <- simd_math_accuracy(mode)
    for (tier in tiers_to_test()) {
      simd_with_impl(tier, {
        got <- suppressWarnings(simd_sincospi(x))
        want <- suppressWarnings(list(sin = simd_sinpi(x), cos = simd_cospi(x)))
        expect_identical(got, want, info = paste(tier, mode))
        expect_identical(simd_sincospi(-3:3), list(sin = simd_sinpi(-3:3), cos = simd_cospi(-3:3)),
          info = tier
        )
      })
    }
    simd_math_accuracy(old)
  }
  if (base_has_platform_sinpi()) {
    expect_close_to(suppressWarnings(simd_sincospi(x)$sin), suppressWarnings(sinpi(x)), 2,
      zeros = FALSE
    )
  }
  expect_identical(simd_sincospi(numeric(0)), list(sin = numeric(0), cos = numeric(0)))
  w <- warnings_of(simd_sincospi(c(Inf, 1, -Inf)))
  expect_identical(w, "NaNs produced")
  expect_error(simd_sincospi(1i), "simd_sincospi\\(\\) does not support 'x' of type complex")
})

test_that("the new math functions keep a simd_vec's pin", {
  for (tier in tiers_to_test()) {
    v <- simd_vec(c(0.5, 3), impl = tier)
    for (f in list(simd_exp2m1, simd_exp10m1, simd_log2p1, simd_log10p1)) {
      r <- f(v)
      expect_true(is_simd_vec(r))
      expect_identical(simd_impl(r), tier)
    }
    s <- simd_sincospi(v)
    expect_true(is_simd_vec(s$sin) && is_simd_vec(s$cos))
    expect_identical(simd_impl(s$cos), tier)
  }
})
