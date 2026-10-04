# simd_sigmoid(), simd_softmax() and simd_log_softmax(): accuracy against
# base R, the stable forms' tails, missing and infinite values, lengths,
# chunked and ALTREP inputs, and the precision modes.

base_softmax <- function(x) exp(x - max(x)) / sum(exp(x - max(x)))
base_log_softmax <- function(x) x - max(x) - log(sum(exp(x - max(x))))

prec_modes <- c("fast", "pairwise", "compensated")

# The largest relative difference of a from b, abs(a - b) / pmax(abs(b),
# floor), over their finite elements; Inf when their infinities differ.
max_rel <- function(a, b, floor = 0) {
  if (!identical(a[!is.finite(b)], b[!is.finite(b)])) {
    return(Inf)
  }
  ok <- is.finite(b)
  max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), floor), 0)
}

# The relative error allowed against base R's formulas for n elements in
# precision mode `mode`: a few units in the last place for exp, the
# division and the subtraction, plus the relative error of the sum of the
# exponentials (positive terms, so the bounds of helper-expect.R are
# relative ones).
softmax_bound <- function(n, mode = simd_precision()) {
  n <- max(n, 1)
  8 * eps + switch(mode,
    fast = n * 8 * eps,
    pairwise = (128 + log2(n)) * eps,
    compensated = 4 * eps
  )
}

softmax_inputs <- function() {
  list(
    small = c(1, 2, 3),
    spread = math_random(1000, -50, 50, seed = 1L),
    offset = 1e4 + math_random(777, -5, 5, seed = 2L),
    tiny = math_random(129, -1e-3, 1e-3, seed = 3L),
    ties = rep(c(0.5, -2, 0.5), 50),
    big = c(1000, 1001),
    huge = c(-1e308, 1e308, 0),
    neg_inf = c(-Inf, 1, -Inf, 3)
  )
}

test_that("simd_sigmoid() matches plogis() above -708 and exp(x) in the lower tail", {
  x <- c(seq(-708, 800, by = 0.01), math_random(2000, -708, 800, seed = 15L), 1e-300, -1e-300)
  for_each_tier(function(tier) {
    expect_close_to(simd_sigmoid(x), plogis(x), 4, paste("sigmoid vs plogis on", tier), x = x)
  })
  expect_tiers_close(simd_sigmoid, x, ulps = 2)
  # Below -709.78 plogis() underflows to 0 (exp(-x) overflows); the stable
  # form keeps the subnormal value, which is exp(x) itself as 1 + exp(x)
  # rounds to 1.
  tail <- c(seq(-746, -708, by = 0.01), -709.78, -720, -744.4)
  for_each_tier(function(tier) {
    expect_close_to(simd_sigmoid(tail), exp(tail), 2, paste("sigmoid tail on", tier), x = tail)
  })
  expect_identical(plogis(-720), 0)
  expect_gt(simd_sigmoid(-720), 0)
})

test_that("simd_sigmoid() has exact limits and keeps missing values silently", {
  x <- c(0, -0, Inf, -Inf, -800, 800, -746, 40, NA, NaN)
  want <- c(0.5, 0.5, 1, 0, 0, 1, 0, 1, NA, NaN)
  expect_tiers_give(want, simd_sigmoid, x)
  expect_tier_warnings(character(), simd_sigmoid, x)
  expect_identical(plogis(x[1:8]), want[1:8])
  ints <- c(-3L, 0L, 1L, NA, .Machine$integer.max, -.Machine$integer.max)
  for_each_tier(function(tier) {
    expect_identical(simd_sigmoid(ints), simd_sigmoid(as.double(ints)), info = tier)
    expect_identical(simd_sigmoid(c(TRUE, FALSE, NA)), simd_sigmoid(c(1, 0, NA)), info = tier)
  })
})

test_that("simd_sigmoid() is monotonic and symmetric about 0", {
  grid <- seq(-800, 800, by = 0.001)
  small <- seq(-1, 1, by = 1 / 1024)
  for_each_tier(function(tier) {
    expect_false(is.unsorted(simd_sigmoid(grid)), info = tier)
    # 1 - s(x) cancels, so the symmetry is only checked for |x| <= 1.
    expect_close_to(1 - simd_sigmoid(small), simd_sigmoid(-small), 4,
      paste("1 - s(x) vs s(-x) on", tier),
      x = small
    )
  })
})

test_that("simd_softmax() and simd_log_softmax() match base R's formulas in every precision mode", {
  inputs <- softmax_inputs()
  for (mode in prec_modes) {
    old <- simd_precision(mode)
    for (name in names(inputs)) {
      x <- inputs[[name]]
      info <- paste(name, mode)
      for_each_tier(function(tier) {
        p <- simd_softmax(x)
        want <- base_softmax(x)
        bound <- softmax_bound(length(x))
        expect_lte(max_rel(p[want > 0], want[want > 0]), bound)
        expect_identical(p == 0, want == 0, info = paste(info, tier))
        lp <- simd_log_softmax(x)
        expect_lte(max_rel(lp, base_log_softmax(x), 1), bound)
      })
    }
    simd_precision(old)
  }
  expect_identical(base_softmax(c(1, 2, 3)), c(
    0.090030573170380462, 0.244728471054797642, 0.665240955774821785
  ))
})

test_that("simd_log_softmax() is stable where log(softmax) is not", {
  x <- c(1000, 1001)
  want <- c(-1.31326168751822281, -0.31326168751822286)
  for_each_tier(function(tier) {
    expect_close_to(simd_log_softmax(x), want, 1, paste("log_softmax on", tier))
    expect_close_to(simd_softmax(x), exp(want), 2, paste("softmax on", tier))
  })
  # Spread beyond exp's range: the naive log of the softmax is -Inf.
  y <- c(0, 800, 1600)
  expect_identical(log(simd_softmax(y))[1], -Inf)
  expect_tiers_give(c(-1600, -800, 0), simd_log_softmax, y)
  # Moderate inputs: the log of the softmax, and the exponentials sum to 1.
  z <- math_random(5000, -20, 20, seed = 4L)
  for_each_tier(function(tier) {
    lp <- simd_log_softmax(z)
    expect_lte(max_rel(lp, log(simd_softmax(z)), 1), 1e-12)
    expect_lte(abs(sum(exp(lp)) - 1), 1e-12)
  })
})

test_that("missing and infinite values give base R's NA and NaN patterns without warnings", {
  cases <- list(
    c(1, Inf, 2), c(1, Inf, Inf), c(-Inf, -Inf), c(-Inf, 1), Inf, -Inf, c(Inf, -Inf),
    c(1, NA, 2), c(1, NaN, 2), NA_real_, NaN, c(-Inf, NaN), c(Inf, NaN), c(Inf, NA)
  )
  for (x in cases) {
    for (pair in list(list(simd_softmax, base_softmax), list(simd_log_softmax, base_log_softmax))) {
      want <- suppressWarnings(pair[[2]](x))
      res <- expect_tier_warnings(character(), pair[[1]], x)
      expect_identical(is.na(res), is.na(want), info = deparse(x))
      expect_identical(is.nan(res), is.nan(want), info = deparse(x))
      expect_identical(res[!is.na(res)], want[!is.na(want)], info = deparse(x))
    }
  }
  # NA wins over NaN wherever they are (base R depends on their order).
  for (x in list(c(NaN, 1, NA), c(NA, NaN), c(Inf, NA, NaN, -Inf), c(1:3, NA))) {
    expect_tiers_give(rep(NA_real_, length(x)), simd_softmax, x)
    expect_tiers_give(rep(NA_real_, length(x)), simd_log_softmax, x)
  }
  expect_tiers_give(rep(NaN, 3), simd_softmax, c(1, NaN, Inf))
  expect_tiers_give(c(0, 1, 0), simd_softmax, c(-Inf, 5, -Inf))
  expect_tiers_give(c(-Inf, 0, -Inf), simd_log_softmax, c(-Inf, 5, -Inf))
})

test_that("results have every length, chunked and ALTREP inputs included, and are bare", {
  for (n in c(edge_lengths(), if (!reduced_lengths()) 2^20 + 7)) {
    x <- math_random(n, -30, 30, seed = n %% 1000)
    expect_length(simd_softmax(x), n)
    expect_length(simd_log_softmax(x), n)
    if (n > 0) {
      for_each_tier(function(tier) {
        bound <- softmax_bound(n)
        expect_lte(max_rel(simd_softmax(x), base_softmax(x)), bound)
        expect_lte(max_rel(simd_log_softmax(x), base_log_softmax(x), 1), bound)
      })
    }
  }
  expect_identical(simd_softmax(double()), double())
  expect_identical(simd_log_softmax(integer()), double())
  expect_tiers_give(1, simd_softmax, 42)
  expect_tiers_give(0, simd_log_softmax, -7L)
  for (x in altrep_inputs(5000)) {
    expect_true(takes_region_path(x))
    for_each_tier(function(tier) {
      expect_identical(simd_softmax(x), simd_softmax(as.double(x) + 0), info = tier)
      expect_identical(simd_log_softmax(x), simd_log_softmax(as.double(x) + 0), info = tier)
    })
  }
  ints <- c(3L, -2L, 0L, 7L, -.Machine$integer.max)
  for_each_tier(function(tier) {
    expect_identical(simd_softmax(ints), simd_softmax(as.double(ints)), info = tier)
    expect_identical(simd_softmax(c(TRUE, FALSE)), simd_softmax(c(1, 0)), info = tier)
    expect_identical(simd_log_softmax(ints), simd_log_softmax(as.double(ints)), info = tier)
  })
  # Matrices and named vectors are normalised as a whole and lose their
  # attributes.
  m <- matrix(c(1, 2, 3, 4), 2)
  expect_null(attributes(simd_softmax(m)))
  expect_identical(simd_softmax(m), simd_softmax(c(1, 2, 3, 4)))
  expect_null(attributes(simd_log_softmax(c(a = 1, b = 2))))
  expect_null(attributes(simd_sigmoid(m)))
})

test_that("the softmax sums to 1 within the precision mode's bound", {
  for (n in c(1e3, 1e5, 1e6)) {
    x <- math_random(n, -10, 10, seed = n %% 977)
    for (mode in prec_modes) {
      old <- simd_precision(mode)
      for_each_tier(function(tier) {
        err <- abs(sum(simd_softmax(x)) - 1)
        bound <- if (mode == "fast") 1e-10 else n * .Machine$double.eps
        expect_lte(err, bound)
      })
      simd_precision(old)
    }
  }
})

test_that("pairwise softmax does not depend on the chunk stride", {
  skip_if_no_subprocess()
  skip_on_cran()
  child <- function() {
    set.seed(11)
    x <- stats::runif(2^20 + 7, -20, 20)
    old <- rsimd::simd_precision("pairwise")
    on.exit(rsimd::simd_precision(old))
    tiers <- rsimd::simd_available()
    stats::setNames(lapply(tiers, function(tier) {
      rsimd::simd_with_impl(tier, list(rsimd::simd_softmax(x), rsimd::simd_log_softmax(x)))
    }), tiers)
  }
  small <- callr::r(child, env = c(callr::rcmd_safe_env(), RSIMD_DEBUG_STRIDE = "4096"))
  default <- child()
  for (tier in intersect(names(default), names(small))) {
    expect_identical(small[[tier]], default[[tier]], info = tier)
  }
})

test_that("non-numeric input is rejected", {
  for (f in list(simd_sigmoid, simd_softmax, simd_log_softmax)) {
    expect_error(f("a"), "must be an atomic vector")
    expect_error(f(as.raw(1:3)), "non-numeric argument to mathematical function")
    expect_error(f(list(1)), "must be an atomic vector")
  }
  expect_error(simd_softmax(1i), "simd_softmax() does not support 'x' of type complex", fixed = TRUE)
  expect_error(simd_log_softmax(1i), "simd_log_softmax() does not support 'x' of type complex",
    fixed = TRUE
  )
  expect_error(simd_sigmoid(1i), "simd_sigmoid() does not support 'x' of type complex", fixed = TRUE)
  skip_if_not_installed("bit64")
  expect_error(simd_softmax(bit64::as.integer64(1)),
    "simd_softmax() does not support 'x' of type integer64",
    fixed = TRUE
  )
})
