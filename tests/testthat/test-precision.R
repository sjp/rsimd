# Precision modes: simd_precision() and the fast, pairwise and compensated
# accumulation schemes, checked through the self-test kernels.

test_that("simd_precision() gets, sets and validates the mode", {
  old <- options(rsimd.precision = "fast")
  on.exit(options(old))
  expect_identical(simd_precision(), "fast")
  expect_invisible(prev <- simd_precision("compensated"))
  expect_identical(prev, "fast")
  expect_identical(simd_precision(), "compensated")
  expect_identical(getOption("rsimd.precision"), "compensated")
  expect_identical(simd_precision("pairwise"), "compensated")
  # Setting the option directly is honoured.
  options(rsimd.precision = "fast")
  expect_identical(simd_precision(), "fast")
  expect_identical(.precision_code(), 0L)
  options(rsimd.precision = "compensated")
  expect_identical(.precision_code(), 2L)

  for (bad in list("bogus", NA_character_, c("fast", "pairwise"), 1, NULL)) {
    expect_error(simd_precision(bad), "precision mode must be one of", info = deparse(bad))
  }
  expect_identical(simd_precision(), "compensated")
  options(rsimd.precision = "bogus")
  expect_error(simd_precision(), "invalid option rsimd.precision")
})

test_that("cancellation: fast loses the 1, compensated keeps it, on every tier", {
  x <- tier_proof_cancel()
  res <- for_each_tier(function(tier) {
    vapply(prec_codes, function(p) fold(x, precision = p), numeric(1))
  })
  for (tier in names(res)) {
    expect_identical(res[[tier]], c(fast = 0, pairwise = 0, compensated = 1), info = tier)
  }
  # With three elements the fast result depends on the lane width.
  res <- for_each_tier(function(tier) {
    vapply(prec_codes, function(p) fold(c(1e16, 1, -1e16), precision = p), numeric(1))
  })
  for (tier in names(res)) {
    expect_true(res[[tier]][["fast"]] %in% c(0, 1), info = tier)
    expect_true(res[[tier]][["pairwise"]] %in% c(0, 1), info = tier)
    expect_identical(res[[tier]][["compensated"]], 1, info = tier)
  }
  # Base R gets 1 from its long double accumulator.
  if (has_wide_long_double()) expect_identical(sum(c(1e16, 1, -1e16)), 1)
})

test_that("compensated mode keeps Inf and NaN results", {
  for_each_tier(function(tier) {
    p <- prec_codes[["compensated"]]
    expect_identical(fold(c(1, Inf, 2), precision = p), Inf, info = tier)
    expect_identical(fold(c(-Inf, 1e308, 1e308), precision = p), -Inf, info = tier)
    expect_identical(fold(c(1e308, 1e308), precision = p), Inf, info = tier)
    expect_true(is.nan(fold(c(Inf, -Inf, 1), precision = p)), info = tier)
  })
})

test_that("modes agree with an accurate sum within their bounds", {
  set.seed(42)
  x <- stats::rnorm(1e6, mean = 1)
  ref <- sum(x)
  res <- for_each_tier(function(tier) {
    vapply(prec_codes, function(p) fold(x, precision = p), numeric(1))
  })
  fast_tol <- length(x) * 2^-52 * sum(abs(x))
  for (tier in names(res)) {
    r <- res[[tier]]
    expect_lte(abs(r[["fast"]] - ref), fast_tol)
    # Tiers agree: pairwise to 2^-50 relative, compensated to 1 ULP.
    expect_lte(abs(r[["pairwise"]] - res$none[["pairwise"]]), 2^-50 * abs(ref))
    expect_lte(abs(r[["compensated"]] - res$none[["compensated"]]), ulp(ref))
    if (has_wide_long_double()) {
      expect_lte(abs(r[["compensated"]] - ref), ulp(ref))
      expect_lte(abs(r[["pairwise"]] - ref), 2^-50 * abs(ref))
    }
  }
  # The other terms follow the same scheme.
  y <- stats::runif(1e5)
  for (term in c("sq", "abs", "xy")) {
    want <- switch(term,
      sq = sum(x[1:1e5]^2),
      abs = sum(abs(x[1:1e5])),
      xy = sum(x[1:1e5] * y)
    )
    for_each_tier(function(tier) {
      for (p in prec_codes) {
        v <- fold(x[1:1e5], if (term == "xy") y, term, precision = p)
        expect_equal(v, want, tolerance = 1e-12, info = paste(term, p, tier))
      }
    })
  }
})

test_that("pairwise does not depend on how the input is chunked", {
  # A compact sequence is read in 4096-element regions, its copy in one
  # 2^20-element chunk; x^2 sums need rounding at this size.
  n <- 2^20 + 7
  alt <- as.double(seq_len(n))
  # Built from a separate sequence: arithmetic on `alt` would expand it.
  mat <- as.double(seq_len(n)) + 0
  expect_identical(.debug_regions(alt)$path, "regions")
  expect_identical(.debug_regions(mat)$path, "contiguous")
  for_each_tier(function(tier) {
    for (term in c("x", "sq")) {
      a <- vapply(prec_codes, function(p) fold(alt, NULL, term, precision = p), numeric(1))
      m <- vapply(prec_codes, function(p) fold(mat, NULL, term, precision = p), numeric(1))
      info <- paste(term, tier)
      expect_identical(a[["pairwise"]], m[["pairwise"]], info = info)
      expect_lte(abs(a[["compensated"]] - m[["compensated"]]), ulp(m[["compensated"]]))
      expect_equal(a[["fast"]], m[["fast"]], tolerance = n * 2^-52)
    }
  })
})

test_that("pairwise is identical under a small chunk stride", {
  skip_if_no_subprocess()
  skip_on_cran()
  child <- function() {
    set.seed(7)
    x <- stats::runif(2^20 + 7, -1, 2)
    codes <- c(fast = 0L, pairwise = 1L, compensated = 2L)
    tiers <- rsimd::simd_available()
    stats::setNames(lapply(tiers, function(tier) {
      rsimd::simd_with_impl(tier, vapply(codes, function(p) {
        rsimd:::.debug_fold(x, NULL, "x", FALSE, TRUE, p)$value
      }, numeric(1)))
    }), tiers)
  }
  env <- function(stride) c(callr::rcmd_safe_env(), RSIMD_DEBUG_STRIDE = stride)
  small <- callr::r(child, env = env("4096"))
  odd <- callr::r(child, env = env("1000")) # rounded up to 1024
  default <- child()
  # (A child process may see fewer tiers, e.g. under an emulator.)
  for (tier in intersect(names(default), names(small))) {
    expect_identical(small[[tier]][["pairwise"]], default[[tier]][["pairwise"]], info = tier)
    expect_identical(odd[[tier]][["pairwise"]], default[[tier]][["pairwise"]], info = tier)
    expect_lte(
      abs(small[[tier]][["compensated"]] - default[[tier]][["compensated"]]),
      ulp(default[[tier]][["compensated"]])
    )
    expect_equal(small[[tier]][["fast"]], default[[tier]][["fast"]], tolerance = 2^-30)
  }
})
