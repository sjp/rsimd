# Self-tests of the test helpers (helper-tiers.R, helper-expect.R,
# helper-cases.R, helper-extended.R).

test_that("tier helpers cover every available tier and end with none", {
  tiers <- tiers_to_test()
  expect_identical(tiers, simd_available())
  expect_identical(tiers[length(tiers)], "none")
  res <- with_each_tier(function() simd_current())
  expect_identical(names(res), tiers)
  expect_identical(unname(vapply(res, as.vector, "")), tiers)
  expect_no_error(skip_if_no_tier("none"))
  expect_condition(skip_if_no_tier("rvv"), "tier rvv unavailable", class = "skip")
})

test_that("expect_simd_equal passes for agreeing tiers and fails for a wrong oracle", {
  expect_success(expect_simd_equal(function(x) sum(x), 1:10))
  skip_if(length(simd_available()) < 2, "only the none tier is available")
  # A result that depends on the tier, as a broken kernel would.
  by_tier <- function(x) if (simd_current() == "none") 1 else 2
  by_tier_dbl <- function(x) if (simd_current() == "none") 1 else 1 + 1e-6
  expect_failure(expect_simd_equal(by_tier, 1), "differs from none")
  expect_failure(expect_simd_equal(by_tier_dbl, c(1, 2), precision = "compensated"), "tier")
  # Missing values must agree exactly, whatever the tolerance.
  na_or_nan <- function(x) if (simd_current() == "none") NA_real_ else NaN
  expect_failure(expect_simd_equal(na_or_nan, 1, tolerance = 1), "is.nan")
  inf_or_big <- function(x) if (simd_current() == "none") Inf else 1e308
  expect_failure(expect_simd_equal(inf_or_big, 1, precision = "fast"), "infinite")
})

test_that("the double bounds scale with the precision mode", {
  x <- c(1e16, 1, -1e16)
  # One unit of difference on a cancelling sum: within the fast and
  # pairwise bounds (they scale with sum(abs(x))), not the compensated one.
  expect_null(double_mismatch(1, 0, mode_bound(x, "fast"), 2^-50))
  expect_null(double_mismatch(1, 0, mode_bound(x, "pairwise"), 2^-50))
  expect_match(double_mismatch(1, 0, mode_bound(x, "compensated"), 2^-50), "difference")
  expect_null(double_mismatch(1 + 2^-52, 1, 0, 2^-50))
  expect_match(double_mismatch(1L, 1, 0, 0), "types differ")
  expect_match(double_mismatch(c(NA, NaN), c(NaN, NA), 0, 0), "is.nan")
  expect_null(double_mismatch(c(NA, NaN), c(NaN, NA), 0, 0, nan = FALSE))
})

test_that("expect_simd_matches_base compares every tier with base R", {
  expect_success(expect_simd_matches_base(function(x) sum(x), sum, c(1, 2, 3)))
  expect_failure(expect_simd_matches_base(function(x) sum(x) + 1, sum, c(1, 2, 3)), "base R")
  # NA and NaN together: base R's pick depends on the order, so only
  # is.na() is compared.
  expect_success(expect_simd_matches_base(function(x) NA_real_, sum, c(NaN, NA)))
})

test_that("edge_lengths contains the width boundaries of every tier", {
  lens <- edge_lengths()
  for (W in lane_widths) {
    expect_true(all(c(W - 1, W, W + 1, 2 * W + 3, 16 * W + 5) %in% lens), info = W)
  }
  expect_true(all(c(0, 1, 4095, 4096, 4097, 1e6) %in% lens))
  expect_false(is.unsorted(lens))
  expect_identical(lens, unique(lens))
  expect_identical(edge_lengths(4), c(0, 1, 3, 4, 5, 11, 69, 4095, 4096, 4097, 1e6))
})

test_that("edge values contain the documented specials", {
  d <- edge_doubles()
  expect_true(all(c(0, 1, -1, Inf, -Inf, pi, 5e-324, .Machine$double.xmax) %in% d))
  expect_identical(sum(is.na(d) & !is.nan(d)), 1L)
  expect_identical(sum(is.nan(d)), 1L)
  expect_true(any(1 / d[!is.na(d) & d == 0] < 0)) # -0
  expect_false(.Machine$double.xmax %in% edge_doubles(xmax = FALSE))
  expect_identical(length(edge_doubles(xmax = FALSE)), length(d) - 1L)
  expect_identical(edge_ints(), c(0L, 1L, -1L, .Machine$integer.max, -.Machine$integer.max, NA))
  expect_identical(edge_lgl(), c(TRUE, FALSE, NA))
})

test_that("altrep_inputs are compact sequences read in regions", {
  for (n in c(2, 10, 5000)) {
    inputs <- altrep_inputs(n)
    expect_named(inputs, c("colon", "seq_len", "double"))
    for (nm in names(inputs)) {
      expect_true(takes_region_path(inputs[[nm]]), info = paste(nm, n))
      expect_identical(length(inputs[[nm]]), as.integer(n))
    }
  }
  expect_false(takes_region_path(c(1L, 2L)))
  expect_error(altrep_inputs(1))
})

test_that("rand_vec is reproducible, typed and leaves the RNG state alone", {
  set.seed(42)
  state <- .Random.seed
  a <- rand_vec("double", 1000, seed = 3)
  expect_identical(.Random.seed, state)
  expect_identical(a, rand_vec("double", 1000, seed = 3))
  expect_false(identical(a, rand_vec("double", 1000, seed = 4)))
  expect_type(a, "double")
  expect_true(any(is.na(a) & !is.nan(a)) && any(is.nan(a)) && any(is.infinite(a)))
  i <- rand_vec("integer", 1000, seed = 3)
  expect_type(i, "integer")
  expect_true(anyNA(i) && all(abs(i[!is.na(i)]) <= 1e6))
  l <- rand_vec("logical", 1000, seed = 3)
  expect_type(l, "logical")
  expect_true(anyNA(l))
  expect_false(anyNA(rand_vec("double", 100, na_frac = 0, nan_frac = 0)))
  expect_error(rand_vec("complex", 10), "unknown type")
})

test_that("skip_unless_extended follows RSIMD_EXTENDED_TESTS", {
  old <- Sys.getenv("RSIMD_EXTENDED_TESTS", unset = NA)
  on.exit(if (is.na(old)) {
    Sys.unsetenv("RSIMD_EXTENDED_TESTS")
  } else {
    Sys.setenv(RSIMD_EXTENDED_TESTS = old)
  })
  Sys.setenv(RSIMD_EXTENDED_TESTS = "false")
  expect_condition(skip_unless_extended(), "RSIMD_EXTENDED_TESTS", class = "skip")
  Sys.setenv(RSIMD_EXTENDED_TESTS = "true")
  expect_no_condition(skip_unless_extended())
})
