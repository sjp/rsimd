# simd_sum(): every tier against the none oracle and against base R.

modes <- c("fast", "pairwise", "compensated")

# Missing values follow the rsimd rule: NA if any NA, else NaN if any NaN.
expected_missing <- function(x) {
  if (any(is.na(x) & !is.nan(x))) NA else if (any(is.nan(x))) NaN else 0
}

test_that("simd_sum matches the oracle on every tier for all edge lengths", {
  for (n in edge_lengths()) {
    seed <- as.integer(n %% 1000) + 1L
    inputs <- list(
      double = rand_vec("double", n, seed = seed),
      double_clean = rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = seed),
      integer = rand_vec("integer", n, seed = seed),
      integer_clean = rand_vec("integer", n, na_frac = 0, seed = seed),
      logical = rand_vec("logical", n, seed = seed)
    )
    for (nm in names(inputs)) {
      x <- inputs[[nm]]
      for (na_rm in c(FALSE, TRUE)) {
        # Precision only affects doubles.
        for (mode in if (is.double(x)) modes else "fast") {
          expect_simd_equal(simd_sum, x, na.rm = na_rm, precision = mode)
        }
      }
      if (!anyNA(x)) expect_simd_equal(simd_sum, x, na_check = FALSE)
    }
  }
})

test_that("simd_sum agrees with base R across edge lengths", {
  for (n in edge_lengths()) {
    seed <- as.integer(n %% 1000) + 1L
    # One kind of missing value at a time, where base R is deterministic.
    xs <- list(
      rand_vec("double", n, nan_frac = 0, seed = seed),
      rand_vec("double", n, na_frac = 0, seed = seed),
      rand_vec("integer", n, seed = seed),
      rand_vec("logical", n, seed = seed)
    )
    for (x in xs) {
      for (na_rm in c(FALSE, TRUE)) {
        for (mode in if (is.double(x)) modes else "fast") {
          expect_simd_matches_base(simd_sum, sum, x, na.rm = na_rm, precision = mode)
        }
      }
    }
  }
})

test_that("edge doubles in every order agree with the oracle and the NA rules", {
  d <- edge_doubles(xmax = FALSE)
  check <- function(x) {
    for (mode in modes) {
      res <- expect_simd_equal(simd_sum, x, precision = mode)
      if (anyNA(x)) {
        expect_identical(is.nan(res$none), is.nan(expected_missing(x)), info = deparse(x))
      } else {
        # Without missing values only Inf + -Inf gives NaN.
        expect_identical(is.nan(res$none), Inf %in% x && -Inf %in% x, info = deparse(x))
      }
      expect_simd_equal(simd_sum, x, na.rm = TRUE, precision = mode)
    }
  }
  for (i in seq_along(d)) {
    for (j in seq_along(d)) check(d[c(i, j)])
  }
  for (k in 1:10) {
    p <- with_seed(k, sample(d))
    check(p)
    check(rev(p))
    # Spread over a longer vector so the values land in different lanes.
    long <- numeric(16 * 32 + 5)
    long[with_seed(k, sort(sample.int(length(long), length(p))))] <- p
    check(long)
  }
})

test_that("NA and NaN in either order give NA on every tier", {
  for (x in list(c(NA, NaN), c(NaN, NA), c(1, NaN, 2, NA), c(NA, 1:40, NaN))) {
    for (mode in modes) {
      res <- with_each_tier(function() {
        old <- simd_precision(mode)
        on.exit(simd_precision(old))
        simd_sum(x)
      })
      for (tier in names(res)) {
        expect_true(is.na(res[[tier]]) && !is.nan(res[[tier]]), info = paste(tier, mode))
      }
    }
  }
  res <- with_each_tier(function() simd_sum(c(1, NaN, Inf)))
  for (tier in names(res)) expect_true(is.nan(res[[tier]]), info = tier)
  res <- with_each_tier(function() simd_sum(c(1, NA, NaN), na.rm = TRUE))
  for (tier in names(res)) expect_identical(res[[tier]], 1, info = tier)
})

test_that("double.xmax: order-independent cases agree on every tier", {
  xmax <- .Machine$double.xmax
  cases <- list(
    list(x = xmax, want = xmax),
    list(x = c(xmax, Inf), want = Inf),
    list(x = c(-xmax, -Inf), want = -Inf),
    list(x = c(xmax, NA), want = NA_real_),
    list(x = c(xmax, NaN), want = NaN),
    list(x = c(Inf, -Inf), want = NaN),
    list(x = c(xmax, -xmax), want = 0)
  )
  for (case in cases) {
    for (mode in modes) {
      res <- expect_simd_equal(simd_sum, case$x, precision = mode)
      for (tier in names(res)) {
        expect_identical(res[[tier]], case$want, info = paste(deparse(case$x), mode, tier))
      }
    }
  }
})

test_that("integer and logical sums follow base R result types", {
  res <- with_each_tier(function() {
    list(
      simd_sum(1:10), simd_sum(c(TRUE, NA)), simd_sum(c(TRUE, FALSE, TRUE)),
      simd_sum(integer(0)), simd_sum(logical(0)), simd_sum(numeric(0)),
      simd_sum(rep(.Machine$integer.max, 3L)), simd_sum(c(.Machine$integer.max, 1L)),
      simd_sum(rep(-.Machine$integer.max, 2L)), simd_sum(c(NA, 1L), na.rm = TRUE),
      simd_sum(NA_integer_, na.rm = TRUE), simd_sum(c(-0, -0))
    )
  })
  want <- list(
    55L, NA_integer_, 2L, 0L, 0L, 0, 6442450941, 2147483648, -4294967294, 1L, 0L, 0
  )
  for (tier in names(res)) expect_identical(res[[tier]], want, info = tier)
  # The same as base R.
  for (x in list(1:10, c(TRUE, NA), integer(0), rep(.Machine$integer.max, 3L))) {
    expect_simd_matches_base(simd_sum, sum, x)
  }
  # Overflow and NA: NA wins without na.rm, and the rest still overflows
  # into a double with it.
  x <- c(rep(.Machine$integer.max, 3L), NA)
  expect_simd_matches_base(simd_sum, sum, x)
  expect_simd_matches_base(simd_sum, sum, x, na.rm = TRUE)
})

test_that("cancellation: compensated gives 1, fast depends only on the lane width", {
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      old <- simd_precision("compensated")
      expect_identical(simd_sum(c(1e16, 1, -1e16)), 1, info = tier)
      expect_identical(simd_sum(tier_proof_cancel()), 1, info = tier)
      simd_precision("fast")
      expect_true(simd_sum(c(1e16, 1, -1e16)) %in% c(0, 1), info = tier)
      expect_identical(simd_sum(tier_proof_cancel()), 0, info = tier)
      simd_precision(old)
    })
  }
})

test_that("compact sequences are summed in regions and equal the expanded vector", {
  expect_true(takes_region_path(1:1e6))
  expect_simd_matches_base(simd_sum, sum, 1:1e6)
  for (n in c(2, 4097, 1e6)) {
    for (x in altrep_inputs(n)) {
      expect_true(takes_region_path(x))
      expanded <- x + if (is.integer(x)) 0L else 0
      expect_false(takes_region_path(expanded))
      for (mode in modes) {
        res <- expect_simd_equal(simd_sum, x, precision = mode)
        res2 <- expect_simd_equal(simd_sum, expanded, precision = mode)
        expect_identical(res, res2, info = paste(typeof(x), n, mode))
      }
    }
  }
})

test_that("compensated sum of 1e6 random doubles is within 1 ULP of the exact sum", {
  skip_if_not(has_wide_long_double(), "no long double accumulator in base R")
  x <- rand_vec("double", 1e6, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = 11)
  # Exact enough: base R accumulates in long double, in order of magnitude.
  exact <- sum(x[order(abs(x))])
  res <- with_each_tier(function() {
    old <- simd_precision("compensated")
    on.exit(simd_precision(old))
    simd_sum(x)
  })
  for (tier in names(res)) {
    expect_lte(abs(res[[tier]] - exact), ulp(exact), label = tier)
  }
})

test_that("the result does not depend on the chunking for pairwise", {
  skip_on_cran()
  skip_if_not_installed("callr")
  x <- rand_vec("double", 2^20 + 7, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = 5)
  child <- function(x, modes) {
    library(rsimd)
    lapply(stats::setNames(nm = simd_available()), function(tier) {
      simd_with_impl(tier, vapply(modes, function(m) {
        old <- simd_precision(m)
        on.exit(simd_precision(old))
        simd_sum(x)
      }, numeric(1)))
    })
  }
  here <- child(x, modes)
  env <- c(callr::rcmd_safe_env(), RSIMD_DEBUG_STRIDE = "4096")
  there <- callr::r(child, list(x, modes), env = env)
  # A child of an emulated R may see fewer tiers.
  for (tier in intersect(names(here), names(there))) {
    expect_identical(there[[tier]][["pairwise"]], here[[tier]][["pairwise"]], info = tier)
    expect_lte(abs(there[[tier]][["compensated"]] - here[[tier]][["compensated"]]),
      ulp(here[[tier]][["compensated"]]),
      label = tier
    )
    expect_lte(abs(there[[tier]][["fast"]] - here[[tier]][["fast"]]),
      mode_bound(x, "fast"),
      label = tier
    )
  }
})

test_that("na_check = FALSE gives the same result on NA-free input", {
  x <- rand_vec("integer", 5000, na_frac = 0, seed = 2)
  d <- rand_vec("double", 5000, na_frac = 0, nan_frac = 0, seed = 2)
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_identical(simd_sum(x, na_check = FALSE), simd_sum(x), info = tier)
      expect_identical(simd_sum(d, na_check = FALSE), simd_sum(d), info = tier)
      old <- options(rsimd.na_check = FALSE)
      expect_identical(simd_sum(x), sum(x), info = tier)
      options(old)
      # A broken promise gives an unspecified value, but no error.
      expect_type(simd_sum(c(1L, NA), na_check = FALSE), "integer")
      expect_type(simd_sum(c(1, NA), na_check = FALSE), "double")
    })
  }
})

test_that("simd_sum validates its arguments", {
  expect_error(simd_sum(as.raw(1)), "invalid 'type' (raw) of argument", fixed = TRUE)
  expect_error(simd_sum("a"), "'x' must be an atomic vector", fixed = TRUE)
  expect_error(simd_sum(factor("a")), "not factor", fixed = TRUE)
  expect_error(simd_sum(list(1)), "not list", fixed = TRUE)
  expect_error(simd_sum(1, na.rm = NA), "'na.rm' must be TRUE or FALSE", fixed = TRUE)
  expect_error(simd_sum(1, na.rm = "yes"), "'na.rm' must be TRUE or FALSE", fixed = TRUE)
  expect_error(simd_sum(1, na_check = 1:2), "'na_check' must be TRUE or FALSE", fixed = TRUE)
})

test_that("attributes of x do not matter", {
  m <- matrix(as.double(1:6), 2, dimnames = list(c("a", "b"), NULL))
  expect_identical(simd_sum(m), 21)
  expect_identical(simd_sum(c(a = 1L, b = 2L)), 3L)
})

test_that("extended: random vectors agree with the oracle and base R", {
  skip_on_cran()
  skip_unless_extended()
  for (k in 1:200) {
    spec <- with_seed(k, list(
      type = sample(c("double", "integer", "logical"), 1),
      n = sample(c(sample.int(300, 1), sample.int(1e5, 1)), 1),
      na_frac = sample(c(0, 0.01, 0.2), 1),
      mode = sample(modes, 1),
      na_rm = sample(c(FALSE, TRUE), 1)
    ))
    x <- rand_vec(spec$type, spec$n, na_frac = spec$na_frac, nan_frac = 0, seed = k)
    expect_simd_equal(simd_sum, x, na.rm = spec$na_rm, precision = spec$mode)
    expect_simd_matches_base(simd_sum, sum, x, na.rm = spec$na_rm, precision = spec$mode)
  }
})
