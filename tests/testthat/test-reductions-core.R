# Core reductions other than simd_sum (which has test-sum.R): prod, mean,
# min, max, range, which_min, which_max, any, all, any_na, count_na,
# which_na. Every tier against the none oracle and against base R.

modes <- c("fast", "pairwise", "compensated")

# Base R's which.min/which.max/any/all/anyNA under the simd_* names.
base_which_na <- function(x) which(is.na(x))
base_count_na <- function(x) as.double(sum(is.na(x)))

# Products of many random doubles over- or underflow depending on the order
# of multiplication, so products use values near 1 (plus specials).
near_one <- function(n, seed, na_frac = 0.05, nan_frac = 0.02) {
  x <- 1 + rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = seed) / 1e4
  with_seed(seed + 1L, {
    x[sample.int(n, stats::rbinom(1, n, nan_frac))] <- NaN
    x[sample.int(n, stats::rbinom(1, n, na_frac))] <- NA
  })
  x
}

# Relative tolerance of a product of n doubles (one rounding per factor).
prod_tol <- function(n) max(n, 1) * 2^-52 * 2

# Expects simd_mean(x, ...) within the precision-mode bound of base mean(),
# on every tier.
expect_mean_matches_base <- function(x, ..., precision) {
  n <- max(sum(!is.na(x)), 1)
  bound <- (mode_bound(x, precision) + length(x) * 2^-52 * abs_mass(x)) / n
  expect_simd_matches_base(simd_mean, mean, x, ...,
    tolerance = NULL, precision = precision
  )
  got <- with_each_tier(function() {
    old <- simd_precision(precision)
    on.exit(simd_precision(old))
    simd_mean(x, ...)
  })
  want <- mean(x, ...)
  for (tier in names(got)) {
    if (is.finite(want)) {
      check_true(abs(got[[tier]] - want) <= bound + 2^-52 * abs(want),
        info = paste(tier, precision, length(x))
      )
    }
  }
}

test_that("every reduction matches the oracle on every tier for all edge lengths", {
  batch_expectations({
    for (n in sweep_lengths()) {
      seed <- as.integer(n %% 1000) + 1L
      d <- rand_vec("double", n, seed = seed)
      i <- rand_vec("integer", n, seed = seed)
      l <- rand_vec("logical", n, seed = seed)
      p <- near_one(n, seed)
      inputs <- list(
        d, rand_vec("double", n, na_frac = 0, nan_frac = 0, seed = seed), i,
        rand_vec("integer", n, na_frac = 0, seed = seed), l
      )
      for (x in sweep_inputs(n, inputs)) {
        for (na_rm in c(FALSE, TRUE)) {
          suppressWarnings({
            expect_simd_identical(simd_min, x, na.rm = na_rm)
            expect_simd_identical(simd_max, x, na.rm = na_rm)
            expect_simd_identical(simd_range, x, na.rm = na_rm)
          })
          expect_simd_identical(simd_any, x != 0, na.rm = na_rm)
          expect_simd_identical(simd_all, x != 0, na.rm = na_rm)
          for (mode in if (is.double(x)) modes else "fast") {
            expect_simd_equal(simd_mean, x, na.rm = na_rm, precision = mode)
          }
          if (!is.double(x)) {
            x25 <- x[seq_len(min(n, 25))]
            expect_simd_equal(simd_prod, x25, na.rm = na_rm, tolerance = prod_tol(25))
          }
        }
        expect_simd_identical(simd_which_min, x)
        expect_simd_identical(simd_which_max, x)
        expect_simd_identical(simd_any_na, x)
        expect_simd_identical(simd_count_na, x)
        expect_simd_identical(simd_which_na, x)
        if (!anyNA(x)) {
          suppressWarnings(expect_simd_identical(simd_range, x, na_check = FALSE))
          expect_simd_equal(simd_mean, x, na_check = FALSE)
        }
      }
      suppressWarnings({
        expect_simd_identical(simd_any, d)
        expect_simd_identical(simd_all, d, na.rm = TRUE)
        expect_simd_identical(simd_any, i)
        expect_simd_identical(simd_all, i)
      })
      for (na_rm in c(FALSE, TRUE)) {
        expect_simd_equal(simd_prod, p, na.rm = na_rm, tolerance = prod_tol(n))
      }
    }
  })
})

test_that("reductions agree with base R across edge lengths", {
  batch_expectations({
    for (n in sweep_lengths()) {
      seed <- as.integer(n %% 1000) + 1L
      # One kind of missing value at a time, where base R is deterministic.
      xs <- list(
        rand_vec("double", n, nan_frac = 0, seed = seed),
        rand_vec("double", n, na_frac = 0, seed = seed),
        rand_vec("integer", n, seed = seed),
        rand_vec("logical", n, seed = seed)
      )
      for (x in sweep_inputs(n, xs)) {
        for (na_rm in c(FALSE, TRUE)) {
          suppressWarnings({
            expect_simd_matches_base(simd_min, min, x, na.rm = na_rm)
            expect_simd_matches_base(simd_max, max, x, na.rm = na_rm)
            expect_simd_matches_base(simd_range, range, x, na.rm = na_rm)
            expect_simd_matches_base(simd_any, any, x != 0, na.rm = na_rm)
            expect_simd_matches_base(simd_all, all, x != 0, na.rm = na_rm)
          })
          for (mode in if (is.double(x)) modes else "fast") {
            expect_mean_matches_base(x, na.rm = na_rm, precision = mode)
          }
        }
        expect_simd_matches_base(simd_which_min, which.min, x)
        expect_simd_matches_base(simd_which_max, which.max, x)
        expect_simd_matches_base(simd_any_na, anyNA, x)
        expect_simd_matches_base(simd_count_na, base_count_na, x)
        expect_simd_matches_base(simd_which_na, base_which_na, x)
      }
      p <- near_one(n, seed, nan_frac = 0)
      for (na_rm in c(FALSE, TRUE)) {
        expect_simd_matches_base(simd_prod, prod, p, na.rm = na_rm, tolerance = prod_tol(n))
      }
      x <- rand_vec("integer", min(n, 30), seed = seed) %/% 1000L
      expect_simd_matches_base(simd_prod, prod, x, na.rm = TRUE, tolerance = prod_tol(30))
    }
  })
})

test_that("the base R behaviour table is reproduced", {
  batch_expectations({
    # Each case: simd call, base call. Results and warnings must be identical
    # on every tier.
    cases <- list(
      quote(prod(integer(0))), quote(prod(1:3)), quote(min(integer(0))),
      quote(max(numeric(0))), quote(range(integer(0))), quote(mean(numeric(0))),
      quote(mean(integer(0))), quote(min(c(0, -0))), quote(min(c(-0, 0))),
      quote(max(c(0, -0))), quote(max(c(-0, 0))), quote(max(c(1, NaN))), quote(max(c(NaN, 1))),
      quote(min(c(NA, NaN))), quote(min(c(NaN, NA))), quote(mean(c(NA, NaN))),
      quote(prod(c(NA, NaN))), quote(any(as.raw(1))), quote(all(as.raw(c(1, 0)))),
      quote(any(c(0, 2.5))), quote(any(c(1L, 0L))), quote(any(c(NaN, 0))),
      quote(any(c(NA, TRUE))), quote(any(c(NA, FALSE))), quote(all(c(NA, TRUE))),
      quote(all(c(NA, FALSE))), quote(any(logical(0))), quote(all(logical(0))),
      quote(any(numeric(0))), quote(all(raw(0))), quote(which.max(c(1, 3, 3))),
      quote(which.max(c(NA, 3, NaN, 3))), quote(which.max(numeric(0))),
      quote(which.max(c(NA, NA))), quote(which.max(c(0, -0))), quote(which.min(c(0, -0))),
      quote(which.min(c(TRUE, FALSE))), quote(anyNA(c(1, NaN))), quote(anyNA(as.raw(1))),
      quote(which.max(as.raw(1:3))), quote(which.min(as.raw(c(3, 1, 1)))),
      quote(which.max(raw(0))), quote(min(c(TRUE, FALSE))),
      quote(range(c(TRUE, NA), na.rm = TRUE)), quote(min(NA_integer_, na.rm = TRUE)),
      quote(range(NA_real_, na.rm = TRUE)), quote(range(c(1, NaN, 3))),
      quote(range(c(NaN, 1, NA))), quote(mean(c(TRUE, FALSE, TRUE))),
      quote(prod(c(TRUE, NA))), quote(mean(c(1L, NA), na.rm = TRUE)),
      quote(which.max(c(-Inf, NA))), quote(which.min(c(Inf, Inf))),
      quote(which.min(c(-Inf, NA, -Inf))), quote(max(c(-Inf, NaN))),
      quote(min(c(Inf, NaN), na.rm = TRUE)), quote(all(c(NaN, 1))),
      quote(all(c(NA_integer_, 0L), na.rm = TRUE)), quote(any(c(NA_real_, 1), na.rm = TRUE)),
      quote(mean(c(.Machine$integer.max, .Machine$integer.max))), quote(prod(c(TRUE, TRUE))),
      quote(min(c(1, NA, NaN))), quote(max(c(NaN, NA, 1))),
      quote(1 / min(c(0, -0, NA), na.rm = TRUE)), quote(sum(c(NA, NaN)))
    )
    simd_name <- c(
      prod = "simd_prod", min = "simd_min", max = "simd_max", range = "simd_range",
      mean = "simd_mean", any = "simd_any", all = "simd_all", which.max = "simd_which_max",
      which.min = "simd_which_min", anyNA = "simd_any_na", sum = "simd_sum"
    )
    to_simd <- function(e) {
      if (is.call(e)) {
        f <- as.character(e[[1L]])
        if (f %in% names(simd_name)) e[[1L]] <- as.name(simd_name[[f]])
        e[-1L] <- lapply(as.list(e[-1L]), to_simd)
      }
      e
    }
    # The value and the warnings of evaluating e in env.
    outcome <- function(e, env) {
      w <- character(0)
      r <- withCallingHandlers(eval(e, env), warning = function(c) {
        w <<- c(w, conditionMessage(c))
        invokeRestart("muffleWarning")
      })
      list(r, w)
    }
    here <- environment()
    for (e in cases) {
      want <- outcome(e, baseenv())
      simd_e <- to_simd(e)
      got <- with_each_tier(function() outcome(simd_e, here))
      for (tier in names(got)) {
        g <- got[[tier]]
        check_identical(g, want, info = paste(tier, deparse1(e)))
        check_identical(g[[1L]], want[[1L]], info = paste(tier, deparse1(e)), num.eq = FALSE)
      }
    }
  })
})

test_that("NA beats NaN in either order for prod, mean, min, max and range", {
  batch_expectations({
    for (x in list(c(NA, NaN), c(NaN, NA), c(1, NaN, 2, NA), c(NA, 1:40, NaN), c(NaN, 1:40, NA))) {
      for (f in list(simd_prod, simd_mean, simd_min, simd_max, simd_range)) {
        res <- with_each_tier(function() f(x))
        for (tier in names(res)) {
          check_true(all(is.na(res[[tier]]) & !is.nan(res[[tier]])), info = tier)
        }
      }
    }
  })
})

test_that("a zero extremum has the sign of the first zero, wherever the zeros are", {
  batch_expectations({
    for (n in c(2, 3, 5, 9, 17, 33, 65, 130, 300)) {
      for (k in 1:12) {
        pos <- with_seed(n * 100 + k, sort(sample.int(n, min(n, 4))))
        x <- with_seed(k, stats::runif(n, 1, 2))
        signs <- with_seed(k + 7, sample(c(0, -0), length(pos), replace = TRUE))
        x[pos] <- signs
        lo <- expect_simd_identical(simd_min, x)
        for (tier in names(lo)) {
          check_identical(lo[[tier]], min(x), info = paste(tier, n, k), num.eq = FALSE)
        }
        y <- -x
        hi <- expect_simd_identical(simd_max, y)
        for (tier in names(hi)) {
          check_identical(hi[[tier]], max(y), info = paste(tier, n, k), num.eq = FALSE)
        }
        r <- expect_simd_identical(simd_range, c(x, y))
        for (tier in names(r)) {
          check_identical(r[[tier]], range(c(x, y)), info = tier, num.eq = FALSE)
        }
        expect_simd_identical(simd_which_min, x)
        expect_simd_matches_base(simd_which_min, which.min, x)
        expect_simd_matches_base(simd_which_max, which.max, y)
      }
    }
  })
})

test_that("prod is exact about zeros, infinities and out-of-range partial products", {
  # Partial products in different lanes overflow and underflow; a zero in
  # one lane meeting an Inf in another must not give NaN.
  x <- with_seed(42L, round(stats::rnorm(2^20 + 1) * 100))
  big <- c(rep(1e300, 40), -2, rep(1e-300, 40))
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_identical(simd_prod(x), 0, info = tier)
      expect_identical(simd_prod(as.integer(x)), 0, info = tier)
      expect_identical(simd_prod_diffs(x, 0L), 0, info = tier)
      expect_identical(simd_prod_sums(as.integer(x), x), 0, info = tier)
      expect_equal(simd_prod(c(1e300, 1e300, 1e-300, 1e-300)), 1, tolerance = 4 * 2^-52,
                   info = tier)
      expect_equal(simd_prod(big), -2, tolerance = prod_tol(81), info = tier)
      expect_equal(simd_prod_sums(big, 0), -2, tolerance = prod_tol(81), info = tier)
      # A zero beats a product that leaves the range of base R's long double too.
      expect_identical(1 / simd_prod(c(rep(1e300, 40), -0, 3)), -Inf, info = tier)
      expect_identical(1 / simd_prod(c(rep(-1e300, 41), 0)), -Inf, info = tier)
      expect_identical(simd_prod(c(rep(1e300, 40), 0, Inf)), NaN, info = tier)
      expect_identical(simd_prod(c(rep(1e300, 40), -Inf)), -Inf, info = tier)
      expect_identical(simd_prod(rep(-1e300, 41)), -Inf, info = tier)
      expect_identical(simd_prod(c(rep(1e-300, 40), 1e300)), 0, info = tier)
      # Missing values decide the result as before, or are removed.
      expect_identical(simd_prod(c(rep(1e300, 40), NaN, 0)), NaN, info = tier)
      expect_identical(simd_prod(c(rep(1e300, 40), NA, 0, NaN)), NA_real_, info = tier)
      expect_identical(simd_prod(c(rep(1e300, 40), NA, 0, NaN), na.rm = TRUE), 0, info = tier)
      expect_identical(simd_prod(c(NA, 0L, -3L), na.rm = TRUE), -0, info = tier)
      expect_identical(simd_prod_diffs(c(rep(1e300, 40), NA, 1), 1, na.rm = TRUE), 0,
                       info = tier)
      # Compact sequences go through the chunked (ALTREP) path.
      expect_identical(simd_prod(-5:100000), -0, info = tier)
    })
  }
})

test_that("which_*: first of ties, NA and NaN ignored, empty gives integer(0)", {
  batch_expectations({
    cases <- list(
      c(1, 3, 3), c(3, 1, 3, 1), c(NA, 3, NaN, 3), c(NaN, NaN), c(NA, NA), numeric(0),
      c(rep(NA, 40), 2, 2), c(rep(5L, 70), 9L, 9L), c(NA_integer_, NA_integer_), integer(0),
      c(FALSE, NA, TRUE, TRUE), c(-Inf, -Inf, NA), c(Inf, NaN, Inf), as.raw(c(0, 255, 255, 0))
    )
    for (x in cases) {
      expect_simd_matches_base(simd_which_min, which.min, x)
      expect_simd_matches_base(simd_which_max, which.max, x)
      expect_simd_identical(simd_which_min, x)
      expect_simd_identical(simd_which_max, x)
    }
    # A tie between lanes, chunks and vector blocks: the first wins.
    for (n in c(7, 64, 1000, 5000)) {
      x <- rep(1, n)
      x[c(n, n %/% 2, 3)] <- 0
      expect_simd_matches_base(simd_which_min, which.min, x)
      expect_simd_matches_base(simd_which_max, which.max, -x)
      xi <- as.integer(x)
      expect_simd_matches_base(simd_which_min, which.min, xi)
    }
  })
})

test_that("compact sequences give the same results as expanded vectors", {
  batch_expectations({
    for (n in c(2, 4097, long_length())) {
      for (x in altrep_inputs(n)) {
        check_true(takes_region_path(x))
        expanded <- x + if (is.integer(x)) 0L else 0
        expect_false(takes_region_path(expanded))
        for (f in list(
          simd_min, simd_max, simd_range, simd_which_min, simd_which_max,
          simd_any_na, simd_count_na, simd_which_na
        )) {
          a <- expect_simd_identical(f, x)
          b <- expect_simd_identical(f, expanded)
          check_identical(a, b)
        }
        for (mode in modes) {
          a <- expect_simd_equal(simd_mean, x, precision = mode)
          b <- expect_simd_equal(simd_mean, expanded, precision = mode)
          check_identical(a, b, info = mode)
        }
        check_identical(simd_which_max(x), which.max(expanded))
        check_identical(simd_mean(x), mean(expanded))
      }
    }
  })
})

test_that("compact sequences of every sign and direction reduce as expanded vectors", {
  fns <- list(
    simd_sum, simd_mean, simd_min, simd_max, simd_range, simd_which_min, simd_which_max,
    simd_min_abs, simd_max_abs, simd_which_min_abs, simd_which_max_abs, simd_sum_sq,
    simd_norm, simd_sum_abs
  )
  seqs <- list(
    1:10, 10:1, -5:5, 5:-5, -3:3, 3:-3, -7:2, 2:-7, -10:-1, -1:-10, 0:5, -5:0, 5:0,
    -50000:49999, -2147483647:-2147480000, 2147480000:2147483647
  )
  seqs <- c(seqs, lapply(seqs, as.double))
  batch_expectations({
    for (x in seqs) {
      expanded <- x + if (is.integer(x)) 0L else 0
      for (k in seq_along(fns)) {
        check_identical(fns[[k]](x), fns[[k]](expanded), info = c(k, x[1], length(x)))
      }
    }
  })
})

test_that("compact sequences are reduced from their endpoints, sums rounded once", {
  skip_if(getRversion() < "4.6.0", "R_altrep_class_name() is new in R 4.6.0")
  # 1:3e9 is a double sequence; reading its 3e9 elements would take seconds.
  x <- 1:3e9
  expect_identical(simd_sum(x), 0x1.f399b14655518p+61)
  expect_identical(simd_sum_sq(x), 0x1.d14a021dcc7b8p+92)
  expect_identical(simd_mean(x), 1500000000.5)
  expect_identical(simd_range(x), c(1, 3e9))
  expect_identical(simd_which_max(x), 3e9)
  expect_identical(simd_which_min_abs(-3e9:5), 3e9 + 1)
  expect_identical(simd_sum(-1e15:1e15), 0)
  expect_identical(simd_sum_abs(-1e15:1e15), 0x1.93e5939a08cf1p+99)
  # Squares near 2^62, whose running sum rounds.
  x <- 2147483000:2147483647
  expect_identical(simd_sum_sq(x), 0x1.43fff995380adp+71)
  expect_identical(simd_norm(x), 0x1.974b1f2c64b9dp+35)
})

test_that("na.rm = TRUE on all-missing input follows base R", {
  for (x in list(NA_real_, c(NA, NaN), rep(NA_integer_, 20), c(NA, NA))) {
    res <- with_each_tier(function() {
      list(
        sum = simd_sum(x, na.rm = TRUE), prod = simd_prod(x, na.rm = TRUE),
        mean = simd_mean(x, na.rm = TRUE),
        min = suppressWarnings(simd_min(x, na.rm = TRUE)),
        max = suppressWarnings(simd_max(x, na.rm = TRUE))
      )
    })
    want <- list(
      sum = if (is.double(x)) 0 else 0L, prod = 1, mean = NaN, min = Inf, max = -Inf
    )
    for (tier in names(res)) expect_identical(res[[tier]], want, info = tier)
  }
})

test_that("na_check = FALSE gives the same results on NA-free input", {
  x <- rand_vec("integer", 5000, na_frac = 0, seed = 2)
  d <- rand_vec("double", 5000, na_frac = 0, nan_frac = 0, seed = 2)
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      for (f in list(simd_min, simd_max, simd_range, simd_mean)) {
        expect_identical(f(x, na_check = FALSE), f(x), info = tier)
        expect_identical(f(d, na_check = FALSE), f(d), info = tier)
      }
      expect_identical(simd_prod(x[1:20], na_check = FALSE), simd_prod(x[1:20]), info = tier)
      # A broken promise gives an unspecified value, but no error.
      expect_type(simd_min(c(1L, NA), na_check = FALSE), "integer")
      expect_type(simd_mean(c(1, NA), na_check = FALSE), "double")
      # na.rm = TRUE always looks for missing values.
      expect_identical(simd_min(c(3L, NA), na.rm = TRUE, na_check = FALSE), 3L, info = tier)
    })
  }
})

test_that("mean: pairwise and compensated refine with a second pass", {
  x <- c(1e16, rep(1, 1000), -1e16 + 2)
  old <- simd_precision("compensated")
  on.exit(simd_precision(old))
  res <- with_each_tier(function() simd_mean(x))
  # The exact mean, 1002 / 1002; base R's is 0 where long double is double.
  for (tier in names(res)) {
    expect_equal(res[[tier]], 1, tolerance = 1e-12, info = tier)
  }
  # Integer means are exact sums divided once.
  skip_if_not(has_wide_long_double(), "no long double accumulator in base R")
  for (x in list(rep(.Machine$integer.max, 7L), c(-5L, 7L, 9L), 1:1e6)) {
    expect_simd_matches_base(simd_mean, mean, x, tolerance = 0)
    res <- with_each_tier(function() simd_mean(x))
    for (tier in names(res)) expect_identical(res[[tier]], mean(x), info = tier)
  }
})

test_that("raw any/all find a zero or non-zero byte at every position", {
  batch_expectations({
    for (n in c(1:70, 127, 128, 129, 1000)) {
      ones <- as.raw(rep(c(1, 128, 255, 129), length.out = n))
      expect_simd_identical(function(x) suppressWarnings(simd_all(x)), ones)
      for (p in unique(c(1, max(n %/% 2, 1), n, seq_len(min(n, 40))))) {
        x <- ones
        x[p] <- as.raw(0)
        res <- with_each_tier(function() suppressWarnings(c(simd_all(x), simd_any(x))))
        want <- c(FALSE, n > 1)
        for (tier in names(res)) check_identical(res[[tier]], want, info = paste(tier, n, p))
        z <- raw(n)
        z[p] <- as.raw(1)
        res <- with_each_tier(function() suppressWarnings(c(simd_all(z), simd_any(z))))
        want <- c(n == 1, TRUE)
        for (tier in names(res)) check_identical(res[[tier]], want, info = paste(tier, n, p))
      }
      res <- with_each_tier(function() suppressWarnings(c(simd_all(ones), simd_any(raw(n)))))
      for (tier in names(res)) check_identical(res[[tier]], c(TRUE, FALSE), info = paste(tier, n))
    }
  })
})

test_that("missing-value queries cover doubles, integers, logicals, complex and raw", {
  batch_expectations({
    for (n in c(1, 2, 3, 7, 33, 129, 4097)) {
      for (k in 1:3) {
        d <- rand_vec("double", n, na_frac = 0.1, nan_frac = 0.1, seed = n + k)
        z <- complex(real = d, imaginary = rev(d))
        for (x in list(d, rand_vec("integer", n, na_frac = 0.2, seed = k), z)) {
          expect_simd_matches_base(simd_any_na, anyNA, x)
          expect_simd_matches_base(simd_count_na, base_count_na, x)
          expect_simd_matches_base(simd_which_na, base_which_na, x)
        }
      }
    }
    check_identical(simd_any_na(as.raw(1:3)), FALSE)
    check_identical(simd_count_na(as.raw(1:3)), 0)
    check_identical(simd_which_na(as.raw(1:3)), integer(0))
    check_identical(simd_which_na(numeric(0)), integer(0))
    check_identical(simd_count_na(integer(0)), 0)
    check_identical(simd_any_na(complex(real = 1, imaginary = NaN)), TRUE)
  })
})

test_that("chunked results equal unchunked ones", {
  skip_on_cran()
  skip_if_no_subprocess()
  n <- 2^20 + 7
  d <- rand_vec("double", n, na_frac = 0.001, nan_frac = 0.001, seed = 9)
  z <- rep(1, n)
  z[c(200000, 3000, 900000)] <- c(0, -0, 0) # the first zero is -0
  i <- rand_vec("integer", n, na_frac = 0.001, seed = 9)
  clean <- rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = 9)
  child <- function(d, z, i, clean, modes) {
    library(rsimd)
    lapply(stats::setNames(nm = simd_available()), function(tier) {
      simd_with_impl(tier, list(
        exact = suppressWarnings(list(
          simd_range(d, na.rm = TRUE), simd_range(z), simd_range(i, na.rm = TRUE),
          simd_which_min(d), simd_which_max(i), simd_which_min(z), simd_any(i > 999000L),
          simd_all(i > -999000L, na.rm = TRUE), simd_any(d > 1e300), simd_any_na(d),
          simd_count_na(d), simd_which_na(i), simd_min(d), simd_range(i),
          simd_mean(i, na.rm = TRUE)
        )),
        mean = vapply(modes, function(m) {
          old <- simd_precision(m)
          on.exit(simd_precision(old))
          simd_mean(clean)
        }, numeric(1)),
        prod = simd_prod(1 + clean / 1e8)
      ))
    })
  }
  here <- child(d, z, i, clean, modes)
  env <- c(callr::rcmd_safe_env(), RSIMD_DEBUG_STRIDE = "128")
  there <- callr::r(child, list(d, z, i, clean, modes), env = env)
  for (tier in intersect(names(here), names(there))) {
    expect_true(identical(there[[tier]]$exact, here[[tier]]$exact, num.eq = FALSE), info = tier)
    expect_identical(1 / here[[tier]]$exact[[2]][1], -Inf, info = tier)
    expect_identical(there[[tier]]$mean[["pairwise"]], here[[tier]]$mean[["pairwise"]])
    expect_lte(abs(there[[tier]]$mean[["compensated"]] - here[[tier]]$mean[["compensated"]]),
      2 * ulp(here[[tier]]$mean[["compensated"]]),
      label = tier
    )
    expect_lte(abs(there[[tier]]$mean[["fast"]] - here[[tier]]$mean[["fast"]]),
      mode_bound(clean, "fast") / n,
      label = tier
    )
    expect_equal(there[[tier]]$prod, here[[tier]]$prod, tolerance = prod_tol(n), info = tier)
  }
  # And the unchunked tiers agree with the oracle.
  for (tier in names(here)) {
    expect_true(identical(here[[tier]]$exact, here[["none"]]$exact, num.eq = FALSE), info = tier)
  }
})

test_that("warnings are emitted once per call, as base R", {
  count_warnings <- function(expr) {
    k <- 0L
    withCallingHandlers(expr, warning = function(w) {
      k <<- k + 1L
      invokeRestart("muffleWarning")
    })
    k
  }
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_identical(count_warnings(simd_min(numeric(0))), 1L, info = tier)
      expect_identical(count_warnings(simd_max(NA_integer_, na.rm = TRUE)), 1L, info = tier)
      expect_identical(count_warnings(simd_range(integer(0))), 2L, info = tier)
      expect_identical(count_warnings(simd_any(as.double(0:5000))), 1L, info = tier)
      expect_identical(count_warnings(simd_all(as.raw(1:200))), 1L, info = tier)
      expect_identical(count_warnings(simd_any(numeric(0))), 0L, info = tier)
      expect_identical(count_warnings(simd_any(1:3)), 0L, info = tier)
    })
  }
  expect_warning(simd_any(c(0, 1)), "coercing argument of type 'double' to logical", fixed = TRUE)
  expect_warning(simd_all(as.raw(1)), "coercing argument of type 'raw' to logical", fixed = TRUE)
})

test_that("complex any/all read values as base R does, with its warning", {
  cplx <- function(r, i) complex(real = r, imaginary = i)
  vals <- list(
    0i, 1i, cplx(1, 0), cplx(-0, -0), NA_complex_, cplx(NaN, 0), cplx(0, NaN),
    cplx(NA, 1), cplx(NaN, 1), cplx(Inf, 0), c(0i, NA_complex_), c(1i, NA_complex_),
    c(NA_complex_, 0i, 1i), complex(0), c(rep(1i, 5000), 0i), c(rep(0i, 5000), 1i)
  )
  for (v in vals) {
    for (rm in c(FALSE, TRUE)) {
      for (f in list(c(simd_any, any), c(simd_all, all))) {
        expect_tiers_give(suppressWarnings(f[[2]](v, na.rm = rm)),
          function(...) suppressWarnings(f[[1]](...)), v, na.rm = rm)
      }
    }
  }
  expect_warning(simd_any(1i), "coercing argument of type 'complex' to logical", fixed = TRUE)
  expect_warning(simd_all(c(0i, 1i)), "coercing argument of type 'complex' to logical",
    fixed = TRUE)
  expect_no_warning(simd_any(complex(0)))
  expect_identical(suppressWarnings(simd_any(simd_vec(c(0i, 1i)))), TRUE)
})

test_that("arguments and types are validated", {
  for (f in c("simd_prod", "simd_mean", "simd_min", "simd_max", "simd_range")) {
    expect_error(get(f)(as.raw(1)), paste0(f, "() does not support 'x' of type raw"),
      fixed = TRUE
    )
  }
  for (f in c("simd_min", "simd_max", "simd_range")) {
    expect_error(get(f)(1i), paste0(f, "() does not support 'x' of type complex"), fixed = TRUE)
  }
  no_complex <- c("simd_which_min", "simd_which_max")
  for (f in no_complex) {
    msg <- tryCatch(get(f)(1i), error = conditionMessage)
    expect_identical(msg, paste0(f, "() does not support 'x' of type complex"))
  }
  for (f in c("simd_min", "simd_mean", "simd_which_max", "simd_any", "simd_count_na")) {
    expect_error(get(f)("a"), "'x' must be an atomic vector", fixed = TRUE)
    expect_error(get(f)(list(1)), "not list", fixed = TRUE)
  }
  expect_error(simd_min(1, na.rm = NA), "'na.rm' must be TRUE or FALSE", fixed = TRUE)
  expect_error(simd_any(TRUE, na.rm = "yes"), "'na.rm' must be TRUE or FALSE", fixed = TRUE)
  expect_error(simd_mean(1, na_check = 1:2), "'na_check' must be TRUE or FALSE", fixed = TRUE)
  skip_if_not_installed("bit64")
  x <- bit64::as.integer64(1:3)
  for (f in c("simd_prod", "simd_mean")) {
    expect_error(get(f)(x), paste0(f, "() does not support 'x' of type integer64 yet"),
      fixed = TRUE
    )
  }
})

# Expects every function with na.rm to error on an extra positional
# argument, naming the user's call, rather than read it as na.rm:
# simd_any(FALSE, TRUE) is not any(FALSE, TRUE).
test_that("na.rm, finite and na_check must be named", {
  one <- c(
    "simd_sum", "simd_prod", "simd_mean", "simd_min", "simd_max", "simd_range",
    "simd_max_abs", "simd_min_abs", "simd_sum_sq", "simd_norm", "simd_sum_abs",
    "simd_var", "simd_sd", "simd_popcount_total", "simd_any", "simd_all", "simd_count"
  )
  two <- c(
    "simd_prod_sums", "simd_prod_diffs", "simd_dot", "simd_dist", "simd_hamming",
    "simd_pmin", "simd_pmax"
  )
  for (f in c(one, two)) {
    fn <- get(f)
    named <- if ("finite" %in% names(formals(fn))) {
      "na.rm, finite and na_check"
    } else if ("na_check" %in% names(formals(fn))) {
      "na.rm and na_check"
    } else {
      "na.rm"
    }
    what <- if (f %in% one) {
      "one vector; combine several with c()"
    } else if (f %in% c("simd_pmin", "simd_pmax")) {
      "two vectors; nest calls for more"
    } else {
      "two vectors"
    }
    args <- if (f %in% one) list(TRUE, TRUE) else list(TRUE, TRUE, TRUE)
    call <- as.call(c(as.name(f), args))
    cnd <- tryCatch(eval(call), error = identity)
    expect_identical(
      conditionMessage(cnd), paste0(f, "() takes ", what, ", and ", named, " must be named"),
      info = f
    )
    expect_identical(conditionCall(cnd), call, info = f)
    expect_error(do.call(fn, c(args[-1L], na = TRUE)), "must be named", info = f)
    expect_no_error(do.call(fn, c(args[-1L], na.rm = TRUE)))
  }
})

test_that("attributes of x do not matter", {
  m <- matrix(c(3, NA, 1, 2), 2, dimnames = list(c("a", "b"), NULL))
  expect_identical(simd_min(m, na.rm = TRUE), 1)
  expect_identical(simd_which_max(m), 1L)
  expect_identical(simd_which_na(m), 2L)
  expect_identical(simd_range(c(a = 1L, b = 2L)), c(1L, 2L))
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
    suppressWarnings({
      expect_simd_identical(simd_range, x, na.rm = spec$na_rm)
      expect_simd_matches_base(simd_range, range, x, na.rm = spec$na_rm)
    })
    expect_simd_matches_base(simd_which_min, which.min, x)
    expect_simd_matches_base(simd_which_max, which.max, x)
    expect_simd_matches_base(simd_which_na, base_which_na, x)
    expect_simd_matches_base(simd_any, any, x > 0, na.rm = spec$na_rm)
    expect_simd_matches_base(simd_all, all, x > 0, na.rm = spec$na_rm)
    expect_mean_matches_base(x, na.rm = spec$na_rm, precision = spec$mode)
    expect_simd_equal(simd_mean, x, na.rm = spec$na_rm, precision = spec$mode)
  }
})

test_that("extended: long vectors", {
  skip_on_cran()
  skip_unless_extended()
  n <- 2^31 + 10
  x <- seq_len(n) # a compact double sequence, never expanded
  expect_true(takes_region_path(x))
  expect_identical(simd_count_na(x), 0)
  expect_false(simd_any_na(x))
  expect_identical(simd_which_na(x), integer(0))
  expect_identical(simd_which_max(x), n)
  expect_identical(simd_which_min(x), 1)
})

test_that("min and max skip NA (a signalling NaN) and NaN with na.rm on every tier", {
  set.seed(11)
  for (n in c(5L, 9L, 17L, 33L, 5000L)) {
    x <- runif(n, -100, 100)
    for (pos in unique(c(1L, n %/% 2L + 1L, n))) {
      for (miss in list(NA_real_, NaN)) {
        y <- x
        y[pos] <- miss
        expect_tiers_give(min(y, na.rm = TRUE), function(v) simd_min(v, na.rm = TRUE), y)
        expect_tiers_give(max(y, na.rm = TRUE), function(v) simd_max(v, na.rm = TRUE), y)
        expect_tiers_give(max(abs(y), na.rm = TRUE), function(v) simd_max_abs(v, na.rm = TRUE), y)
        expect_tiers_give(which.min(y), simd_which_min, y)
        expect_tiers_give(min(y), simd_min, y)
      }
    }
  }
})

test_that("min and max find NA and NaN anywhere, and Inf + -Inf is not taken for NaN", {
  # The vector tiers test the sum of each four vectors for NaN, which
  # Inf + -Inf in a lane (or finite values overflowing into it) also give.
  # Base R function, rsimd function, whether both take na.rm.
  fs <- list(
    list(min, simd_min, TRUE), list(max, simd_max, TRUE), list(range, simd_range, TRUE),
    list(function(v, ...) max(abs(v), ...), simd_max_abs, TRUE),
    list(function(v, ...) min(abs(v), ...), simd_min_abs, TRUE),
    list(which.min, simd_which_min, FALSE), list(which.max, simd_which_max, FALSE),
    list(function(v) which.min(abs(v)), simd_which_min_abs, FALSE),
    list(function(v) which.max(abs(v)), simd_which_max_abs, FALSE)
  )
  expect_all <- function(y) {
    for (f in fs) {
      expect_tiers_give(f[[1]](y), f[[2]], y)
      if (f[[3]]) {
        expect_tiers_give(f[[1]](y, na.rm = TRUE), function(v) f[[2]](v, na.rm = TRUE), y)
      }
    }
  }
  for (n in c(64L, 101L, 3e4L + 3L)) {
    inf <- with_seed(n, sample(c(Inf, -Inf, 1e308, -1e308, runif(4)), n, TRUE))
    expect_all(inf)
    for (f in fs[1:3]) {
      expect_tiers_give(f[[1]](inf), function(v) f[[2]](v, na_check = FALSE), inf)
    }
    # A NaN or an NA in each of the first four vectors of a block, at a
    # block edge, in a later block and in the tail.
    for (pos in unique(pmin(c(1:33, 4096:4097, 20000L, n - 0:9), n))) {
      for (miss in list(NaN, NA_real_)) {
        y <- inf
        y[pos] <- miss
        expect_all(y)
        # NA wins over a NaN before or after it.
        y[if (pos > 1L) 1L else n] <- if (is.nan(miss)) NA_real_ else NaN
        expect_tiers_give(NA_real_, simd_min, y)
      }
    }
  }
})

test_that("finite = TRUE drops infinities and missing values on every tier", {
  fin <- function(v) v[is.finite(v)]
  # Base R function of the finite values, rsimd function.
  fs <- list(
    list(min, simd_min), list(max, simd_max), list(range, simd_range),
    list(function(v) max(abs(v)), simd_max_abs),
    list(function(v) min(abs(v)), simd_min_abs)
  )
  for (n in c(6L, 37L, 3e4L + 3L)) {
    x <- with_seed(n, runif(n, -100, 100))
    edges <- unique(pmin(c(1L, n, 1024:1025, n - 1L), n))
    # Infinities at the ends, at a block edge and in the tail; an infinity
    # of each sign alone; missing values among them; no infinity at all.
    cases <- list(
      replace(x, edges, rep_len(c(Inf, -Inf), length(edges))),
      replace(x, n %/% 2L + 1L, Inf), replace(x, 1L, -Inf),
      replace(x, unique(c(2L, n)), c(NA, NaN)),
      replace(x, unique(c(1L, 3L, n)), c(Inf, NA, -Inf)), x
    )
    for (y in cases) {
      for (f in fs) {
        expect_tiers_give(f[[1]](fin(y)), function(v) f[[2]](v, finite = TRUE), y)
        # finite implies na.rm, as in base R.
        expect_tiers_give(f[[1]](fin(y)), function(v) f[[2]](v, na.rm = FALSE, finite = TRUE), y)
      }
      expect_tiers_give(range(y, finite = TRUE), function(v) simd_range(v, finite = TRUE), y)
    }
  }
  # Nothing finite left: base R's values and warnings.
  for (y in list(c(Inf, -Inf, NA), c(Inf, Inf), numeric())) {
    expect_identical(
      suppressWarnings(simd_range(y, finite = TRUE)),
      suppressWarnings(range(y, finite = TRUE))
    )
  }
  expect_warning(simd_max(c(Inf, NaN), finite = TRUE), "no non-missing arguments to max")
  expect_warning(simd_min_abs(-Inf, finite = TRUE), "no non-missing arguments to min")
  # Integer and logical vectors hold no infinities: finite is na.rm.
  expect_identical(simd_range(c(4L, NA, -2L), finite = TRUE), c(-2L, 4L))
  expect_identical(simd_max(c(TRUE, NA), finite = TRUE), 1L)
  expect_identical(simd_min_abs(c(-5L, NA, 3L), finite = TRUE), 3L)
  expect_identical(simd_range(1:10, finite = TRUE), c(1L, 10L))
  # Complex magnitudes: a modulus that is infinite (or missing) goes.
  z <- complex(real = c(3, Inf, NA, 0, -Inf), imaginary = c(4, 0, 1, 1, NaN))
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_identical(simd_max_abs(z, finite = TRUE), 5)
      expect_identical(simd_min_abs(z, finite = TRUE), 1)
    })
  }
  expect_error(simd_range(1, finite = NA), "'finite' must be TRUE or FALSE", fixed = TRUE)
  expect_error(simd_max(1, finite = "yes"), "'finite' must be TRUE or FALSE", fixed = TRUE)
  expect_error(simd_range(1, 2), "na.rm, finite and na_check must be named", fixed = TRUE)
  skip_if_not_installed("bit64")
  x <- bit64::as.integer64(c(7, NA, -3))
  expect_identical(simd_range(x, finite = TRUE), bit64::as.integer64(c(-3, 7)))
  expect_identical(simd_max_abs(x, finite = TRUE), bit64::as.integer64(7))
})

test_that("reductions without na.rm stop at an NA and still give NA on every tier", {
  n <- 3e4
  x <- runif(n)
  xi <- sample.int(1000L, n, TRUE)
  for (pos in c(1, 4096, 4097, n)) {
    y <- x
    y[pos] <- NA
    yi <- xi
    yi[pos] <- NA
    for (f in list(simd_sum, simd_mean, simd_min, simd_max, simd_sum_sq, simd_prod)) {
      expect_tiers_give(NA_real_, f, y)
    }
    expect_tiers_give(NA_integer_, simd_sum, yi)
    expect_tiers_give(NA_integer_, simd_min, yi)
    expect_tiers_give(NA_real_, simd_mean, yi)
    expect_tiers_give(NA_real_, simd_prod, yi)
    for (f in list(simd_prod_sums, simd_prod_diffs, simd_dot, simd_dist, simd_cosine)) {
      expect_tiers_give(NA_real_, f, y, x)
      expect_tiers_give(NA_real_, f, x, y)
      expect_tiers_give(NA_real_, f, yi, x)
    }
    expect_tiers_give(NA_real_, simd_prod_sums, y, 0.5)
    expect_tiers_give(NA_real_, simd_prod_diffs, 0.5, yi)
    # NA wins over a NaN before it.
    y[1] <- NaN
    if (pos > 1) {
      for (f in list(simd_sum, simd_prod)) expect_tiers_give(NA_real_, f, y)
      for (f in list(simd_prod_sums, simd_dot, simd_dist, simd_cosine)) {
        expect_tiers_give(NA_real_, f, y, x)
      }
    }
  }
  # A broadcast NA, and a NaN alone, read to the end.
  expect_tiers_give(NA_real_, simd_prod_sums, x, NA_real_)
  expect_tiers_give(NA_real_, simd_prod_diffs, NA_integer_, x)
  y <- x
  y[5000] <- NaN
  for (f in list(simd_prod, simd_sum)) expect_tiers_give(NaN, f, y)
  for (f in list(simd_prod_sums, simd_dot, simd_dist)) expect_tiers_give(NaN, f, y, x)
})

test_that("NaN in rescanned blocks still lets a later NA decide, and alone gives NaN", {
  # Blocks with a NaN are rescanned for NA as they end; the final scan
  # covers only the vectors after the last block, so the NA positions are
  # in the tail, the last block, at a block edge and in a later block.
  n <- 3e4 + 3
  x <- runif(n)
  nan_at <- list(every_64 = seq(1, n - 40, by = 64), first = 1, block_2 = 4100)
  old <- options(rsimd.precision = "fast")
  on.exit(options(old))
  for (prec in c("fast", "compensated", "pairwise")) {
    options(rsimd.precision = prec)
    for (where in names(nan_at)) {
      y <- x
      y[nan_at[[where]]] <- NaN
      for (f in list(simd_sum, simd_mean, simd_prod, simd_min, simd_max, simd_sum_sq)) {
        expect_tiers_give(NaN, f, y)
      }
      for (f in list(simd_prod_sums, simd_prod_diffs, simd_dot, simd_dist)) {
        expect_tiers_give(NaN, f, y, x)
        expect_tiers_give(NaN, f, x, y)
      }
      expect_tiers_give(NaN, simd_prod_sums, y, 0.5)
      for (pos in c(n, n - 2, n - 100, 8193, 20000)) {
        z <- y
        z[pos] <- NA
        for (f in list(simd_sum, simd_mean, simd_prod, simd_min, simd_max, simd_sum_sq)) {
          expect_tiers_give(NA_real_, f, z)
        }
        for (f in list(simd_prod_sums, simd_prod_diffs, simd_dot, simd_dist)) {
          expect_tiers_give(NA_real_, f, z, x)
        }
        # The NaN in x, the NA in y.
        w <- x
        w[pos] <- NA
        for (f in list(simd_prod_sums, simd_dot, simd_dist)) {
          expect_tiers_give(NA_real_, f, y, w)
        }
        expect_tiers_give(NA_real_, simd_prod_diffs, z, 0.5)
      }
    }
  }
  # With na.rm the NaN are still counted out. Whole numbers make the sums
  # exact, so in fast mode (no refinement) the mean is the exact sum over
  # the count rounded once; base R's mean() is that only with long double.
  options(rsimd.precision = "fast")
  y <- with_seed(n, as.double(sample.int(1000L, n, TRUE)))
  y[nan_at$every_64] <- NaN
  expect_tiers_give(
    sum(y, na.rm = TRUE) / sum(!is.na(y)), function(v) simd_mean(v, na.rm = TRUE), y
  )
  expect_tiers_give(min(y, na.rm = TRUE), function(v) simd_min(v, na.rm = TRUE), y)
  expect_tiers_give(which.max(y), simd_which_max, y)
})

test_that("an NA past the first chunk, or before it, gives NA on every tier", {
  skip_unless_extended()
  n <- 2^21 + 37
  x <- runif(n, 0.9999, 1.0001)
  xi <- sample(c(-1L, 1L), n, TRUE)
  for (pos in c(1, 2^20 + 1, n)) {
    y <- x
    y[pos] <- NA
    yi <- xi
    yi[pos] <- NA
    expect_tiers_give(NA_real_, simd_prod, y)
    expect_tiers_give(NA_real_, simd_prod, yi)
    for (f in list(simd_prod_sums, simd_prod_diffs, simd_dot, simd_dist, simd_cosine)) {
      expect_tiers_give(NA_real_, f, y, x)
      expect_tiers_give(NA_real_, f, x, yi)
    }
  }
})

test_that("mean, var and sd with na.rm give the same bits as without when nothing is removed", {
  # The passes after the first skip the NA handling when the first one
  # removed nothing, and keep it when it removed something.
  fs <- list(list(mean, simd_mean), list(stats::var, simd_var), list(stats::sd, simd_sd))
  for (precision in c("fast", "pairwise", "compensated")) {
    old <- simd_precision(precision)
    on.exit(simd_precision(old))
    for (n in c(7L, 101L, 3e4L + 3L)) {
      x <- with_seed(n, 1e6 + runif(n))
      xi <- with_seed(n, sample.int(1e6, n, TRUE))
      for (v in list(x, xi)) {
        for (f in fs) {
          plain <- with_each_tier(function() f[[2]](v))
          expect_identical(with_each_tier(function() f[[2]](v, na.rm = TRUE)), plain)
          y <- v
          y[c(1L, n %/% 2L, n)] <- NA
          if (is.double(v)) y[2L] <- NaN
          # Removing elements moves the others to other lanes, so the
          # result can differ from that of the kept elements in the last bit.
          want <- f[[1]](y, na.rm = TRUE)
          for (got in with_each_tier(function() f[[2]](y, na.rm = TRUE))) {
            expect_equal(got, want, tolerance = 1e-12)
          }
        }
      }
    }
    simd_precision(old)
  }
})
