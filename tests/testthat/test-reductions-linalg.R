# sum_sq, sum_abs, norm, dot, dist, cosine, var and sd: every tier against
# the none oracle and against base R.

modes <- c("fast", "pairwise", "compensated")

# Base R equivalents. The pair functions drop a pair when either element is
# missing, under na.rm.
base_sum_sq <- function(x, na.rm = FALSE) sum(as.double(x)^2, na.rm = na.rm)
base_sum_abs <- function(x, na.rm = FALSE) sum(abs(x), na.rm = na.rm)
base_norm <- function(x, na.rm = FALSE) sqrt(base_sum_sq(x, na.rm = na.rm))
pair_ok <- function(x, y, na.rm) if (na.rm) !is.na(x) & !is.na(y) else rep(TRUE, length(x))
base_dot <- function(x, y, na.rm = FALSE) {
  ok <- pair_ok(x, y, na.rm)
  sum(as.double(x[ok]) * as.double(y[ok]))
}
base_dist <- function(x, y, na.rm = FALSE) {
  ok <- pair_ok(x, y, na.rm)
  sqrt(sum((as.double(x[ok]) - as.double(y[ok]))^2))
}
base_cosine <- function(x, y) {
  x <- as.double(x)
  y <- as.double(y)
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}
base_var <- function(x, na.rm = FALSE) stats::var(as.double(x), na.rm = na.rm)
base_sd <- function(x, na.rm = FALSE) stats::sd(as.double(x), na.rm = na.rm)

# Expects f(x, ...) on every tier within `rel` relative of `mass` (an
# absolute bound rel * mass) of the none tier's result, with missing and
# infinite values identical.
expect_fused_equal <- function(f, x, ..., mass, rel, precision = "fast") {
  old <- simd_precision(precision)
  on.exit(simd_precision(old))
  res <- with_each_tier(function() f(x, ...))
  for (tier in setdiff(names(res), "none")) {
    problem <- double_mismatch(as.double(res[[tier]]), as.double(res[["none"]]), rel * mass, rel)
    simd_expect(is.null(problem), sprintf(
      "tier %s differs from none (%s mode, n = %d): %s", tier, precision, length(x), problem
    ))
    check_identical(typeof(res[[tier]]), typeof(res[["none"]]))
  }
  invisible(res)
}

# The finite part of a double vector.
finite_mass <- function(v) sum(abs(v[is.finite(v)]))

test_that("fused reductions match the oracle on every tier for all edge lengths", {
  batch_expectations({
    for (n in sweep_lengths()) {
      seed <- as.integer(n %% 1000) + 7L
      d <- rand_vec("double", n, seed = seed)
      d2 <- rand_vec("double", n, seed = seed + 1L)
      clean <- rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = seed)
      i <- rand_vec("integer", n, seed = seed)
      i2 <- rand_vec("integer", n, seed = seed + 2L)
      l <- rand_vec("logical", n, seed = seed)
      rel <- 4 * max(n, 1) * 2^-52 + 2^-50
      for (x in sweep_inputs(n, list(d, clean, i, l))) {
        xd <- as.double(x)
        for (na_rm in c(FALSE, TRUE)) {
          for (mode in modes) {
            expect_fused_equal(simd_sum_sq, x,
              na.rm = na_rm,
              mass = finite_mass(xd^2), rel = rel, precision = mode
            )
            expect_fused_equal(simd_norm, x,
              na.rm = na_rm,
              mass = sqrt(finite_mass(xd^2)), rel = rel, precision = mode
            )
            if (is.double(x)) {
              expect_fused_equal(simd_sum_abs, x,
                na.rm = na_rm,
                mass = finite_mass(xd), rel = rel, precision = mode
              )
            }
            expect_fused_equal(simd_var, x,
              na.rm = na_rm,
              mass = 0, rel = 8 * rel, precision = mode
            )
            expect_fused_equal(simd_sd, x,
              na.rm = na_rm,
              mass = 0, rel = 8 * rel, precision = mode
            )
          }
          if (!is.double(x)) expect_simd_identical(simd_sum_abs, x, na.rm = na_rm)
        }
      }
      pairs <- list(
        list(d, d2), list(clean, d), list(i, i2), list(d, i), list(i, d), list(l, d),
        list(l, i)
      )
      for (p in sweep_inputs(n, pairs)) {
        x <- p[[1]]
        y <- p[[2]]
        xy <- as.double(x) * as.double(y)
        sq <- (as.double(x) - as.double(y))^2
        for (mode in modes) {
          for (na_rm in c(FALSE, TRUE)) {
            expect_fused_equal(simd_dot, x, y,
              na.rm = na_rm,
              mass = finite_mass(xy), rel = rel, precision = mode
            )
            expect_fused_equal(simd_dist, x, y,
              na.rm = na_rm,
              mass = sqrt(finite_mass(sq)), rel = rel, precision = mode
            )
          }
          expect_fused_equal(simd_cosine, x, y, mass = 0, rel = 8 * rel, precision = mode)
        }
      }
    }
  })
})

test_that("fused reductions agree with base R within 1e-12 on random data", {
  batch_expectations({
    for (n in c(1, 2, 3, 5, 17, 100, 1001, 4097, 1e5)) {
      seed <- as.integer(n %% 1000) + 11L
      # Only NA as missing value, where base R is deterministic.
      d <- rand_vec("double", n, nan_frac = 0, inf_frac = 0, seed = seed)
      d2 <- rand_vec("double", n, nan_frac = 0, inf_frac = 0, seed = seed + 1L)
      pos <- abs(rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = seed))
      i <- rand_vec("integer", n, seed = seed)
      l <- rand_vec("logical", n, seed = seed)
      for (mode in modes) {
        for (x in list(d, pos, i, l)) {
          for (na_rm in c(FALSE, TRUE)) {
            expect_simd_matches_base(simd_sum_sq, base_sum_sq, x,
              na.rm = na_rm, tolerance = 1e-12, precision = mode
            )
            expect_simd_matches_base(simd_sum_abs, base_sum_abs, x,
              na.rm = na_rm, tolerance = 1e-12, precision = mode
            )
            expect_simd_matches_base(simd_norm, base_norm, x,
              na.rm = na_rm, tolerance = 1e-12, precision = mode
            )
            expect_simd_matches_base(simd_var, base_var, x,
              na.rm = na_rm, tolerance = 1e-12, precision = mode
            )
            expect_simd_matches_base(simd_sd, base_sd, x,
              na.rm = na_rm, tolerance = 1e-12, precision = mode
            )
          }
        }
        for (p in list(list(pos, abs(d2)), list(i, pos), list(l, i))) {
          for (na_rm in c(FALSE, TRUE)) {
            expect_simd_matches_base(simd_dot, base_dot, p[[1]], p[[2]],
              na.rm = na_rm, tolerance = 1e-12, precision = mode
            )
          }
        }
        # dist and cosine are well conditioned on any data.
        for (p in list(list(d, d2), list(i, d2), list(l, i))) {
          for (na_rm in c(FALSE, TRUE)) {
            expect_simd_matches_base(simd_dist, base_dist, p[[1]], p[[2]],
              na.rm = na_rm, tolerance = 1e-12, precision = mode
            )
          }
          ok <- !is.na(p[[1]]) & !is.na(p[[2]])
          expect_simd_matches_base(simd_cosine, base_cosine, p[[1]][ok], p[[2]][ok],
            tolerance = 1e-12, precision = mode
          )
        }
      }
    }
  })
})

test_that("products of small integers are exact on every tier, fused or not", {
  batch_expectations({
    for (n in c(3, 31, 1000, 4097)) {
      x <- as.double(with_seed(n, sample(-1000:1000, n, replace = TRUE)))
      y <- as.double(with_seed(n + 1, sample(-1000:1000, n, replace = TRUE)))
      for (mode in modes) {
        old <- simd_precision(mode)
        expect_simd_identical(simd_sum_sq, x)
        expect_simd_identical(simd_dot, x, y)
        expect_simd_identical(simd_dist, x, y)
        expect_simd_identical(simd_sum_abs, x)
        check_identical(simd_dot(x, y), sum(x * y))
        check_identical(simd_sum_sq(as.integer(x)), sum(x^2))
        simd_precision(old)
      }
    }
  })
})

test_that("missing values: NA beats NaN, var and sd give NA, na.rm drops pairs", {
  batch_expectations({
    cases <- list(c(1, NA, NaN), c(NaN, 1, NA), c(NaN, 2), c(NA, 2), c(1:40, NaN, NA))
    for (x in cases) {
      want_na <- anyNA(x) && any(is.na(x) & !is.nan(x))
      res <- with_each_tier(function() {
        c(
          simd_sum_sq(x), simd_sum_abs(x), simd_norm(x), simd_dot(x, x), simd_dist(x, 0 * x),
          simd_cosine(x, x + 1)
        )
      })
      for (tier in names(res)) {
        check_true(all(is.na(res[[tier]])), info = tier)
        check_identical(any(is.nan(res[[tier]])), !want_na, info = tier)
      }
      res <- with_each_tier(function() c(simd_var(x), simd_sd(x)))
      for (tier in names(res)) {
        check_identical(res[[tier]], c(NA_real_, NA_real_), info = tier)
      }
    }
    # NA in x or in y removes the pair.
    x <- c(1, NA, 3, 4, NaN)
    y <- c(2, 5, NA, 1, 1)
    for (tier in tiers_to_test()) {
      simd_with_impl(tier, {
        check_identical(simd_dot(x, y, na.rm = TRUE), 6)
        check_identical(simd_dist(x, y, na.rm = TRUE), sqrt(10))
        check_identical(simd_sum_sq(x, na.rm = TRUE), 26)
        check_identical(simd_var(c(1, 2, NaN, NA), na.rm = TRUE), 0.5)
        check_identical(simd_dot(c(1L, NA), c(2L, 3L), na.rm = TRUE), 2)
        check_identical(simd_sum_abs(c(-2L, NA), na.rm = TRUE), 2L)
        check_identical(simd_sum_abs(c(-2L, NA)), NA_integer_)
        check_identical(simd_dot(c(NA, NA), c(1, 2), na.rm = TRUE), 0)
      })
    }
  })
})

test_that("na.rm drops pairs whose term is NaN (Inf * 0, Inf - Inf), as base R does", {
  for (mode in c("fast", "pairwise", "compensated")) {
    old <- simd_precision(mode)
    for (n in c(2L, 7L, 33L, 5000L)) {
      set.seed(n)
      x <- round(rnorm(n) * 8)
      y <- round(rnorm(n) * 8)
      x[c(1L, n)] <- c(Inf, -Inf)
      y[c(1L, n)] <- c(0, -Inf)
      if (n > 2L) x[2L] <- NA
      if (n > 7L) y[n - 3L] <- NaN
      want_dot <- sum(x * y, na.rm = TRUE)
      want_dist <- sqrt(sum((x - y)^2, na.rm = TRUE))
      res <- with_each_tier(function() {
        c(
          simd_dot(x, y, na.rm = TRUE), simd_dist(x, y, na.rm = TRUE),
          simd_dot(simd_vec(x), simd_vec(y), na.rm = TRUE),
          simd_dot(x, y, na.rm = TRUE, na_check = FALSE),
          simd_dist(x, y, na.rm = TRUE, na_check = FALSE)
        )
      })
      for (tier in names(res)) {
        check_identical(res[[tier]], c(want_dot, want_dist, want_dot, want_dot, want_dist),
                        info = sprintf("%s %s n = %d", tier, mode, n))
      }
    }
    simd_precision(old)
  }
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      check_identical(simd_dot(c(Inf, 2), c(0L, 3L), na.rm = TRUE), 6)
      check_identical(simd_dot(c(0L, 2L), c(Inf, 3), na.rm = TRUE), 6)
      check_identical(simd_dist(c(Inf, 1), c(Inf, 4), na.rm = TRUE), 3)
      # NA-free input is still scanned for NaN terms under na.rm.
      check_identical(simd_dot(simd_vec(c(Inf, 1)), simd_vec(c(0, 1)), na.rm = TRUE), 1)
      # Without na.rm the term stays: NaN, or NA when an element is NA.
      check_identical(simd_dot(c(Inf, 1), c(0, 1)), NaN)
      check_identical(simd_dist(c(Inf, 1), c(Inf, 1)), NaN)
      check_identical(simd_dot(c(Inf, NA, 1), c(0, 1, 1)), NA_real_)
    })
  }
})

test_that("edge values: empty, zero vectors, Inf, overflow and short input", {
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_identical(simd_sum_sq(numeric(0)), 0)
      expect_identical(simd_sum_sq(integer(0)), 0)
      expect_identical(simd_sum_abs(integer(0)), 0L)
      expect_identical(simd_sum_abs(logical(0)), 0L)
      expect_identical(simd_sum_abs(numeric(0)), 0)
      expect_identical(simd_norm(numeric(0)), 0)
      expect_identical(simd_dot(numeric(0), integer(0)), 0)
      expect_identical(simd_dist(numeric(0), numeric(0)), 0)
      expect_identical(simd_cosine(numeric(0), numeric(0)), NaN)
      expect_identical(simd_cosine(c(0, 0), c(1, 2)), NaN)
      expect_identical(simd_cosine(c(Inf, 1), c(1, 1)), NaN)
      expect_identical(simd_cosine(c(Inf, 1), c(0, 1)), NaN)
      expect_equal(simd_cosine(c(1, 2), c(2, 4)), 1, tolerance = 4 * 2^-52)
      expect_identical(simd_cosine(1:2, c(-2, 1)), 0)
      expect_identical(simd_norm(c(3L, 4L)), 5)
      expect_identical(simd_norm(c(-3, 4)), 5)
      # Naive sums of squares overflow for huge values.
      expect_identical(simd_norm(c(1e200, 1e200)), Inf)
      expect_identical(simd_sum_sq(c(-Inf, 1)), Inf)
      expect_identical(simd_var(1), NA_real_)
      expect_identical(simd_var(numeric(0)), NA_real_)
      expect_identical(simd_var(c(1, NA), na.rm = TRUE), NA_real_)
      expect_identical(simd_var(TRUE), NA_real_)
      expect_identical(simd_var(c(TRUE, FALSE)), 0.5)
      expect_identical(simd_var(c(1, Inf)), NaN)
      expect_identical(simd_sd(c(5, 5, 5)), 0)
      # The integer sum of |x| becomes double when it overflows, as sum().
      big <- rep(-.Machine$integer.max, 3L)
      expect_identical(simd_sum_abs(big), sum(abs(big)))
      expect_type(simd_sum_abs(big), "double")
      expect_identical(simd_sum_abs(c(-.Machine$integer.max, 0L)), .Machine$integer.max)
    })
  }
})

test_that("mean, var and sd rescale a sum of finite data that overflows, on every tier and length", {
  old <- simd_precision()
  on.exit(simd_precision(old))
  M <- .Machine$double.xmax
  for (mode in c("fast", "pairwise", "compensated")) {
    simd_precision(mode)
    batch_expectations({
      # The double sum overflows, so the mean is recomputed from x scaled by a
      # power of 2, as base R's long double mean gives it.
      for (n in 1:65) {
        x <- rep(1e308, n)
        expect_tiers_give(1e308, simd_mean, x)
        expect_tiers_give(-1e308, simd_mean, -x)
        expect_tiers_give(1e308, simd_mean, c(x, NA, NaN), na.rm = TRUE)
        expect_tiers_give(M, simd_mean, rep(M, n))
        if (n >= 2L) {
          want <- c(0, 0)
          expect_tiers_give(want, function(x) c(simd_var(x), simd_sd(x)), x)
          expect_tiers_give(0, simd_var, -x)
          expect_tiers_give(0, simd_var, c(x, NA), na.rm = TRUE)
        }
      }
      expect_tiers_give(NA_real_, simd_mean, c(1e308, 1e308, NA))
      expect_tiers_give(1.5 * 2^1023, simd_mean, 2^1023 * c(1.5, 1.75, 1.25, 1.5))
      # Partial sums that overflow both ways (in every lane, on every tier)
      # cancel once scaled.
      expect_tiers_give(2^1000 / 129, simd_mean, c(rep(2^1023, 64), rep(-2^1023, 64), 2^1000))
      # A deviation from the mean can still overflow: (M + M) - mean.
      expect_tiers_give(Inf, simd_var, c(M, -M, M, -M))
      # An infinite element decides the mean, whatever the finite ones give;
      # its deviation is Inf - Inf, so var and sd are NaN.
      expect_tiers_give(Inf, simd_mean, c(Inf, -1e308, -1e308, -1e308))
      expect_tiers_give(-Inf, simd_mean, c(1e308, 1e308, -Inf, 1e308))
      expect_tiers_give(NaN, simd_mean, c(Inf, 1e308, 1e308, -Inf))
      for (n in 2:17) {
        expect_tiers_give(NaN, simd_var, c(rep(1e308, n), Inf))
        expect_tiers_give(NaN, simd_sd, c(Inf, rep(1, n), -Inf))
      }
    })
  }
})

test_that("the base R table rows for var and sd are reproduced", {
  cases <- list(
    quote(var(1)), quote(var(numeric(0))), quote(var(c(1, NA))),
    quote(var(c(1, NA), na.rm = TRUE)), quote(sd(c(1, 2, 3, 4))), quote(var(1:3)),
    quote(var(c(1, NaN))), quote(var(c(NaN, NA))), quote(sd(c(TRUE, FALSE, TRUE)))
  )
  for (e in cases) {
    want <- eval(e, asNamespace("stats"))
    simd_e <- e
    simd_e[[1L]] <- as.name(paste0("simd_", as.character(e[[1L]])))
    got <- with_each_tier(function() eval(simd_e))
    for (tier in names(got)) {
      # Base R's var works in long double, so finite values may differ in the
      # last bit.
      if (is.finite(want)) {
        expect_equal(got[[tier]], want, tolerance = 4 * 2^-52, info = paste(tier, deparse1(e)))
      } else {
        expect_identical(got[[tier]], want, info = paste(tier, deparse1(e)))
      }
    }
  }
})

test_that("var is two-pass: no cancellation for a large mean", {
  x <- 1e9 + with_seed(5, stats::rnorm(1e5))
  for (mode in modes) {
    expect_simd_matches_base(simd_var, stats::var, x, tolerance = 1e-9, precision = mode)
  }
  xi <- 1e8L + with_seed(6, sample.int(100L, 1e5, replace = TRUE))
  expect_simd_matches_base(simd_var, base_var, xi, tolerance = 1e-12)
})

test_that("var refines the mean in every mode, as base R does", {
  # The fast-mode sum of 1e15 + {0, 1} gives a mean 0.5 off, which without
  # base R's refinement pass made the variance twice its true value.
  for (mean in c(1e15, 1e12)) {
    x <- mean + with_seed(2, sample(c(0, 1), 1e5, replace = TRUE))
    xna <- replace(x, c(1, 500, 7777), c(NA, NaN, NA))
    for (mode in modes) {
      for (f in list(c(simd_var, stats::var), c(simd_sd, stats::sd))) {
        expect_simd_matches_base(f[[1]], f[[2]], x, tolerance = 1e-12, precision = mode)
        expect_simd_matches_base(f[[1]], f[[2]], xna, na.rm = TRUE,
                                 tolerance = 1e-12, precision = mode)
      }
    }
  }
})

test_that("na_check = FALSE gives the same results on NA-free input", {
  d <- rand_vec("double", 5000, na_frac = 0, nan_frac = 0, seed = 3)
  i <- rand_vec("integer", 5000, na_frac = 0, seed = 3)
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      for (x in list(d, i)) {
        for (f in list(simd_sum_sq, simd_sum_abs, simd_norm, simd_var, simd_sd)) {
          expect_identical(f(x, na_check = FALSE), f(x), info = tier)
        }
        for (f in list(simd_dot, simd_dist, simd_cosine)) {
          expect_identical(f(x, d, na_check = FALSE), f(x, d), info = tier)
        }
      }
      expect_type(simd_dot(c(1, NA), c(1, 2), na_check = FALSE), "double")
    })
  }
})

test_that("compact sequences give the same results as expanded vectors", {
  batch_expectations({
    for (n in c(2, 4097, 1e5)) {
      for (x in altrep_inputs(n)) {
        check_true(takes_region_path(x))
        expanded <- x + if (is.integer(x)) 0L else 0
        for (f in list(simd_sum_sq, simd_sum_abs, simd_norm, simd_var, simd_sd)) {
          check_identical(f(x), f(expanded))
        }
        check_identical(simd_dot(x, rev(x)), simd_dot(expanded, rev(x)))
        check_identical(simd_dist(x, expanded), 0)
        check_identical(simd_cosine(x, expanded), simd_cosine(expanded, expanded))
      }
    }
  })
})

test_that("var and sd span chunks: forced small chunks match base R", {
  skip_on_cran()
  skip_if_no_subprocess()
  n <- 2^20 + 7
  d <- 1e6 + rand_vec("double", n, na_frac = 0.001, nan_frac = 0, inf_frac = 0, seed = 4)
  i <- rand_vec("integer", n, na_frac = 0.001, seed = 4)
  y <- rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = 5)
  child <- function(d, i, y) {
    library(rsimd)
    lapply(stats::setNames(nm = simd_available()), function(tier) {
      simd_with_impl(tier, lapply(c("fast", "pairwise", "compensated"), function(m) {
        old <- simd_precision(m)
        on.exit(simd_precision(old))
        c(
          simd_var(d, na.rm = TRUE), simd_sd(i, na.rm = TRUE), simd_var(i),
          simd_dot(d, y, na.rm = TRUE), simd_dist(i, y, na.rm = TRUE), simd_cosine(y, y + 1),
          simd_sum_sq(i, na.rm = TRUE), simd_sum_abs(d, na.rm = TRUE)
        )
      }))
    })
  }
  ok <- !is.na(d)
  want <- c(
    stats::var(d, na.rm = TRUE), stats::sd(i, na.rm = TRUE), NA, sum(d[ok] * y[ok]),
    base_dist(i, y, na.rm = TRUE), base_cosine(y, y + 1), base_sum_sq(i, na.rm = TRUE),
    sum(abs(d), na.rm = TRUE)
  )
  env <- c(callr::rcmd_safe_env(), RSIMD_DEBUG_STRIDE = "128")
  there <- callr::r(child, list(d, i, y), env = env)
  here <- child(d, i, y)
  # Without a long double accumulator, base R's sums of 2^20 terms are only
  # within about n * eps.
  tol <- if (has_wide_long_double()) 1e-12 else n * eps
  for (tier in names(there)) {
    for (k in 1:3) {
      expect_equal(there[[tier]][[k]], want, tolerance = tol, info = paste(tier, k))
    }
    # Pairwise sums are independent of the chunking.
    if (tier %in% names(here)) expect_identical(there[[tier]][[2]], here[[tier]][[2]])
  }
})

test_that("argument errors", {
  expect_error(simd_dot(1:3, 1:2), "lengths of 'x' \\(3\\) and 'y' \\(2\\) must be equal")
  expect_error(simd_dot(1:3, 1), "must be equal")
  expect_error(simd_dist(1, 1:3), "must be equal")
  expect_error(simd_cosine(1:2, 1:3), "must be equal")
  for (f in list(simd_sum_sq, simd_sum_abs, simd_norm, simd_var, simd_sd)) {
    expect_error(f(as.raw(1:3)), "invalid 'type' (raw) of argument", fixed = TRUE)
    expect_error(f(c(1i, 2i)), "does not support 'x' of type complex")
    expect_error(f(letters), "must be an atomic vector")
    expect_error(f(1:3, na.rm = NA), "'na.rm' must be TRUE or FALSE")
  }
  expect_error(simd_dot(1:3, as.raw(1:3)), "invalid 'type' (raw) of argument", fixed = TRUE)
  expect_error(simd_dot(1:2, c(1i, 2i)), "does not support 'y' of type complex")
  expect_error(simd_dist(factor(1:2), 1:2), "must be an atomic vector")
  skip_if_not_installed("bit64")
  expect_error(simd_dot(bit64::as.integer64(1:2), 1:2), "does not support 'x' of type integer64$")
  expect_error(simd_var(bit64::as.integer64(1:2)), "does not support 'x' of type integer64$")
})

test_that("results are bare doubles (or integers for integer sum_abs)", {
  x <- c(a = 1, b = -2)
  expect_identical(simd_sum_sq(x), 5)
  expect_identical(simd_sum_abs(c(a = -1L)), 1L)
  expect_null(attributes(simd_dot(x, x)))
  expect_null(attributes(simd_var(structure(1:3, dim = 3L))))
})
