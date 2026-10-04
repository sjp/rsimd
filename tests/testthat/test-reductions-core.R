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
      expect_lte(abs(got[[tier]] - want), bound + 2^-52 * abs(want),
        label = paste(tier, precision, length(x))
      )
    }
  }
}

test_that("every reduction matches the oracle on every tier for all edge lengths", {
  for (n in edge_lengths()) {
    seed <- as.integer(n %% 1000) + 1L
    d <- rand_vec("double", n, seed = seed)
    i <- rand_vec("integer", n, seed = seed)
    l <- rand_vec("logical", n, seed = seed)
    p <- near_one(n, seed)
    inputs <- list(
      d, rand_vec("double", n, na_frac = 0, nan_frac = 0, seed = seed), i,
      rand_vec("integer", n, na_frac = 0, seed = seed), l
    )
    for (x in inputs) {
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

test_that("reductions agree with base R across edge lengths", {
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

test_that("the base R behaviour table is reproduced", {
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
      expect_identical(g, want, info = paste(tier, deparse1(e)))
      expect_true(identical(g[[1L]], want[[1L]], num.eq = FALSE), info = paste(tier, deparse1(e)))
    }
  }
})

test_that("NA beats NaN in either order for prod, mean, min, max and range", {
  for (x in list(c(NA, NaN), c(NaN, NA), c(1, NaN, 2, NA), c(NA, 1:40, NaN), c(NaN, 1:40, NA))) {
    for (f in list(simd_prod, simd_mean, simd_min, simd_max, simd_range)) {
      res <- with_each_tier(function() f(x))
      for (tier in names(res)) {
        expect_true(all(is.na(res[[tier]]) & !is.nan(res[[tier]])), info = tier)
      }
    }
  }
})

test_that("a zero extremum has the sign of the first zero, wherever the zeros are", {
  for (n in c(2, 3, 5, 9, 17, 33, 65, 130, 300)) {
    for (k in 1:12) {
      pos <- with_seed(n * 100 + k, sort(sample.int(n, min(n, 4))))
      x <- with_seed(k, stats::runif(n, 1, 2))
      signs <- with_seed(k + 7, sample(c(0, -0), length(pos), replace = TRUE))
      x[pos] <- signs
      lo <- expect_simd_identical(simd_min, x)
      for (tier in names(lo)) {
        expect_true(identical(lo[[tier]], min(x), num.eq = FALSE), info = paste(tier, n, k))
      }
      y <- -x
      hi <- expect_simd_identical(simd_max, y)
      for (tier in names(hi)) {
        expect_true(identical(hi[[tier]], max(y), num.eq = FALSE), info = paste(tier, n, k))
      }
      r <- expect_simd_identical(simd_range, c(x, y))
      for (tier in names(r)) {
        expect_true(identical(r[[tier]], range(c(x, y)), num.eq = FALSE), info = tier)
      }
      expect_simd_identical(simd_which_min, x)
      expect_simd_matches_base(simd_which_min, which.min, x)
      expect_simd_matches_base(simd_which_max, which.max, y)
    }
  }
})

test_that("which_*: first of ties, NA and NaN ignored, empty gives integer(0)", {
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

test_that("compact sequences give the same results as expanded vectors", {
  for (n in c(2, 4097, 1e6)) {
    for (x in altrep_inputs(n)) {
      expect_true(takes_region_path(x))
      expanded <- x + if (is.integer(x)) 0L else 0
      expect_false(takes_region_path(expanded))
      for (f in list(
        simd_min, simd_max, simd_range, simd_which_min, simd_which_max,
        simd_any_na, simd_count_na, simd_which_na
      )) {
        a <- expect_simd_identical(f, x)
        b <- expect_simd_identical(f, expanded)
        expect_identical(a, b)
      }
      for (mode in modes) {
        a <- expect_simd_equal(simd_mean, x, precision = mode)
        b <- expect_simd_equal(simd_mean, expanded, precision = mode)
        expect_identical(a, b, info = mode)
      }
      expect_identical(simd_which_max(x), which.max(expanded))
      expect_identical(simd_mean(x), mean(expanded))
    }
  }
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
  for (tier in names(res)) {
    expect_equal(res[[tier]], mean(x), tolerance = 1e-12, info = tier)
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
  for (n in c(1:70, 127, 128, 129, 1000)) {
    ones <- as.raw(rep(c(1, 128, 255, 129), length.out = n))
    expect_simd_identical(function(x) suppressWarnings(simd_all(x)), ones)
    for (p in unique(c(1, max(n %/% 2, 1), n, seq_len(min(n, 40))))) {
      x <- ones
      x[p] <- as.raw(0)
      res <- with_each_tier(function() suppressWarnings(c(simd_all(x), simd_any(x))))
      want <- c(FALSE, n > 1)
      for (tier in names(res)) expect_identical(res[[tier]], want, info = paste(tier, n, p))
      z <- raw(n)
      z[p] <- as.raw(1)
      res <- with_each_tier(function() suppressWarnings(c(simd_all(z), simd_any(z))))
      want <- c(n == 1, TRUE)
      for (tier in names(res)) expect_identical(res[[tier]], want, info = paste(tier, n, p))
    }
    res <- with_each_tier(function() suppressWarnings(c(simd_all(ones), simd_any(raw(n)))))
    for (tier in names(res)) expect_identical(res[[tier]], c(TRUE, FALSE), info = paste(tier, n))
  }
})

test_that("missing-value queries cover doubles, integers, logicals, complex and raw", {
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
  expect_identical(simd_any_na(as.raw(1:3)), FALSE)
  expect_identical(simd_count_na(as.raw(1:3)), 0)
  expect_identical(simd_which_na(as.raw(1:3)), integer(0))
  expect_identical(simd_which_na(numeric(0)), integer(0))
  expect_identical(simd_count_na(integer(0)), 0)
  expect_identical(simd_any_na(complex(real = 1, imaginary = NaN)), TRUE)
})

test_that("chunked results equal unchunked ones", {
  skip_on_cran()
  skip_if_not_installed("callr")
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

test_that("arguments and types are validated", {
  for (f in c("simd_prod", "simd_mean", "simd_min", "simd_max", "simd_range")) {
    expect_error(get(f)(as.raw(1)), "invalid 'type' (raw) of argument", fixed = TRUE)
  }
  for (f in c("simd_min", "simd_max", "simd_range")) {
    expect_error(get(f)(1i), "invalid 'type' (complex) of argument", fixed = TRUE)
  }
  no_complex <- c(
    "simd_prod", "simd_mean", "simd_which_min", "simd_which_max", "simd_any", "simd_all"
  )
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
