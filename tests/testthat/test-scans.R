# cumsum, cumprod, cummin and cummax: every tier against the none oracle
# and against base R.

scans <- list(
  cumsum = list(simd = simd_cumsum, base = cumsum),
  cumprod = list(simd = simd_cumprod, base = cumprod),
  cummin = list(simd = simd_cummin, base = cummin),
  cummax = list(simd = simd_cummax, base = cummax)
)

# Expects a and b (double vectors) to have their missing values, NaN and
# infinities in the same places, and the other elements within `bound`
# (a vector, per element) of each other.
expect_scan_close <- function(a, b, bound, info = NULL) {
  check_identical(is.na(a), is.na(b), info = info)
  check_identical(is.nan(a), is.nan(b), info = info)
  inf <- is.infinite(a) | is.infinite(b)
  check_identical(a[inf], b[inf], info = info)
  ok <- !is.na(a) & !inf
  bad <- which(!(abs(a[ok] - b[ok]) <= bound[ok]))
  simd_expect(length(bad) == 0L, sprintf(
    "%s: element %d is %.17g, expected %.17g (bound %.3g)", paste(info, collapse = " "),
    which(ok)[bad[1]], a[ok][bad[1]], b[ok][bad[1]], bound[ok][bad[1]]
  ))
}

# Error bound of a reordered running sum of x: (i + 8) * eps times the
# running sum of |x| (finite elements; the bound stops mattering after the
# first missing or infinite value, whose position must match exactly).
cumsum_bound <- function(x) {
  a <- abs(as.double(x))
  a[!is.finite(a)] <- 0
  (seq_along(a) + 8) * 2^-52 * cumsum(a)
}

# Products of doubles near 1: relative bound 2 * (i + 8) * eps.
cumprod_bound <- function(want) 2 * (seq_along(want) + 8) * 2^-52 * abs(want)

near_one <- function(n, seed, na_frac = 0.002, nan_frac = 0.002) {
  x <- 1 + rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = seed) / 1e4
  with_seed(seed + 1L, {
    x[sample.int(n, stats::rbinom(1, n, nan_frac))] <- NaN
    x[sample.int(n, stats::rbinom(1, n, na_frac))] <- NA
  })
  x
}

test_that("the base R behaviour table is reproduced", {
  cases <- list(
    quote(cumsum(c(.Machine$integer.max, 1L, 2L))), quote(cumsum(c(1L, NA, 3L))),
    quote(cumsum(c(1, NA, 3))), quote(cumsum(c(1, NaN, NA))), quote(cumsum(c(1, NA, NaN))),
    quote(cummax(c(1, NaN, 3))), quote(cummax(c(1, NA, 3))), quote(cummin(c(3L, NA, 1L))),
    quote(cumprod(1:3)), quote(cumsum(c(TRUE, TRUE))), quote(cumsum(numeric(0))),
    quote(cumsum(integer(0))), quote(cumprod(logical(0))), quote(cummax(integer(0))),
    quote(cummax(c(1, NaN, NA))), quote(cummin(c(NaN, NA))), quote(cumprod(c(1, NA, NaN))),
    quote(cumprod(c(NA, 2L))), quote(cumsum(c(NA, .Machine$integer.max, 1L))),
    quote(cumsum(c(-.Machine$integer.max, -1L))), quote(cummin(c(TRUE, FALSE))),
    quote(cummax(c(TRUE, NA))), quote(cumsum(c(1, Inf, -Inf, 2))), quote(cumsum(c(Inf, NA))),
    quote(cumsum(c(Inf, -Inf, NA))), quote(cumprod(c(1e308, 10, NA))),
    quote(cummax(c(0, -0))), quote(cummax(c(-0, 0))), quote(cummin(c(0, -0))),
    quote(cummin(c(-0, 0))), quote(cumsum(c(-0, -0))), quote(cumprod(c(-1, 0))),
    quote(cumsum(c(1.5, 2.5, -4))), quote(cumprod(c(2L, -3L, NA, 4L)))
  )
  to_simd <- function(e) {
    e[[1L]] <- as.name(paste0("simd_", as.character(e[[1L]])))
    e
  }
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
      expect_true(identical(got[[tier]], want, num.eq = FALSE), info = paste(tier, deparse1(e)))
    }
  }
})

test_that("integer scans and cummin/cummax are bit-identical to none and base R", {
  batch_expectations({
    for (n in sweep_lengths()) {
      seed <- as.integer(n %% 1000) + 3L
      d <- rand_vec("double", n, na_frac = 0.002, nan_frac = 0.002, seed = seed)
      inputs <- list(
        rand_vec("integer", n, na_frac = 0, seed = seed) %/% 1000L,
        rand_vec("integer", n, na_frac = 0.002, seed = seed) %/% 1000L,
        rand_vec("logical", n, na_frac = 0, seed = seed),
        rand_vec("logical", n, na_frac = 0.002, seed = seed),
        rand_vec("double", n, na_frac = 0, nan_frac = 0, seed = seed), d,
        round(rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = seed))
      )
      for (x in sweep_inputs(n, inputs)) {
        for (f in c("cummin", "cummax", if (!is.double(x)) "cumsum")) {
          res <- expect_simd_identical(scans[[f]]$simd, x)
          want <- scans[[f]]$base(x)
          for (tier in names(res)) {
            check_identical(res[[tier]], want, info = paste(tier, f, n), num.eq = FALSE)
          }
        }
      }
    }
  })
})

test_that("double cumsum and cumprod are within tolerance of none and base R", {
  batch_expectations({
    for (n in sweep_lengths()) {
      seed <- as.integer(n %% 1000) + 5L
      for (x in sweep_inputs(n, list(
        rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = seed),
        rand_vec("double", n, na_frac = 0.001, nan_frac = 0.001, inf_frac = 0.001, seed = seed)
      ))) {
        res <- with_each_tier(function() simd_cumsum(x))
        for (tier in names(res)) {
          expect_scan_close(res[[tier]], res[["none"]], cumsum_bound(x), info = c(tier, n))
          expect_scan_close(res[[tier]], cumsum(x), cumsum_bound(x), info = c(tier, n, "base"))
        }
      }
      p <- near_one(n, seed)
      want <- cumprod(p)
      res <- with_each_tier(function() simd_cumprod(p))
      for (tier in names(res)) {
        expect_scan_close(res[[tier]], res[["none"]], cumprod_bound(want), info = c(tier, n))
        expect_scan_close(res[[tier]], want, cumprod_bound(want), info = c(tier, n, "base"))
      }
    }
  })
})

test_that("double cumsum is within 1e-12 relative of base R on random data", {
  for (n in c(10, 1000, 1e5)) {
    x <- abs(rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = 8))
    for (mode in c("fast", "pairwise", "compensated")) {
      expect_simd_matches_base(simd_cumsum, cumsum, x, tolerance = 1e-12, precision = mode)
    }
  }
})

test_that("compensated cumsum is the same sequential sum on every tier", {
  old <- simd_precision("compensated")
  on.exit(simd_precision(old))
  for (n in c(1, 7, 100, 4097)) {
    x <- rand_vec("double", n, na_frac = 0.01, nan_frac = 0.01, seed = n)
    expect_simd_identical(simd_cumsum, x)
  }
  # Compensation recovers what a plain running sum loses.
  x <- c(1e16, 1, 1, -1e16, 1)
  res <- with_each_tier(function() simd_cumsum(x))
  for (tier in names(res)) expect_identical(res[[tier]], c(1e16, 1e16 + 1, 1e16 + 2, 2, 3))
  simd_precision("fast")
  expect_identical(simd_with_impl("none", simd_cumsum(x))[4:5], c(0, 1))
})

test_that("NA and NaN are placed as in base R wherever they occur", {
  batch_expectations({
    base <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
    for (n in c(5, 12, 40)) {
      x0 <- rep_len(base, n)
      for (a in seq(1, n, by = 3)) {
        for (b in unique(c(a, a + 1, n))) {
          for (first in list(NA, NaN)) {
            x <- x0
            x[a] <- first
            x[b] <- if (is.nan(first)) NA else NaN
            for (f in names(scans)) {
              want <- scans[[f]]$base(x)
              res <- with_each_tier(function() scans[[f]]$simd(x))
              for (tier in names(res)) {
                check_identical(is.na(res[[tier]]), is.na(want), info = paste(tier, f, n, a, b))
                check_identical(is.nan(res[[tier]]), is.nan(want), info = paste(tier, f, n, a, b))
              }
            }
          }
        }
      }
    }
  })
})

test_that("integer cumsum overflows at the right element with one warning", {
  batch_expectations({
    for (n in c(3, 9, 40, 300)) {
      for (pos in unique(pmax(pmin(c(2, 3, 5, n %/% 2, n), n), 2))) {
        x <- rep(1L, n)
        x[pos - 1L] <- .Machine$integer.max - as.integer(pos) + 2L
        for (sign in c(1L, -1L)) {
          y <- sign * x
          want <- withCallingHandlers(cumsum(y),
            warning = function(w) invokeRestart("muffleWarning")
          )
          for (tier in tiers_to_test()) {
            w <- character(0)
            got <- withCallingHandlers(simd_with_impl(tier, simd_cumsum(y)), warning = function(c) {
              w <<- c(w, conditionMessage(c))
              invokeRestart("muffleWarning")
            })
            check_identical(got, want, info = paste(tier, n, pos, sign))
            check_identical(w, "integer overflow in 'cumsum'; use 'cumsum(as.numeric(.))'")
          }
        }
      }
    }
    # Partial sums that leave the int32 range inside a vector while the
    # running total stays in range are not an overflow.
    x <- c(-2000000000L, rep(c(2000000000L, -2000000000L), 50))
    expect_simd_identical(simd_cumsum, x)
    check_identical(simd_cumsum(x), cumsum(x))
    # An NA before the overflow gives no warning.
    y <- c(1L, NA, .Machine$integer.max, 5L)
    expect_silent(check_identical(simd_cumsum(y), c(1L, rep(NA, 3))))
    check_identical(suppressWarnings(simd_cumsum(1:70000))[65535:65536], c(2147450880L, NA))
  })
})

test_that("compact sequences scan like expanded vectors", {
  batch_expectations({
    for (n in c(2, 4097, 65535)) {
      for (x in altrep_inputs(n)) {
        check_true(takes_region_path(x))
        expanded <- x + if (is.integer(x)) 0L else 0
        for (f in names(scans)) {
          # Double products are reassociated differently by each tier.
          if (f == "cumprod") {
            a <- with_each_tier(function() simd_cumprod(x))
            b <- with_each_tier(function() simd_cumprod(expanded))
          } else {
            a <- expect_simd_identical(scans[[f]]$simd, x)
            b <- expect_simd_identical(scans[[f]]$simd, expanded)
          }
          check_identical(a, b, info = f)
        }
        check_identical(simd_cumsum(x), cumsum(expanded))
      }
    }
  })
})

test_that("scans carry across chunk boundaries", {
  n <- 2^20 + 7
  i <- rand_vec("integer", n, na_frac = 0, seed = 12) %/% 100000L
  d <- rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = 12)
  expect_identical(simd_cumsum(i), cumsum(i))
  expect_identical(simd_cummax(d), cummax(d))
  expect_identical(simd_cummin(i), cummin(i))
  res <- with_each_tier(function() simd_cumsum(d))
  for (tier in names(res)) expect_scan_close(res[[tier]], cumsum(d), cumsum_bound(d), info = tier)
  skip_on_cran()
  skip_if_no_subprocess()
  child <- function(i, d) {
    library(rsimd)
    lapply(stats::setNames(nm = simd_available()), function(tier) {
      simd_with_impl(tier, list(
        simd_cumsum(i), simd_cummax(d), simd_cumsum(d), simd_cumprod(1 + d / 1e9),
        simd_cumsum(c(d[1:5000], NaN, d[1:5000], NA, 1))
      ))
    })
  }
  here <- child(i, d)
  env <- c(callr::rcmd_safe_env(), RSIMD_DEBUG_STRIDE = "128")
  there <- callr::r(child, list(i, d), env = env)
  for (tier in intersect(names(here), names(there))) {
    # Chunks start at multiples of every vector width, so even the double
    # scans are unchanged.
    expect_true(identical(there[[tier]], here[[tier]], num.eq = FALSE), info = tier)
  }
})

test_that("argument errors and attributes", {
  for (f in scans) {
    expect_error(f$simd(as.raw(1:3)), "invalid 'type' (raw) of argument", fixed = TRUE)
    if (identical(f$simd, simd_cummin) || identical(f$simd, simd_cummax)) {
      expect_error(f$simd(c(1i, 2i)), "does not support 'x' of type complex$")
    }
    expect_error(f$simd(letters), "must be an atomic vector")
    expect_identical(f$simd(c(a = 1, b = 2)), unname(f$base(c(a = 1, b = 2))))
  }
  expect_null(attributes(simd_cumsum(matrix(1:4, 2))))
  skip_if_not_installed("bit64")
  expect_error(simd_cumprod(bit64::as.integer64(1:2)), "type integer64$")
})
