# simd_which, simd_count, the magnitude reductions (max_abs, min_abs,
# which_max_abs, which_min_abs) and the products of sums and differences.
# Every tier against base R: exactly, except the products, which round in
# a different order.

# Products of many random doubles over- or underflow depending on the order
# of multiplication, so products use values near 1 (plus specials in x).
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

# The first index of the maximum (minimum) of an integer64 vector, NA
# ignored: bit64 has no which.max().
which_i64 <- function(a, f) {
  if (all(is.na(a))) {
    return(integer(0))
  }
  which(a == f(a, na.rm = TRUE))[1L]
}

test_that("simd_which and simd_count are which() and sum() of a logical vector", {
  for (n in edge_lengths()) {
    seed <- as.integer(n %% 1000) + 1L
    for (x in list(
      rand_vec("logical", n, seed = seed), rand_vec("logical", n, na_frac = 0, seed = seed),
      rep_len(c(TRUE, NA), n), rep_len(FALSE, n), rep_len(TRUE, n)
    )) {
      expect_tiers_give(which(x), simd_which, x)
      expect_tiers_give(sum(x), simd_count, x)
      expect_tiers_give(sum(x, na.rm = TRUE), simd_count, x, na.rm = TRUE)
    }
  }
  expect_identical(simd_which(logical(0)), integer(0))
  expect_identical(simd_count(logical(0)), 0L)
  expect_identical(simd_count(NA), NA_integer_)
  expect_identical(simd_count(NA, na.rm = TRUE), 0L)
  expect_identical(simd_which(c(a = TRUE, b = FALSE)), 1L)
})

test_that("simd_which and simd_count take only logical vectors", {
  for (x in list(1, 1L, as.raw(1), 1i)) {
    expect_error(simd_which(x), "simd_which\\(\\) does not support 'x' of type ")
    expect_error(simd_count(x), "simd_count\\(\\) does not support 'x' of type ")
  }
  if (has_bit64()) {
    expect_error(simd_which(simd_as_integer64(1)), "of type integer64")
  }
})

test_that("simd_count skips the NA scan for an NA-free simd_vec, and simd_which reads one", {
  x <- rep_len(c(TRUE, FALSE, FALSE), 1001)
  v <- simd_vec(x)
  expect_identical(simd_count(v), sum(x))
  expect_identical(simd_which(v), which(x))
  expect_identical(simd_count(simd_vec(c(x, NA))), NA_integer_)
})

test_that("the magnitude reductions are max, min, which.max and which.min of abs(x)", {
  base_max_abs <- function(x, ...) max(abs(x), ...)
  base_min_abs <- function(x, ...) min(abs(x), ...)
  for (n in edge_lengths()) {
    seed <- as.integer(n %% 1000) + 1L
    d <- rand_vec("double", n, seed = seed)
    inputs <- list(
      d, rand_vec("double", n, na_frac = 0, nan_frac = 0, seed = seed),
      rand_vec("integer", n, seed = seed), rand_vec("integer", n, na_frac = 0, seed = seed),
      rand_vec("logical", n, seed = seed), -abs(rand_vec("double", n, na_frac = 0, seed = seed))
    )
    for (x in inputs) {
      for (na_rm in c(FALSE, TRUE)) {
        want <- suppressWarnings(list(
          max = base_max_abs(x, na.rm = na_rm), min = base_min_abs(x, na.rm = na_rm)
        ))
        # NA wins over NaN in rsimd, whatever the order (base R depends on it).
        if (!na_rm && anyNA(x) && any(is.na(x) & !is.nan(x))) {
          want <- lapply(want, function(w) w[NA_integer_])
        }
        quiet_max <- function(...) suppressWarnings(simd_max_abs(...))
        quiet_min <- function(...) suppressWarnings(simd_min_abs(...))
        expect_tiers_give(want$max, quiet_max, x, na.rm = na_rm)
        expect_tiers_give(want$min, quiet_min, x, na.rm = na_rm)
      }
      expect_tiers_give(which.max(abs(x)), simd_which_max_abs, x)
      expect_tiers_give(which.min(abs(x)), simd_which_min_abs, x)
    }
  }
})

test_that("the magnitude reductions follow base R at the edges", {
  # Ties go to the first index; -0 and 0 tie; the signed element is x[i].
  x <- c(-3, 3, -0, 0, 2, -3)
  expect_tiers_give(1L, simd_which_max_abs, x)
  expect_tiers_give(3L, simd_which_min_abs, x)
  expect_identical(x[simd_which_max_abs(x)], -3)
  # The magnitude is never negative, and a zero one is +0.
  expect_identical(1 / simd_min_abs(c(-0, 5)), Inf)
  expect_identical(simd_max_abs(c(-7L, 3L)), 7L)
  expect_identical(simd_max_abs(c(TRUE, FALSE)), 1L)
  # Missing values: NA unless na.rm, NA over NaN; which_* skips them.
  expect_identical(simd_max_abs(c(1, NaN, -2)), NaN)
  expect_identical(simd_max_abs(c(1, NaN, NA)), NA_real_)
  expect_identical(simd_max_abs(c(1, NaN, -2), na.rm = TRUE), 2)
  expect_identical(simd_max_abs(c(1L, NA)), NA_integer_)
  expect_identical(simd_which_max_abs(c(NA, -5, NaN, 4)), 2L)
  expect_identical(simd_which_max_abs(c(NA, NaN)), integer(0))
  expect_identical(simd_which_min_abs(integer(0)), integer(0))
  # Empty input, as max(abs(x)) and min(abs(x)).
  expect_warning(
    expect_identical(simd_max_abs(numeric(0)), -Inf),
    "no non-missing arguments to max; returning -Inf"
  )
  expect_warning(
    expect_identical(simd_min_abs(c(NA, NA), na.rm = TRUE), Inf),
    "no non-missing arguments to min; returning Inf"
  )
  expect_warning(expect_identical(simd_max_abs(integer(0)), -Inf), "returning -Inf")
  # |NA_integer_| stays NA (it is not 2^31), and |.Machine$integer.max| fits.
  m <- .Machine$integer.max
  expect_tiers_give(m, simd_max_abs, c(-m, NA, 5L), na.rm = TRUE)
  expect_tiers_give(1L, simd_which_max_abs, c(-m, NA, m))
  expect_tiers_give(5L, simd_min_abs, c(-m, NA, 5L), na.rm = TRUE)
  expect_error(simd_max_abs(as.raw(1)), "simd_max_abs\\(\\) does not support 'x' of type raw")
  expect_error(simd_which_min_abs(as.raw(1)), "does not support 'x' of type raw")
})

test_that("the magnitude reductions of complex numbers use the moduli of simd_abs()", {
  for (n in edge_lengths()) {
    seed <- as.integer(n %% 1000) + 1L
    z <- complex(
      real = rand_vec("double", n, seed = seed),
      imaginary = rand_vec("double", n, na_frac = 0, nan_frac = 0, seed = seed + 1L)
    )
    for (na_rm in c(FALSE, TRUE)) {
      res <- with_each_tier(function() {
        list(
          max = suppressWarnings(simd_max_abs(z, na.rm = na_rm)),
          min = suppressWarnings(simd_min_abs(z, na.rm = na_rm)),
          max_ref = suppressWarnings(simd_max(simd_abs(z), na.rm = na_rm)),
          min_ref = suppressWarnings(simd_min(simd_abs(z), na.rm = na_rm)),
          which_max = simd_which_max_abs(z), which_min = simd_which_min_abs(z),
          which_max_ref = simd_which_max(simd_abs(z)), which_min_ref = simd_which_min(simd_abs(z))
        )
      })
      for (tier in names(res)) {
        r <- res[[tier]]
        expect_identical(r$max, r$max_ref, info = paste(tier, n))
        expect_identical(r$min, r$min_ref, info = paste(tier, n))
        expect_identical(r$which_max, r$which_max_ref, info = paste(tier, n))
        expect_identical(r$which_min, r$which_min_ref, info = paste(tier, n))
      }
    }
    # Within an ULP of base R's Mod().
    if (n > 0) {
      want <- suppressWarnings(max(Mod(z), na.rm = TRUE))
      got <- suppressWarnings(simd_max_abs(z, na.rm = TRUE))
      if (is.finite(want)) expect_lte(abs(got - want), 2^-52 * want)
    }
  }
  z <- complex(real = c(3, 0, -1, NA), imaginary = c(4, -6, 1, 0))
  expect_identical(simd_max_abs(z), NA_real_)
  expect_identical(simd_max_abs(z, na.rm = TRUE), 6)
  expect_identical(simd_which_max_abs(z), 2L)
  expect_identical(simd_which_min_abs(z), 3L)
  # Mod(Inf + NaN i) is Inf, as in base R.
  expect_identical(simd_max_abs(complex(real = Inf, imaginary = NaN)), Inf)
  expect_warning(expect_identical(simd_max_abs(complex(0)), -Inf), "returning -Inf")
})

test_that("the magnitude reductions of integer64 are integer64", {
  skip_if_not(has_bit64())
  x <- i64_dec(c("-9223372036854775807", "5", NA, "-3"))
  expect_tiers_i64(i64_dec("9223372036854775807"), simd_max_abs, x, na.rm = TRUE)
  expect_tiers_i64(i64_dec("3"), simd_min_abs, x, na.rm = TRUE)
  expect_tiers_i64(i64_dec(NA), simd_max_abs, x)
  expect_tiers_give(1L, simd_which_max_abs, x)
  expect_tiers_give(4L, simd_which_min_abs, x)
  for (n in edge_lengths()) {
    seed <- as.integer(n %% 1000) + 1L
    y <- rand_i64(n, seed = seed)
    a <- abs(y)
    expect_tiers_i64(suppressWarnings(max(a, na.rm = TRUE)),
      function(...) suppressWarnings(simd_max_abs(...)), y,
      na.rm = TRUE
    )
    expect_tiers_give(which_i64(a, max), simd_which_max_abs, y)
    expect_tiers_give(which_i64(a, min), simd_which_min_abs, y)
  }
  expect_warning(
    expect_i64(
      simd_max_abs(simd_as_integer64(c(NA, NA)), na.rm = TRUE), i64_dec("-9223372036854775807")
    ),
    "no non-NA value"
  )
})

test_that("simd_prod_sums and simd_prod_diffs are prod(x + y) and prod(x - y)", {
  for (n in edge_lengths()) {
    if (n > 5000) next
    seed <- as.integer(n %% 1000) + 1L
    x <- near_one(n, seed)
    y <- rand_vec("double", n, na_frac = 0, nan_frac = 0, inf_frac = 0, seed = seed) / 1e6
    yi <- rand_vec("integer", n, seed = seed) %% 3L
    for (na_rm in c(FALSE, TRUE)) {
      for (mode in c("fast", "compensated")) {
        expect_simd_matches_base(simd_prod_sums, function(x, y, ...) prod(x + y, ...), x, y,
          na.rm = na_rm, tolerance = prod_tol(n), precision = mode
        )
        expect_simd_matches_base(simd_prod_diffs, function(x, y, ...) prod(x - y, ...), x, y,
          na.rm = na_rm, tolerance = prod_tol(n), precision = mode
        )
      }
      # Integer and logical operands, and broadcasts.
      xi <- as.integer(round(x))
      expect_simd_matches_base(simd_prod_sums, function(x, y, ...) prod(x + y, ...), xi, yi,
        na.rm = na_rm, tolerance = prod_tol(n)
      )
      expect_simd_matches_base(simd_prod_diffs, function(x, y, ...) prod(x - y, ...), x, 1e-6,
        na.rm = na_rm, tolerance = prod_tol(n)
      )
      expect_simd_matches_base(simd_prod_sums, function(x, y, ...) prod(x + y, ...), 1, y,
        na.rm = na_rm, tolerance = prod_tol(n)
      )
    }
  }
})

test_that("the products of sums follow simd_prod's rules for missing values and types", {
  # A pair is missing when its sum is: NA over NaN, and na.rm removes it,
  # Inf - Inf included.
  expect_tiers_give(NA_real_, simd_prod_sums, c(1, NA, NaN, 2), c(1, 1, 1, 1))
  expect_tiers_give(NaN, simd_prod_sums, c(1, NaN, 2), 1)
  expect_tiers_give(NaN, simd_prod_sums, c(Inf, 2), c(-Inf, 1))
  expect_tiers_give(3, simd_prod_sums, c(Inf, 2), c(-Inf, 1), na.rm = TRUE)
  expect_tiers_give(3, simd_prod_diffs, c(Inf, 2, NA), c(Inf, -1, 1), na.rm = TRUE)
  expect_tiers_give(NA_real_, simd_prod_diffs, c(5L, NA), 1L)
  expect_identical(simd_prod_sums(numeric(0), numeric(0)), 1)
  expect_identical(simd_prod_sums(numeric(0), 1), 1)
  expect_identical(simd_prod_sums(c(TRUE, TRUE), TRUE), 4)
  # Integer sums are computed in double, so they do not overflow.
  m <- .Machine$integer.max
  expect_identical(simd_prod_sums(m, 1L), m + 1)
  expect_identical(simd_prod_diffs(-m, 1L), -m - 1)
  expect_error(
    simd_prod_sums(1:3, 1:2),
    "lengths of 'x' \\(3\\) and 'y' \\(2\\) must be equal or one of them must be 1"
  )
  if (has_bit64()) {
    expect_error(
      simd_prod_sums(simd_as_integer64(1), 1),
      "simd_prod_sums\\(\\) does not support 'x' of type integer64"
    )
  }
})

test_that("the products of sums of complex numbers are prod(x + y) of complex", {
  z <- complex(real = c(1, 0.5, -2), imaginary = c(1, -1, 0.25))
  w <- complex(real = c(0.5, 1, 1), imaginary = c(0, 2, -1))
  for (mode in c("fast", "pairwise", "compensated")) {
    old <- simd_precision(mode)
    expect_equal(simd_prod_sums(z, w), prod(z + w), tolerance = 1e-14)
    expect_equal(simd_prod_diffs(z, w), prod(z - w), tolerance = 1e-14)
    expect_identical(simd_prod_sums(z, w), simd_prod(simd_add(z, w)))
    expect_equal(simd_prod_sums(z, 2), prod(z + 2), tolerance = 1e-14)
    expect_equal(simd_prod_diffs(1L, z), prod(1L - z), tolerance = 1e-14)
    simd_precision(old)
  }
  expect_identical(simd_prod_sums(c(z, NA), 1), complex(real = NA, imaginary = NA))
  expect_equal(simd_prod_sums(c(z, NA), 1, na.rm = TRUE), prod(z + 1), tolerance = 1e-14)
  # Long vectors go through the blocks of sums.
  zz <- complex(real = 1 + seq_len(3000) / 1e5, imaginary = 1 / seq_len(3000) / 1e3)
  expect_equal(simd_prod_sums(zz, 1e-6), prod(zz + 1e-6), tolerance = 1e-11)
})

test_that("the new reductions take simd_vec operands and return plain results", {
  v <- simd_vec(c(-3, 1, 2))
  expect_identical(simd_max_abs(v), 3)
  expect_identical(simd_which_min_abs(v), 2L)
  expect_identical(simd_prod_sums(v, simd_vec(c(1, 1, 1))), -12)
  expect_identical(simd_count(simd_vec(c(TRUE, NA))), NA_integer_)
})
