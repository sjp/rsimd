# Complex numbers: add, sub, neg, conj, re, im, sum and the predicates,
# on every tier, against base R.

modes <- c("fast", "pairwise", "compensated")

# Complex numbers whose parts are independent random doubles with NA, NaN,
# +-Inf (as rand_vec() makes them) and signed zeros.
rand_cplx <- function(n, seed = 1L, ...) {
  re <- rand_vec("double", n, seed = seed, ...)
  im <- rand_vec("double", n, seed = seed + 1000L, ...)
  with_seed(seed, {
    k <- stats::rbinom(2, n, 0.05)
    re[sample.int(n, k[1])] <- rep_len(c(0, -0), k[1])
    im[sample.int(n, k[2])] <- rep_len(c(-0, 0), k[2])
  })
  complex(real = re, imaginary = im)
}

# Complex lengths around the number of complex values per vector (half the
# 64-bit lanes: 1 to 16), the four-vector blocks of the folds, the region
# size and the sum block.
cplx_lengths <- c(0:18, 31:35, 63:67, 511:513, 4095:4097)

# x op y for doubles part by part, with NA_real_ wherever an operand's part
# is NA: the rsimd rule for missing values (base R may give NaN for
# NaN + NA on some CPUs).
part_op <- function(op, x, y) {
  complex(
    real = na_merged(op(Re(x), Re(y)), Re(x), Re(y)),
    imaginary = na_merged(op(Im(x), Im(y)), Im(x), Im(y))
  )
}

test_that("add, sub, neg, conj, re and im match base R exactly on every tier", {
  for (n in cplx_lengths) {
    x <- rand_cplx(n, seed = n + 1L)
    y <- rand_cplx(n, seed = n + 2L)
    expect_tiers_give(part_op(`+`, x, y), simd_add, x, y)
    expect_tiers_give(part_op(`-`, x, y), simd_sub, x, y)
    expect_tiers_give(-x, simd_neg, x)
    expect_tiers_give(Conj(x), simd_conj, x)
    expect_tiers_give(Re(x), simd_re, x)
    expect_tiers_give(Im(x), simd_im, x)
  }
})

test_that("without NaN, add and sub equal base R's own result", {
  for (n in c(1, 5, 17, 100)) {
    x <- rand_cplx(n, seed = n, nan_frac = 0)
    y <- rand_cplx(n, seed = n + 7L, nan_frac = 0)
    expect_tiers_give(x + y, simd_add, x, y)
    expect_tiers_give(x - y, simd_sub, x, y)
  }
})

test_that("missing parts behave as in base R", {
  na_1 <- complex(real = NA, imaginary = 1)
  one_nan <- complex(real = 1, imaginary = NaN)
  expect_identical(simd_is_na(c(na_1, one_nan, 1 + 1i)), c(TRUE, TRUE, FALSE))
  expect_identical(simd_is_nan(c(NA_complex_, na_1, one_nan)), c(FALSE, FALSE, TRUE))
  expect_identical(simd_is_nan(complex(real = NaN, imaginary = NA)), TRUE)
  expect_identical(is.na(simd_im(simd_conj(complex(real = 1, imaginary = NA)))), TRUE)
  expect_identical(simd_re(simd_neg(complex(real = 1, imaginary = NA))), -1)
  expect_identical(simd_re(na_1), NA_real_)
  expect_identical(simd_im(na_1), 1)
  expect_identical(simd_im(simd_conj(na_1)), -1)
  expect_identical(is.na(simd_re(simd_conj(na_1))), TRUE)
  expect_identical(is.na(simd_re(simd_neg(na_1))), TRUE)
  expect_identical(simd_im(simd_neg(na_1)), -1)
  # Conj and neg flip the sign of a zero part.
  expect_identical(1 / simd_im(simd_conj(complex(real = 1, imaginary = 0))), -Inf)
  expect_identical(1 / simd_im(simd_neg(complex(real = 1, imaginary = 0))), -Inf)
  expect_identical(simd_neg(1 + 2i), -1 - 2i)
  expect_identical(simd_sub(1 + 2i, 3 + 4i), -2 - 2i)
})

test_that("mixed operands are promoted to complex and length-1 operands broadcast", {
  expect_identical(simd_add(c(1, 2), 1i), c(1 + 1i, 2 + 1i))
  for (n in c(1, 2, 3, 9, 33, 600)) {
    z <- rand_cplx(n, seed = n + 3L)
    d <- rand_vec("double", n, seed = n + 4L)
    i <- rand_vec("integer", n, seed = n + 5L)
    l <- rand_vec("logical", n, seed = n + 6L)
    for (w in list(d, i, l)) {
      expect_tiers_give(part_op(`+`, z, as.complex(w)), simd_add, z, w)
      expect_tiers_give(part_op(`-`, as.complex(w), z), simd_sub, w, z)
    }
    s <- complex(real = 2.5, imaginary = -1)
    expect_tiers_give(part_op(`+`, z, s), simd_add, z, s)
    expect_tiers_give(part_op(`-`, s, z), simd_sub, s, z)
    expect_tiers_give(part_op(`+`, z, 3 + 0i), simd_add, z, 3)
    expect_tiers_give(part_op(`-`, 2 + 0i, z), simd_sub, 2L, z)
    expect_tiers_give(part_op(`+`, z, NA_complex_), simd_add, z, NA_complex_)
    expect_tiers_give(part_op(`+`, as.complex(d), s), simd_add, d, s)
  }
  # integer and logical NA become NA + 0i, as with as.complex().
  expect_identical(simd_im(simd_add(1i, NA_integer_)), 1)
  expect_identical(simd_im(simd_add(NA, 1i)), 1)
  expect_identical(simd_add(complex(0), 1), complex(0))
  expect_identical(simd_add(1i, complex(0)), complex(0))
  expect_error(simd_add(c(1i, 2i), 1:3),
    "lengths of 'x' (2) and 'y' (3) must be equal or one of them must be 1",
    fixed = TRUE
  )
})

test_that("na_check = FALSE gives the same result without missing values", {
  x <- rand_cplx(100, seed = 9L, na_frac = 0, nan_frac = 0)
  y <- rand_cplx(100, seed = 10L, na_frac = 0, nan_frac = 0)
  expect_tiers_give(x + y, simd_add, x, y, na_check = FALSE)
  expect_tiers_give(x - y, simd_sub, x, y, na_check = FALSE)
})

# NA_real_ if x has an NA, else NaN if it has a NaN, else NULL: the rsimd
# rule for a sum with na.rm = FALSE (base R's choice depends on the order).
part_missing <- function(x) {
  if (any(is.na(x) & !is.nan(x))) NA_real_ else if (any(is.nan(x))) NaN
}

# Expects simd_sum(z, na.rm) in every precision mode and tier to agree,
# part by part, with the none tier within the mode's bound and with base
# R's sum within that bound plus base R's own rounding, missing values
# following part_missing().
expect_complex_sum <- function(z, na.rm = FALSE) {
  keep <- if (na.rm) !is.na(z) else rep(TRUE, length(z))
  base <- sum(z, na.rm = na.rm)
  base_eps <- if (has_wide_long_double()) 2^-63 else eps
  old <- simd_precision()
  on.exit(simd_precision(old))
  problems <- character(0)
  for (mode in modes) {
    simd_precision(mode)
    res <- with_each_tier(function() simd_sum(z, na.rm = na.rm))
    for (part in c("real", "imaginary")) {
      f <- if (part == "real") Re else Im
      x <- f(z)[keep]
      want <- part_missing(x)
      if (is.null(want)) want <- f(base)
      bound <- mode_bound(x, mode)
      for (tier in names(res)) {
        got <- res[[tier]]
        problem <- if (!is.complex(got) || length(got) != 1L) {
          "not a single complex number"
        } else {
          double_mismatch(f(got), f(res[["none"]]), bound, 2^-50)
        }
        if (is.null(problem)) {
          problem <- double_mismatch(f(got), want,
            bound + max(length(x), 1) * base_eps * abs_mass(x), max(1e-12, 2^-50)
          )
        }
        if (!is.null(problem)) {
          problems <- c(problems, sprintf(
            "tier %s, %s mode, %s part, n = %d, na.rm = %s: %s", tier, mode, part, length(z),
            na.rm, problem
          ))
        }
      }
    }
  }
  testthat::expect(length(problems) == 0L, paste(problems, collapse = "\n"))
}

test_that("sum matches the oracle and base R part by part in every mode", {
  for (n in cplx_lengths) {
    seed <- n + 11L
    zs <- list(
      rand_cplx(n, seed = seed),
      rand_cplx(n, seed = seed, na_frac = 0, nan_frac = 0, inf_frac = 0)
    )
    for (z in zs) {
      expect_complex_sum(z)
      expect_complex_sum(z, na.rm = TRUE)
    }
  }
})

test_that("each part of the sum is the sum of that part as a double vector", {
  z <- rand_cplx(10000, seed = 5L, na_frac = 0, nan_frac = 0, inf_frac = 0)
  old <- simd_precision("pairwise")
  on.exit(simd_precision(old))
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      s <- simd_sum(z)
      expect_identical(Re(s), simd_sum(Re(z)), info = tier)
      expect_identical(Im(s), simd_sum(Im(z)), info = tier)
    })
  }
})

test_that("sum propagates missing parts separately and na.rm drops whole elements", {
  na_1 <- complex(real = NA, imaginary = 1)
  one_nan <- complex(real = 1, imaginary = NaN)
  for (mode in modes) {
    old <- simd_precision(mode)
    expect_tiers_give(complex(real = NA, imaginary = 2), simd_sum, c(1 + 2i, NA))
    expect_tiers_give(complex(real = 2, imaginary = NaN), simd_sum, c(1 + 2i, one_nan))
    expect_tiers_give(complex(real = NA, imaginary = NaN), simd_sum, c(one_nan, na_1))
    expect_tiers_give(4 + 3i, simd_sum, c(1 + 2i, NA, 3 + 1i, one_nan), na.rm = TRUE)
    expect_tiers_give(0 + 0i, simd_sum, complex(0))
    expect_tiers_give(0 + 0i, simd_sum, c(na_1, one_nan), na.rm = TRUE)
    expect_tiers_give(NA_complex_, simd_sum, NA_complex_)
    simd_precision(old)
  }
  # Compensated mode keeps what fast mode may cancel, as base R's long
  # double does.
  old <- simd_precision("compensated")
  on.exit(simd_precision(old))
  expect_tiers_give(1 + 1i, simd_sum, c(1e16 + 0i, 1 + 1i, -1e16 + 0i))
})

test_that("predicates match base R on every tier", {
  preds <- list(
    na = list(simd_is_na, simd_is_na_any, simd_is_na_all, is.na),
    nan = list(simd_is_nan, simd_is_nan_any, simd_is_nan_all, is.nan),
    finite = list(simd_is_finite, simd_is_finite_any, simd_is_finite_all, is.finite),
    infinite = list(simd_is_infinite, simd_is_infinite_any, simd_is_infinite_all, is.infinite)
  )
  sp <- c(NA, NaN, Inf, -Inf, 0, 1)
  grid <- all_pairs(sp)
  specials <- complex(real = grid$x, imaginary = grid$y)
  for (n in cplx_lengths) {
    zs <- list(rand_cplx(n, seed = n + 21L), rep(NA_complex_, n), rep(1i, n))
    for (z in zs) {
      for (p in preds) {
        want <- p[[4]](z)
        expect_tiers_give(want, p[[1]], z)
        expect_tiers_give(any(want), p[[2]], z)
        expect_tiers_give(all(want), p[[3]], z)
      }
    }
  }
  for (p in preds) expect_tiers_give(p[[4]](specials), p[[1]], specials)
})

test_that("missing-value queries cover complex lengths around the vector width", {
  for (n in cplx_lengths) {
    z <- rand_cplx(n, seed = n + 31L)
    expect_tiers_give(anyNA(z), simd_any_na, z)
    expect_tiers_give(as.double(sum(is.na(z))), simd_count_na, z)
    expect_tiers_give(which(is.na(z)), simd_which_na, z)
  }
})

test_that("conj, re and im of double, integer and logical input match base R", {
  for (x in list(c(1.5, -0, NA, NaN, Inf), c(1L, NA, -3L), c(TRUE, NA, FALSE), double(0))) {
    expect_tiers_give(Conj(x), simd_conj, x)
    expect_tiers_give(Re(x), simd_re, x)
    expect_tiers_give(Im(x), simd_im, x)
  }
  expect_null(attributes(simd_re(c(a = 1 + 1i))))
  expect_null(attributes(simd_conj(c(a = 1.5))))
  expect_error(simd_re(as.raw(1)), "non-numeric argument to function", fixed = TRUE)
  expect_error(simd_re("a"), "'z' must be an atomic vector")
  expect_error(simd_conj(structure(0, class = "integer64")),
    "simd_conj() does not support 'z' of type integer64",
    fixed = TRUE
  )
})

test_that("chunk boundaries give the unchunked result", {
  n <- 2^20 + 7
  z <- complex(real = as.double(seq_len(n)), imaginary = -as.double(seq_len(n)) / 3)
  expect_identical(simd_add(z, 1i), z + 1i)
  expect_identical(simd_sub(z, z), complex(real = rep(0, n), imaginary = rep(0, n)))
  expect_identical(simd_conj(z), Conj(z))
  expect_identical(simd_im(z), Im(z))
  expect_identical(simd_is_na(z), logical(n))
  expect_complex_sum(z)
  expect_identical(simd_re(as.complex(seq_len(5000))), as.double(seq_len(5000)))
})

test_that("a small interrupt stride gives the same results", {
  skip_on_cran()
  skip_if_no_subprocess()
  z <- rand_cplx(5000, seed = 41L)
  w <- rand_cplx(5000, seed = 42L)
  child <- function(z, w) {
    library(rsimd)
    lapply(stats::setNames(nm = simd_available()), function(t) {
      simd_with_impl(t, {
        old <- simd_precision("pairwise")
        on.exit(simd_precision(old))
        list(
          add = simd_add(z, w), neg = simd_neg(z), conj = simd_conj(z), re = simd_re(z),
          na = simd_is_na(z), sum = simd_sum(z, na.rm = TRUE), count = simd_count_na(z)
        )
      })
    })
  }
  here <- child(z, w)
  env <- c(callr::rcmd_safe_env(), RSIMD_DEBUG_STRIDE = "128")
  there <- callr::r(child, list(z, w), env = env)
  for (tier in intersect(names(here), names(there))) {
    for (nm in names(here[[tier]])) {
      expect_true(same_values(there[[tier]][[nm]], here[[tier]][[nm]]), info = paste(tier, nm))
    }
  }
})

test_that("other functions reject complex input", {
  expect_error(simd_mul(1i, 1i), "simd_mul() does not support 'x' of type complex", fixed = TRUE)
  expect_error(simd_exp(1i), "simd_exp() does not support 'x' of type complex", fixed = TRUE)
  expect_error(simd_max(1i), "invalid 'type' (complex) of argument", fixed = TRUE)
  expect_error(simd_is_zero(1i), "simd_is_zero() does not support 'x' of type complex",
    fixed = TRUE
  )
  expect_error(simd_add(1i, structure(0, class = "integer64")),
    "simd_add() cannot combine complex and integer64 operands",
    fixed = TRUE
  )
  expect_error(simd_add(1i, as.raw(1)), "non-numeric argument to binary operator", fixed = TRUE)
  expect_error(simd_sub(as.raw(1), 1i), "non-numeric argument to binary operator", fixed = TRUE)
})
