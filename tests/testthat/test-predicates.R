# Predicates and their _any/_all forms: identical to base R on every tier,
# NA kept apart from NaN for every NaN payload.

# Base R versions of the predicates (negative: sign bit set, not NaN).
base_pred <- list(
  na = is.na,
  nan = function(x) if (is.double(x)) is.nan(x) else rep(FALSE, length(x)),
  finite = is.finite,
  infinite = is.infinite,
  negative = function(x) {
    if (is.double(x)) !is.na(x) & (x < 0 | (x == 0 & 1 / x < 0)) else !is.na(x) & x < 0
  },
  zero = function(x) !is.na(x) & x == 0
)
simd_pred <- function(op, form = "") get(paste0("simd_is_", op, form))

test_that("double predicates match base R for every NaN payload, zero and denormal", {
  x <- pred_doubles()
  for (op in names(base_pred)) {
    expected <- base_pred[[op]](x)
    expect_tiers_give(expected, simd_pred(op), x)
    expect_tiers_give(any(expected), simd_pred(op, "_any"), x)
    expect_tiers_give(all(expected), simd_pred(op, "_all"), x)
  }
  # NA and NaN are told apart whatever the payload.
  p <- nan_payloads()
  expect_tiers_give(is.na(p) & !is.nan(p), function(x) simd_is_na(x) & !simd_is_nan(x), p)
  expect_identical(simd_is_nan(NA_real_), FALSE)
  expect_identical(simd_is_na(NaN), TRUE)
  expect_identical(
    simd_is_negative(c(-0, -Inf, NaN, NA, -5e-324)),
    c(TRUE, TRUE, FALSE, FALSE, TRUE)
  )
  expect_identical(simd_is_zero(c(-0, 0, 5e-324)), c(TRUE, TRUE, FALSE))
})

test_that("integer and logical predicates match base R; NA only at INT_MIN", {
  ints <- c(0L, 1L, -1L, .Machine$integer.max, -.Machine$integer.max, NA, 7L, NA)
  lgls <- c(TRUE, FALSE, NA, FALSE)
  for (x in list(ints, lgls)) {
    for (op in names(base_pred)) {
      expected <- base_pred[[op]](x)
      expect_tiers_give(expected, simd_pred(op), x)
      expect_tiers_give(any(expected), simd_pred(op, "_any"), x)
      expect_tiers_give(all(expected), simd_pred(op, "_all"), x)
    }
  }
  expect_identical(simd_is_na(c(1L, NA, 3L)), c(FALSE, TRUE, FALSE))
  expect_identical(simd_is_finite(NA_integer_), FALSE)
  expect_identical(simd_is_infinite(NA), FALSE)
  expect_identical(simd_is_nan(1L), FALSE)
})

test_that("raw is never missing, finite or infinite, as in base R", {
  x <- as.raw(c(0, 1, 255))
  for (op in c("na", "nan", "finite", "infinite")) {
    expect_tiers_give(base_pred[[op]](x), simd_pred(op), x)
    expect_identical(simd_pred(op, "_any")(x), FALSE)
    expect_identical(simd_pred(op, "_all")(x), FALSE)
  }
  expect_identical(simd_is_finite(as.raw(0)), is.finite(as.raw(0)))
  expect_error(simd_is_zero(x), "simd_is_zero() does not support 'x' of type raw", fixed = TRUE)
  expect_error(simd_is_negative_any(x), "does not support 'x' of type raw", fixed = TRUE)
})

test_that("_any is FALSE and _all is TRUE for empty input", {
  for (op in names(base_pred)) {
    for (x in list(double(0), integer(0), logical(0))) {
      expect_identical(simd_pred(op, "_any")(x), FALSE)
      expect_identical(simd_pred(op, "_all")(x), TRUE)
      expect_identical(simd_pred(op)(x), logical(0))
    }
  }
  expect_identical(simd_is_na_all(raw(0)), TRUE)
})

test_that("predicates match base R on random vectors of edge lengths", {
  for (n in edge_lengths()) {
    for (type in c("double", "integer")) {
      x <- rand_vec(type, n, seed = n + 3L)
      for (op in names(base_pred)) {
        expected <- base_pred[[op]](x)
        expect_tiers_give(expected, simd_pred(op), x)
        expect_tiers_give(any(expected), simd_pred(op, "_any"), x)
        expect_tiers_give(all(expected), simd_pred(op, "_all"), x)
      }
    }
  }
})

test_that("early exit does not change _any/_all when the last element decides", {
  n <- 2^20 + 7
  fin <- rep(1.5, n)
  last_na <- fin
  last_na[n] <- NA
  last_zero <- rep(3L, n)
  last_zero[n] <- 0L
  expect_tiers_give(FALSE, simd_is_na_any, fin)
  expect_tiers_give(TRUE, simd_is_na_any, last_na)
  expect_tiers_give(TRUE, simd_is_finite_all, fin)
  expect_tiers_give(FALSE, simd_is_finite_all, last_na)
  expect_tiers_give(TRUE, simd_is_zero_any, last_zero)
  expect_tiers_give(FALSE, simd_is_zero_all, last_zero)
  expect_tiers_give(FALSE, simd_is_zero_any, rep(3L, n))
  expect_tiers_give(TRUE, simd_is_negative_all, -fin)
  expect_tiers_give(is.na(last_na), simd_is_na, last_na)
})

test_that("compact sequences are read without expanding them", {
  for (x in altrep_inputs(5000)) {
    expect_true(takes_region_path(x))
    expect_tiers_give(rep(FALSE, 5000), simd_is_na, x)
    expect_tiers_give(TRUE, simd_is_finite_all, x)
    expect_tiers_give(FALSE, simd_is_zero_any, x)
  }
})

test_that("attributes are dropped and unsupported types rejected", {
  x <- c(a = 1, b = NA)
  expect_identical(simd_is_na(x), c(FALSE, TRUE))
  expect_identical(simd_is_na(matrix(c(1L, NA), 1)), c(FALSE, TRUE))
  expect_error(simd_is_negative(1i), "simd_is_negative() does not support 'x' of type complex",
    fixed = TRUE
  )
  expect_error(simd_is_zero_any(1i), "does not support 'x' of type complex$")
  expect_error(simd_is_zero("a"), "must be an atomic vector")
  expect_error(simd_is_na(factor("a")), "must be an atomic vector")
})
