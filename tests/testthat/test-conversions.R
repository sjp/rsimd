# Type conversions: checked mode identical to base R, warnings included;
# saturating and truncating against exact references; on every tier.

int_msg <- "NAs introduced by coercion to integer range"
raw_msg <- "out-of-range values treated as 0 in coercion to raw"

# Exact references of the saturating and truncating modes.
sat_int <- function(x) {
  t <- trunc(x)
  as.integer(ifelse(is.na(x), NA, pmin(pmax(t, -.Machine$integer.max), .Machine$integer.max)))
}
trunc_int <- function(x) {
  t <- trunc(x)
  w <- t - 2^32 * floor(t / 2^32)
  w[!is.finite(x)] <- NA
  from_u32(w)
}
sat_raw <- function(x) {
  t <- trunc(as.double(x))
  as.raw(ifelse(is.na(t), 0, pmin(pmax(t, 0), 255)))
}
trunc_raw <- function(x) {
  t <- trunc(as.double(x))
  w <- t - 256 * floor(t / 256)
  as.raw(ifelse(is.finite(t), w, 0))
}

conv_doubles <- function() {
  c(
    pred_doubles(), 2147483647.5, -2147483647.5, 2147483647.9, -2147483648.5, 2147483648,
    -2147483648, 2^32 + 5, -(2^32 + 5), 1e19, -1e19, 2.9, -2.9, 1954, 255.9, 256.5, -0.5,
    -0.99, -1, 2.5, 3e9, 2^31 + 256, 2^63
  )
}

test_that("every row of base R's coercion table is reproduced with its warnings", {
  expect_tiers_like_base(as.integer, simd_as_integer, c(2^31, -2^31, 2.9, -2.9, NaN, Inf, 1954))
  expect_tiers_like_base(as.integer, simd_as_integer, -2147483647.5)
  expect_tiers_like_base(as.integer, simd_as_integer, 2147483647.9)
  expect_tiers_like_base(as.raw, simd_as_raw, c(256, -1, 2.7, NA, 255))
  expect_tiers_like_base(as.logical, simd_as_logical, c(0, 2.5, NaN, NA, -Inf))
  expect_tiers_like_base(as.logical, simd_as_logical, as.raw(c(0, 2)))
  expect_tiers_like_base(as.integer, simd_as_integer, as.raw(255))
  expect_tiers_like_base(as.integer, simd_as_integer, TRUE)
  # Through the integer range first, as base R: two warnings, in order.
  expect_tiers_like_base(as.raw, simd_as_raw, -Inf)
  expect_tiers_like_base(as.raw, simd_as_raw, c(2^31, NaN, 300))
})

test_that("checked conversions of doubles match base R exactly", {
  x <- conv_doubles()
  expect_tiers_like_base(as.integer, simd_as_integer, x)
  expect_tiers_like_base(as.raw, simd_as_raw, x)
  expect_tiers_like_base(as.logical, simd_as_logical, x)
  expect_tiers_like_base(as.integer, simd_as_integer, c(1.5, -0, 5e-324, 2^53))
  # One warning per call, however many elements trigger it.
  expect_tiers_like_base(as.integer, simd_as_integer, rep(c(1e10, 1), 5000))
  expect_tiers_like_base(as.raw, simd_as_raw, rep(c(1e10, 300, 1), 5000))
})

test_that("the acceptance boundaries hold in all three modes", {
  x <- c(-2147483648.5, -2147483647.5)
  expect_tier_warnings(int_msg, simd_as_integer, x)
  expect_identical(suppressWarnings(simd_as_integer(x)), c(NA, -2147483647L))
  expect_tiers_give(c(-2147483647L, -2147483647L), simd_as_integer, x, "saturating")
  expect_tiers_give(c(NA, -2147483647L), simd_as_integer, x, "truncating")
})

test_that("saturating and truncating modes match the exact references, silently", {
  x <- conv_doubles()
  expect_tier_warnings(character(0), simd_as_integer, x, "saturating")
  expect_tier_warnings(character(0), simd_as_integer, x, "truncating")
  expect_tier_warnings(character(0), simd_as_raw, x, "saturating")
  expect_tier_warnings(character(0), simd_as_raw, x, "truncating")
  expect_tiers_give(sat_int(x), simd_as_integer, x, "saturating")
  expect_tiers_give(trunc_int(x), simd_as_integer, x, "truncating")
  expect_tiers_give(sat_raw(x), simd_as_raw, x, "saturating")
  expect_tiers_give(trunc_raw(x), simd_as_raw, x, "truncating")
  expect_identical(simd_as_integer(c(2^31, 2^32 + 5, Inf), "truncating"), c(NA, 5L, NA))
  expect_identical(simd_as_raw(c(255, 256, -1, -0.5), "truncating"), as.raw(c(255, 0, 255, 0)))
  expect_identical(
    simd_as_raw(c(255, 256, -1, NaN, Inf), "saturating"),
    as.raw(c(255, 255, 0, 0, 255))
  )
})

test_that("integer and logical to raw follow each mode", {
  x <- c(0L, 1L, 255L, 256L, -1L, 1000L, -300L, NA, .Machine$integer.max, -.Machine$integer.max)
  expect_tiers_like_base(as.raw, simd_as_raw, x)
  expect_tiers_like_base(as.raw, simd_as_raw, c(TRUE, FALSE, NA))
  expect_tiers_give(sat_raw(x), simd_as_raw, x, "saturating")
  expect_tiers_give(trunc_raw(x), simd_as_raw, x, "truncating")
  expect_tiers_give(as.raw(c(1, 0, 0)), simd_as_raw, c(TRUE, FALSE, NA), "truncating")
})

test_that("conversions to double and logical, and from raw, match base R", {
  for (x in list(
    c(0L, 1L, -1L, NA, .Machine$integer.max), c(TRUE, FALSE, NA), as.raw(0:255),
    conv_doubles()
  )) {
    expect_tiers_like_base(as.double, simd_as_double, x)
    expect_tiers_like_base(as.logical, simd_as_logical, x)
    expect_tiers_like_base(as.integer, simd_as_integer, x)
  }
  expect_identical(simd_as_integer(as.raw(255), "saturating"), 255L)
})

test_that("input of the target type is copied without attributes", {
  expect_identical(simd_as_integer(c(a = 1L, b = NA)), c(1L, NA))
  expect_identical(simd_as_double(matrix(c(1.5, NA), 1)), c(1.5, NA))
  expect_identical(simd_as_logical(c(x = NA)), NA)
  expect_identical(simd_as_raw(as.raw(1:3)), as.raw(1:3))
  expect_identical(simd_as_integer(c(TRUE, NA)), c(1L, NA))
})

test_that("conversions match base R on random vectors of edge lengths", {
  for (n in sweep_lengths()) {
    x <- rand_vec("double", n, seed = n + 31L) * 1e6
    expect_tiers_like_base(as.integer, simd_as_integer, x)
    expect_tiers_give(trunc_int(x), simd_as_integer, x, "truncating")
    expect_tiers_give(as.logical(x), simd_as_logical, x)
    i <- rand_vec("integer", n, seed = n + 32L) %/% 3000L
    expect_tiers_like_base(as.raw, simd_as_raw, i)
    expect_tiers_give(as.double(i), simd_as_double, i)
    r <- as.raw(with_seed(n, sample.int(256, n, TRUE) - 1L))
    expect_tiers_give(as.double(r), simd_as_double, r)
  }
  for (x in altrep_inputs(5000)) {
    expect_tiers_give(as.integer(x), simd_as_integer, x)
    expect_tiers_give(as.double(x), simd_as_double, x)
  }
})

test_that("argument and type errors", {
  expect_error(simd_as_integer(1, mode = "wrap"),
    "'mode' must be one of \"checked\", \"saturating\" or \"truncating\"",
    fixed = TRUE
  )
  expect_error(simd_as_integer(1, mode = ""), "'mode' must be one of")
  expect_error(simd_as_integer(1, mode = NA_character_), "'mode' must be one of")
  expect_error(simd_as_integer(2.5, c("checked", "saturating")), "'mode' must be a single string")
  expect_error(simd_as_integer64(1, mode = 1), "'mode' must be a single string")
  expect_error(simd_as_raw(1, mode = character(0)), "'mode' must be a single string")
  e <- tryCatch(simd_as_integer(1, mode = "wrap"), error = identity)
  expect_identical(conditionCall(e), quote(simd_as_integer(1, mode = "wrap")))
  expect_identical(simd_as_integer(1e10, "sat"), .Machine$integer.max)
  expect_identical(
    simd_as_integer(2.5, c("checked", "saturating", "truncating")),
    simd_as_integer(2.5)
  )
  expect_error(simd_as_integer(1i), "simd_as_integer() does not support 'x' of type complex",
    fixed = TRUE
  )
  expect_error(simd_as_raw("1"), "must be an atomic vector")
})
