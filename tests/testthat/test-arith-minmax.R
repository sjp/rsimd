# pmin, pmax and clamp: bit-identical to base R's pmin/pmax on every tier,
# including which missing value wins and the sign of zero.

edge_mm <- function() c(0, -0, 1, -1, 2.5, -Inf, Inf, NA, NaN, 5e-324, 1e308)

test_that("pmin and pmax match base R for every pair of edge values", {
  p <- all_pairs(edge_mm())
  for (na_rm in c(FALSE, TRUE)) {
    expect_tiers_give(pmin(p$x, p$y, na.rm = na_rm), simd_pmin, p$x, p$y, na.rm = na_rm)
    expect_tiers_give(pmax(p$x, p$y, na.rm = na_rm), simd_pmax, p$x, p$y, na.rm = na_rm)
  }
  expect_tiers_give(pmin(p$x, p$y, na.rm = TRUE), simd_pmin_num, p$x, p$y)
  expect_tiers_give(pmax(p$x, p$y, na.rm = TRUE), simd_pmax_num, p$x, p$y)
})

test_that("base R's NA, NaN and signed-zero rules reproduce", {
  kind <- function(v) ifelse(is.na(v) & !is.nan(v), "NA", ifelse(is.nan(v), "NaN", "num"))
  expect_identical(kind(simd_pmax(NaN, NA)), "NA")
  expect_identical(kind(simd_pmax(NA, NaN)), "NaN")
  expect_identical(kind(simd_pmax(c(1, NaN), c(NA, 2))), c("NA", "NaN"))
  expect_identical(simd_pmax(c(1, NaN), c(NA, 2), na.rm = TRUE), c(1, 2))
  expect_identical(simd_pmax_num(NaN, 1), 1)
  expect_identical(kind(simd_pmax_num(NaN, NA)), "NA")
  expect_identical(1 / simd_pmin(c(0, -0), c(-0, 0)), c(Inf, -Inf))
  expect_identical(simd_pmax(1:3, 2L), pmax(1:3, 2L))
  expect_type(simd_pmax(1:3, 2), "double")
  expect_type(simd_pmin(TRUE, FALSE), "integer")
})

test_that("integer pmin and pmax match base R", {
  vals <- c(0L, 1L, -1L, .Machine$integer.max, -.Machine$integer.max, NA)
  p <- all_pairs(vals)
  for (na_rm in c(FALSE, TRUE)) {
    expect_tiers_give(pmin(p$x, p$y, na.rm = na_rm), simd_pmin, p$x, p$y, na.rm = na_rm)
    expect_tiers_give(pmax(p$x, p$y, na.rm = na_rm), simd_pmax, p$x, p$y, na.rm = na_rm)
  }
  x <- rand_vec("integer", 1003, seed = 51L)
  y <- rand_vec("integer", 1003, seed = 52L)
  expect_tiers_give(pmax(x, y), simd_pmax, x, y)
  expect_tiers_give(pmin(x, y, na.rm = TRUE), simd_pmin_num, x, y)
  expect_tiers_give(pmax(x, 7L), simd_pmax, x, 7L)
  expect_tiers_give(pmin(TRUE, x), simd_pmin, TRUE, x)
})

test_that("pmin and pmax of mixed and random doubles match base R", {
  for (n in c(1, 3, 4, 5, 8, 9, 17, 4097)) {
    x <- rand_vec("double", n, seed = n + 60L)
    y <- rand_vec("double", n, seed = n + 61L)
    expect_tiers_give(pmin(x, y), simd_pmin, x, y)
    expect_tiers_give(pmax(x, y, na.rm = TRUE), simd_pmax, x, y, na.rm = TRUE)
  }
  i <- c(1L, NA, -3L, 4L)
  d <- c(NaN, 2, NA, 4)
  expect_tiers_give(pmax(i, d), simd_pmax, i, d)
  expect_tiers_give(pmin(d, i, na.rm = TRUE), simd_pmin, d, i, na.rm = TRUE)
  expect_tiers_give(pmin(i, 2.5), simd_pmin, i, 2.5)
})

test_that("clamp is pmin(pmax(x, lo), hi)", {
  x <- c(edge_mm(), -5, 5)
  for (lohi in list(c(-1, 1), c(0, 0), c(-Inf, Inf), c(-0, 0))) {
    expect_tiers_give(pmin(pmax(x, lohi[1]), lohi[2]), simd_clamp, x, lohi[1], lohi[2])
  }
  lo <- c(-1, NA, 0, NaN, -2, 0)
  hi <- c(1, 2, NA, 3, -1, NaN)
  x <- c(5, 1, 1, 1, NA, -3)
  expect_tiers_give(pmin(pmax(x, lo), hi), simd_clamp, x, lo, hi)
  xi <- c(-10L, 0L, 10L, NA)
  expect_tiers_give(pmin(pmax(xi, -3L), 3L), simd_clamp, xi, -3L, 3L)
  expect_tiers_give(pmin(pmax(xi, c(-1L, NA, 1L, 0L)), 5L), simd_clamp, xi, c(-1L, NA, 1L, 0L), 5L)
  expect_tiers_give(pmin(pmax(xi, -3), 3.5), simd_clamp, xi, -3, 3.5)
  expect_type(simd_clamp(1:3, 1L, 2L), "integer")
  expect_tiers_give(c(1, 1, 2), simd_clamp, 0:2, 1, 2)
})

test_that("clamp errors when lo > hi anywhere", {
  msg <- "'lo' must not be greater than 'hi'"
  expect_error(simd_clamp(1:3, 3L, 1L), msg, fixed = TRUE)
  expect_error(simd_clamp(1:3, c(0, 0, 5), 4), msg, fixed = TRUE)
  expect_error(simd_clamp(c(1, 2), 0L, c(1L, -1L)), msg, fixed = TRUE)
  for (tier in simd_available()) {
    expect_error(simd_with_impl(tier, simd_clamp(seq_len(37), c(rep(0, 36), 2), 1)), msg,
      fixed = TRUE, info = tier
    )
  }
  # Missing bounds are not compared.
  expect_identical(simd_clamp(c(1, 2), c(NA, 0), c(0, NaN)), c(NA, NaN))
  expect_identical(simd_clamp(c(1L, 2L), c(NA, 0L), c(0L, NA)), c(NA_integer_, NA_integer_))
})

test_that("bad arguments are rejected", {
  expect_error(simd_pmin(1, 2, na.rm = NA), "'na.rm' must be TRUE or FALSE")
  expect_error(simd_pmax(1:3, 1:2), "lengths of 'x' (3) and 'y' (2)", fixed = TRUE)
  expect_error(simd_clamp(1:3, 1:2, 3), "'lo' (2)", fixed = TRUE)
  expect_error(simd_pmax(as.raw(1), as.raw(2)), "non-numeric argument to binary operator")
  expect_error(simd_pmin(1i, 1), "simd_pmin() does not support 'x' of type complex",
    fixed = TRUE
  )
})

test_that("chunk boundaries and ALTREP inputs give the oracle's result", {
  n <- 2^20 + 7
  x <- seq_len(n)
  expect_identical(simd_pmax(x, 500000L), pmax(x, 500000L))
  expect_simd_identical(simd_pmin, as.double(x), 777.5)
  expect_identical(simd_clamp(x, 10L, 1000000L), pmin(pmax(x, 10L), 1000000L))
})
