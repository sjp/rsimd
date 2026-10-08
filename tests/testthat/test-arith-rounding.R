# floor, ceiling, trunc and round: bit-identical to base R on every tier,
# including -0, values of 2^52 and more, NaN and infinities.

edge_round <- function() {
  big <- 2^52
  c(
    0, -0, 0.5, -0.5, 1.5, -1.5, 2.5, -2.5, 0.4, -0.4, 0.6, -0.6, 0.49999999999999994,
    -0.49999999999999994, 1 - 2^-53, -(1 - 2^-53), big - 0.5, -(big - 0.5), big, -big,
    big + 1, 2^53 + 2, 1e300, -1e300, 5e-324, -5e-324, Inf, -Inf, NA, NaN, 2.675, 1234.5, pi
  )
}

test_that("floor, ceiling, trunc and round(x) match base R bit for bit", {
  x <- edge_round()
  expect_tiers_give(floor(x), simd_floor, x)
  expect_tiers_give(ceiling(x), simd_ceiling, x)
  expect_tiers_give(trunc(x), simd_trunc, x)
  expect_tiers_give(round(x), simd_round, x)
  expect_identical(1 / simd_round(c(-0.5, -0.4, 0.5)), c(-Inf, -Inf, Inf))
  expect_identical(1 / simd_ceiling(-0.5), -Inf)
  expect_identical(1 / simd_trunc(-0.5), -Inf)
  expect_identical(simd_round(c(0.5, 1.5, 2.5, -2.5)), c(0, 2, 2, -2))
})

test_that("rounding random doubles of edge lengths matches base R", {
  for (n in c(1, 3, 4, 5, 8, 9, 17, 4097)) {
    x <- rand_vec("double", n, seed = n + 70L) * 3
    x[seq(1, n, by = 3)] <- round(x[seq(1, n, by = 3)]) + 0.5
    expect_tiers_give(floor(x), simd_floor, x)
    expect_tiers_give(ceiling(x), simd_ceiling, x)
    expect_tiers_give(trunc(x), simd_trunc, x)
    expect_tiers_give(round(x), simd_round, x)
  }
})

test_that("integer input gives a double result", {
  x <- c(-3L, 0L, NA, 7L)
  expect_tiers_give(floor(x), simd_floor, x)
  expect_tiers_give(round(x), simd_round, x)
  expect_tiers_give(as.double(c(TRUE, NA)), simd_trunc, c(TRUE, NA))
  expect_tiers_give(round(x, -1), simd_round, x, -1)
})

test_that("round with digits uses base R's algorithm", {
  expect_identical(simd_round(2.675, 2), 2.67)
  expect_identical(simd_round(1234.5, -1), 1230)
  x <- with_seed(81L, stats::runif(1e5, -1e4, 1e4))
  for (d in c(2, -1, 5L)) expect_identical(simd_round(x, d), round(x, d), info = d)
  x <- c(NA, NaN, Inf, -0.125, 0)
  expect_true(same_values(simd_round(x, 2), round(x, 2)))
  expect_true(same_values(simd_round(c(1.5, NaN), NA), round(c(1.5, NaN), NA)))
  expect_true(same_values(simd_round(1.5, NaN), round(1.5, NaN)))
  expect_identical(simd_round(1.2345, 2.6), round(1.2345, 2.6))
})

test_that("digits must be a single number", {
  expect_error(simd_round(1, 1:2), "'digits' must be a single number")
  expect_error(simd_round(1, "2"), "'digits' must be a single number")
  expect_error(simd_round(1.26, structure(5e-324, class = "integer64")), "'digits' must be a single number")
  expect_identical(simd_round(2.55, TRUE), round(2.55, TRUE))
  expect_identical(simd_round(2.567, TRUE), round(2.567, TRUE))
  expect_identical(simd_round(2.55, FALSE), round(2.55, FALSE))
  expect_identical(simd_round(c(1.5, NaN), NA), round(c(1.5, NaN), NA))
  expect_identical(simd_round(1.26, simd_vec(1)), round(1.26, 1))
  expect_identical(as.double(simd_round(simd_vec(1.26), TRUE)), 1.3)
  expect_error(simd_round(as.raw(1)), "non-numeric argument to mathematical function")
  expect_error(simd_round(as.raw(1), 2), "non-numeric argument to mathematical function")
  expect_error(simd_floor(1i), "simd_floor() does not support 'x' of type complex", fixed = TRUE)
})

test_that("chunk boundaries and ALTREP inputs give the oracle's result", {
  n <- 2^20 + 7
  x <- seq_len(n) / 4
  expect_simd_identical(simd_round, x)
  expect_identical(simd_floor(seq_len(n)), as.double(seq_len(n)))
  expect_identical(simd_round(seq_len(5000) / 8, 1), round(seq_len(5000) / 8, 1))
})
