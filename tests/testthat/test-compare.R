# Elementwise comparisons: identical() to base R's operators for every type
# pair, NA positions included, on every tier.

cmp_ops <- list(
  simd_eq = `==`, simd_ne = `!=`, simd_lt = `<`, simd_le = `<=`,
  simd_gt = `>`, simd_ge = `>=`
)

cmp_values <- list(
  double = c(0, -0, 1, -1, 1.5, 200, 255.5, 2147483647, -2147483647, 1e300, Inf, -Inf, NA, NaN),
  integer = c(0L, 1L, -1L, 200L, 255L, .Machine$integer.max, -.Machine$integer.max, NA),
  logical = c(TRUE, FALSE, NA),
  raw = as.raw(c(0, 1, 200, 255))
)

test_that("every type pair matches base R, NA and NaN giving NA", {
  for (tx in names(cmp_values)) {
    for (ty in names(cmp_values)) {
      p <- all_pairs(cmp_values[[tx]], cmp_values[[ty]])
      for (name in names(cmp_ops)) {
        expected <- cmp_ops[[name]](p$x, p$y)
        expect_tiers_give(expected, get(name), p$x, p$y)
      }
    }
  }
  expect_identical(simd_eq(-0, 0), TRUE)
  expect_identical(simd_eq(as.raw(1), 1L), TRUE)
  expect_identical(simd_eq(TRUE, as.raw(200)), TRUE)
  expect_identical(simd_lt(as.raw(200), as.raw(3)), FALSE)
  expect_identical(simd_ne(NaN, NaN), NA)
})

test_that("a length-1 operand is broadcast on either side", {
  x <- c(-2, NA, 0, 3.5, NaN, 7)
  for (name in names(cmp_ops)) {
    f <- get(name)
    expect_tiers_give(cmp_ops[[name]](x, 0), f, x, 0)
    expect_tiers_give(cmp_ops[[name]](0, x), f, 0, x)
    expect_tiers_give(cmp_ops[[name]](3L, 1:9), f, 3L, 1:9)
    expect_tiers_give(cmp_ops[[name]](NA, x), f, NA, x)
    expect_tiers_give(cmp_ops[[name]](as.raw(1:9), as.raw(4)), f, as.raw(1:9), as.raw(4))
  }
  expect_error(simd_eq(1:3, 1:2), "lengths of 'x' (3) and 'y' (2)", fixed = TRUE)
  expect_identical(simd_lt(integer(0), 1L), logical(0))
})

test_that("comparisons match base R on random vectors of edge lengths", {
  for (n in edge_lengths()) {
    xd <- rand_vec("double", n, seed = n + 5L)
    yd <- rand_vec("double", n, seed = n + 6L)
    xi <- rand_vec("integer", n, seed = n + 7L)
    yi <- rand_vec("integer", n, seed = n + 8L)
    xd[seq_len(n) %% 3 == 0] <- yd[seq_len(n) %% 3 == 0] # ties
    for (name in names(cmp_ops)) {
      f <- get(name)
      expect_tiers_give(cmp_ops[[name]](xd, yd), f, xd, yd)
      expect_tiers_give(cmp_ops[[name]](xi, yi), f, xi, yi)
      expect_tiers_give(cmp_ops[[name]](xi, yd), f, xi, yd)
    }
  }
})

test_that("attributes are dropped and unsupported types rejected", {
  expect_identical(simd_eq(c(a = 1, b = 2), 1), c(TRUE, FALSE))
  expect_identical(simd_eq(1i, 1), FALSE)
  expect_error(simd_lt(1, 1i), "invalid comparison with complex values", fixed = TRUE)
  expect_error(simd_gt("a", 1), "must be an atomic vector")
})
