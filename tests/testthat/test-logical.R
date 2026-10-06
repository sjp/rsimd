# Three-valued logic: identical to base R's &, |, xor and ! on every tier.

logic_ops <- list(simd_and = `&`, simd_or = `|`, simd_xor = xor)

test_that("the full truth tables match base R", {
  p <- all_pairs(c(TRUE, FALSE, NA))
  for (name in names(logic_ops)) {
    expect_tiers_give(logic_ops[[name]](p$x, p$y), get(name), p$x, p$y)
  }
  expect_tiers_give(!c(TRUE, FALSE, NA), simd_not, c(TRUE, FALSE, NA))
  expect_identical(simd_and(NA, FALSE), FALSE)
  expect_identical(simd_and(NA, TRUE), NA)
  expect_identical(simd_or(NA, TRUE), TRUE)
  expect_identical(simd_or(NA, FALSE), NA)
  expect_identical(simd_xor(NA, TRUE), NA)
  expect_identical(simd_not(NA), NA)
})

test_that("integer and double operands are read as logical values", {
  vals <- list(
    logical = c(TRUE, FALSE, NA),
    integer = c(0L, 2L, -1L, NA, .Machine$integer.max),
    double = c(0, -0, 2.5, NA, NaN, Inf, 5e-324)
  )
  for (a in vals) {
    for (b in vals) {
      p <- all_pairs(a, b)
      for (name in names(logic_ops)) {
        expect_tiers_give(logic_ops[[name]](p$x, p$y), get(name), p$x, p$y)
      }
    }
    expect_tiers_give(!a, simd_not, a)
  }
  expect_identical(class(simd_and(1L, 2L)), "logical")
})

test_that("a length-1 operand is broadcast on either side", {
  x <- c(TRUE, FALSE, NA, TRUE, FALSE, NA, TRUE)
  for (name in names(logic_ops)) {
    f <- get(name)
    for (s in list(TRUE, FALSE, NA, 0L, NaN)) {
      expect_tiers_give(logic_ops[[name]](x, s), f, x, s)
      expect_tiers_give(logic_ops[[name]](s, x), f, s, x)
    }
  }
  expect_error(simd_and(c(TRUE, FALSE), c(TRUE, FALSE, NA)), "must be equal or one of them")
})

test_that("raw operands combine bytewise, for all pairs of bytes", {
  v <- as.raw(0:255)
  p <- all_pairs(v)
  for (name in names(logic_ops)) {
    expect_tiers_give(logic_ops[[name]](p$x, p$y), get(name), p$x, p$y)
  }
  expect_tiers_give(!v, simd_not, v)
  expect_identical(simd_and(as.raw(0xF0), as.raw(0x3C)), as.raw(0x30))
  expect_identical(simd_or(as.raw(0xF0), as.raw(0x3C)), as.raw(0xfc))
  expect_identical(simd_xor(as.raw(0xF0), as.raw(0x3C)), as.raw(0xcc))
  expect_identical(simd_not(as.raw(0xF0)), as.raw(0x0f))
})

test_that("mixing raw with another type errors with base R's message", {
  msg <- tryCatch(TRUE & as.raw(1), error = conditionMessage)
  expect_error(simd_and(TRUE, as.raw(1)), msg, fixed = TRUE)
  expect_error(simd_or(as.raw(1), 1L), msg, fixed = TRUE)
  expect_error(simd_xor(as.raw(1), 2.5), msg, fixed = TRUE)
  expect_error(simd_and(1i, TRUE), "simd_and() does not support 'x' of type complex", fixed = TRUE)
})

test_that("logic matches base R on random vectors of edge lengths", {
  for (n in sweep_lengths()) {
    x <- rand_vec("logical", n, seed = n + 11L)
    y <- rand_vec("logical", n, seed = n + 12L)
    d <- rand_vec("double", n, seed = n + 13L)
    for (name in names(logic_ops)) {
      expect_tiers_give(logic_ops[[name]](x, y), get(name), x, y)
      expect_tiers_give(logic_ops[[name]](d, y), get(name), d, y)
    }
    expect_tiers_give(!d, simd_not, d)
  }
  expect_identical(simd_and(c(a = TRUE), TRUE), TRUE)
})
