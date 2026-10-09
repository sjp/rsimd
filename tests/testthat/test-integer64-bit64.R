# integer64 results compared with bit64's own methods, including the texts
# of the warnings that the documentation quotes. Runs only with bit64.

skip_if_not_installed("bit64")
suppressPackageStartupMessages(library(bit64))

# Random integer64 values in (-2^40, 2^40) with NA: sums of many of them,
# and their products with small values, stay in range, so bit64's
# step-by-step overflow checks agree with rsimd's exact ones.
rand_b64 <- function(n, seed, na_frac = 0.05, bound = 2^40) {
  with_seed(seed, {
    x <- as.integer64(round(stats::runif(n, -bound, bound)))
    x[stats::runif(n) < na_frac] <- NA
    x
  })
}

# The value and warnings of f(...) on every tier match g(...) (bit64).
expect_like_bit64 <- function(g, f, ...) {
  want <- value_and_warnings(g(...))
  res <- with_each_tier(function() value_and_warnings(f(...)))
  for (tier in names(res)) {
    expect_identical(res[[tier]]$value, want$value, info = paste("tier", tier))
    expect_identical(res[[tier]]$warnings, want$warnings, info = paste("tier", tier))
  }
}

test_that("reductions match bit64", {
  for (n in c(0, 1, 7, 33, 1000, 2^16 + 3)) {
    x <- rand_b64(n, seed = n + 1)
    expect_like_bit64(sum, simd_sum, x)
    expect_like_bit64(function(x, na.rm) sum(x, na.rm = na.rm), simd_sum, x, na.rm = TRUE)
    expect_like_bit64(min, simd_min, x)
    expect_like_bit64(function(x, na.rm) max(x, na.rm = na.rm), simd_max, x, na.rm = TRUE)
    expect_like_bit64(function(x, na.rm) range(x, na.rm = na.rm), simd_range, x, na.rm = TRUE)
    expect_like_bit64(anyNA, simd_any_na, x)
    expect_like_bit64(function(x) which(is.na(x)), simd_which_na, x)
    expect_like_bit64(function(x) as.logical(is.na(x)), simd_is_na, x)
    expect_like_bit64(function(x, na.rm) any(x != 0, na.rm = na.rm), simd_any, x, na.rm = TRUE)
    y <- x[!is.na(x)]
    # bit64 has no which.min method; the first position of the minimum.
    expect_like_bit64(
      function(x) if (length(y)) which(x == min(y))[1L] else integer(0),
      simd_which_min, x
    )
    expect_like_bit64(
      function(x) if (length(y)) which(x == max(y))[1L] else integer(0),
      simd_which_max, x
    )
    small <- rand_b64(n, seed = n + 2, na_frac = 0.001, bound = 2^20)
    expect_like_bit64(cumsum, simd_cumsum, small)
    expect_like_bit64(cummin, simd_cummin, small)
    expect_like_bit64(cummax, simd_cummax, small)
  }
  x <- rand_b64(1e5, seed = 3, na_frac = 0)
  expect_like_bit64(sum, simd_sum, x)
})

test_that("simd_sum of 10^7 integer64 values matches bit64", {
  skip_unless_extended()
  x <- rand_b64(1e7, seed = 4, na_frac = 0)
  expect_like_bit64(sum, simd_sum, x)
})

test_that("overflow, division by zero and empty input warn as bit64 does", {
  mx <- lim.integer64()[2L]
  one <- as.integer64(1)
  expect_like_bit64(sum, simd_sum, c(mx, one))
  expect_like_bit64(function(x, y) x + y, simd_add, mx, one)
  expect_like_bit64(function(x, y) x - y, simd_sub, -mx, one)
  expect_like_bit64(function(x, y) x * y, simd_mul, as.integer64(3037000500), as.integer64(3037000500))
  expect_like_bit64(function(x, y) x * y, simd_mul, as.integer64(3037000499), as.integer64(3037000499))
  expect_like_bit64(cumsum, simd_cumsum, c(mx, one, one))
  expect_like_bit64(function(x, y) x %/% y, simd_idiv, as.integer64(c(7, 8)), as.integer64(c(0, 3)))
  expect_like_bit64(function(x, y) x %% y, simd_mod, as.integer64(c(7, -8)), as.integer64(c(0, 3)))
  expect_like_bit64(min, simd_min, integer64())
  expect_like_bit64(max, simd_max, integer64())
  expect_like_bit64(range, simd_range, integer64())
  expect_like_bit64(function(x, na.rm) min(x, na.rm = na.rm), simd_min, as.integer64(NA), na.rm = TRUE)
  # Where the step-by-step check differs: bit64 overflows, rsimd does not.
  expect_warning(b <- sum(c(mx, one, -one)), "NAs produced by integer64 overflow")
  expect_true(is.na(b))
  expect_identical(simd_sum(c(mx, one, -one)), mx)
})

test_that("elementwise arithmetic matches bit64", {
  ops <- list(
    list(`+`, simd_add), list(`-`, simd_sub), list(`*`, simd_mul), list(`%/%`, simd_idiv),
    list(`%%`, simd_mod)
  )
  for (n in c(1, 5, 64, 1001)) {
    x <- rand_b64(n, seed = n + 10, bound = 2^31)
    y <- rand_b64(n, seed = n + 11, bound = 2^31)
    y[!is.na(y) & y == 0] <- as.integer64(3)
    for (op in ops) {
      expect_like_bit64(op[[1L]], op[[2L]], x, y)
      expect_like_bit64(op[[1L]], op[[2L]], x, y[1L])
      expect_like_bit64(op[[1L]], op[[2L]], x, 3L)
    }
    big <- rand_b64(n, seed = n + 12, bound = 2^62)
    for (op in ops[1:3]) expect_like_bit64(op[[1L]], op[[2L]], big, big)
    expect_like_bit64(function(x, y) as.double(x) / as.double(y), simd_div, x, y)
    expect_like_bit64(abs, simd_abs, x)
    expect_like_bit64(sign, simd_sign, x)
    expect_like_bit64(function(x) -x, simd_neg, x)
  }
})

# Operands for %/% and %%: blocks of the kernel (runs of 256) with all
# operands within +-2^51 are divided in double lanes, others in scalar code,
# and after eight of those in a row the rest is.
intdiv_cases <- function(n, seed) {
  p51 <- as.integer64(2)^51L
  edges <- c(p51 - 1L, -p51, -p51 + 1L, p51, -p51 - 1L, lim.integer64(), as.integer64(c(0, 1, -1, 7, -7, NA)))
  small <- rand_b64(n, seed = seed, bound = 2^51 - 1)
  big <- rand_b64(n, seed = seed + 1, bound = 2^62)
  run <- (seq_len(n) - 1) %/% 256
  mix <- function(pick) {
    x <- small
    x[pick] <- big[pick]
    x
  }
  edged <- small
  edged[seq(1, n, by = 37)] <- rep_len(edges, length(seq(1, n, by = 37)))
  list(
    small = small, edged = edged, big = big, big_first = mix(run < 9), alternating = mix(run %% 2 == 0)
  )
}
intdiv_divisors <- function() {
  p51 <- as.integer64(2)^51L
  c(
    as.integer64(c(3, -7, 1000, 1, -1, 2, 6700417, -1000000007, 2^40, 0, NA)), p51 - 1L, -p51, p51,
    as.integer64(2)^62L, lim.integer64()
  )
}

test_that("%/% and %% match bit64 in and out of the double-lane range", {
  n <- 13 * 256 + 77
  y <- rand_b64(n, seed = 32, bound = 1e6)
  y[seq(5, n, by = 101)] <- as.integer64(0)
  yi <- as.integer(y)
  ops <- list(list(`%/%`, simd_idiv), list(`%%`, simd_mod))
  for (x in intdiv_cases(n, seed = 30)) {
    for (op in ops) {
      expect_like_bit64(op[[1L]], op[[2L]], x, y)
      expect_like_bit64(op[[1L]], op[[2L]], x, yi)
      expect_like_bit64(op[[1L]], op[[2L]], y, x)
      expect_like_bit64(op[[1L]], op[[2L]], x[1L], y)
      expect_like_bit64(op[[1L]], op[[2L]], x, -1000L)
      for (d in as.list(intdiv_divisors())) expect_like_bit64(op[[1L]], op[[2L]], x, d)
    }
  }
})

test_that("%/% and %% match bit64 on long random input (extended)", {
  skip_unless_extended()
  n <- 2^20 + 1000
  ops <- list(list(`%/%`, simd_idiv), list(`%%`, simd_mod))
  for (seed in 1:3) {
    y <- rand_b64(n, seed = 40 + seed, bound = 2^(10 * seed))
    y[y == 0] <- as.integer64(-3)
    for (x in intdiv_cases(n, seed = 50 + seed)) {
      for (op in ops) {
        expect_like_bit64(op[[1L]], op[[2L]], x, y)
        d <- with_seed(seed, rand_b64(20, seed = 60 + seed, na_frac = 0, bound = 2^(20 * seed)))
        for (k in seq_along(d)) expect_like_bit64(op[[1L]], op[[2L]], x, d[k])
      }
    }
  }
})

test_that("comparisons and predicates match bit64", {
  x <- rand_b64(500, seed = 20, bound = 50)
  y <- rand_b64(500, seed = 21, bound = 50)
  cmp <- list(
    list(`==`, simd_eq), list(`!=`, simd_ne), list(`<`, simd_lt), list(`<=`, simd_le),
    list(`>`, simd_gt), list(`>=`, simd_ge)
  )
  for (op in cmp) {
    expect_like_bit64(op[[1L]], op[[2L]], x, y)
    expect_like_bit64(op[[1L]], op[[2L]], x, 3L)
  }
  expect_like_bit64(function(x) as.logical(!is.na(x) & x < 0), simd_is_negative, x)
  expect_like_bit64(function(x) as.logical(!is.na(x) & x == 0), simd_is_zero, x)
  expect_like_bit64(function(x) as.logical(!is.na(x) & x %% 2L == 0L), simd_is_even, x)
  expect_like_bit64(function(x) as.logical(!is.na(x) & x %% 2L == 1L), simd_is_odd, x)
  expect_like_bit64(function(x) as.logical(!is.na(x) & x != 0L), simd_is_normal, x)
  p2 <- as.integer64(c(1, 2, 3, 4, 0, -4, NA, 1024, 1023)) * as.integer64(2)^40L
  expect_like_bit64(function(x) c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE),
    simd_is_pow2, p2)
})

test_that("Hamming distances match bit64's sum(x != y)", {
  x <- rand_b64(1000, seed = 22, bound = 3)
  y <- rand_b64(1000, seed = 23, bound = 3)
  for (rm in c(FALSE, TRUE)) {
    expect_like_bit64(function(x, y) as.double(sum(x != y, na.rm = rm)),
      function(x, y) simd_hamming(x, y, na.rm = rm), x, y)
    expect_like_bit64(function(x, y) as.double(sum(x != y, na.rm = rm)),
      function(x, y) simd_hamming(x, y, na.rm = rm), x, 2L)
  }
})

test_that("conversions match bit64", {
  x <- c(rand_b64(200, seed = 30, bound = 2^62), as.integer64(c(NA, 0, 1, -1)))
  expect_like_bit64(as.double, simd_as_double, x)
  expect_like_bit64(as.double, simd_as_double, rand_b64(50, seed = 31, bound = 2^52))
  expect_like_bit64(as.integer, simd_as_integer, x)
  expect_like_bit64(as.integer, simd_as_integer, rand_b64(50, seed = 32, bound = 2^31))
  expect_like_bit64(as.logical, simd_as_logical, x)
  d <- c(with_seed(33, stats::runif(300, -2^63, 2^63)), 2^63, -2^63, 9.5, -9.5, NA, 0)
  expect_like_bit64(as.integer64, simd_as_integer64, d)
  expect_like_bit64(as.integer64, simd_as_integer64, c(2^63, -2^63, 9.5, NaN))
  expect_like_bit64(as.integer64, simd_as_integer64, c(1L, NA, -7L))
  expect_like_bit64(as.integer64, simd_as_integer64, c(TRUE, NA, FALSE))
})

test_that("integer64 mixed with double converts the integer64 operand", {
  # bit64 converts the double to integer64 instead; rsimd keeps the double.
  expect_identical(as.integer64(1) + 1.5, as.integer64(2))
  expect_warning(r <- simd_add(as.integer64(1), 1.5), "integer64 coerced to double", fixed = TRUE)
  expect_identical(r, 2.5)
})

test_that("integer64 results print and combine with bit64", {
  x <- simd_add(as.integer64(c(1, NA)), 1L)
  expect_s3_class(x, "integer64")
  expect_identical(as.character(x), c("2", NA))
  expect_identical(simd_as_integer64("1" == "1"), as.integer64(1))
})
