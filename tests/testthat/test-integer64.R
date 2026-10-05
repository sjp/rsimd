# integer64 support, tested without bit64: inputs are built from their bit
# patterns (helper-int64.R), every tier is compared with the none tier bit
# for bit, and with exact expectations.

ovf_msg <- "NAs produced by integer64 overflow"
lens <- c(0, 1, 2, 3, 7, 8, 9, 31, 32, 33, 67, 130, 1000, 4097)

test_that("the integer64 test helpers are exact", {
  x <- i64(c(0, 1, -1, 2^53 + 2, -2^62, NA))
  expect_identical(i64_str(x), c("0", "1", "-1", "9007199254740994", "-4611686018427387904", NA))
  s <- c("9223372036854775807", "-9223372036854775807", "9007199254740993", "-4294967297")
  expect_identical(i64_str(i64_dec(s)), s)
  expect_i64(i64_dec("-1"), i64(-1))
  expect_i64(i64_from_bits(i64_bits(i64_dec(s))), i64_dec(s))
  expect_identical(i64_value(i64_dec("9007199254740993")), 2^53)
})

# ---- Reductions ---------------------------------------------------------------

test_that("simd_sum of integer64 is exact and integer64", {
  expect_tiers_i64(i64(0), simd_sum, i64(numeric(0)))
  expect_tiers_i64(na64, simd_sum, i64(c(1, 2, NA, 4)))
  expect_tiers_i64(i64(7), simd_sum, i64(c(1, 2, NA, 4)), na.rm = TRUE)
  expect_tiers_i64(i64(0), simd_sum, i64(c(NA, NA)), na.rm = TRUE)
  for (n in lens[lens > 0]) {
    v <- round(with_seed(n, stats::runif(n, -2^40, 2^40)))
    expect_tiers_i64(i64(sum(v)), simd_sum, i64(v))
  }
})

test_that("simd_sum overflows only when the exact total leaves int64", {
  one <- i64(1)
  # The running total leaves int64 but the total does not.
  expect_tiers_i64(max64, simd_sum, i64c(max64, one, i64(-1)))
  lo <- i64_dec("-9223372036854775807")
  expect_tiers_i64(i64(0), simd_sum, i64c(max64, max64, lo, lo))
  # Many wraps in a random order.
  v <- c(rep(c("9223372036854775807", "-9223372036854775807"), 500), "5")
  v <- v[with_seed(3L, sample(length(v)))]
  expect_tiers_i64(i64(5), simd_sum, i64_dec(v))
  # Totals outside int64, and INT64_MIN (the NA pattern): NA, one warning.
  expect_tiers_i64(na64, simd_sum, i64c(max64, one), msgs = ovf_msg)
  expect_tiers_i64(na64, simd_sum, i64_dec(c("-9223372036854775807", "-1")), msgs = ovf_msg)
  expect_tiers_i64(na64, simd_sum, i64_dec(rep("9223372036854775807", 5)), msgs = ovf_msg)
  expect_tiers_i64(
    i64_dec("-9223372036854775807"), simd_sum,
    i64_dec(c("-9223372036854775807", "-1", "1"))
  )
  # A missing value wins over an overflow, without a warning.
  expect_tiers_i64(na64, simd_sum, i64c(max64, one, na64))
})

test_that("simd_sum of integer64 agrees across tiers on random input", {
  for (n in lens) {
    for (kind in c("mix", "full", "small")) {
      x <- rand_i64(n, kind, seed = n + 7L)
      expect_simd_identical(function(x) suppressWarnings(simd_sum(x)), x)
      expect_simd_identical(function(x) suppressWarnings(simd_sum(x, na.rm = TRUE)), x)
    }
    x <- rand_i64(n, "small", na_frac = 0, seed = n)
    expect_simd_identical(simd_sum, x, na_check = FALSE)
  }
  # Several chunks.
  v <- round(with_seed(11L, stats::runif(2^20 + 7, -2^31, 2^31)))
  expect_tiers_i64(i64(sum(v)), simd_sum, i64(v))
  x <- rand_i64(2^20 + 7, "full", seed = 12L)
  expect_simd_identical(function(x) suppressWarnings(simd_sum(x, na.rm = TRUE)), x)
})

test_that("simd_sum of 10^7 integer64 values is exact", {
  skip_unless_extended()
  v <- round(with_seed(13L, stats::runif(1e7, -2^29, 2^29)))
  x <- i64(v)
  expect_tiers_i64(i64(sum(v)), simd_sum, x)
})

test_that("min, max and range of integer64", {
  x <- i64(c(3, NA, 1, -7, 1))
  expect_tiers_i64(na64, simd_min, x)
  expect_tiers_i64(i64(-7), simd_min, x, na.rm = TRUE)
  expect_tiers_i64(i64(3), simd_max, x, na.rm = TRUE)
  expect_tiers_i64(i64(c(-7, 3)), simd_range, x, na.rm = TRUE)
  expect_tiers_i64(i64(c(NA, NA)), simd_range, x)
  ends <- i64_dec(c("9223372036854775807", "-9223372036854775807", "0"))
  expect_tiers_i64(i64_at(ends, 2:1), simd_range, ends)
  # Empty input, as bit64: the extreme values with a warning.
  lo <- i64_dec("-9223372036854775807")
  hi_msg <- "no non-NA value, returning the highest possible integer64 value +9223372036854775807"
  lo_msg <- "no non-NA value, returning the lowest possible integer64 value -9223372036854775807"
  expect_tiers_i64(max64, simd_min, i64(numeric(0)), msgs = hi_msg)
  expect_tiers_i64(lo, simd_max, i64(numeric(0)), msgs = lo_msg)
  expect_tiers_i64(max64, simd_min, i64(c(NA, NA)), na.rm = TRUE, msgs = hi_msg)
  expect_tiers_i64(i64c(max64, lo), simd_range, i64(NA),
    na.rm = TRUE,
    msgs = "no non-NA value, returning c(+9223372036854775807, -9223372036854775807)"
  )
  for (n in lens) {
    x <- rand_i64(n, seed = n + 1L)
    for (f in list(simd_min, simd_max, simd_range)) {
      expect_simd_identical(function(x, ...) suppressWarnings(f(x, ...)), x)
      expect_simd_identical(function(x, ...) suppressWarnings(f(x, ...)), x, na.rm = TRUE)
    }
  }
})

test_that("which_min and which_max of integer64 ignore NA", {
  x <- i64(c(NA, 3, 1, -7, -7, 9, 9))
  expect_tiers_give(4L, simd_which_min, x)
  expect_tiers_give(6L, simd_which_max, x)
  expect_tiers_give(integer(0), simd_which_min, i64(c(NA, NA)))
  expect_tiers_give(integer(0), simd_which_max, i64(numeric(0)))
  for (n in lens) {
    x <- rand_i64(n, seed = n + 2L)
    expect_simd_identical(simd_which_min, x)
    expect_simd_identical(simd_which_max, x)
  }
})

test_that("missing values of integer64", {
  x <- i64(c(1, NA, 0, NA, -5))
  expect_tiers_give(TRUE, simd_any_na, x)
  expect_tiers_give(2, simd_count_na, x)
  expect_tiers_give(c(2L, 4L), simd_which_na, x)
  expect_tiers_give(c(FALSE, TRUE, FALSE, TRUE, FALSE), simd_is_na, x)
  expect_tiers_give(FALSE, simd_any_na, max64)
  # INT64_MIN + 1 is a number.
  expect_tiers_give(FALSE, simd_is_na, i64_dec("-9223372036854775807"))
  for (n in lens) {
    x <- rand_i64(n, seed = n + 3L, na_frac = 0.2)
    for (f in list(
      simd_any_na, simd_count_na, simd_which_na, simd_is_na, simd_is_na_any,
      simd_is_na_all
    )) {
      expect_simd_identical(f, x)
    }
  }
})

test_that("predicates of integer64", {
  x <- i64(c(-3, 0, NA, 5))
  expect_tiers_give(c(FALSE, FALSE, FALSE, FALSE), simd_is_nan, x)
  expect_tiers_give(c(TRUE, TRUE, FALSE, TRUE), simd_is_finite, x)
  expect_tiers_give(c(FALSE, FALSE, FALSE, FALSE), simd_is_infinite, x)
  expect_tiers_give(c(TRUE, FALSE, FALSE, FALSE), simd_is_negative, x)
  expect_tiers_give(c(FALSE, TRUE, FALSE, FALSE), simd_is_zero, x)
  expect_tiers_give(TRUE, simd_is_negative_any, x)
  expect_tiers_give(FALSE, simd_is_finite_all, x)
  expect_tiers_give(TRUE, simd_is_finite_all, i64(numeric(0)))
  expect_tiers_give(FALSE, simd_is_zero_any, i64(numeric(0)))
  for (n in lens) {
    x <- rand_i64(n, seed = n + 4L)
    for (f in list(
      simd_is_finite, simd_is_negative, simd_is_zero, simd_is_negative_any,
      simd_is_negative_all, simd_is_zero_any, simd_is_finite_all
    )) {
      expect_simd_identical(f, x)
    }
  }
})

test_that("number classes of integer64 follow their values; NA is in none", {
  is_na <- function(x) i64_str(x) %in% NA
  ref <- function(x) {
    b <- i64_bits(x)
    na <- is_na(x)
    nz <- colSums(b) > 0
    list(
      normal = nz & !na, subnormal = rep(FALSE, length(x)), whole = !na,
      even = b[1, ] == 0 & !na, odd = b[1, ] == 1,
      pow2 = colSums(b) == 1 & b[64, ] == 0
    )
  }
  x <- i64_dec(c(
    "0", "1", "-1", "2", "-2", NA, "9223372036854775807", "-9223372036854775807",
    "4611686018427387904", "4294967296", "4294967297", "6", "7", "-4611686018427387904"
  ))
  r <- ref(x)
  expect_identical(r$pow2, c(FALSE, TRUE, FALSE, TRUE, rep(FALSE, 4), TRUE, TRUE, rep(FALSE, 4)))
  for (op in names(r)) {
    f <- get(paste0("simd_is_", op))
    expect_tiers_give(r[[op]], f, x)
    expect_tiers_give(any(r[[op]]), get(paste0("simd_is_", op, "_any")), x)
    expect_tiers_give(all(r[[op]]), get(paste0("simd_is_", op, "_all")), x)
  }
  for (n in lens) {
    x <- rand_i64(n, seed = n + 11L)
    r <- ref(x)
    for (op in names(r)) {
      expect_tiers_give(r[[op]], get(paste0("simd_is_", op)), x)
      expect_simd_identical(get(paste0("simd_is_", op, "_any")), x)
      expect_simd_identical(get(paste0("simd_is_", op, "_all")), x)
    }
  }
})

test_that("Hamming distances of integer64", {
  # Pairs differ when their bit patterns do; NA pairs are missing.
  # The length-1 broadcast rule: a length-0 operand gives no pairs.
  pairs <- function(x, y) if (min(length(x), length(y)) == 0) 0 else max(length(x), length(y))
  ref <- function(x, y, na.rm = FALSE) {
    n <- pairs(x, y)
    x <- i64_at(x, rep_len(seq_along(x), n))
    y <- i64_at(y, rep_len(seq_along(y), n))
    miss <- i64_str(x) %in% NA | i64_str(y) %in% NA
    if (any(miss) && !na.rm) {
      return(NA_real_)
    }
    as.double(sum((colSums(i64_bits(x) != i64_bits(y)) > 0)[!miss]))
  }
  bits <- function(x, y) {
    n <- pairs(x, y)
    as.double(sum(i64_bits(i64_at(x, rep_len(seq_along(x), n))) !=
      i64_bits(i64_at(y, rep_len(seq_along(y), n)))))
  }
  for (n in lens) {
    x <- rand_i64(n, kind = "small", seed = n + 12L)
    y <- i64_at(x, rev(seq_len(n)))
    if (n > 2) y <- i64c(i64_at(x, 1:2), i64_at(y, 3:n))
    for (rm in c(FALSE, TRUE)) {
      expect_tiers_give(ref(x, y, rm), simd_hamming, x, y, na.rm = rm)
      expect_tiers_give(ref(x, i64(1), rm), simd_hamming, x, i64(1), na.rm = rm)
    }
    expect_tiers_give(bits(x, y), simd_hamming_bits, x, y)
    expect_tiers_give(bits(x, max64), simd_hamming_bits, x, max64)
  }
  # integer and logical operands are compared as integer64, without a warning.
  x <- i64_dec(c("1", "4294967298", NA, "-3"))
  expect_tiers_give(1, simd_hamming, x, c(1L, 2L, 5L, -3L), na.rm = TRUE)
  expect_tiers_give(NA_real_, simd_hamming, x, c(1L, 2L, 5L, -3L))
  expect_tiers_give(1, simd_hamming, c(FALSE, NA, TRUE), i64(c(1, 1, 1)), na.rm = TRUE)
  expect_tiers_give(64, simd_hamming_bits, i64(-1), i64(0))
  expect_tiers_give(1, simd_hamming_bits, na64, i64(0))
  # With a double or complex operand they are converted to double first.
  expect_warning(
    expect_identical(simd_hamming(x, c(1, 4294967298, 0, -3.5), na.rm = TRUE), 1),
    "integer64 coerced to double"
  )
  expect_warning(
    expect_identical(simd_hamming(x, c(1, 2, 3, -3) + 0i, na.rm = TRUE), 1),
    "integer64 coerced to double"
  )
})

test_that("any and all of integer64 read non-zero as TRUE, without a warning", {
  x <- i64(c(0, 5, NA))
  expect_tier_warnings(character(0), simd_any, x)
  expect_tiers_give(TRUE, simd_any, x)
  expect_tiers_give(FALSE, simd_all, x)
  expect_tiers_give(NA, simd_any, i64(c(0, NA)))
  expect_tiers_give(FALSE, simd_any, i64(c(0, NA)), na.rm = TRUE)
  expect_tiers_give(NA, simd_all, i64(c(1, NA)))
  expect_tiers_give(TRUE, simd_all, i64(c(1, NA)), na.rm = TRUE)
  for (n in lens) {
    x <- rand_i64(n, seed = n + 5L, na_frac = 0.01)
    expect_simd_identical(simd_any, x)
    expect_simd_identical(simd_all, x, na.rm = TRUE)
    z <- i64(rep(0, n))
    expect_simd_identical(simd_any, z)
  }
})

test_that("cumsum, cummin and cummax of integer64", {
  x <- i64(c(5, -3, 7, NA, 2))
  expect_tiers_i64(i64(c(5, 2, 9, NA, NA)), simd_cumsum, x)
  expect_tiers_i64(i64(c(5, -3, -3, NA, NA)), simd_cummin, x)
  expect_tiers_i64(i64(c(5, 5, 7, NA, NA)), simd_cummax, x)
  expect_tiers_i64(i64(numeric(0)), simd_cumsum, i64(numeric(0)))
  expect_tiers_i64(i64c(max64, na64, na64), simd_cumsum, i64c(max64, i64(1), i64(-1)),
    msgs = ovf_msg
  )
  # The extreme values are reached without overflow.
  ends <- i64_dec(c("9223372036854775807", "-9223372036854775807"))
  expect_tiers_i64(ends, simd_cummin, ends)
  expect_tiers_i64(i64c(max64, max64), simd_cummax, ends)
  for (n in lens) {
    x <- rand_i64(n, "small", seed = n + 6L, na_frac = 0.001)
    for (f in list(simd_cumsum, simd_cummin, simd_cummax)) expect_simd_identical(f, x)
  }
})

# ---- Arithmetic --------------------------------------------------------------

test_that("checked integer64 arithmetic gives NA and one warning on overflow", {
  one <- i64(1)
  lo <- i64_dec("-9223372036854775807")
  expect_tiers_i64(i64(c(6, -2, NA)), simd_add, i64(c(5, -3, NA)), one)
  expect_tiers_i64(i64c(na64, i64(2)), simd_add, i64c(max64, one), one, msgs = ovf_msg)
  expect_tiers_i64(na64, simd_sub, lo, one, msgs = ovf_msg)
  expect_tiers_i64(i64_dec("-9223372036854775806"), simd_sub, lo, i64(-1))
  expect_tiers_i64(na64, simd_mul, i64(3037000500), i64(3037000500), msgs = ovf_msg)
  expect_tiers_i64(i64_dec("9223372030926249001"), simd_mul, i64(3037000499), i64(3037000499))
  expect_tiers_i64(i64_dec("-9223372030926249001"), simd_mul, i64(3037000499), i64(-3037000499))
  # Products that overflow only in the high half, or are INT64_MIN.
  expect_tiers_i64(na64, simd_mul, i64(2^32), i64(2^32), msgs = ovf_msg)
  expect_tiers_i64(na64, simd_mul, i64(-2^31), i64(2^32), msgs = ovf_msg)
  expect_tiers_i64(i64(2^62), simd_mul, i64(2^31), i64(2^31))
  expect_tiers_i64(i64_dec("-4611686014132420609"), simd_mul, i64(2^31 - 1), i64(-(2^31 - 1)))
  # An NA operand gives NA without a warning.
  expect_tiers_i64(i64(c(NA, NA)), simd_mul, i64(c(NA, 2^40)), i64(c(2^40, NA)))
  # Integer and logical operands are taken as integer64.
  expect_tiers_i64(i64(c(6, NA, 2)), simd_add, i64(c(5, 1, 1)), c(1L, NA, 1L))
  expect_tiers_i64(i64(c(2, 1)), simd_add, c(TRUE, FALSE), i64(1))
  expect_tiers_i64(i64(c(-(2^31 - 1), NA)), simd_mul, c(-.Machine$integer.max, NA), i64(1))
})

test_that("wrapping integer64 arithmetic wraps", {
  lo <- i64_dec("-9223372036854775807")
  expect_tiers_i64(i64_dec("-9223372036709301616"), simd_mul_wrap, i64(3037000500), i64(3037000500))
  expect_tiers_i64(lo, simd_add_wrap, max64, i64(2))
  expect_tiers_i64(max64, simd_sub_wrap, lo, i64(2))
  expect_tiers_i64(i64(0), simd_mul_wrap, i64(2^32), i64(2^32))
  expect_tiers_i64(i64(c(NA, NA)), simd_add_wrap, i64(c(NA, 1)), i64(c(1, NA)))
  expect_error(simd_add_wrap(i64(1), 1.5), "does not support 'y' of type double")
})

test_that("neg, abs and sign of integer64", {
  x <- i64c(i64(c(-5, 0, 7, NA)), max64, i64_dec("-9223372036854775807"))
  expect_tiers_i64(i64c(i64(c(5, 0, -7, NA)), i64_dec("-9223372036854775807"), max64), simd_neg, x)
  expect_tiers_i64(i64c(i64(c(5, 0, 7, NA)), max64, max64), simd_abs, x)
  expect_tiers_i64(i64(c(-1, 0, 1, NA, 1, -1)), simd_sign, x)
  expect_tiers_i64(simd_neg(x), simd_neg_wrap, x)
  expect_tiers_i64(simd_abs(x), simd_abs_wrap, x)
})

test_that("%/% and %% of integer64 floor, and a zero divisor warns", {
  x <- c(7, -7, 7, -7, 0, 2^31 - 5)
  y <- c(2, 2, -2, -2, 3, -17)
  expect_tiers_i64(i64(x %/% y), simd_idiv, i64(x), i64(y))
  expect_tiers_i64(i64(x %% y), simd_mod, i64(x), i64(y))
  big <- i64_dec(c("9223372036854775807", "-9223372036854775807"))
  expect_tiers_i64(i64_dec(c("922337203685477580", "-922337203685477581")), simd_idiv, big, i64(10))
  expect_tiers_i64(i64(c(7, 3)), simd_mod, big, i64(10))
  expect_tiers_i64(i64_dec(c("-9223372036854775807", "9223372036854775807")), simd_idiv, big, i64(-1))
  dz <- "NAs produced due to division by zero"
  expect_tiers_i64(i64(c(NA, 2)), simd_idiv, i64(c(5, 4)), i64(c(0, 2)), msgs = dz)
  expect_tiers_i64(i64(c(NA, NA)), simd_mod, i64(c(5, 4)), 0L, msgs = dz)
  expect_tiers_i64(i64(c(NA, NA)), simd_idiv, i64(c(NA, 4)), i64(c(0, NA)))
})

test_that("pmin, pmax and clamp of integer64", {
  x <- i64(c(1, NA, 5, NA))
  y <- i64(c(0, 3, NA, NA))
  expect_tiers_i64(i64(c(0, NA, NA, NA)), simd_pmin, x, y)
  expect_tiers_i64(i64(c(1, NA, NA, NA)), simd_pmax, x, y)
  expect_tiers_i64(i64(c(0, 3, 5, NA)), simd_pmin, x, y, na.rm = TRUE)
  expect_tiers_i64(i64(c(1, 3, 5, NA)), simd_pmax_num, x, y)
  expect_tiers_i64(i64(c(1, NA, 2, NA)), simd_pmin, x, 2L)
  v <- i64(c(-5, 0, 5, NA, 9))
  expect_tiers_i64(i64(c(-1, 0, 3, NA, 3)), simd_clamp, v, i64(-1), 3L)
  expect_error(simd_clamp(v, i64(4), i64(3)), "'lo' must not be greater than 'hi'")
})

test_that("mul_add and add_mul of integer64 count an overflow once", {
  expect_tiers_i64(i64(c(7, NA)), simd_mul_add, i64(c(3, 2^40)), i64(c(2, 2^40)), 1L,
    msgs = ovf_msg
  )
  expect_tiers_i64(i64(c(12, NA)), simd_add_mul, i64(c(3, NA)), 1L, i64(3))
  expect_tiers_i64(na64, simd_add_mul, max64, i64(1), i64(1), msgs = ovf_msg)
})

test_that("integer64 arithmetic agrees across tiers on random input", {
  ops <- list(
    simd_add, simd_sub, simd_mul, simd_add_wrap, simd_sub_wrap, simd_mul_wrap, simd_idiv,
    simd_mod, simd_pmin, simd_pmax, simd_pmin_num, simd_pmax_num
  )
  for (n in c(1, 3, 8, 33, 130, 4097)) {
    for (kind in c("mix", "small")) {
      x <- rand_i64(n, kind, seed = n + 20L)
      y <- rand_i64(n, kind, seed = n + 21L)
      yi <- with_seed(n, sample(c(-3:3, NA, .Machine$integer.max), n, replace = TRUE))
      for (f in ops) {
        g <- function(x, y) suppressWarnings(f(x, y))
        expect_simd_identical(g, x, y)
        expect_simd_identical(g, x, yi)
        expect_simd_identical(g, x, i64_at(y, 1))
        expect_simd_identical(g, i64_at(x, 1), y)
      }
      for (f in list(simd_neg, simd_abs, simd_sign)) expect_simd_identical(f, x)
      z <- rand_i64(n, kind, seed = n + 22L)
      for (f in list(simd_mul_add, simd_add_mul)) {
        expect_simd_identical(function(x, y, z) suppressWarnings(f(x, y, z)), x, y, z)
      }
      expect_simd_identical(function(x) simd_clamp(x, i64(-2^40), .Machine$integer.max), x)
    }
  }
  x <- rand_i64(2^20 + 7, "mix", seed = 30L)
  expect_simd_identical(function(x) suppressWarnings(simd_mul(x, x)), x)
  expect_simd_identical(function(x) suppressWarnings(simd_add(x, x, na_check = FALSE)), x)
})

test_that("integer64 mixed with double converts to double with a warning", {
  co <- "integer64 coerced to double"
  expect_warning(r <- simd_add(i64(c(1, NA)), 1.5), co, fixed = TRUE)
  expect_identical(r, c(2.5, NA))
  expect_warning(r <- simd_mul(2, i64(3)), co, fixed = TRUE)
  expect_identical(r, 6)
  expect_warning(r <- simd_pmin(i64(c(1, 5)), 2), co, fixed = TRUE)
  expect_identical(r, c(1, 2))
  expect_warning(r <- simd_clamp(i64(9), 0, i64(4)), co, fixed = TRUE)
  expect_identical(r, 4)
  # Precision is lost in the conversion (and only the coercion warns).
  expect_identical(warnings_of(r <- simd_eq(i64_dec("9007199254740993"), 2^53)), co)
  expect_identical(r, TRUE)
  expect_identical(warnings_of(simd_add(i64(2^60), 0.5)), co)
  # The warning names the function.
  w <- tryCatch(simd_sub(i64(1), 1), warning = function(w) w)
  expect_identical(conditionCall(w)[[1L]], as.name("simd_sub"))
  # simd_div gives double; it warns only with a double operand.
  expect_identical(warnings_of(r <- simd_div(i64(c(7, NA)), i64(2))), character(0))
  expect_identical(r, c(3.5, NA))
  expect_identical(simd_div(i64(7), 2L), 3.5)
  expect_warning(r <- simd_div(i64(7), 2), co, fixed = TRUE)
  expect_identical(r, 3.5)
})

test_that("other integer64 combinations and functions are rejected", {
  x <- i64(1)
  expect_error(simd_add(x, as.raw(1)), "non-numeric argument to binary operator", fixed = TRUE)
  expect_error(simd_add(x, 1i), "cannot combine complex and integer64", fixed = TRUE)
  for (f in c("simd_fma", "simd_lerp")) {
    expect_error(get(f)(x, 1, 1), paste0(f, "() does not support 'x' of type integer64"),
      fixed = TRUE
    )
  }
  for (f in c(
    "simd_sqrt", "simd_recip", "simd_floor", "simd_round", "simd_exp", "simd_var",
    "simd_sum_sq", "simd_cumprod"
  )) {
    expect_error(get(f)(x), paste0(f, "() does not support 'x' of type integer64"), fixed = TRUE)
  }
  expect_error(simd_copysign(x, 1), "does not support 'x' of type integer64", fixed = TRUE)
  expect_error(simd_prod(x), "simd_prod() does not support 'x' of type integer64 yet", fixed = TRUE)
  expect_error(simd_mean(x), "simd_mean() does not support 'x' of type integer64 yet", fixed = TRUE)
  expect_error(simd_and(x, TRUE), "simd_and() does not support 'x' of type integer64", fixed = TRUE)
  expect_error(simd_bit_and(x, 1), "does not support 'y' of type double", fixed = TRUE)
})

# ---- Comparisons ---------------------------------------------------------------

test_that("comparisons of integer64", {
  x <- i64_dec(c("9223372036854775807", "-5", "9007199254740993", "NA", "3"))
  y <- i64_dec(c("9223372036854775806", "-5", "9007199254740992", "1", "NA"))
  expect_tiers_give(c(FALSE, TRUE, FALSE, NA, NA), simd_eq, x, y)
  expect_tiers_give(c(TRUE, FALSE, TRUE, NA, NA), simd_ne, x, y)
  expect_tiers_give(c(FALSE, FALSE, FALSE, NA, NA), simd_lt, x, y)
  expect_tiers_give(c(FALSE, TRUE, FALSE, NA, NA), simd_le, x, y)
  expect_tiers_give(c(TRUE, FALSE, TRUE, NA, NA), simd_gt, x, y)
  expect_tiers_give(c(TRUE, TRUE, TRUE, NA, NA), simd_ge, x, y)
  expect_tiers_give(c(TRUE, FALSE, NA), simd_lt, i64(c(-2^40, 2^40, 0)), c(1L, 1L, NA))
  expect_tiers_give(c(FALSE, TRUE), simd_eq, as.raw(c(1, 3)), i64(3))
  for (n in lens) {
    x <- rand_i64(n, seed = n + 40L)
    y <- rand_i64(n, seed = n + 41L)
    yi <- with_seed(n, sample(c(-3:3, NA), n, replace = TRUE))
    for (f in list(simd_eq, simd_ne, simd_lt, simd_le, simd_gt, simd_ge)) {
      expect_simd_identical(f, x, y)
      expect_simd_identical(f, x, yi)
      expect_simd_identical(f, x, x)
      expect_simd_identical(f, i64_at(y, 1), x)
    }
  }
})

# ---- Bitwise ops ---------------------------------------------------------------

# Bitwise references on the bit matrix (least significant bit first).
ref_shift <- function(x, k, how) {
  b <- i64_bits(x)
  out <- apply(b, 2, function(v) {
    switch(how,
      shl = c(rep(0L, k), v[seq_len(64 - k)]),
      shr = c(v[seq_len(64 - k) + k], rep(0L, k)),
      sar = c(v[seq_len(64 - k) + k], rep(v[64], k)),
      rotl = c(v[seq_len(k) + 64 - k], v[seq_len(64 - k)]),
      rotr = c(v[seq_len(64 - k) + k], v[seq_len(k)])
    )
  })
  r <- i64_from_bits(out)
  na <- is.na(i64_str(x))
  structure(ifelse(na, unclass(na64), unclass(r)), class = "integer64")
}
ref_count <- function(x, how) {
  b <- i64_bits(x)
  r <- apply(b, 2, function(v) {
    switch(how,
      popcount = sum(v),
      lzcnt = if (any(v == 1L)) 64L - max(which(v == 1L)) else 64L,
      tzcnt = if (any(v == 1L)) min(which(v == 1L)) - 1L else 64L
    )
  })
  r <- as.integer(r)
  r[is.na(i64_str(x))] <- NA
  r
}
ref_bitop <- function(x, y, op) {
  r <- i64_from_bytes(op(i64_bytes(x), i64_bytes(y)))
  na <- is.na(i64_str(x)) | is.na(i64_str(y))
  structure(ifelse(na, unclass(na64), unclass(r)), class = "integer64")
}

test_that("bitwise ops of integer64 work on the bit patterns", {
  x <- rand_i64(40, "full", seed = 50L, na_frac = 0.1)
  y <- rand_i64(40, "full", seed = 51L, na_frac = 0.1)
  expect_tiers_i64(ref_bitop(x, y, `&`), simd_bit_and, x, y)
  expect_tiers_i64(ref_bitop(x, y, `|`), simd_bit_or, x, y)
  expect_tiers_i64(ref_bitop(x, y, xor), simd_bit_xor, x, y)
  expect_tiers_i64(ref_bitop(x, x, function(a, b) !a), simd_bit_not, x)
  yi <- c(-1L, 255L)
  expect_tiers_i64(
    ref_bitop(i64(c(2^40 + 300, 2^40 + 300)), i64(yi), `&`), simd_bit_and,
    i64(c(2^40 + 300, 2^40 + 300)), yi
  )
  for (k in c(0L, 1L, 31L, 32L, 33L, 62L, 63L)) {
    for (how in c("shl", "shr", "sar", "rotl", "rotr")) {
      expect_tiers_i64(ref_shift(x, k, how), get(paste0("simd_", how)), x, k)
    }
  }
  # Results with the NA pattern are NA; counts outside 0..63 give NA.
  expect_tiers_i64(na64, simd_shl, i64(1), 63)
  expect_tiers_i64(i64(c(NA, NA)), simd_shl, i64(c(1, 2)), 64)
  expect_tiers_i64(i64(c(NA, NA)), simd_shr, i64(c(1, 2)), -1)
  expect_tiers_i64(i64(c(NA, NA)), simd_sar, i64(c(1, 2)), NA)
  expect_tiers_i64(simd_rotl(x, 1), simd_rotl, x, 65)
  expect_tiers_i64(simd_rotr(x, 1), simd_rotl, x, -1)
  for (how in c("popcount", "lzcnt", "tzcnt")) {
    expect_tiers_give(ref_count(x, how), get(paste0("simd_", how)), x)
  }
  edges <- i64_dec(c("0", "-1", "1", "9223372036854775807", "-9223372036854775807", "NA"))
  expect_tiers_give(c(0L, 64L, 1L, 63L, 2L, NA), simd_popcount, edges)
  expect_tiers_give(c(64L, 0L, 63L, 1L, 0L, NA), simd_lzcnt, edges)
  expect_tiers_give(c(64L, 0L, 0L, 0L, 0L, NA), simd_tzcnt, edges)
  expect_tiers_give(as.double(sum(ref_count(x, "popcount"), na.rm = TRUE)),
    simd_popcount_total, x,
    na.rm = TRUE
  )
  expect_tiers_give(NA_real_, simd_popcount_total, x)
})

test_that("integer64 bitwise ops agree across tiers on random input", {
  for (n in lens) {
    x <- rand_i64(n, "full", seed = n + 60L)
    y <- rand_i64(n, "full", seed = n + 61L)
    for (f in list(simd_bit_and, simd_bit_or, simd_bit_xor)) {
      expect_simd_identical(f, x, y)
      expect_simd_identical(f, x, i64_at(y, 1))
    }
    for (f in list(simd_bit_not, simd_popcount, simd_lzcnt, simd_tzcnt)) {
      expect_simd_identical(f, x)
      expect_simd_identical(f, x, na_check = FALSE)
    }
    for (f in list(simd_shl, simd_shr, simd_sar, simd_rotl, simd_rotr)) {
      expect_simd_identical(f, x, 13)
    }
    expect_simd_identical(simd_popcount_total, x, na.rm = TRUE)
  }
})

# ---- Conversions ---------------------------------------------------------------

test_that("simd_as_double of integer64 rounds to nearest and warns beyond 2^53", {
  prec <- "integer precision lost while converting to double"
  x <- i64_dec(c(
    "9007199254740993", "9007199254740995", "-9007199254740993", "9223372036854775807",
    "-9223372036854775807", "1152921504606846977", "NA", "0", "-1"
  ))
  expect_tier_warnings(prec, simd_as_double, x)
  expect_tiers_give(i64_value(x), function(x) suppressWarnings(simd_as_double(x)), x)
  expect_identical(suppressWarnings(simd_as_double(i64_dec("9007199254740993"))), 2^53)
  expect_identical(suppressWarnings(simd_as_double(i64_dec("9007199254740995"))), 2^53 + 4)
  expect_tier_warnings(character(0), simd_as_double, i64(c(2^53 - 1, -(2^53 - 1), NA)))
  expect_tier_warnings(prec, simd_as_double, i64(2^53))
  expect_tier_warnings(prec, simd_as_double, i64(-2^53))
  for (n in lens) {
    x <- rand_i64(n, seed = n + 70L)
    expect_simd_identical(function(x) suppressWarnings(simd_as_double(x)), x)
    expect_identical(suppressWarnings(simd_as_double(x)), i64_value(x))
  }
})

test_that("simd_as_integer64 converts doubles, integers, logicals and raw", {
  expect_tier_warnings(ovf_msg, simd_as_integer64, c(2^63, -2^63, 9.5, NaN))
  expect_tiers_i64(i64(c(NA, NA, 9, NA)), simd_as_integer64, c(2^63, -2^63, 9.5, NaN),
    msgs = ovf_msg
  )
  expect_tiers_i64(i64(c(-9, NA, 2^62)), simd_as_integer64, c(-9.5, NA, 2^62))
  expect_tiers_i64(i64_dec("-9223372036854774784"), simd_as_integer64, -2^63 + 1024)
  expect_tiers_i64(i64(NA), simd_as_integer64, Inf, msgs = ovf_msg)
  sat <- c(Inf, -Inf, NaN, 1e300, -1e300, 5.7, -2^63)
  expect_tiers_i64(
    i64c(
      max64, i64_dec("-9223372036854775807"), na64, max64,
      i64_dec("-9223372036854775807"), i64(5), i64_dec("-9223372036854775807")
    ),
    simd_as_integer64, sat,
    mode = "saturating"
  )
  tr <- c(2^64 + 2^12, 2^63, 2^63 + 2048, -2^63 - 4096, Inf, NaN, -7.9)
  expect_tiers_i64(i64_dec(c(
    "4096", NA, "-9223372036854773760", "9223372036854771712", NA, NA,
    "-7"
  )), simd_as_integer64, tr, mode = "truncating")
  expect_tiers_i64(i64(c(1, NA, -5)), simd_as_integer64, c(1L, NA, -5L))
  expect_tiers_i64(i64(c(1, 0, NA)), simd_as_integer64, c(TRUE, FALSE, NA))
  expect_tiers_i64(i64(c(0, 255)), simd_as_integer64, as.raw(c(0, 255)))
  x <- i64c(max64, na64)
  expect_tiers_i64(x, simd_as_integer64, x)
  expect_tiers_i64(x, simd_as_integer64, structure(unclass(x),
    names = c("a", "b"),
    class = "integer64"
  ))
  expect_error(simd_as_integer64(1i), "does not support 'x' of type complex", fixed = TRUE)
  expect_error(simd_as_integer64("1"), "'x' must be an atomic vector")
  for (n in lens) {
    d <- with_seed(n, c(stats::runif(n, -2^64, 2^64), NA, NaN, Inf, -Inf, 0.5))
    for (mode in c("checked", "saturating", "truncating")) {
      expect_simd_identical(function(x) suppressWarnings(simd_as_integer64(x, mode)), d)
    }
    ii <- with_seed(n, sample(c(-5:5, NA, .Machine$integer.max), n, replace = TRUE))
    expect_simd_identical(simd_as_integer64, ii)
  }
})

test_that("integer64 converts to integer, logical and raw", {
  iovf <- "NAs produced by integer overflow"
  x <- i64(c(2^31 - 1, -(2^31 - 1), 2^31, -2^31, NA, 2^40 + 5))
  expect_tier_warnings(iovf, simd_as_integer, x)
  expect_tiers_give(
    c(.Machine$integer.max, -.Machine$integer.max, NA, NA, NA, NA),
    function(x) suppressWarnings(simd_as_integer(x)), x
  )
  expect_tiers_give(
    c(
      .Machine$integer.max, -.Machine$integer.max, .Machine$integer.max,
      -.Machine$integer.max, NA, .Machine$integer.max
    ),
    simd_as_integer, x,
    mode = "saturating"
  )
  expect_tiers_give(c(.Machine$integer.max, -.Machine$integer.max, NA, NA, NA, 5L),
    simd_as_integer, x,
    mode = "truncating"
  )
  expect_tiers_give(-1L, simd_as_integer, i64(2^32 - 1), mode = "truncating")
  expect_tiers_give(c(TRUE, FALSE, NA, TRUE), simd_as_logical, i64(c(-3, 0, NA, 2^60)))
  rw <- i64(c(0, 255, 256, -1, NA))
  expect_tier_warnings("out-of-range values treated as 0 in coercion to raw", simd_as_raw, rw)
  expect_tiers_give(as.raw(c(0, 255, 0, 0, 0)), function(x) suppressWarnings(simd_as_raw(x)), rw)
  expect_tiers_give(as.raw(c(0, 255, 255, 0, 0)), simd_as_raw, rw, mode = "saturating")
  expect_tiers_give(as.raw(c(0, 255, 0, 255, 0)), simd_as_raw, rw, mode = "truncating")
  for (n in lens) {
    x <- rand_i64(n, seed = n + 80L)
    for (mode in c("checked", "saturating", "truncating")) {
      expect_simd_identical(function(x) suppressWarnings(simd_as_integer(x, mode)), x)
      expect_simd_identical(function(x) suppressWarnings(simd_as_raw(x, mode)), x)
    }
    expect_simd_identical(simd_as_logical, x)
  }
})

# ---- Attributes ------------------------------------------------------------------

test_that("integer64 results keep only the class", {
  x <- structure(unclass(i64(c(1, 2))), names = c("a", "b"), dim = 2L, class = "integer64")
  expect_i64(simd_add(x, 1L), i64(c(2, 3)))
  expect_i64(simd_sum(x), i64(3))
  expect_i64(simd_cumsum(x), i64(c(1, 3)))
  expect_i64(simd_bit_not(x), i64(c(-2, -3)))
  expect_identical(simd_eq(x, 1L), c(TRUE, FALSE))
})
