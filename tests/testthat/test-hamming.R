# Hamming distances: simd_hamming() equals sum(x != y) (with na.rm) for every
# type and promotion; simd_hamming_bits() is the population count of x ^ y.

# Base R's count, as a double.
ham_ref <- function(x, y, na.rm = FALSE) as.double(sum(x != y, na.rm = na.rm))

# The number of differing bits of two integer vectors, NA by its pattern.
bits_ref <- function(x, y) {
  as.double(sum(as.integer(xor(intToBits(x), intToBits(y)))))
}

expect_hamming <- function(x, y) {
  for (rm in c(FALSE, TRUE)) {
    expect_tiers_give(ham_ref(x, y, rm), simd_hamming, x, y, na.rm = rm)
    expect_tiers_give(ham_ref(y, x, rm), simd_hamming, y, x, na.rm = rm)
  }
}

test_that("simd_hamming matches sum(x != y) for every type, with NA and NaN", {
  x <- c(1, 2, NA, NaN, 0, -0, Inf, 3, NaN, 5)
  y <- c(1, 3, 1, 1, -0, 0, Inf, NA, NaN, 6)
  expect_identical(simd_hamming(x, y), NA_real_)
  expect_identical(simd_hamming(x, y, na.rm = TRUE), 2)
  expect_identical(simd_hamming(c(NaN, 1), c(NaN, 2)), NA_real_)
  expect_identical(simd_hamming(numeric(0), numeric(0)), 0)
  expect_identical(simd_hamming(1, numeric(0)), 0)
  expect_hamming(x, y)
  expect_hamming(c(1L, NA, 3L, 4L), c(1L, 2L, NA, 5L))
  expect_hamming(c(TRUE, NA, FALSE), c(TRUE, TRUE, TRUE))
  expect_hamming(c(1L, NA, 3L, 4L), c(1, 2, NaN, 4.5))
  expect_hamming(as.raw(c(0, 1, 255)), as.raw(c(0, 2, 255)))
  expect_hamming(c(1 + 1i, NA, complex(real = 1, imaginary = NaN), 2i), c(1 + 2i, 1, 1, 2i))
  expect_hamming(c(1 + 0i, 2 + 0i, 3 + 1i), c(1, NA, 3))
})

test_that("operands are promoted as base R's != promotes them", {
  r <- as.raw(c(0, 1, 2))
  expect_identical(simd_hamming(r, c(0L, 1L, 3L)), ham_ref(r, c(0L, 1L, 3L)))
  expect_identical(simd_hamming(r, c(FALSE, TRUE, TRUE)), ham_ref(r, c(FALSE, TRUE, TRUE)))
  expect_identical(simd_hamming(c(0.5, 1), r[1:2]), ham_ref(c(0.5, 1), r[1:2]))
  expect_identical(simd_hamming(r, 1 + 0i), ham_ref(r, 1 + 0i))
  expect_identical(simd_hamming(c(TRUE, FALSE), c(1i, 0i)), ham_ref(c(TRUE, FALSE), c(1i, 0i)))
})

test_that("length-1 operands broadcast and other length mismatches are errors", {
  x <- c(1, 2, 2, NA)
  expect_hamming(x, 2)
  expect_hamming(c(1L, 2L, 2L), 2L)
  expect_hamming(c(1 + 1i, 2 + 0i), 2 + 0i)
  expect_hamming(as.raw(1:3), as.raw(2))
  expect_error(simd_hamming(1:3, 1:2), "lengths of 'x' (3) and 'y' (2) must be equal", fixed = TRUE)
  expect_error(simd_hamming(c(1i, 2i), 1:3 + 0i), "must be equal or one of them must be 1")
  expect_error(simd_hamming_bits(1:3, 1:2), "must be equal or one of them must be 1")
})

test_that("simd_hamming matches base R on random vectors of edge lengths", {
  for (n in sweep_lengths()) {
    a <- with_seed(n, sample(c(1:3, NA, NaN, -0, 0), n, replace = TRUE))
    b <- with_seed(n + 1L, sample(c(1:3, NA, 0), n, replace = TRUE))
    clean_a <- with_seed(n + 2L, sample(c(1:3, -0), n, replace = TRUE))
    for (y in sweep_inputs(n, list(b, clean_a[rev(seq_len(n))]))) {
      expect_hamming(a, y)
      expect_hamming(clean_a, y)
      expect_hamming(as.integer(a), as.integer(y))
      expect_hamming(as.integer(clean_a), y)
      expect_hamming(complex(real = a, imaginary = rev(a)), complex(real = y, imaginary = 1))
    }
    ra <- as.raw(with_seed(n + 3L, sample(0:3, n, replace = TRUE)))
    expect_hamming(ra, rev(ra))
  }
})

test_that("an NA early in a long vector gives NA; na.rm counts the rest", {
  n <- 2^20 + 7
  x <- rep(c(1, 2), length.out = n)
  y <- rep(1, n)
  x[3] <- NA
  expect_tiers_give(NA_real_, simd_hamming, x, y)
  expect_tiers_give(ham_ref(x, y, TRUE), simd_hamming, x, y, na.rm = TRUE)
  xi <- as.integer(x)
  expect_tiers_give(NA_real_, simd_hamming, xi, 1L)
  expect_tiers_give(ham_ref(xi, 1L, TRUE), simd_hamming, xi, 1L, na.rm = TRUE)
})

test_that("compact sequences are compared without expanding them", {
  for (x in altrep_inputs(5000)) {
    expect_true(takes_region_path(x))
    expect_tiers_give(ham_ref(x, rev(x)), simd_hamming, x, rev(x))
    expect_tiers_give(4999, simd_hamming, x, 1)
  }
  expect_tiers_give(bits_ref(1:5000, 5000:1), simd_hamming_bits, 1:5000, 5000:1)
})

test_that("simd_hamming_bits counts differing bits of integer, logical and raw", {
  x <- c(0L, 1L, -1L, NA, .Machine$integer.max, 12345L)
  y <- c(0L, 2L, 0L, 0L, NA, -12345L)
  expect_tiers_give(bits_ref(x, y), simd_hamming_bits, x, y)
  expect_tiers_give(32, simd_hamming_bits, -1L, 0L)
  expect_tiers_give(1, simd_hamming_bits, NA_integer_, 0L)
  expect_tiers_give(bits_ref(c(TRUE, NA), c(FALSE, FALSE)), simd_hamming_bits, c(TRUE, NA), FALSE)
  expect_tiers_give(bits_ref(c(TRUE, NA), 3L), simd_hamming_bits, c(TRUE, NA), 3L)
  r <- as.raw(c(0, 255, 15, 1))
  expect_tiers_give(20, simd_hamming_bits, r, as.raw(c(255, 0, 0, 1)))
  expect_tiers_give(13, simd_hamming_bits, r, as.raw(0))
  expect_identical(simd_hamming_bits(integer(0), integer(0)), 0)
  for (n in sweep_lengths()) {
    a <- with_seed(n, sample(c(-50:50, NA, .Machine$integer.max), n, replace = TRUE))
    b <- with_seed(n + 1L, sample(c(-50:50, NA), n, replace = TRUE))
    expect_tiers_give(bits_ref(a, b), simd_hamming_bits, a, b)
    expect_tiers_give(bits_ref(a, 77L), simd_hamming_bits, a, 77L)
    ra <- as.raw(with_seed(n + 2L, sample(0:255, n, replace = TRUE)))
    expect_tiers_give(bits_ref(as.integer(ra), as.integer(rev(ra))), simd_hamming_bits, ra, rev(ra))
  }
})

test_that("simd_hamming_bits rejects other types and mixed widths", {
  expect_error(simd_hamming_bits(1, 2), "simd_hamming_bits() does not support 'x' of type double",
    fixed = TRUE
  )
  expect_error(simd_hamming_bits(1L, 1i), "does not support 'y' of type complex", fixed = TRUE)
  expect_error(simd_hamming_bits(1L, as.raw(1)),
    "simd_hamming_bits() needs 'x' and 'y' of the same width: integer and raw",
    fixed = TRUE
  )
  expect_error(simd_hamming_bits(1L, i64(1)), "same width: integer and integer64", fixed = TRUE)
  expect_error(simd_hamming("a", "b"), "must be an atomic vector")
  expect_error(simd_hamming(1, 2, na.rm = NA), "'na.rm' must be TRUE or FALSE", fixed = TRUE)
})
