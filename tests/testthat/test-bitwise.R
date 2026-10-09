# Bitwise ops, shifts, rotates and bit counts: integers against base R's
# bitw* family and exact references, raw bytes exhaustively, on every tier.

imax <- .Machine$integer.max

# Integers whose bit patterns cover the boundaries, and random ones.
boundary_bits <- function() {
  c(
    0L, 1L, -1L, imax, -imax, 0x7FFF0000L, 0x0000FFFFL, 2L, -2L, 0x40000000L, 0x55555555L,
    -0x55555556L, with_seed(7L, sample.int(imax, 12) * sample(c(-1L, 1L), 12, TRUE)), NA
  )
}

# References on unsigned values in doubles.
sar_ref <- function(x, k) as.integer(floor(as.double(x) / 2^k))
rotl_ref <- function(x, k) {
  u <- u32(x)
  if (k == 0) {
    return(x)
  }
  from_u32((u * 2^k) %% 2^32 + u %/% 2^(32 - k))
}
rotr_ref <- function(x, k) if (k == 0) x else rotl_ref(x, 32 - k)

test_that("and, or, xor and not match base R's bitw* functions", {
  p <- all_pairs(boundary_bits())
  expect_tiers_give(bitwAnd(p$x, p$y), simd_bit_and, p$x, p$y)
  expect_tiers_give(bitwOr(p$x, p$y), simd_bit_or, p$x, p$y)
  expect_tiers_give(bitwXor(p$x, p$y), simd_bit_xor, p$x, p$y)
  x <- boundary_bits()
  expect_tiers_give(bitwNot(x), simd_bit_not, x)
  expect_identical(simd_bit_and(5L, 3L), 1L)
  expect_identical(simd_bit_not(0L), -1L)
  expect_identical(simd_bit_and(NA_integer_, 1L), NA_integer_)
  # A result with the bit pattern 0x80000000 is NA, as in base R.
  expect_identical(simd_bit_not(imax), bitwNot(imax))
  expect_identical(simd_bit_xor(-imax, 1L), NA_integer_)
  expect_identical(simd_bit_xor(-imax, 1L), bitwXor(-imax, 1L))
})

test_that("shifts by every count match bitwShiftL, bitwShiftR and floor division", {
  x <- boundary_bits()
  for (k in 0:31) {
    expect_tiers_give(bitwShiftL(x, k), simd_shl, x, k)
    expect_tiers_give(bitwShiftR(x, k), simd_shr, x, k)
    expect_tiers_give(sar_ref(x, k), simd_sar, x, k)
  }
  for (k in c(-1, 32, 100, NA)) {
    expect_tiers_give(rep(NA_integer_, length(x)), simd_shl, x, k)
    expect_tiers_give(rep(NA_integer_, length(x)), simd_sar, x, k)
  }
  expect_identical(simd_shl(1L, 31L), NA_integer_)
  expect_identical(simd_shl(3L, 30L), bitwShiftL(3L, 30L))
  expect_identical(simd_shr(-1L, 1L), 2147483647L)
  expect_identical(simd_sar(-8L, 1L), -4L)
  expect_identical(simd_shl(1L, 2.7), 4L)
  expect_identical(simd_shl(TRUE, 3L), 8L)
  # Counts on and off the fast path for plain numbers in range.
  expect_identical(simd_shl(1L, 30.9), bitwShiftL(1L, 30L))
  expect_identical(simd_shl(1L, 31.5), NA_integer_)
  expect_identical(simd_shl(1L, -0.5), 1L)
  expect_identical(simd_shl(1L, c(k = 2L)), 4L)
  expect_identical(simd_shl(as.raw(1), 8), as.raw(0))
  expect_identical(simd_shr(as.raw(0x80), 7.9), as.raw(1))
  expect_error(simd_shl(1L, factor("a")), "'n' must be a single number", fixed = TRUE)
})

test_that("rotates by every count match the reference, modulo 32", {
  x <- boundary_bits()
  for (k in 0:31) {
    expect_tiers_give(rotl_ref(x, k), simd_rotl, x, k)
    expect_tiers_give(rotr_ref(x, k), simd_rotr, x, k)
  }
  expect_identical(simd_rotl(x, 33), simd_rotl(x, 1))
  expect_identical(simd_rotl(x, -1), simd_rotr(x, 1))
  expect_identical(simd_rotl(1L, 31L), NA_integer_)
  expect_error(simd_rotl(1L, NA), "'n' must be a finite number", fixed = TRUE)
  expect_identical(simd_rotl(x, 31.9), simd_rotl(x, 31L))
  expect_identical(simd_rotl(x, -0.5), x)
  expect_identical(simd_rotl(as.raw(0x81), 7.5), as.raw(0xc0))
  expect_identical(simd_rotl(as.raw(0x81), 9), as.raw(0x03))
  expect_error(simd_rotl(1L, Inf), "'n' must be a finite number", fixed = TRUE)
})

test_that("bit counts match the references on random and boundary integers", {
  rand <- with_seed(11L, sample.int(2^31 - 1, 1e5, TRUE) * sample(c(-1L, 1L), 1e5, TRUE))
  x <- c(boundary_bits(), rand)
  expect_tiers_give(popcount_ref(x), simd_popcount, x)
  expect_tiers_give(lzcnt_ref(x), simd_lzcnt, x)
  expect_tiers_give(tzcnt_ref(x), simd_tzcnt, x)
  expect_identical(simd_popcount(c(-1L, 0L)), c(32L, 0L))
  expect_identical(simd_lzcnt(0L), 32L)
  expect_identical(simd_tzcnt(0L), 32L)
  expect_identical(simd_popcount(c(TRUE, FALSE, NA)), c(1L, 0L, NA))
})

test_that("popcount_total sums the counts, NA unless na.rm", {
  x <- boundary_bits()
  ok <- x[!is.na(x)]
  expect_tiers_give(NA_real_, simd_popcount_total, x)
  total <- as.double(sum(popcount_ref(ok)))
  expect_tiers_give(total, simd_popcount_total, x, na.rm = TRUE)
  expect_tiers_give(total, simd_popcount_total, ok)
  expect_identical(simd_popcount_total(integer(0)), 0)
  expect_identical(simd_popcount_total(as.raw(0:255)), 1024)
  big <- rep(-1L, 2^20 + 7)
  expect_tiers_give(32 * (2^20 + 7), simd_popcount_total, big)
  big[length(big)] <- NA
  expect_tiers_give(NA_real_, simd_popcount_total, big)
})

test_that("na_check = FALSE reads NA as the bit pattern 0x80000000", {
  expect_identical(simd_popcount(NA_integer_, na_check = FALSE), 1L)
  expect_identical(simd_lzcnt(NA_integer_, na_check = FALSE), 0L)
  expect_identical(simd_shr(NA_integer_, 31L, na_check = FALSE), 1L)
  expect_identical(simd_popcount_total(c(NA, 1L), na_check = FALSE), 2)
  expect_identical(simd_popcount_total(c(NA, 1L), na.rm = TRUE, na_check = FALSE), 1)
})

test_that("every raw op and count matches base R or the reference, exhaustively", {
  v <- as.raw(0:255)
  iv <- 0:255
  for (k in 0:8) {
    expect_tiers_give(rawShift(v, k), simd_shl, v, k)
    expect_tiers_give(rawShift(v, -k), simd_shr, v, k)
  }
  for (k in 0:7) {
    rot <- as.raw((iv * 2^k) %% 256 + iv %/% 2^(8 - k))
    expect_tiers_give(rot, simd_rotl, v, k)
    expect_tiers_give(rot, simd_rotr, v, 8 - k)
  }
  p <- all_pairs(v)
  expect_tiers_give(p$x & p$y, simd_bit_and, p$x, p$y)
  expect_tiers_give(p$x | p$y, simd_bit_or, p$x, p$y)
  expect_tiers_give(xor(p$x, p$y), simd_bit_xor, p$x, p$y)
  expect_tiers_give(!v, simd_bit_not, v)
  bits <- matrix(as.integer(rawToBits(v)), 8)
  expect_tiers_give(as.integer(colSums(bits)), simd_popcount, v)
  expect_tiers_give(as.integer(ifelse(iv == 0, 8, 7 - floor(log2(pmax(iv, 1))))), simd_lzcnt, v)
  low_zeros <- apply(bits, 2, function(b) if (any(b == 1)) which(b == 1)[1] - 1 else 8)
  expect_tiers_give(as.integer(low_zeros), simd_tzcnt, v)
  expect_identical(simd_shl(as.raw(0x81), 1), as.raw(0x02))
  expect_identical(simd_shr(as.raw(0x81), 1), as.raw(0x40))
})

test_that("ops match on random vectors of edge lengths and compact sequences", {
  for (n in sweep_lengths()) {
    x <- rand_vec("integer", n, seed = n + 21L) * 1000L
    y <- rand_vec("integer", n, seed = n + 22L)
    expect_tiers_give(bitwAnd(x, y), simd_bit_and, x, y)
    expect_tiers_give(bitwXor(x, 5L), simd_bit_xor, x, 5L)
    expect_tiers_give(bitwShiftL(x, 3L), simd_shl, x, 3L)
    expect_tiers_give(rotl_ref(x, 7), simd_rotl, x, 7)
    expect_tiers_give(popcount_ref(x), simd_popcount, x)
    r <- as.raw(with_seed(n, sample.int(256, n, TRUE) - 1L))
    expect_tiers_give(rawShift(r, 3), simd_shl, r, 3)
  }
  x <- seq_len(5000)
  expect_true(takes_region_path(x))
  expect_tiers_give(bitwAnd(x, 255L), simd_bit_and, x, 255L)
  expect_tiers_give(popcount_ref(x), simd_popcount, x)
})

test_that("argument and type errors", {
  expect_error(simd_bit_and(1L, as.raw(1)), "simd_bit_and() cannot combine integer and raw operands",
    fixed = TRUE
  )
  expect_error(simd_bit_and(1.5, 1L), "simd_bit_and() does not support 'x' of type double",
    fixed = TRUE
  )
  expect_error(simd_shl(2, 1L), "simd_shl() does not support 'x' of type double", fixed = TRUE)
  expect_error(simd_popcount_total(2), "does not support 'x' of type double", fixed = TRUE)
  expect_error(simd_sar(as.raw(1), 1L), "simd_sar() does not support 'x' of type raw", fixed = TRUE)
  expect_error(simd_shl(as.raw(1), 9), "argument 'n' must be a small integer", fixed = TRUE)
  expect_error(simd_shl(as.raw(1), -1), "argument 'n' must be a small integer", fixed = TRUE)
  expect_error(simd_shr(as.raw(1), NA), "argument 'n' must be a small integer", fixed = TRUE)
  expect_error(simd_shl(1L, 1:2), "'n' must be a single number", fixed = TRUE)
  expect_error(simd_shl(1L, "1"), "'n' must be a single number", fixed = TRUE)
  expect_error(simd_popcount(1i), "simd_popcount() does not support 'x' of type complex",
    fixed = TRUE
  )
  expect_identical(simd_bit_and(c(a = 1L), 1L), 1L)
})
