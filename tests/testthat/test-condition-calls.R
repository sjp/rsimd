# Warnings and errors from the C side, and the R side's integer64
# coercion warning, name the user's call, not an internal helper.

# The call of the first warning or error that `expr` signals.
condition_call <- function(expr) {
  tryCatch(expr, warning = conditionCall, error = conditionCall)
}

test_that("warnings name the user's call", {
  expect_identical(
    condition_call(simd_add(.Machine$integer.max, 1L)),
    quote(simd_add(.Machine$integer.max, 1L))
  )
  expect_identical(condition_call(simd_sqrt(-1)), quote(simd_sqrt(-1)))
  expect_identical(condition_call(simd_log(-1)), quote(simd_log(-1)))
  expect_identical(condition_call(simd_sincos(-Inf)), quote(simd_sincos(-Inf)))
  expect_identical(condition_call(simd_as_integer(1e10)), quote(simd_as_integer(1e10)))
  expect_identical(condition_call(simd_as_raw(300)), quote(simd_as_raw(300)))
  expect_identical(condition_call(simd_min(numeric())), quote(simd_min(numeric())))
  expect_identical(
    condition_call(simd_cumsum(c(.Machine$integer.max, 1L))),
    quote(simd_cumsum(c(.Machine$integer.max, 1L)))
  )
  x <- simd_vec(-1)
  expect_identical(condition_call(simd_log(x)), quote(simd_log(x)))
})

test_that("errors from the C side name the user's call", {
  expect_identical(condition_call(simd_fma(1:3, 1:2, 1)), quote(simd_fma(1:3, 1:2, 1)))
  expect_identical(condition_call(simd_hypot(1:3, 1:2)), quote(simd_hypot(1:3, 1:2)))
  expect_identical(condition_call(simd_scaleb(1:3, 1:2)), quote(simd_scaleb(1:3, 1:2)))
  expect_identical(condition_call(simd_dot(1:3, 1:2)), quote(simd_dot(1:3, 1:2)))
  expect_identical(condition_call(simd_prod_sums(1:3, 1:2)), quote(simd_prod_sums(1:3, 1:2)))
  expect_identical(condition_call(simd_lt(1:3, 1:2)), quote(simd_lt(1:3, 1:2)))
  expect_identical(
    condition_call(simd_and(1:3, c(TRUE, FALSE))),
    quote(simd_and(1:3, c(TRUE, FALSE)))
  )
  expect_identical(condition_call(simd_bit_and(1:3, 1:2)), quote(simd_bit_and(1:3, 1:2)))
  expect_identical(condition_call(simd_hamming(1:3, 1:2)), quote(simd_hamming(1:3, 1:2)))
  expect_identical(condition_call(simd_is_na(list(1))), quote(simd_is_na(list(1))))
  expect_identical(condition_call(simd_which_max(list(1))), quote(simd_which_max(list(1))))
  expect_identical(condition_call(simd_conj(list(1))), quote(simd_conj(list(1))))
  expect_identical(condition_call(simd_softmax(list(1))), quote(simd_softmax(list(1))))
  expect_identical(condition_call(simd_sum(list(1))), quote(simd_sum(list(1))))
})

test_that("integer64 warnings name the user's call", {
  skip_if_not_installed("bit64")
  z <- bit64::as.integer64(2)^62
  expect_identical(condition_call(simd_mul(z, z)), quote(simd_mul(z, z)))
  expect_identical(condition_call(simd_sum(c(z, z))), quote(simd_sum(c(z, z))))
  expect_identical(condition_call(simd_add(z, 1.5)), quote(simd_add(z, 1.5)))
  expect_identical(condition_call(simd_eq(z, 1.5)), quote(simd_eq(z, 1.5)))
  expect_identical(condition_call(simd_hamming(z, 1.5)), quote(simd_hamming(z, 1.5)))
})
