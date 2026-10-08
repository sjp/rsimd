# Warnings and errors from the C and R sides, and the R side's integer64
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

test_that("errors from the R side name the user's call", {
  f <- function(y) simd_pmin(y, 1, na.rm = NA)
  expect_identical(condition_call(f(1)), quote(simd_pmin(y, 1, na.rm = NA)))
  f <- function(y) simd_shl(y, "a")
  expect_identical(condition_call(f(1L)), quote(simd_shl(y, "a")))
  f <- function(y) simd_var(y)
  expect_identical(condition_call(f(1i)), quote(simd_var(y)))
  expect_identical(condition_call(simd_add(Sys.Date(), 1)), quote(simd_add(Sys.Date(), 1)))
  expect_identical(condition_call(simd_sum(1, 2)), quote(simd_sum(1, 2)))
  expect_identical(condition_call(simd_round(1, "a")), quote(simd_round(1, "a")))
  expect_identical(condition_call(simd_rotl(1L, Inf)), quote(simd_rotl(1L, Inf)))
  expect_identical(condition_call(simd_ldexp(1, 0.5)), quote(simd_ldexp(1, 0.5)))
  expect_identical(condition_call(simd_arg("a")), quote(simd_arg("a")))
  expect_identical(condition_call(simd_lt(1i, 1)), quote(simd_lt(1i, 1)))
  expect_identical(
    condition_call(simd_hamming_bits(1L, as.raw(1))),
    quote(simd_hamming_bits(1L, as.raw(1)))
  )
  expect_identical(condition_call(simd_use("bogus")), quote(simd_use("bogus")))
  expect_identical(
    condition_call(simd_vec(1, check_na = NA)),
    quote(simd_vec(1, check_na = NA))
  )
  # Through simd_with_impl() and lapply(), the call into rsimd is named.
  expect_identical(
    condition_call(simd_with_impl("none", simd_shl(1L, "a"))),
    quote(simd_shl(1L, "a"))
  )
  expect_identical(
    condition_call(lapply(1, simd_round, digits = "a")),
    quote(FUN(X[[i]], ...))
  )
  # An invalid rsimd.impl option names the call that found it.
  old_opts <- options(rsimd.impl = "bogus")
  on.exit(options(old_opts))
  expect_identical(condition_call(simd_sum(1)), quote(simd_sum(1)))
})

test_that("errors from simd_vec methods name the method's call, not an inner helper", {
  v <- simd_vec(c(1, 2))
  expect_identical(condition_call(v + "a"), quote(Ops.simd_vec(v, "a")))
  expect_identical(condition_call(+simd_vec(as.raw(1))), quote(Ops.simd_vec(simd_vec(as.raw(1)))))
})

test_that("R/ raises no call-less errors", {
  skip_on_cran()
  dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(dir), "needs a source checkout")
  # Only the load-time warnings in aaa-onload.R, which have no user call.
  files <- setdiff(list.files(dir, "[.]R$", full.names = TRUE), file.path(dir, "aaa-onload.R"))
  hits <- unlist(lapply(files, function(f) {
    lines <- readLines(f)
    paste0(basename(f), ":", grep("call. = FALSE", lines, fixed = TRUE))
  }))
  hits <- hits[!grepl(":$", hits)]
  expect_identical(hits, character(), label = paste(hits, collapse = ", "))
})
