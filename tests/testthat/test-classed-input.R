# Classed data other than integer64 and simd_vec is rejected by every
# export, in every data argument, with one message.

# Exports that take no data argument, or only a simd_vec.
no_data <- c(
  "simd_available", "simd_complex_variants", "simd_cpu_features", "simd_current",
  "simd_math_accuracy", "simd_na_check", "simd_options", "simd_precision", "simd_tiers",
  "simd_use", "simd_build_info",
  "simd_with_impl",
  "simd_impl", "simd_impl<-", "simd_na_free", "is_simd_vec"
)

# Scalar parameters, which keep their own validation.
scalar_args <- list(
  simd_shl = "n", simd_shr = "n", simd_sar = "n", simd_rotl = "n", simd_rotr = "n"
)

classed <- list(
  Date = as.Date("2026-01-01") + 0:2,
  POSIXct = as.POSIXct("2026-01-01", tz = "UTC") + 0:2,
  difftime = as.difftime(c(1, 2, 3), units = "hours"),
  foo = structure(1:3, class = "foo")
)

test_that("every export rejects classed data in every data argument", {
  exports <- sort(setdiff(getNamespaceExports("rsimd"), no_data))
  ns <- asNamespace("rsimd")
  checked <- 0L
  batch_expectations(for (f in exports) {
    fun <- get(f, envir = ns)
    fm <- formals(fun)
    required <- setdiff(names(fm)[vapply(fm, identical, NA, quote(expr = ))], "...")
    if (!length(required)) required <- names(fm)[[1L]]
    data_args <- setdiff(required, scalar_args[[f]])
    for (pos in data_args) {
      for (cls in names(classed)) {
        args <- lapply(required, function(a) {
          if (a == pos) classed[[cls]] else if (a %in% scalar_args[[f]]) 1L else 1:3
        })
        names(args) <- required
        msg <- paste0(
          "'", pos, "' must be an atomic vector \\(double, integer, logical, raw, complex or ",
          "integer64\\), not ", cls
        )
        check_error(do.call(fun, args), msg, info = paste0(f, "(", pos, " = <", cls, ">)"))
        checked <- checked + 1L
      }
    }
  })
  # Every export with data was reached.
  expect_gt(checked, 4L * 150L)
})

test_that("the simd_vec constructor and the C side give the same message", {
  msg <- "'x' must be an atomic vector \\(double, integer, logical, raw, complex or integer64\\),"
  for (cls in names(classed)) {
    expect_error(simd_vec(classed[[cls]]), paste(msg, "not", cls))
    expect_error(as_simd_vec(classed[[cls]]), paste(msg, "not", cls))
    expect_error(simd_sum(classed[[cls]]), paste(msg, "not", cls))
    expect_error(simd_unwrap(classed[[cls]]), paste(msg, "not", cls))
  }
  expect_error(simd_sum(factor("a")), paste(msg, "not factor"))
  expect_error(simd_sum(list(1)), paste(msg, "not list"))
  expect_error(simd_vec(list(1)), paste(msg, "not list"))
  expect_error(simd_sum("a"), paste(msg, "not character"))
  expect_error(simd_vec("a"), paste(msg, "not character"))
})

test_that("integer64, simd_vec and unclassed attributes are still taken", {
  expect_identical(simd_sum(matrix(1:4, 2)), 10L)
  expect_identical(simd_sum(c(a = 1, b = 2)), 3)
  expect_identical(simd_sum(simd_vec(c(1, 2))), 3)
  if (has_bit64()) expect_identical(simd_add(simd_as_integer64(1), 1L), simd_as_integer64(2))
})
