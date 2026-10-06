# Tests of the R vector access layer in src/rvec.c, driven through the
# internal debug routines of R/debug.R.

type_msg <- function(arg, what) {
  paste0(
    "'", arg, "' must be an atomic vector (double, integer, logical, raw, complex or ",
    "integer64), not ", what
  )
}

len_msg <- function(nx, ny) {
  paste0(
    "lengths of 'x' (", nx, ") and 'y' (", ny, ") must be equal or one of them must be 1"
  )
}

integer64 <- function(x) structure(x, class = "integer64")

test_that("element types are classified", {
  type <- function(x) .debug_regions(x)$type
  expect_identical(type(c(1, 2)), "double")
  expect_identical(type(c(1L, 2L)), "integer")
  expect_identical(type(c(TRUE, NA)), "logical")
  expect_identical(type(as.raw(1:3)), "raw")
  expect_identical(type(c(1i, 2)), "complex")
  # integer64 is recognised by its class, without bit64.
  expect_identical(type(integer64(0)), "integer64")
  expect_identical(type(matrix(1:4, 2)), "integer")
  expect_identical(type(structure(1, class = "Date")), "double")
})

test_that("unsupported inputs error with the class or type", {
  expect_error(.debug_regions(factor("a")), type_msg("x", "factor"), fixed = TRUE)
  expect_error(.debug_regions(data.frame(a = 1)), type_msg("x", "data.frame"), fixed = TRUE)
  expect_error(.debug_regions(list(1)), type_msg("x", "list"), fixed = TRUE)
  expect_error(.debug_regions("a"), type_msg("x", "character"), fixed = TRUE)
  expect_error(.debug_regions(NULL), type_msg("x", "NULL"), fixed = TRUE)
  expect_error(.debug_regions(sum), type_msg("x", "builtin"), fixed = TRUE)
  expect_error(.debug_regions(asS4(1)), "'x' must be an atomic vector", fixed = TRUE)
  expect_error(.debug_bin(1, factor("a")), type_msg("y", "factor"), fixed = TRUE)
})

test_that("errors report the R call that made the .Call", {
  err <- tryCatch(.debug_bin(c(1, 2), c(1, 2, 3)), error = identity)
  expect_identical(conditionMessage(err), len_msg(2, 3))
  expect_identical(conditionCall(err), quote(.debug_bin(c(1, 2), c(1, 2, 3))))
})

test_that("compact sequences are read in regions without materialising", {
  # R makes no compact sequence of length 1; that case is contiguous below.
  for (n in c(2, 4095, 4096, 4097, 2^20 - 1, 2^20, 2^20 + 1)) {
    expected <- n * (n + 1) / 2
    for (x in list(1:n, seq_len(n), as.double(1:n))) {
      r <- .debug_regions(x)
      expect_identical(r$path, "regions", info = paste(typeof(x), n))
      expect_identical(r$regions, ceiling(n / 4096), info = paste(typeof(x), n))
      expect_identical(r$sum, expected, info = paste(typeof(x), n))
      expect_identical(.debug_copy(x), x[seq_len(n)], info = paste(typeof(x), n))
      # Still compact afterwards: reading did not expand it.
      expect_identical(.debug_regions(x)$path, "regions", info = paste(typeof(x), n))
    }
  }
})

test_that("ordinary vectors take the contiguous path", {
  for (n in c(0, 1, 4095, 4096, 4097, 2^20 - 1, 2^20, 2^20 + 1)) {
    expected <- n * (n + 1) / 2
    for (x in list(c(seq_len(n)), rev(rev(seq_len(n))), c(as.double(seq_len(n))))) {
      r <- .debug_regions(x)
      expect_identical(r$path, "contiguous", info = paste(typeof(x), n))
      expect_identical(r$regions, ceiling(n / 2^20), info = paste(typeof(x), n))
      expect_identical(r$sum, expected, info = paste(typeof(x), n))
    }
  }
  # rev() of a sequence is an ordinary vector, read like any other.
  r <- .debug_regions(rev(1:10))
  expect_identical(r$path, "contiguous")
  expect_identical(r$sum, 55)
})

test_that("empty inputs visit no regions", {
  for (x in list(integer(0), double(0), seq_len(0), logical(0), raw(0), complex(0))) {
    r <- .debug_regions(x)
    expect_identical(r$regions, 0, info = typeof(x))
    expect_identical(r$sum, 0, info = typeof(x))
    expect_identical(.debug_copy(x), x, info = typeof(x))
  }
})

test_that("every element type is read and copied exactly", {
  inputs <- list(
    c(1.5, NA, -0, Inf, NaN),
    c(1L, NA, -5L, .Machine$integer.max),
    c(TRUE, NA, FALSE),
    as.raw(c(0, 1, 255)),
    c(1 + 2i, NA, 3 - 4i)
  )
  for (x in inputs) expect_identical(.debug_copy(x), x, info = typeof(x))
  expect_identical(.debug_regions(c(1L, NA, 3L))[c("sum", "na")], list(sum = 4, na = 1))
  expect_identical(.debug_regions(as.raw(c(1, 255)))$sum, 256)
  # Complex values are read as interleaved real and imaginary parts.
  expect_identical(.debug_regions(c(1 + 2i, 3 + 4i))$sum, 10)
  # integer64 values 5, NA (INT64_MIN) and -2, from their little-endian bytes.
  bytes <- c(5, rep(0, 7), rep(0, 7), 0x80, 0xfe, rep(0xff, 7))
  bits <- integer64(readBin(as.raw(bytes), "double", n = 3, size = 8, endian = "little"))
  expect_identical(.debug_regions(bits)[c("sum", "na")], list(sum = 3, na = 1))
})

test_that("the known-NA-free hint comes from ALTREP, raw and scalars", {
  expect_true(.debug_regions(1:10)$no_na_hint)
  expect_true(.debug_regions(as.double(1:10))$no_na_hint)
  # sort() returns a wrapper ALTREP that records it has no NA. Recent R clears
  # that record once a data pointer is taken, so each check gets a fresh one.
  wrapped <- .debug_regions(sort(c(3, 1, 2)))
  expect_identical(wrapped$path, "contiguous")
  expect_true(wrapped$no_na_hint)
  expect_true(.debug_regions(sort(c(3L, 1L)))$no_na_hint)
  expect_true(.debug_regions(as.raw(1))$no_na_hint)
  expect_false(.debug_regions(c(1, 2))$no_na_hint)
  expect_false(.debug_regions(c(1L, 2L))$no_na_hint)
  expect_false(.debug_regions(integer64(c(0, 1)))$no_na_hint)
  # A length-1 operand is checked directly.
  expect_true(.debug_regions(2)$no_na_hint)
  expect_true(.debug_regions(integer64(0))$no_na_hint)
  expect_false(.debug_regions(NA_real_)$no_na_hint)
  expect_false(.debug_regions(NaN)$no_na_hint)
  expect_false(.debug_regions(NA_integer_)$no_na_hint)
  expect_false(.debug_regions(NA)$no_na_hint)
  expect_false(.debug_regions(complex(real = 1, imaginary = NaN))$no_na_hint)
})

test_that("na_check is cleared for inputs known to be NA-free", {
  plain <- c(1, 2)
  expect_identical(.debug_opts(plain), list(na_rm = FALSE, na_check = TRUE, precision = 0L))
  expect_true(.debug_opts(plain, na_check = TRUE)$na_check)
  expect_false(.debug_opts(plain, na_check = FALSE)$na_check)
  expect_false(.debug_opts(1:10, na_check = TRUE)$na_check)
  expect_false(.debug_opts(sort(c(2, 1)), na_check = TRUE)$na_check)
  expect_true(.debug_opts(plain, na_rm = TRUE)$na_rm)
  expect_identical(.debug_opts(plain, precision = 2L)$precision, 2L)
  expect_identical(.debug_opts(plain, precision = 1)$precision, 1L)
})

test_that("argument errors name the argument", {
  expect_error(.debug_opts(1, na_rm = NA), "'na.rm' must be TRUE or FALSE", fixed = TRUE)
  expect_error(.debug_opts(1, na_rm = c(TRUE, FALSE)), "'na.rm' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(.debug_opts(1, na_check = "yes"), "'na_check' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(.debug_opts(1, precision = 1.5), "'precision' must be a single integer",
    fixed = TRUE
  )
  expect_error(.debug_opts(1, precision = 3L), "invalid precision code 3", fixed = TRUE)
  expect_error(.debug_finish(1, "double", 1), "'op' must be a single string", fixed = TRUE)
  expect_error(.debug_finish("sum", "double", NA), "'n' must be a single number", fixed = TRUE)
})

test_that("binary operands follow the length-1 broadcast rule", {
  x <- c(1, 2, 3, 4, 5)
  r <- .debug_bin(x, x)
  expect_identical(r$value, x + x)
  expect_false(r$x_scalar || r$y_scalar)

  r <- .debug_bin(x, 10)
  expect_identical(r$value, x + 10)
  expect_true(r$y_scalar)
  expect_false(r$x_scalar)

  r <- .debug_bin(10, x)
  expect_identical(r$value, 10 + x)
  expect_true(r$x_scalar)
  expect_false(r$y_scalar)

  # A single element on both sides is not a broadcast.
  r <- .debug_bin(1, 2)
  expect_identical(r$value, 3)
  expect_false(r$x_scalar || r$y_scalar)

  for (args in list(list(numeric(0), 1), list(1, numeric(0)), list(numeric(0), numeric(0)))) {
    r <- do.call(.debug_bin, args)
    expect_identical(r$value, numeric(0))
    expect_identical(r$n, 0)
    expect_identical(r$regions, 0)
  }

  expect_error(.debug_bin(x, c(1, 2)), len_msg(5, 2), fixed = TRUE)
  expect_error(.debug_bin(c(1, 2), x), len_msg(2, 5), fixed = TRUE)
  expect_error(.debug_bin(numeric(0), c(1, 2, 3)), len_msg(0, 3), fixed = TRUE)
  expect_error(.debug_bin(c(1, 2, 3), numeric(0)), len_msg(3, 0), fixed = TRUE)
})

test_that("binary chunking mixes ALTREP, contiguous and scalar operands", {
  for (n in c(4095, 4097, 2^20 + 1)) {
    alt <- as.double(seq_len(n))
    flat <- c(as.double(seq_len(n))) * 2
    expected <- 3 * as.double(seq_len(n))
    r <- .debug_bin(alt, flat)
    expect_identical(r$value, expected, info = n)
    # The ALTREP operand caps chunks at its region size.
    expect_identical(r$regions, ceiling(n / 4096), info = n)
    expect_identical(.debug_bin(flat, alt)$value, expected, info = n)
    expect_identical(.debug_bin(alt, 0.5)$value, alt + 0.5, info = n)
    expect_identical(.debug_bin(seq_len(n), 2L)$value, seq_len(n) + 2, info = n)
    expect_identical(.debug_bin(2L, seq_len(n))$value, seq_len(n) + 2, info = n)
    expect_identical(.debug_bin(flat, flat)$regions, ceiling(n / 2^20), info = n)
  }
  expect_identical(.debug_bin(c(1L, NA, 3L), TRUE)$value, c(2, NA, 4))
})

test_that("results are bare except for the integer64 and simd_vec classes", {
  expect_identical(.debug_copy(c(a = 1, b = 2)), c(1, 2))
  expect_identical(.debug_copy(matrix(1:4, 2, dimnames = list(c("a", "b"), NULL))), 1:4)
  expect_identical(.debug_copy(structure(1, class = "Date")), 1)
  expect_identical(.debug_copy(structure(1, foo = "bar")), 1)

  i64 <- structure(c(a = 0), class = "integer64")
  expect_identical(.debug_copy(i64), integer64(0))
  # Any class that inherits from "integer64" (as the access layer reads
  # it) gives a plain "integer64" result, a simd_vec one included.
  expect_identical(.debug_copy(structure(0, class = c("integer64", "foo"))), integer64(0))
  sv64 <- structure(c(1, 2), class = c("simd_vec", "integer64"), rsimd_impl = "none")
  expect_identical(.debug_copy(sv64), sv64)

  # A forged flag (no token) is unknown; a flagged result gets a token.
  sv <- structure(c(1, 2), class = "simd_vec", rsimd_impl = "none", rsimd_na_free = TRUE)
  expect_identical(
    .debug_copy(sv),
    structure(c(1, 2), class = "simd_vec", rsimd_impl = "none")
  )
  expect_null(simd_na_free(.debug_copy(sv, no_na = TRUE)))
  flagged <- .debug_copy(simd_vec(c(1, 2), impl = "none", check_na = TRUE), no_na = TRUE)
  expect_true(all.equal(flagged, simd_vec(c(1, 2), impl = "none", check_na = TRUE)))
})

test_that("integer64 class is kept when bit64 is loaded", {
  skip_if_not_installed("bit64")
  x <- bit64::as.integer64(c(1, -5, NA))
  out <- .debug_copy(x)
  expect_s3_class(out, "integer64")
  expect_identical(as.character(out), c("1", "-5", NA))
})

test_that("operand pairs are promoted following base R", {
  same <- c(1, 2)
  expect_identical(.promote_pair(same, 3), list(x = same, y = 3))
  expect_identical(.promote_pair(1L, 2.5), list(x = 1, y = 2.5))
  expect_identical(.promote_pair(2.5, TRUE), list(x = 2.5, y = 1))
  expect_identical(.promote_pair(TRUE, 2L), list(x = 1L, y = 2L))
  # Logicals stay logical; kernels read their int storage.
  expect_identical(.promote_pair(TRUE, NA), list(x = TRUE, y = NA))
  expect_identical(.promote_pair(1i, 2L), list(x = 1i, y = 2 + 0i))
  expect_identical(.promote_pair(as.raw(1), as.raw(2)), list(x = as.raw(1), y = as.raw(2)))

  expect_error(.promote_pair(as.raw(1), 1),
    "cannot combine 'x' of type raw with 'y' of type double",
    fixed = TRUE
  )
  expect_error(.promote_pair(1, as.raw(1)),
    "cannot combine 'x' of type double with 'y' of type raw",
    fixed = TRUE
  )
  expect_error(.promote_pair(1i, integer64(0)),
    "cannot combine 'x' of type complex with 'y' of type integer64",
    fixed = TRUE
  )
  expect_error(.promote_pair(factor("a"), 1), type_msg("x", "factor"), fixed = TRUE)
  expect_error(.promote_pair(1, "a"), type_msg("y", "character"), fixed = TRUE)
})

test_that("integer64 promotes to integer64, or to double with a warning", {
  skip_if_not_installed("bit64")
  big <- bit64::as.integer64("9007199254740993")
  out <- .promote_pair(big, 2L)
  expect_identical(out$x, big)
  expect_identical(out$y, bit64::as.integer64(2L))
  expect_identical(.promote_pair(TRUE, big)$x, bit64::as.integer64(1L))

  caller <- function(x, y) .promote_pair(x, y)
  w <- tryCatch(caller(big, 0.5), warning = identity)
  expect_identical(conditionMessage(w), "integer64 coerced to double")
  expect_identical(conditionCall(w), quote(caller(big, 0.5)))
  out <- suppressWarnings(caller(big, 0.5))
  expect_identical(out, list(x = 9007199254740992, y = 0.5))
  expect_warning(
    .promote_pair(1, bit64::as.integer64(7L)), "integer64 coerced to double",
    fixed = TRUE
  )
})

test_that("reductions start from their identity", {
  # Finishing an empty min or max warns; only the identity matters here.
  init <- function(op) suppressWarnings(.debug_finish(op, "double", 0))$init
  sum0 <- init("sum")
  expect_identical(sum0$f64, 0)
  expect_identical(sum0$comp, 0)
  expect_identical(sum0$i64, 0)
  expect_identical(sum0$idx, -1)
  expect_identical(sum0$count, 0)
  expect_false(any(unlist(sum0[c("saw_na", "saw_nan", "overflow", "any_true", "any_false")])))
  expect_identical(init("prod")$f64, 1)
  for (op in c("min", "which_min")) {
    expect_identical(init(op)$f64, Inf)
    expect_identical(init(op)$i64, 2^63)
  }
  for (op in c("max", "which_max")) {
    expect_identical(init(op)$f64, -Inf)
    expect_identical(init(op)$i64, -2^63)
  }
})

test_that("finished reductions have base R's result types", {
  # Every element survived unless the test says otherwise.
  fin <- function(op, type, fields = list(), n = 10, precision = 0L) {
    if (is.null(fields$count)) fields$count <- n
    .debug_finish(op, type, n, fields, precision)$value
  }
  # sum
  expect_identical(fin("sum", "integer", list(i64 = 55)), 55L)
  expect_identical(fin("sum", "logical", list(i64 = 2)), 2L)
  expect_identical(fin("sum", "integer", list(i64 = .Machine$integer.max)), .Machine$integer.max)
  expect_identical(fin("sum", "integer", list(i64 = -.Machine$integer.max)), -.Machine$integer.max)
  expect_identical(fin("sum", "integer", list(i64 = 3 * .Machine$integer.max)), 6442450941)
  expect_identical(fin("sum", "integer", list(i64 = -2^31)), -2^31)
  expect_identical(fin("sum", "double", list(f64 = 1.5, comp = 0.25)), 1.5)
  expect_identical(fin("sum", "double", list(f64 = 1.5, comp = 0.25), precision = 2L), 1.75)
  expect_identical(fin("sum", "integer64", list(i64 = 5)), integer64(2.5e-323))
  # always double
  for (op in c("prod", "mean", "sum_sq", "sum_abs", "dot", "norm", "dist", "cosine", "var", "sd")) {
    expect_identical(fin(op, "integer", list(f64 = 2)), 2, info = op)
    expect_identical(fin(op, "double", list(f64 = 2)), 2, info = op)
  }
  # min and max
  expect_identical(fin("min", "integer", list(i64 = -3)), -3L)
  expect_identical(fin("max", "logical", list(i64 = 1)), 1L)
  expect_identical(fin("max", "double", list(f64 = 2.5)), 2.5)
  expect_identical(fin("min", "integer64", list(i64 = 5)), integer64(2.5e-323))
  # which_*: 1-based, double for long vectors
  long <- 2^31 + 10
  expect_identical(fin("which_min", "double", list(idx = 4)), 5L)
  expect_identical(fin("which_max", "integer", list(idx = 0)), 1L)
  expect_identical(fin("which_min", "double", list(idx = 4), n = long), 5)
  expect_identical(fin("which_max", "double", list(idx = 2^31 + 1), n = long), 2^31 + 2)
  expect_identical(fin("which_min", "double"), integer(0))
  expect_identical(fin("which_min", "double", n = long), numeric(0))
  # logical results
  expect_identical(fin("any", "logical", list(any_true = 1)), TRUE)
  expect_identical(fin("any", "logical"), FALSE)
  expect_identical(fin("all", "logical", list(any_false = 1)), FALSE)
  expect_identical(fin("all", "logical"), TRUE)
  expect_identical(fin("any_na", "double", list(saw_nan = 1)), TRUE)
  expect_identical(fin("any_na", "integer", list(saw_na = 1)), TRUE)
  expect_identical(fin("any_na", "double"), FALSE)
  # count_na
  expect_identical(fin("count_na", "double", list(i64 = 3)), 3)
  expect_identical(fin("count_na", "double", list(i64 = 3), n = long), 3)
  # unsupported types
  expect_error(fin("sum", "raw"), "invalid 'type' (raw) of argument", fixed = TRUE)
  expect_error(fin("sum", "complex"), "invalid 'type' (complex) of argument", fixed = TRUE)
  expect_error(fin("mean", "integer64"), "invalid 'type' (integer64) of argument", fixed = TRUE)
  expect_error(fin("bogus", "double"), "unknown reduction 'bogus'", fixed = TRUE)
  expect_error(fin("sum", "bogus"), "unknown type 'bogus'", fixed = TRUE)
  expect_error(fin("sum", "double", list(bogus = 1)), "unknown field 'bogus'", fixed = TRUE)
})

test_that("RSIMD_DEBUG_STRIDE sets the chunk stride at load", {
  skip_if_no_subprocess()
  skip_on_cran()
  child <- function() {
    warnings <- character()
    r <- withCallingHandlers(
      rsimd:::.debug_regions(c(as.double(1:10000))),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(
      regions = r, warnings = warnings,
      alt = rsimd:::.debug_regions(as.double(1:10000)),
      bin = rsimd:::.debug_bin(as.double(1:10000), c(as.double(1:10000)))
    )
  }
  env <- function(stride) c(callr::rcmd_safe_env(), RSIMD_DEBUG_STRIDE = stride)

  res <- callr::r(child, env = env("4096"))
  expect_identical(res$regions$stride, 4096)
  expect_identical(res$regions$regions, 3)
  expect_identical(res$regions$sum, 50005000)
  expect_identical(res$bin$value, 2 * as.double(1:10000))
  expect_length(res$warnings, 0)

  # A stride below the region size also caps ALTREP regions. Strides are
  # rounded up to a multiple of 128 (the pairwise summation leaf).
  res <- callr::r(child, env = env("1000"))
  expect_identical(res$regions$stride, 1024)
  expect_identical(res$alt$regions, 10)
  expect_identical(res$alt$sum, 50005000)
  expect_identical(res$bin$regions, 10)
  expect_identical(res$bin$value, 2 * as.double(1:10000))

  for (bad in c("abc", "0", "-5", "1.5")) {
    res <- callr::r(child, env = env(bad))
    expect_identical(res$regions$stride, 2^20, info = bad)
    expect_match(res$warnings, "RSIMD_DEBUG_STRIDE must be a positive integer", fixed = TRUE)
  }
})

test_that("reading a compact sequence does not allocate it", {
  skip_on_cran()
  mb_used <- function() sum(gc(full = TRUE)[, 2L])
  for (x in list(seq_len(1e7), as.double(seq_len(1e7)))) {
    before <- mb_used()
    r <- .debug_regions(x)
    after <- mb_used()
    expect_identical(r$path, "regions")
    expect_lt(after - before, 1)
    expect_identical(.debug_regions(x)$path, "regions")
  }
})

test_that("long vectors are read with R_xlen_t indices", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("RSIMD_EXTENDED_TESTS"), "true"),
    "RSIMD_EXTENDED_TESTS is not true"
  )
  skip_if(.Machine$sizeof.pointer < 8, "no long vectors on 32-bit platforms")
  n <- 2^31 + 10
  x <- seq_len(n) # a compact double sequence: long, but takes no memory
  r <- .debug_regions(x)
  expect_identical(r$path, "regions")
  expect_identical(r$regions, ceiling(n / 4096))
  # A plain running sum of 2^31 doubles drifts by about 1e-9 relative.
  expect_equal(r$sum, n * (n + 1) / 2, tolerance = 1e-6)
})
