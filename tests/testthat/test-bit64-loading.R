# rsimd loads bit64 before it hands out an integer64 object, so that
# bit64's methods are registered; without bit64 installed, only a call that
# would make an integer64 result errors. Run in child processes, where bit64
# starts out unloaded.

test_that("creating integer64 loads bit64 when it is installed but not loaded", {
  skip_if_not_installed("bit64")
  skip_if_not_installed("callr")
  child <- function() {
    library(rsimd)
    # An integer64 vector made without bit64, as readRDS() would give it.
    raw64 <- structure(unclass(simd_as_integer64(c(3, -1))), class = "integer64")
    paths <- list(
      as_integer64 = function() simd_as_integer64(1),
      add = function() simd_add(raw64, 1L),
      sum = function() simd_sum(raw64),
      promote = function() simd_add(1L, raw64),
      simd_vec = function() simd_vec(raw64),
      unwrap = function() simd_unwrap(simd_vec(raw64))
    )
    loaded <- vapply(paths, function(f) {
      # Unloading also checks that rsimd forgets it had loaded bit64.
      unloadNamespace("bit64")
      stopifnot(!isNamespaceLoaded("bit64"))
      f()
      isNamespaceLoaded("bit64")
    }, NA)
    unloadNamespace("bit64")
    v <- simd_vec(simd_as_integer64(c(3, -1, 2, 2^40)))
    df <- data.frame(v = v)
    list(
      loaded = loaded,
      sort = as.character(sort(v)),
      unique = as.character(unique(v)),
      df_class = class(df$v),
      df = format(df$v)
    )
  }
  res <- callr::r(child)
  expect_true(all(res$loaded), label = paste(names(res$loaded)[!res$loaded], collapse = ", "))
  expect_identical(res$sort, c("-1", "2", "3", "1099511627776"))
  expect_identical(res$unique, c("3", "-1", "2", "1099511627776"))
  expect_identical(res$df_class, "integer64")
  expect_identical(trimws(res$df), c("3", "-1", "2", "1099511627776"))
})

test_that("without bit64 only integer64 results error, with an install hint", {
  skip_if_not_installed("callr")
  skip_on_cran()
  # A library holding rsimd but not bit64.
  lib <- tempfile("lib")
  dir.create(lib)
  link <- file.path(lib, "rsimd")
  if (!file.symlink(find.package("rsimd"), link)) skip("cannot symlink rsimd")
  child <- function() {
    library(rsimd)
    # integer64 3, from its bits: simd_as_integer64() itself needs bit64.
    raw64 <- structure(readBin(as.raw(c(3, 0, 0, 0, 0, 0, 0, 0)), "double"), class = "integer64")
    try_msg <- function(expr) {
      tryCatch(
        {
          expr
          "ok"
        },
        error = conditionMessage
      )
    }
    list(
      installed = requireNamespace("bit64", quietly = TRUE),
      as_integer64 = try_msg(simd_as_integer64(1)),
      add = try_msg(simd_add(raw64, 1L)),
      sum = try_msg(simd_sum(raw64)),
      simd_vec = try_msg(simd_vec(raw64)),
      as_double = simd_as_double(raw64),
      eq = simd_eq(raw64, raw64),
      which_max = simd_which_max(raw64)
    )
  }
  res <- callr::r(child, libpath = c(lib, .Library))
  unlink(link)
  unlink(lib, recursive = TRUE)
  if (res$installed) skip("bit64 is in R's system library")
  hint <- "integer64 results need package 'bit64'; install it with install.packages(\"bit64\")"
  expect_identical(res$as_integer64, hint)
  expect_identical(res$add, hint)
  expect_identical(res$sum, hint)
  expect_identical(res$simd_vec, hint)
  expect_identical(res$as_double, 3)
  expect_identical(res$eq, TRUE)
  expect_identical(res$which_max, 1L)
})
