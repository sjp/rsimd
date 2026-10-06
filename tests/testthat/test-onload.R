test_that(".onLoad sets default options", {
  expect_type(getOption("rsimd.impl"), "character")
  expect_identical(getOption("rsimd.precision"), "fast")
  expect_identical(getOption("rsimd.na_check"), TRUE)
})

test_that(".onLoad takes rsimd.impl from the option, then RSIMD_IMPL, then auto", {
  onload <- get(".onLoad", envir = asNamespace("rsimd"))
  old_impl <- simd_current()
  old_opts <- options(rsimd.impl = NULL, rsimd.precision = NULL, rsimd.na_check = NULL)
  old_env <- Sys.getenv("RSIMD_IMPL", unset = NA)
  on.exit({
    options(old_opts)
    if (is.na(old_env)) Sys.unsetenv("RSIMD_IMPL") else Sys.setenv(RSIMD_IMPL = old_env)
    simd_use(attr(old_impl, "requested"))
  })

  Sys.unsetenv("RSIMD_IMPL")
  onload(NULL, "rsimd")
  expect_identical(getOption("rsimd.impl"), "auto")
  expect_identical(attr(simd_current(), "requested"), "auto")

  options(rsimd.impl = NULL)
  Sys.setenv(RSIMD_IMPL = "none")
  onload(NULL, "rsimd")
  expect_identical(getOption("rsimd.impl"), "none")
  expect_identical(as.vector(simd_current()), "none")

  # An option set before loading wins over the environment variable.
  options(rsimd.impl = "auto")
  onload(NULL, "rsimd")
  expect_identical(getOption("rsimd.impl"), "auto")
  expect_identical(as.vector(simd_current()), simd_available()[1])
})

test_that(".onLoad warns about an unusable rsimd.impl and falls back to auto", {
  onload <- get(".onLoad", envir = asNamespace("rsimd"))
  old_impl <- simd_current()
  old_opts <- options(rsimd.impl = NULL)
  old_env <- Sys.getenv("RSIMD_IMPL", unset = NA)
  on.exit({
    options(old_opts)
    if (is.na(old_env)) Sys.unsetenv("RSIMD_IMPL") else Sys.setenv(RSIMD_IMPL = old_env)
    simd_use(attr(old_impl, "requested"))
  })

  Sys.setenv(RSIMD_IMPL = "bogus")
  expect_warning(onload(NULL, "rsimd"), "RSIMD_IMPL: unknown implementation 'bogus'")
  expect_identical(getOption("rsimd.impl"), "auto")

  Sys.unsetenv("RSIMD_IMPL")
  options(rsimd.impl = "rvv")
  expect_warning(onload(NULL, "rsimd"), "option rsimd.impl: implementation 'rvv' is not available")
  expect_identical(getOption("rsimd.impl"), "auto")
  expect_identical(attr(simd_current(), "requested"), "auto")
})

test_that(".onLoad keeps user-set precision and NA-check options", {
  onload <- get(".onLoad", envir = asNamespace("rsimd"))
  old_opts <- options(rsimd.precision = "compensated", rsimd.na_check = FALSE)
  on.exit(options(old_opts))

  onload(NULL, "rsimd")
  expect_identical(getOption("rsimd.precision"), "compensated")
  expect_identical(getOption("rsimd.na_check"), FALSE)
})

test_that("unloading the namespace unloads the shared library", {
  skip_if_not_installed("callr")
  child <- function() {
    loadNamespace("rsimd")
    loaded <- "rsimd" %in% names(getLoadedDLLs())
    unloadNamespace("rsimd")
    unloaded <- !("rsimd" %in% names(getLoadedDLLs()))
    # Loading it again works.
    reloaded <- rsimd::simd_sum(c(1, 2, 3))
    c(loaded = loaded, unloaded = unloaded, reloaded = reloaded == 6)
  }
  res <- callr::r(child)
  expect_true(all(res), label = paste(names(res)[!res], collapse = ", "))
})
