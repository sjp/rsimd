test_that(".onLoad sets default options", {
  expect_type(getOption("rsimd.impl"), "character")
  expect_identical(getOption("rsimd.precision"), "fast")
  expect_identical(getOption("rsimd.na_check"), TRUE)
})

test_that(".onLoad takes rsimd.impl from the option, then RSIMD_IMPL, then auto", {
  onload <- get(".onLoad", envir = asNamespace("rsimd"))
  old_impl <- getOption("rsimd.impl")
  old_opts <- options(rsimd.impl = NULL, rsimd.precision = NULL, rsimd.na_check = NULL)
  old_env <- Sys.getenv("RSIMD_IMPL", unset = NA)
  on.exit({
    options(old_opts)
    if (is.na(old_env)) Sys.unsetenv("RSIMD_IMPL") else Sys.setenv(RSIMD_IMPL = old_env)
    simd_use(old_impl)
  })

  Sys.unsetenv("RSIMD_IMPL")
  onload(NULL, "rsimd")
  expect_identical(getOption("rsimd.impl"), "auto")
  expect_identical(.impl_state$requested, "auto")

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
  old_impl <- getOption("rsimd.impl")
  old_opts <- options(rsimd.impl = NULL)
  old_env <- Sys.getenv("RSIMD_IMPL", unset = NA)
  on.exit({
    options(old_opts)
    if (is.na(old_env)) Sys.unsetenv("RSIMD_IMPL") else Sys.setenv(RSIMD_IMPL = old_env)
    simd_use(old_impl)
  })

  Sys.setenv(RSIMD_IMPL = "bogus")
  expect_warning(onload(NULL, "rsimd"), "RSIMD_IMPL: unknown implementation 'bogus'")
  expect_identical(getOption("rsimd.impl"), "auto")

  Sys.unsetenv("RSIMD_IMPL")
  options(rsimd.impl = "rvv")
  expect_warning(onload(NULL, "rsimd"), "option rsimd.impl: implementation 'rvv' is not available")
  expect_identical(getOption("rsimd.impl"), "auto")
  expect_identical(.impl_state$requested, "auto")
})

test_that(".onLoad keeps user-set precision and NA-check options", {
  onload <- get(".onLoad", envir = asNamespace("rsimd"))
  old_opts <- options(rsimd.precision = "compensated", rsimd.na_check = FALSE)
  on.exit(options(old_opts))

  onload(NULL, "rsimd")
  expect_identical(getOption("rsimd.precision"), "compensated")
  expect_identical(getOption("rsimd.na_check"), FALSE)
})

test_that(".onLoad warns about invalid option values and uses the defaults", {
  onload <- get(".onLoad", envir = asNamespace("rsimd"))
  old_opts <- options(
    rsimd.precision = "bogus", rsimd.math_accuracy = "FAST", rsimd.na_check = "yes"
  )
  on.exit(options(old_opts))
  warnings <- character()
  withCallingHandlers(onload(NULL, "rsimd"), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(warnings, c(
    paste(
      "option rsimd.precision: precision mode must be one of \"fast\", \"pairwise\",",
      "\"compensated\"; using \"fast\""
    ),
    paste(
      "option rsimd.math_accuracy: math accuracy mode must be one of \"accurate\",",
      "\"fast\"; using \"accurate\""
    ),
    "option rsimd.na_check: must be TRUE or FALSE; using TRUE"
  ))
  expect_identical(getOption("rsimd.precision"), "fast")
  expect_identical(getOption("rsimd.math_accuracy"), "accurate")
  expect_identical(getOption("rsimd.na_check"), TRUE)
})

test_that("invalid options set before library(rsimd) warn and fall back", {
  skip_if_no_subprocess()
  skip_on_cran()
  child <- function(opts) {
    options(opts)
    warnings <- character()
    withCallingHandlers(library(rsimd), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(
      warnings = warnings,
      options = options()[c("rsimd.impl", "rsimd.precision", "rsimd.math_accuracy", "rsimd.na_check")],
      sum = simd_sum(c(1, 2))
    )
  }
  cases <- list(
    list("rsimd.impl", "bogus", "auto"),
    list("rsimd.precision", "bogus", "fast"),
    list("rsimd.math_accuracy", "bogus", "accurate"),
    list("rsimd.na_check", "bogus", TRUE)
  )
  for (case in cases) {
    res <- callr::r(child, list(structure(list(case[[2]]), names = case[[1]])))
    expect_length(res$warnings, 1)
    expect_match(res$warnings, paste0("^option ", case[[1]], ": "), info = case[[1]])
    expect_identical(res$options[[case[[1]]]], case[[3]], info = case[[1]])
    expect_identical(res$sum, 3, info = case[[1]])
  }
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
