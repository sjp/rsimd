test_that(".onLoad sets default options", {
  expect_type(getOption("rsimd.impl"), "character")
  expect_identical(getOption("rsimd.precision"), "fast")
  expect_identical(getOption("rsimd.na_check"), TRUE)
})

test_that(".onLoad defaults rsimd.impl to auto and honours RSIMD_IMPL", {
  onload <- get(".onLoad", envir = asNamespace("rsimd"))
  old_opts <- options(rsimd.impl = NULL, rsimd.precision = NULL, rsimd.na_check = NULL)
  old_env <- Sys.getenv("RSIMD_IMPL", unset = NA)
  on.exit({
    options(old_opts)
    if (is.na(old_env)) Sys.unsetenv("RSIMD_IMPL") else Sys.setenv(RSIMD_IMPL = old_env)
  })

  Sys.unsetenv("RSIMD_IMPL")
  onload(NULL, "rsimd")
  expect_identical(getOption("rsimd.impl"), "auto")

  options(rsimd.impl = "none")
  onload(NULL, "rsimd")
  expect_identical(getOption("rsimd.impl"), "none")

  Sys.setenv(RSIMD_IMPL = "neon")
  onload(NULL, "rsimd")
  expect_identical(getOption("rsimd.impl"), "neon")
})

test_that(".onLoad keeps user-set precision and NA-check options", {
  onload <- get(".onLoad", envir = asNamespace("rsimd"))
  old_opts <- options(rsimd.precision = "compensated", rsimd.na_check = FALSE)
  on.exit(options(old_opts))

  onload(NULL, "rsimd")
  expect_identical(getOption("rsimd.precision"), "compensated")
  expect_identical(getOption("rsimd.na_check"), FALSE)
})
