# Implementation selection: simd_tiers(), simd_available(), simd_current(),
# simd_use(), simd_with_impl() and the dispatcher's fill-down.

preference <- function(arch) {
  switch(arch,
    x86_64 = ,
    i686 = c("avx512", "avx2", "sse2", "none"),
    aarch64 = ,
    armv7 = c("sve2", "sve", "neon", "none"),
    "none"
  )
}

test_that("simd_tiers() lists every tier id in enum order", {
  expect_identical(
    simd_tiers(),
    c("none", "sse2", "avx2", "avx512", "neon", "sve", "sve2", "rvv", "wasm128")
  )
  expect_identical(simd_tiers(), names(simd_cpu_tiers()))
})

test_that("simd_available() is the compiled, supported tiers in preference order", {
  avail <- simd_available()
  compiled <- as.vector(simd_compiled_tiers())
  cpu <- names(which(simd_cpu_tiers()))
  pref <- preference(simd_cpu_features()$arch)
  expect_identical(avail, pref[pref %in% compiled & pref %in% cpu])
  expect_identical(avail[length(avail)], "none")
  expect_false(any(c("rvv", "wasm128") %in% avail))
})

test_that("simd_available() matches the hardware on known machines", {
  f <- simd_cpu_features()
  compiled <- simd_compiled_tiers()
  has <- f$features
  neon_only <- f$arch == "aarch64" && has[["neon"]] && !has[["sve"]]
  avx2_only <- f$arch == "x86_64" && has[["avx2"]] && has[["fma"]] && !has[["avx512f"]]
  if (neon_only && "neon" %in% compiled) {
    expect_identical(simd_available(), c("neon", "none"))
  }
  if (avx2_only && all(c("sse2", "avx2") %in% compiled)) {
    expect_identical(simd_available(), c("avx2", "sse2", "none"))
  }
  expect_true(TRUE)
})

test_that("simd_use() selects none and auto", {
  old_impl <- attr(simd_current(), "requested")
  on.exit(simd_use(old_impl), add = TRUE)
  expect_silent(simd_use("none"))
  expect_identical(as.vector(simd_current()), "none")
  expect_identical(attr(simd_current(), "requested"), "none")
  expect_identical(getOption("rsimd.impl"), "none")

  expect_identical(simd_use("auto"), "none")
  expect_identical(as.vector(simd_current()), simd_available()[1])
  expect_identical(attr(simd_current(), "requested"), "auto")

  for (tier in simd_available()) {
    simd_use(tier)
    expect_identical(as.vector(simd_current()), tier)
  }
})

test_that("simd_use() returns the previous request invisibly", {
  old_impl <- attr(simd_current(), "requested")
  on.exit(simd_use(old_impl), add = TRUE)
  simd_use("auto")
  expect_invisible(simd_use("none"))
  expect_identical(simd_use("auto"), "none")
  expect_identical(simd_use("auto"), "auto")
})

test_that("simd_use() rejects an unavailable tier and keeps the selection", {
  old_impl <- attr(simd_current(), "requested")
  on.exit(simd_use(old_impl), add = TRUE)
  simd_use("none")
  avail <- paste(simd_available(), collapse = ", ")
  unavailable <- setdiff(simd_tiers(), simd_available())
  for (tier in unavailable) {
    expect_error(
      simd_use(tier),
      paste0("implementation '", tier, "' is not available on this machine; available: ", avail),
      fixed = TRUE
    )
  }
  expect_identical(as.vector(simd_current()), "none")
  expect_identical(getOption("rsimd.impl"), "none")
  if (!"avx512" %in% simd_available()) {
    expect_error(simd_use("avx512"), "'avx512' is not available on this machine", fixed = TRUE)
  }
})

test_that("simd_use() rejects unknown names, listing the tiers", {
  old_impl <- attr(simd_current(), "requested")
  on.exit(simd_use(old_impl), add = TRUE)
  simd_use("none")
  expect_error(
    simd_use("bogus"),
    paste0(
      "unknown implementation 'bogus'; use \"auto\" or one of: ",
      paste(simd_tiers(), collapse = ", ")
    ),
    fixed = TRUE
  )
  expect_error(simd_use("AUTO"), "unknown implementation")
  expect_error(simd_use(NA_character_), "single string")
  expect_error(simd_use(c("none", "auto")), "single string")
  expect_error(simd_use(1), "single string")
  expect_error(simd_use(NULL), "single string")
  expect_identical(as.vector(simd_current()), "none")
})

test_that("simd_with_impl() selects temporarily and restores", {
  old_impl <- attr(simd_current(), "requested")
  on.exit(simd_use(old_impl), add = TRUE)
  simd_use("auto")
  expect_identical(as.vector(simd_with_impl("none", simd_current())), "none")
  expect_identical(attr(simd_current(), "requested"), "auto")
  expect_identical(getOption("rsimd.impl"), "auto")

  expect_error(simd_with_impl("none", stop("boom")), "boom")
  expect_identical(attr(simd_current(), "requested"), "auto")
  expect_identical(as.vector(simd_current()), simd_available()[1])

  expect_error(simd_with_impl("bogus", 1), "unknown implementation")
  expect_identical(attr(simd_current(), "requested"), "auto")

  simd_use("none")
  nested <- simd_with_impl("auto", simd_with_impl("none", simd_current()))
  expect_identical(as.vector(nested), "none")
  expect_identical(attr(simd_current(), "requested"), "none")
})

test_that("setting the rsimd.impl option directly is honoured", {
  old_impl <- attr(simd_current(), "requested")
  on.exit(simd_use(old_impl), add = TRUE)
  old <- options(rsimd.impl = "auto")
  on.exit(options(old), add = TRUE)
  simd_use("auto")

  options(rsimd.impl = "none")
  expect_identical(as.vector(simd_current()), "none")
  expect_identical(attr(simd_current(), "requested"), "none")

  options(rsimd.impl = "auto")
  expect_identical(as.vector(simd_current()), simd_available()[1])

  # simd_with_impl() restores a value set through the option.
  options(rsimd.impl = "none")
  simd_with_impl("auto", NULL)
  expect_identical(getOption("rsimd.impl"), "none")
  expect_identical(as.vector(simd_current()), "none")

  options(rsimd.impl = "bogus")
  expect_error(simd_current(), "unknown implementation 'bogus'", fixed = TRUE)
  options(rsimd.impl = "rvv")
  expect_error(simd_current(), "implementation 'rvv' is not available on this machine")
  options(rsimd.impl = "none")
  expect_identical(as.vector(simd_current()), "none")
})

test_that("RSIMD_IMPL and rsimd.impl initialise a fresh session", {
  skip_if_no_subprocess()
  skip_on_cran()
  child <- function() {
    warnings <- character()
    current <- withCallingHandlers(
      rsimd::simd_current(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(
      current = as.vector(current), requested = attr(current, "requested"),
      available = rsimd::simd_available(), warnings = warnings
    )
  }
  env <- function(impl) c(callr::rcmd_safe_env(), RSIMD_IMPL = impl)

  res <- callr::r(child, env = env("none"))
  expect_identical(res$current, "none")
  expect_length(res$warnings, 0)

  res <- callr::r(child, env = env("bogus"))
  expect_identical(res$current, res$available[1])
  expect_identical(res$requested, "auto")
  expect_length(res$warnings, 1)
  expect_match(res$warnings, "RSIMD_IMPL: unknown implementation 'bogus'", fixed = TRUE)

  res <- callr::r(child, env = env("rvv"))
  expect_identical(res$current, res$available[1])
  expect_match(res$warnings, "implementation 'rvv' is not available", fixed = TRUE)

  # An option set before loading wins over the environment variable.
  res <- callr::r(function(child) {
    options(rsimd.impl = "none")
    child()
  }, list(child), env = env("auto"))
  expect_identical(res$current, "none")
  expect_identical(res$requested, "none")
})

test_that("every slot of an available tier resolves to an available tier", {
  for (tier in simd_available()) {
    k <- simd_kernel_tiers(tier)
    expect_type(k, "character")
    expect_identical(names(k)[1:2], c("tier_name", "fill_probe"))
    expect_false(anyNA(names(k)) || anyDuplicated(names(k)) > 0)
    expect_true(all(k %in% simd_available()))
    # A slot never runs a kernel of a better tier than the one selected.
    avail <- simd_available()
    expect_true(all(match(k, avail) >= match(tier, avail)))
  }
  k <- simd_kernel_tiers("none")
  expect_identical(unname(k), rep("none", length(k)))
})

test_that("an empty slot is filled from below: fill_probe always runs none", {
  hole <- attr(simd_compiled_tiers(), "test_hole")
  res <- for_each_tier(function(tier) {
    list(reported = simd_kernel_tiers(), called = simd_probe_slots())
  })
  for (tier in names(res)) {
    expect_identical(res[[tier]]$reported[["fill_probe"]], "none")
    called <- res[[tier]]$called
    expect_identical(called, res[[tier]]$reported[names(called)])
    if (hole != "tier_name") {
      expect_identical(res[[tier]]$reported[["tier_name"]], tier)
    }
  }
})

test_that("a RSIMD_TEST_HOLE slot is filled from the next tier down", {
  hole <- attr(simd_compiled_tiers(), "test_hole")
  skip_if(hole == "", "package not built with RSIMD_TEST_HOLE")
  # The hole is in every tier but none, so the search passes over every
  # lower tier, which lacks the kernel too, and ends at none.
  res <- for_each_tier(function(tier) {
    list(reported = simd_kernel_tiers()[[hole]], called = simd_probe_slots()[[hole]])
  })
  for (tier in names(res)) {
    expect_identical(res[[tier]]$reported, "none")
    expect_identical(res[[tier]]$called, "none")
  }
})

test_that("simd_kernel_tiers() defaults to the active tier and checks its argument", {
  old_impl <- attr(simd_current(), "requested")
  on.exit(simd_use(old_impl), add = TRUE)
  simd_use("none")
  expect_identical(simd_kernel_tiers(), simd_kernel_tiers("none"))
  expect_error(simd_kernel_tiers("rvv"), "not available")
  expect_error(simd_kernel_tiers("bogus"), "not available")
  expect_error(simd_kernel_tiers(1), "single string")
})
