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

test_that("simd_current() returns a plain tier id", {
  cur <- simd_current()
  expect_identical(attributes(cur), NULL)
  expect_true(cur %in% simd_available())
  expect_identical(capture.output(print(cur)), paste0('[1] "', cur, '"'))
})

test_that("simd_use() selects none and auto", {
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  expect_silent(simd_use("none"))
  expect_identical(as.vector(simd_current()), "none")
  expect_identical(.impl_state$requested, "none")
  expect_identical(getOption("rsimd.impl"), "none")

  expect_identical(simd_use("auto"), "none")
  expect_identical(as.vector(simd_current()), simd_available()[1])
  expect_identical(.impl_state$requested, "auto")

  for (tier in simd_available()) {
    simd_use(tier)
    expect_identical(as.vector(simd_current()), tier)
  }
})

test_that("simd_use() returns the previous request invisibly", {
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  simd_use("auto")
  expect_invisible(simd_use("none"))
  expect_identical(simd_use("auto"), "none")
  expect_identical(simd_use("auto"), "auto")
})

test_that("simd_use() without impl returns the current request", {
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  simd_use("auto")
  expect_visible(simd_use())
  expect_identical(simd_use(), "auto")
  simd_use("none")
  expect_identical(simd_use(), "none")
  # A value set through the option is the request too.
  options(rsimd.impl = "auto")
  expect_identical(simd_use(), "auto")
  expect_identical(as.vector(simd_current()), simd_available()[1])
  # An invalid option is reported as by the compute functions.
  options(rsimd.impl = "bogus")
  expect_error(simd_use(), "invalid option rsimd.impl: unknown implementation 'bogus'",
    fixed = TRUE
  )
  expect_error(simd_use(), "; reset it with simd_use()", fixed = TRUE)
  expect_identical(simd_use("none"), "auto")
  expect_identical(simd_use(), "none")
})

test_that("simd_options() reports every option", {
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  old <- options(
    rsimd.precision = "pairwise", rsimd.math_accuracy = "fast", rsimd.na_check = FALSE
  )
  on.exit(options(old), add = TRUE)
  simd_use("auto")
  opts <- simd_options()
  expect_s3_class(opts, "simd_options")
  expect_identical(unclass(opts), list(
    impl = "auto", current = simd_available()[1], precision = "pairwise",
    math_accuracy = "fast", na_check = FALSE
  ))
  expect_output(print(opts), "impl: +auto \\(running [a-z0-9]+\\)")
  expect_output(print(opts), "precision: +pairwise")
  expect_output(print(opts), "math_accuracy: +fast")
  expect_output(print(opts), "na_check: +FALSE")
  expect_invisible(print(opts))
  # A tier asked for by name is shown once.
  simd_use("none")
  expect_identical(simd_options()$current, "none")
  expect_output(print(simd_options()), "impl: +none\n")
  # An invalid option is an error, as for the function that reads it.
  options(rsimd.precision = "bogus")
  expect_error(simd_options(), "invalid option rsimd.precision")
  options(rsimd.precision = "fast", rsimd.na_check = "yes")
  expect_error(simd_options(), "; reset it with simd_na_check()", fixed = TRUE)
})

test_that("simd_use() rejects an unavailable tier and keeps the selection", {
  old_impl <- getOption("rsimd.impl")
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
  old_impl <- getOption("rsimd.impl")
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
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  simd_use("auto")
  expect_identical(as.vector(simd_with_impl("none", simd_current())), "none")
  expect_identical(.impl_state$requested, "auto")
  expect_identical(getOption("rsimd.impl"), "auto")

  expect_error(simd_with_impl("none", stop("boom")), "boom")
  expect_identical(.impl_state$requested, "auto")
  expect_identical(as.vector(simd_current()), simd_available()[1])

  expect_error(simd_with_impl("bogus", 1), "unknown implementation")
  expect_identical(.impl_state$requested, "auto")

  simd_use("none")
  nested <- simd_with_impl("auto", simd_with_impl("none", simd_current()))
  expect_identical(as.vector(nested), "none")
  expect_identical(.impl_state$requested, "none")
})

test_that("setting the rsimd.impl option directly is honoured", {
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  old <- options(rsimd.impl = "auto")
  on.exit(options(old), add = TRUE)
  simd_use("auto")

  options(rsimd.impl = "none")
  expect_identical(as.vector(simd_current()), "none")
  expect_identical(.impl_state$requested, "none")

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

test_that("compute functions honour rsimd.impl set directly, without simd_current()", {
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  old <- options(rsimd.impl = "auto")
  on.exit(options(old), add = TRUE)
  simd_use("auto")
  best <- simd_available()[1]

  # The C side syncs at the start of each call (.debug_active() is one).
  options(rsimd.impl = "none")
  expect_identical(.debug_active(1), "none")
  options(rsimd.impl = "auto")
  expect_identical(.debug_active(1), best)

  # Back to the value last synced, after a simd_use() in between.
  options(rsimd.impl = "none")
  expect_identical(.debug_active(1), "none")
  simd_use("auto")
  options(rsimd.impl = "none")
  expect_identical(.debug_active(1), "none")

  # An unset option means "auto".
  options(rsimd.impl = NULL)
  expect_identical(.debug_active(1), best)

  # A value that cannot be selected is an error naming the option, on every
  # call until it is fixed; the selection is unchanged.
  simd_use("none")
  options(rsimd.impl = "bogus")
  expect_error(simd_sum(1), "invalid option rsimd.impl: unknown implementation 'bogus'",
    fixed = TRUE
  )
  expect_error(simd_add(1, 2), "invalid option rsimd.impl: unknown implementation 'bogus'",
    fixed = TRUE
  )
  options(rsimd.impl = 1)
  expect_error(simd_sum(1), "invalid option rsimd.impl: 'impl' must be a single string",
    fixed = TRUE
  )
  options(rsimd.impl = "rvv")
  expect_error(simd_sum(1), "invalid option rsimd.impl: implementation 'rvv' is not available")
  # The error says how to recover, and simd_use() does.
  expect_error(simd_sum(1), "; reset it with simd_use()", fixed = TRUE)
  expect_identical(simd_use("auto"), "none")
  expect_identical(simd_sum(c(1, 2)), 3)
  expect_identical(.debug_active(1), best)
})

test_that("an option set in a calling handler of an rsimd call is honoured inside it", {
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  old <- options(rsimd.impl = "auto")
  on.exit(options(old), add = TRUE)
  simd_use("auto")
  best <- simd_available()[1]

  # A warning handler runs once the call has finished.
  inside <- NULL
  res <- withCallingHandlers(simd_as_integer(c(1e10, 2)), warning = function(w) {
    options(rsimd.impl = "none")
    inside <<- .debug_active(1)
    invokeRestart("muffleWarning")
  })
  expect_identical(res, c(NA, 2L))
  expect_identical(inside, "none")
  expect_identical(.debug_active(1), "none")

  # A handler of the invalid-option error runs after the sync, so a call
  # it makes syncs again.
  options(rsimd.impl = "auto")
  expect_identical(.debug_active(1), best)
  options(rsimd.impl = "bogus")
  inside <- NULL
  expect_error(
    withCallingHandlers(simd_sum(1), error = function(e) {
      options(rsimd.impl = "none")
      inside <<- list(.debug_active(1), simd_sum(c(1, 2)))
    }),
    "invalid option rsimd.impl: unknown implementation 'bogus'"
  )
  expect_identical(inside, list("none", 3))
  expect_identical(.debug_active(1), "none")
  # The error still names the user's call.
  options(rsimd.impl = "bogus")
  err <- tryCatch(simd_sum(1), error = identity)
  expect_identical(conditionCall(err), quote(simd_sum(1)))
  options(rsimd.impl = "auto")
  expect_identical(.debug_active(1), best)
})

test_that("an interrupt between chunks leaves the next call's state clean", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("callr")
  rs <- callr::r_session$new()
  on.exit(rs$close(), add = TRUE)
  rs$call(function() {
    tiers <- rsimd::simd_available()
    # A pinned operand switches the table for the call; the interrupt
    # leaves the call with it still switched.
    x <- rsimd::simd_vec(as.double(seq_len(5e7)), impl = tiers[[length(tiers)]])
    tryCatch(repeat rsimd::simd_sum(x), interrupt = function(e) "interrupted")
  })
  Sys.sleep(1)
  rs$interrupt()
  expect_identical(rs$poll_process(10000), "ready")
  expect_identical(rs$read()$result, "interrupted")
  res <- rs$run(function() {
    list(
      selected = rsimd::simd_current(), active = rsimd:::.debug_active(1),
      sum = rsimd::simd_sum(c(1, 2, NA), na.rm = TRUE),
      warning = tryCatch(rsimd::simd_as_integer(1e10), warning = conditionMessage)
    )
  })
  expect_identical(res$active, res$selected)
  expect_identical(res$sum, 3)
  expect_identical(res$warning, "NAs introduced by coercion to integer range")
})

test_that("forked workers inherit the selection; callr sessions start from RSIMD_IMPL", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("callr")
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  child <- function() list(rsimd::simd_current(), rsimd:::.debug_active(1))

  simd_use("none")
  res <- parallel::mclapply(1:2, function(i) child(), mc.cores = 2)
  expect_identical(res, rep(list(list("none", "none")), 2))
  # A selection made in the worker stays there.
  res <- parallel::mclapply(1:2, function(i) {
    simd_use("auto")
    child()
  }, mc.cores = 2)
  best <- simd_available()[1]
  expect_identical(res, rep(list(list(best, best)), 2))
  expect_identical(.debug_active(1), "none")

  # A callr session is a fresh process: RSIMD_IMPL, not the parent's
  # selection.
  env <- function(impl) c(callr::rcmd_safe_env(), RSIMD_IMPL = impl)
  expect_identical(callr::r(child, env = env("")), list(best, best))
  simd_use("auto")
  expect_identical(callr::r(child, env = env("none")), list("none", "none"))
})

test_that("simd_with_impl() restores the selection when a compute call errors", {
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  simd_use("auto")
  expect_error(simd_with_impl("none", simd_add(1:3, 1:2)), "lengths")
  expect_identical(.debug_active(1), simd_available()[1])
  expect_identical(getOption("rsimd.impl"), "auto")
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
      current = current, requested = getOption("rsimd.impl"),
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
  old_impl <- getOption("rsimd.impl")
  on.exit(simd_use(old_impl), add = TRUE)
  simd_use("none")
  expect_identical(simd_kernel_tiers(), simd_kernel_tiers("none"))
  expect_error(simd_kernel_tiers("rvv"), "not available")
  expect_error(simd_kernel_tiers("bogus"), "not available")
  expect_error(simd_kernel_tiers(1), "single string")
})
