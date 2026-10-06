# Fast mode (option rsimd.math_accuracy = "fast"): the setter, the 3.5-ULP
# bounds against the none tier, and missing values, special values, zero
# signs and warnings identical to accurate mode.

# Runs expr with option rsimd.math_accuracy set to mode.
with_accuracy <- function(mode, expr) {
  old <- options(rsimd.math_accuracy = mode)
  on.exit(options(old))
  force(expr)
}

# The functions that use a 3.5-ULP SLEEF variant in fast mode; the others
# give the same results in both modes.
fast_math1 <- c(
  "log", "log2", "cbrt", "sin", "cos", "tan", "asin", "acos", "atan", "sinpi", "cospi",
  "tanpi", "sinh", "cosh", "tanh"
)

# Every function's result in both modes on every tier, as list(accurate =
# , fast = ) of lists named by tier, each with the warnings it gave.
both_modes <- function(f, ...) {
  run <- function() {
    w <- character(0)
    value <- withCallingHandlers(f(...), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = w)
  }
  list(
    accurate = with_accuracy("accurate", with_each_tier(run)),
    fast = with_accuracy("fast", with_each_tier(run))
  )
}

# Expects the elements of fast where accurate is missing or infinite, or
# the none tier's result is exactly zero, to be bit-identical to
# accurate's (NA vs NaN and the sign of zero included), and both modes to
# give the same warnings. (A subnormal result may round to a zero of
# either sign in either mode: SLEEF's tan_u10 gives +0 for -5e-324.)
expect_same_specials <- function(res, label) {
  none <- res$accurate[["none"]]$value
  for (tier in names(res$accurate)) {
    a <- res$accurate[[tier]]
    f <- res$fast[[tier]]
    info <- paste(label, "on", tier)
    sp <- !is.finite(a$value) | none %in% 0
    expect_true(identical(f$value[sp], a$value[sp], num.eq = FALSE), info = info)
    expect_identical(f$warnings, a$warnings, info = info)
  }
}

test_that("simd_math_accuracy() gets, sets and validates the mode", {
  old <- options(rsimd.math_accuracy = NULL)
  on.exit(options(old))
  expect_identical(simd_math_accuracy(), "accurate")
  expect_invisible(prev <- simd_math_accuracy("fast"))
  expect_identical(prev, "accurate")
  expect_identical(simd_math_accuracy(), "fast")
  expect_identical(getOption("rsimd.math_accuracy"), "fast")
  expect_identical(simd_math_accuracy("accurate"), "fast")
  # Setting the option directly is honoured.
  options(rsimd.math_accuracy = "fast")
  expect_identical(simd_math_accuracy(), "fast")

  for (bad in list("bogus", NA_character_, c("fast", "accurate"), 1, NULL, "FAST")) {
    expect_error(simd_math_accuracy(bad), "math accuracy mode must be one of",
      info = deparse(bad)
    )
  }
  expect_identical(simd_math_accuracy(), "fast")
  options(rsimd.math_accuracy = "bogus")
  expect_error(simd_math_accuracy(), "invalid option rsimd.math_accuracy")
  expect_error(simd_sin(1), "invalid option rsimd.math_accuracy")
  expect_error(simd_sincos(1), "invalid option rsimd.math_accuracy")
  expect_error(simd_exp(1i), "invalid option rsimd.math_accuracy")
  expect_error(simd_max_abs(1i), "invalid option rsimd.math_accuracy")
  # Functions the mode does not affect do not read it.
  expect_identical(simd_next_up(0), 4.9406564584124654e-324)
})

test_that("fast unary functions are within 3.5 ULP of SLEEF's bound of the none tier", {
  tab <- math1_table()
  for (name in names(tab)) {
    spec <- tab[[name]]
    x <- math_inputs(spec[[3]])
    # The none tier (libm) is the reference in both modes; its own error
    # (the accurate bound less SLEEF's 1 ULP) is added to the 3.5 ULP.
    none_bound <- if (is.null(spec$none)) 2 else spec$none
    bound <- if (name %in% fast_math1) none_bound + 2.5 else none_bound
    with_accuracy("fast", expect_tiers_close(spec[[1]], x,
      ulps = bound, label = name, x = x, abs = trig_ref_err(x, name)
    ))
  }
})

test_that("fast binary functions, log to a base and sincos are within their bounds", {
  x <- c(math_specials(), math_random(1500, -300, 300, TRUE))
  y <- c(rev(math_specials()), math_random(1500, -300, 300, TRUE, seed = 7L))
  px <- suppressWarnings(abs(x) %% 30)
  py <- suppressWarnings(y %% 8)
  with_accuracy("fast", {
    expect_tiers_close(simd_atan2, y, x, ulps = 4.5, label = "atan2")
    expect_tiers_close(simd_hypot, x, y, ulps = 4.5, label = "hypot")
    expect_tiers_close(simd_pow, px, py, ulps = 2, label = "pow")
    expect_tiers_close(simd_log, abs(x), 3, ulps = 4.5, label = "log base 3")
    s <- math_random(3000, -1e3, 1e3)
    expect_tiers_close(function(v) simd_sincos(v)$sin, s,
      ulps = 4.5, label = "sincos sin", abs = trig_ref_err(s, "sin")
    )
    expect_tiers_close(function(v) simd_sincos(v)$cos, s,
      ulps = 4.5, label = "sincos cos", abs = trig_ref_err(s, "cos")
    )
  })
})

test_that("fast mode keeps missing values, infinities, zero signs and warnings", {
  tab <- math1_table()
  for (name in names(tab)) {
    x <- c(math_specials(), math_inputs(tab[[name]][[3]]))
    expect_same_specials(both_modes(tab[[name]][[1]], x), name)
  }
  x <- c(math_specials(), -2, 3, NA, NaN)
  y <- rev(x)
  for (name in c("pow", "atan2", "hypot")) {
    f <- get(paste0("simd_", name))
    expect_same_specials(both_modes(f, x, y), name)
    expect_same_specials(both_modes(f, x, 2), paste(name, "scalar y"))
  }
  expect_same_specials(both_modes(simd_log, x, 3), "log base 3")
  expect_same_specials(both_modes(simd_log, x, -1), "log base -1")
  expect_same_specials(both_modes(function(v) simd_sincos(v)$sin, x), "sincos sin")
  expect_same_specials(both_modes(function(v) simd_sincos(v)$cos, x), "sincos cos")
  # Integer and logical input, NA included.
  expect_same_specials(both_modes(simd_log, c(-3:3, NA)), "log integer")
  expect_same_specials(both_modes(simd_asin, c(TRUE, FALSE, NA)), "asin logical")
})

test_that("fast mode warns exactly where accurate mode does", {
  cases <- list(
    log = -1, log2 = -1, sin = Inf, cos = -Inf, tan = Inf, asin = 2, acos = -2,
    sinpi = Inf, cospi = Inf, tanpi = 0.5, cbrt = NaN, atan = NA
  )
  for (name in names(cases)) {
    f <- get(paste0("simd_", name))
    warns <- !name %in% c("cbrt", "atan")
    for (mode in c("accurate", "fast")) {
      for (tier in tiers_to_test()) {
        info <- paste(name, mode, tier)
        with_accuracy(mode, simd_with_impl(tier, {
          if (warns) {
            expect_warning(f(c(1, cases[[name]])), "NaNs produced", info = info)
          } else {
            expect_no_warning(f(c(1, cases[[name]])))
          }
        }))
      }
    }
  }
})

test_that("functions without a fast variant are identical in both modes", {
  tab <- math1_table()
  for (name in setdiff(names(tab), fast_math1)) {
    x <- math_inputs(tab[[name]][[3]])
    res <- both_modes(tab[[name]][[1]], x)
    expect_identical(res$fast, res$accurate, info = name)
  }
  x <- c(math_specials(), math_random(1500, -300, 300, TRUE))
  res <- both_modes(simd_pow, suppressWarnings(abs(x) %% 30), suppressWarnings(x %% 8))
  expect_identical(res$fast, res$accurate, info = "pow")
  res <- both_modes(simd_log, abs(x), 10)
  expect_identical(res$fast, res$accurate, info = "log base 10")
})

test_that("the ML helpers give the same results in both modes", {
  x <- c(math_random(5000, -50, 50), -800, 800, 0, -0, NA)
  for (f in list(simd_sigmoid, simd_softmax, simd_log_softmax)) {
    res <- both_modes(f, x)
    expect_identical(res$fast, res$accurate)
    res <- both_modes(f, x[!is.na(x)])
    expect_identical(res$fast, res$accurate)
  }
})

test_that("fast mode makes no difference on the none tier", {
  tab <- math1_table()
  for (name in names(tab)) {
    x <- math_inputs(tab[[name]][[3]])
    f <- tab[[name]][[1]]
    a <- simd_with_impl("none", suppressWarnings(f(x)))
    expect_identical(with_accuracy("fast", simd_with_impl("none", suppressWarnings(f(x)))), a,
      info = name
    )
  }
})

test_that("accurate mode is the default", {
  tab <- math1_table()
  for (name in names(tab)) {
    x <- math_inputs(tab[[name]][[3]])
    f <- tab[[name]][[1]]
    unset <- with_accuracy(NULL, with_each_tier(function() suppressWarnings(f(x))))
    set <- with_accuracy("accurate", with_each_tier(function() suppressWarnings(f(x))))
    expect_identical(set, unset, info = name)
  }
})

test_that("simd_vec operands follow the global mode, which they do not pin", {
  x <- math_random(1000, -10, 10)
  for (tier in tiers_to_test()) {
    v <- simd_vec(x, impl = tier)
    want <- with_accuracy("fast", simd_with_impl(tier, simd_sin(x)))
    got <- with_accuracy("fast", sin(v))
    expect_s3_class(got, "simd_vec")
    expect_identical(simd_impl(got), tier)
    expect_identical(.sv_strip(got), want, info = tier)
  }
})
