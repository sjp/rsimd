# The bundled SLEEF inline headers: which tiers configure built with them,
# and the SLEEF exp of the selftest_sleef_exp slot on every tier. Elementary
# function semantics are tested with the math functions themselves.

split_list <- function(x) {
  if (identical(x, "")) character() else strsplit(x, ",", fixed = TRUE)[[1]]
}

sleef_tiers <- function() split_list(attr(simd_compiled_tiers(), "sleef"))

# The tier whose kernel fills the slot when `tier` is selected: the tier
# itself if it was built with SLEEF, else the next tier down that was, else
# none.
expected_owner <- function(tier) {
  avail <- simd_available()
  below <- avail[seq(match(tier, avail), length(avail))]
  c(intersect(below, sleef_tiers()), "none")[[1]]
}

# Largest error of y against libm's exp(x) in units of the last place of
# the reference (at least the smallest subnormal), over the elements where
# the reference is finite and not zero; the others must match exactly (NaN
# for NaN). SLEEF's exp is within 1 ULP of the exact value and libm's
# close to 0.5, so they can differ by up to about 1.5 ULP.
max_ulp_error <- function(y, x) {
  ref <- exp(x)
  finite <- is.finite(ref) & ref != 0
  expect_identical(is.nan(y[!finite]), is.nan(ref[!finite]))
  expect_identical(y[!finite & !is.nan(ref)], ref[!finite & !is.nan(ref)])
  if (!any(finite)) {
    return(0)
  }
  max(abs(y[finite] - ref[finite]) / pmax(ulp(ref[finite]), 5e-324))
}

test_that("SLEEF is only reported for compiled tiers that have a header", {
  tiers <- simd_compiled_tiers()
  sleef <- sleef_tiers()
  expect_true(all(sleef %in% split_list(attr(tiers, "configured"))))
  # none uses libm; armv7 NEON has no double-precision vectors.
  expect_false("none" %in% sleef)
  if (simd_cpu_features()$arch == "armv7") {
    expect_length(sleef, 0)
  }
})

test_that("the SLEEF slot runs the tier's own kernel or fills down", {
  res <- for_each_tier(function(tier) {
    list(reported = simd_kernel_tiers()[["selftest_sleef_exp"]], expected = expected_owner(tier))
  })
  for (tier in names(res)) {
    expect_identical(res[[tier]]$reported, res[[tier]]$expected, info = tier)
  }
})

test_that("SLEEF exp is within 2 ULP of libm on every tier", {
  x <- c(
    0, -0, 1, -1, 0.5, -0.5, 1e-300, -1e-300, 5e-324, 1e-8, 20, -20, 100, -100,
    700, 709.78, -708, -745, -746, 710, Inf, -Inf, NaN,
    seq(-30, 30, length.out = 241), with_seed(13L, stats::runif(500, -740, 709))
  )
  res <- for_each_tier(function(tier) .debug_sleef_exp(x))
  for (tier in names(res)) {
    expect_lte(max_ulp_error(res[[tier]], x), 2, label = paste("max ULP error on", tier))
    # exp(0) and exp(-0) are exactly 1 and the infinities are exact.
    expect_identical(res[[tier]][c(1, 2)], c(1, 1), info = tier)
  }
})

test_that("SLEEF exp handles every length, including the vector tails", {
  for (n in edge_lengths()) {
    x <- with_seed(n, stats::runif(n, -50, 50))
    res <- for_each_tier(function(tier) .debug_sleef_exp(x))
    for (tier in names(res)) {
      expect_length(res[[tier]], n)
      expect_lte(max_ulp_error(res[[tier]], x), 2, label = paste("n =", n, "on", tier))
    }
  }
})

test_that("a build without SLEEF takes exp from none on every tier", {
  skip_if(length(sleef_tiers()) > 0, "package built with SLEEF")
  res <- for_each_tier(function(tier) simd_kernel_tiers()[["selftest_sleef_exp"]])
  expect_true(all(unlist(res) == "none"))
  x <- c(-1, 0, 1, 2.5)
  for (tier in names(res)) {
    expect_identical(simd_with_impl(tier, .debug_sleef_exp(x)), exp(x), info = tier)
  }
})

test_that(".debug_sleef_exp() checks its argument", {
  expect_error(.debug_sleef_exp(1L), "double")
  expect_identical(.debug_sleef_exp(double()), double())
})
