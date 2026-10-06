# The bundled SLEEF inline headers: which tiers configure built with them,
# and which tier's elementary-function kernels each tier uses. Elementary
# function semantics are tested in test-math-*.R.

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

test_that("the math slots run the tier's own kernels or fill down", {
  res <- for_each_tier(function(tier) {
    list(
      reported = simd_kernel_tiers()[c("math1_f64", "math2_f64", "sincos_f64", "softmax_f64")],
      expected = expected_owner(tier)
    )
  })
  for (tier in names(res)) {
    expect_true(all(res[[tier]]$reported == res[[tier]]$expected), info = tier)
  }
})

test_that("SLEEF exp is within 2 ULP of libm on every tier, at every length", {
  x <- c(
    0, -0, 1, -1, 0.5, -0.5, 1e-300, -1e-300, 5e-324, 1e-8, 20, -20, 100, -100,
    700, 709.78, -708, -745, -746, 710, Inf, -Inf, NaN,
    seq(-30, 30, length.out = 241), with_seed(13L, stats::runif(500, -740, 709))
  )
  ref <- expect_tiers_close(simd_exp, x, ulps = 2)
  expect_identical(ref, exp(x))
  for (n in sweep_lengths()) {
    expect_tiers_close(simd_exp, with_seed(n, stats::runif(n, -50, 50)), ulps = 2)
  }
})

test_that("a build without SLEEF takes the math kernels from none on every tier", {
  skip_if(length(sleef_tiers()) > 0, "package built with SLEEF")
  res <- for_each_tier(function(tier) simd_kernel_tiers()[c("math1_f64", "softmax_f64")])
  expect_true(all(unlist(res) == "none"))
  x <- c(-1, 0, 1, 2.5)
  for (tier in names(res)) {
    expect_identical(simd_with_impl(tier, simd_exp(x)), exp(x), info = tier)
    expect_identical(simd_with_impl(tier, simd_softmax(x)), simd_with_impl("none", simd_softmax(x)),
      info = tier
    )
  }
})
