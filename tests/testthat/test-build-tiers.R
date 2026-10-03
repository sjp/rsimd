split_tiers <- function(x) {
  if (identical(x, "")) character() else strsplit(x, ",", fixed = TRUE)[[1]]
}

test_that("every configured tier object is linked and names itself", {
  tiers <- simd_compiled_tiers()
  expect_type(tiers, "character")
  expect_identical(as.vector(tiers), split_tiers(attr(tiers, "configured")))
  expect_true("none" %in% tiers)
  expect_true(all(tiers %in% names(simd_cpu_tiers())))
  # simd_tiers() order, which is also the enum order.
  expect_identical(as.vector(tiers), intersect(names(simd_cpu_tiers()), tiers))
})

test_that("the platform baseline tier is compiled unless disabled", {
  tiers <- simd_compiled_tiers()
  disabled <- split_tiers(attr(tiers, "disabled"))
  baseline <- switch(simd_cpu_features()$arch, x86_64 = "sse2", aarch64 = "neon", NULL)
  skip_if(is.null(baseline), "no baseline tier on this architecture")
  expect_true(baseline %in% tiers || baseline %in% disabled)
})

test_that("reserved tiers are never compiled", {
  expect_false(any(c("rvv", "wasm128") %in% simd_compiled_tiers()))
})

test_that("tiers for another architecture are never compiled", {
  tiers <- simd_compiled_tiers()
  arch <- simd_cpu_features()$arch
  if (arch %in% c("x86_64", "i686")) {
    expect_false(any(c("neon", "sve", "sve2") %in% tiers))
  }
  if (arch %in% c("aarch64", "armv7")) {
    expect_false(any(c("sse2", "avx2", "avx512") %in% tiers))
  }
})
