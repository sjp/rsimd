test_that("simd_version() returns a single string", {
  v <- simd_version()
  expect_type(v, "character")
  expect_length(v, 1L)
  expect_false(is.na(v))
})
