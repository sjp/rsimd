# Elementwise ops with every shape of operand: double or integer, vector
# or length 1, in every position. The vector tiers read double vectors and
# scalars with one loop and int32 vectors with another, so each mix of
# shapes is checked against the none tier at lengths with and without a
# tail on every tier.

shape_lengths <- function() c(2, 3, 9, 67)

# Operand k of shape `shape` at length n: a double or int vector with
# special values among ordinary ones, or a double or int scalar (one of a
# few, picked by k, NA among them).
shape_operand <- function(shape, n, k) {
  switch(shape,
    dv = with_seed(n * 7L + k, {
      v <- round(stats::runif(n, -50, 50), 2)
      i <- sample.int(n, n %/% 3)
      v[i] <- sample(edge_doubles(xmax = FALSE), length(i), TRUE)
      v
    }),
    iv = with_seed(n * 7L + k, {
      v <- sample(-60:60, n, TRUE)
      v[sample.int(n, n %/% 4)] <- NA
      v
    }),
    ds = c(2.5, -0, NA, -3, NaN, Inf)[(k %% 6) + 1L],
    is = c(3L, NA, -7L, 0L)[(k %% 4) + 1L]
  )
}

# Expects f to give the none tier's result for every mix of `nargs`
# operand shapes at every length. One expectation naming each case that
# differs.
expect_shapes_like_none <- function(f, nargs, name) {
  shapes <- expand.grid(rep(list(c("dv", "ds", "iv", "is")), nargs), stringsAsFactors = FALSE)
  g <- function(...) suppressWarnings(f(...))
  problems <- character(0)
  for (n in shape_lengths()) {
    for (s in seq_len(nrow(shapes))) {
      args <- lapply(seq_len(nargs), function(k) shape_operand(shapes[s, k], n, s + k))
      expected <- simd_with_impl("none", do.call(g, args))
      res <- with_each_tier(function() do.call(g, args))
      for (tier in names(res)) {
        if (!same_values(res[[tier]], expected)) {
          problems <- c(problems, sprintf(
            "%s(%s) at n = %d, tier %s: %s", name, paste(shapes[s, ], collapse = ", "), n, tier,
            first_difference(res[[tier]], expected)
          ))
        }
      }
    }
  }
  simd_expect(length(problems) == 0L, paste(utils::head(problems, 10L), collapse = "\n"))
}

test_that("binary arithmetic gives the none tier's result for every operand shape", {
  for (name in c(
    "simd_add", "simd_sub", "simd_mul", "simd_div", "simd_idiv", "simd_mod", "simd_pmin",
    "simd_pmax", "simd_pmin_num", "simd_pmax_num", "simd_copysign"
  )) {
    expect_shapes_like_none(get(name), 2L, name)
  }
})

test_that("comparisons and logical ops give the none tier's result for every operand shape", {
  for (name in c(
    "simd_eq", "simd_ne", "simd_lt", "simd_le", "simd_gt", "simd_ge", "simd_and", "simd_or",
    "simd_xor"
  )) {
    expect_shapes_like_none(get(name), 2L, name)
  }
})

test_that("three-operand ops give the none tier's result for every operand shape", {
  for (name in c("simd_fma", "simd_mul_add", "simd_add_mul", "simd_lerp")) {
    expect_shapes_like_none(get(name), 3L, name)
  }
  # lo <= hi wherever both are numbers, so that no tier warns.
  clamp <- function(x, lo, hi) simd_clamp(x, -abs(lo) - 1, abs(hi) + 1)
  expect_shapes_like_none(clamp, 3L, "simd_clamp")
})

test_that("a scalar operand gives what the same value repeated gives, on either side", {
  x <- shape_operand("dv", 67, 1L)
  for (s in list(2.5, -3L, NA, NA_integer_, NaN)) {
    r <- rep(s, length(x))
    expect_tiers_give(simd_with_impl("none", simd_sub(r, x)), simd_sub, s, x)
    expect_tiers_give(simd_with_impl("none", simd_div(x, r)), simd_div, x, s)
    expect_tiers_give(simd_with_impl("none", simd_lt(r, x)), simd_lt, s, x)
    expect_tiers_give(simd_with_impl("none", simd_fma(x, r, r)), simd_fma, x, s, s)
    expect_tiers_give(simd_with_impl("none", simd_fma(r, x, r)), simd_fma, s, x, s)
  }
  expect_tiers_give(na_merged(1 - x, x), simd_sub, 1, x)
  expect_tiers_give(na_merged(2L * x + 1, x), simd_fma, 2L, x, 1)
})
