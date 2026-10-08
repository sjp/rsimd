# Complex multiply, divide, Mod, Arg, prod, mean, cumsum, cumprod and the
# comparisons, on every tier, against base R.

# Complex numbers with independent random parts: rand_vec() doubles (NA,
# NaN and +-Inf included unless the fractions are 0) and some signed zeros;
# `wide` spreads the magnitudes over the whole exponent range.
rand_z <- function(n, seed, wide = FALSE, ...) {
  part <- function(s) {
    x <- rand_vec("double", n, seed = s, ...)
    if (wide) {
      x <- with_seed(s + 7L, x * 10^stats::runif(n, -300, 300))
    }
    with_seed(s + 13L, {
      k <- stats::rbinom(1, n, 0.05)
      x[sample.int(n, k)] <- rep_len(c(0, -0), k)
    })
    x
  }
  complex(real = part(seed), imaginary = part(seed + 1000L))
}

# Every combination of `values` for the four parts of x and y.
grid_pairs <- function(values) {
  g <- expand.grid(a = values, b = values, c = values, d = values)
  list(x = complex(real = g$a, imaginary = g$b), y = complex(real = g$c, imaginary = g$d))
}

lengths_c <- c(0:18, 31:35, 63:67, 511:513, 4095:4097)
annex_g <- c(1.5, 0, -0, Inf, -Inf, NaN, NA)
extremes <- c(1e300, -1e300, 1e-300, -1e-300, 5e-324, .Machine$double.xmax, 2.5)

# fma(x, y, z) exactly (the none tier calls C's fma).
fma <- function(x, y, z) simd_with_impl("none", simd_fma(x, y, z))

# x * y in the rounding variants vr and vi (0 unfused, 1 FMA1, 2 FMA2).
cmul_ref <- function(z, w, vr, vi) {
  a <- Re(z)
  b <- Im(z)
  c <- Re(w)
  d <- Im(w)
  re <- switch(vr + 1L,
    a * c - b * d,
    fma(a, c, -(b * d)),
    fma(-b, d, a * c)
  )
  im <- switch(vi + 1L,
    a * d + b * c,
    fma(b, c, a * d),
    fma(a, d, b * c)
  )
  complex(real = re, imaginary = im)
}

# x / y as libgcc's __divdc3 (GCC 12 and later) computes it before its
# recovery of infinities, fused or not.
cdiv_ref <- function(z, w, fused) {
  f2 <- function(x, y, z) if (fused) fma(x, y, z) else x * y + z
  a <- Re(z)
  b <- Im(z)
  c <- Re(w)
  d <- Im(w)
  dmin <- .Machine$double.xmin
  rmax2 <- .Machine$double.xmax / 2 * .Machine$double.eps
  sw <- abs(c) < abs(d)
  big <- ifelse(sw, d, c)
  small <- ifelse(sw, c, d)
  abig <- abs(big)
  tiny <- abig < rmax2 & ((abs(a) < dmin & abs(b) < rmax2) | (abs(b) < dmin & abs(a) < rmax2))
  f <- ifelse(abig >= .Machine$double.xmax / 2, 0.5,
    ifelse(abig < .Machine$double.eps | tiny, 1 / .Machine$double.eps, 1)
  )
  a <- a * f
  b <- b * f
  big <- big * f
  small <- small * f
  u <- ifelse(sw, a, b)
  v <- ifelse(sw, b, a)
  s <- ifelse(sw, b, -a)
  t <- ifelse(sw, -a, b)
  ratio <- small / big
  denom <- f2(small, ratio, big)
  alt <- !(abs(ratio) > dmin)
  x <- ifelse(alt, f2(small, u / big, v), f2(u, ratio, v))
  y <- ifelse(alt, f2(small, s / big, t), f2(s, ratio, t))
  complex(real = x / denom, imaginary = y / denom)
}

# Runs code with the complex arithmetic variants set to `codes`.
with_variants <- function(codes, code) {
  old <- .debug_c128_variants(codes)
  on.exit(.debug_c128_variants(old))
  force(code)
}

test_that("the load-time probe picks a variant for each operation", {
  v <- simd_complex_variants()
  expect_identical(names(v), c("mul", "div", "cumprod"))
  expect_true(all(nzchar(v)))
  # GCC on aarch64 contracts by default, and libgcc 12+ divides.
  if (R.version$arch == "aarch64" && Sys.info()[["sysname"]] == "Linux") {
    expect_identical(v[["mul"]], "fma(a, c, -bd), fma(b, c, ad)")
    expect_identical(v[["div"]], "libgcc, fused")
    expect_identical(v[["cumprod"]], "fma(-b, d, ac), fma(b, c, ad)")
  }
})

test_that("mul and div are identical to base R on every tier", {
  for (n in lengths_c) {
    x <- rand_z(n, seed = n + 1L)
    y <- rand_z(n, seed = n + 2L)
    expect_tiers_give(x * y, simd_mul, x, y)
    expect_tiers_give(x / y, simd_div, x, y)
  }
  x <- rand_z(5000, seed = 3L, wide = TRUE)
  y <- rand_z(5000, seed = 4L, wide = TRUE)
  expect_tiers_give(x * y, simd_mul, x, y)
  # Where the operands hold both NA and NaN, whether base R's quotient keeps
  # the NA or the NaN depends on the C library's division (macOS's differs
  # from libgcc's), so those elements are compared only as missing values.
  parts <- cbind(Re(x), Im(x), Re(y), Im(y))
  mixed <- rowSums(is.na(parts) & !is.nan(parts)) > 0 & rowSums(is.nan(parts)) > 0
  expect_tiers_give(x[!mixed] / y[!mixed], simd_div, x[!mixed], y[!mixed])
  for (tier in tiers_to_test()) {
    q <- simd_with_impl(tier, simd_div(x[mixed], y[mixed]))
    expect_identical(is.na(Re(q)) | is.na(Im(q)), rep(TRUE, sum(mixed)), info = tier)
  }
})

test_that("Annex G special values and extreme magnitudes match base R", {
  for (values in list(annex_g, extremes, c(0, -0, Inf, NA, 1e-300, 1e300))) {
    p <- grid_pairs(values)
    expect_tiers_give(p$x * p$y, simd_mul, p$x, p$y)
    expect_tiers_give(p$x / p$y, simd_div, p$x, p$y)
  }
  expect_identical(simd_div(1, 0 + 0i), 1 / (0 + 0i))
  expect_identical(
    simd_mul(complex(real = Inf, imaginary = NaN), 1 + 1i),
    complex(real = Inf, imaginary = NaN) * (1 + 1i)
  )
})

test_that("scalars broadcast on either side, and other types are promoted", {
  x <- rand_z(1001, seed = 5L)
  for (s in list(2 - 3i, complex(real = NA, imaginary = 1), complex(real = 0, imaginary = -0))) {
    expect_tiers_give(x * s, simd_mul, x, s)
    expect_tiers_give(s * x, simd_mul, s, x)
    expect_tiers_give(x / s, simd_div, x, s)
    expect_tiers_give(s / x, simd_div, s, x)
  }
  d <- rand_vec("double", 1001, seed = 6L)
  expect_tiers_give(x * d, simd_mul, x, d)
  expect_tiers_give(d / x, simd_div, d, x)
  expect_identical(simd_div(x, 3L), x / 3L)
  expect_identical(simd_mul(TRUE, x), TRUE * x)
  expect_identical(simd_mul(complex(0), 1i), complex(0))
  expect_error(simd_mul(1i, structure(0, class = "integer64")),
    "simd_mul() cannot combine complex and integer64 operands",
    fixed = TRUE
  )
  expect_error(simd_div(as.raw(1), 1i), "non-numeric argument to binary operator", fixed = TRUE)
})

test_that("every rounding variant's kernels compute its formula", {
  x <- rand_z(1003, seed = 7L, na_frac = 0, nan_frac = 0, inf_frac = 0)
  y <- rand_z(1003, seed = 8L, na_frac = 0, nan_frac = 0, inf_frac = 0)
  for (vr in 0:2) {
    for (vi in 0:2) {
      expected <- cmul_ref(x, y, vr, vi)
      with_variants(c(vr, vi, 1L, 1L, 1L), {
        expect_tiers_give(expected, simd_mul, x, y)
        expect_tiers_give(expected[1:7], simd_mul, x[1:7], y[1:7])
      })
    }
  }
  # Division over the scaling ranges (no NaN results, so no fix-up).
  xw <- rand_z(3000, seed = 9L, wide = TRUE, na_frac = 0, nan_frac = 0, inf_frac = 0)
  yw <- rand_z(3000, seed = 10L, wide = TRUE, na_frac = 0, nan_frac = 0, inf_frac = 0)
  xw <- c(xw, 1 + 0i, 0 + 1i, 1e-310 + 2i, 1e308 + 1e308i)
  yw <- c(yw, 3 + 0i, 0 - 2i, 5e-320 + 1i, 1e308 - 1e300i)
  # Zero divisors, and NaN parts from the formula, are recovered by base
  # R's operator, not the formula. Base R's quotient alone does not tell:
  # it is not libgcc's when clang built R (macOS).
  no_nan <- function(z) !is.nan(Re(z)) & !is.nan(Im(z))
  keep <- no_nan(xw / yw) & no_nan(cdiv_ref(xw, yw, TRUE)) & no_nan(cdiv_ref(xw, yw, FALSE)) &
    yw != 0
  xw <- xw[keep]
  yw <- yw[keep]
  fused <- cdiv_ref(xw, yw, TRUE)
  if (simd_complex_variants()[["div"]] == "libgcc, fused") {
    expect_true(same_values(fused, xw / yw))
  }
  with_variants(c(1L, 1L, 1L, 1L, 1L), expect_tiers_give(fused, simd_div, xw, yw))
  with_variants(c(0L, 0L, 2L, 0L, 0L), {
    expect_tiers_give(cdiv_ref(xw, yw, FALSE), simd_div, xw, yw)
  })
  # Base R's operators alone.
  with_variants(c(-1L, -1L, 0L, -1L, -1L), {
    p <- grid_pairs(annex_g)
    expect_tiers_give(p$x * p$y, simd_mul, p$x, p$y)
    expect_tiers_give(p$x / p$y, simd_div, p$x, p$y)
  })
})

test_that("cumsum and cumprod are identical to base R, NA/NaN fix-up included", {
  for (n in c(0:9, 17, 100, 1000)) {
    x <- rand_z(n, seed = n + 20L)
    expect_tiers_give(cumsum(x), simd_cumsum, x)
    expect_tiers_give(cumprod(x), simd_cumprod, x)
  }
  z <- complex(modulus = 1 + (1:3000) * 1e-7, argument = 1:3000)
  expect_tiers_give(cumprod(z), simd_cumprod, z)
  expect_tiers_give(cumsum(z), simd_cumsum, z)
  na_1 <- complex(real = NA, imaginary = 1)
  nan_1 <- complex(real = NaN, imaginary = 1)
  for (v in list(
    c(1i, nan_1, 2, na_1, 3), c(1i, na_1, nan_1, 3), c(Inf + 0i, -Inf + 0i, 1i),
    c(1 + 1i, complex(real = 1, imaginary = NaN), 2), c(Inf * 1i, 0, nan_1)
  )) {
    expect_tiers_give(cumsum(v), simd_cumsum, v)
    expect_tiers_give(cumprod(v), simd_cumprod, v)
  }
  for (codes in list(c(1L, 1L, 1L, 0L, 0L), c(1L, 1L, 1L, 1L, 2L), c(1L, 1L, 1L, -1L, -1L))) {
    with_variants(codes, {
      out <- simd_cumprod(z[1:200])
      p <- 1 + 0i
      ref <- complex(200)
      for (i in 1:200) {
        p <- if (codes[4] < 0) z[i] * p else cmul_ref(z[i], p, codes[4], codes[5])
        ref[i] <- p
      }
      if (codes[4] >= 0) expect_true(same_values(out, ref), info = toString(codes))
    })
  }
})

test_that("prod follows the precision mode", {
  old <- simd_precision()
  on.exit(simd_precision(old))
  z <- complex(modulus = stats::runif(20000, 0.999, 1.001), argument = stats::runif(20000, -3, 3))
  base <- prod(z)
  # Base R multiplies in long double where it has one; without, its own
  # product is only within about 1e-11.
  tight <- if (has_wide_long_double()) 2^-50 else 1e-11
  for (mode in c("fast", "pairwise", "compensated")) {
    simd_precision(mode)
    res <- with_each_tier(function() simd_prod(z))
    for (tier in names(res)) {
      expect_lt(Mod(res[[tier]] - base) / Mod(base), if (mode == "compensated") tight else 1e-11)
    }
    # pairwise and compensated multiply in an order that does not depend on
    # the tier.
    if (mode != "fast") {
      expect_true(all(vapply(res, same_values, TRUE, res[["none"]])), info = mode)
    }
    expect_identical(simd_prod(complex(0)), 1 + 0i)
    expect_identical(simd_prod(c(1 + 1i, NA, NaN)), NA_complex_)
    expect_true(is.nan(Re(simd_prod(c(1 + 1i, NaN)))) && is.nan(Im(simd_prod(c(1 + 1i, NaN)))))
    expect_identical(simd_prod(c(2i, NA, 3), na.rm = TRUE), prod(c(2i, NA, 3), na.rm = TRUE))
    # Every multiply is base R's x * y, Annex G recovery included (base
    # prod multiplies in long double without it).
    one <- 1 + 0i
    simd_with_impl("none", {
      expect_identical(simd_prod(c(Inf + 0i, 1i)), one * (Inf + 0i) * 1i)
      if (mode != "compensated") expect_identical(simd_prod(z[1:3]), one * z[1] * z[2] * z[3])
    })
  }
})

test_that("mean is the sum divided by n, with sum's missing-value rules", {
  z <- rand_z(5001, seed = 30L, na_frac = 0, nan_frac = 0, inf_frac = 0)
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      s <- simd_sum(z)
      expect_identical(simd_mean(z), complex(real = Re(s) / 5001, imaginary = Im(s) / 5001))
    })
  }
  expect_lt(Mod(simd_mean(z) - mean(z)), 1e-12)
  expect_identical(
    simd_mean(c(1 + 3i, complex(real = 2, imaginary = NA))),
    mean(c(1 + 3i, complex(real = 2, imaginary = NA)))
  )
  expect_identical(simd_mean(c(1 + 3i, NA), na.rm = TRUE), 1 + 3i)
  expect_identical(simd_mean(complex(0)), mean(complex(0)))
})

test_that("eq and ne match base R, missing values included", {
  p <- grid_pairs(annex_g)
  expect_tiers_give(p$x == p$y, simd_eq, p$x, p$y)
  expect_tiers_give(p$x != p$y, simd_ne, p$x, p$y)
  for (n in lengths_c) {
    x <- rand_z(n, seed = n + 40L)
    y <- x
    y[seq_len(n) %% 3 == 0] <- 1i
    expect_tiers_give(x == y, simd_eq, x, y)
    expect_tiers_give(x != y, simd_ne, x, y)
  }
  x <- rand_z(100, seed = 41L)
  expect_identical(simd_eq(x, 1.5), x == 1.5)
  expect_identical(simd_ne(2L, x), 2L != x)
  expect_identical(simd_eq(c(1i, 0i), FALSE), c(1i, 0i) == FALSE)
  expect_identical(simd_eq(1i, 1i), TRUE)
  expect_warning(simd_eq(1i, structure(0, class = "integer64")), "integer64 coerced to double")
  for (f in list(simd_lt, simd_le, simd_gt, simd_ge)) {
    expect_error(f(1i, 1), "invalid comparison with complex values", fixed = TRUE)
    expect_error(f(2, 1i), "invalid comparison with complex values", fixed = TRUE)
  }
})

test_that("abs (Mod) and Arg match base R", {
  old <- simd_math_accuracy()
  on.exit(simd_math_accuracy(old))
  p <- grid_pairs(c(annex_g, 1e300, 5e-324))
  z <- c(p$x, rand_z(5000, seed = 50L, wide = TRUE))
  for (accuracy in c("accurate", "fast")) {
    simd_math_accuracy(accuracy)
    simd_with_impl("none", {
      expect_identical(simd_abs(z), Mod(z))
      expect_identical(simd_arg(z), Arg(z))
    })
    # Elements with a non-finite part are libm's on every tier.
    special <- !is.finite(Re(z)) | !is.finite(Im(z))
    for (tier in tiers_to_test()) {
      simd_with_impl(tier, {
        m <- simd_abs(z)
        a <- simd_arg(z)
        expect_identical(m[special], Mod(z)[special], info = tier)
        expect_identical(a[special], Arg(z)[special], info = tier)
        expect_null(math_mismatch(m, Mod(z), if (accuracy == "fast") 4 else 1), label = tier)
        expect_null(math_mismatch(a, Arg(z), if (accuracy == "fast") 4 else 2), label = tier)
      })
    }
  }
  expect_identical(simd_abs(c(-1.5, 2)), c(1.5, 2))
  expect_identical(simd_arg(c(-1, 0, 2, NA)), Arg(c(-1, 0, 2, NA)))
  expect_identical(simd_arg(c(-1L, 1L)), Arg(c(-1L, 1L)))
  expect_identical(simd_arg(c(TRUE, NA)), Arg(c(TRUE, NA)))
  expect_identical(simd_arg(complex(0)), double(0))
  expect_error(simd_arg("a"), "non-numeric argument to function", fixed = TRUE)
  expect_error(simd_arg(as.raw(1)), "non-numeric argument to function", fixed = TRUE)
  expect_error(simd_arg(structure(0, class = "integer64")),
    "simd_arg() does not support 'z' of type integer64",
    fixed = TRUE
  )
})

test_that("simd_vec * and / use the kernels", {
  z <- c(1 + 2i, -3i, NA)
  w <- c(2 - 1i, 1i, 1)
  expect_identical(unclass(simd_vec(z) * w), z * w)
  expect_identical(unclass(simd_vec(z) / w), z / w)
  kernel <- function(x, y, ...) "kernel"
  local_mocked_bindings(simd_mul = kernel, simd_div = kernel)
  expect_identical(simd_vec(z) * w, "kernel")
  expect_identical(w / simd_vec(z), "kernel")
})

test_that("chunk boundaries give the unchunked results", {
  skip_on_cran()
  skip_if_no_subprocess()
  n <- 2^20 + 7
  z <- complex(real = as.double(seq_len(n)) / 7, imaginary = -as.double(seq_len(n)) / 3)
  w <- complex(modulus = 1 + 1e-9 * seq_len(n), argument = seq_len(n))
  z[c(5, 300000)] <- complex(real = NA, imaginary = 1)
  w[1000001] <- NaN
  child <- function(z, w) {
    library(rsimd)
    lapply(stats::setNames(nm = simd_available()), function(t) {
      simd_with_impl(t, {
        old <- simd_precision("pairwise")
        on.exit(simd_precision(old))
        list(
          mul = simd_mul(z, w), div = simd_div(z, w), cumsum = simd_cumsum(z),
          cumprod = simd_cumprod(w), prod = simd_prod(w), eq = simd_eq(z, w),
          mod = simd_abs(w), mean = simd_mean(z, na.rm = TRUE)
        )
      })
    })
  }
  here <- child(z, w)
  expect_identical(here[["none"]]$mul, z * w)
  expect_identical(here[["none"]]$div, z / w)
  expect_identical(here[["none"]]$cumprod, cumprod(w))
  env <- c(callr::rcmd_safe_env(), RSIMD_DEBUG_STRIDE = "128")
  there <- callr::r(child, list(z, w), env = env)
  for (tier in intersect(names(here), names(there))) {
    for (nm in names(here[[tier]])) {
      expect_true(same_values(there[[tier]][[nm]], here[[tier]][[nm]]), info = paste(tier, nm))
    }
  }
})
