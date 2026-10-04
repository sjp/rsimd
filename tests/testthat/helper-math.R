# Helpers for the elementary-function tests.

# The distance between doubles a and b in units in the last place (0 for
# two NAs, two NaNs, equal values and +0 vs -0; Inf for different missing
# kinds or a missing value against a number).
ulp_dist <- function(a, b) .ulp_dist(as.double(a), as.double(b))

# The 1/x sign of each zero of x: 1 for +0, -1 for -0, NA elsewhere.
zero_sign <- function(x) ifelse(!is.na(x) & x == 0, sign(1 / x), NA)

# NULL if a and b agree: the same missing kinds (NA, NaN), infinities and
# (unless zeros = FALSE) signs of zeros, and every other element within `ulps` units in the last
# place, or within `abs` where both are below `abs_below` in magnitude; else
# a description of the first difference.
math_mismatch <- function(a, b, ulps, x = NULL, zeros = TRUE, abs = 0, abs_below = 0) {
  if (!is.double(a) || length(a) != length(b)) return("types or lengths differ")
  kind <- function(v) ifelse(is.na(v) & !is.nan(v), "NA", ifelse(is.nan(v), "NaN", ""))
  bad <- kind(a) != kind(b) | xor(is.infinite(a), is.infinite(b))
  both_zero <- !is.na(a) & !is.na(b) & a == 0 & b == 0
  if (zeros) bad <- bad | (both_zero & zero_sign(a) != zero_sign(b))
  num <- !is.na(a) & !is.na(b) & !is.infinite(a) & !is.infinite(b)
  d <- ulp_dist(a, b)
  small <- num & pmax(abs(a), abs(b)) < abs_below & abs(a - b) <= abs
  bad <- bad | (num & d > ulps & !small) | (is.infinite(a) & a != b)
  bad[is.na(bad)] <- TRUE
  if (!any(bad)) return(NULL)
  i <- which(bad)[1L]
  sprintf(
    "%d of %d differ; first at %s: %a vs %a (%g ULP, bound %g)", sum(bad), length(bad),
    if (is.null(x)) i else sprintf("x = %a", x[i]), a[i], b[i], d[i], ulps
  )
}

# Expects f(...) on every tier to match the none tier within `ulps`
# (math_mismatch), and returns the none tier's value.
expect_tiers_close <- function(f, ..., ulps = 2, label = "", x = NULL) {
  res <- with_each_tier(function() suppressWarnings(f(...)))
  ref <- res[["none"]]
  problems <- character(0)
  for (tier in setdiff(names(res), "none")) {
    m <- math_mismatch(res[[tier]], ref, ulps, x)
    if (!is.null(m)) problems <- c(problems, sprintf("%s on %s: %s", label, tier, m))
  }
  testthat::expect(length(problems) == 0L, paste(problems, collapse = "\n"))
  invisible(ref)
}

# Expects `got` to match base R's `want` within `ulps` (math_mismatch).
expect_close_to <- function(got, want, ulps, label = "", ...) {
  m <- math_mismatch(got, want, ulps, ...)
  testthat::expect(is.null(m), sprintf("%s: %s", label, m))
}

# Does base R's sinpi() come from the platform's math library (glibc 2.41
# or later, macOS) rather than R's own fallback? The fallback returns +0
# for sinpi(-0) and is inaccurate near the zeros, so accuracy comparisons
# with base R's sinpi and cospi are only made when it is the platform's.
base_has_platform_sinpi <- function() 1 / sinpi(-0) < 0

# R 4.6.0 returns +0 instead of -0 for these functions of -0 (base R's
# small-argument shortcuts); rsimd keeps C's -0.
r46_plus_zero <- c(
  "sin", "tan", "asin", "atan", "sinh", "tanh", "asinh", "atanh", "expm1", "log1p"
)

# Special values for every function.
math_specials <- function() {
  c(
    0, -0, Inf, -Inf, NaN, NA, .Machine$double.xmin, 5e-324, -5e-324, .Machine$double.xmax,
    -.Machine$double.xmax, 1, -1, 0.5, -0.5, 2, -2, 1e-300, 1e300, -1e300, 1e-8, -1e-8
  )
}

# Seeded random values: n uniform on [lo, hi], or signs times 10^u for u
# uniform on [lo, hi] when log_scale.
math_random <- function(n, lo, hi, log_scale = FALSE, signed = TRUE, seed = 20261003L) {
  with_seed(seed, {
    if (log_scale) {
      v <- 10^stats::runif(n, lo, hi)
      if (signed) v <- v * sample(c(-1, 1), n, replace = TRUE)
      v
    } else {
      stats::runif(n, lo, hi)
    }
  })
}

# The test inputs of a function over its domain `dom`, list(lo, hi,
# log_scale): a grid and random values (about 2e3 in all by default, 1e5
# random values and a fine grid in extended runs) and subnormals.
math_inputs <- function(dom, extended = FALSE) {
  n_grid <- if (extended) 2e4 else 1000
  n_rand <- if (extended) 1e5 else 1000
  grid <- if (isTRUE(dom$log_scale)) {
    g <- 10^seq(dom$lo, dom$hi, length.out = n_grid)
    c(g, -g)
  } else {
    seq(dom$lo, dom$hi, length.out = n_grid)
  }
  denorm <- if (extended) {
    2^-(1022:1074) * rep(c(1, 1.5, -1, -1.75), each = 53)
  } else {
    c(5e-324, -5e-324, 1e-310, -1e-310, 2^-1023)
  }
  c(math_specials(), grid, math_random(n_rand, dom$lo, dom$hi, isTRUE(dom$log_scale)), denorm)
}

dom <- function(lo, hi, log_scale = FALSE) list(lo = lo, hi = hi, log_scale = log_scale)

# The unary functions: the rsimd function, its base R reference (NULL
# when base R has no accurate equivalent), the domain of the test inputs,
# and the ULP bounds against none and base R.
math1_table <- function() {
  # glibc's cbrt is up to 3 ULP from the exact value on aarch64, and base R
  # has no cbrt: compare with the none tier (libm) at 4 ULP and with
  # sign(x) * abs(x)^(1/3) loosely. glibc 2.43's sinpi is up to 2 ULP off.
  list(
    exp = list(simd_exp, exp, dom(-750, 710)),
    exp2 = list(simd_exp2, function(x) 2^x, dom(-1080, 1025)),
    exp10 = list(simd_exp10, function(x) 10^x, dom(-325, 309)),
    expm1 = list(simd_expm1, expm1, dom(-40, 710)),
    log = list(simd_log, log, dom(-310, 308, TRUE)),
    log2 = list(simd_log2, log2, dom(-310, 308, TRUE)),
    log10 = list(simd_log10, log10, dom(-310, 308, TRUE)),
    log1p = list(simd_log1p, log1p, dom(-1.5, 1e3)),
    cbrt = list(simd_cbrt, NULL, dom(-310, 308, TRUE), none = 4),
    sin = list(simd_sin, sin, dom(-1e3, 1e3)),
    cos = list(simd_cos, cos, dom(-1e3, 1e3)),
    tan = list(simd_tan, tan, dom(-1e3, 1e3)),
    asin = list(simd_asin, asin, dom(-1.05, 1.05)),
    acos = list(simd_acos, acos, dom(-1.05, 1.05)),
    atan = list(simd_atan, atan, dom(-300, 300, TRUE)),
    sinpi = list(simd_sinpi, sinpi, dom(-1e3, 1e3), base = 3, platform = TRUE),
    cospi = list(simd_cospi, cospi, dom(-1e3, 1e3), base = 3, platform = TRUE),
    tanpi = list(simd_tanpi, NULL, dom(-1e3, 1e3)),
    sinh = list(simd_sinh, sinh, dom(-712, 712)),
    cosh = list(simd_cosh, cosh, dom(-712, 712)),
    tanh = list(simd_tanh, tanh, dom(-25, 25)),
    asinh = list(simd_asinh, asinh, dom(-300, 308, TRUE)),
    acosh = list(simd_acosh, acosh, dom(-300, 308, TRUE)),
    atanh = list(simd_atanh, atanh, dom(-1.05, 1.05))
  )
}
