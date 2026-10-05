# Complex elementary functions (sqrt, exp, log, the trigonometric and
# hyperbolic functions and their inverses, powers, log with a base, atan2):
# the none tier is base R exactly, the other tiers within the stated
# bounds, with base R's special values, branch cuts and whole-number powers.

cm_fns <- c(
  "sqrt", "exp", "log", "sin", "cos", "tan", "sinh", "cosh", "tanh", "asin", "acos", "atan",
  "asinh", "acosh", "atanh"
)

# rsimd's bounds in ULPs (accurate, fast) as ?simd_exp states them, plus 6
# for the C library's own error: the tests compare with the none tier.
cm_bound <- list(
  sqrt = c(2, 3), exp = c(2.5, 3), log = c(3, 5), sin = c(3, 4), cos = c(3, 4),
  tan = c(4.5, 6), sinh = c(3, 4), cosh = c(3, 4), tanh = c(5, 6), asin = c(3.5, 4.5),
  acos = c(3.5, 4.5), atan = c(3, 3.5), asinh = c(3.5, 4.5), acosh = c(3.5, 4.5),
  atanh = c(3, 3.5), pow = c(3, 5), logb = c(4.5, 6.5), atan2 = c(3.5, 4)
)
cm_tol <- function(f, fast) cm_bound[[f]][[fast + 1L]] + 6

# Is the C library's complex arithmetic accurate (glibc, macOS)? Windows'
# is not: its clog, casin, catan, ctan and the rest lose many bits near the
# unit circle and the axes, and its cexp, csin and ccos reduce large
# arguments as inaccurately as its sin and cos. Base R's functions, and so
# the none tier, are then no reference for accuracy, and the SIMD tiers are
# compared with the first SIMD tier instead, whose accuracy the other
# platforms check.
libm_complex_accurate <- function() {
  z <- complex(real = 1 + 2^-30, imaginary = 2^-30)
  libm_reduces_trig() && abs(Re(log(z)) / (0.5 * log1p(2^-29 + 2^-59)) - 1) < 1e-12
}

# The tier, among `tiers`, that the other tiers' accuracy is checked against.
cm_ref_tier <- function(tiers) {
  simd <- setdiff(tiers, "none")
  if (libm_complex_accurate() || !length(simd)) "none" else simd[[1L]]
}

# How much log(z, w) and atan2(z, w) magnify the errors of the logarithms
# and arc tangent they are built from (1 where well conditioned): the
# relative error of log(v) is about 1 / |log(v)| ULP, and atan(u) has
# condition number |u| / (|1 + u^2| |atan(u)|).
cm_logb_scale <- function(z, w) {
  s <- 1 + 1 / Mod(log(w)) + 1 / Mod(log(z))
  s[!is.finite(s)] <- 1
  s
}
cm_atan2_scale <- function(z, w) {
  u <- z / w
  s <- 1 + suppressWarnings(Mod(u) / (Mod(1 + u^2) * Mod(atan(u))))
  s[!is.finite(s)] <- 1
  s
}

# Complex test inputs: moderate parts, the whole exponent range, near the
# unit circle, near +-1 and +-i, on and near the axes with signed zeros,
# and every combination of special parts.
cm_inputs <- function(n = 2000L, seed = 30L) {
  with_seed(seed, {
    r <- function(k) stats::rnorm(k, sd = 2)
    tiny <- function(k) stats::rnorm(k) * 10^-stats::runif(k, 0, 18)
    t <- stats::runif(n, -pi, pi)
    s <- sample(c(-1, 1), n, TRUE)
    c(
      complex(real = r(n), imaginary = r(n)),
      complex(
        real = r(n) * 10^stats::runif(n, -300, 300),
        imaginary = r(n) * 10^stats::runif(n, -300, 300)
      ),
      complex(modulus = 1 + tiny(n), argument = t),
      complex(real = s + tiny(n), imaginary = tiny(n)),
      complex(real = tiny(n), imaginary = s + tiny(n)),
      complex(real = r(n), imaginary = sample(c(0, -0), n, TRUE)),
      complex(real = sample(c(0, -0), n, TRUE), imaginary = r(n)),
      complex(real = r(n) * 400, imaginary = r(n) * 400),
      cm_specials()
    )
  })
}

cm_special_values <- c(
  0, -0, 0.5, -0.5, 1, -1, 2, -2, 30, -30, 710, -710, 1e300, -1e300, 1e-300, 5e-324, Inf, -Inf,
  NaN, NA
)
cm_specials <- function(v = cm_special_values) {
  g <- expand.grid(r = v, i = v)
  complex(real = g$r, imaginary = g$i)
}

bit_identical <- function(a, b) identical(a, b, num.eq = FALSE)

# TRUE where the sign bit of x is set (-0 included).
signbit <- function(x) ifelse(x == 0, 1 / x < 0, x < 0)

# The ULP distance of each part of a from b, or of |b| for a part much
# smaller than |b|, whichever is smaller.
cm_ulps <- function(a, b) {
  m <- Mod(b)
  ulp_m <- 2^(floor(log2(pmax(m, .Machine$double.xmin))) - 52)
  part <- function(p) pmin(ulp_dist(p(a), p(b)), abs(p(a) - p(b)) / ulp_m)
  pmax(part(Re), part(Im))
}

# NULL if a and b have the same missing kinds, infinities and (with zeros)
# signs of zero parts in each part, and are within ulps (a vector, per
# element) elsewhere; else a description.
cm_mismatch <- function(a, b, ulps, z = NULL, zeros = TRUE) {
  kind <- function(v) ifelse(is.na(v) & !is.nan(v), "NA", ifelse(is.nan(v), "NaN", ""))
  bad <- rep(FALSE, length(b))
  for (p in list(Re, Im)) {
    x <- p(a)
    y <- p(b)
    bad <- bad | kind(x) != kind(y) | xor(is.infinite(x), is.infinite(y)) |
      (is.infinite(x) & is.infinite(y) & x != y)
    if (zeros) bad <- bad | (!is.na(x) & !is.na(y) & x == 0 & y == 0 & zero_sign(x) != zero_sign(y))
  }
  num <- is.finite(Re(a)) & is.finite(Im(a)) & is.finite(Re(b)) & is.finite(Im(b))
  e <- cm_ulps(a, b)
  bad <- bad | (num & e > ulps)
  bad[is.na(bad)] <- TRUE
  if (!any(bad)) {
    return(NULL)
  }
  i <- which(bad)[1L]
  sprintf(
    "%d of %d differ; first at %s: %s vs %s (%g ULP, bound %g)", sum(bad), length(bad),
    if (is.null(z)) i else format(z[i], digits = 17), format(a[i], digits = 17),
    format(b[i], digits = 17), e[i], rep_len(ulps, length(b))[i]
  )
}

# The warnings an expression gives, and its value.
warnings_of <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

simd_fn <- function(f) get(paste0("simd_", f), envir = asNamespace("rsimd"))

with_accuracy <- function(acc, code) {
  old <- simd_math_accuracy(acc)
  on.exit(simd_math_accuracy(old))
  code
}

test_that("the none tier is identical to base R, warnings included", {
  z <- cm_inputs()
  for (f in setdiff(cm_fns, "acosh")) {
    b <- warnings_of(get(f, baseenv())(z))
    for (acc in c("accurate", "fast")) {
      r <- with_accuracy(acc, simd_with_impl("none", warnings_of(simd_fn(f)(z))))
      expect_true(bit_identical(r$value, b$value), label = paste(f, acc))
      expect_identical(r$warnings, b$warnings, label = paste(f, acc, "warnings"))
    }
  }
  expect_identical(
    warnings_of(simd_exp(complex(real = Inf, imaginary = Inf)))$warnings,
    "NaNs produced in function \"exp\""
  )
  expect_identical(warnings_of(simd_sqrt(1i))$warnings, character(0))
})

test_that("acosh is C99's principal value, not base R's acos(z) * i", {
  z <- cm_inputs(500L)
  ok <- !is.na(z)
  ref <- cm_ref_tier(tiers_to_test())
  for (tier in tiers_to_test()) {
    a <- simd_with_impl(tier, simd_acosh(z))
    # Re >= 0, and the imaginary part has the sign of Im(z).
    expect_true(all(Re(a[ok]) >= 0 | is.nan(Re(a[ok]))), label = tier)
    s <- !is.nan(Im(a)) & ok
    expect_identical(signbit(Im(a[s])), signbit(Im(z[s])), label = tier)
    # Base R's value negated where its real part is negative.
    b <- suppressWarnings(acosh(z))
    fix <- !is.na(b) & (Re(b) < 0)
    b[fix] <- -b[fix]
    if (tier != "none" && ref != "none") b <- simd_with_impl(ref, simd_acosh(z))
    off_axis <- ok & Im(z) != 0 & is.finite(Re(z)) & is.finite(Im(z))
    expect_null(cm_mismatch(a[off_axis], b[off_axis], 12, z[off_axis]), label = tier)
  }
  m0 <- -0
  w <- simd_acosh(complex(real = c(2, 2, -2, -2, 0.5, 0.5), imaginary = c(0, m0, 0, m0, 0, m0)))
  expect_equal(Re(w), c(rep(acosh(2), 4), 0, 0))
  expect_identical(zero_sign(Im(w[1:2])), c(1, -1))
  expect_equal(Im(w[3:6]), c(pi, -pi, acos(0.5), -acos(0.5)))
  expect_equal(simd_acosh(2 - 1i), Conj(simd_acosh(2 + 1i)))
  expect_equal(Re(simd_acosh(2 - 1i)), Re(acosh(2 + 1i)))
})

test_that("SIMD tiers are within the bounds and give base R's special values", {
  z <- cm_inputs()
  for (acc in c("accurate", "fast")) {
    fast <- acc == "fast"
    for (f in cm_fns) {
      res <- with_accuracy(acc, with_each_tier(function() suppressWarnings(simd_fn(f)(z))))
      ref_tier <- cm_ref_tier(names(res))
      ref <- res[[ref_tier]]
      for (tier in setdiff(names(res), c("none", ref_tier))) {
        expect_null(cm_mismatch(res[[tier]], ref, cm_tol(f, fast), z), label = paste(f, acc, tier))
      }
    }
  }
})

test_that("special values and branch cuts are base R's exactly on every tier", {
  z <- c(cm_specials(), cm_specials(c(0, -0, 1, -1, 2, -2, 1.5, -1.5, 25, -25)))
  ref_tier <- cm_ref_tier(tiers_to_test())
  for (f in cm_fns) {
    want <- suppressWarnings(simd_with_impl("none", simd_fn(f)(z)))
    exact <- !is.finite(Re(z)) | !is.finite(Im(z)) | !is.finite(Re(want)) | !is.finite(Im(want))
    close_to <- suppressWarnings(simd_with_impl(ref_tier, simd_fn(f)(z)))
    for (tier in tiers_to_test()) {
      got <- suppressWarnings(simd_with_impl(tier, simd_fn(f)(z)))
      expect_true(bit_identical(got[exact], want[exact]), label = paste(f, tier))
      if (tier == "none" && ref_tier != "none") next
      expect_null(cm_mismatch(got, close_to, cm_tol(f, FALSE), z), label = paste(f, tier))
    }
  }
})

test_that("whole-number powers are identical to base R on every tier", {
  z <- cm_inputs(300L)
  k <- with_seed(31L, sample(c(-70:70, 65536, -65536, 65535), length(z), TRUE))
  set <- c(-65536, -3, -2, -1, 0, 1, 2, 3, 17, 65536)
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_true(bit_identical(simd_pow(z, k), z^k), label = tier)
      expect_true(bit_identical(simd_pow(z, as.complex(k)), z^as.complex(k)), label = tier)
      for (p in set) expect_true(bit_identical(simd_pow(z, p), z^p), label = paste(tier, p))
      expect_true(bit_identical(simd_pow(1 + 1i, k), (1 + 1i)^k), label = tier)
    })
  }
})

test_that("general powers, log with a base and atan2 are within the bounds", {
  z <- cm_inputs(400L)
  w <- with_seed(32L, sample(c(
    z, complex(real = stats::rnorm(500), imaginary = stats::rnorm(500)),
    stats::rnorm(200), 0.5, -0.5, 1 / 3
  ), length(z), TRUE))
  for (acc in c("accurate", "fast")) {
    fast <- acc == "fast"
    with_accuracy(acc, {
      none <- simd_with_impl("none", list(
        pow = simd_pow(z, w), logb = suppressWarnings(simd_log(z, w)),
        atan2 = suppressWarnings(simd_atan2(z, w))
      ))
      scale <- 1 + Mod(w * log(z))
      scale[!is.finite(scale)] <- 1
      if (libm_complex_accurate()) {
        expect_true(bit_identical(none$pow, z^w))
      } else {
        # Windows: the polar form of base R's power calls pow, sin and cos,
        # which R.dll and this package's DLL need not take from the same
        # library, so the none tier matches base R only within the bound.
        # Nor need they keep the NA payload through those calls: for an NA
        # power one gives NA and the other NaN.
        nan_for_na <- function(v) {
          na <- function(x) ifelse(is.na(x), NaN, x)
          complex(real = na(Re(v)), imaginary = na(Im(v)))
        }
        tol <- cm_tol("pow", fast) * scale
        m <- cm_mismatch(nan_for_na(none$pow), nan_for_na(z^w), tol, z, zeros = FALSE)
        expect_null(m, label = paste("none pow", acc))
      }
      expect_true(bit_identical(none$logb, suppressWarnings(log(z, w))))
      expect_true(bit_identical(none$atan2, suppressWarnings(atan2(z, w))))
      ref_tier <- cm_ref_tier(tiers_to_test())
      ref <- if (ref_tier == "none") {
        none
      } else {
        simd_with_impl(ref_tier, list(
          pow = simd_pow(z, w), logb = suppressWarnings(simd_log(z, w)),
          atan2 = suppressWarnings(simd_atan2(z, w))
        ))
      }
      for (tier in setdiff(tiers_to_test(), c("none", ref_tier))) {
        simd_with_impl(tier, {
          # General powers: the bound grows with |w log z|; zero parts may
          # differ in sign.
          got <- simd_pow(z, w)
          tol <- cm_tol("pow", fast) * scale
          expect_null(cm_mismatch(got, ref$pow, tol, z, zeros = FALSE),
            label = paste("pow", acc, tier)
          )
          got <- suppressWarnings(simd_log(z, w))
          expect_null(cm_mismatch(got, ref$logb, cm_tol("logb", fast) * cm_logb_scale(z, w), z),
            label = paste("logb", acc, tier)
          )
          got <- suppressWarnings(simd_atan2(z, w))
          expect_null(cm_mismatch(got, ref$atan2, cm_tol("atan2", fast) * cm_atan2_scale(z, w), z),
            label = paste("atan2", acc, tier)
          )
        })
      }
    })
  }
  # Special cases of base R's power: 0^w, z^0, z^1.
  zero <- complex(real = c(0, -0, 0), imaginary = c(0, 0, -0))
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_true(bit_identical(simd_pow(zero, 2.5), zero^2.5), label = tier)
      expect_true(bit_identical(simd_pow(zero, -1), zero^-1), label = tier)
      expect_true(bit_identical(simd_pow(zero, 1i), zero^1i), label = tier)
      expect_true(bit_identical(simd_pow(z, 0), z^0), label = tier)
      expect_true(bit_identical(simd_pow(z, 1), z^1), label = tier)
    })
  }
})

test_that("log2, log10, log(z, base) and atan2 follow base R's rules", {
  z <- cm_inputs(200L)
  simd_with_impl("none", {
    expect_true(bit_identical(suppressWarnings(simd_log2(z)), suppressWarnings(log2(z))))
    expect_true(bit_identical(suppressWarnings(simd_log10(z)), suppressWarnings(log10(z))))
    expect_true(bit_identical(suppressWarnings(simd_log(z, 2)), suppressWarnings(log(z, 2))))
    expect_true(bit_identical(suppressWarnings(simd_log(2, z)), suppressWarnings(log(2, z))))
    expect_true(bit_identical(suppressWarnings(simd_atan2(z, 2)), suppressWarnings(atan2(z, 2))))
    expect_true(bit_identical(suppressWarnings(simd_atan2(2, z)), suppressWarnings(atan2(2, z))))
  })
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_true(
        bit_identical(suppressWarnings(simd_atan2(0i, 0i)), suppressWarnings(atan2(0i, 0i))),
        label = tier
      )
      expect_identical(
        warnings_of(simd_atan2(0i, 0i))$warnings,
        warnings_of(atan2(0i, 0i))$warnings
      )
      expect_identical(
        warnings_of(simd_log2(complex(real = NaN, imaginary = 1)))$warnings,
        warnings_of(log2(complex(real = NaN, imaginary = 1)))$warnings
      )
      expect_true(bit_identical(simd_log(NA_complex_, NA_complex_), log(NA_complex_, NA_complex_)))
      expect_true(bit_identical(
        suppressWarnings(simd_log(c(NA, 1 + 1i), NA)),
        suppressWarnings(log(c(NA, 1 + 1i), NA))
      ), label = tier)
    })
  }
  expect_identical(
    warnings_of(simd_log(-1 + 0i, 0i))$warnings,
    warnings_of(log(-1 + 0i, 0i))$warnings
  )
})

test_that("operands are promoted to complex, and other types are rejected as base R does", {
  simd_with_impl("none", {
    expect_true(bit_identical(simd_pow(2, 0.5i), 2^0.5i))
    expect_true(bit_identical(simd_pow(1:3, 1i), (1:3)^1i))
    expect_true(bit_identical(simd_pow(1i, 1:3), (1i)^(1:3)))
    expect_true(bit_identical(simd_pow(TRUE, 1i), TRUE^1i))
    expect_true(bit_identical(simd_log(1i, 2L), log(1i, 2L)))
    expect_true(bit_identical(simd_atan2(1, 1i), atan2(1, 1i)))
  })
  expect_identical(simd_pow(complex(0), 1i), complex(0))
  expect_identical(simd_sqrt(complex(0)), complex(0))
  expect_error(simd_pow(1i, as.raw(1)), "non-numeric argument to binary operator", fixed = TRUE)
  expect_error(simd_log(1i, "a"), "non-numeric argument to mathematical function", fixed = TRUE)
  expect_error(simd_pow(1i, 1:2 + 0i)[1], NA)
  expect_error(simd_pow(c(1i, 2i, 3i), 1:2))
  # Functions base R does not define on complex numbers.
  for (f in c("expm1", "log1p", "cbrt", "exp2", "exp10", "sinpi", "cospi", "tanpi")) {
    expect_error(simd_fn(f)(1i), "does not support 'x' of type complex", fixed = TRUE)
  }
  expect_error(simd_hypot(1i, 1), "does not support 'x' of type complex", fixed = TRUE)
  skip_if_not_installed("bit64")
  expect_error(simd_pow(1i, bit64::as.integer64(2)), "cannot combine complex and integer64",
    fixed = TRUE
  )
})

test_that("tails and lengths below the vector width match whole vectors", {
  z <- cm_inputs(40L)
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      for (f in c("sqrt", "exp", "log", "tan", "asin", "atanh")) {
        whole <- suppressWarnings(simd_fn(f)(z))
        for (n in c(1:9, 31:35)) {
          expect_true(bit_identical(suppressWarnings(simd_fn(f)(z[seq_len(n)])), whole[seq_len(n)]),
            label = paste(f, tier, n)
          )
        }
      }
    })
  }
})

test_that("simd_vec ^ and the Math group use the kernels", {
  z <- cm_inputs(50L)
  best <- simd_available()[[1L]]
  v <- simd_vec(z, impl = best)
  for (f in c("sqrt", "exp", "log", "log2", "log10", "sin", "acosh", "atanh")) {
    got <- suppressWarnings(get(f, baseenv())(v))
    expect_true(is_simd_vec(got), label = f)
    want <- suppressWarnings(simd_with_impl(best, simd_fn(f)(z)))
    expect_true(bit_identical(as.vector(unclass(got)), want),
      label = f
    )
  }
  expect_true(bit_identical(
    as.vector(unclass(suppressWarnings(log(v, 3i)))),
    suppressWarnings(simd_with_impl(best, simd_log(z, 3i)))
  ))
  kernel <- function(x, y, ...) "kernel"
  local_mocked_bindings(simd_pow = kernel)
  expect_identical(v^2, "kernel")
  expect_identical(2^simd_vec(1i), "kernel")
})

test_that("chunk boundaries give the unchunked results", {
  skip_on_cran()
  skip_if_no_subprocess()
  n <- 2^20 + 7
  z <- complex(real = sin(seq_len(n)) * 3, imaginary = cos(seq_len(n) / 7) * 2)
  z[c(5, 300000)] <- complex(real = NA, imaginary = 1)
  child <- function(z) {
    library(rsimd)
    lapply(stats::setNames(nm = simd_available()), function(t) {
      simd_with_impl(t, list(
        sqrt = simd_sqrt(z), exp = simd_exp(z), log = simd_log(z), asin = simd_asin(z),
        tanh = simd_tanh(z), pow = simd_pow(z, 1.5 - 0.5i), ipow = simd_pow(z, 3),
        atan2 = simd_atan2(z, rev(z))
      ))
    })
  }
  here <- child(z)
  env <- c(callr::rcmd_safe_env(), RSIMD_DEBUG_STRIDE = "128")
  there <- callr::r(child, list(z), env = env)
  for (tier in intersect(names(here), names(there))) {
    for (nm in names(here[[tier]])) {
      expect_true(bit_identical(there[[tier]][[nm]], here[[tier]][[nm]]), label = paste(tier, nm))
    }
  }
})

test_that("without a matching multiply or divide variant, powers stay exact", {
  # As on platforms where no variant reproduces base R (macOS division,
  # i686): base R's operators lane by lane.
  z <- cm_inputs(200L)
  k <- with_seed(33L, sample(-30:30, length(z), TRUE))
  w <- with_seed(34L, complex(real = stats::rnorm(length(z)), imaginary = stats::rnorm(length(z))))
  old <- .debug_c128_variants(c(-1L, -1L, 0L, -1L, -1L))
  on.exit(.debug_c128_variants(old))
  # The operators are then the package compiler's; when that is not the
  # compiler that built R (clang against a GCC-built R), they round
  # differently from base R's, which the load-time probe would have noticed.
  same_ops <- simd_with_impl(
    "none",
    bit_identical(simd_mul(z, w), z * w) && bit_identical(simd_div(1, z), 1 / z)
  )
  if (!same_ops) skip("this compiler's complex operators are not base R's")
  ref_tier <- cm_ref_tier(tiers_to_test())
  ref <- simd_with_impl(ref_tier, list(
    logb = suppressWarnings(simd_log(z, w)), atan2 = suppressWarnings(simd_atan2(z, w))
  ))
  logb_tol <- cm_tol("logb", FALSE) * cm_logb_scale(z, w)
  atan2_tol <- cm_tol("atan2", FALSE) * cm_atan2_scale(z, w)
  for (tier in tiers_to_test()) {
    simd_with_impl(tier, {
      expect_true(bit_identical(simd_pow(z, k), z^k), label = tier)
      expect_true(bit_identical(simd_pow(z, -3), z^-3), label = tier)
    })
    if (tier == "none" && ref_tier != "none") next
    simd_with_impl(tier, {
      expect_null(
        cm_mismatch(suppressWarnings(simd_log(z, w)), ref$logb, logb_tol, z),
        label = tier
      )
      expect_null(
        cm_mismatch(suppressWarnings(simd_atan2(z, w)), ref$atan2, atan2_tol, z),
        label = tier
      )
    })
  }
})
