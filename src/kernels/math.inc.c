/* Elementary functions: out[i] = f(x[i]) (math1_f64), f(x[i], y[i])
 * (math2_f64), and sin and cos together (sincos_f64), with the op codes of
 * kernel_types.h.
 *
 * Missing values follow base R. A unary function returns a NaN input as it
 * is (payload kept, so NA stays NA), as R's math1 does; atan2 and hypot
 * give NA if either operand is NA and otherwise NaN if either is NaN, as
 * R's math2 does; pow follows R's ^ (rsimd_pow_f64). A NaN result from
 * operands that were not NaN sets RSIMD_EW_NAN_PRODUCED ("NaNs produced"),
 * except for pow, since base R's ^ never warns.
 *
 * The none tier calls libm, as base R does, except for sinpi, cospi and
 * tanpi, which reduce the argument exactly and then call libm's sin and
 * cos on |t| <= 1/4 or tan on |t| < 1/2 (base R uses the platform's sinpi when there is one and
 * a less accurate fallback otherwise), and exp2 and exp10, which are base
 * R's 2^x and 10^x. It is the reference the SIMD tiers are tested
 * against. The SIMD tiers call SLEEF (the RSIMD_SLEEF_* wrappers of
 * common.inc.h) with the same reductions and special cases, and recompute
 * with libm the rare lanes outside the range where SLEEF is accurate:
 * sinh and cosh for |x| > 709, asinh and acosh for |x| > 1e154, remainder
 * with |a| > 2^1000 |b| (SLEEF's quotient overflows), and pow with an
 * infinite operand (base R's rules there are not C's). In accurate
 * mode a tier may instead run the none tier's loops for the ops where libm was
 * measured faster than SLEEF on that tier (RSIMD_MATH1_LIBM and
 * RSIMD_MATH2_LIBM: log, log2, logb, cosh, asinh, acosh and pow on neon,
 * sve and sve2),
 * and in fast mode for those of them SLEEF has no fast variant of. A SIMD tier
 * built without SLEEF leaves these slots empty, so they are filled from
 * the next tier down, ultimately none.
 *
 * With RSIMD_MATH_FAST set in the op code (fast mode) the SIMD tiers call
 * the rsimd_sleef_<f>_fast wrappers instead, with the same reductions,
 * special cases and fallbacks, so only the last bits of finite results
 * change. The kernels test the bit once per call: each runs one of two
 * instantiations of its loops. The none tier ignores it.
 *
 * Signed zeros follow C and IEEE 754: f(-0) = -0 for every odd function
 * (base R 4.6 returns +0 for some of them). sinpi of an integer is a zero
 * with the sign of the argument; tanpi of an integer is +0, as in base R.
 */

#define RSIMD_MATH_PI 3.141592653589793
/* pi - RSIMD_MATH_PI, the low part of pi in double-double. */
#define RSIMD_MATH_PI_LO 1.2246467991473532e-16
/* log(2), log(10), 1 / log(2) and 1 / log(10) as double-doubles. */
#define RSIMD_MATH_LN2 0x1.62e42fefa39efp-1
#define RSIMD_MATH_LN2_LO 0x1.abc9e3b39803fp-56
#define RSIMD_MATH_LN10 0x1.26bb1bbb55516p+1
#define RSIMD_MATH_LN10_LO -0x1.f48ad494ea3e9p-53
#define RSIMD_MATH_INV_LN2 0x1.71547652b82fep+0
#define RSIMD_MATH_INV_LN2_LO 0x1.777d0ffda0d24p-56
#define RSIMD_MATH_INV_LN10 0x1.bcb7b1526e50ep-2
#define RSIMD_MATH_INV_LN10_LO 0x1.95355baaafad3p-57

/* ---- scalar functions ----------------------------------------------------
   Compiled in every tier: the none kernels use them, the SIMD tiers use
   rsimd_rpow_f64 for fallback lanes, and the vector layer test compares
   the SIMD kernels with them. */

/* Base R's x ^ y (R_pow in R's arithmetic.c), including its special cases
   y == 2, y == 3 or 4 with |x| <= 11, x == 0, and its rules for infinite
   operands, which differ from C's pow: (-Inf)^0.5 and (-2)^Inf are NaN. */
static inline double rsimd_rpow_f64(double x, double y) {
  if (y == 2.0) return x * x;
  if (x == 1.0 || y == 0.0) return 1.0;
  if (x == 0.0) {
    if (y > 0.0) return 0.0;
    if (y < 0.0) return INFINITY;
    return y;
  }
  if (x >= -11.0 && x <= 11.0) {
    if (y == 4.0) return x * x * x * x;
    if (y == 3.0) return x * x * x;
  }
  if (isfinite(x) && isfinite(y)) {
#if defined(_WIN64) && defined(__MINGW64_VERSION_MAJOR) && __MINGW64_VERSION_MAJOR >= 3
    /* R's USE_POWL_IN_R_POW: Mingw-w64's pow is inaccurate, so 64-bit
       Windows builds of R use powl. */
    return (double) powl(x, y);
#else
    return pow(x, y);
#endif
  }
  if (isnan(x) || isnan(y)) return x + y;
  if (isinf(x)) {
    if (x > 0) return y < 0.0 ? 0.0 : INFINITY;
    if (isfinite(y) && y == floor(y)) return y < 0.0 ? 0.0 : (fmod(y, 2.0) != 0 ? x : -x);
  }
  if (isinf(y) && x >= 0) {
    if (y > 0) return x >= 1 ? INFINITY : 0.0;
    return x < 1 ? INFINITY : 0.0;
  }
  return NAN;
}

/* x ^ y with the missing-value rule of the package's elementwise ops: 1
   when x is 1 or y is 0 (whatever the other operand), else NA if either
   operand is NA. */
static inline double rsimd_pow_f64(double x, double y) {
  if (x == 1.0 || y == 0.0) return 1.0;
  return rsimd_na_merge_f64(rsimd_rpow_f64(x, y), x, y);
}

/* atan2 and hypot with R's math2 missing-value rule. */
static inline double rsimd_math2_na_f64(double r, double x, double y) {
  if (rsimd_is_na_f64(x) || rsimd_is_na_f64(y)) return rsimd_na_real();
  if (isnan(x) || isnan(y)) return NAN;
  return r;
}

/* sin(pi t) and cos(pi t) for |t| <= 1/4, with pi t as the double-double
   hi + lo, so that libm's sin and cos and the final addition are the only
   rounding errors. */
static inline double rsimd_sinpi_small(double t) {
  double hi = RSIMD_MATH_PI * t, lo = rsimd_fma(RSIMD_MATH_PI, t, -hi) + RSIMD_MATH_PI_LO * t;
  return sin(hi) + cos(hi) * lo;
}
static inline double rsimd_cospi_small(double t) {
  double hi = RSIMD_MATH_PI * t, lo = rsimd_fma(RSIMD_MATH_PI, t, -hi) + RSIMD_MATH_PI_LO * t;
  return cos(hi) - sin(hi) * lo;
}

/* tan(pi t) for 0 < t < 1/2, corrected to first order for the low part
   of pi t (within about 1 ULP, measured). */
static inline double rsimd_tanpi_small(double t) {
  double hi = RSIMD_MATH_PI * t, lo = rsimd_fma(RSIMD_MATH_PI, t, -hi) + RSIMD_MATH_PI_LO * t;
  double th = tan(hi);
  return th + (1.0 + th * th) * lo;
}

/* sin(pi x): x reduced exactly modulo 2 into (-1, 1] as base R does, exact
   at the integers (a zero with the sign of x) and half-integers. */
static inline double rsimd_sinpi_f64(double x) {
  double r, a, s;
  if (isnan(x)) return x;
  if (isinf(x)) return NAN;
  r = fmod(x, 2.0);
  if (r <= -1.0) {
    r += 2.0;
  } else if (r > 1.0) {
    r -= 2.0;
  }
  if (r == 0.0 || r == 1.0) return copysign(0.0, x);
  a = fabs(r);
  if (a == 0.5) {
    s = 1.0;
  } else {
    if (a > 0.5) a = 1.0 - a;
    s = a <= 0.25 ? rsimd_sinpi_small(a) : rsimd_cospi_small(0.5 - a);
  }
  return r < 0 ? -s : s;
}

/* cos(pi x): |x| reduced exactly modulo 2, exact at the integers and
   half-integers (+0 there). */
static inline double rsimd_cospi_f64(double x) {
  double a, c;
  int neg;
  if (isnan(x)) return x;
  if (isinf(x)) return NAN;
  a = fmod(fabs(x), 2.0);
  if (a > 1.0) a = 2.0 - a;
  if (a == 0.5) return 0.0;
  if (a == 0.0) return 1.0;
  if (a == 1.0) return -1.0;
  neg = a > 0.5;
  if (neg) a = 1.0 - a;
  c = a <= 0.25 ? rsimd_cospi_small(a) : rsimd_sinpi_small(0.5 - a);
  return neg ? -c : c;
}

/* tan(pi x) as base R's tanpi: x reduced exactly modulo 1 into
   (-1/2, 1/2], +0 at the integers, NaN at the half-integers, exactly 1 and
   -1 at 1/4 and -1/4, otherwise tan(pi |r|) with the sign of r. */
static inline double rsimd_tanpi_f64(double x) {
  double r, a, t;
  if (isnan(x)) return x;
  if (isinf(x)) return NAN;
  r = fmod(x, 1.0);
  if (r <= -0.5) {
    r += 1.0;
  } else if (r > 0.5) {
    r -= 1.0;
  }
  if (r == 0.0) return 0.0;
  if (r == 0.5) return NAN;
  if (r == 0.25) return 1.0;
  if (r == -0.25) return -1.0;
  a = fabs(r);
  t = rsimd_tanpi_small(a);
  return r < 0 ? -t : t;
}

/* The logistic function, in the form that is accurate in both tails:
   1 / (1 + e) for x >= 0, e / (1 + e) below, with e = exp(-|x|). */
static inline double rsimd_sigmoid_f64(double x) {
  double e = exp(-fabs(x));
  return (x >= 0 ? 1.0 : e) / (1.0 + e);
}

/* C23's exp2m1(x) = 2^x - 1 (ten = 0) and exp10m1(x) = 10^x - 1 (ten =
   1), which C99 lacks. For t = x log(b), split exactly into th + tl (th =
   x c rounded, with c log(b) rounded, and tl its rounding error plus x
   times the low part of log(b)), expm1(t) is e + tl (1 + e) to first order,
   with e = expm1(th): accurate within expm1's error plus a rounding. Where
   b^x - 1 is exact (whole x in [-53, 53] for 2^x, in [1, 15] for 10^x) it
   is b^x - 1, from base R's b^x; above 54 (17) it is b^x itself, rounded,
   and below -1100 (-330) it is -1. Below 2^-1000 in magnitude th is
   correctly rounded, and keeps the sign of a zero. */
static inline double rsimd_expm1_base_f64(double x, int ten) {
  const double c = ten ? RSIMD_MATH_LN10 : RSIMD_MATH_LN2, cl = ten ? RSIMD_MATH_LN10_LO : RSIMD_MATH_LN2_LO;
  const double b = ten ? 10.0 : 2.0;
  double th, tl, e;
  if (isnan(x)) return x;
  if (x > (ten ? 17 : 54)) return rsimd_rpow_f64(b, x);
  if (x < (ten ? -330 : -1100)) return -1.0;
  if (x != 0 && x == trunc(x) && x <= (ten ? 15 : 53) && x >= (ten ? 1 : -53)) {
    return rsimd_rpow_f64(b, x) - 1.0;
  }
  th = x * c;
  if (fabs(x) < 0x1p-1000) return th;
  tl = rsimd_fma(x, c, -th) + x * cl;
  e = expm1(th);
  return rsimd_fma(tl, 1.0 + e, e);
}

/* C23's log2p1(x) = log2(1 + x) (ten = 0) and log10p1(x) = log10(1 + x)
   (ten = 1): log2(u) or log10(u) when u = 1 + x is exact (exact at the
   powers of 2 or 10), otherwise log1p(x) / log(b) with the division a
   multiplication by 1 / log(b) in double-double, within log1p's error plus
   a rounding. Zeros are returned as they are. */
static inline double rsimd_log1p_base_f64(double x, int ten) {
  double u = 1.0 + x, l;
  if (isnan(x) || x == 0) return x;
  if (u - 1.0 == x) return ten ? log10(u) : log2(u);
  l = log1p(x);
  return ten ? rsimd_fma(l, RSIMD_MATH_INV_LN10, l * RSIMD_MATH_INV_LN10_LO)
             : rsimd_fma(l, RSIMD_MATH_INV_LN2, l * RSIMD_MATH_INV_LN2_LO);
}

/* IEEE nextUp and nextDown (C23 nextup and nextdown, which C99 lacks). */
static inline double rsimd_next_up_f64(double x) { return nextafter(x, INFINITY); }
static inline double rsimd_next_down_f64(double x) { return nextafter(x, -INFINITY); }

/* x * 2^n for a whole n (C's scalbn: exact, or rounded once to a
   subnormal or infinite result). n beyond +-2200 changes nothing. */
static inline double rsimd_scaleb_f64(double x, double n) {
  return scalbn(x, n > 2200 ? 2200 : n < -2200 ? -2200 : (int) n);
}

/* r^an as the double-double hi + lo, for a whole an in [1, 512], by
   binary powering with exact products (fma gives each product's rounding
   error). The caller keeps r^an and its intermediate powers normal. */
static inline void rsimd_ipow_dd(double r, int an, double *hi, double *lo) {
  double ph = 1.0, pl = 0.0, bh = r, bl = 0.0, h, l;
  for (;;) {
    if (an & 1) {
      h = ph * bh;
      l = rsimd_fma(ph, bh, -h) + (ph * bl + pl * bh);
      ph = h + l;
      pl = l - (ph - h);
    }
    an >>= 1;
    if (an == 0) break;
    h = bh * bh;
    l = rsimd_fma(bh, bh, -h) + 2.0 * bh * bl;
    bh = h + l;
    bl = l - (bh - h);
  }
  *hi = ph;
  *lo = pl;
}

/* The real n-th root of x for a whole n, as C23's rootn. n = 0, and a
   negative x (-Inf too) with an even n, give NaN; zeros and infinities
   give a zero or an infinity, with the sign of x for an odd n; n = 1 is x,
   n = -1 is 1 / x and n = 2 is sqrt(x). Otherwise |x| = m * 2^(k |n|)
   exactly, with m in [1, 2^|n|) when |n| <= 512 (k = 0 above), and the
   root is r = exp2(log2(m) / n) after one Newton step, r - r c with
   c = ((r^n - m) / r^n) / n, scaled by 2^k (2^-k for a negative n). The
   step roughly squares the relative error of r, so r needs only about 30
   correct bits. For
   |n| <= 512, r^|n| is computed in double-double (rsimd_ipow_dd), which
   makes c almost exact and the root correctly rounded in all but rare
   ties (measured: within 0.5 ULP); for |n| > 512, where m can be
   subnormal or near overflow, c is ((q - 1) / q) / n with
   q = r^h (r^(n - h) / m) and h = trunc(n / 2), which neither overflows
   nor underflows. Exact for exact roots. The SIMD tiers run the same
   steps with SLEEF's functions (rsimd_math_rootn), its 3.5-ULP log2 for
   the first r. */
static inline double rsimd_rootn_f64(double x, double n) {
  double ax, an, k, m, r, c;
  int odd;
  if (isnan(x)) return x;
  if (n == 0) return NAN;
  if (n == 1) return x;
  if (n == -1) return 1.0 / x;
  odd = fmod(n, 2.0) != 0;
  if (x < 0 && !odd) return NAN;
  if (x == 0 || isinf(x)) {
    double mag = (x == 0) == (n > 0) ? 0.0 : INFINITY;
    return odd ? copysign(mag, x) : mag;
  }
  if (n == 2) return sqrt(x);
  ax = fabs(x);
  an = fabs(n);
  k = an <= 512 ? floor(ilogb(ax) / an) : 0;
  m = scalbn(ax, (int) (-k * an));
  r = exp2(log2(m) / n);
  if (an <= 512) {
    double h, l;
    rsimd_ipow_dd(r, (int) an, &h, &l);
    /* r^n = h + l for n > 0, 1 / (h + l) for n < 0. */
    c = (n > 0 ? ((h - m) + l) / h : rsimd_fma(-m, h, 1.0) - m * l) / n;
  } else {
    double h = trunc(n / 2), q = pow(r, h) * (pow(r, n - h) / m);
    c = ((q - 1) / q) / n;
  }
  r = r - r * c;
  r = scalbn(r, (int) (n > 0 ? k : -k));
  return odd ? copysign(r, x) : r;
}

/* f(x) for unary op code `op`; p is the divisor of LOGB. */
static inline double rsimd_math1_f64(int op, double x, double p) {
  switch (op) {
  case RSIMD_MATH_EXP: return exp(x);
  case RSIMD_MATH_EXP2: return rsimd_rpow_f64(2.0, x);
  case RSIMD_MATH_EXP10: return rsimd_rpow_f64(10.0, x);
  case RSIMD_MATH_EXPM1: return expm1(x);
  case RSIMD_MATH_LOG: return log(x);
  case RSIMD_MATH_LOG2: return log2(x);
  case RSIMD_MATH_LOG10: return log10(x);
  case RSIMD_MATH_LOG1P: return log1p(x);
  case RSIMD_MATH_LOGB: return log(x) / p;
  case RSIMD_MATH_CBRT: return cbrt(x);
  case RSIMD_MATH_SIN: return sin(x);
  case RSIMD_MATH_COS: return cos(x);
  case RSIMD_MATH_TAN: return tan(x);
  case RSIMD_MATH_ASIN: return asin(x);
  case RSIMD_MATH_ACOS: return acos(x);
  case RSIMD_MATH_ATAN: return atan(x);
  case RSIMD_MATH_SINPI: return rsimd_sinpi_f64(x);
  case RSIMD_MATH_COSPI: return rsimd_cospi_f64(x);
  case RSIMD_MATH_TANPI: return rsimd_tanpi_f64(x);
  case RSIMD_MATH_SINH: return sinh(x);
  case RSIMD_MATH_COSH: return cosh(x);
  case RSIMD_MATH_TANH: return tanh(x);
  case RSIMD_MATH_ASINH: return asinh(x);
  case RSIMD_MATH_ACOSH: return acosh(x);
  case RSIMD_MATH_ATANH: return atanh(x);
  case RSIMD_MATH_SIGMOID: return rsimd_sigmoid_f64(x);
  case RSIMD_MATH_EXP2M1: return rsimd_expm1_base_f64(x, 0);
  case RSIMD_MATH_EXP10M1: return rsimd_expm1_base_f64(x, 1);
  case RSIMD_MATH_LOG2P1: return rsimd_log1p_base_f64(x, 0);
  case RSIMD_MATH_LOG10P1: return rsimd_log1p_base_f64(x, 1);
  case RSIMD_MATH_NEXT_UP: return rsimd_next_up_f64(x);
  case RSIMD_MATH_NEXT_DOWN: return rsimd_next_down_f64(x);
  case RSIMD_MATH_RSQRT:
  case RSIMD_MATH_RSQRT_APPROX: return 1.0 / sqrt(x);
  case RSIMD_MATH_RECIP_APPROX: return 1.0 / x;
  default: return NAN;
  }
}

#if !RSIMD_TIER_IS(none) && !defined(RSIMD_NO_F64_SIMD)

/* ---- vector helpers that need no SLEEF ---------------------------------- */

/* 2^k for whole k in [-1022, 1023], built from the bits: k + 2^52 + 2^51
   holds k in its low bits. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_pow2i(rsimd_vf64 k) {
  const rsimd_vf64 magic = rsimd_vf64_set1(6755399441055744.0);
  rsimd_vi64 i = rsimd_vi64_sub(rsimd_vf64_as_vi64(rsimd_vf64_add(k, magic)),
                                rsimd_vf64_as_vi64(magic));
  return rsimd_vi64_as_vf64(rsimd_vi64_sll(rsimd_vi64_add(i, rsimd_vi64_set1(1023)), 52));
}

/* rsimd_scaleb_f64() lane by lane, for whole n (musl's scalbn): factors of
   2^1023 or 2^-969 (2^-1022 2^53, which keeps a second rounding away from
   the subnormal range) bring n into [-1022, 1023], then one multiply by
   2^n rounds once. Lanes where n is NaN are unspecified. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_scaleb(rsimd_vf64 x, rsimd_vf64 n) {
  const rsimd_vf64 one = rsimd_vf64_set1(1.0), zero = rsimd_vf64_zero(),
                   hi = rsimd_vf64_set1(1023.0), lo = rsimd_vf64_set1(-1022.0);
  int j;
  n = rsimd_vf64_min(rsimd_vf64_max(n, rsimd_vf64_set1(-2200.0)), rsimd_vf64_set1(2200.0));
  for (j = 0; j < 2; j++) {
    rsimd_mf64 m = rsimd_vf64_cmp_gt(n, hi);
    if (rsimd_mf64_any(m)) {
      x = rsimd_vf64_mul(x, rsimd_vf64_blend(one, rsimd_vf64_set1(0x1p1023), m));
      n = rsimd_vf64_sub(n, rsimd_vf64_blend(zero, hi, m));
    }
    m = rsimd_vf64_cmp_lt(n, lo);
    if (rsimd_mf64_any(m)) {
      x = rsimd_vf64_mul(x, rsimd_vf64_blend(one, rsimd_vf64_set1(0x1p-969), m));
      n = rsimd_vf64_add(n, rsimd_vf64_blend(zero, rsimd_vf64_set1(969.0), m));
    }
  }
  n = rsimd_vf64_min(rsimd_vf64_max(n, lo), hi);
  return rsimd_vf64_mul(x, rsimd_math_pow2i(n));
}

/* The unbiased exponent of each lane of ax (finite and > 0) as a double:
   C's ilogb, exact for subnormals, which are scaled by 2^54 first. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_ilogb(rsimd_vf64 ax) {
  const rsimd_vf64 two52 = rsimd_vf64_set1(4503599627370496.0);
  rsimd_mf64 sub = rsimd_vf64_cmp_lt(ax, rsimd_vf64_set1(0x1p-1022));
  rsimd_vf64 e;
  ax = rsimd_vf64_blend(ax, rsimd_vf64_mul(ax, rsimd_vf64_set1(0x1p54)), sub);
  /* The biased exponent E, as the double 2^52 + E minus 2^52. */
  e = rsimd_vf64_sub(rsimd_vi64_as_vf64(rsimd_vi64_or(rsimd_vi64_srl(rsimd_vf64_as_vi64(ax), 52),
                                                      rsimd_vf64_as_vi64(two52))),
                     two52);
  return rsimd_vf64_sub(e, rsimd_vf64_blend(rsimd_vf64_set1(1023.0), rsimd_vf64_set1(1077.0), sub));
}

/* rsimd_next_up_f64() lane by lane, on the bits: one more for a positive
   number, one less (a smaller magnitude) for a negative one, the smallest
   subnormal for a zero, and +Inf and NaN stay. A NaN is put back
   explicitly: stepping the bits of one with mantissa 1 gives an
   infinity, which RSIMD_MATH1_FINISH would not restore. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_next_up(rsimd_vf64 a) {
  const rsimd_vf64 zero = rsimd_vf64_zero();
  rsimd_vi64 d = rsimd_vi64_blend(rsimd_vi64_set1(-1), rsimd_vi64_set1(1),
                                  rsimd_mf64_to_mi64(rsimd_vf64_cmp_ge(a, zero)));
  rsimd_vf64 r = rsimd_vi64_as_vf64(rsimd_vi64_add(rsimd_vf64_as_vi64(a), d));
  r = rsimd_vf64_blend(r, rsimd_vf64_set1(4.9406564584124654e-324), rsimd_vf64_cmp_eq(a, zero));
  /* +Inf and NaN: the lanes where a < +Inf is false. */
  return rsimd_vf64_blend(r, a, rsimd_mf64_not(rsimd_vf64_cmp_lt(a, rsimd_vf64_set1(INFINITY))));
}

#endif

/* ---- the libm kernels ----------------------------------------------------
   math1_f64 and math2_f64 element by element with libm, as base R: the
   none tier's kernels, and what a SIMD tier runs for the ops it hands to
   libm (RSIMD_MATH1_LIBM, RSIMD_MATH2_LIBM below). */

/* Element i of operand k as a double (an int32 NA as NA_real_). */
static inline double rsimd_math_get(const void *p, int flags, int k, R_xlen_t i) {
  R_xlen_t j = (flags & RSIMD_EW_SCALAR(k)) ? 0 : i;
  if (flags & RSIMD_EW_I32(k)) {
    int v = ((const int *) p)[j];
    return v == RSIMD_NA_I32 ? rsimd_na_real() : (double) v;
  }
  return ((const double *) p)[j];
}

/* Runs out[i] = f(x[i]) over the chunk with math1's rule; f is an
   expression in a. Double input (flags 0) has its own loop. */
#define RSIMD_MATH1_NONE_ELT(get, f)                                             \
  for (i = 0; i < n; i++) {                                                      \
    double a = (get), r;                                                         \
    if (isnan(a)) {                                                              \
      r = a;                                                                     \
    } else {                                                                     \
      r = (f);                                                                   \
      if (isnan(r)) st = RSIMD_EW_NAN_PRODUCED;                                  \
    }                                                                            \
    out[i] = r;                                                                  \
  }
#define RSIMD_MATH1_NONE_LOOP(f)                                                 \
  if (flags == 0) {                                                              \
    RSIMD_MATH1_NONE_ELT(((const double *) x)[i], f)                             \
  } else {                                                                       \
    RSIMD_MATH1_NONE_ELT(rsimd_math_get(x, flags, 0, i), f)                      \
  }

/* math1_f64 for op (without RSIMD_MATH_FAST) with libm. */
static inline int rsimd_math1_libm(int op, const void *x, R_xlen_t n, int flags, double p,
                                   double *out) {
  int st = 0;
  R_xlen_t i;
  /* The common functions, and those a SIMD tier may hand to libm, get
     their own loops; the others go through rsimd_math1_f64(). */
  switch (op) {
  case RSIMD_MATH_EXP: RSIMD_MATH1_NONE_LOOP(exp(a)) break;
  case RSIMD_MATH_LOG: RSIMD_MATH1_NONE_LOOP(log(a)) break;
  case RSIMD_MATH_LOG2: RSIMD_MATH1_NONE_LOOP(log2(a)) break;
  case RSIMD_MATH_LOGB: RSIMD_MATH1_NONE_LOOP(log(a) / p) break;
  case RSIMD_MATH_SIN: RSIMD_MATH1_NONE_LOOP(sin(a)) break;
  case RSIMD_MATH_COS: RSIMD_MATH1_NONE_LOOP(cos(a)) break;
  case RSIMD_MATH_COSH: RSIMD_MATH1_NONE_LOOP(cosh(a)) break;
  case RSIMD_MATH_TANH: RSIMD_MATH1_NONE_LOOP(tanh(a)) break;
  case RSIMD_MATH_ASINH: RSIMD_MATH1_NONE_LOOP(asinh(a)) break;
  case RSIMD_MATH_ACOSH: RSIMD_MATH1_NONE_LOOP(acosh(a)) break;
  case RSIMD_MATH_SIGMOID: RSIMD_MATH1_NONE_LOOP(rsimd_sigmoid_f64(a)) break;
  case RSIMD_MATH_NEXT_UP: RSIMD_MATH1_NONE_LOOP(rsimd_next_up_f64(a)) break;
  case RSIMD_MATH_NEXT_DOWN: RSIMD_MATH1_NONE_LOOP(rsimd_next_down_f64(a)) break;
  case RSIMD_MATH_RSQRT:
  case RSIMD_MATH_RSQRT_APPROX: RSIMD_MATH1_NONE_LOOP(1.0 / sqrt(a)) break;
  case RSIMD_MATH_RECIP_APPROX: RSIMD_MATH1_NONE_LOOP(1.0 / a) break;
  default: RSIMD_MATH1_NONE_LOOP(rsimd_math1_f64(op, a, p)) break;
  }
  return st;
}

#undef RSIMD_MATH1_NONE_LOOP
#undef RSIMD_MATH1_NONE_ELT

/* math2_f64 for op (without RSIMD_MATH_FAST) with libm. */
static inline int rsimd_math2_libm(int op, const void *x, const void *y, R_xlen_t n, int flags,
                                   double *out) {
  int st = 0;
  R_xlen_t i;
  if (op == RSIMD_MATH_POW) {
    /* Base R's ^ never warns. */
    for (i = 0; i < n; i++) {
      out[i] = rsimd_pow_f64(rsimd_math_get(x, flags, 0, i), rsimd_math_get(y, flags, 1, i));
    }
    return 0;
  }
  for (i = 0; i < n; i++) {
    double a = rsimd_math_get(x, flags, 0, i), b = rsimd_math_get(y, flags, 1, i), r;
    switch (op) {
    case RSIMD_MATH_ATAN2: r = rsimd_math2_na_f64(atan2(a, b), a, b); break;
    case RSIMD_MATH_HYPOT: r = rsimd_math2_na_f64(hypot(a, b), a, b); break;
    case RSIMD_MATH_NEXTAFTER: r = rsimd_math2_na_f64(nextafter(a, b), a, b); break;
    case RSIMD_MATH_REMAINDER: r = rsimd_math2_na_f64(remainder(a, b), a, b); break;
    case RSIMD_MATH_SCALEB:
      r = rsimd_math2_na_f64(isnan(b) ? b : rsimd_scaleb_f64(a, b), a, b);
      break;
    case RSIMD_MATH_ROOTN:
      r = rsimd_math2_na_f64(isnan(b) ? b : rsimd_rootn_f64(a, b), a, b);
      break;
    default: r = NAN; break;
    }
    if (isnan(r) && !isnan(a) && !isnan(b)) st = RSIMD_EW_NAN_PRODUCED;
    out[i] = r;
  }
  return st;
}

#if RSIMD_TIER_IS(none)

int RSIMD_KERNEL(math1_f64)(int op, const void *x, R_xlen_t n, int flags, double p, double *out);
int RSIMD_KERNEL(math1_f64)(int op, const void *x, R_xlen_t n, int flags, double p, double *out) {
  return rsimd_math1_libm(op & ~RSIMD_MATH_FAST, x, n, flags, p, out);
}

int RSIMD_KERNEL(math2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                            double *out);
int RSIMD_KERNEL(math2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                            double *out) {
  return rsimd_math2_libm(op & ~RSIMD_MATH_FAST, x, y, n, flags, out);
}

int RSIMD_KERNEL(sincos_f64)(int op, const void *x, R_xlen_t n, int flags, double *s,
                             double *c);
int RSIMD_KERNEL(sincos_f64)(int op, const void *x, R_xlen_t n, int flags, double *s,
                             double *c) {
  int st = 0, pi = (op & ~RSIMD_MATH_FAST) == RSIMD_MATH_SINPI;
  R_xlen_t i;
  /* Separate loops: in one, GCC fuses sin(a) and cos(a) into sincos(),
     which glibc does not always round as sin() (as base R calls). */
  for (i = 0; i < n; i++) {
    double a = rsimd_math_get(x, flags, 0, i);
    s[i] = isnan(a) ? a : pi ? rsimd_sinpi_f64(a) : sin(a);
    if (!isnan(a) && isnan(s[i])) st = RSIMD_EW_NAN_PRODUCED;
  }
  for (i = 0; i < n; i++) {
    double a = rsimd_math_get(x, flags, 0, i);
    c[i] = isnan(a) ? a : pi ? rsimd_cospi_f64(a) : cos(a);
    if (!isnan(a) && isnan(c[i])) st = RSIMD_EW_NAN_PRODUCED;
  }
  return st;
}

#elif defined(RSIMD_HAVE_SLEEF)

/* Lanes of the vector r (of `lanes` active lanes) where `m` is set,
   recomputed from the operands a (and b) with the scalar function f1 (or
   f2): the fallback for lanes outside SLEEF's accurate range. */
static inline rsimd_vf64 rsimd_math_fix1(rsimd_vf64 r, rsimd_vf64 a, rsimd_mf64 m, int lanes,
                                         double (*f1)(double)) {
  double ba[RSIMD_MAX_LANES_64], br[RSIMD_MAX_LANES_64], bm[RSIMD_MAX_LANES_64];
  int j;
  rsimd_vf64_storeu(ba, a);
  rsimd_vf64_storeu(br, r);
  rsimd_vf64_storeu(bm, rsimd_vf64_blend(rsimd_vf64_zero(), rsimd_vf64_set1(1.0), m));
  for (j = 0; j < lanes; j++) {
    if (bm[j] != 0) br[j] = f1(ba[j]);
  }
  return rsimd_vf64_loadu(br);
}
static inline rsimd_vf64 rsimd_math_fix2(rsimd_vf64 r, rsimd_vf64 a, rsimd_vf64 b, rsimd_mf64 m,
                                         int lanes, double (*f2)(double, double)) {
  double ba[RSIMD_MAX_LANES_64], bb[RSIMD_MAX_LANES_64], br[RSIMD_MAX_LANES_64],
    bm[RSIMD_MAX_LANES_64];
  int j;
  rsimd_vf64_storeu(ba, a);
  rsimd_vf64_storeu(bb, b);
  rsimd_vf64_storeu(br, r);
  rsimd_vf64_storeu(bm, rsimd_vf64_blend(rsimd_vf64_zero(), rsimd_vf64_set1(1.0), m));
  for (j = 0; j < lanes; j++) {
    if (bm[j] != 0) br[j] = f2(ba[j], bb[j]);
  }
  return rsimd_vf64_loadu(br);
}

/* Lanes where |a| > t. */
RSIMD_ALWAYS_INLINE rsimd_mf64 rsimd_math_abs_gt(rsimd_vf64 a, double t) {
  return rsimd_vf64_cmp_gt(rsimd_vf64_abs(a), rsimd_vf64_set1(t));
}

/* f(a) for the libm function f, recomputed in the lanes where |a| > t. */
#define RSIMD_MATH_FIXED(r, a, t, f)                                             \
  do {                                                                           \
    rsimd_mf64 big_ = rsimd_math_abs_gt(a, t);                                   \
    if (rsimd_mf64_any(big_)) r = rsimd_math_fix1(r, a, big_, lanes, f);         \
  } while (0)

/* rsimd_sleef_<f>(args), or rsimd_sleef_<f>_fast(args) when `fast` (a
   constant in each instantiation) is set. */
#define RSIMD_SLEEF_CALL(f, ...)                                                 \
  (fast ? rsimd_sleef_##f##_fast(__VA_ARGS__) : rsimd_sleef_##f(__VA_ARGS__))

/* A zero with the sign of each lane of a. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_signed_zero(rsimd_vf64 a) {
  return rsimd_vf64_and(rsimd_vf64_set1(-0.0), a);
}

/* rsimd_sinpi_f64() lane by lane (none tier): x - 2 trunc(x / 2) is
   fmod(x, 2) up to the sign of a zero result, exactly. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_sinpi(rsimd_vf64 x, int fast) {
  const rsimd_vf64 one = rsimd_vf64_set1(1.0), two = rsimd_vf64_set1(2.0);
  rsimd_vf64 r = rsimd_vf64_sub(x, rsimd_vf64_mul(two, rsimd_vf64_trunc(rsimd_vf64_mul(x, rsimd_vf64_set1(0.5)))));
  rsimd_vf64 s;
  r = rsimd_vf64_blend(r, rsimd_vf64_add(r, two), rsimd_vf64_cmp_le(r, rsimd_vf64_set1(-1.0)));
  r = rsimd_vf64_blend(r, rsimd_vf64_sub(r, two), rsimd_vf64_cmp_gt(r, one));
  s = RSIMD_SLEEF_CALL(sinpi, r);
  return rsimd_vf64_blend(
    s, rsimd_math_signed_zero(x),
    rsimd_mf64_or(rsimd_vf64_cmp_eq(r, rsimd_vf64_zero()), rsimd_vf64_cmp_eq(r, one)));
}

/* rsimd_cospi_f64() lane by lane. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_cospi(rsimd_vf64 x, int fast) {
  const rsimd_vf64 two = rsimd_vf64_set1(2.0), half = rsimd_vf64_set1(0.5);
  rsimd_vf64 a = rsimd_vf64_abs(x);
  a = rsimd_vf64_sub(a, rsimd_vf64_mul(two, rsimd_vf64_trunc(rsimd_vf64_mul(a, half))));
  a = rsimd_vf64_blend(a, rsimd_vf64_sub(two, a), rsimd_vf64_cmp_gt(a, rsimd_vf64_set1(1.0)));
  /* SLEEF is exact at 0 and 1 but gives -0 at 1/2. */
  return rsimd_vf64_blend(RSIMD_SLEEF_CALL(cospi, a), rsimd_vf64_zero(),
                          rsimd_vf64_cmp_eq(a, half));
}

/* rsimd_tanpi_f64() lane by lane. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_tanpi(rsimd_vf64 x, int fast) {
  const rsimd_vf64 one = rsimd_vf64_set1(1.0), half = rsimd_vf64_set1(0.5),
                   quarter = rsimd_vf64_set1(0.25), zero = rsimd_vf64_zero();
  rsimd_vf64 r = rsimd_vf64_sub(x, rsimd_vf64_trunc(x)), s, c, t;
  r = rsimd_vf64_blend(r, rsimd_vf64_add(r, one), rsimd_vf64_cmp_le(r, rsimd_vf64_neg(half)));
  r = rsimd_vf64_blend(r, rsimd_vf64_sub(r, one), rsimd_vf64_cmp_gt(r, half));
  if (fast) {
    rsimd_sleef_sincospi_fast(r, &s, &c);
  } else {
    rsimd_sleef_sincospi(r, &s, &c);
  }
  t = rsimd_vf64_div(s, c);
  t = rsimd_vf64_blend(t, zero, rsimd_vf64_cmp_eq(r, zero));
  t = rsimd_vf64_blend(t, rsimd_vf64_set1(NAN), rsimd_vf64_cmp_eq(r, half));
  t = rsimd_vf64_blend(t, one, rsimd_vf64_cmp_eq(r, quarter));
  return rsimd_vf64_blend(t, rsimd_vf64_neg(one), rsimd_vf64_cmp_eq(r, rsimd_vf64_neg(quarter)));
}

/* rsimd_pow_f64() lane by lane: SLEEF's pow for finite operands, base R's
   special cases blended over it in reverse order of precedence, and the
   lanes with an infinite operand recomputed by the scalar rules. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_pow(rsimd_vf64 a, rsimd_vf64 b, int lanes) {
  const rsimd_vf64 zero = rsimd_vf64_zero(), one = rsimd_vf64_set1(1.0),
                   inf = rsimd_vf64_set1(INFINITY);
  rsimd_vf64 r = rsimd_sleef_pow(a, b), a3 = rsimd_vf64_mul(rsimd_vf64_mul(a, a), a);
  rsimd_mf64 m, az = rsimd_vf64_cmp_eq(a, zero);
  rsimd_mf64 small = rsimd_vf64_cmp_le(rsimd_vf64_abs(a), rsimd_vf64_set1(11.0));
  r = rsimd_vf64_blend(r, a3, rsimd_mf64_and(small, rsimd_vf64_cmp_eq(b, rsimd_vf64_set1(3.0))));
  r = rsimd_vf64_blend(r, rsimd_vf64_mul(a3, a),
                       rsimd_mf64_and(small, rsimd_vf64_cmp_eq(b, rsimd_vf64_set1(4.0))));
  r = rsimd_vf64_blend(r, zero, rsimd_mf64_and(az, rsimd_vf64_cmp_gt(b, zero)));
  r = rsimd_vf64_blend(r, inf, rsimd_mf64_and(az, rsimd_vf64_cmp_lt(b, zero)));
  m = rsimd_mf64_or(rsimd_vf64_cmp_eq(rsimd_vf64_abs(a), inf),
                    rsimd_vf64_cmp_eq(rsimd_vf64_abs(b), inf));
  if (rsimd_mf64_any(m)) r = rsimd_math_fix2(r, a, b, m, lanes, rsimd_rpow_f64);
  r = rsimd_vf64_blend(r, rsimd_vf64_mul(a, a), rsimd_vf64_cmp_eq(b, rsimd_vf64_set1(2.0)));
  m = rsimd_mf64_or(rsimd_vf64_is_nan(a), rsimd_vf64_is_nan(b));
  if (rsimd_mf64_any(m)) {
    r = rsimd_vf64_blend(r, rsimd_vf64_set1(NAN), m);
    r = rsimd_vf64_na_merge(r, a, b);
  }
  return rsimd_vf64_blend(r, one,
                          rsimd_mf64_or(rsimd_vf64_cmp_eq(a, one), rsimd_vf64_cmp_eq(b, zero)));
}

/* The vector form of rsimd_sigmoid_f64(). */
static inline rsimd_vf64 rsimd_math_sigmoid(rsimd_vf64 a) {
  const rsimd_vf64 one = rsimd_vf64_set1(1.0);
  rsimd_vf64 e = rsimd_sleef_exp(rsimd_vf64_neg(rsimd_vf64_abs(a)));
  rsimd_vf64 num = rsimd_vf64_blend(e, one, rsimd_vf64_cmp_ge(a, rsimd_vf64_zero()));
  return rsimd_vf64_div(num, rsimd_vf64_add(one, e));
}

/* Runs `expr`, which sets the rsimd_vf64 r from a (and b), over the
   chunk, with the loads of the elementwise kernels (arith.inc.c); `lanes`
   is the number of active lanes. */
#define RSIMD_MATH_LOOP(nargs, expr)                                             \
  do {                                                                           \
    ptrdiff_t i = 0;                                                             \
    for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {                       \
      const int lanes = (int) RSIMD_LANES_64;                                    \
      rsimd_vf64 a = rsimd_ew_ld(x, flags, 0, bc0, i, 1), b = bc1, r;            \
      if ((nargs) > 1) b = rsimd_ew_ld(y, flags, 1, bc1, i, 1);                  \
      (void) b;                                                                  \
      (void) lanes;                                                              \
      expr;                                                                      \
      rsimd_vf64_storeu(out + i, r);                                             \
    }                                                                            \
    if (i < n) {                                                                 \
      rsimd_p64 pg = rsimd_p64_while(i, n);                                      \
      const int lanes = rsimd_p64_count(pg);                                     \
      rsimd_vf64 a = rsimd_ew_ld_p(x, flags, 0, bc0, i, pg, 1), b = bc1, r;      \
      if ((nargs) > 1) b = rsimd_ew_ld_p(y, flags, 1, bc1, i, pg, 1);            \
      (void) b;                                                                  \
      (void) lanes;                                                              \
      expr;                                                                      \
      rsimd_vf64_storeu_p(pg, out + i, r);                                       \
    }                                                                            \
  } while (0)

/* Finishes r = f(a) with math1's rule: a NaN input is returned as it is,
   and a NaN result from any other input sets the status. Every function
   gives NaN for a NaN input, so only a vector with a NaN result needs
   either. */
#define RSIMD_MATH1_FINISH()                                                     \
  do {                                                                           \
    rsimd_mf64 rn_ = rsimd_vf64_is_nan(r);                                       \
    if (rsimd_mf64_any(rn_)) {                                                   \
      rsimd_mf64 in_ = rsimd_vf64_is_nan(a);                                     \
      if (rsimd_mf64_any(rsimd_mf64_andnot(in_, rn_))) st = RSIMD_EW_NAN_PRODUCED; \
      r = rsimd_vf64_blend(r, a, in_);                                           \
    }                                                                            \
  } while (0)

#define RSIMD_MATH1_RESULT(e)                                                    \
  do {                                                                           \
    r = (e);                                                                     \
    RSIMD_MATH1_FINISH();                                                        \
  } while (0)

/* As RSIMD_MATH1_RESULT, with the lanes where |a| > t recomputed by the
   libm function f. */
#define RSIMD_MATH1_FIXED(e, t, f)                                               \
  do {                                                                           \
    r = (e);                                                                     \
    RSIMD_MATH_FIXED(r, a, t, f);                                                \
    RSIMD_MATH1_FINISH();                                                        \
  } while (0)

/* For the functions with f(x) ~ x at 0 (sin, tan, asin, atan, sinh,
   tanh, asinh, atanh, expm1, log1p): x itself where |x| < 2^-1000, which
   is f(x) correctly rounded. SLEEF loses the smallest subnormals and the
   sign of their zero result on some tiers. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_tiny(rsimd_vf64 r, rsimd_vf64 a) {
  rsimd_mf64 tiny = rsimd_vf64_cmp_lt(rsimd_vf64_abs(a), rsimd_vf64_set1(0x1p-1000));
  return rsimd_mf64_any(tiny) ? rsimd_vf64_blend(r, a, tiny) : r;
}

/* math1_f64 with `fast` constant: rsimd_math1_run(..., 0) is the
   accurate kernel and rsimd_math1_run(..., 1) the fast one. */
RSIMD_ALWAYS_INLINE int rsimd_math1_run(int op, const void *x, R_xlen_t n, int flags, double p,
                                        double *out, int fast) {
  const void *y = NULL;
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, 1), bc1 = rsimd_vf64_zero();
  const rsimd_vf64 vp = rsimd_vf64_set1(p);
  int st = 0;
  (void) y;
#define RSIMD_MATH1_CASE(OP, f)                                                  \
  case OP: RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(RSIMD_SLEEF_CALL(f, a))); break;
#define RSIMD_MATH1_CASE_ODD(OP, f)                                              \
  case OP:                                                                       \
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_tiny(RSIMD_SLEEF_CALL(f, a), a))); \
    break;
  switch (op) {
    RSIMD_MATH1_CASE(RSIMD_MATH_EXP, exp)
    RSIMD_MATH1_CASE(RSIMD_MATH_EXP2, exp2)
    RSIMD_MATH1_CASE(RSIMD_MATH_EXP10, exp10)
    RSIMD_MATH1_CASE_ODD(RSIMD_MATH_EXPM1, expm1)
    RSIMD_MATH1_CASE(RSIMD_MATH_LOG, log)
    RSIMD_MATH1_CASE(RSIMD_MATH_LOG2, log2)
    RSIMD_MATH1_CASE(RSIMD_MATH_LOG10, log10)
    RSIMD_MATH1_CASE_ODD(RSIMD_MATH_LOG1P, log1p)
  case RSIMD_MATH_LOGB:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_vf64_div(RSIMD_SLEEF_CALL(log, a), vp)));
    break;
    RSIMD_MATH1_CASE(RSIMD_MATH_CBRT, cbrt)
    RSIMD_MATH1_CASE_ODD(RSIMD_MATH_SIN, sin)
    RSIMD_MATH1_CASE(RSIMD_MATH_COS, cos)
    RSIMD_MATH1_CASE_ODD(RSIMD_MATH_TAN, tan)
    RSIMD_MATH1_CASE_ODD(RSIMD_MATH_ASIN, asin)
    RSIMD_MATH1_CASE(RSIMD_MATH_ACOS, acos)
    RSIMD_MATH1_CASE_ODD(RSIMD_MATH_ATAN, atan)
  case RSIMD_MATH_SINPI:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_sinpi(a, fast)));
    break;
  case RSIMD_MATH_COSPI:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_cospi(a, fast)));
    break;
  case RSIMD_MATH_TANPI:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_tanpi(a, fast)));
    break;
  case RSIMD_MATH_SINH:
    RSIMD_MATH_LOOP(
      1, RSIMD_MATH1_FIXED(rsimd_math_tiny(RSIMD_SLEEF_CALL(sinh, a), a), 709.0, sinh));
    break;
  case RSIMD_MATH_COSH:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_FIXED(RSIMD_SLEEF_CALL(cosh, a), 709.0, cosh));
    break;
    RSIMD_MATH1_CASE_ODD(RSIMD_MATH_TANH, tanh)
  case RSIMD_MATH_ASINH:
    RSIMD_MATH_LOOP(
      1, RSIMD_MATH1_FIXED(rsimd_math_tiny(RSIMD_SLEEF_CALL(asinh, a), a), 1e154, asinh));
    break;
  case RSIMD_MATH_ACOSH:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_FIXED(RSIMD_SLEEF_CALL(acosh, a), 1e154, acosh));
    break;
    RSIMD_MATH1_CASE_ODD(RSIMD_MATH_ATANH, atanh)
  /* Only exp, which has no fast variant. */
  case RSIMD_MATH_SIGMOID: RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_sigmoid(a))); break;
  default: break;
  }
#undef RSIMD_MATH1_CASE
#undef RSIMD_MATH1_CASE_ODD
  return st;
}

/* 1 / a and 1 / sqrt(a) within 2^-22 (rsimd_vf64_recip_approx and
   rsqrt_approx), exact in the lanes outside their range: zeros,
   infinities, NaN, subnormals, |a| >= 2^1022 and, for rsqrt, negative
   numbers. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_recip_approx(rsimd_vf64 a) {
#ifdef RSIMD_RECIP_APPROX_EXACT
  /* The layer's recip_approx is the exact quotient in every lane. */
  return rsimd_vf64_recip_approx(a);
#else
  rsimd_vf64 ax = rsimd_vf64_abs(a), r = rsimd_vf64_recip_approx(a);
  rsimd_mf64 ok = rsimd_mf64_and(rsimd_vf64_cmp_ge(ax, rsimd_vf64_set1(0x1p-1022)),
                                 rsimd_vf64_cmp_lt(ax, rsimd_vf64_set1(0x1p1022)));
  if (!rsimd_mf64_all(ok)) r = rsimd_vf64_blend(rsimd_vf64_div(rsimd_vf64_set1(1.0), a), r, ok);
  return r;
#endif
}
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_rsqrt_approx(rsimd_vf64 a) {
  rsimd_vf64 r = rsimd_vf64_rsqrt_approx(a);
  rsimd_mf64 ok = rsimd_mf64_and(rsimd_vf64_cmp_ge(a, rsimd_vf64_set1(0x1p-1022)),
                                 rsimd_vf64_cmp_lt(a, rsimd_vf64_set1(INFINITY)));
  if (!rsimd_mf64_all(ok)) {
    r = rsimd_vf64_blend(rsimd_vf64_div(rsimd_vf64_set1(1.0), rsimd_vf64_sqrt(a)), r, ok);
  }
  return r;
}

/* rsimd_expm1_base_f64() lane by lane, with SLEEF's expm1, exp2 and
   exp10 (the 1-ULP ones, which are exact at whole numbers). */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_expm1_base(rsimd_vf64 a, const int ten) {
  const rsimd_vf64 one = rsimd_vf64_set1(1.0), c = rsimd_vf64_set1(ten ? RSIMD_MATH_LN10 : RSIMD_MATH_LN2);
  rsimd_vf64 th = rsimd_vf64_mul(a, c), tl, e, r;
  rsimd_mf64 m;
  tl = rsimd_vf64_add(rsimd_vf64_fma(a, c, rsimd_vf64_neg(th)),
                      rsimd_vf64_mul(a, rsimd_vf64_set1(ten ? RSIMD_MATH_LN10_LO : RSIMD_MATH_LN2_LO)));
  e = rsimd_sleef_expm1(th);
  r = rsimd_vf64_fma(tl, rsimd_vf64_add(one, e), e);
  r = rsimd_vf64_blend(r, th, rsimd_vf64_cmp_lt(rsimd_vf64_abs(a), rsimd_vf64_set1(0x1p-1000)));
  m = rsimd_vf64_cmp_gt(a, rsimd_vf64_set1(ten ? 17.0 : 54.0));
  if (rsimd_mf64_any(m)) r = rsimd_vf64_blend(r, ten ? rsimd_sleef_exp10(a) : rsimd_sleef_exp2(a), m);
  r = rsimd_vf64_blend(r, rsimd_vf64_neg(one), rsimd_vf64_cmp_lt(a, rsimd_vf64_set1(ten ? -330.0 : -1100.0)));
  /* Whole numbers with b^x - 1 exact. */
  m = rsimd_mf64_and(rsimd_mf64_and(rsimd_vf64_cmp_eq(a, rsimd_vf64_trunc(a)),
                                    rsimd_vf64_cmp_ne(a, rsimd_vf64_zero())),
                     rsimd_mf64_and(rsimd_vf64_cmp_le(a, rsimd_vf64_set1(ten ? 15.0 : 53.0)),
                                    rsimd_vf64_cmp_ge(a, rsimd_vf64_set1(ten ? 1.0 : -53.0))));
  if (rsimd_mf64_any(m)) {
    r = rsimd_vf64_blend(r, rsimd_vf64_sub(ten ? rsimd_sleef_exp10(a) : rsimd_sleef_exp2(a), one), m);
  }
  return r;
}

/* rsimd_log1p_base_f64() lane by lane, with SLEEF's log1p, log2 and
   log10 (exact at the powers of 2 and 10); the second is computed only
   for vectors with a lane where 1 + x is exact. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_log1p_base(rsimd_vf64 a, const int ten) {
  const rsimd_vf64 one = rsimd_vf64_set1(1.0);
  rsimd_vf64 u = rsimd_vf64_add(one, a), l, r;
  rsimd_mf64 exact = rsimd_mf64_and(rsimd_vf64_cmp_eq(rsimd_vf64_sub(u, one), a),
                                    rsimd_vf64_cmp_ne(a, rsimd_vf64_zero()));
  if (rsimd_mf64_all(exact)) return ten ? rsimd_sleef_log10(u) : rsimd_sleef_log2(u);
  l = rsimd_math_tiny(rsimd_sleef_log1p(a), a);
  r = rsimd_vf64_fma(l, rsimd_vf64_set1(ten ? RSIMD_MATH_INV_LN10 : RSIMD_MATH_INV_LN2),
                     rsimd_vf64_mul(l, rsimd_vf64_set1(ten ? RSIMD_MATH_INV_LN10_LO : RSIMD_MATH_INV_LN2_LO)));
  if (rsimd_mf64_any(exact)) {
    r = rsimd_vf64_blend(r, ten ? rsimd_sleef_log10(u) : rsimd_sleef_log2(u), exact);
  }
  return r;
}

/* The math1_f64 ops without a fast variant (EXP2M1 .. LOG10P1, which
   SLEEF has no 3.5-ULP expm1 and log1p for, and NEXT_UP .. RSQRT_APPROX):
   one instantiation. */
static int rsimd_math1_extra(int op, const void *x, R_xlen_t n, int flags, double *out) {
  const void *y = NULL;
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, 1), bc1 = rsimd_vf64_zero();
  const rsimd_vf64 one = rsimd_vf64_set1(1.0);
  int st = 0;
  (void) y;
  switch (op) {
  case RSIMD_MATH_EXP2M1: RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_expm1_base(a, 0))); break;
  case RSIMD_MATH_EXP10M1:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_expm1_base(a, 1)));
    break;
  case RSIMD_MATH_LOG2P1: RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_log1p_base(a, 0))); break;
  case RSIMD_MATH_LOG10P1:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_log1p_base(a, 1)));
    break;
  case RSIMD_MATH_NEXT_UP: RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_next_up(a))); break;
  case RSIMD_MATH_NEXT_DOWN:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_vf64_neg(rsimd_math_next_up(rsimd_vf64_neg(a)))));
    break;
  case RSIMD_MATH_RSQRT:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_vf64_div(one, rsimd_vf64_sqrt(a))));
    break;
  case RSIMD_MATH_RECIP_APPROX:
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_recip_approx(a)));
    break;
  case RSIMD_MATH_RSQRT_APPROX:
#if RSIMD_TIER_IS(neon)
    /* On Apple cores 1 / sqrt(a) is faster than the estimate and its two
       Newton steps (and exact). */
    if (rsimd_apple_core) {
      RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_vf64_div(one, rsimd_vf64_sqrt(a))));
      break;
    }
#endif
    RSIMD_MATH_LOOP(1, RSIMD_MATH1_RESULT(rsimd_math_rsqrt_approx(a)));
    break;
  default: break;
  }
  return st;
}

/* ---- ops handed to the C math library -----------------------------------
   The ops for which this tier runs the none tier's loops
   (rsimd_math1_libm, rsimd_math2_libm), and so libm as base R does,
   because libm was measured faster on this tier than SLEEF's 1-ULP vector
   function (issue 036, amending D13): bit k of RSIMD_MATH1_LIBM is
   math1_f64 op k, of RSIMD_MATH2_LIBM math2_f64 op k. A function belongs
   here only if libm wins on every core measured: on neon, Apple's (glibc
   and macOS's libm) and Neoverse N2's; on sve and sve2, which are 128-bit
   there, Neoverse N2's (glibc). Windows is not measured. LOGB is log(x) /
   log(base), so it follows LOG. Fast mode uses SLEEF's 3.5-ULP functions:
   only the ops SLEEF has none of (RSIMD_MATH*_LIBM_FAST) stay on libm
   there, since without one fast mode is accurate mode. Configure's
   RSIMD_NO_MATH_LIBM=1 empties the lists, to time SLEEF against libm. */
#if (RSIMD_TIER_IS(neon) || RSIMD_TIER_IS(sve) || RSIMD_TIER_IS(sve2)) && !defined(_WIN32) && \
  !defined(RSIMD_NO_MATH_LIBM)
#define RSIMD_MATH_BIT(op) ((uint64_t) 1 << (op))
#define RSIMD_MATH1_LIBM_FAST (RSIMD_MATH_BIT(RSIMD_MATH_ASINH) | RSIMD_MATH_BIT(RSIMD_MATH_ACOSH))
#define RSIMD_MATH1_LIBM                                                         \
  (RSIMD_MATH1_LIBM_FAST | RSIMD_MATH_BIT(RSIMD_MATH_LOG) | RSIMD_MATH_BIT(RSIMD_MATH_LOG2) | \
   RSIMD_MATH_BIT(RSIMD_MATH_LOGB) | RSIMD_MATH_BIT(RSIMD_MATH_COSH))
#define RSIMD_MATH2_LIBM RSIMD_MATH_BIT(RSIMD_MATH_POW)
#define RSIMD_MATH2_LIBM_FAST RSIMD_MATH2_LIBM
#else
#define RSIMD_MATH1_LIBM ((uint64_t) 0)
#define RSIMD_MATH1_LIBM_FAST ((uint64_t) 0)
#define RSIMD_MATH2_LIBM ((uint64_t) 0)
#define RSIMD_MATH2_LIBM_FAST ((uint64_t) 0)
#endif

int RSIMD_KERNEL(math1_f64)(int op, const void *x, R_xlen_t n, int flags, double p, double *out);
int RSIMD_KERNEL(math1_f64)(int op, const void *x, R_xlen_t n, int flags, double p, double *out) {
  if ((op & RSIMD_MATH_FAST ? RSIMD_MATH1_LIBM_FAST : RSIMD_MATH1_LIBM) >>
        (op & ~RSIMD_MATH_FAST) & 1u) {
    return rsimd_math1_libm(op & ~RSIMD_MATH_FAST, x, n, flags, p, out);
  }
  if ((op & ~RSIMD_MATH_FAST) >= RSIMD_MATH_EXP2M1) {
    return rsimd_math1_extra(op & ~RSIMD_MATH_FAST, x, n, flags, out);
  }
  if (op & RSIMD_MATH_FAST) return rsimd_math1_run(op & ~RSIMD_MATH_FAST, x, n, flags, p, out, 1);
  return rsimd_math1_run(op, x, n, flags, p, out, 0);
}

/* r = f(a, b) with math2's rule: NA if either operand is NA, else NaN if
   either is NaN; a NaN result from operands that are not NaN sets the
   status. */
#define RSIMD_MATH2_RESULT(e)                                                    \
  do {                                                                           \
    rsimd_mf64 rn_, in_;                                                         \
    r = (e);                                                                     \
    rn_ = rsimd_vf64_is_nan(r);                                                  \
    in_ = rsimd_mf64_or(rsimd_vf64_is_nan(a), rsimd_vf64_is_nan(b));             \
    if (rsimd_mf64_any(rsimd_mf64_or(rn_, in_))) {                               \
      if (rsimd_mf64_any(rsimd_mf64_andnot(in_, rn_))) st = RSIMD_EW_NAN_PRODUCED; \
      r = rsimd_vf64_na_merge(rsimd_vf64_blend(r, rsimd_vf64_set1(NAN), in_), a, b); \
    }                                                                            \
  } while (0)

/* math2_f64 with `fast` constant (pow has no fast variant). */
RSIMD_ALWAYS_INLINE int rsimd_math2_run(int op, const void *x, const void *y, R_xlen_t n,
                                        int flags, double *out, int fast) {
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, 1), bc1 = rsimd_ew_bcast(y, flags, 1, 1);
  int st = 0;
  switch (op) {
  case RSIMD_MATH_POW: RSIMD_MATH_LOOP(2, r = rsimd_math_pow(a, b, lanes)); break;
  case RSIMD_MATH_ATAN2:
    RSIMD_MATH_LOOP(2, RSIMD_MATH2_RESULT(RSIMD_SLEEF_CALL(atan2, a, b)));
    break;
  case RSIMD_MATH_HYPOT:
    RSIMD_MATH_LOOP(2, RSIMD_MATH2_RESULT(RSIMD_SLEEF_CALL(hypot, a, b)));
    break;
  default: break;
  }
  return st;
}

/* rsimd_ipow_dd() lane by lane, for whole an in [0, 512] (lanes with
   an = 0 give 1): binary powering until every lane's exponent is used. */
RSIMD_ALWAYS_INLINE void rsimd_math_ipow_dd(rsimd_vf64 r, rsimd_vf64 an, rsimd_vf64 *hi,
                                            rsimd_vf64 *lo) {
  const rsimd_vf64 zero = rsimd_vf64_zero(), half = rsimd_vf64_set1(0.5),
                   two = rsimd_vf64_set1(2.0);
  rsimd_vf64 ph = rsimd_vf64_set1(1.0), pl = zero, bh = r, bl = zero, h, l, nh, d;
  for (;;) {
    d = rsimd_vf64_floor(rsimd_vf64_mul(an, half));
    {
      rsimd_mf64 odd = rsimd_vf64_cmp_ne(an, rsimd_vf64_mul(d, two));
      if (rsimd_mf64_any(odd)) {
        h = rsimd_vf64_mul(ph, bh);
        l = rsimd_vf64_add(rsimd_vf64_fma(ph, bh, rsimd_vf64_neg(h)),
                           rsimd_vf64_add(rsimd_vf64_mul(ph, bl), rsimd_vf64_mul(pl, bh)));
        nh = rsimd_vf64_add(h, l);
        pl = rsimd_vf64_blend(pl, rsimd_vf64_sub(l, rsimd_vf64_sub(nh, h)), odd);
        ph = rsimd_vf64_blend(ph, nh, odd);
      }
    }
    an = d;
    if (!rsimd_mf64_any(rsimd_vf64_cmp_gt(an, zero))) break;
    h = rsimd_vf64_mul(bh, bh);
    l = rsimd_vf64_add(rsimd_vf64_fma(bh, bh, rsimd_vf64_neg(h)),
                       rsimd_vf64_mul(two, rsimd_vf64_mul(bh, bl)));
    bh = rsimd_vf64_add(h, l);
    bl = rsimd_vf64_sub(l, rsimd_vf64_sub(bh, h));
  }
  *hi = ph;
  *lo = pl;
}

/* rsimd_rootn_f64() lane by lane, with SLEEF's functions; b holds whole
   numbers. The general steps run in every lane and the special cases are
   blended over them. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_rootn(rsimd_vf64 a, rsimd_vf64 b) {
  const rsimd_vf64 zero = rsimd_vf64_zero(), one = rsimd_vf64_set1(1.0),
                   half = rsimd_vf64_set1(0.5), inf = rsimd_vf64_set1(INFINITY);
  rsimd_vf64 ax = rsimd_vf64_abs(a), an = rsimd_vf64_abs(b), k, m, r, c = zero;
  rsimd_mf64 small = rsimd_vf64_cmp_le(an, rsimd_vf64_set1(512.0)), odd, z, msk;
  k = rsimd_vf64_blend(zero, rsimd_vf64_floor(rsimd_vf64_div(rsimd_math_ilogb(ax), an)), small);
  m = rsimd_math_scaleb(ax, rsimd_vf64_neg(rsimd_vf64_mul(k, an)));
  r = rsimd_sleef_exp2_fast(rsimd_vf64_div(rsimd_sleef_log2_fast(m), b));
  if (rsimd_mf64_any(small)) {
    rsimd_vf64 h, l, cp, cn;
    rsimd_math_ipow_dd(r, rsimd_vf64_blend(zero, an, small), &h, &l);
    cp = rsimd_vf64_div(rsimd_vf64_add(rsimd_vf64_sub(h, m), l), h);
    cn = rsimd_vf64_sub(rsimd_vf64_fma(rsimd_vf64_neg(m), h, one), rsimd_vf64_mul(m, l));
    c = rsimd_vf64_div(rsimd_vf64_blend(cp, cn, rsimd_vf64_cmp_lt(b, zero)), b);
  }
  if (!rsimd_mf64_all(small)) {
    rsimd_vf64 h = rsimd_vf64_trunc(rsimd_vf64_mul(b, half)), q;
    q = rsimd_vf64_mul(rsimd_sleef_pow(r, h),
                       rsimd_vf64_div(rsimd_sleef_pow(r, rsimd_vf64_sub(b, h)), m));
    c = rsimd_vf64_blend(rsimd_vf64_div(rsimd_vf64_div(rsimd_vf64_sub(q, one), q), b), c, small);
  }
  r = rsimd_vf64_sub(r, rsimd_vf64_mul(r, c));
  r = rsimd_vf64_mul(
    r, rsimd_math_pow2i(rsimd_vf64_blend(k, rsimd_vf64_neg(k), rsimd_vf64_cmp_lt(b, zero))));
  /* b - 2 floor(b / 2) is exact for whole b. */
  odd = rsimd_vf64_cmp_ne(
    rsimd_vf64_sub(b, rsimd_vf64_mul(rsimd_vf64_set1(2.0), rsimd_vf64_floor(rsimd_vf64_mul(b, half)))),
    zero);
  r = rsimd_vf64_blend(r, rsimd_vf64_or(r, rsimd_math_signed_zero(a)), odd);
  msk = rsimd_vf64_cmp_eq(b, rsimd_vf64_set1(2.0));
  if (rsimd_mf64_any(msk)) r = rsimd_vf64_blend(r, rsimd_vf64_sqrt(a), msk);
  /* Zeros and infinities: 0 or Inf, signed like a for an odd b. */
  z = rsimd_vf64_cmp_eq(a, zero);
  msk = rsimd_mf64_or(z, rsimd_vf64_cmp_eq(ax, inf));
  if (rsimd_mf64_any(msk)) {
    rsimd_mf64 pos = rsimd_vf64_cmp_gt(b, zero);
    rsimd_mf64 to0 = rsimd_mf64_or(rsimd_mf64_and(z, pos), rsimd_mf64_andnot(rsimd_mf64_or(z, pos),
                                                                              msk));
    rsimd_vf64 mag = rsimd_vf64_blend(inf, zero, to0);
    r = rsimd_vf64_blend(r, rsimd_vf64_blend(mag, rsimd_vf64_or(mag, rsimd_math_signed_zero(a)), odd),
                         msk);
  }
  msk = rsimd_mf64_andnot(odd, rsimd_vf64_cmp_lt(a, zero));
  if (rsimd_mf64_any(msk)) r = rsimd_vf64_blend(r, rsimd_vf64_set1(NAN), msk);
  msk = rsimd_vf64_cmp_eq(an, one);
  if (rsimd_mf64_any(msk)) {
    r = rsimd_vf64_blend(r, rsimd_vf64_blend(a, rsimd_vf64_div(one, a), rsimd_vf64_cmp_lt(b, zero)),
                         msk);
  }
  return rsimd_vf64_blend(r, rsimd_vf64_set1(NAN), rsimd_vf64_cmp_eq(b, zero));
}

/* SLEEF's remainder, which is exact while |a / b| < 2^1000 but gives NaN
   once its quotient overflows (measured: exact up to exponents 1020
   apart): those lanes are recomputed by libm. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_remainder(rsimd_vf64 a, rsimd_vf64 b, int lanes) {
  rsimd_vf64 r = rsimd_sleef_remainder(a, b), ax = rsimd_vf64_abs(a);
  rsimd_mf64 big =
    rsimd_mf64_and(rsimd_vf64_cmp_gt(ax, rsimd_vf64_mul(rsimd_vf64_abs(b), rsimd_vf64_set1(0x1p1000))),
                   rsimd_vf64_cmp_lt(ax, rsimd_vf64_set1(INFINITY)));
  if (rsimd_mf64_any(big)) r = rsimd_math_fix2(r, a, b, big, lanes, remainder);
  return r;
}

/* The math2_f64 ops without a fast variant (NEXTAFTER .. ROOTN): one
   instantiation. SLEEF's nextafter is exact, and so is its remainder in
   the range rsimd_math_remainder() leaves it. */
static int rsimd_math2_extra(int op, const void *x, const void *y, R_xlen_t n, int flags,
                             double *out) {
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, 1), bc1 = rsimd_ew_bcast(y, flags, 1, 1);
  int st = 0;
  switch (op) {
  case RSIMD_MATH_NEXTAFTER:
    RSIMD_MATH_LOOP(2, RSIMD_MATH2_RESULT(rsimd_sleef_nextafter(a, b)));
    break;
  case RSIMD_MATH_REMAINDER:
    RSIMD_MATH_LOOP(2, RSIMD_MATH2_RESULT(rsimd_math_remainder(a, b, lanes)));
    break;
  case RSIMD_MATH_SCALEB: RSIMD_MATH_LOOP(2, RSIMD_MATH2_RESULT(rsimd_math_scaleb(a, b))); break;
  case RSIMD_MATH_ROOTN: RSIMD_MATH_LOOP(2, RSIMD_MATH2_RESULT(rsimd_math_rootn(a, b))); break;
  default: break;
  }
  return st;
}

int RSIMD_KERNEL(math2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                            double *out);
int RSIMD_KERNEL(math2_f64)(int op, const void *x, const void *y, R_xlen_t n, int flags,
                            double *out) {
  if ((op & RSIMD_MATH_FAST ? RSIMD_MATH2_LIBM_FAST : RSIMD_MATH2_LIBM) >>
        (op & ~RSIMD_MATH_FAST) & 1u) {
    return rsimd_math2_libm(op & ~RSIMD_MATH_FAST, x, y, n, flags, out);
  }
  if ((op & ~RSIMD_MATH_FAST) >= RSIMD_MATH_NEXTAFTER) {
    return rsimd_math2_extra(op & ~RSIMD_MATH_FAST, x, y, n, flags, out);
  }
  if (op & RSIMD_MATH_FAST) return rsimd_math2_run(op & ~RSIMD_MATH_FAST, x, y, n, flags, out, 1);
  return rsimd_math2_run(op, x, y, n, flags, out, 0);
}

/* rsimd_math_sinpi() and rsimd_math_cospi() of a, with one call of
   SLEEF's sincospi on the reduction of sinpi: the cosine of r is that of
   |r|, which is cospi's reduction, and SLEEF's sincospi computes the same
   as its sinpi and cospi, bit for bit. */
RSIMD_ALWAYS_INLINE void rsimd_math_sincospi(rsimd_vf64 x, rsimd_vf64 *vs, rsimd_vf64 *vc,
                                             int fast) {
  const rsimd_vf64 one = rsimd_vf64_set1(1.0), two = rsimd_vf64_set1(2.0);
  rsimd_vf64 r = rsimd_vf64_sub(x, rsimd_vf64_mul(two, rsimd_vf64_trunc(rsimd_vf64_mul(x, rsimd_vf64_set1(0.5)))));
  rsimd_vf64 s, c;
  r = rsimd_vf64_blend(r, rsimd_vf64_add(r, two), rsimd_vf64_cmp_le(r, rsimd_vf64_set1(-1.0)));
  r = rsimd_vf64_blend(r, rsimd_vf64_sub(r, two), rsimd_vf64_cmp_gt(r, one));
  if (fast) {
    rsimd_sleef_sincospi_fast(r, &s, &c);
  } else {
    rsimd_sleef_sincospi(r, &s, &c);
  }
  *vs = rsimd_vf64_blend(
    s, rsimd_math_signed_zero(x),
    rsimd_mf64_or(rsimd_vf64_cmp_eq(r, rsimd_vf64_zero()), rsimd_vf64_cmp_eq(r, one)));
  *vc = rsimd_vf64_blend(c, rsimd_vf64_zero(), rsimd_vf64_cmp_eq(rsimd_vf64_abs(r), rsimd_vf64_set1(0.5)));
}

/* sincos_f64 with `fast` and `pi` (op RSIMD_MATH_SINPI) constant. */
RSIMD_ALWAYS_INLINE int rsimd_sincos_run(const void *x, R_xlen_t n, int flags, double *s,
                                         double *c, int fast, int pi) {
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, 1);
  int st = 0;
  ptrdiff_t i = 0;
/* sin and cos (sinpi and cospi) of a into vs and vc, with math1's rule for
   each. */
#define RSIMD_MATH_SINCOS(a, vs, vc)                                             \
  do {                                                                           \
    rsimd_mf64 in_ = rsimd_vf64_is_nan(a);                                       \
    if (pi) {                                                                    \
      rsimd_math_sincospi(a, &vs, &vc, fast);                                    \
    } else if (fast) {                                                           \
      rsimd_sleef_sincos_fast(a, &vs, &vc);                                      \
    } else {                                                                     \
      rsimd_sleef_sincos(a, &vs, &vc);                                           \
    }                                                                            \
    if (rsimd_mf64_any(rsimd_mf64_andnot(                                        \
          in_, rsimd_mf64_or(rsimd_vf64_is_nan(vs), rsimd_vf64_is_nan(vc))))) {  \
      st = RSIMD_EW_NAN_PRODUCED;                                                \
    }                                                                            \
    vs = pi ? rsimd_vf64_blend(vs, a, in_) : rsimd_math_tiny(rsimd_vf64_blend(vs, a, in_), a); \
    vc = rsimd_vf64_blend(vc, a, in_);                                           \
  } while (0)
  for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
    rsimd_vf64 a = rsimd_ew_ld(x, flags, 0, bc0, i, 1), vs, vc;
    RSIMD_MATH_SINCOS(a, vs, vc);
    rsimd_vf64_storeu(s + i, vs);
    rsimd_vf64_storeu(c + i, vc);
  }
  if (i < n) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    rsimd_vf64 a = rsimd_ew_ld_p(x, flags, 0, bc0, i, pg, 1), vs, vc;
    RSIMD_MATH_SINCOS(a, vs, vc);
    rsimd_vf64_storeu_p(pg, s + i, vs);
    rsimd_vf64_storeu_p(pg, c + i, vc);
  }
#undef RSIMD_MATH_SINCOS
  return st;
}

int RSIMD_KERNEL(sincos_f64)(int op, const void *x, R_xlen_t n, int flags, double *s,
                             double *c);
int RSIMD_KERNEL(sincos_f64)(int op, const void *x, R_xlen_t n, int flags, double *s,
                             double *c) {
  if ((op & ~RSIMD_MATH_FAST) == RSIMD_MATH_SINPI) {
    if (op & RSIMD_MATH_FAST) return rsimd_sincos_run(x, n, flags, s, c, 1, 1);
    return rsimd_sincos_run(x, n, flags, s, c, 0, 1);
  }
  if (op & RSIMD_MATH_FAST) return rsimd_sincos_run(x, n, flags, s, c, 1, 0);
  return rsimd_sincos_run(x, n, flags, s, c, 0, 0);
}

#undef RSIMD_MATH_LOOP
#undef RSIMD_MATH1_FINISH
#undef RSIMD_MATH1_RESULT
#undef RSIMD_MATH1_FIXED
#undef RSIMD_MATH2_RESULT
#undef RSIMD_MATH_FIXED
#undef RSIMD_SLEEF_CALL
#undef RSIMD_MATH1_LIBM
#undef RSIMD_MATH1_LIBM_FAST
#undef RSIMD_MATH2_LIBM
#undef RSIMD_MATH2_LIBM_FAST
#undef RSIMD_MATH_BIT

#else /* a SIMD tier without SLEEF: the slots are filled from below */
#define RSIMD_SKIP_math1_f64 1
#define RSIMD_SKIP_math2_f64 1
#define RSIMD_SKIP_sincos_f64 1
#endif

/* ---- ilogb --------------------------------------------------------------
   C's ilogb as an int, NA (INT_MIN) for zero, an infinity, NaN or NA,
   where C returns its sentinels. Needs no SLEEF. */
#if RSIMD_TIER_IS(none)

void RSIMD_KERNEL(ilogb_f64)(const void *x, R_xlen_t n, int flags, int *out);
void RSIMD_KERNEL(ilogb_f64)(const void *x, R_xlen_t n, int flags, int *out) {
  R_xlen_t i;
  for (i = 0; i < n; i++) {
    double a = rsimd_math_get(x, flags, 0, i);
    out[i] = (isnan(a) || isinf(a) || a == 0) ? RSIMD_NA_I32 : ilogb(a);
  }
}

#elif !defined(RSIMD_NO_F64_SIMD)

/* The exponent as a double, NA lanes as INT_MIN (exactly a double), stored
   by the truncating conversion. */
RSIMD_ALWAYS_INLINE rsimd_vf64 rsimd_math_ilogb_na(rsimd_vf64 a) {
  rsimd_vf64 ax = rsimd_vf64_abs(a);
  rsimd_mf64 ok = rsimd_mf64_and(rsimd_vf64_cmp_gt(ax, rsimd_vf64_zero()),
                                 rsimd_vf64_cmp_lt(ax, rsimd_vf64_set1(INFINITY)));
  return rsimd_vf64_blend(rsimd_vf64_set1(-2147483648.0), rsimd_math_ilogb(ax), ok);
}

void RSIMD_KERNEL(ilogb_f64)(const void *x, R_xlen_t n, int flags, int *out);
void RSIMD_KERNEL(ilogb_f64)(const void *x, R_xlen_t n, int flags, int *out) {
  const rsimd_vf64 bc0 = rsimd_ew_bcast(x, flags, 0, 1);
  ptrdiff_t i = 0;
  for (; i + RSIMD_LANES_64 <= n; i += RSIMD_LANES_64) {
    rsimd_vf64_storeu_i32(out + i, rsimd_math_ilogb_na(rsimd_ew_ld(x, flags, 0, bc0, i, 1)));
  }
  if (i < n) {
    rsimd_p64 pg = rsimd_p64_while(i, n);
    rsimd_vf64_storeu_i32_p(pg, out + i,
                            rsimd_math_ilogb_na(rsimd_ew_ld_p(x, flags, 0, bc0, i, pg, 1)));
  }
}

#else /* RSIMD_NO_F64_SIMD: 32-bit ARM, the none tier does doubles */
#define RSIMD_SKIP_ilogb_f64 1
#endif

#undef RSIMD_MATH_PI
#undef RSIMD_MATH_PI_LO
#undef RSIMD_MATH_LN2
#undef RSIMD_MATH_LN2_LO
#undef RSIMD_MATH_LN10
#undef RSIMD_MATH_LN10_LO
#undef RSIMD_MATH_INV_LN2
#undef RSIMD_MATH_INV_LN2_LO
#undef RSIMD_MATH_INV_LN10
#undef RSIMD_MATH_INV_LN10_LO
